#!/bin/sh

# $1 is contigs $2 is reads file; extend contig combine abundance information

#bin=/homes/kang/bin/hybrid/version2/script2
bin=/prj/metastrain/bin/git/StrainXpress/script2

#1.filter contigs with length <150bp
perl -ane'if(/^>/){if($. == 1){print"$F[0]\n";}else{print"\n$F[0]\n";};}else{chomp;s/[^ATCG]/C/g; print;}' $1 > contigs_b.fa

perl -e'open A, $ARGV[0] or die;my$i=0; $/=">";<A>;while(<A>){chomp;my@a=split/\n/;next if length($a[1])<150;$i+=1;print ">$i\n$a[1]\n";}' contigs_b.fa > contigs_b.fasta


#convert fasta to fastq
perl -e '$/=">";open I, $ARGV[0] or die;<I>;while(<I>){chomp;my@a=split/\n/;my$seq=join("",@a[1..$#a]);my$len=length($seq);print "@";print "$a[0]\n";print "$seq\n";print "+\n";print "=" x $len;print "\n";}' contigs_b.fasta > contigs_b.fastq

contig_counts=`cat contigs_b.fasta |grep "^>"|wc -l`

python $bin/build_abundance_between_contigs_v3.py -con contigs_b.fasta  -reads $2 -out sfoverlaps.out -p 0.90


#convert rust-overlap format to savage format
python $bin/sfo2overlaps.py --in sfoverlaps.out --out sfoverlap.out.savage --num_singles $contig_counts --num_pairs 0

mkdir -p fastq 
cp contigs_b.fastq ./fastq/singles.fastq

python $bin/pipeline_per_stage.py --no_error_correction --remove_branches true  --stage b --min_overlap_len 500 --min_overlap_perc 0 --edge_threshold 0.995 --overlaps sfoverlap.out.savage --fastq ./fastq --max_tip_len 1000 --num_threads 20

python $bin/fastq2fasta.py singles.fastq final_result.fa

rm contigs* graph* p* 




