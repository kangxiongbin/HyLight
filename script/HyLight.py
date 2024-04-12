#!/usr/bin/env python

import sys
import os
import ast
import json
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import shutil
 
from utils import *
from toolkits import *
from logging import handlers

__author__ = "Xiongbin Kang; Xiao Luo;"
__version__ = "1.0.1"
__license__ = "GPL"
__date__ = 'May, 2022'

usage = """Haplotype-aware de novo assembly of metagenome from hybrid sequencing data (long reads, short reads)"""

def main():
    parser = ArgumentParser(prog='python HyLight.py', description=usage, epilog='Built by: {}'.
                            format(__author__), formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--short_reads', dest='short_reads', type=str, required=False,
                        help="input short reads file in FASTQ format")
    parser.add_argument('-l', '--long_reads', dest='long_reads', type=str, required=True,
                        help="input long reads file in FASTQ format")
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, required=False, default='./',
                        help="output directory")
    parser.add_argument('-t', '--threads', dest='threads', type=int, required=False, default=20,
                        help="number of threads. Default is 20.")
    parser.add_argument('--corrected', dest='corrected', action='store_true', required=False,
                        help="The long and short reads had been corrected.")
    parser.add_argument('--low_q', dest='low_quality', action='store_true', required=False,
                        help="when the short reads is low quality and the long reads cannot be corrected perfectly.")
    parser.add_argument('--nsplit', dest='nsplit', type=int, required=False, default=60,
                        help="number of splitted input fasta/fastq files. Default is 60.")
    parser.add_argument('--min_identity', dest='min_identity', type=float, required=False, default=0.95,
                        help="min identity for filtering overlaps")
    parser.add_argument('--min_ovlp_len', dest='min_ovlp_len', type=int, required=False, default=3000,
                        help="min overlap length between long reads.")
    parser.add_argument('--size', dest = 'size', default = 15000, type = int, 
                        help = "The maximum size of the cluster for short reads. Default is 15000. ")
    parser.add_argument('--max_tip_len', dest='max_tip_len', type=int, required=False, default=10000,
                        help="max length to be removed as tips. Default is 10000.")
    parser.add_argument('--insert_size', dest = 'insert_size', default = 450, type = int, 
                        help = "The length of insert size of short reads. Default is 450. ")
    parser.add_argument('--average_read_len', dest = 'average_read_len', default = 250, type = int, 
                        help = "The average length of short reads. Default is 250. ")
                        
    parser.add_argument('--version', '-v', action='version', version='%(prog)s version: 1.0.1', help='show the version')
    
    args = parser.parse_args()
    
    if len(sys.argv[1:])==0:
        print(usage)
        parser.exit()


    global bin, len_over, max_cluster_size, fastx_files, miniasm
    
    bin = os.path.split(os.path.realpath(__file__))[0]

    miniasm = bin[:-7]+"/tools/miniasm/miniasm"

    short_reads = os.path.abspath(args.short_reads) 
    long_reads = os.path.abspath(args.long_reads)
    id = "HiStrain"
    (outdir, nsplit, threads, len_over, iden) = (args.outdir, args.nsplit, args.threads, args.min_ovlp_len, args.min_identity)
    (size, insert_size, average_read_len, corrected, max_tip_len) = (args.size, args.insert_size, args.average_read_len, args.corrected, args.max_tip_len)
    
    log = Logger('assembly_long_reads.log', level='debug')
    execute("mkdir -p %s" %outdir)
    execute("mkdir -p %s/tmp/" %outdir)
    outdir1 = outdir + "/tmp/"

#    log.logger.info('correct short reads with bfc')

    if corrected:
    
        infile = long_reads
        short_reads = short_reads

    else:

        cor_short = os.popen("bfc -s 3g -t%s %s 2> /dev/null" %(threads, short_reads))

#        cor_short = os.popen("cat %s" %(short_reads))

        short_reads = outdir1 + "/cor_short_reads.fq"

        with open(short_reads, "w") as outfile:
            for i, line in enumerate(cor_short):
                if (i + 1) % 4 == 1:
                    outfile.write(f"{line.split()[0]}\n")
                else:
                    outfile.write(line)
                
        log.logger.info('correcting long reads with short reads')

        fmlrc2 = shutil.which('fmlrc2') 
        ropebwt2 = shutil.which('ropebwt2')
        convert = shutil.which('fmlrc2-convert')

        cmd1 = "cat {} | awk 'NR % 4 == 2' | sort | tr NT TN | {} -LR | tr NT TN | {} {}/comp_msbwt.npy;".format(short_reads, ropebwt2, convert, outdir1)
        cmd2 = "cd %s; %s -t 30 -m 2  comp_msbwt.npy %s fmlrc1.fasta ; %s -t 30 -m 2 comp_msbwt.npy fmlrc1.fasta fmlrc2.fasta ; %s -t 30 -m 2 comp_msbwt.npy fmlrc2.fasta fmlrc3.fasta ; rm fmlrc1.fasta fmlrc2.fasta;" %(outdir1, fmlrc2, long_reads,fmlrc2, fmlrc2)

        execute(cmd1)
        execute(cmd2)

        infile = outdir1 + "/fmlrc3.fasta"
    
    log.logger.info('generate overlap information among long reads')
    
    model1 = fq_or_fa(infile)
    fastx_files = filter_non_atcg(infile, outdir, model1)
    infile = fastx_files

    # split reads and generate long reads overlap file

    log.logger.info('computing all-vs-all long read overlaps and then filter wrong overlap by SNPs and sort overlap by overlap score.')
    
    new_out_dir = outdir + "/2.overlap"
    os.system("rm -rf {}".format(new_out_dir))
   
    execute("mkdir -p %s" % new_out_dir)
    out_file = new_out_dir + "/s1_s1.paf"

    overlap = split_reads2(fastx_files, fastx_files, nsplit, outdir1, out_file, bin, threads = threads, len_over = 6000, mc = 2, iden = iden, long = True)

    log.logger.info("assemble long reads with clean graph")
    long_con = outdir1+"/contigs1.fa"
    gfa = outdir1+"/contigs1.gfa"
    
    if(args.low_quality):
        execute("cd %s; %s -d %s -n 3 -e 1 -c 3 -f %s %s> %s;"%(outdir1, miniasm, max_tip_len, infile, overlap, gfa))
        
    else:
        execute("cd %s; %s -d %s -n 1 -e 1 -c 1 -f %s %s> %s;"%(outdir1, miniasm, max_tip_len, infile, overlap, gfa))
 
    gfa2fa(gfa, long_con)

    log.logger.info("Start assemble short reads")

    log.logger.info("generate overlap between reads and temporary reference")
    
    ov_long_ref = outdir1 + "/ov_long_ref.paf"
    overlap2 = split_reads2(infile, long_con, nsplit, outdir1, ov_long_ref, bin, threads = threads, len_over = len_over, mc = 2, iden = iden, long = True)

    p1 = outdir1 + "/polish1.fa"
    execute("cd %s; racon --no-trimming -u -t 30 %s %s %s > %s;" %(outdir1, infile, ov_long_ref, long_con, p1))
    
    p2 = outdir1 + "/polish2.fa"
    
    ti = 0
    
    while(ti < 2):
        if(args.low_quality):

            con_file = os.popen("cat %s" %(p1))

        else:

            # utilize the remain long reads to generate contigs
            remain_long = pick_up(overlap2, outdir1, infile)
            ov_long_remain = outdir1 + "/ov_long_remain.paf"
            overlap3 = split_reads2(remain_long, remain_long, nsplit, outdir1, ov_long_remain, bin, threads = threads, len_over = len_over, mc = 2, iden = iden, long = True)

            remain_gfa = outdir1+"/remain.gfa"
            execute("cd %s;%s -d %s -n 1 -e 1 -c 1 -f %s %s> %s;"%(outdir1, miniasm, max_tip_len, infile, ov_long_remain, remain_gfa))

            if os.path.isfile(remain_gfa) and os.path.getsize(remain_gfa) == 0:
                break

            remain_con = outdir1+"/remain_con.fa"
            gfa2fa(remain_gfa, remain_con) 
    
            ov_long_ref2 = outdir1 + "/ov_long_ref2.paf"
            overlap3 = split_reads2(infile, remain_con, nsplit, outdir1, ov_long_ref2, bin, threads = threads, len_over = len_over, mc = 2, iden = iden, long = True)
    
            execute("cd %s; racon --no-trimming -u -t 30 %s %s %s >> %s; cat %s >> %s" %(outdir1, infile, ov_long_ref2, remain_con, p2, overlap3, overlap2))
        
            ti = ti + 1

    con_file = os.popen("cat %s %s" %(p1, p2))
    long_con2 = outdir1+"/long_con_polished.fa" 
    
    with open(long_con2,'w') as longr_con:
        for num, line in enumerate(con_file):
                if num % 2 == 1:
                    longr_con.write(line)
                else:
                    new_id1 = str(int(num /2))
                    new_id = ">longr_con_" + new_id1 + "\n"
                    longr_con.write(new_id)

    # pick up short reads for remaining strains
    ov_short = outdir1 + "/shortr1.paf"
    overlap4 = split_reads2(short_reads, long_con2, nsplit, outdir1, ov_short, bin, threads = threads, len_over = 70, mc = 3, iden = iden)

    long_con3 = outdir+"/long_con_polished.fa" 
    execute("cd %s; racon --no-trimming -u -t 30 %s %s %s > %s;" %(outdir, short_reads, ov_short, long_con2, long_con3))

    remain_short = pick_up(ov_short, outdir1, short_reads)
    ov_short2 = outdir1 + "/shortr2.paf"
    overlap5 = split_reads2(short_reads, remain_short, nsplit, outdir1, ov_short2, bin, threads = threads, len_over = 70, mc = 3, iden = iden)


    # assemble short reads

    folder_name = "/fq_"+ str(size)
    outdir2 = outdir1 + "/fq_"+ str(size)
    
    cmd_fq_name = "python %s/get_readnames.py %s %s/readnames.txt" %(bin, short_reads, outdir1) # get the name of reads
    execute(cmd_fq_name)

    cmd_cluster1 = "cd %s; python %s/bin_pointer_limited_filechunks_shortpath2.py %s readnames.txt %s %s %s" %(outdir1, bin, ov_short2, size, id, threads)
    cmd_cluster2 = "cd %s; python %s/getclusters.py %s_max%s_final %s ;" %(outdir1, bin, id, size, threads)
    cmd_cluster3 = "cd %s; python %s/get_fq_cluster.py %s_max%s_final_clusters_grouped.json %s %s/%s ;" %(outdir1, bin, id, size, short_reads, outdir1, folder_name)
    cmd_rm = "cd %s; rm -rf Chunkfile*; rm %s_max%s_final_clustersizes.json %s_max%s_final_clusters_unchained.json %s_max%s_final_clusters.json;" %(outdir1, id, size, id, size, id, size)
    
    execute(cmd_cluster1)
    execute(cmd_cluster2)
    execute(cmd_cluster3)
    execute(cmd_rm)

    fq_id = os.popen("ls %s" %(outdir2))

    for i in fq_id:
        i = i.strip()
        fq1 = "%s/%s/%s.1.fq" %(outdir2, i, i)
        fq2 = "%s/%s/%s.2.fq" %(outdir2, i, i)
        read_folder = "%s/%s" %(outdir2, i)

        cmd_polyte = "cd %s; python %s/polyte.tune_params.py -p1 %s -p2 %s  -m 50 -m_EC 60 -t 1 \
        --hap_cov 10 --insert_size  %s --stddev 27  --mismatch_rate 0 --min_clique_size 2 \
        --average_read_len %s --edge1 0.93 --edge2 1.0 --no_EC > log.txt 2>&1" \
        %(read_folder, bin, fq1, fq2, insert_size, average_read_len)

        with open('%s/cmd_polyte.sh' %outdir1, 'a') as fa:
            fa.write(cmd_polyte)
            fa.write("\n")
    
    cmd_polyte = "cat %s/cmd_polyte.sh | xargs -i -P %s bash -c \"{}\";" %(outdir1, threads)
    execute(cmd_polyte)
    
    fq_id = os.popen("ls %s" %(outdir2))

    for i in fq_id:
        i = i.strip()
        read_folder = "%s/%s/" %(outdir2, i)
        file_contigs = read_folder + "contigs.fasta"
        
        if not os.path.isfile(file_contigs):
            print("You need to rerun polyte in %s or decrease the size of cluster and rerun the whole steps" %read_folder)

    # collect all contigs from individual group.
    fname = "%s/all.contigs_%s.fasta" %(outdir1, size)
    cmd_merge_contigs = "cat %s/*/contigs.fasta > %s" %(outdir2, fname)
    execute(cmd_merge_contigs)

    # extend the length of short reads contigs
    short_con = outdir + "/short_stageb.fa"
    extend_con(fname, outdir1, short_con, bin, threads = 30, len_c = 500000)

    # extend the length of short and long reads contigs
    all_con = outdir + "/all_contigs.fa"

    if file_exists_and_nonempty(short_con):
    
        execute(f"cat {short_con} {long_con3} > {all_con}")
    else:

        execute(f"cat {fname} {long_con3} > {all_con}")

#    execute("cd %s; cat long_con_polished.fa short_stageb.fa > %s" %(outdir, all_con))

    out_file = outdir + "/final_contigs.fa"
    extend_con(all_con, outdir1, out_file, bin, threads = 30)
 
def extend_con(input_con, outdir1, out_file, bin, threads = 30, len_c = 50000000):

    #1.filter contigs with length < 150bp
    conb = '%s/contigs_b.fastq' %outdir1
    if os.path.exists(conb):
        execute("rm %s" %conb)

    if os.path.exists(input_con):
        with open(input_con, 'r') as f:
            n = 0
            for i in f:
                i = i.strip()
                if re.search(r"^>",i):
                    con_id = i
                else:
                    lenq = len(i)
                    if lenq > 150:
                        n +=1
                        seq = ("@"+str(n), i, "+", "="*lenq)
                        con_fq = "\n".join(seq)

                        with open(conb, 'a') as con_fq_w:
                            con_fq_w.write(con_fq)
                            con_fq_w.write("\n")                        
                            
    contig_counts = n
    execute("cd %s; mkdir -p stageb;" %outdir1)
    cmd_get_overlap = "cd %s/stageb; minimap2 -t %s --sr -X -c -k 21 -w 11 -s 60 -m 30 -n 2 -r 0 \
     -A 4 -B 2 --end-bonus=100 ../contigs_b.fastq ../contigs_b.fastq | python %s/filter_trans_ovlp_inline_v3.py \
      -len 90 -iden 0.99 -oh 2 -sfo > sfoverlaps.out;" %(outdir1, threads, bin)
    execute(cmd_get_overlap)
    
    # convert minimap format to savage format
    cmd_convert = "cd %s/stageb; python %s/sfo2overlaps.py --in sfoverlaps.out \
    --out sfoverlap.out.savage --num_singles %s --num_pairs 0; mkdir -p fastq;\
    cp ../contigs_b.fastq ./fastq/singles.fastq;" %(outdir1, bin, contig_counts)
    execute(cmd_convert)
    
    cmd_extend = "cd %s/stageb; python %s/pipeline_per_stage.py --no_error_correction --remove_branches true \
     --stage b --min_overlap_len 300 --min_overlap_perc 0 --edge_threshold 1 --overlaps ./sfoverlap.out.savage \
     --fastq ./fastq --max_tip_len 1000 --len_c %s --num_threads %s; python %s/fastq2fasta.py ./singles.fastq \
     %s;" %(outdir1, bin, len_c, threads, bin, out_file)
    execute(cmd_extend)

    return(out_file)

def gfa2fa(gfa, fa):

    with open(gfa, 'r') as gfa_file, open(fa, 'w') as output_file:
    
        for line in gfa_file:
            fields = line.split()
            if fields[0] == "S":
                header = fields[1]
                sequence = fields[2]
                output_file.write(f'>{header}\n{sequence}\n')

def file_exists_and_nonempty(filename):

    if os.path.isfile(filename) and os.path.getsize(filename) > 0:

        return True
    else:
        return False

def pick_up(ovlap, outdir, fq):
    
    d = {}
    sub = "sub"+str(time.time())[-3:]
    with open(ovlap) as f:
        for line in f:
            key = line.split( )
            k1 = key[0].split("/")
            k2 = key[5].split("/")
            d[k1[0]] = 1
            d[k2[0]] = 1

    out_put = outdir + "/" + sub + "_remain.fq"
    if os.path.exists(out_put):
        execute("rm -r %s" %out_put)

    mode = fq_or_fa(fq)
    nu = 4 if mode == 'fastq' else 2

    with open(fq,'r') as file:

        for num, line in enumerate(file):
            if num % nu == 0:
                kk1 = line.strip().split("/")
                if kk1[0][1:] in d:
                    p1 = 0
                else:
                    p1 = 1
            if p1 > 0:
                with open(out_put,'a') as fw:
                    fw.write(line)
    return out_put

if __name__ == '__main__':
    sys.exit(main())
