import sys
import os
import ast
import json
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import shutil

from parallel import *
from assembly import *
from utils import *
from toolkits import *
from clustering import limitedClusterReads, split_infiles_by_cluster
from logging import handlers

__author__ = "Xiongbin; Xiao Luo"
__version__ = "1.0.1"
__license__ = "GPL"
__date__ = 'November, 2023'

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
                        help="number of threads")
    parser.add_argument('-f', '--fast', dest='fast', action='store_true', required=False,
                        help="use fast parameters")
    parser.add_argument('-te', '--trim_ends', dest='trim_ends', action='store_true', required=False,
                        help="Trim the end of contigs")
    parser.add_argument('--overlap', dest='overlap', type=str, required=False, help="input file in PAF format")
    parser.add_argument('--corrected', dest='corrected', action='store_true', required=False, help="The long reads were corrected")
    parser.add_argument('--nsplit', dest='nsplit', type=int, required=False, default=60,
                        help="number of splitted input fasta/fastq files")
    parser.add_argument('--min_identity', dest='min_identity', type=float, required=False, default=0.95,
                        help="min identity for filtering overlaps")
    parser.add_argument('--min_read_len', dest='min_read_len', type=int, required=False, default=1000,
                        help="min read length for processing")
    parser.add_argument('--min_sread_len', dest='min_sread_len', type=int, required=False, default=3000,
                        help="min seed read length")
    parser.add_argument('--min_ovlp_len', dest='min_ovlp_len', type=int, required=False, default=3000,
                        help="min overlap length for super reads construction")
    parser.add_argument('--max_oh', dest='max_oh', type=int, required=False, default=40,
                        help="max overhang length")
    parser.add_argument('--sp_oh', dest='sp_oh', type=int, required=False, default=10,
                        help="max overhang length")
    parser.add_argument('--sp_min_identity', dest='sp_min_identity', type=float, required=False, default=0.995,
                        help="super reads min identity for filtering overlaps")
    parser.add_argument('--sp_min_ovlplen', dest='sp_min_ovlplen', type=int, required=False, default=3000,
                        help="min overlap length for super reads construction")
    parser.add_argument('--max_tip_len', dest='max_tip_len', type=int, required=False, default=1000,
                        help="max length to be removed as tips")
    parser.add_argument('--min_cluster_size', dest='min_cluster_size', type=int, required=False, default=2,
                        help="min size of read clusters")
    parser.add_argument('--level', dest='level', type=int, required=False, default=1,
                        help="iteration times of clustering read")
    parser.add_argument('--rm_trans', dest='rm_trans', type=int, required=False, default=1,
                        help="choose to (0) keep all edges, (1) remove transitive edges, (2) to remove double transitive edges")          
    parser.add_argument('--rename', dest='rename', type=ast.literal_eval, required=False, default=False,
                        help="rename read name or not, should be either True or False")
    parser.add_argument('--sort_by_len', dest='sort_by_len', type=ast.literal_eval, required=False, default=False,
                        help="sort super reads by read length or not(by number of overlaps), should be either True or False")
    parser.add_argument('--max_cluster_size', dest='max_cluster_size', type=int, required=False, default=1000,
                        help="max cluster size")
    parser.add_argument('--limited_times', dest='limited_times', type=int, required=False, default=1,
                        help="max times used for read")
    parser.add_argument('--max_ovlps', dest='max_ovlps', type=int, required=False, default=1000,
                        help="max number of overlaps for read")
    parser.add_argument('--min_cov', dest='min_cov', type=float, required=False, default=2.0,
                        help="min coverage for trimming consensus")
    parser.add_argument('--version', '-v', action='version', version='%(prog)s version: 1.0.1', help='show the version')
    
    args = parser.parse_args()
    
    if len(sys.argv[1:])==0:
        print(usage)
        parser.exit()

    global bin, len_over, max_cluster_size, fastx_files, min_cov, trim_ends
    
    bin = os.path.split(os.path.realpath(__file__))[0]
    
    overlap = args.overlap
    corrected = args.corrected
    short_reads = os.path.abspath(args.short_reads) 
    long_reads = os.path.abspath(args.long_reads)
#    short_reads = args.short_reads
#    long_reads = args.long_reads
    outdir = args.outdir
    nsplit = args.nsplit
    threads = args.threads 
    len_over = args.min_ovlp_len
    iden = args.min_identity 
    sort_by_len = args.sort_by_len
    min_sread_len = args.min_sread_len
    max_cluster_size = args.max_cluster_size
    level = args.level
    limited_times = args.limited_times
    max_ovlps = args.max_ovlps
    max_tip_len = args.max_tip_len
    rm_trans = args.rm_trans
    min_cluster_size = args.min_cluster_size
    trim_ends = args.trim_ends
    
    log = Logger('assembly_long_reads.log', level='debug')
    execute("mkdir -p %s" %outdir)

    log.logger.info('correct short reads with bfc')


    if corrected:
    
        infile = long_reads
        short_reads = short_reads

    else:

        cor_short = os.popen("bfc -s 3g -t%s %s 2> /dev/null" %(threads, short_reads))
        short_reads = outdir + "/cor_short_reads.fq"

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

        cmd1 = "cat {} | awk 'NR % 4 == 2' | sort | tr NT TN | {} -LR | tr NT TN | {} {}/comp_msbwt.npy;".format(short_reads, ropebwt2, convert, outdir)
        cmd2 = "cd %s; %s -t 30 -m 2  comp_msbwt.npy %s fmlrc1.fasta ; %s -t 30 -m 2 comp_msbwt.npy fmlrc1.fasta fmlrc2.fasta ; %s -t 30 -m 2 comp_msbwt.npy fmlrc2.fasta fmlrc3.fasta ; rm fmlrc1.fasta fmlrc2.fasta;" %(outdir, fmlrc2, long_reads,fmlrc2, fmlrc2)

        execute("mkdir -p %s" %outdir)
        execute(cmd1)
        execute(cmd2)

        infile = outdir + "/fmlrc3.fasta"
    
    log.logger.info('generate overlap information among long reads')
    
    model1 = fq_or_fa(infile)
    fastx_files = filter_non_atcg(infile, outdir, model1)
    infile = fastx_files
        
    # split reads and generate overlap file
    print('computing all-vs-all read overlaps...')

    log.logger.info('filter wrong overlap by SNPs and sort overlap by overlap score.')
    
    new_out_dir = outdir + "/2.overlap"
    os.system("rm -rf {}".format(new_out_dir))
   
    execute("mkdir -p %s" % new_out_dir)
    out_file = new_out_dir + "/s1_s1.paf"

    overlap = split_reads(fastx_files, fastx_files, nsplit, outdir, out_file, bin, threads = threads, len_over = 3000, mc = 2, iden = iden, long = True)

    # separate reads into different clusters
    print("Assign long reads into different clusters")
    fa = outdir+"/contigs1.fa"
    gfa = outdir+"/contigs1.gfa"
    execute("cd %s;miniasm -d 10000 -n 1 -e 1 -c 1 -f %s %s> %s;"%(outdir, infile, overlap, gfa))
    gfa2fa(gfa, fa)

    ref = outdir+"/ref_tmp.fa"
    with open(ref,'w') as ref1:
        with open(fa,'r') as con:
            for num, line in enumerate(con):
                if num % 2 == 1:
                    ref1.write(line)
                else:
                    new_id1 = str(int(num /2))
                    new_id = ">ref_" + new_id1 + "\n"
                    ref1.write(new_id)

    log.logger.info("Start assemble short reads")

    print("generate overlap between reads and temporary reference")
    
    ov_long_ref = outdir + "/ov_long_ref.paf"
    overlap2 = split_reads(infile, ref, nsplit, outdir, ov_long_ref, bin, threads = threads, len_over = 3000, mc = 2, iden = iden, long = True)
    remain_long = pick_up(overlap2, outdir, infile)

    ov_long_remain = outdir + "/ov_long_remain.paf"
    overlap3 = split_reads(infile, remain_long, nsplit, outdir, ov_long_remain, bin, threads = threads, len_over = 3000, mc = 2, iden = iden, long = True)
    
    remain_gfa = outdir+"/remain.gfa"
    execute("cd %s;miniasm -d 10000 -n 1 -e 1 -c 1 -f %s %s> %s;"%(outdir, infile, ov_long_remain, remain_gfa))
    remain_con = outdir+"/remain_con.fa"
    gfa2fa(remain_gfa, remain_con)

    # polish the long reads contigs
    ov_long_ref2 = outdir + "/ov_long_ref2.paf"
    overlap3 = split_reads(infile, remain_con, nsplit, outdir, ov_long_ref2, bin, threads = threads, len_over = 3000, mc = 2, iden = iden, long = True)
    
    p1 = outdir + "/polish1.fa"
    p2 = outdir + "/polish2.fa"
    execute("cd %s; racon -t 30 %s %s %s > polish1.fa; racon -t 30 %s %s %s > polish2.fa;" %(outdir, infile, ov_long_ref, ref, infile, ov_long_ref2, remain_con))
    
    long_con = outdir+"/long_con.fa"
    con_file = os.popen("cat %s %s" %(p1, p2))

    with open(long_con,'w') as longr_con:
        for num, line in enumerate(con_file):
                if num % 2 == 1:
                    longr_con.write(line)
                else:
                    new_id1 = str(int(num /2))
                    new_id = ">longr_con_" + new_id1 + "\n"
                    longr_con.write(new_id)

    execute("rm %s %s %s %s %s %s %s %s %s %s"%(gfa, fa, ref, ov_long_ref, ov_long_remain, remain_gfa, remain_con, ov_long_ref2, p1, p2))

    # pick up short reads
    ov_short = outdir + "/shortr1.paf"
    overlap4 = split_reads(short_reads, long_con, nsplit, outdir, ov_short, bin, threads = threads, len_over = 70, mc = 5, iden = iden)
    remain_short = pick_up(ov_short, outdir, short_reads)
    
    ov_short2 = outdir + "/shortr2.paf"
    overlap5 = split_reads(short_reads, remain_short, nsplit, outdir, ov_short2, bin, threads = threads, len_over = 70, mc = 5, iden = iden)

#    short_paf = generate_and_filter_overlap(short_reads, ref, nsplit, outdir, bin, threads = 30, len_over = 70, iden = iden)
#    remain_short = pick_up(short_paf, outdir, short_reads)
#    short_paf_remain = generate_and_filter_overlap(short_reads, remain_short, nsplit, outdir, bin, threads = 30, len_over = 70, iden = 0.9)
#    execute("rm %s %s" %(remain_long, remain_short))


def gfa2fa(gfa, fa):

    with open(gfa, 'r') as gfa_file, open(fa, 'w') as output_file:
    
        for line in gfa_file:
            fields = line.split()
            if fields[0] == "S":
                header = fields[1]
                sequence = fields[2]
                output_file.write(f'>{header}\n{sequence}\n')
                    


def get_fq4cluster(tar_clusters, file, outdir, read2cluster):
    for cluster in tar_clusters:
        outdir4c=outdir + "/" + cluster + "/"
        os.system("mkdir -p " + outdir4c)
        locals()["fq1w%s"%cluster] = open(outdir4c + cluster + '.1.fq', 'w')
        locals()["fq2w%s"%cluster] = open(outdir4c + cluster + '.2.fq', 'w')
    #
    pre_record = []
    pre_read = ''
    pre_header = ''
    with open(file, 'r') as fq:

        for i, line in enumerate(fq):
            if i%1000000 == 0: print('this is the: '+str(i)+' for 100w lines')
            if i % 4 == 0:
                if (i and (pre_read in read2cluster)) and (read2cluster[pre_read] in tar_clusters):
                    if re.search(r'/1$', pre_header):
                        locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
                    else:
                        locals()["fq2w%s"%read2cluster[pre_read]].writelines(pre_record)
                pre_record = []
                pre_record.append(line)
                pre_read = re.split('[@/]', line)[-2]
                pre_header = line
            else:
                pre_record.append(line)
                
        if (i and (pre_read in read2cluster)) and (read2cluster[pre_read] in tar_clusters):
            if re.search(r'/1$', pre_header):
                locals()["fq1w%s"%read2cluster[pre_read]].writelines(pre_record)
            else:
                locals()["fq2w%s"%read2cluster[pre_read]].writelines(pre_record)
    fq.close()
    for cluster in tar_clusters:
        locals()["fq1w%s"%cluster].close()
        locals()["fq2w%s"%cluster].close()

def get_merge_id(v1, v2, merged_id):
    cluster_id = ':'.join(sorted([v1, v2]))
    if cluster_id in merged_id:
        merged_id[cluster_id] = merged_id[cluster_id] + 1
    else:
        merged_id[cluster_id] = 1

def transfer_chain(key, cluster_remain_fq, clusters, clustersize):
    
    if key in cluster_remain_fq:
        keys = [k for k,v in cluster_remain_fq.items() if v == cluster_remain_fq[key]]

        for k in keys:
            clusters[k] = clusters[key]
            clustersize[clusters[key]] = clustersize[clusters[key]] +1
            del cluster_remain_fq[k]

def replace_id(k0, k1, clustersize, replace, maxsize): # k0 > k1
        
        kt = k0
        if k0 in replace:
            k0 = replace[k0]
        
        k2 = k1
        if k1 in replace:
            k1 = replace[k1]
            
        if(k0 != k1):
            
            if clustersize[k0] +  clustersize[k1] < maxsize:
            
                if clustersize[k0] < clustersize[k1]:
                    (k0, k1) = (k1, k0)
                replace[k2] = k0
                replace[k1] = k0

                if k1 in replace.values():
                    k3 = [k for k,v in replace.items() if v == k1]
                    for k4 in k3:
                        replace[k4] = k0
                    
                if k2 in replace.values():
                    k3 = [k for k,v in replace.items() if v == k2]
                    for k4 in k3:
                        replace[k4] = k0

                clustersize[k0] = clustersize[k0] + clustersize[k1]

def get_contigs4cluster(tar_clusters, file, outdir, contigs2cluster):
    for cluster in tar_clusters:
        outdir4c=outdir + "/" + cluster + "/"
        os.system("mkdir -p " + outdir4c)
        locals()["con1w%s"%cluster] = open(outdir4c + cluster + '.fa', 'w')

    pre_record = []
    pre_read = ''
    with open(file, 'r') as fa:

        for i, line in enumerate(fa):
            if i%1000000 == 0: print('this is the: '+str(i)+' for 100w lines')
            if i % 2 == 0:
                if (i and (pre_read in contigs2cluster)) and (str(contigs2cluster[pre_read]) in tar_clusters):
                    locals()["con1w%s"%contigs2cluster[pre_read]].writelines(pre_record)
                pre_record = []
                pre_record.append(line)
                line = line.strip('\n')
                pre_read = re.split('[>]', line)[-1]
            else:
                pre_record.append(line)

        if (i and (pre_read in contigs2cluster)) and (str(contigs2cluster[pre_read]) in tar_clusters):
            locals()["con1w%s"%contigs2cluster[pre_read]].writelines(pre_record)

    fa.close()
    for cluster in tar_clusters:
        locals()["con1w%s"%cluster].close()

def get_longsr4cluster2(tar_clusters, file, outdir, reads2cluster):
    for cluster in tar_clusters:
        outdir4c=outdir + "/" + cluster + "/"
        os.system("mkdir -p " + outdir4c)
        locals()["con1w%s"%cluster] = open(outdir4c + cluster + '_long.fa', 'w')

    pre_record = []
    pre_read = ''
    with open(file, 'r') as fq:

        for i, line in enumerate(fq):
            if i%1000000 == 0: print('this is the: '+str(i)+' for 100w lines')
            if i % 2 == 0:
                if (i and (pre_read in reads2cluster)) and (str(reads2cluster[pre_read]) in tar_clusters):
                    locals()["con1w%s"%reads2cluster[pre_read]].writelines(pre_record)
                pre_record = []
                pre_record.append(line)
                line = line.strip('\n')
                pre_read = re.split('[>]', line)[-1]
            else:
                pre_record.append(line)

        if (i and (pre_read in reads2cluster)) and (str(reads2cluster[pre_read]) in tar_clusters):
            locals()["con1w%s"%reads2cluster[pre_read]].writelines(pre_record)

    fq.close()
    for cluster in tar_clusters:
        locals()["con1w%s"%cluster].close()
        
def get_longsr4cluster(tar_clusters, file, outdir, reads2cluster):
    for cluster in tar_clusters:
        outdir4c=outdir + "/" + cluster + "/"
        os.system("mkdir -p " + outdir4c)
        locals()["con1w%s"%cluster] = open(outdir4c + cluster + '_long.fq', 'w')

    pre_record = []
    pre_read = ''
    with open(file, 'r') as fq:

        for i, line in enumerate(fq):
            if i%1000000 == 0: print('this is the: '+str(i)+' for 100w lines')
            if i % 4 == 0:
                if (i and (pre_read in reads2cluster)) and (str(reads2cluster[pre_read]) in tar_clusters):
                    locals()["con1w%s"%reads2cluster[pre_read]].writelines(pre_record)
                pre_record = []
                pre_record.append(line)
                line = line.strip('\n')
                pre_read = re.split('[@]', line)[-1]
            else:
                pre_record.append(line)

        if (i and (pre_read in reads2cluster)) and (str(reads2cluster[pre_read]) in tar_clusters):
            locals()["con1w%s"%reads2cluster[pre_read]].writelines(pre_record)

    fq.close()
    for cluster in tar_clusters:
        locals()["con1w%s"%cluster].close()

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
