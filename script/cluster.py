#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
from itertools import islice

__author__ = "Xiongbin Kang, Luo Xiao"


usage = """%prog -fq <input fq file>

This program is used to cluster the reads that stem from identical species.

# need to install panda

"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-fq', dest = 'fq', type = str, help = "The input fq file.")
    parser.add_argument('-t', dest = 'threads', default = 10, type = int, help = 'The number of threads when run strainxpress. Default is 10.')
    parser.add_argument('-size', dest = 'size', default = 15000, type = int, help = "The maximum size of the cluster. Default is 15000. ")
    parser.add_argument('-insert_size', dest = 'insert_size', default = 450, type = int, help = "The length of insert size of reads. Default is 450. ")
    parser.add_argument('-average_read_len', dest = 'average_read_len', default = 250, type = int, help = "The average length of reads. Default is 250. ")
    parser.add_argument('-split_nu', dest = 'split_nu', default = 8, type = int, help = "Split the fq file into several files. Default is 8.")
    parser.add_argument('-fast', dest = 'fast', action='store_true', help = "When fq file is big or have some high coverage bacteria, suggest use the fast model. Defualt is don't use fast.")
    args = parser.parse_args()
    
    (fq, threads, size, insert_size, average_read_len, split_nu) = (args.fq, args.threads, args.size, args.insert_size, args.average_read_len, args.split_nu)
    
    if not (fq):
        print("Specify min over length and min identical.")
        parser.print_help()
        sys.exit()

    id = "strainxpress"
    folder = os.getcwd()
    bin = os.path.split(os.path.realpath(__file__))[0]  
    k = os.path.split(os.path.realpath(__file__))[0]
    nu = int(os.popen("wc -l " + fq ).readline().split()[0])
    nu_sub = int((nu/(8*split_nu)+1)*8)
    folder_name = "fq_"+ str(size)

    split_line = "split {} -l {} -d -a 2 sub".format(fq, nu_sub)
    execute(split_line)

    cmd = ["minimap2 -t %s -c --sr -X -k 21 -w 11 -s 60 -m 30 -n 2 -r 0 -A 4 -B 2 --end-bonus=100 %s sub0%s 2> /dev/null | python %s/filter_trans_ovlp_inline_v3.py -len 30 -iden 0.90 -oh 1 > sub0%s.map " %(threads, fq, i, bin,i) for i in range(0, split_nu )]

    cmd2 = "\n".join(cmd)
    with open('cmd_overlap.sh', 'w') as fa:
        fa.write(cmd2)
        fa.write("\n")

    cmd_minimap = "cat cmd_overlap.sh | xargs -i -P %s bash -c \"{}\";" %threads
    execute(cmd_minimap) # run the minimap2 and get the overlap file

    if args.fast:
        execute("cat sub*.map > all_reads_sort.map")
    else:
        execute("for X in sub*.map; do sort  -k3 -nr < $X > sorted-$X; done;")
        execute("sort -k3 -nr -m sorted-sub*.map > all_reads_sort.map;")

    execute("rm *sub*;")    
    cmd_fq_name = "python %s/get_readnames.py %s readnames.txt" %(bin, fq) # get the name of reads
    execute(cmd_fq_name)

    cmd_cluster1 =  "python %s/bin_pointer_limited_filechunks_shortpath.py all_reads_sort.map readnames.txt %s %s %s" %(bin, size, id, threads)
    cmd_rm = "rm -rf Chunkfile*; rm %s_max%s_final_clustersizes.json %s_max%s_final_clusters_unchained.json %s_max%s_final_clusters.json" %(id, size, id, size, id, size)
    cmd_cluster2 = "python %s/getclusters.py %s_max%s_final %s" %(bin, id, size, threads)
    cmd_cluster3 = "python %s/get_fq_cluster.py %s_max%s_final_clusters_grouped.json %s %s/%s" %(bin, id, size, fq, folder, folder_name)
    execute(cmd_cluster1)
    execute(cmd_cluster2)
    execute(cmd_cluster3)
    execute(cmd_rm)
    

def execute(cmd):
    te = os.system(cmd + " 2>out.txt")
    if te:
        with open("output.txt","r") as file:
            print("Don't execute the command: %s " %cmd, end='')
            print(file.read())

        



if __name__ == '__main__':
        sys.exit(main())
