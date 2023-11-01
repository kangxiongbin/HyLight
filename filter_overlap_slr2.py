#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
# from itertools import islice

__author__ = "Xiongbin Kang"


usage = """%prog python filter_overlap_slr.py -r <reads >  -c <contigs> -min_mis <int/mini mismatch> -t <int/threads>

 This program is used reads to polish contigs (refine racon)

"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-r', dest = 'reads', required=True, type=str, help='sequence')
    parser.add_argument('-c', dest = 'con', required=True, type=str, help='contigs')
    parser.add_argument('-t', dest = 'threads', required=False, type=int, default = 10, help='the number of threads')
    parser.add_argument('-pl', dest = 'platform', required=False, type=str, default = "pb", help='sequencing platform: pb/ont')
    parser.add_argument('-min_len', dest = 'min_len', required=False, type=int, default = 3000, help='the min length of overlap region')
    parser.add_argument('-min_mis', dest = 'min_mis', required=False, type=float, default = 0.02, help='the min mismatch in the overlap region')
    parser.add_argument('-long_reads', dest = 'long_reads', action='store_true', help='Long reads activte the parameter')
    
    args = parser.parse_args()

    if not (args.reads and args.con):
        print("Please input reads and contigs")
        parser.print_help()
        sys.exit()

    bin = os.path.split(os.path.realpath(__file__))[0]
    con = args.con
    reads = args.reads
    threads = args.threads
    min_len = args.min_len
    min_mis = args.min_mis
    
    platform = "pb" if args.platform == "pb" else "ont"

    overlap_long = "%s_tmp_overlap1.paf" %(con)
    overlap = "%s_tmp_overlap2.paf" %(con)

    if args.long_reads:
    
        execute("minimap2 -t %s -L --eqx -cx ava-%s -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 %s %s 2> /dev/null 1> %s" %(threads, platform, con, reads, overlap_long))
        execute("minimap2 -t %s -L --eqx -cx ava-%s -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 %s %s 2> /dev/null 1>> %s" %(threads, platform, reads, con, overlap_long))
    else:
    
        execute("minimap2 -t %s  -L --eqx -c --sr -DP --no-long-join -k 21 -w 11 -s 60 -m 30 -n 2 -A 4 -B 2 --end-bonus=100 %s %s 2> /dev/null 1> %s" %(threads, con, reads, overlap_long))
        
    paf_ovlp1 = os.popen("cat %s" %overlap_long)
    execute("rm %s"%overlap_long)

    dedup = {}
    with open (overlap, "w") as olvp_out:
            
            for line in paf_ovlp1:
                len_sp = line.split('\t')
                
                if int(len_sp[9]) < min_len:
                    continue

                mi = sum_before_X(len_sp[-1])
                mis = int(mi) / int(len_sp[9])

                if mis > min_mis:
                    continue

                k1 = ':'.join(sorted([len_sp[0], len_sp[5]]))
                if k1 in dedup:
                    continue 
                else:
                    dedup[k1] = 1

                [qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length] = len_sp[:11]
                
                mlen=(int(qlen)+int(slen))/2
                score = format((0.4*(int(length)/mlen)+0.6*(1 - mis)), '.3f')
                line2 ="\t".join([qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length, str(score), "\n"])

                olvp_out.write(line2)

def sum_before_X(string):
    digits_sum = 0
    for i in range(len(string)):
        if string[i] == "X":
            digits_sum += int(string[i-1])
    return digits_sum
    

def execute(cmd):
    te = os.system(cmd + " 2>out.txt")
    if te:
        with open("out.txt","r") as file:
            print("Don't execute the command: %s " %cmd)
            print(file.read())
      
if __name__ == '__main__':
        sys.exit(main())
