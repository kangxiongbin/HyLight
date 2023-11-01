#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser
from collections import defaultdict
from itertools import islice

__author__ = "Xiongbin Kang"


usage = """python merge_paired_reads.py -r1 <The input read1: forward reads> -r2 <The input read2: reverse reads> -o <The output file>

 This program is used to forward and reverse reads (illumina)
 
"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-r1', dest = 'r1', type=str, help='The input read1: forward reads')
    parser.add_argument('-r2', dest = 'r2', type=str, help='The input read2: reverse reads')
    parser.add_argument('-o', dest = 'out', type=str, help='The output file')
    args = parser.parse_args()
    
    d = defaultdict(list)

    if len(sys.argv[1:])==0:
        print(usage)
        parser.exit()

    with open(args.r1) as f:
    
        for i, line in enumerate(f):
        
            if i % 4 == 0:
                pre_record = []
                pre_record.append(line)
                key = line.split("/")[0]
            else:
                pre_record.append(line)

            if i % 4 == 3:
                d[key[1:]] = pre_record
                
    with open(args.out,'w') as fw:
    
        with open(args.r2) as r:
        
            for i, line in enumerate(r):
                if i % 4 == 0 :
                    kk1 = line.split("/")[0]
                    if kk1[1:] in d:
                        p1 = 1
                        f1 = "".join(d[kk1[1:]])
                        fw.write(f1)
                    else:
                        p1 = 0
                        
                if p1 > 0:
                    fw.write(line)


def execute(cmd):
    te = os.system(cmd + " 2>out.txt")
    if te:
        with open("output.txt","r") as file:
            print("Don't execute the command: %s " %cmd, end='')
            print(file.read())
            

if __name__ == '__main__':
        sys.exit(main())
