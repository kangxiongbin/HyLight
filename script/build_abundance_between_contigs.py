#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
from utils import *
from toolkits import *

__author__ = "Xiongbin Kang, Luo Xiao"


usage = """%prog python build_abundance_between_contigs.py -reads <reads >  -con <contigs>

 This program is used to filter read overlaps(i.e.,duplicate or internal overlaps) and improve the result of extending
 
"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-reads', dest = 'reads', required=True, type=str, help='short reads')
    parser.add_argument('-con', dest = 'con', required=True, type=str, help='require the format of contigs must be fastq')
    parser.add_argument('-out', dest = 'out', type=str, default = "sfoverlaps.out", help='output file')
    parser.add_argument('-iden', dest = 'iden', type=float, default = 0.95, help='Identity within overlap region')
    parser.add_argument('-nsplit', dest='nsplit', type=int, required=False, default=60,
                        help="number of splitted input fasta/fastq files. Default is 60.")
                        
    args = parser.parse_args()

    if not (args.reads and args.con):
        print("Please input reads and contigs")
        parser.print_help()
        sys.exit()

    con = args.con
    iden = args.iden
    reads = args.reads
    nsplit = args.nsplit
    outdir1 = os.getcwd()
    bin = os.path.split(os.path.realpath(__file__))[0]
    
    cr = outdir1 + "/cr.fq"
    os.system("cat %s %s > %s" %(reads, con, cr))

    rc_ovlp = outdir1 + "/rc_ovlp.paf"
    reads_cons_ovlp = split_reads(cr, con, nsplit, outdir1, rc_ovlp, bin, threads = 30, len_over = 70, mc = 4, iden = iden)

    rc_ovlp_popen = os.popen("cat %s" %rc_ovlp)
    with open (args.out, "w") as out:
        for line1 in rc_ovlp_popen:

            [qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length] = line1.split('\t')[:11]
            (qlen, qstart, qend, slen, sstart, send, matchcount, length) = map(int, [qlen, qstart, qend, slen, sstart, send, matchcount, length])
            
            if(qseqid.isdigit() and sseqid.isdigit()):
                line = minimap22sfo(qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length)
                out.write(line + "\n")

def minimap22sfo(qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length):
            ori = 'N' if qori == '+' else 'I'
            
            if ori == 'N':
                OHA = qstart - sstart
                OHB = slen - sstart - (qlen - qstart)
            else:
                OHA = qstart - (slen - send)
                OHB = send - (qlen - qstart)
                
            if OHA >= 0:
                OLA = min(qlen - OHA, slen)
            else:
                OLA = min(slen + OHA, qlen)
                
            OLB = OLA
            mismatch = length - matchcount
            
            sfo_line = '\t'.join([qseqid, sseqid, ori, str(OHA), str(OHB), str(OLA), str(OLB), str(mismatch)])
            return sfo_line

if __name__ == '__main__':
        sys.exit(main())
