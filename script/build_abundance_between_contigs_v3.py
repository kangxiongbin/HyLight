#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
# from itertools import islice

__author__ = "Xiongbin Kang, Luo Xiao"


usage = """%prog python build_abundance_between_contigs.py -reads <reads >  -con <contigs> -p <percentage>

 This program is used to filter read overlaps(i.e.,duplicate or internal overlaps) and improve the result of extending
 
"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-reads', dest = 'reads', required=True, type=str, help='sequence')
    parser.add_argument('-con', dest = 'con', required=True, type=str, help='contigs')
    parser.add_argument('-out', dest = 'out', type=str, default = "sfoverlaps.out", help='output file')
    parser.add_argument('-iden', dest = 'iden', type=float, default = 0.95, help='Identity within overlap region')
    parser.add_argument('-p', dest = 'perc', type = float, default = 0.8, help='the threshold for picking up major overlap')

    args = parser.parse_args()

    if not (args.reads and args.con):
        print("Please input reads and contigs")
        parser.print_help()
        sys.exit()

    bin = os.path.split(os.path.realpath(__file__))[0]
    con = args.con
    iden = args.iden
    reads = args.reads
    perc = args.perc

    overlap_contigs = os.popen("minimap2 -t 30 -cx ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 %s %s 2> /dev/null \
    | python %s/filter_trans_ovlp_inline.py -len 70 -iden %s " %(con, con, bin, iden))

    os.system("minimap2 -t 30 -L --eqx -ax map-pb -cx ava-pb -Hk19 -Xw5 -m100 -g10000 --max-chain-skip 25 %s %s \
    2> /dev/null | cut -f1-9 > sam.file" %(con, con))
    
    os.system("minimap2 -t 30 -L --eqx -c --sr -DP --no-long-join -k 21 -w 11 -s 60 -m 30 -n 2 -A 4 -B 2 --end-bonus=100 25 %s %s \
    2> /dev/null | cut -f1-9 >> sam.file" %(con, reads))

    sam_file = os.popen("cat sam.file")
    os.system("rm sam.file")

    (snp, map_po, start_po, all_snp_read) = prpare_mutation(sam_file)
    mutation = mutation_re(snp, start_po, map_po, th = 5)

    # count overlap bewteen contigs for preparing to continues extend those contigs
    ovlcon = {} # all nodes in the contigs
    ovlconf = {} # the forward nodes of a contigs
    novlconf = defaultdict(list) # the number of forward nodes
    mutation2 = {}

    for line1 in overlap_contigs:
            line = line1.split( )
            
            if len(line) < 5 :
                continue

            fkey = ":".join(sorted([line[0], line[5]]))

            if fkey in mutation:
                fre = mutation[fkey]/int(line[5])
                
                if fre > 0.0025:
                    continue

            if line[1] == line[3]:
            
                key = ":".join([line[0], line[5]])

                if line[0] in ovlconf:
                
                    ovlconf[line[0]].append(key)
                    novlconf[line[0]] += 1
                else:
                    ovlconf[line[0]] = [key]
                    novlconf[line[0]] = 1

            else:
        
                key = ":".join([line[5], line[0]])
            
                if line[5] in ovlconf:
                
                    ovlconf[line[5]].append(key)
                    novlconf[line[5]] += 1
                else:
                    ovlconf[line[5]] = [key]
                    novlconf[line[5]] = 1

            ovlcon[key] = 1

  
    # filter the overlap file for next extend step
    sfoverlaps = os.popen("minimap2 -t 30 -c -x map-pb  -X -k 21 %s %s 2> /dev/null \
    | python %s/filter_trans_ovlp_inline.py -len 40 -iden %s -sfo " %(con, con, bin, iden))
	
    with open (args.out, "w") as out:
        for line in sfoverlaps:
            keys = line.split("\t")
            k1 = ":".join([keys[0], keys[1]])
            k2 = ":".join([keys[1], keys[0]])
            if (k1 in ovlcon) or (k2 in ovlcon):
                out.write(line)
             
def plus_dict(key, dict):
    
    dict[key] = dict[key] + 1 if key in dict else 1
    
    return(dict)

def plus_value(key, value, dict):
    
    if key in dict:
        dict[key].append(value)
                    
    else:
        dict[key] = [value]
    
    return(dict)
    
    
def prpare_mutation(sam_file):

        snp = defaultdict(list) # snp position and times
        map_po = defaultdict(list) # snp position and the read that support the snp
        all_snp_read = defaultdict(list) # all SNPs in one read (read id and snp id include snp position)
        start_po = defaultdict(list) # the start position of match reads
        deduplication = defaultdict(list) # delete the duplication of read combination

        ref_len = {}
        for line in sam_file:
            
            len_t = line.split('\t')

            if len(len_t) > 4:
                break
                
            ref = len_t[1].split(':')[1]
            ref_length = len_t[2].split(':')[1]
            ref_len[ref] = int(ref_length)
            
        for line in sam_file:
            len_t = line.split('\t')
            
            if len(len_t) < 6:
                continue

            [qseqid, flag, refid, refposi, quality, cigar] = line.split('\t')[:6]

            if qseqid == refid:
                continue
                
            if cigar == "*":
                continue

            k1 = ':'.join(sorted([qseqid, refid]))
            
            if k1 in deduplication:
                continue        
            else:
                deduplication[k1]
            
            ty1 =  re.split("\d+", cigar) # the array contains all type of mutation
            po1 = re.split("\D", cigar) # # the array contains the number of mutation (how many based in a mutation)
            
            if  int(flag) == 256:
                ty = ty1
                ty.pop(0)
            else:
                ty = ty1[::-1]
            
            if int(flag) == 256:
                po = po1

            else:

                po = po1[::-1]
                po.pop(0)
            
            pos1 = int(0) 
            
            if int(flag) == 256:
                pos2 = int(refposi) - 1
                
            else:
                pos2 = ref_len[refid] - int(refposi) + 1
                
                for t in range(len(ty)):
                    
                    if ty[t] == "=":
                        pos2 -= int(po[t])
                        
                    elif ty[t] == "D":
                        pos2 -= int(po[t])
                        
                    elif ty[t] == "X":
                        pos2 -= int(po[t])
                
            start_po = plus_value(refid, refposi, start_po)
            start_po = plus_value(qseqid, po[0], start_po) if ty[0] == "S" else plus_value(qseqid, "1", start_po)
                
            for t in range(len(ty)):
                
                if ty[t] == "=":
                    pos1 += int(po[t])
                    pos2 += int(po[t])
                    
                elif ty[t] == "S":
                    pos1 += int(po[t])

                elif ty[t] == "I":
                    pos1 += int(po[t])
                    
                elif ty[t] == "D":
                    pos2 += int(po[t])
                    
                elif ty[t] == "X":

                    pos1 += int(po[t])
                    pos2 += int(po[t])

                    key1 = ":".join([qseqid, str(pos1)])
                    key2 = ":".join([refid, str(pos2)])

                    snp = plus_dict(key1, snp)
                    
                    snp = plus_dict(key2, snp)
                    
                    map_po = plus_value(key1, refid, map_po)
                    
                    map_po = plus_value(key2, qseqid, map_po)

                    all_snp_read = plus_value(qseqid, key1, all_snp_read)

                    all_snp_read = plus_value(refid, key2, all_snp_read)
                    
        return(snp, map_po, start_po, all_snp_read)


def mutation_re(snp, start_po, map_po, th = 5):

    mutation = {}
    mutation_test = {}
    for k,v in snp.items():
        
        if int(v) < 3: # at least three reads support SNP
                continue

        (read, position) = k.split(":")
        position = int(position)
        
        con = 0
        for p in start_po[read]:
            
            if int(p) < position:
                con += 1
                
        c = con - int(v)  
        if c > int(th): # how many reads support the query reads base type
            
            for ke in map_po[k]:
                
                fkey = ':'.join(sorted([ke, read]))
                
                mutation = plus_dict(fkey, mutation)
#                mutation_test = plus_value(fkey, k, mutation_test)
 
    return(mutation)

           
if __name__ == '__main__':
        sys.exit(main())
