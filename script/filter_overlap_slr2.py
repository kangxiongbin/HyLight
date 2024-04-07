#!/usr/bin/env python

import sys, os, re
from argparse import ArgumentParser
from collections import defaultdict
# from itertools import islice

__author__ = "Xiongbin Kang"


usage = """%prog python filter_overlap_slr.py -r <reads >  -c <contigs> -mc <int/mini coverage> -t <int/threads>

 This program is used reads to polish contigs (refine racon)

"""

def main():
    parser = ArgumentParser(description = usage)
    parser.add_argument('-r', dest = 'reads', required=True, type=str, help='sequence')
    parser.add_argument('-c', dest = 'con', required=True, type=str, help='contigs')
    parser.add_argument('-mc', dest = 'mc', required=False, type=int, default = 5, help='the mini coverage of one mutation (snp)')
    parser.add_argument('-t', dest = 'threads', required=False, type=int, default = 10, help='the number of threads')
    parser.add_argument('-len', dest = 'min_ovlp_len', required=False, default = 70, type = int, help='the shortest length of overlap')
    parser.add_argument('-thre', dest = 'threshold', required=False, default = 0.0025, type = float, help='allow the number of mutations in the overlap region between two reads')
    parser.add_argument('-iden', dest = 'iden', required=False, default = 0.95, type = float, help='the minimum identity in the overlap region between two reads')
    parser.add_argument('-oh', dest = 'min_o', required=False, default = 4, type = int, help='allow the number of over hang in the overlap region')
    parser.add_argument('-long_reads', dest = 'long_reads', action='store_true', help='Long reads activte the parameter')
    args = parser.parse_args()

    if not (args.reads and args.con):
        print("Please input reads and contigs")
        parser.print_help()
        sys.exit()

    bin = os.path.split(os.path.realpath(__file__))[0]
    con = args.con
    reads = args.reads
    mc = args.mc
    min_o = args.min_o
    threads = args.threads
    min_ovlp_len = args.min_ovlp_len
    threshold = args.threshold
    iden = args.iden

    overlap = "%s_tmp_overlap2.paf" %(con)
    overlap_sort = "%s_tmp_overlap3.paf" %(con)
    overlap2 = "%s_tmp_overlap4.paf" %(con)

    if args.long_reads:

        execute("minimap2 -N 40 -t %s -L --eqx -cx ava-pb -Hk19  -m100 -g10000 --max-chain-skip 25 %s %s 2> /dev/null | python %s/filter_trans_ovlp_inline_v4.py -len 30  -oh 3 > %s " %(threads, con, reads,  bin, overlap))

    else:
    
        execute("minimap2 -t %s  -L --eqx -c --sr -DP --no-long-join -k 21 -w 11 -s 60 -m 30 -n 2 -A 4 -B 2 --end-bonus=100 %s %s 2> /dev/null | python %s/filter_trans_ovlp_inline_v4.py -len 30  -oh 3 > %s " %(threads, con, reads, bin, overlap))
		
    execute("sort -nk7 -k8 -k9 -k5 %s > %s; rm %s; " %(overlap, overlap_sort, overlap))
    
    paf_file2 = os.popen("cat %s" %overlap_sort)
    
    if args.long_reads:
    
        (snp, map_po, start_po) = prpare_mutation2(paf_file2)
        
    else:
    
        (snp, map_po, start_po) = prpare_mutation(paf_file2)
 
    start_po_sorted = defaultdict(list) # read start and end
    for k in start_po.keys():
        start_po_sorted[k] = sorted(start_po[k], key = lambda x: (x[0], x[1]))

    mutation = mutation_re(snp, start_po_sorted, map_po, mc = mc)
	
    mutation2 = defaultdict(list)

    paf_file2 = os.popen("cat %s" %overlap_sort)
    dedup = {}
    
    with open (overlap2, "w") as out:

        for line1 in paf_file2:
            len_sp = line1.split('\t')

            [qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length] = len_sp[:11]
            (qlen, qstart, qend, slen, sstart, send, matchcount, length) = map(int, [qlen, qstart, qend, slen, sstart, send, matchcount, length])
            
            fkey = ":".join(sorted([qseqid, sseqid]))
            
            if args.long_reads:
                
                if fkey in mutation:
                    fre = mutation[fkey]/matchcount
                    if fre > threshold:
                    
                        continue
            else:
            
                if fkey in mutation:
                    continue                

            if qseqid == sseqid:
                continue

            if matchcount < min_ovlp_len:
                continue

            fkey = ":".join(sorted([qseqid, sseqid]))
            if fkey in mutation :
                fre = mutation[fkey]/matchcount
                if fre > threshold:
                    continue
                
            mlen=(qlen+slen)/2

            # filter overlap with overhang parameter
            if qori == '-':
                sstart2 = slen - send  # end in original seq
                send2 = slen - sstart  # start in original seq
            
                overhang = min(qstart, sstart2) + min(qlen - qend, slen - send2)
                maplen = max(qend - qstart, send2 - sstart2)

            else:
                # from algorithm 5 in minimap paper
                overhang = min(qstart, sstart) + min(qlen - qend, slen - send)
                maplen = max(qend - qstart, send - sstart)

            if overhang > min(min_o, maplen * 0.8):
                # internal match
                continue
            
            if fkey in dedup:
                continue
            else:
                dedup[fkey] = 1

            mi = sum_before_X(len_sp[-1])
            mis = int(mi)/ matchcount
            
            
            score = format((0.4*(matchcount/mlen)+0.6*(matchcount/length)), '.4f')
            score2 = format((1 - mis), '.4f')
            score3 = format((matchcount/length), '.4f')
            
            if(float(score2) < iden):
                continue
                
            (qlen, qstart, qend, slen, sstart, send, matchcount, length, score, score2, score3) = map(str, [qlen, qstart, qend, slen, sstart, send, matchcount, length, score, score2, score3])

            line2 ="\t".join([qseqid, qlen, qstart, qend, qori, sseqid, slen, sstart, send, matchcount, length, score, score2, score3,"\n"]) 
            out.write(line2)

    execute("rm %s;"%overlap_sort)

def sum_before_X(string):
    digits_sum = 0
    for i in range(len(string)):
        if string[i] == "X":
            digits_sum += int(string[i-1])
    return digits_sum
    

def pick_remain(polished_fa, con, remain_fa):

    d = defaultdict(list)
    
    with open(polished_fa) as f:
        i = 0
        for line in f:
            if i % 2 == 0 :
                line.replace('\n','')
                key = line.split( )[0]
                d[key] = 1

            i = i + 1
        
    with open(remain_fa,'a') as fw:
        with open(con) as q:
    
            i = 0
            for line2 in q:
                if i % 2 == 0 :
                    line2.replace('\n','')
                    key2 = line2.split( )[0]
                    if key2 in d:
                        p1 = 0
                    else:
                        p1 = 1

                if p1 > 0:
                    fw.write(line2)
                i = i + 1
    

def plus_dict(key, dict):
    
    dict[key] = dict[key] + 1 if key in dict else 1
    
    return(dict)

def plus_value(key, value, dict):
    
    if key in dict:
        dict[key].append(value)
                    
    else:
        dict[key] = [value]
    
    return(dict)

def plus_value2(key, value1, value2, dict):
    
    if key in dict:
        dict[key].append([value1, value2])
                    
    else:
        dict[key] = [[value1, value2]]
    
    return(dict)
    
def execute(cmd):
    te = os.system(cmd + " 2>out.txt")
    if te:
        with open("out.txt","r") as file:
            print("Don't execute the command: %s " %cmd)
            print(file.read())
        
def prpare_mutation(paf_file):

        snp = defaultdict(list) # snp position and times
        map_po = defaultdict(list) # snp position and the read that support the snp
        start_po = defaultdict(list) # the start and end position of match reads
        deduplication = defaultdict(list) # delete the duplication of read combination

        ref_len = {}

        for line in paf_file:
            len_sp = line.split('\t')
            
            if len(len_sp) < 6:
                continue

            qseqid = len_sp[0]
            flag = len_sp[4]
            refid = len_sp[5]
            refposi = len_sp[7]
            refposi2 = len_sp[8]
            cigar1 = len_sp[-1]

            if qseqid == refid:
                continue
                
            if cigar1 == "*":
                continue
                
            cigar = cigar1.strip()
            ty =  re.split("\d+", cigar) # the array contains all type of mutation
            po1 = re.split("\D", cigar) # # the array contains the number of mutation (how many based in a mutation)

            ty.pop(0)
            po = po1[5:-1]

            pos2 = int(refposi)

            start_po = plus_value2(refid, int(refposi), int(refposi2), start_po)
                
            for t in range(len(ty)):
                
                if ty[t] == "=":
                    pos2 += int(po[t])
 
                elif ty[t] == "D":
                    pos2 += int(po[t])
                    
                elif ty[t] == "X":

                    pos2 += int(po[t])

                    key2 = ":".join([refid, str(pos2)])

                    snp = plus_dict(key2, snp)

                    map_po = plus_value(key2, qseqid, map_po)
                    
                    
        return(snp, map_po, start_po)

def prpare_mutation2(paf_file):

        snp = defaultdict(list) # snp position and times
        map_po = defaultdict(list) # snp position and the read that support the snp
        start_po = defaultdict(list) # the start and end position of match reads
        deduplication = defaultdict(list) # delete the duplication of read combination

        ref_len = {}

        for line in paf_file:
            len_sp = line.split('\t')
            
            if len(len_sp) < 6:
                continue
                
            qseqid = len_sp[0]
            qseqlen = len_sp[1]
            qseqposi = len_sp[2]
            qseqposi2 = len_sp[3]
            flag = len_sp[4]
            refid = len_sp[5]
            refposi = len_sp[7]
            refposi2 = len_sp[8]
            cigar1 = len_sp[-1]
            cigar = cigar1.strip()
            
            if qseqid == refid:
                continue
                
            if cigar == "*":
                continue
                
            k1 = ':'.join(sorted([qseqid, refid]))
            
            if k1 in deduplication:
                continue        
            else:
                deduplication[k1] = 1
            
            ty =  re.split("\d+", cigar) # the array contains all type of mutation
            po1 = re.split("\D", cigar) # # the array contains the number of mutation (how many based in a mutation)

            ty.pop(0)
            po = po1[5:-1]
            
            pos1 = int(qseqposi) if flag == "+" else int(qseqlen) - int(qseqposi2)
            
            pos2 = int(refposi)

            start_po = plus_value2(refid, int(refposi), int(refposi2), start_po)
            start_po = plus_value2(qseqid, int(qseqposi), int(qseqposi2), start_po)

            for t in range(len(ty)):
                
                if ty[t] == "=":
                    pos1 += int(po[t])
                    pos2 += int(po[t])

                elif ty[t] == "I":
                    pos1 += int(po[t])
                    
                elif ty[t] == "D":
                    pos2 += int(po[t])
                    
                elif ty[t] == "X":

                    pos1 += int(po[t])
                    pos2 += int(po[t])

                    key1 = ":".join([qseqid, str(pos1)]) if flag == "+" else ":".join([qseqid, str(int(qseqlen) - pos1 + 1)])
                    key2 = ":".join([refid, str(pos2)])

                    snp = plus_dict(key1, snp)
                    snp = plus_dict(key2, snp)
                    
                    map_po = plus_value(key1, refid, map_po)
                    map_po = plus_value(key2, qseqid, map_po)                    

        return(snp, map_po, start_po)


def mutation_re(snp, start_po, map_po, mc = 5):

    mutation = {}

    for k,v in snp.items():
        
        if v < mc: # at least three reads support SNP
                continue

        (read, position) = k.split(":")
        position = int(position)
        
        con = 0
        for p in start_po[read]:
            
            if position > p[1]:
                continue

            elif position > p[0] and position < p[1] :
                con += 1

            elif position < p[0]:
                break
				
        c = con - v 

        if c < mc: # how many reads support the query reads base type
             continue
            
        for ke in map_po[k]:
                
              fkey = ':'.join(sorted([ke, read]))
                
              mutation = plus_dict(fkey, mutation)
 
    return(mutation)


           
if __name__ == '__main__':
        sys.exit(main())