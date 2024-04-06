import sys, re, os, copy, gzip, time
from toolkits import *
from multiprocessing import Pool
from collections import defaultdict
from itertools import combinations, islice
from fnmatch import fnmatch
# log = Logger('whatshapdenovo.log', level='debug')

def split_reads(fa, ref, nsplit, out_dir, out_file, bin, threads = 30, len_over = 3000, mc = 2, iden = 0.95, long = False):

    # split contigs file
    nu = int(os.popen("wc -l " + ref ).readline().split()[0])
    nu_sub = int(nu/(8*nsplit)+1)*8
    sub = "sub"+str(time.time())[-3:]
    split_line = "cd %s; split %s -l %s -d -a 5 %s;" %(out_dir, ref, nu_sub, sub)
    execute(split_line)
    sub_overlap = [name for name in os.listdir(out_dir) if fnmatch(name, sub+"*")]
    th = 3 + int(threads/nsplit)
    
    if long:

        cmd = ["cd %s; python %s/filter_overlap_slr.py -r %s -c %s -t %s -len %s -mc %s -iden %s -long_reads; sort -k12  -T /tmp -nr %s_tmp_overlap4.paf > %s_sorted_overlap.paf;" %(out_dir, bin, fa, i, th, len_over, mc, iden, i, i) for i in sub_overlap]

    else:

        cmd = ["cd %s; python %s/filter_overlap_slr.py -r %s -c %s -t %s -len %s -mc %s -iden %s; sort -k12  -T /tmp -nr %s_tmp_overlap4.paf > %s_sorted_overlap.paf; " %(out_dir, bin, fa, i, th, len_over, mc, iden, i, i) for i in sub_overlap]

    cmd2 = "\n".join(cmd)
    with open('%s/cmd_overlap.sh' %out_dir, 'w') as fa:
        fa.write(cmd2)
        fa.write("\n")
        
    cmd_minimap = "cat %s/cmd_overlap.sh | xargs -i -P %s bash -c \"{}\";" %(out_dir, threads)
    
    execute(cmd_minimap) # run the minimap2 and get the overlap file

    execute("cd %s; sort -T /tmp -k12 -nr --parallel %s %s*_sorted_overlap.paf > %s; rm %s*;" %(out_dir, threads, sub, out_file, sub))

    return(out_file)

def split_reads2(fa, ref, nsplit, out_dir, out_file, bin, threads = 30, len_over = 3000, mc = 2, iden = 0.95, long = False):

    # split contigs file
    nu = int(os.popen("wc -l " + ref ).readline().split()[0])
    nu_sub = int(nu/(8*nsplit)+1)*8
    sub = "sub"+str(time.time())[-3:]
    split_line = "cd %s; split %s -l %s -d -a 5 %s;" %(out_dir, ref, nu_sub, sub)
    execute(split_line)
    sub_overlap = [name for name in os.listdir(out_dir) if fnmatch(name, sub+"*")]
    th = 3 + int(threads/nsplit)
    
    if long:

        cmd = ["cd %s; python %s/filter_overlap_slr2.py -r %s -c %s -t %s -len %s -mc %s -iden %s -long_reads; sort -k12  -T /tmp -nr %s_tmp_overlap4.paf > %s_sorted_overlap.paf;" %(out_dir, bin, fa, i, th, len_over, mc, iden, i, i) for i in sub_overlap]

    else:

        cmd = ["cd %s; python %s/filter_overlap_slr2.py -r %s -c %s -t %s -len %s -mc %s -iden %s; sort -k12  -T /tmp -nr %s_tmp_overlap4.paf > %s_sorted_overlap.paf; " %(out_dir, bin, fa, i, th, len_over, mc, iden, i, i) for i in sub_overlap]

    cmd2 = "\n".join(cmd)
    with open('%s/cmd_overlap.sh' %out_dir, 'w') as fa:
        fa.write(cmd2)
        fa.write("\n")
        
    cmd_minimap = "cat %s/cmd_overlap.sh | xargs -i -P %s bash -c \"{}\";" %(out_dir, threads)
    
    execute(cmd_minimap) # run the minimap2 and get the overlap file

    execute("cd %s; sort -T /tmp -k12 -nr --parallel %s %s*_sorted_overlap.paf > %s; rm %s*;" %(out_dir, threads, sub, out_file, sub))

    return(out_file)

def execute(cmd):
    te = os.system(cmd + " 2>output.txt")
    if te != 0:
        with open("output.txt", "r") as file:
            print("Error executing the command: %s" % cmd)
            print(file.read())
        sys.exit(1)

def filter_non_atcg(fq, out_dir, model):

    """Filter out non actg base"""

    i = 0
    new_out_dir = out_dir + "/1.split_fastx"
    execute("mkdir -p %s" %new_out_dir)
    new_out_fa = new_out_dir+"/s1.fa"

    with open(new_out_fa, 'w') as out1:
    
        if model == "fastq":
            with open(fq,'r') as file:
                for num, line in enumerate(file):
                    if num % 4 == 1:
                        newline = re.sub(r'[^ATGCN\n]', "N", line.upper())
                        out1.write(newline)
                    elif num % 4 == 0 and line.startswith('@'):
                        
                        new_id1 = line.strip().split(" ")[0]
                        new_id = ">" + new_id1[1:] + "\n"
                        out1.write(new_id)
        else:
            with open(fq,'r') as file:
                for num, line in enumerate(file):
                    if num % 2 == 1:
                        newline = re.sub(r'[^ATGCN\n]', "N", line.upper())
                        out1.write(newline)
                    else:
                        new_id1 = line.strip().split(" ")[0]
                        new_id = new_id1 + "\n"
                        out1.write(new_id)           
        
    return(new_out_fa)

## testing functions now##
if __name__ == '__main__':
    # preprocess_fastx(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])

    cluster_reads([sys.argv[1]], sys.argv[2], bool(sys.argv[3]), int(sys.argv[4]))
    # cluster_reads(ovlp_files, outdir, sort_by_len, k):
