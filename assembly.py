import os
import pysam
import networkx as nx

from utils import *
from overlap2graph import ovlp2graph
from toolkits import *
from consensus import get_err_consensus
from correction import *

# from correction import *
# log = Logger('whatshapdenovo.log', level='debug')


######################################################
################ For Raw Reads Assembly  #############
######################################################

def get_superead4(param):
    '''
    get the super reads for a cluster
    '''
    [i, outdir,  max_tip_len, rm_trans] = param
    i = str(i)
    outdir = outdir + '/c' + i

#    clog = Logger(outdir + '/' + i + '.log', level='debug')

    fasta = outdir + '/' + i + '_long.fa'
    paf = outdir + '/' + i + '.paf'
    os.system('ln -fs {} {}'.format(i + '_long.fa', outdir + '/' + i + '.corrected.fa'))

    if (os.path.getsize(fasta) == 0) or (os.path.getsize(paf) == 0):
        clog.logger.error("cluster {} is empty, skipping.}".format(i))
        return

    # get ad-hoc reference for calling variants
    ref = assemble_raw_reads(i, fasta, paf, outdir, max_tip_len, rm_trans)
    
    #polish the raw consensus
    
    ref = polish_seq2(i, ref, fasta, outdir)  # TODO: necessary when HiFi ??
    
    hap = "1"
    ref2 = scan_seq_by_depth(i, hap, ref, fasta, outdir)
   
#    os.system("cd %s; rm graph* overlaps* singles.fastq digraph.txt nonedge_overlaps.txt reads.id_map *renamed.fa  *corrected.fa *tmp.fa;" %outdir)
    os.system("cd %s; rm graph* overlaps* singles.fastq digraph.txt nonedge_overlaps.txt reads.id_map *corrected.fa *tmp.fa;" %outdir)
    
    return(ref2)

def get_superead3(param):
    '''
    get the super reads for a cluster
    '''
    [i, outdir,  max_tip_len, rm_trans] = param
    i = str(i)
    outdir = outdir + '/c' + i

    clog = Logger(outdir + '/' + i + '.log', level='debug')

    fasta = outdir + '/' + i + '_long.fa'
    paf = outdir + '/' + i + '.paf'
    os.system('ln -fs {} {}'.format(i + '_long.fa', outdir + '/' + i + '.corrected.fa'))

    if (os.path.getsize(fasta) == 0) or (os.path.getsize(paf) == 0):
        clog.logger.error("cluster {} is empty, skipping.}".format(i))
        return

    # get ad-hoc reference for calling variants
    ref = assemble_raw_reads(i, fasta, paf, outdir, max_tip_len, rm_trans)
    
    #polish the raw consensus
    
    ref = polish_seq2(i, ref, fasta, outdir)  # TODO: necessary when HiFi ??
    
#    os.system("cd %s; rm graph* overlaps* singles.fastq digraph.txt nonedge_overlaps.txt reads.id_map *renamed.fa  *corrected.fa;" %outdir)
    
    os.system("cd %s; rm graph* overlaps* singles.fastq digraph.txt nonedge_overlaps.txt reads.id_map *corrected.fa *tmp.fa;" %outdir)
    
    return

def assemble_raw_reads(id, fasta, paf, outdir, max_tip_len, rm_trans):
    # rename reads id, but this is not necessary when debug finished
    fasta2 = outdir + '/' + str(id) + '.renamed.fa'
    paf2 = outdir + '/' + str(id) + '.renamed.paf'
    id_map_file = rename_fa(fasta, fasta2, outdir)
    rename_paf(paf, paf2, id_map_file)
    digraph_file = ovlp2graph(fasta2, paf2, rm_trans, threads=1, remove_inclusions='false', rm_tips='true',
                              min_read_len=0, max_tip_len=max_tip_len)
    ovlp2record = read_paf(paf2)
    read2seq = get_read2seq(fasta2, 'fasta')
    G = construct_digraph(digraph_file)
    nodes = longest_path(G)

    ref_file = get_ctg_from_simple_path(id, nodes, ovlp2record, read2seq, outdir, as_file=True)
    return ref_file

def rename_fa(fasta, fasta2, outdir, min_len=0):
    id_map_file = outdir + '/reads.id_map'
    i = -1
    out = []
    id_map = []
    old_id = ''
    new_id = ''
    with open(fasta) as fr:
        for line in fr:
            if line.startswith('>'):
                old_id = line[1:].split()[0]
            else:
                if len(line.strip()) >= int(min_len):
                    i += 1
                    new_id = str(i)
                    out.append('>' + new_id + '\n' + line)
                    id_map.append(old_id + '\t' + new_id)
    with open(fasta2, 'w') as fw:
        fw.write(''.join(out))
    with open(id_map_file, 'w') as fw:
        fw.write('\n'.join(id_map))
    return id_map_file


def rename_paf(paf, paf2, id_map_file, min_ovlp_len=0):
    old2new = {}
    with open(id_map_file) as fr:
        for line in fr:
            a = line.strip().split()
            old2new[a[0]] = a[1]
    out = []
    with open(paf) as fr:
        for line in fr:
            a = line.strip().split()
            if (a[0] in old2new) and (a[5] in old2new) and (int(a[10]) >= min_ovlp_len):
                a[0] = old2new[a[0]]
                a[5] = old2new[a[5]]
                out.append('\t'.join(a))
    with open(paf2, 'w') as fw:
        fw.write('\n'.join(out))
    return

def construct_digraph(digraph_file):
    edges = []
    with open(digraph_file) as fr:
        for line in fr:
            a, b = line.strip().split()
            edges.append((a, b))
    G = nx.DiGraph()
    G.add_edges_from(edges)
    # Gr = nx.dag.transitive_reduction(G) #already removed using varialquasispecies program
    return G

def longest_path(G):
    '''
    used for merging raw reads into superead
    '''
    path = nx.dag_longest_path(G)  # ['1','3','2']
    return path

def get_ctg_from_simple_path(id, nodes, ovlp2record, read2seq, outdir, as_file=True):
    '''
    get supereads or contigs from a simple path[i.e., longest path/unitig path]
    '''
    # ovlp2record = read_paf("in.paf")  # OR read overlaps from paf file
    # ovlp2record = parse_paf_str(ovlp_str)
    # read2seq = get_read2seq(fasta, 'fasta')

    ctg = get_err_consensus(nodes, ovlp2record, read2seq)

    ref_file = outdir + "/" + str(id) + ".ref.fa"
    if as_file:
        with open(ref_file, 'w') as fw:
            fw.write(">" + str(id) + "\n")
            fw.write(ctg + "\n")
            return ref_file
    else:
        return ctg


def polish_seq2(i, ref, reads_fa, outdir):
    '''
    polish sequence using raw or corrected reads
    :param i: the i pivot read
    :param ref: fasta which needs to be polished
    :param reads_fa: reads used to polish
    :return a polished fasta file
    '''
    prefix = str(i)
    polish_paf = outdir + '/' + prefix + '.polish.paf'
    polished_fa = outdir + '/' + prefix + '.ref.polished.fa'
    tmp_fa = outdir + '/' + prefix + '.tmp.fa'
    os.system("cp {} {}".format(ref, tmp_fa))
    logfile = "{}/{}.log".format(outdir, i)

    os.system("minimap2 --secondary=no -x map-pb -c -t 1 {} {} 2>/dev/null |cut -f 1-12 >{}" .format(tmp_fa, reads_fa, polish_paf))
    try:
        os.system("racon  -t 1 {} {} {} >{} 2>{}".format(reads_fa, polish_paf, tmp_fa, polished_fa, logfile))
    except RuntimeError as reason:
        print("Error occurs when running Racon:\n" +
            "The reason is {}\n".format(reason) +
            "Maybe the cluster is not homogeneous, Racon can not support this!\n" +
            "Try compiling Racon on each machine individually if possible!",
            file=open("{}/{}.log".format(outdir, i), 'a'))
    os.system("cp {} {}".format(polished_fa, tmp_fa))  # the input and output fasta are identical now

    return polished_fa
    

if __name__ == '__main__':
    fasta, outdir, threads, min_read_len, min_ovlp_len, min_identity, o, r, max_tip_len = sys.argv[1:]
    assemble_supereads(fasta, outdir, threads, min_read_len, min_ovlp_len, min_identity, o, r, max_tip_len)
