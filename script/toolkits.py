import sys
import os 
import logging
from logging import handlers


def fq_or_fa(file):
    # with open(file) as fr:
        # s = fr.readline()[0]
    s = os.popen("less {}|head -1".format(file)).read()[0]
    mode = ''
    if s == '>':
        mode = 'fasta'
    elif s == '@':
        mode = 'fastq'
    else:
        raise Exception("invalid input file, must be FASTA/FASTQ format.", file)
    return mode
    
def fq2fa2(outdir, i):
    
    fq_p = outdir + "/%s_long.fq"%i
    fa_p = outdir + "/long.fa"
    with open(fa_p, "w" ) as fa_file:

        with open(fq_p) as r2:
            i = 1
            for line in r2:
                if i % 4 == 1:

                    line_new = line[1:]
                    fa_file.write('>'+line_new)

                elif i % 4 == 2:
                    newline = re.sub(r'[^ATGCN\n]', "N", line)
                    fa_file.write(newline)
                    
                i += 1 
                
def fq2fa(fq,fa):
    out=[]
    with open(fq,'r') as fr:
        c=0
        for line in fr:
            c+=1
            if (c%4) == 1:
                out.append('>' + line[1:])
            elif (c%4) == 2:
                out.append(line)
    with open(fa,'w') as fw:
        fw.write(''.join(out))
    return

def fa2fq(fa,fq):
    seq=''
    out=[]
    n_reads=0
    with open(fa,'r') as fr:
        for line in fr:
            if line.startswith('>'):
                n_reads+=1
                if seq:
                    score='I'*len(seq)
                    out.append(seq+'\n+\n'+score+'\n')
                    seq=''
                out.append('@'+line[1:])
            else:
                seq+=line.strip()
    #the last one
    score='I'*len(seq)
    out.append(seq+'\n+\n'+score+'\n')
    with open(fq,'w') as fw:
        fw.write(''.join(out))
    return n_reads

def get_read2seq(file, mode='fastq'):
    read2seq = {}
    read = ''
    seq = ''
    if mode == 'fastq':
        with open(file) as fr:
            for i, line in enumerate(fr):
                if i % 4 == 0:
                    read = line.strip().strip('@')
                elif i % 4 == 1:
                    seq = line.strip()
                    read2seq[read] = seq
                else:
                    continue
    elif mode == 'fasta':
        with open(file) as fr:
            for i, line in enumerate(fr):
                if i % 2 == 0:
                    read = line.strip().strip('>')
                elif i % 2 == 1:
                    seq = line.strip()
                    read2seq[read] = seq
                else:
                    continue
    else:
        print("Error: unknown mode, only fastq or fasta permitted.")
        sys.exit(1)
    return read2seq

def get_fasta(i, nodes, outdir,read2seq):
    '''
    :param nodes: reads ID list
    :return: fasta file
    '''
    prefix = str(i)
    fa = outdir + '/' + prefix + '.fa'
    with open(fa, 'w') as fw:
        for node in nodes:
            fw.write(">" + node + "\n" + read2seq[node] + "\n")
    return fa



def get_paf(i, reads, outdir,ovlp2record):
    '''
    generate paf file from given read ids
    '''
    # print(reads)
    paf = outdir + '/' + str(i) + '.paf'
    out = []
    for i in range(len(reads) - 1):
        for j in range(i + 1, len(reads)):
            key = reads[i] + ':' + reads[j]
            if key in ovlp2record:
                out.append('\t'.join(ovlp2record[key]))
            # else:
            # print('{} does not exist in overlap file!'.format(key))
    # print(list(ovlp2record.keys()))
    # print('paf out:')
    # print('\n'.join(out))
    with open(paf, 'w') as fw:
        fw.write('\n'.join(out))
    return paf


class Logger(object):
    level_relations = {
        'debug': logging.DEBUG,
        'info': logging.INFO,
        'warning': logging.WARNING,
        'error': logging.ERROR,
        'crit': logging.CRITICAL
    }

    def __init__(self, filename, level='info', when='D', backCount=3,
                 fmt='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'):
        self.logger = logging.getLogger(filename)
        format_str = logging.Formatter(fmt)
        self.logger.setLevel(self.level_relations.get(level))
        sh = logging.StreamHandler()
        sh.setFormatter(format_str)
        th = handlers.TimedRotatingFileHandler(filename=filename, when=when, backupCount=backCount, encoding='utf-8')
        th.setFormatter(format_str)
        self.logger.addHandler(sh)
        self.logger.addHandler(th)


def read_paf(paf):
    ovlp2record = {}  # store the record of v2w and w2v
    with open(paf, 'r') as fr:
        for line in fr:
            a = line.strip().split()
            node1 = a[0]
            node2 = a[5]
            key = node1 + ":" + node2
            key2 = node2 + ":" + node1
            ovlp2record[key] = a
            ovlp2record[key2] = ovlp2record[key]
    return ovlp2record
