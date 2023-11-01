import json
import sys
import os
from itertools import combinations, islice


# cluster reads
def limitedClusterReads(ovlp_files, outdir, sort_by_len, min_sread_len, min_cluster_size, k, limited_times, max_ovlps,
                        max_cluster_size):  # TODO: sort overlap file by overlap length,10,11 at first !!
    '''
    assign reads to limited clusters
    '''
    read2neigb = {}  # read->[read_len,neighbor reads]
    for ovlp_file in ovlp_files:
        print(ovlp_file)
        with open(ovlp_file, 'r') as fr:
            for line in fr:
                a = line.split()
                qname, qlen, tname, tlen = a[0], int(a[1]), a[5], int(a[6])
                if qname in read2neigb and len(read2neigb[qname][1]) <= max_ovlps:
                    read2neigb[qname][1].append(tname)
                else:
                    read2neigb[qname] = [qlen, [tname]]

                if tname in read2neigb and len(read2neigb[qname][1]) <= max_ovlps:
                    read2neigb[tname][1].append(qname)
                else:
                    read2neigb[tname] = [tlen, [qname]]

    print('reading overlap done...')
    used_reads = {read: 0 for read in read2neigb.keys()}
    clusters_file = outdir + "/clustered_reads.list"
    fw = open(clusters_file, 'w')
    if sort_by_len:
        # TODO: Really need to sort?
        for item in sorted(read2neigb.items(), key=lambda d: d[1][0], reverse=True):
            sread = item[0]
            if used_reads[sread] or read2neigb[sread][0] < min_sread_len:
                continue
            clustered_reads = {}
            clustered_reads[sread] = 1

            last_reads = {}
            for i in range(k):
                if i == 0:
                    for r in read2neigb[sread][1]:
                        if used_reads[r] >= limited_times:
                            continue
                        last_reads[r] = 1
                        clustered_reads[r] = 1
                else:
                    tmp_reads = {}
                    if len(last_reads):
                        for sread_tmp in last_reads.keys():
                            if used_reads[sread_tmp] >= limited_times:
                                continue
                            for r in read2neigb[sread_tmp][1]:
                                if used_reads[r] >= limited_times:
                                    continue
                                clustered_reads[r] = 1
                                tmp_reads[r] = 1
                            if len(clustered_reads) > max_cluster_size:
                                break
                        last_reads = tmp_reads
                if len(clustered_reads) > max_cluster_size:
                    break
            if len(clustered_reads) < min_cluster_size:
                continue

            # mark as used reads
            for read in clustered_reads.keys():
                used_reads[read] += 1
            fw.write(' '.join([r for r in clustered_reads.keys()]) + '\n')

    else:
        for sread in read2neigb.keys():
            if used_reads[sread] or read2neigb[sread][0] < min_sread_len:
                continue
            clustered_reads = {}
            clustered_reads[sread] = 1

            last_reads = {}
            for i in range(k):
                if i == 0:
                    for r in read2neigb[sread][1]:
                        if used_reads[r] >= limited_times:
                            continue
                        last_reads[r] = 1
                        clustered_reads[r] = 1
                else:
                    tmp_reads = {}
                    if len(last_reads):
                        for sread_tmp in last_reads.keys():
                            if used_reads[sread_tmp] >= limited_times:
                                continue
                            for r in read2neigb[sread_tmp][1]:
                                if used_reads[r] >= limited_times:
                                    continue
                                clustered_reads[r] = 1
                                tmp_reads[r] = 1
                            if len(clustered_reads) > max_cluster_size:
                                break
                        last_reads = tmp_reads
                if len(clustered_reads) > max_cluster_size:
                    break
            if len(clustered_reads) < min_cluster_size:
                continue

            # mark as used reads
            for read in clustered_reads.keys():
                used_reads[read] += 1
            fw.write(' '.join([r for r in clustered_reads.keys()]) + '\n')

    del read2neigb, used_reads
    return clusters_file

def split_infiles_by_cluster(fastx_files, ovlp_files, clusters_file, outdir, threads):
    read2clusters = {}
    with open(clusters_file) as fr:
        cluster_id = 0
        for line in fr:
            cluster_id += 1
            os.system("mkdir -p {}/c{}".format(outdir, cluster_id))
            for read in line.strip().split():
                if read in read2clusters:
                    read2clusters[read] = read2clusters[read] + ' ' + str(cluster_id)
                else:
                    read2clusters[read] = str(cluster_id)
    num_clusters = cluster_id

    # generate fasta for each cluster
    cutoff = 1000000
    num_reads = 0
    cluster2fa = {}

#    for fastx_file in fastx_files:
    with open(fastx_files) as fr:
            while True:
                read_item = list(islice(fr, 2))
                if not read_item:
                    break
                num_reads += 1
                read = read_item[0][1:].strip()
                if read in read2clusters:
                    for cluster_id in set(read2clusters[read].split()):
                        if cluster_id in cluster2fa:
                            cluster2fa[cluster_id].append(''.join(read_item))
                        else:
                            cluster2fa[cluster_id] = [''.join(read_item)]

                if num_reads > cutoff:
                    for cluster_id in cluster2fa.keys():
                        out_fa = outdir + "/c" + cluster_id + "/" + cluster_id + "_long.fa"
                        with open(out_fa, 'a') as fw:
                            fw.write(''.join(cluster2fa[cluster_id]))
                    num_reads = 0
                    cluster2fa = {}
    # write the last part
    for cluster_id in cluster2fa.keys():
        out_fa = outdir + "/c" + cluster_id + "/" + cluster_id + "_long.fa"
        with open(out_fa, 'a') as fw:
            fw.write(''.join(cluster2fa[cluster_id]))
    cluster2fa = {}

    # generate paf for each cluster
    cutoff = 50000000
    num_ovlps = 0
    k = 1000
    cluster2ovlps = {}
    for ovlp_file in ovlp_files:
        with open(ovlp_file) as fr:
            while True:
                next_n_lines = list(islice(fr, k))
                if not next_n_lines:
                    break
                num_ovlps += k
                for line in next_n_lines:
                    a = line.split()
                    qname, tname = a[0], a[5]
                    if (qname in read2clusters) and (tname in read2clusters):
                        for cluster_id in set(read2clusters[qname].split()).intersection(
                                set(read2clusters[tname].split())):
                            if cluster_id in cluster2ovlps:
                                cluster2ovlps[cluster_id].append(line)
                            else:
                                cluster2ovlps[cluster_id] = [line]
                if num_ovlps > cutoff:
                    for cluster_id in cluster2ovlps.keys():
                        out_paf = outdir + "/c" + cluster_id + "/" + cluster_id + ".paf"
                        with open(out_paf, 'a') as fw:
                            fw.write(''.join(cluster2ovlps[cluster_id]))
                    num_ovlps = 0
                    cluster2ovlps = {}
    # write the last part
    for cluster_id in cluster2ovlps.keys():
        out_paf = outdir + "/c" + cluster_id + "/" + cluster_id + ".paf"
        with open(out_paf, 'a') as fw:
            fw.write(''.join(cluster2ovlps[cluster_id]))

    cluster2ovlps = {}

    return num_clusters


if __name__ == '__main__':
    ovlp_files, outdir, sort_by_len, min_cluster_size, k = sys.argv[1:]
    max_ovlps = 90
    max_cluster_size = 30000000000000000
    limited_times = 10000000000000
    limitedClusterReads([ovlp_files], outdir, False, int(min_cluster_size), int(k), int(limited_times), max_ovlps,
                        max_cluster_size)
