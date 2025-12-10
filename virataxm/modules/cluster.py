import os
import subprocess
import pandas as pd
import numpy as np
from itertools import combinations
from sklearn.cluster import HDBSCAN
import time


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="cluster")

    parser.add_argument("-o", "--output_dir", type=str, default="virataxm",
                        help="The directory of virataxm classification result files. Default: virataxm")

    parser.add_argument('-s', '--threads', type = int, default = 8,
                        help = "The number of CPU threads. Default: 8")


# 2 clustering by class, input with class label. cutoff, seq, protein, diamond file. output with cluster result in this class.
def clustering_by_class(df_c, c, cutoff, cluster_path):
    def run_mcl_clustering(class_path, inflation=2):
        output_prefix="mcl_output"

        # command
        input_file = f"{class_path}/{output_prefix}.abc"
        mci_file = f"{class_path}/{output_prefix}.mci"
        tab_file = f"{class_path}/{output_prefix}.tab"
        cluster_file = f"{class_path}/{output_prefix}.out"

        commands = [
            f"mcxload -abc {input_file} --stream-mirror -write-tab {tab_file} -o {mci_file}",
            f"mcl {mci_file} -I {inflation} -o {cluster_file}",
            f"mcxdump -icl {cluster_file} -tabr {tab_file} -o {class_path}/{output_prefix}_clusters.txt"
        ]

        for cmd in commands:
            # subprocess.run(cmd, shell=True, check=True)
            subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # parse the clustering result
        clusters = {}
        with open(f"{class_path}/{output_prefix}_clusters.txt") as f:
            for cluster_id, line in enumerate(f, 1):
                seqs = line.strip().split('\t')
                clusters[cluster_id] = seqs

        # remove tempotary file
        for f in [input_file, mci_file, tab_file, cluster_file]:
            if os.path.exists(f):
                os.remove(f)

        return clusters

    class_path = f"{cluster_path}/{c}"
    prot_file_c = f"{class_path}/prot.fa"
    diamond_db = f"{class_path}/prot.dmnd"
    diamond_out = f"{class_path}/diamond.tab"

    Seq_acc = df_c.index.tolist()

    # build_edge()
    dict_query_target = {}  # {query: {target: {gene_name : best_score}}}
    dict_best_score = {}
    pair_set = set()  # to avoid duplicate pairs
    with open(diamond_out, "r") as f:
        for line in f:
            prot_query, prot_target, idt ,_ ,_ ,_ ,_ ,_ ,_, _ , evalue, score = line.strip().split()
            query, target = prot_query.rsplit("_", 1)[0], prot_target.rsplit("_", 1)[0]
            idt, evalue, score = float(idt), float(evalue), float(score)
            dict_best_score[prot_query] = max(score, dict_best_score.get(prot_query, 0))

            if query not in dict_query_target:
                dict_query_target[query] = {}
            if target not in dict_query_target[query]:
                dict_query_target[query][target] = {}

            if prot_query not in dict_query_target[query][target]:
                dict_query_target[query][target][prot_query] = float(score)
                tmp_tuple = (query, target) if query < target else (target, query)
                pair_set.add(tmp_tuple)
            else:
                dict_query_target[query][target][prot_query] = max(float(score), dict_query_target[query][target][prot_query])

    # turn into mcl required file
    print(f"\t\tBuilding edges ({len(pair_set)} pairs)... ")
    abc_file = f"{class_path}/mcl_output.abc"
    with open(abc_file, 'w') as f:
        # df_sim = pd.DataFrame(0, index=Seq_acc, columns=Seq_acc, dtype=float)
        # df_score = pd.DataFrame(0, index=Seq_acc, columns=Seq_acc, dtype=int)

        # record the self alignment score / self-loop
        df_seq_self = {}
        for seq in Seq_acc:
            # df_sim.loc[seq, seq] = 1
            # self alignment
            if seq not in dict_query_target or seq not in dict_query_target[seq]:
                df_seq_self[seq] = 1
            else:
                df_seq_self[seq] = sum(dict_query_target[seq][seq].values())
            f.write(f"{seq}\t{seq}\t1\n")
            # df_score.loc[seq, seq] = df_seq_self[seq]


        # record the pairwise alignment score
        # pbar = tqdm(total=len(pair_set), desc="Calculating similarity scores")
        dict_pair_sim = {}
        for (seq_a, seq_b) in pair_set:
        # for (seq_a, seq_b) in combinations(Seq_acc, 2):
        #     pbar.update(1)
            score_a_b = score_b_a = 0
            # numerator
            if seq_a in dict_query_target and seq_b in dict_query_target[seq_a]:
                score_a_b = sum(dict_query_target[seq_a][seq_b].values())
            if seq_b in dict_query_target and seq_a in dict_query_target[seq_b]:
                score_b_a = sum(dict_query_target[seq_b][seq_a].values())

            if df_seq_self[seq_a] != 0 and df_seq_self[seq_b] != 0:
                sim = 0.5 * (score_a_b + score_b_a) / ((df_seq_self[seq_a] * df_seq_self[seq_b])**0.5)
                dict_pair_sim[(seq_a, seq_b)] = sim
                if sim > cutoff:  # 只写入正相似度
                    f.write(f"{seq_a}\t{seq_b}\t{sim}\n")


    # run mcl clustering
    print(f"\t\tRunning MCL clustering of Class - {c}... ")
    start_time = time.perf_counter()
    I, p = 2, 1
    # # df_sim = df_sim.applymap(lambda x: x ** p)

    clusters = run_mcl_clustering(class_path, inflation=I)
    df_cluster = pd.DataFrame(clusters.items(), columns=["Clust_id", "Seqs"])
    df_cluster["Clust_size"] = df_cluster["Seqs"].map(lambda x: len(x))
    # df_cluster.to_csv(f"{class_path}/cluster_result.csv", header=True, index=True)    # store the class-specific result

    df_cluster.index = df_cluster.index.map(lambda x: f"{c}_{x}" if not str(x).startswith(f"{c}_") else x)

    dict_seq2cluster = {seq: rid for rid, row in df_cluster.iterrows() for seq in row["Seqs"] }

    # double check each cluster for single genus confidence
    df_cluster["is_single_genus"] = None
    is_single_genus_clust = []
    not_single_genus_clust = []
    for rid, row in df_cluster.iterrows():
        seq_list = row["Seqs"]
        if len(seq_list)>=3:
            # check_confidence
            df_sim = pd.DataFrame(0, index=seq_list, columns=seq_list, dtype=float)
            np.fill_diagonal(df_sim.values, 1)
            # print(df_sim.values)
            for (seq_a, seq_b) in combinations(seq_list, 2):
                seq_a, seq_b = (seq_a, seq_b) if seq_a < seq_b else (seq_b, seq_a)
                df_sim.loc[seq_a, seq_b] = dict_pair_sim.get((seq_a, seq_b), 0)
                df_sim.loc[seq_b, seq_a] = df_sim.loc[seq_a, seq_b]
                # print((seq_a, seq_b), dict_pair_sim.get((seq_a, seq_b), 0))

            dist_mtx = 1 - df_sim.values
            dist_mtx = np.clip(dist_mtx, 0, 1)
            mask_num = df_sim[df_sim < cutoff].count().sum() / 2

            if mask_num > 0:
                clust_model = HDBSCAN(min_cluster_size=2, metric='precomputed', cluster_selection_method='eom')
                clusters = clust_model.fit_predict(dist_mtx)
                if 1 in clusters and len(seq_list) > 2:  # (-1 in clusters or 1 in clusters) and len(seqs) > 2:
                    not_single_genus_clust.append(rid)
                    continue
            is_single_genus_clust.append(rid)
        elif len(seq_list)==2:
            is_single_genus_clust.append(rid)

    df_cluster.loc[not_single_genus_clust, "is_single_genus"] = "low_confidence"
    df_cluster.loc[is_single_genus_clust, "is_single_genus"] = "high_confidence"
    # len(seq)=1 >> singleton
    df_cluster.to_csv(f"{class_path}/cluster_result.csv", header=True, index=True)

    end_time = time.perf_counter()
    print(f"\t\t\t clsutering take {end_time - start_time} seconds")
    return df_cluster, dict_seq2cluster



# 1 partition the classfied result and fasta file by class
def partition_by_class(df_unG, output_dir, threads):
    cluster_path = f"{output_dir}/class"
    os.makedirs(cluster_path, exist_ok=True)

    id2class = df_unG["pred|Class"].to_dict()
    for c in df_unG["pred|Class"].unique():
        class_path = f"{cluster_path}/{c}"
        os.makedirs(class_path, exist_ok=True)

        df_c = df_unG[df_unG["pred|Class"] == c]
        df_c.to_csv(f"{class_path}/df.csv")

    # partition seq file into class/{c}/seq.fasta
    fasta_file = f"{output_dir}/filtered_contigs.fa"
    file_handles = {}
    with open(fasta_file) as f:
        current_class = None
        current_handle = None
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].split()[0]
                if seq_id in id2class:
                    current_class = id2class[seq_id]
                    if current_class not in file_handles:
                        fh = open(os.path.join(f"{cluster_path}/{current_class}", f"seq.fasta"), "w")
                        file_handles[current_class] = fh
                    current_handle = file_handles[current_class]
                    current_handle.write(line)
                else:
                    current_class = None
                    current_handle = None
            else:
                if current_handle:
                    current_handle.write(line)
    for fh in file_handles.values():
        fh.close()

    # partition protein file into class/{c}/prot.fa
    prot_file = f"{output_dir}/prot.fa"
    file_handles = {}
    with open(prot_file) as f:
        current_class = None
        current_handle = None
        for line in f:
            if line.startswith(">"):
                # extract seq acc of the protein
                seq_id = line[1:].split()[0].rsplit("_", 1)[0]  # prot_id = {seq_id}_{protN}
                if seq_id in id2class:
                    current_class = id2class[seq_id]
                    if current_class not in file_handles:
                        fh = open(os.path.join(f"{cluster_path}/{current_class}", f"prot.fa"), "w")
                        file_handles[current_class] = fh
                    current_handle = file_handles[current_class]
                    current_handle.write(line)
                else:
                    # if seq not in df, then continue
                    current_class = None
                    current_handle = None
            else:
                if current_handle:
                    current_handle.write(line)
    for fh in file_handles.values():
        fh.close()

    # make diamond db + self alignment
    for c in df_unG["pred|Class"].unique():
        class_path = f"{cluster_path}/{c}"
        prot_file_c = f"{class_path}/prot.fa"
        diamond_db = f"{class_path}/prot.dmnd"

        if not os.path.exists(prot_file_c) or os.path.getsize(prot_file_c) == 0:
            # build empty output to avoid errors in downstream steps
            open(f"{class_path}/prot.fa", "w").close()
            open(f"{class_path}/prot.dmnd", "w").close()
            open(f"{class_path}/diamond.tab", "w").close()
            continue

        cmd = f"diamond makedb --in {prot_file_c} -d {diamond_db}"
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


        diamond_out = f"{class_path}/diamond.tab"
        cmd = f'diamond blastp -q {prot_file_c} ' \
                      f'-d {diamond_db} --threads {threads} ' \
                      f'-o {diamond_out} --top 90 --sensitive'
        subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)   #stdout=subprocess.DEVNULL,

    return cluster_path


def _cluster(output_dir, threads):

    threads = threads

    package_path = str(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    df_class_cutoff = pd.read_csv(f"{package_path}/genus_minimum_score_cutoff.csv", index_col=0)

    df_vtm = pd.read_csv(f"{output_dir}/vtm_result.csv", index_col=0)
    df_vtm["pred|Class"].fillna("Others", inplace=True)
    # df_vtm["pred|Genus"]=None  # for testing purpose
    df_ungenus = df_vtm[df_vtm["pred|Genus"].isna()]

    # 1) read, load, and partition
    print("1 - Partitioning sequences by class...")
    cluster_path = partition_by_class(df_ungenus, output_dir, threads)
    print("\tPartitioning done.")

    # 2) loop over each class, read the predefined cutoff, clustering, subcluster/single genus check
    print("2 - Clustering sequences by class...")
    df_cluster_all = pd.DataFrame()
    dict_seq2cluster = {}
    for c in df_ungenus["pred|Class"].unique():
        df_c = df_ungenus[df_ungenus["pred|Class"] == c]
        print(f"\tHandling class : {c} with {len(df_c)} sequences")
        genus_score_cutoff = df_class_cutoff.loc[c, "genus_minimum_score_cutoff"] if c in df_class_cutoff.index else 0.2
        df_cluster_c, dict_seq2cluster_c = clustering_by_class(df_c, c, genus_score_cutoff , cluster_path)
        df_cluster_all = pd.concat([df_cluster_all, df_cluster_c], axis=0)
        dict_seq2cluster.update(dict_seq2cluster_c)

    df_cluster_all["Clust_id"] = [i for i in range(len(df_cluster_all))]
    df_cluster_all.to_csv(f"{output_dir}/cluster_result.csv", header=True, index=True)

    # 3) merge classify and cluster result
    df_final = df_vtm.copy()
    df_final.loc[list(dict_seq2cluster.keys()), "pred|Genus"] = list(dict_seq2cluster.values())
    df_final.to_csv(f"{output_dir}/vtm_result_clustered.csv")


def main(args):
    output_dir = args["output_dir"]
    if not os.path.exists(output_dir):
        raise ValueError(f"The ViraTaxM result folder {output_dir} does not exist.")
    if not os.path.exists(output_dir + "/vtm_result.csv"):
        raise ValueError(f"The ViraTaxM result file {output_dir}/vtm_result.csv does not exist.")

    threads = args["threads"]

    _cluster(output_dir, threads)


if __name__ == "__main__":

    package_path = str(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    output_dir = f"{package_path}/test/virataxm"
    if not os.path.exists(output_dir):
        raise ValueError(f"The ViraTaxM result folder {output_dir} does not exist.")
    if not os.path.exists(output_dir + "/vtm_result.csv"):
        raise ValueError(f"The ViraTaxM result file {output_dir}/vtm_result.csv does not exist.")

    threads = 8

    _cluster(output_dir, threads)