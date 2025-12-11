#coding=utf-8

import os
import subprocess

from Bio import SeqIO
import pandas as pd
from collections import Counter


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="classify")
    # parser.set_defaults(description = "Classify the query viruses till genus level")

    required = parser.add_argument_group("required arguments")
    required.add_argument('-i', '--input', type = str, help = "The sequences of query viruses (.fasta)")
    ###### should be path
    required.add_argument('-o', '--output', type = str, default = ".",
                        help="The output path to store the classification results, will output to {output_path}/virataxm. Default: current working directory")


    optional = parser.add_argument_group("optional arguments")
    optional.add_argument('-r', '--ref', type=str, default=None,
                          help="The ViraTaxM reference database path. Default: built-in reference database")
    optional.add_argument('-g', '--gini', type = float, default = 0.44,
                        help = "The Gene Gini impurity cutoff to determine whether a genus assignment is confident. Range: (0, 1); Default: 0.44")
    optional.add_argument('-l', '--len', type=int, default=500,
                        help="The minimal length cutoff for the query sequences (bp). Default: 500 bp")
    optional.add_argument('--top', type = float, default = 10,
                        help = "Record protein alignments within this percentage range of top alignment score. Range: [0, 100]; Default: 10")
    optional.add_argument('-s', '--threads', type = int, default = 8,
                        help = "The number of CPU threads. Default: 8")

    # HIGH / LOW cutoff


def gini_impurity(input_list):
    if len(input_list) == 0:
        return 0
    counter = Counter(input_list)
    total_count = sum(counter.values())
    normalized_counts = {item: count / total_count for item, count in counter.items()}
    gini = 1 - sum(count ** 2 for count in normalized_counts.values())
    return gini


def get_topkGenus(input_dict, k=5):
    if len(input_dict) == 0:
        return []
    # from large to small
    topk_pairs = sorted(input_dict.items(), key=lambda x: x[1], reverse=True)[:k]  # [tuple]

    return topk_pairs

# 4 prediction strategey
def pred(df, df_taxon2lineage, high, low, gini_cutoff):
    def get_LCA(G_list):
        G_list_new = [x for x in G_list if pd.notna(x)]
        if len(G_list_new) == 0:
            return None
        elif len(G_list_new) == 1:
            return G_list_new[0]
        else:
            # 现在只用了是否hit，可以增加 ratio cutoff 的判断
            for rank in ["Family", "Order", "Class", "Phylum", "Kingdom"]:
                LCA = set()
                for G in G_list_new:
                    LCA.add(df_taxon2lineage.loc[G, rank])
                if len(LCA) == 1:
                    tmp = LCA.pop()
                    if pd.notna(tmp):
                        return tmp
        return None

    df["output|Genus"] = False
    df.loc[(df["G2_ratio"] < high) & (df["G2_ratio"] >= low) & (df["prot_G1_Gini"]<gini_cutoff), "output|Genus"] = True
    df.loc[(df["G2_ratio"] <low), "output|Genus"] = True
    df.loc[(df["G2"].isna()) & (df["G1_score"]>0), "output|Genus"] = True

    df["pred|final"] = df.apply(lambda x: x["G1"] if x["output|Genus"] else get_LCA(x[["G1", "G2", "G3"]].to_list()), axis=1)
    df["pred|rank"] = df["pred|final"].map(lambda x: df_taxon2lineage.loc[x, "Rank"] if pd.notna(x) else None)
    df.loc[~(df["G1"].isna()) & (df["pred|rank"].isna()), "pred|rank"] = "Higher than Kingdom"
    for rank in ["Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]:
        df[f"pred|{rank}"] = df["pred|final"].map(lambda x: df_taxon2lineage.loc[x, rank] if pd.notna(x) else None)
    df.drop(columns=["output|Genus", "pred|final"], inplace=True)

    return df

# 3 aggregate protein level information to sequence level
def aggregate_prot_to_seq(filtered_seq_list, dict_seq_len, df_prot_stat, dict_prot_genus_score):
    # 先 count : acc   prot_num	hit_num	prot_len	hit_len hit_len_ratio	tot_score	#score_vec
    df_seq_stat = pd.DataFrame(index=filtered_seq_list, columns=["seq_len", "prot_num", "hit_num", "prot_len", "prot_hit_len", "hit_len_ratio",
                                                            "tot_score",]) #"score_vec", "topkGenus"])
    df_seq_stat[["prot_num", "hit_num", "prot_len", "prot_hit_len", "hit_len_ratio"]] = 0
    df_seq_stat["tot_score"] = 0.0
    # df_seq_stat["score_vec"] = df_seq_stat["score_vec"].apply(lambda x: [], axis=1)

    df_prot_stat["seq"] = df_prot_stat.index.map(lambda x:x.rsplit("_", 1)[0])
    df_prot_stat["aligned_state"] = df_prot_stat["sum_cnt"].map(lambda x:x>0)
    df_protGroupByContig = df_prot_stat.groupby("seq")
    for seq, group in df_protGroupByContig: # seq:str, group:df
        df_seq_stat.loc[seq, "prot_num"] = len(group)  #df_protGroupByContig["seq"].size() #size
        df_seq_stat.loc[seq, "hit_num"] = group["aligned_state"].sum()   #df_protGroupByContig["aligned_state"].sum()
        df_seq_stat.loc[seq, "prot_len"] = group["len"].sum()  #df_protGroupByContig["len"].sum()  # size
        df_seq_stat.loc[seq, "prot_hit_len"] = group["hit_len"].sum()    #df_protGroupByContig["hit_len"].sum()
        df_seq_stat.loc[seq, "tot_score"] = group["bitscore"].sum()   #df_protGroupByContig["bitscore"].sum()

    df_seq_stat["hit_len_ratio"] = df_seq_stat.apply(lambda x:x["prot_hit_len"] / x["prot_len"] if x["prot_len"] > 0 else 0, axis=1)#df_seq_stat["prot_hit_len"] / df_seq_stat["prot_len"]

    # got ranking;
    dict_seq_genus_score = {seq: {} for seq in filtered_seq_list}
    for prot, dict_genus_score in dict_prot_genus_score.items():
        seq = prot.rsplit("_", 1)[0]
        for genus, score in dict_genus_score.items():
            if genus not in dict_seq_genus_score[seq]:
                dict_seq_genus_score[seq][genus] = score
            else:
                dict_seq_genus_score[seq][genus] += score

    df_seq_stat["topkGenus"] = df_seq_stat.index.map(lambda x: get_topkGenus(dict_seq_genus_score[x], k=3))

    for i in range(0, 3):
        df_seq_stat[f"G{i + 1}"] = df_seq_stat["topkGenus"].map(lambda x: x[i][0] if len(x) > i else None)
        df_seq_stat[f"G{i + 1}_score"] = df_seq_stat["topkGenus"].map(lambda x: x[i][1] if len(x) > i else None)
        df_seq_stat[f"G{i + 1}_ratio"] = df_seq_stat.apply(lambda x: x[f"G{i + 1}_score"] / x["G1_score"] if x["G1_score"]!=0 and pd.notna(x[f"G{i + 1}_score"]) else 0, axis=1)

    # prot_G1_list  prot_G1_Gini  #G1_hit_prot_num
    df_seq_stat["prot_G1_list"] = df_seq_stat.apply(lambda x: [], axis=1)
    df_seq_stat["prot_G1_score_list"] = df_seq_stat.apply(lambda x: {}, axis=1)
    for prot_acc, row in df_prot_stat.iterrows():
        seq, G1 = row["seq"], row["G1"]
        if row["aligned_state"]:
            df_seq_stat.loc[seq, "prot_G1_list"].append(G1)
            if G1 not in df_seq_stat.loc[seq, "prot_G1_score_list"]:
                df_seq_stat.loc[seq, "prot_G1_score_list"][G1] = 0
            df_seq_stat.loc[seq, "prot_G1_score_list"][G1] += row["G1_score"]

    df_seq_stat["prot_G1_Gini"] = df_seq_stat["prot_G1_list"].map(lambda x: gini_impurity(x))
    # df_seq_stat["prot_G1_score_Gini"] = df_seq_stat["prot_G1_score_list"].map(lambda x: gini_score_impurtiy(x))
    df_seq_stat["seq_len"] = df_seq_stat.index.map(lambda x: dict_seq_len[x] if x in dict_seq_len else 0)
    return df_seq_stat


    # generate cnt file using data-stream for faster processing
    # file0 = open(f"{output_dir}/prot_cnt_all.csv", "w")
    # file0.write("accession,prot_num,hit_num,prot_len,hit_len,hit_len_ratio,tot_score,score_vec\n")

# 2 parse the diamond blastp result and generate the taxonomic profile for proteins
def handle_prot(ref_label, prot_file, diamond_out, top_cutoff):
    """
    :param ref_label: the reference label file
    :param prot_file:
    :param diamond_out:
    :param cutoff: the minimum score cutoff to consider a hit as aligned
    :return: the taxonomic scoring profile of proteins and associated information
    """

    df_train = pd.read_csv(ref_label, index_col=0, header=0)
    df_diamond = pd.read_csv(diamond_out, index_col=False, header=None, delimiter="\t")

    prot_acc, prot_len = [], []
    with open(prot_file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            prot_acc.append(record.id)
            prot_len.append(len(record.seq))

    # record the diamond result
    df_diamond.columns = ["query", "subject", "identity", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    df_diamond["query_hit_len"] = df_diamond.apply(lambda x:abs(x["qend"] - x["qstart"]), axis=1)
    df_diamond["target_genus"] = df_diamond["subject"].map(lambda x: df_train.loc[x.rsplit("_", 1)[0], "Genus"])
    # df_diamond["query_genus"] = df_diamond["query"].map(lambda x: df_val.loc[x.rsplit("_", 1)[0], "Genus"])
    df_GroupByQuery = df_diamond.groupby("query")

    # apply the cutoff to the diamond result and group by query
    # max * percent (0.98)  # 只保留0.98的 alignment
    cutoff_percent = 1 - top_cutoff/100   # 10 > 90%
    df_prot_cutoff = df_GroupByQuery["bitscore"].max().to_frame()
    df_prot_cutoff.columns = ["bitscore"]
    df_prot_cutoff["cutoff"] = df_prot_cutoff["bitscore"] * cutoff_percent
    df_diamond["cutoff"] = df_diamond['query'].map(lambda x: df_prot_cutoff.loc[x, "cutoff"])
    df_diamond = df_diamond.loc[df_diamond["bitscore"]>=df_diamond["cutoff"]]

    df_GroupByQuery = df_diamond.groupby("query")  # renew the group above the cutoff
    df_prot_best = df_diamond.loc[df_GroupByQuery["bitscore"].idxmax()]
    df_prot_best.set_index("query", inplace=True)
    df_prot_genus = df_GroupByQuery["target_genus"].value_counts().to_frame()
    df_prot_genus.columns = ["count"]

    # get basic stat for each aligned protein
    df_prot_stat = pd.DataFrame(index=prot_acc, columns=["hit_genus_num", "hit_genus", "cnt", "sum_cnt"])
    df_prot_stat["hit_genus"] = df_prot_stat["hit_genus"].apply(lambda x: [])
    df_prot_stat["cnt"] = df_prot_stat["cnt"].apply(lambda x: [])
    df_prot_stat["len"] = prot_len

    for index, row in df_prot_genus.iterrows():
        # index[0]: query_acc; index[1]: target_genus;
        df_prot_stat.loc[index[0], "hit_genus"].append(index[1])
        df_prot_stat.loc[index[0], "cnt"].append(row["count"])
    df_prot_stat["sum_cnt"] = df_prot_stat["cnt"].apply(lambda x: sum(x))
    df_prot_stat["hit_genus_num"] = df_prot_stat["hit_genus"].apply(lambda x: len(x))
    # some query might not be aligned, and thus not exist in df_prot_best.index
    common_indices = df_prot_stat.index.intersection(df_prot_best.index)
    df_prot_stat[["evalue", "bitscore", "hit_idt", "hit_len"]] = [100.0, 0.0 ,0.0 ,0]    # default to a big value, means insignificant alignment
    df_prot_stat.loc[common_indices, "evalue"] = df_prot_best.loc[common_indices, "evalue"]
    df_prot_stat.loc[common_indices, "bitscore"] = df_prot_best.loc[common_indices, "bitscore"]
    df_prot_stat.loc[common_indices, "hit_idt"] = df_prot_best.loc[common_indices, "identity"]
    df_prot_stat.loc[common_indices, "hit_len"] = df_prot_best.loc[common_indices, "query_hit_len"]

    df_prot_stat["hit_ratio"] = df_prot_stat["hit_len"] / df_prot_stat["len"]

    # record the alignment dict of each protein, preparing for the ranking
    dict_prot_genus_score = {prot:{} for prot in prot_acc}
    df_GroupByQueryNGenus = df_diamond.groupby(["query", "target_genus"])
    df_GroupByQueryNGenus_score = df_GroupByQueryNGenus["bitscore"].max()
    # index is a tuple (query, target_genus), value is the max bitscore
    for index, score in df_GroupByQueryNGenus_score.items():
        dict_prot_genus_score[index[0]][index[1]] = score
    # print(dict_prot_genus_score)

    # structure [topkGenus, G1~5, G1~5_score]
    # 怎么从  dict_prot_genus_score >> G1 G1_score (prot 不需要 ratio)?
    # df_blast_score = pd.DataFrame(index=prot_acc, columns=["topkGenus"])
    df_prot_stat["topkGenus"] = None

    df_prot_stat["topkGenus"] = df_prot_stat.index.map(lambda x:get_topkGenus(dict_prot_genus_score[x], k=5))

    for i in range(0, 5):   # G1 ~ Gk
        df_prot_stat[f"G{i + 1}"] = df_prot_stat["topkGenus"].map(lambda x: x[i][0] if len(x) > i else None)
        df_prot_stat[f"G{i + 1}_score"] = df_prot_stat["topkGenus"].map(lambda x: x[i][1] if len(x) > i else None)


    return df_prot_stat, dict_prot_genus_score

# 1 filter by length, translate with prodigal, run diamond blastp against ref DB
def preprocessing_data(input_data, ref_path, output_path, top_cutoff, length, threads):

    # Check folders
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
    output_dir = output_path + "/virataxm"
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Filter short contigs
    filtered_seq_list = []
    rec = []
    dict_seq_len = dict()
    for record in SeqIO.parse(input_data, 'fasta'):
        if len(record.seq) > length:
            rec.append(record)
            filtered_seq_list.append(record.id)
            dict_seq_len[record.id] = len(record.seq)
    if not os.path.exists(f'{output_dir}/filtered_contigs.fa'):
        SeqIO.write(rec, f'{output_dir}/filtered_contigs.fa', 'fasta')

    # Prodigal translation
    prodigal_out_protein = f"{output_dir}/prot.fa"
    prodigal_cmd = f'prodigal -i {output_dir}/filtered_contigs.fa ' \
                   f'-a {prodigal_out_protein} -p meta'

    if not os.path.exists(prodigal_out_protein):
        print("\tProdigal translation...")
        _ = subprocess.check_call(prodigal_cmd, shell=True,
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("\tProtein file generated.")
    else:
        print(f"\tProtein file already exists, skipping prodigal step.")

    database = f"{ref_path}/prot_db"

    # running alignment
    top_c = top_cutoff + 1
    diamond_out_fp = f"{output_dir}/diamond_top{top_c}.tab"
    diamond_cmd = f'diamond blastp -q {output_dir}/prot.fa ' \
                  f'-d {database}.dmnd --threads {threads} ' \
                  f'-o {diamond_out_fp} --top {top_c} --sensitive'  # 10 + 1 = 11
    print("\tRunning Diamond...")
    if not os.path.exists(diamond_out_fp):
        # print("\n\n" + "{:-^80}".format("Diamond BLASTp"))
        print("\tAligning proteins against reference database...")
        _ = subprocess.check_call(diamond_cmd, shell=True,
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print("\tDiamond alignment finished.")
    else:
        print(f"\tDiamond output file already exists, skipping diamond step.")

    return output_dir, filtered_seq_list, prodigal_out_protein, diamond_out_fp, dict_seq_len


def _classify(input_data, ref_path, output_path, gini_cutoff, top_cutoff, length, threads):

    # ref files
    df_taxon2lineage = pd.read_csv(f"{ref_path}/taxon_lineage.csv", index_col=0, header=0)
    ref_label = f"{ref_path}/ref.csv"


    print("1 - Preprocessing ... (length filter / translation / DB construction / protein alignment)")
    output_dir, filtered_seq_list, prot_file, diamond_out, dict_seq_len =  \
        preprocessing_data(input_data, ref_path, output_path, top_cutoff, length, threads)

    print("2 - Parsing proteins... (read and parse the protein alignment result, generate protein taxonomic info)")
    #   1、alignment information of best hit
    #   2、how many genus got hit in the significant range
    #   3、top 1～k，name + score + ratio
    df_prot_stat, dict_prot_genus_score = handle_prot(ref_label, prot_file, diamond_out, top_cutoff=top_cutoff)
    df_prot_stat.to_csv(f"{output_dir}/vtm_prot.csv", header=True, index=True)

    # summarize protein info to contig level
    print("3 - Parsing contigs... (read protein taxonomic info and generate the contig taxonmic info)")
    df_seq_stat = aggregate_prot_to_seq(filtered_seq_list, dict_seq_len, df_prot_stat, dict_prot_genus_score)


    # make decision : filering + LCA
    print("4 - Making decision... (classify contigs)")
    df_seq_stat = pred(df_seq_stat, df_taxon2lineage=df_taxon2lineage, high=0.9, low=0.1, gini_cutoff=gini_cutoff)


    df_seq_stat.to_csv(f"{output_dir}/vtm_result.csv", header=True, index=True)

    return output_dir


def main(args):
    # check_cmd(args)


    input, output_path, ref_path = args["input"], args["output"], args["ref"]
    if ref_path == None:
        package_path = str(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        ref_path = f"{package_path}/ref"

    gini_cutoff = args["gini"]
    threads = args["threads"]
    length_cutoff = args["len"]
    top_cutoff = args["top"]

    output_dir = _classify(input_data=input, #input
                  ref_path=ref_path,   #ref
                  output_path=output_path, #out
                  length=length_cutoff,
                  top_cutoff=top_cutoff,
                  threads=threads,
                  gini_cutoff=gini_cutoff)

if __name__ == "__main__":
    package_path = str(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


    # top cutoff > x > 1-x/100
    output_dir = _classify(input_data=f"{package_path}/test/seq_test.fasta",
                            ref_path = f"{package_path}/ref",
                            output_path=f"{package_path}/test",
                            length=200,
                            top_cutoff=10,
                            threads=8,
                            gini_cutoff=0.44)
