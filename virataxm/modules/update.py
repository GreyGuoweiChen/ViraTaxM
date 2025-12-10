#coding=utf-8

import os
import subprocess
from Bio import SeqIO
import pandas as pd
import csv
import zipfile

# num_classes = 1326
PACKAGE_PATH = str(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def fetch_arguments(parser):
    parser.set_defaults(func=main)
    parser.set_defaults(program="update")

    parser.add_argument('--auto', action='store_true', help="Automatically build from default reference database, (no other arguments needed).")

    parser.add_argument('-i', '--input', type = str, help = "The viral reference sequences (.fasta)")
    parser.add_argument('-t', '--taxa', type = str, help = "The taxonomic lineage of reference sequences. (.csv)")
    parser.add_argument('-o', '--output', default = "./vtm_ref", type = str,
                        help="The reference path，including the taxonomic files and genomic files.")

    parser.add_argument('-f', action='store_true', help="Overwrite the existing database files, which takes around an hour to go.")



def check_cmd(args):
    if args["auto"]:
        if args["input"] is not None or args["taxa"] is not None:
            raise Exception("When using --auto, do not provide -i/--input or -t/--taxa arguments.")
        return
    else:
        file_seq = args["input"]
        if not os.path.exists(file_seq):
            raise Exception("The input fasta file do not exist.")

        file_taxa = args["taxa"]
        if not os.path.exists(file_taxa):
            raise Exception("The taxonomic information is required.")

        try:    # check seq number of fasta and file_taxa
            num = 0
            with open(file_seq) as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    num +=1
            df_taxa = pd.read_csv(file_taxa, index_col=0, header=0)
            if num != len(df_taxa):
                raise Exception(f"The number of input sequences do not equal to the number of taxonomic information\n"
                                f"Fasta file includes {num}\n"
                                f"Taxonomic file includes {len(df_taxa)} records")
        except Exception as e:
            print(e)
            exit()

    output_dir = args["output"]
    if os.path.exists(output_dir):
        raise Exception(f"The output directory ({output_dir}) already exists.")



def _update(fasta, label, output_dir, threads, should_overwrite):
    print(f"Updating the reference database...")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Prodigal translation
    prodigal_out_protein = f"{output_dir}/prot.faa"
    pdg_annotation = f"{output_dir}/prot.gbk"
    if should_overwrite or not os.path.exists(prodigal_out_protein):
        print(f"\tTranslating contigs to proteins...")
        prodigal_cmd = f'prodigal -i {fasta} -a {prodigal_out_protein} -o {pdg_annotation} -p meta'
        _ = subprocess.check_call(prodigal_cmd, shell=True,)
                                  # stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.path.exists(pdg_annotation):
            os.remove(pdg_annotation)
        print("\tTranslation finished.")
    else:
        print("\tWarning: The reference files already exist. Use -f to overwrite it.")
        return output_dir
    # else:
    #     print("\tThe protein reference already exist.")

    # build diamond database
    diamond_db = f"{output_dir}/prot_db"
    diamond_cmd = f'diamond makedb --in {prodigal_out_protein} -d {diamond_db}'
    print("\tBuilding Diamond database...")
    _ = subprocess.check_call(diamond_cmd, shell=True,
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print("\tDiamond database built.")

    ########## taxon lineage + count ##########
    df = pd.read_csv(label, index_col=0, header=0)
    rank_list = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    df_taxon_lineage = pd.DataFrame(columns=["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Count", "Rank"])
    for i, rank in enumerate(rank_list):
        df_rank = df[rank].value_counts(ascending=False).to_frame().reset_index()   # dropna
        df_rank.columns = [rank, "Count"]
        # print(df_rank)
        df_rank_taxon = pd.DataFrame(index=df_rank[rank].unique(),
                                columns=["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"])
        df_group = df.groupby([rank])
        for j in range(0, i+1):
            rank_j = rank_list[j]
            df_rank_taxon[rank_j] = df_group[rank_j].first()
        df_rank_taxon["Count"] = df_rank["Count"].tolist()
        df_rank_taxon["Rank"] = rank
        df_taxon_lineage = pd.concat([df_taxon_lineage, df_rank_taxon], axis=0)
    df_taxon_lineage.to_csv(f"{output_dir}/taxon_lineage.csv", index=True)
    print(f"\tTaxonomic lineage file generated as ref_path/taxon_lineage.csv.")

    # # record the protein taxonomic lineage
    out_file = f"{output_dir}/prot_lineage.csv"
    header = [None, "length", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]
    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for record in SeqIO.parse(prodigal_out_protein, "fasta"):
            prot_id = record.id
            prot_len = len(record.seq)
            seq_acc = prot_id.rsplit("_", 1)[0]

            lineage = df.loc[seq_acc, ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]].fillna("").tolist()
            writer.writerow([prot_id, prot_len] + lineage)

    print(f"\tProtein taxonomic lineage file generated as ref_path/prot_lineage.csv")

    return output_dir



def main(args):
    check_cmd(args)

    if args["auto"]:
        input_fasta = f"{PACKAGE_PATH}/ref/ref.fasta"
        input_taxa = f"{PACKAGE_PATH}/ref/ref.csv"
        output_dir = f"{PACKAGE_PATH}/ref"
        if args["auto"] and (args["input"] or args["taxa"]):
            raise IOError("--auto cannot be used together with --input/--taxa")

        output_dir_zip = f"{PACKAGE_PATH}/ref.zip"
        if os.path.exists(output_dir_zip) and not os.path.exists(output_dir):
            print("\tUnzip the built-in reference database ...")
            with zipfile.ZipFile(output_dir_zip, 'r') as zip_ref:
                zip_ref.extractall(PACKAGE_PATH)
            print("\tThe reference database is set")
            return
    else:
        input_fasta, input_taxa, output_dir = args["input"], args["taxa"], args["output"]
    should_overwrite = args["f"]
    os.makedirs(output_dir, exist_ok=True)

    _update(fasta=input_fasta,
            label=input_taxa,
            output_dir= output_dir,
            threads =32,
            should_overwrite=should_overwrite)

# if run directly
if __name__ == "__main__":
    input_fasta = f"{PACKAGE_PATH}/ref/ref.fasta"
    input_taxa = f"{PACKAGE_PATH}/ref/ref.csv"
    output_dir = f"{PACKAGE_PATH}/ref"
    _update(fasta=input_fasta,
            label=input_taxa,
            output_dir=output_dir,
            threads =32)