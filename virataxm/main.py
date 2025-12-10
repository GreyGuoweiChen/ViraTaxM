import argparse
import sys
from .modules import update, classify, cluster
#import virataxm

def main():
    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f"""ViraTaxM: a Viral Taxonomic classification pipeline for Metagenomic sequences
    https://github.com/GreyGuoweiChen/VirGenus.git

    usage: virataxm <program> [options]

    programs:
        update              Rebuild reference database using updated data
        classify            Classify the query viruses till genus level
        cluster             For novel sequences at the genus level, provide a potential genus assignment
        predict             Run both classify and cluster modules sequentially  """,
    )

    subparsers = parser.add_subparsers(help=argparse.SUPPRESS)

    # update
    update_parser = subparsers.add_parser(
        "update",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Rebuild reference database using updated data
    \nusage: virataxm update auto
    \n   or: virataxm update -i <input.fasta> -t <lineage.csv> -d <database_path>""",
    )
    update.fetch_arguments(update_parser)
    # 提供两种模式，一种是 virataxm update auto，不用接后续参数，直接更新到默认路径；
    # 另一种是提供参数，指定输入输出路径。

    # classify
    classify_parser = subparsers.add_parser(
        "classify",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Classify the query viruses till genus level
    \nusage: virataxm classify -i <input.fasta> -o <output_path> -d <custom_database_path> [options] -""",
    )
    classify.fetch_arguments(classify_parser)


    # cluster
    cluster_parser = subparsers.add_parser(
        "cluster",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""For novel sequences at the genus level, provide a potential genus assignment
    \nusage: virataxm cluster -o <output_path> [options]""",
    )
    cluster.fetch_arguments(cluster_parser)


    # predict
    predict_parser = subparsers.add_parser(
        "predict",
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Run both classify and cluster modules sequentially
        \nusage: virataxm predict -i <input.fasta> -o <output_path> -d <custom_database_path> [options]""",
    )

    # 直接复用 classify 的参数定义
    classify.fetch_arguments(predict_parser)

    def run_predict(args):
        print("--------------------classification--------------------")
        classify.main(args)

        # print(args)
        print("--------------------clustering--------------------")
        cluster_args = {
            "output_dir": args["output"] + "/virataxm",
            "threads": args["threads"]
            # 如果 cluster 还有其他参数，可以在这里补充
        }
        cluster.main(cluster_args)

    predict_parser.set_defaults(func=run_predict)
    predict_parser.set_defaults(program="predict")



    # sys.argv[0] >> the first parameter of a command line
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "update":
            update_parser.print_help()
            sys.exit(0)
        if sys.argv[1] == "classify":
            classify_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "cluster":
            cluster_parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "predict":
            predict_parser.print_help()
            sys.exit(0)
        else:
            parser.print_help()
            sys.exit(0)

    print(f"You are runnning virataxm's '{sys.argv[1]}' program.")
    print(f"\t with command : {' '.join(sys.argv)}")
    args = vars(parser.parse_args())
    args["func"](args)

#main()
