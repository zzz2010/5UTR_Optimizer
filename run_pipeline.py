import sys,os
import argparse
def get_cmd_parser():
    main_parser = argparse.ArgumentParser()

    main_parser.add_argument('--help',
                        help='run_pipeline.py [task]  [task parameters]')

    main_parser.add_argument('--output_dir',
                        default='./output',
                        help='the output directory')

    subparsers = main_parser.add_subparsers(help='feature_extraction | model_build | model_eval | sequence_generation', dest='Pass task name and task args')
    feature_extraction_parser = subparsers.add_parser(name="feature_extraction",
                                                    help="feature extraction step, convert fasta file to sequence feature format")

    feature_extraction_parser.add_argument("--input_fasta",help="5utr sequences in fasta format")


    model_build_parser = subparsers.add_parser(name="model_build",
                                                    help="build machine learning model from given sequence features and label data ")

    model_build_parser.add_argument("--input_prefix",help="output prefix from feature_extraction step")


    return main_parser

if __name__ == '__main__':
    parser=get_cmd_parser()
    args = parser.parse_args()
    print(args)