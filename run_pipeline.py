import sys,os
import argparse
def get_cmd_parser():
    main_parser = argparse.ArgumentParser(add_help=False)

    main_parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s 1.0', help="Show program's version number and exit.")
    main_parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                        help='run_pipeline.py [task]  [task parameters]',)



    subparsers = main_parser.add_subparsers(help='feature_extract | model_build | model_eval | sequence_generate | full_pipe', dest='taskname')
    feature_extraction_parser = subparsers.add_parser(name="feature_extract",
                                                    help="""feature extraction step, convert fasta file to sequence feature format.
                                                     Example: feature_extract --input_fasta data/gencode_v17_5utr_15bpcds.fa --output_dir output/""")

    feature_extraction_parser.add_argument("--input_fasta",required=True,help="5utr sequences in fasta format")

    feature_extraction_parser.add_argument('--output_dir',
                        default='./output',
                        help='the output directory')


    model_build_parser = subparsers.add_parser(name="model_build",
                                                    help="""build machine learning model from given sequence features and label data 
                                                    Example: model_build --prefix output/input.fa --annotation_file data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt --min_rna_rpkm 5 --min_riboseq_rpkm 0.1 --model 1 --out output/muscle_randomforest.model
                                                    """)

    model_build_parser.add_argument("-p", "--prefix",required=True,help="sequence feature files prefix")
    model_build_parser.add_argument("-k", "--min_rna_rpkm", help="minimal RNASEQ RPKM for the gene included in training. In paper, we used 5 for muscle cell")
    model_build_parser.add_argument("-r", "--min_riboseq_rpkm", help="minimal RiboSEQ RPKM for the gene included in training. In paper, we used 0.1 for muscle cell")
    model_build_parser.add_argument("-a", "--annotation_file", required=True,help="transcript annotation file with gene id and RPKM information for rnaseq and riboseq")
    model_build_parser.add_argument("-m", "--model" , default=1, help="1: randomforest, 2: glmnet, 3: regression tree, 4: SVM  [default= %default]")
    model_build_parser.add_argument("-j", "--n_cpus", default=16, help="number of cpus [default= %default]")
    model_build_parser.add_argument("-o", "--out", default="trained.model", help="output file name [default= %default]")




    model_eval_parser = subparsers.add_parser(name="model_eval",
                                                    help="eval the built ML model using cross validation")


    model_eval_parser.add_argument("-p", "--prefix",required=True,help="sequence feature files prefix")
    model_eval_parser.add_argument("-k", "--min_rna_rpkm", help="minimal RNASEQ RPKM for the gene included in training. In paper, we used 5 for muscle cell")
    model_eval_parser.add_argument("-r", "--min_riboseq_rpkm", help="minimal RiboSEQ RPKM for the gene included in training. In paper, we used 0.1 for muscle cell")
    model_eval_parser.add_argument("-a", "--annotation_file", required=True,help="transcript annotation file with gene id and RPKM information for rnaseq and riboseq")
    model_eval_parser.add_argument("-m", "--modellist" , default="1,2,3,4", help="1: randomforest, 2: glmnet, 3: regression tree, 4: SVM")
    model_eval_parser.add_argument("-j", "--n_cpus", default=16, help="number of cpus [default= %default]")
    model_eval_parser.add_argument("-o", "--out", default="evaluation.pdf", help="output evaluation pdf file name [default= %default]")

    sequence_generate_parser = subparsers.add_parser(name="sequence_generate",
                                                    help="""generate synthetic UTR sequences with optimized Riob level or transcription efficiency based on the pretraned ML model\
                                                    Example: sequence_generate --prefix output/input.fa --annotation_file data/df_counts_and_len.TE_sorted.Muscle.with_annot.txt --min_rna_rpkm 5 --min_riboseq_rpkm 0.1 -t ribo  -m output/muscle_randomforest.model  -o output/
                                                    """)

    sequence_generate_parser.add_argument("-p", "--prefix",required=True,help="sequence feature files prefix")
    sequence_generate_parser.add_argument("-k", "--min_rna_rpkm", help="minimal RNASEQ RPKM for the gene included in training. In paper, we used 5 for muscle cell")
    sequence_generate_parser.add_argument("-r", "--min_riboseq_rpkm", help="minimal RiboSEQ RPKM for the gene included in training. In paper, we used 0.1 for muscle cell")
    sequence_generate_parser.add_argument("-a", "--annotation_file", required=True,help="transcript annotation file with gene id and RPKM information for rnaseq and riboseq")
    sequence_generate_parser.add_argument("-m", "--model_file" ,required=True,   help="path to the saved model file")
    sequence_generate_parser.add_argument("-j", "--n_cpus", default=16, help="number of cpus [default= %default]")
    sequence_generate_parser.add_argument("-o", "--out", default="./output", help="output folder path [default= %default]")
    sequence_generate_parser.add_argument("-t", "--optimize_target", default="ribo",
                                          help="optimize_target: ribo | te  [default= %default]")
    sequence_generate_parser.add_argument("-n", "--n_total", default=3585,
                                          help="total number of synthetic sequences")

    ##### full pipeline section ######
    fullpipe_parser = subparsers.add_parser(name="full_pipe",
                                                    help="run full pipeline from start to the end: feature_extract -> model_build -> model_eval -> sequence_generate ")

    fullpipe_parser.add_argument("-i","--input_fasta", required=True, help="5utr sequences in fasta format")

    fullpipe_parser.add_argument("-o",'--output_dir',
                        default='./output',
                        help='the output directory')

    fullpipe_parser.add_argument("-k", "--min_rna_rpkm", help="minimal RNASEQ RPKM for the gene included in training. In paper, we used 5 for muscle cell")
    fullpipe_parser.add_argument("-r", "--min_riboseq_rpkm", help="minimal RiboSEQ RPKM for the gene included in training. In paper, we used 0.1 for muscle cell")
    fullpipe_parser.add_argument("-a", "--annotation_file", help="transcript annotation file with gene id and RPKM information for rnaseq and riboseq")
    fullpipe_parser.add_argument("-m", "--model" , default=1, help="1: randomforest, 2: glmnet, 3: regression tree, 4: SVM  [default= %default]")
    fullpipe_parser.add_argument("-j", "--n_cpus", default=16, help="number of cpus [default= %default]")
    fullpipe_parser.add_argument("-n", "--n_total", default=3585,
                                          help="total number of synthetic sequences")



    return main_parser

def oss(cmd):
    print("executing:", cmd)
    os.system(cmd)

if __name__ == '__main__':
    parser=get_cmd_parser()
    args = parser.parse_args()
    if args.taskname=="feature_extract" or args.taskname=="full_pipe":
        os.makedirs(args.output_dir,exist_ok=True)
        cmd="python FeatureExtraction_final.py %s %s"%(args.input_fasta,args.output_dir)
        oss(cmd)

        if args.taskname=="full_pipe":
            args.prefix=args.output_dir+"/input.fa"



    if args.taskname=="model_build" or args.taskname=="full_pipe":
        if args.taskname=="full_pipe":
            args.out=args.output_dir+"/trained.model"
        cmd="Rscript buildModel_final.R -p %s -a %s -k %s -r %s -m %s -o %s"%(args.prefix,
                                                                              args.annotation_file,
                                                                            args.min_rna_rpkm,
                                                                             args.min_riboseq_rpkm,
                                                                              args.model,args.out  )
        oss(cmd)


    if args.taskname=="sequence_generate" or args.taskname=="full_pipe":
        for seed in range(int(args.n_total)):
            if args.taskname=="full_pipe":
                args.model_file=args.output_dir+"/trained.model"
                for optimize_target in ['ribo','te']:
                    cmd = "Rscript evolutionDesign.R -p %s -a %s -t %s -k %s -r %s -m %s  -o %s" % (args.prefix,
                                                                                                args.annotation_file,
                                                                                                optimize_target,
                                                                                                args.min_rna_rpkm,
                                                                                                args.min_riboseq_rpkm,
                                                                                                args.model_file, args.out)
                    oss(cmd)

            else:
                cmd="Rscript evolutionDesign.R -p %s -a %s -t %s -k %s -r %s -m %s  -o %s"%(args.prefix,
                                                                                      args.annotation_file, args.optimize_target,
                                                                                    args.min_rna_rpkm,
                                                                                     args.min_riboseq_rpkm,
                                                                                      args.model_file,args.out  )
                oss(cmd)
