#!/usr/bin/python2.7

import logging
import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser(prog = "joinRes",
                                     description = "Join eQTL prediction results.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-a', '--matrix-anova', help = "MatrixEQTL Anova result file.",
                        required = True, dest = 'anova_file')
    parser.add_argument('-l', '--matrix-linear', help = "MatrixEQTL Linear result file.",
                        required = True, dest = 'linear_file')
    parser.add_argument('-m', '--mRMR', help = "mRMR result file.",
                        required = True, dest = 'mrmr_file')
    parser.add_argument('-e', '--eqtl-em', help = "R/QTL EM result file.",
                        required = True, dest = 'em_file')
    parser.add_argument('-k', '--eqtl-hk', help = "R/QTL HK result file.",
                        required = True, dest = 'hk_file')
    parser.add_argument('-t', '--truth', help = "True predictions.",
                        required = True, dest = 'truth_file')
    parser.add_argument('-o', '--out-file', help = "CSV output file.",
                        required = True, dest = 'out_file')
    parser.add_argument('-v', '--verbose',
                        help='increase output verbosity',
                        action='count', default=0)
    args = parser.parse_args()

    if args.verbose == 0:
        log_level = logging.INFO
    elif args.verbose == 1:                                                                                             
        log_level = logging.DEBUG
    else:                                                                                                               
        log_level = logging.DEBUG                                                                                       
                                                                                                                        
    logging.basicConfig(level=log_level,
                        format='%(levelname)-8s [%(asctime)s]  %(message)s',                                            
                        datefmt="%y%m%d %H%M%S")                                                                        
                                                                                                                        
    logging.info("Program Started")
	# Check input files
    if not os.path.isfile(args.em_file):
        logging.error("File " + args.em_file + " not found!")
        exit(1)
    if not os.path.isfile(args.hk_file):
        logging.error("File " + args.hk_file + " not found!")
        exit(1)
    if not os.path.isfile(args.mrmr_file):
        logging.error("File " + args.mrmr_file + " not found!")
        exit(1)
    if not os.path.isfile(args.anova_file):
        logging.error("File " + args.anova_file + " not found!")
        exit(1)
    if not os.path.isfile(args.linear_file):
        logging.error("File " + args.linear_file + " not found!")
        exit(1)

    # R/QTL - EM predictions
    em_dict = {}
    em_start_col = 3
    with open(args.em_file, "r") as em_data:
        em_genes = em_data.readline().rstrip().split(",")[em_start_col:]
        em_snps = set()
        logging.info("EM genes: '%d'", len(em_genes))
        for em_pred in em_data.readlines():
            preds = em_pred.rstrip().split(",")
            em_snps.add(preds[0])
            for p in range(em_start_col, len(preds)):
                em_dict[(preds[0], em_genes[p - em_start_col])] = preds[p]
        logging.info("EM snps: '%d'", len(em_snps))
    logging.info("EM predictions: '%d'", len(em_dict))
    # R/QTL - HK predictions
    hk_dict = {}
    hk_start_col = 3
    with open(args.hk_file, "r") as hk_data:
        hk_genes = hk_data.readline().rstrip().split(",")[hk_start_col:]
        hk_snps = set()
        logging.info("HK genes: '%d'", len(hk_genes))
        for hk_pred in hk_data.readlines():
            preds = hk_pred.rstrip().split(",")
            hk_snps.add(preds[0])
            for p in range(hk_start_col, len(preds)):
                hk_dict[(preds[0], hk_genes[p - hk_start_col])] = preds[p]
        logging.info("HK snps: '%d'", len(hk_snps))
    logging.info("HK predictions: '%d'", len(hk_dict))
    # mRMR predictions
    mrmr_dict = {}
    with open(args.mrmr_file, "r") as mrmr_data:
        h_mrmr = mrmr_data.readline()
        mrmr_genes = set()
        mrmr_snps = set()
        for mrmr_pred in mrmr_data.readlines():
            preds = mrmr_pred.rstrip().split("\t")
            mrmr_genes.add(preds[2])
            mrmr_snps.add(preds[0])
            mrmr_dict[(preds[0], preds[2])] = preds[3:]
        logging.info("mRMR genes: '%d'", len(mrmr_genes))
        logging.info("mRMR snps: '%d'", len(mrmr_snps))
    logging.info("mRMR predictions: '%d'", len(mrmr_dict))
    # Matrix eQTL - Anova predictions
    anova_dict = {}
    with open(args.anova_file, "r") as anova_data:
        h_anova = anova_data.readline()
        anova_genes = set()
        anova_snps = set()
        for anova_pred in anova_data.readlines():
            preds = anova_pred.rstrip().split("\t")
            anova_genes.add(preds[1])
            anova_snps.add(preds[0])
            anova_dict[(preds[0], preds[1])] = preds[2:]
        logging.info("Matrix eQTL Anova genes: '%d'", len(anova_genes))
        logging.info("Matrix eQTL Anova snps: '%d'", len(anova_snps))
    logging.info("Matrix eQTL Anova predictions: '%d'", len(anova_dict))
    # Matrix EQTL - Linear predictions
    linear_dict = {}
    with open(args.linear_file, "r") as linear_data:
        h_linear = linear_data.readline()
        linear_genes = set()
        linear_snps = set()
        for linear_pred in linear_data.readlines():
            preds = linear_pred.rstrip().split("\t")
            linear_genes.add(preds[1])
            linear_snps.add(preds[0])
            linear_dict[(preds[0], preds[1])] = preds[2:]
        logging.info("Matrix eQTL Linear genes: '%d'", len(linear_genes))
        logging.info("Matrix eQTL Linear snps: '%d'", len(linear_snps))
    logging.info("Matrix eQTL Linear predictions: '%d'", len(linear_dict))
    # True predictions
    true_dict = {}
    with open(args.truth_file, "r") as truth:
        for t in truth.readlines():
            t_vals = t.rstrip().split("\t")
            t_s = "snp_" + t_vals[0]
            t_g = "gene_" + t_vals[1]
            logging.debug("snp: %s -- gene: %s", t_s, t_g)
            true_dict[(t_s, t_g)] = 1
    logging.info("True predictions: '%d'", len(true_dict))
    
    f_h = ["Snp", "Gene", "RQTL-EM", "RQTL-HK", "mRMR-Score", "mRMR-Accuracy",
           "MatrixEQTL-Anova-F-test", "MatrixEQTL-Anova-p-value", "MatrixEQTL-Anova-FDR",
           "MatrixEQTL-Linear-beta", "MatrixEQTL-Linear-t-test", "MatrixEQTL-Linear-p-value",
           "MatrixEQTL-Linear-FDR", "Truth"]
    with open(args.out_file, "w") as out:
        out.write(",".join(f_h) + "\n")
        for i in em_dict:
            out.write(",".join(i) + "," + em_dict[i] + "," + hk_dict[i])
            mRMR_out = ["NA", "NA"]
            mRMR_out = ["0", "NA"]
            if i in mrmr_dict:
                mRMR_out = mrmr_dict[i]
            out.write("," + ",".join(mRMR_out))
            anova_out = ["NA", "NA", "NA"]
            anova_out = ["0", "1", "1"]
            if i in anova_dict:
                anova_out = anova_dict[i]
            out.write("," + ",".join(anova_out))
            linear_out = ["NA", "NA", "NA", "NA"]
            linear_out = ["0", "0", "1", "1"]
            if i in linear_dict:
                linear_out = linear_dict[i]
            out.write("," + ",".join(linear_out))
            tr_val = "0"
            if i in true_dict:
                tr_val = "1"
            out.write("," + tr_val)
            out.write("\n")
    logging.info("Program Finished")


if __name__ == "__main__":                                                                                              
    main()  
