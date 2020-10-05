import socket
if "cluster" in socket.gethostname():
  shell.prefix("module load R; ")

configfile: "config/config.yaml"
results_dir = config["results_dir"]

rule all:
  input:
    #expand(["{results_dir}/gtex_v8/download/GTEx_Analysis_v8_eQTL_covariates.tar.gz",
    #        "{results_dir}/gtex_v8/download/GTEx_Analysis_v8_eQTL_expression_matrices.tar",
    #        "{results_dir}/gtex_v8/download/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"],
    #        results_dir = results_dir)
    #
    #expand("{results_dir}/gtex_v8/extracted/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt", results_dir = config['results_dir'], tissue = config['tissues'])
    #
    expand("{results_dir}/gtex_v8/data/normalized/{tissue}.normalized.txt", results_dir = config['results_dir'], tissue = config['tissues']),
    expand("{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt", results_dir = config['results_dir'], tissue = config['tissues']),
    expand("{results_dir}/gtex_v8/data/raw_tpm/{tissue}.txt", results_dir = config['results_dir'], tissue = config['tissues'])
include: "rules/download_gtex.smk"
