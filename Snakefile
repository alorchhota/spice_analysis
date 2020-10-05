import socket
if "cluster" in socket.gethostname():
  shell.prefix("module load R; ")

configfile: "config/config.yaml"
results_dir = config["results_dir"]

rule all:
  input:
    expand("{results_dir}/gtex_v8/data/corrected/{tissue}.corrected.{gene_selection}.{n_genes}.txt", results_dir = config['results_dir'], tissue = config['tissues'], gene_selection=config['gene_selection'], n_genes=config['n_genes'])

include: "rules/download_gtex.smk"
include: "rules/data_correction.smk"
