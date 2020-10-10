shell.prefix("module load R/4.0.2; ")  # load modules in marcc. may produce error message in other systems.

configfile: "config/config.yaml"
results_dir = config["results_dir"]


rule all:
  input:
    expand("{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_network.rds", results_dir = config['results_dir'], tissue = config['validation_tissues'], correction_label = config['validation_correction_labels'], gene_selection=config['validation_gene_selection'], n_genes=config['validation_n_genes'], method = config['validation_methods']),
    expand("{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_ppi.rds", results_dir = config['results_dir'], tissue = config['validation_tissues'], correction_label = config['validation_correction_labels'], gene_selection=config['validation_gene_selection'], n_genes=config['validation_n_genes'])

include: "rules/download_gtex.smk"
include: "rules/download_resources.smk"
include: "rules/data_correction.smk"
include: "rules/run_networks.smk"
include: "rules/string_ppi.smk"

