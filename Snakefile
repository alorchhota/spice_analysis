shell.prefix("module load R/4.0.2; ")  # load modules in marcc. may produce error message in other systems.

configfile: "config/config.yaml"
results_dir = config["results_dir"]

rule all:
  input:
    expand(
      [ # network run
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_network.rds",
      # string ppi
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_ppi.rds",
      # string ppi auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_ppi_auc.rds",
      # string ppi hub auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_ppi_hub_auc.rds",
      # string ppi spearman cor
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_ppi_spearman_cor.rds",
      # string ppi precision
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_ppi_precision.rds",
      # shared pathway auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_shared_pathway_auc_{pathway}.rds",
      #pathway enrichment
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_pathway_enrichment_{pathway}.rds",
      # aggregate evaluations for validation
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/aggregated_evaluations_validation.rds",
      "{results_dir}/gtex_v8/aggregated/all_evaluations_validation.rds"],
      #
        results_dir = config['results_dir'], 
        tissue = config['validation_tissues'], 
        correction_label = config['validation_correction_labels'], 
        gene_selection=config['validation_gene_selection'], 
        n_genes=config['validation_n_genes'], 
        method = config['validation_methods'],
        pathway = config['pathways'])

include: "rules/download_gtex.smk"
include: "rules/download_resources.smk"
include: "rules/data_correction.smk"
include: "rules/run_networks.smk"
include: "rules/string_ppi.smk"
include: "rules/msigdb_genesets.smk"
include: "rules/eval_interactions.smk"
include: "rules/eval_pathways.smk"
include: "rules/aggregate_evals.smk"
