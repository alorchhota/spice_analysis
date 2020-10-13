rule aggregate_validation_eval_per_dataset:
  input:
    expand(
      [# string ppi auc
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
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_pathway_enrichment_{pathway}.rds"],
        results_dir = "{results_dir}", 
        tissue = "{tissue}", 
        correction_label = "{correction_label}", 
        gene_selection = "{gene_selection}", 
        n_genes = "{n_genes}", 
        method = config['validation_methods'],
        pathway = config['pathways'])
  params:
    dir = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}",
    methods = ",".join([str(s) for s in config['validation_methods']])
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/aggregated_evaluations_validation.rds"
  log: 
    "{results_dir}/gtex_v8/logs/aggregate_validation_eval_per_dataset/{tissue}/{correction_label}/{gene_selection}/{n_genes}/aggregate_validation_eval_per_dataset.log"
  shell:
    """
    Rscript src/aggregate_evaluations_per_dataset.R \
      --dir "{params.dir}" \
      --methods "{params.methods}" \
      --o "{output}" \
      2>&1 | tee {log}
    """
  
