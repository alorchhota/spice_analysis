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
  

rule aggregate_all_validation_evals:
  input:
    expand(
      ["{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/aggregated_evaluations_validation.rds"],
        results_dir = "{results_dir}",
        tissue = config['validation_tissues'],
        correction_label = config['validation_correction_labels'],
        gene_selection = config['validation_gene_selection'],
        n_genes = config['validation_n_genes'])
  params:
    dir = "{results_dir}/gtex_v8/results",
    tissues = ",".join([str(s) for s in config['validation_tissues']]),
    correction_labels = ",".join([str(s) for s in config['validation_correction_labels']]),
    gene_selections = ",".join([str(s) for s in config['validation_gene_selection']]),
    n_genes = ",".join([str(s) for s in config['validation_n_genes']]),
    fn = "aggregated_evaluations_validation.rds"
  output:
    "{results_dir}/gtex_v8/aggregated/all_evaluations_validation.rds"
  log: 
    "{results_dir}/gtex_v8/logs/aggregate_all_validation_evals/aggregate_all_validation_evals.log"
  shell:
    """
    Rscript src/aggregate_aggregated_evaluations.R \
      --dir "{params.dir}" \
      --tissue "{params.tissues}" \
      --correction "{params.correction_labels}" \
      --gene_selection "{params.gene_selections}" \
      --n_genes "{params.n_genes}" \
      --fn "{params.fn}" \
      --o "{output}" \
      2>&1 | tee {log}
    """
  
