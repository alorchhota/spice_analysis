def get_chip_files_for_file_aggregate_evals(wildcards, methods):
  tis = wildcards.get('tissue')
  chipseq_files = []
  if tis in chipseq_tissue_map and chipseq_tissue_map[tis] != "":
    chipseq_files = expand([# chipseq interaction auc
                            "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_auc.rds",
                            # chipseq interaction hub auc
                            "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_hub_auc.rds",
                            # chipseq interaction spearman cor
                            "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_spearman_cor.rds",
                            # chipseq interaction precision
                            "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_precision.rds"],
                            results_dir = wildcards.results_dir,
                            chipseq_tissue = tis,
                            correction_label = wildcards.correction_label,
                            gene_selection = wildcards.gene_selection,
                            n_genes = wildcards.n_genes,
                            method = methods)
  return chipseq_files

def get_chip_files_for_rule_aggregate_validation_eval_per_dataset(wildcards):
  return get_chip_files_for_file_aggregate_evals(wildcards, config['validation_methods'])

def get_chip_files_for_rule_aggregate_test_eval_per_dataset(wildcards):
  return get_chip_files_for_file_aggregate_evals(wildcards, config['test_methods'])

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
      # string_exp ppi auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_auc.rds",
      # string_exp ppi hub auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_hub_auc.rds",
      # string_exp ppi spearman cor
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_spearman_cor.rds",
      # string_exp ppi precision
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_precision.rds",
      # kegg interaction auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_auc.rds",
      # kegg interaction hub auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_hub_auc.rds",
      # kegg interaction spearman cor
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_spearman_cor.rds",
      # kegg interaction precision
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_precision.rds",
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
        pathway = config['pathways']),
    chipseq_files = get_chip_files_for_rule_aggregate_validation_eval_per_dataset
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
  
rule aggregate_test_eval_per_dataset:
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
      # string_kegg ppi auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_kegg_ppi_auc.rds",
      # string_kegg ppi hub auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_kegg_ppi_hub_auc.rds",
      # string_kegg ppi spearman cor
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_kegg_ppi_spearman_cor.rds",
      # string_kegg ppi precision
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_kegg_ppi_precision.rds",
      # string_exp ppi auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_auc.rds",
      # string_exp ppi hub auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_hub_auc.rds",
      # string_exp ppi spearman cor
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_spearman_cor.rds",
      # string_exp ppi precision
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_precision.rds",
      # kegg interaction auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_auc.rds",
      # kegg interaction hub auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_hub_auc.rds",
      # kegg interaction spearman cor
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_spearman_cor.rds",
      # kegg interaction precision
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_kegg_interaction_precision.rds",
      # inweb ppi auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_inweb_ppi_auc.rds",
      # inweb ppi hub auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_inweb_ppi_hub_auc.rds",
      # inweb ppi spearman cor
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_inweb_ppi_spearman_cor.rds",
      # inweb ppi precision
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_inweb_ppi_precision.rds",
      # shared pathway auc
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_shared_pathway_auc_{pathway}.rds",
      #pathway enrichment
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_pathway_enrichment_{pathway}.rds"],
        results_dir = "{results_dir}", 
        tissue = "{tissue}", 
        correction_label = "{correction_label}", 
        gene_selection = "{gene_selection}", 
        n_genes = "{n_genes}", 
        method = config['test_methods'],
        pathway = config['pathways']),
    chipseq_files = get_chip_files_for_rule_aggregate_test_eval_per_dataset
  params:
    dir = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}",
    methods = ",".join([str(s) for s in config['test_methods']])
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/aggregated_evaluations_test.rds"
  log: 
    "{results_dir}/gtex_v8/logs/aggregate_test_eval_per_dataset/{tissue}/{correction_label}/{gene_selection}/{n_genes}/aggregate_test_eval_per_dataset.log"
  shell:
    """
    Rscript src/aggregate_evaluations_per_dataset.R \
      --dir "{params.dir}" \
      --methods "{params.methods}" \
      --o "{output}" \
      2>&1 | tee {log}
    """
  

rule aggregate_all_test_evals:
  input:
    expand(
      ["{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/aggregated_evaluations_test.rds"],
        results_dir = "{results_dir}",
        tissue = config['test_tissues'],
        correction_label = config['test_correction_labels'],
        gene_selection = config['test_gene_selection'],
        n_genes = config['test_n_genes'])
  params:
    dir = "{results_dir}/gtex_v8/results",
    tissues = ",".join([str(s) for s in config['test_tissues']]),
    correction_labels = ",".join([str(s) for s in config['test_correction_labels']]),
    gene_selections = ",".join([str(s) for s in config['test_gene_selection']]),
    n_genes = ",".join([str(s) for s in config['test_n_genes']]),
    fn = "aggregated_evaluations_test.rds"
  output:
    "{results_dir}/gtex_v8/aggregated/all_evaluations_test.rds"
  log: 
    "{results_dir}/gtex_v8/logs/aggregate_all_test_evals/aggregate_all_test_evals.log"
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
