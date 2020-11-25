shell.prefix("module load R/4.0.2; ")  # load modules in marcc. may produce error message in other systems.

configfile: "config/config.yaml"
include: "config/derived_config.smk"

rule all:
  input:
    expand(
      [
        # aggregate evaluations for validation
        "{results_dir}/gtex_v8/aggregated/all_evaluations_validation.rds",
        # aggregate evaluations for test
        "{results_dir}/gtex_v8/aggregated/all_evaluations_test.rds"
      ],
      results_dir = config['results_dir'])

rule all_validation:
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
      # chipseq interaction auc
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_auc.rds",
      # chipseq interaction hub auc
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_hub_auc.rds",
      # chipseq interaction spearman cor
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_spearman_cor.rds",
      # chipseq interaction precision
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_precision.rds",
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
        chipseq_tissue = validation_chipseq_tissues,
        correction_label = config['validation_correction_labels'], 
        gene_selection=config['validation_gene_selection'], 
        n_genes=config['validation_n_genes'], 
        method = config['validation_methods'],
        pathway = config['pathways'])

rule all_test:
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
      # string_kegg ppi
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_kegg_ppi.rds",
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
      # chipseq interaction auc
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_auc.rds",
      # chipseq interaction hub auc
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_hub_auc.rds",
      # chipseq interaction spearman cor
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_spearman_cor.rds",
      # chipseq interaction precision
      "{results_dir}/gtex_v8/results/{chipseq_tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_precision.rds",
      # inweb ppi
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/inweb_ppi.rds",
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
      "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_pathway_enrichment_{pathway}.rds",
      # aggregate evaluations for test
      "{results_dir}/gtex_v8/aggregated/all_evaluations_test.rds" 
      ],
      #
        results_dir = config['results_dir'],
        tissue = config['test_tissues'],
        chipseq_tissue = test_chipseq_tissues,
        correction_label = config['test_correction_labels'],
        gene_selection=config['test_gene_selection'],
        n_genes=config['test_n_genes'],
        method = config['test_methods'],
        pathway = config['pathways'])


include: "rules/resource_allocation.smk"
include: "rules/download_gtex.smk"
include: "rules/download_resources.smk"
include: "rules/data_correction.smk"
include: "rules/run_networks.smk"
include: "rules/string_ppi.smk"
include: "rules/inweb_ppi.smk"
include: "rules/msigdb_genesets.smk"
include: "rules/eval_interactions.smk"
include: "rules/eval_string_kegg_interactions.smk"
include: "rules/eval_string_exp_interactions.smk"
include: "rules/eval_inweb_interactions.smk"
include: "rules/eval_kegg_interactions.smk"
include: "rules/eval_chipseq_interactions.smk"
include: "rules/eval_pathways.smk"
include: "rules/eval_trans_eqtl.smk"
include: "rules/aggregate_evals.smk"

