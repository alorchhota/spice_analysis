def get_chip_file_for_rule_get_chipseq_interactions(wildcards):
  tis = wildcards.get('tissue')
  if tis not in chipseq_tissue_map or chipseq_tissue_map[tis] == "":
    raise ValueError("rule get_chipseq_interactions: no chip-seq mapped to " + tis + "tissue.")
  chip_type = chipseq_tissue_map[tis]
  chip_fn = wildcards.results_dir + "/shared_data/chipseq/per_tissue/chipseq_" + chip_type + ".rds"
  return chip_fn

rule split_chipseq_interactions:
  input:
    chipseq_file = "{results_dir}/shared_data/chipseq/chipseq_interactions.bed"
  params:
    out_pfx = "{results_dir}/shared_data/chipseq/per_tissue/chipseq"
  output:
    expand("{results_dir}/shared_data/chipseq/per_tissue/chipseq_{chip_type}.rds",
      results_dir = "{results_dir}",
      chip_type = chipseq_tissue_types)
  log: 
    "{results_dir}/logs/split_chipseq_interactions/split_chipseq_interactions.log"
  shell:
    """
    Rscript src/split_chipseq_per_tissue_type.R \
      --chip "{input.chipseq_file}" \
      --o "{params.out_pfx}" \
      2>&1 | tee {log}
    """

rule get_chipseq_interactions:
  input:
    gene_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/genes.txt",
    chipseq_file = get_chip_file_for_rule_get_chipseq_interactions
  params:
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    interaction_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction.rds"
  log: 
    "{results_dir}/gtex_v8/logs/get_chipseq_interactions/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.log"
  shell:
    """
    Rscript src/prepare_chipseq_interaction_matrix.R \
      --gene "{input.gene_file}" \
      --chip "{input.chipseq_file}" \
      --o "{output.interaction_file}" \
      2>&1 | tee {log}
    """

rule eval_chipseq_interaction_auc:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction.rds"
  params:
    curve = config['eval_opts']['auc']['curve'],
    max = config['eval_opts']['auc']['max'],
    min = config['eval_opts']['auc']['min'],
    rand = config['eval_opts']['auc']['rand'],
    dg = config['eval_opts']['auc']['dg'],
    na = config['eval_opts']['weight']['na_handle'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_auc.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_chipseq_interaction_auc/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction_auc_{method}.log"
  shell:
    """
    Rscript src/eval_interaction_auc.R \
      --net "{input.absnet}" \
      --known "{input.known}" \
      --curve {params.curve} \
      --max {params.max} \
      --min {params.min} \
      --rand {params.rand} \
      --dg {params.dg} \
      --na "{params.na}" \
      --neg "{params.neg}" \
      --o "{output}" \
      2>&1 | tee {log}
    """

rule eval_chipseq_interaction_hub_auc:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction.rds"
  params:
    curve = config['eval_opts']['auc']['curve'],
    max = config['eval_opts']['auc']['max'],
    min = config['eval_opts']['auc']['min'],
    rand = config['eval_opts']['auc']['rand'],
    dg = config['eval_opts']['auc']['dg'],
    na = config['eval_opts']['weight']['na_handle'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_hub_auc.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_chipseq_interaction_hub_auc/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction_hub_auc_{method}.log"
  shell:
    """
    Rscript src/eval_interaction_hub_auc.R \
      --net "{input.absnet}" \
      --known "{input.known}" \
      --curve {params.curve} \
      --max {params.max} \
      --min {params.min} \
      --rand {params.rand} \
      --dg {params.dg} \
      --na "{params.na}" \
      --neg "{params.neg}" \
      --o "{output}" \
      2>&1 | tee {log}
    """

rule eval_chipseq_interaction_spearman_cor:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction.rds"
  params:
    method = "spearman",
    na = config['eval_opts']['weight']['na_handle'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_spearman_cor.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_chipseq_interaction_spearman_cor/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction_spearman_cor_{method}.log"
  shell:
    """
    Rscript src/eval_interaction_cor.R \
      --net "{input.absnet}" \
      --known "{input.known}" \
      --method "{params.method}" \
      --na "{params.na}" \
      --neg "{params.neg}" \
      --o "{output}" \
      2>&1 | tee {log}
    """

rule eval_chipseq_interaction_precision:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction.rds"
  params:
    threshold = ",".join([str(s) for s in config['eval_opts']['precision']['threshold']]),
    top = ",".join([str(s) for s in config['eval_opts']['precision']['top']]),
    na = config['eval_opts']['weight']['na_handle'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_chipseq_interaction_precision.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_chipseq_interaction_precision/{tissue}/{correction_label}/{gene_selection}/{n_genes}/chipseq_interaction_precision_{method}.log"
  shell:
    """
    Rscript src/eval_interaction_precision_at_top.R \
      --net "{input.absnet}" \
      --known "{input.known}" \
      --threshold "{params.threshold}" \
      --top "{params.top}" \
      --na "{params.na}" \
      --neg "{params.neg}" \
      --o "{output}" \
      2>&1 | tee {log}
    """
