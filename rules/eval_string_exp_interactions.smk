rule get_string_exp_ppi:
  input:
    gene_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/genes.txt",
    string_link_file = "{results_dir}/shared_data/string.9606/9606.protein.links.detailed.v11.0.txt.gz",
    string_info_file = "{results_dir}/shared_data/string.9606/9606.protein.info.v11.0.txt.gz"
  params:
    directed = "FALSE",
    norm = "TRUE",
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    ppi_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_ppi.rds"
  log: 
    "{results_dir}/gtex_v8/logs/get_string_exp_ppi/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.log"
  shell:
    """
    Rscript src/prepare_string_experimental_matrix.R \
      --gene "{input.gene_file}" \
      --link "{input.string_link_file}" \
      --info "{input.string_info_file}" \
      --directed "{params.directed}" \
      --norm "{params.norm}" \
      --o "{output.ppi_file}" \
      2>&1 | tee {log}
    """

rule eval_string_exp_ppi_auc:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_ppi.rds"
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
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_auc.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_string_exp_ppi_auc/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_auc_{method}.log"
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

rule eval_string_exp_ppi_hub_auc:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_ppi.rds"
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
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_hub_auc.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_string_exp_ppi_hub_auc/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_hub_auc_{method}.log"
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

rule eval_string_exp_ppi_spearman_cor:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_ppi.rds"
  params:
    method = "spearman",
    na = config['eval_opts']['weight']['na_handle'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_spearman_cor.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_string_exp_ppi_spearman_cor/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_spearman_cor_{method}.log"
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

rule eval_string_exp_ppi_precision:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_ppi.rds"
  params:
    threshold = ",".join([str(s) for s in config['eval_opts']['precision']['threshold']]),
    top = ",".join([str(s) for s in config['eval_opts']['precision']['top']]),
    na = config['eval_opts']['weight']['na_handle'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_precision.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_string_exp_ppi_precision/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_precision_{method}.log"
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

rule eval_string_exp_ppi_geodesic_rankdiff:
  input:
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    known = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_ppi.rds"
  params:
    threshold = config['eval_opts']['geodesic_rankdiff']['string_threshold'],
    d = ",".join([str(s) for s in config['eval_opts']['geodesic_rankdiff']['d']])
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_string_exp_ppi_geodesic_rankdiff.rds"
  log:
    "{results_dir}/gtex_v8/logs/eval_string_exp_ppi_geodesic_rankdiff/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_exp_ppi_geodesic_rankdiff_{method}.log"
  shell:
    """
    Rscript src/eval_geodesic_rankdiff.R \
      --net "{input.absnet}" \
      --known "{input.known}" \
      --threshold {params.threshold} \
      --d "{params.d}" \
      --o "{output}" \
      2>&1 | tee {log}
    """
