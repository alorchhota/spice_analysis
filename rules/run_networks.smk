def select_n_thread(wildcards):
  ng = int(wildcards.n_genes)
  met = wildcards.method
  #
  if ng <= 1500:
    th = 1
    if 'genie3' in met:
      th = 5
    elif 'spice' in met:
      if '_ME_' in met:
        th = 5
  elif ng <= 5000:
    th = 4
    if 'genie3' in met:
      th = 6
    elif 'spice' in met:
      th = 6
  else:
    th = 1
  #
  return th

def select_runtime(wildcards):
  ng = int(wildcards.n_genes)
  met = wildcards.method
  #
  if ng <= 1500:
    rt = 30     # 3p min
    if 'genie3' in met:
      rt = 4 * 60
    elif 'glasso' in met:
      rt = 2 * 60
    elif 'spice' in met:
      if '_ME_' in met:
        rt = 7 * 60
  elif ng <= 5000:
    rt = 4 * 60
    if 'genie3' in met:
      rt = 16 * 60
    elif 'glasso' in met:
      rt = 16 * 60
    elif 'spice' in met:
      rt = 4 * 60
  else:
    rt = 4 * 60
  #
  return rt


rule run_gtex_network:
  input:
    expr = "{results_dir}/gtex_v8/data/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.txt",
    annot = config["gene_annot"],
    cross = "{results_dir}/shared_data/hg38_cross_mappability.txt"
  params:
    method = "{method}",
    cross_th = config["crossmap_threshold"],
    thread = lambda wildcards, output: select_n_thread(wildcards),
    runtime = lambda wildcards, output: select_runtime(wildcards),
    outdir = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}"
  output:
    net="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_network.rds",
    absnet="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    time="{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_time.rds"
  log: 
    "{results_dir}/gtex_v8/logs/run_gtex_network/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}.log"
  shell:
    """
    Rscript src/run_network.R \
      --expr "{input.expr}" \
      --method "{params.method}" \
      --thread {params.thread} \
      --annot "{input.annot}" \
      --cross "{input.cross}" \
      --cross_th {params.cross_th} \
      --old_glasso TRUE \
      --o "{params.outdir}" 2>&1 | tee {log}
    """
