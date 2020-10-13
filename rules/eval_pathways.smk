def eval_shared_pathway_n_thread(wildcards):
  ng = int(wildcards.n_genes)
  #
  if ng <= 1500:
    th = 1
  elif ng <= 5000:
    th = 2
  else:
    th = 1
  #
  return th


rule eval_shared_pathway:
  input:
    absnet = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    pathway = "{results_dir}/shared_data/msigdb/{pathway}_genesets.rds"
  params:
    curve = config['eval_opts']['auc']['curve'],
    max = config['eval_opts']['auc']['max'],
    min = config['eval_opts']['auc']['min'],
    rand = config['eval_opts']['auc']['rand'],
    dg = config['eval_opts']['auc']['dg'],
    na_rm = config['eval_opts']['weight']['na_rm'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: eval_shared_pathway_n_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_shared_pathway_auc_{pathway}.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_shared_pathway/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_shared_pathway_auc_{pathway}.log"
  shell:
    """
    Rscript src/eval_shared_pathway_auc.R \
      --net "{input.absnet}" \
      --pathway "{input.pathway}" \
      --curve {params.curve} \
      --max {params.max} \
      --min {params.min} \
      --rand {params.rand} \
      --dg {params.dg} \
      --na_rm "{params.na_rm}" \
      --neg "{params.neg}" \
      --o "{output}" \
      2>&1 | tee {log}
    """

rule eval_pathway_enrichment:
  input:
    absnet = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    pathway = "{results_dir}/shared_data/msigdb/{pathway}_genesets.rds"
  params:
    min_gene = config['eval_opts']['pathway_enrichment']['min_gene'],
    max_gene = config['eval_opts']['pathway_enrichment']['max_gene'],
    iter = config['eval_opts']['pathway_enrichment']['iter'],
    seed = config['random_seed'],
    na_rm = config['eval_opts']['weight']['na_rm'],
    neg = config['eval_opts']['weight']['neg_handle'],
    thread = lambda wildcards, output: eval_shared_pathway_n_thread(wildcards)
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_pathway_enrichment_{pathway}.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_pathway_enrichment/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_pathway_enrichment_{pathway}.log"
  shell:
    """
    Rscript src/eval_pathway_enrichment.R \
      --net "{input.absnet}" \
      --pathway "{input.pathway}" \
      --min_gene {params.min_gene} \
      --max_gene {params.max_gene} \
      --iter {params.iter} \
      --seed {params.seed} \
      --na_rm "{params.na_rm}" \
      --neg "{params.neg}" \
      --o "{output}" \
      2>&1 | tee {log}
    """
