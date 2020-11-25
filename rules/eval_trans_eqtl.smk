rule process_snps_nearby_genes_for_trans_eqtl:
  input:
    annot = config["gene_annot"],
    genotypes = expand("{geno_pfx}{geno_chr}{geno_sfx}",
      geno_pfx = config['genotype']['prefix'],
      geno_chr = config['genotype']['chromosomes'],
      geno_sfx = config['genotype']['suffix'])
  params:
    geno_pfx = config['genotype']['prefix'],
    geno_sfx = config['genotype']['suffix'],
    geno_chr = ",".join([str(s) for s in config["genotype"]["chromosomes"]]),
    d = config['eval_opts']['trans_eqtl']['max_cis_dist'],
    thread = 3
  output:
    "{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_snps_nearby_gene/snps_nearby_gene.rds"
  log: 
    "{results_dir}/logs/process_snps_nearby_genes_for_trans_eqtl/process_snps_nearby_genes_for_trans_eqtl.log"
  shell:
    """
    Rscript src/process_snps_nearby_genes.R \
      --geno_pfx "{params.geno_pfx}" \
      --geno_sfx "{params.geno_sfx}" \
      --geno_chr "{params.geno_chr}" \
      --gene_annot "{input.annot}" \
      --d {params.d} \
      --o "{output}" \
      2>&1 | tee {log}
    """

rule eval_trans_eqtl_n_egenes:
  input:
    annot = config["gene_annot"],
    genotypes = expand("{geno_pfx}{geno_chr}{geno_sfx}",
      geno_pfx = config['genotype']['prefix'],
      geno_chr = config['genotype']['chromosomes'],
      geno_sfx = config['genotype']['suffix']),
    net = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_absnet.rds",
    expr = "{results_dir}/gtex_v8/data/normalized/{tissue}.txt",
    cov = "{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt",
    gene_snp = "{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_snps_nearby_gene/snps_nearby_gene.rds",
    crossmap = "{results_dir}/shared_data/hg38_cross_mappability.txt"
  params:
    top = config['eval_opts']['trans_eqtl']['top_edge'],
    geno_pfx = config['genotype']['prefix'],
    geno_sfx = config['genotype']['suffix'],
    d = config['eval_opts']['trans_eqtl']['max_crossmap_d'],
    max_snp = config['eval_opts']['trans_eqtl']['max_snp_per_run'],
    repo = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/trans_eqtl_tests_repo.rds",
    thread = 8,
    runtime = 300
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_trans_eqtl_n_egenes.rds"
  log: 
    "{results_dir}/gtex_v8/logs/eval_trans_eqtl_n_egenes/{tissue}/{correction_label}/{gene_selection}/{n_genes}/{method}_eval_trans_eqtl_n_egenes.log"
  shell:
    """
    Rscript src/eval_trans_eqtl.R \
      --net "{input.net}" \
      --top {params.top} \
      --expr "{input.expr}" \
      --cov "{input.cov}" \
      --geno_pfx "{params.geno_pfx}" \
      --geno_sfx "{params.geno_sfx}" \
      --gene_snp "{input.gene_snp}" \
      --annot "{input.annot}" \
      --crossmap "{input.crossmap}" \
      --d {params.d} \
      --repo "{params.repo}" \
      --max_snp {params.max_snp} \
      --o "{output}" \
      2>&1 | tee {log}
    """

