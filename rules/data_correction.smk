rule rm_cov_from_expr:
  input:
    expr="{results_dir}/gtex_v8/data/normalized/{tissue}.txt",
    cov="{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt"
  output:
    "{results_dir}/gtex_v8/data/covariates_removed/{tissue}.txt"
  group: "correct_data"
  shell:
    """
    Rscript src/rm_cov_from_expr.R \
      --expr "{input.expr}" \
      --cov "{input.cov}" \
      --id 1 \
      --st 2 \
      --o "{output}"
    """

rule process_corrected_expr:
  input:
    expr="{results_dir}/gtex_v8/data/covariates_removed/{tissue}.txt",
    tpm="{results_dir}/gtex_v8/data/raw_tpm/{tissue}.txt",
    annot=config["gene_annot"]
  output:
    "{results_dir}/gtex_v8/data/corrected/{tissue}.corrected.{gene_selection}.{n_genes}.txt"
  shell:
    """
    Rscript src/process_corrected_data.R \
      --expr "{input.expr}" \
      --count "{input.tpm}" \
      --annot "{input.annot}" \
      --n {wildcards.n_genes} \
      --select "{wildcards.gene_selection}" \
      --qnorm TRUE \
      --o {output}
    """
