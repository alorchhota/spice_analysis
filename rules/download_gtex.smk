rule download_gtex_expr_matrices:
  params:
    gtex_expr_mat_url="https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar"
  output: "{results_dir}/gtex_v8/data/download/GTEx_Analysis_v8_eQTL_expression_matrices.tar" # temp("")
  group: "gtex_expr_download"
  log: 
    "{results_dir}/gtex_v8/logs/download_gtex_expr_matrices/download_gtex_expr_matrices.log"
  shell:
    """
    curl {params.gtex_expr_mat_url} -o {output}
    """

rule extract_gtex_expr_matrices:
  input: expand("{results_dir}/gtex_v8/data/download/GTEx_Analysis_v8_eQTL_expression_matrices.tar", results_dir = config['results_dir'])
  output: expand("{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz", results_dir = config['results_dir'], tissue = config['tissues'])
  group: "gtex_expr_download"
  shell:
    """
    tar -xvf {input} -C "{results_dir}/gtex_v8/data"
    """
    
rule extract_gtex_expr_matrix:
  input: "{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz"
  output: "{results_dir}/gtex_v8/data/normalized/{tissue}.txt"
  group: "gtex_expr_download"
  log: 
    "{results_dir}/gtex_v8/logs/extract_gtex_expr_matrix/{tissue}.log"
  shell:
    """
    gzip -cd "{input}" | cut -f4- > "{output}"
    """

  
rule download_gtex_tpm:
  params:
    gtex_tpm_url="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"
  output: "{results_dir}/gtex_v8/data/download/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz" # temp
  group: "gtex_tpm_download"
  log: 
    "{results_dir}/gtex_v8/logs/download_gtex_tpm/download_gtex_tpm.log"
  shell:
    """
    curl {params.gtex_tpm_url} -o {output}
    """

rule download_gtex_sample_annot:
  params:
    gtex_sample_annot_url="https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
  output: "{results_dir}/gtex_v8/data/download/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt" # temp
  group: "gtex_tpm_download"
  log: 
    "{results_dir}/gtex_v8/logs/download_gtex_sample_annot/download_gtex_sample_annot.log"
  shell:
    """
    curl {params.gtex_sample_annot_url} -o {output}
    """

    
rule extract_gtex_tpm:
  input: 
    tpm=expand("{results_dir}/gtex_v8/data/download/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", results_dir = config['results_dir']),
    annot=expand("{results_dir}/gtex_v8/data/download/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", results_dir = config['results_dir'])
  output: expand("{results_dir}/gtex_v8/data/raw_tpm/{tissue}.txt", results_dir = config['results_dir'], tissue = config['tissues'])
  group: "gtex_tpm_download"
  shell:
    """
    python src/tissuewise_categorize_tpm_data.py \
      -rpkm "{input.tpm}" \
      -ann "{input.annot}" \
      -st 3 \
      -skiprow 2 \
      -o "{results_dir}/gtex_v8/data/raw_tpm"
    """
    
    
rule download_gtex_cov:
  params:
    gtex_cov_url="https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz"
  output: "{results_dir}/gtex_v8/data/download/GTEx_Analysis_v8_eQTL_covariates.tar.gz"
  group: "gtex_cov_download"
  log: 
    "{results_dir}/gtex_v8/logs/download_gtex_cov/download_gtex_cov.log"
  shell:
    """
    curl {params.gtex_cov_url} -o {output}
    """


rule extract_gtex_cov:
  input: expand("{results_dir}/gtex_v8/data/download/GTEx_Analysis_v8_eQTL_covariates.tar.gz", results_dir = config['results_dir'])
  output: expand("{results_dir}/gtex_v8/data/GTEx_Analysis_v8_eQTL_covariates/{tissue}.v8.covariates.txt", results_dir = config['results_dir'], tissue = config['tissues'])
  group: "gtex_cov_download"
  shell:
    """
    tar -zxvf {input} -C "{results_dir}/gtex_v8/data"
    """
    
