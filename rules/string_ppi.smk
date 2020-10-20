# prepare dummy interactions to download the database into local directory
rule download_string_ppi:
  params:
    version = config['string']['version'],
    gene_file = "{results_dir}/shared_data/dummy_genes_for_string_download.txt",
    string_dir = "{results_dir}/shared_data/string",
  output:
    ppi_file = "{results_dir}/shared_data/dependency_ppi_{species}.rds"
  log: 
    "{results_dir}/logs/download_string_ppi/download_string_ppi_{species}.log"
  shell:
    """
    printf "TP53\nRBM3\nSF3\nLIM12\nATM\nTMEM160\nBCL2L1\nMDM2\nEGFR\nCD96\nKEAP1\nSRSF1\nTSEN2" > "{params.gene_file}"
    mkdir -p "{params.string_dir}"
    Rscript src/prepare_string_ppi_matrix.R \
      --gene "{params.gene_file}" \
      --version "{params.version}" \
      --species "{wildcards.species}"\
      --dir "{params.string_dir}" \
      --o "{output.ppi_file}" \
      2>&1 | tee {log}
    rm "{params.gene_file}"
    """

rule generate_gene_file:
  input:
    "{results_dir}/gtex_v8/data/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.txt"
  output:
    "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/genes.txt"
  log:
    "{results_dir}/gtex_v8/logs/generate_gene_file/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.generate_gene_file.log"
  shell:
    """
    cut -f1 "{input}" | tail -n +2 > "{output}"
    """
    
rule get_string_ppi:
  input:
    gene_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/genes.txt",
    dependency_ppi_file = expand("{results_dir}/shared_data/dependency_ppi_{species}.rds", results_dir = "{results_dir}", species=config['string']['species'])
  params:
    version = config['string']['version'],
    species = config['string']['species'],
    directed = "FALSE",
    score_threshold = 0,
    string_dir = "{results_dir}/shared_data/string",
    norm = "TRUE",
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    ppi_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_ppi.rds"
  log: 
    "{results_dir}/gtex_v8/logs/get_string_ppi/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.log"
  shell:
    """
    Rscript src/prepare_string_ppi_matrix.R \
      --gene "{input.gene_file}" \
      --version "{params.version}" \
      --species "{params.species}"\
      --dir "{params.string_dir}" \
      --o "{output.ppi_file}" \
      2>&1 | tee {log}
    """

rule get_string_kegg_ppi:
  input:
    gene_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/genes.txt",
    kegg_file = "{results_dir}/shared_data/msigdb/kegg_genesets.rds",
    dependency_ppi_file = expand("{results_dir}/shared_data/dependency_ppi_{species}.rds", results_dir = "{results_dir}", species=config['string']['species'])
  params:
    version = config['string']['version'],
    species = config['string']['species'],
    directed = "FALSE",
    score_threshold = 0,
    string_dir = "{results_dir}/shared_data/string",
    norm = "TRUE",
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    ppi_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/string_kegg_ppi.rds"
  log: 
    "{results_dir}/gtex_v8/logs/get_string_kegg_ppi/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.log"
  shell:
    """
    Rscript src/prepare_string_ppi_matrix.R \
      --gene "{input.gene_file}" \
      --version "{params.version}" \
      --directed "{params.directed}" \
      --species "{params.species}" \
      --threshold "{params.score_threshold}" \
      --dir "{params.string_dir}" \
      --pathway "{input.kegg_file}" \
      --norm {params.norm} \
      --o "{output.ppi_file}" \
      2>&1 | tee {log}
    """
