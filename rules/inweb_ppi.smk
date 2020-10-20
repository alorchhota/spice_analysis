rule get_inweb_ppi:
  input:
    gene_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/genes.txt",
    inweb_file = "{results_dir}/shared_data/inweb/human_interactions_with_evidence.tab"
  params:
    thread = lambda wildcards, output: get_matrix_size_dependent_default_thread(wildcards)
  output:
    ppi_file = "{results_dir}/gtex_v8/results/{tissue}/{correction_label}/{gene_selection}/{n_genes}/inweb_ppi.rds"
  log: 
    "{results_dir}/gtex_v8/logs/get_inweb_ppi/{correction_label}/{tissue}.{correction_label}.{gene_selection}.{n_genes}.log"
  shell:
    """
    Rscript src/prepare_inweb_ppi_matrix.R \
      --gene "{input.gene_file}" \
      --inweb "{input.inweb_file}" \
      --o "{output.ppi_file}" \
      2>&1 | tee {log}
    """
