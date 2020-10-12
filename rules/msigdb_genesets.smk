rule get_msigdb_pathway:
  params:
    species =  config['msigdb_species'],
    cat = lambda wildcards: config["pathways"][wildcards.pathway]["cat"],
    subcat = lambda wildcards: config["pathways"][wildcards.pathway]["subcat"]
  output:
    pathway_file = "{results_dir}/shared_data/msigdb/{pathway}_genesets.rds"
  log: 
    "{results_dir}/logs/get_msigdb_pathway/{pathway}.log"
  shell:
    """
    Rscript src/prepare_pathway_geneset.R \
      --species "{params.species}" \
      --cat "{params.cat}" \
      --subcat "{params.subcat}"\
      --o "{output.pathway_file}" \
      2>&1 | tee "{log}"
    """
    