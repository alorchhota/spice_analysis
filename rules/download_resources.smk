rule download_crossmap_hg38:
  params:
    crossmap_url="https://ndownloader.figshare.com/files/13514741"
  output: 
    "{results_dir}/shared_data/hg38_cross_mappability.txt.gz"
  group: 
    "resource_download"
  log: 
    "{results_dir}/gtex_v8/logs/download_crossmap/download_crossmap_hg38.log"
  shell:
    """
    curl {params.crossmap_url} -o {output}
    """

rule extract_crossmap_hg38:
  input: 
    "{results_dir}/shared_data/hg38_cross_mappability.txt.gz"
  output:
    "{results_dir}/shared_data/hg38_cross_mappability.txt"
  group: 
    "resource_download"
  log: 
    "{results_dir}/gtex_v8/logs/extract_crossmap_hg38/extract_crossmap_hg38.log"
  shell:
    """
    gzip -cd "{input}" > "{output}"
    """