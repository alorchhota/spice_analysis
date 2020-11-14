results_dir = config["results_dir"]
validation_tissues = config["validation_tissues"]
test_tissues = config["test_tissues"]
chipseq_tissue_map = config["chipseq_tissue_map"]

validation_chipseq_tissues = [t for t in validation_tissues 
                              if t in chipseq_tissue_map 
                              and chipseq_tissue_map[t] != ""]

test_chipseq_tissues = [t for t in test_tissues 
                        if t in chipseq_tissue_map 
                        and chipseq_tissue_map[t] != ""]

chipseq_tissue_types = [chipseq_tissue_map[t] for t in chipseq_tissue_map 
                        if chipseq_tissue_map[t] != ""]
chipseq_tissue_types = sorted(list(set(chipseq_tissue_types)))
