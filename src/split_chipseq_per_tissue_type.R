### this script splits surya's chip-seq interactions by the nearest tissue types.
library(ioutil)
library(argparser)

args <- arg_parser("program");
args <- add_argument(args, "--chip", help="file with chipseq interactions prepared by surya", default="/work-zfs/abattle4/surya/datasets/encode/for_genemodels/genepromoter_tfbinding_sites_new_with_nearestTissuegroup.bed")
args <- add_argument(args, "--o", help="splitted file prefix (rds)", default="results/per_tissue/chipseq")

### parse args
argv = parse_args(args)
chipseq_fn = argv$chip
processed_chipseq_pfx = argv$o

### read chipseq by surya
chipseq_df = read_df(chipseq_fn, header = T, row.names = F)
dim(chipseq_df) # [1] 14663683       23
head(chipseq_df, n=2)

### split data per tissue type
nearest_tissue_types = sort(unique(chipseq_df$Nearest_Tissue_Type))
for(tt in nearest_tissue_types){
  chipseq_per_tissue_df = chipseq_df[ chipseq_df$Nearest_Tissue_Type == tt, , drop = F ]
  file_tt = gsub(pattern = "[()_ ]", replacement = " ", x = tt)
  file_tt = gsub(pattern = " +", replacement = "_", x = trimws(file_tt))
  fn = sprintf("%s_%s.rds", processed_chipseq_pfx, file_tt)
  saveRDS(chipseq_per_tissue_df, file = fn)
}
