library(argparser)
library(ioutil)
library(reshape2)

args <- arg_parser("program");
args <- add_argument(args, "--gene", help="file with gene names", default="/work-zfs/abattle4/ashis/progres/spice_anlysis/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/5000/genes.txt")
args <- add_argument(args, "--chip", help="chipseq interaction file", default="/work-zfs/abattle4/ashis/progres/spice_anlysis/shared_data/chipseq/per_tissue/chipseq_Blood.rds")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/chipmat.rds")

### parse args
argv = parse_args(args)
gene_fn = argv$gene
chipseq_fn = argv$chip
out_fn = argv$o

### get gene names
gene_df = read_df(gene_fn, header = F, row.names = F, check.names = F)
genes = as.character(gene_df[,1])

### read chipseq info
chipseq_df = readRDS(chipseq_fn)

### construct chipseq matrix
chipseq_df = chipseq_df[(chipseq_df$Target %in% genes) &
                          (chipseq_df$gene_name %in% genes), , drop = T]
chip_signal_val_mat = suppressWarnings(
  acast(
    data = chipseq_df,
    formula = Target ~ gene_name,
    value.var = "signal_val",
    fun.aggregate = max,
    na.rm = T
  )
)
chip_signal_val_mat[is.infinite(chip_signal_val_mat)] = 0
# convert to 0/1 matrix
chip_signal_val_mat[chip_signal_val_mat > 0] = 1

### create a matrix with all genes
chip_mat = matrix(
  NA,
  nrow = length(genes),
  ncol = length(genes),
  dimnames = list(genes, genes)
)
chip_mat[rownames(chip_signal_val_mat), colnames(chip_signal_val_mat)] = chip_signal_val_mat
chip_mat = pmax(chip_mat, t(chip_mat), na.rm = T)
chip_mat[is.na(chip_mat)] = 0

### save
saveRDS(chip_mat, file = out_fn)
