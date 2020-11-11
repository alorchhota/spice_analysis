library(argparser)
library(ioutil)
library(reshape2)

args <- arg_parser("program");
args <- add_argument(args, "--gene", help="file with gene names", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/genes.txt")
args <- add_argument(args, "--kegg", help="file with kegg interactions", default="/Users/ashissaha/github/spice_analysis/results/spice_results/shared_data/kegg/kegg_interactions_unique.txt")
args <- add_argument(args, "--directed", help="Directed edge?", default=FALSE)
args <- add_argument(args, "--o", help="Output file (rds)", default="results/kegg.rds")

### parse args
argv = parse_args(args)
gene_fn = argv$gene
kegg_fn = argv$kegg
directed = as.logical(argv$directed)
out_fn = argv$o

### get gene names
gene_df = read_df(gene_fn, header = F, row.names = F, check.names = F)
genes = as.character(gene_df[,1])

### get kegg interactions
kegg_df = read_df(kegg_fn, header = T, row.names = F)
kegg_df = kegg_df[(kegg_df$src %in% genes) & (kegg_df$dest %in% genes), , drop = F]
kegg_df$score = 1

### covert to a matrix
score_mat = suppressWarnings(
  reshape2::acast(
    data = kegg_df,
    formula = src ~ dest,
    value.var = "score",
    fun.aggregate = max,
    na.rm = T
  )
)
score_mat[!is.finite(score_mat)] = 0

kegg_mat = matrix(0, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
kegg_mat[rownames(score_mat), colnames(score_mat)] = score_mat
if (!directed) 
  kegg_mat = pmax(kegg_mat, t(kegg_mat), na.rm = T)

### save
saveRDS(kegg_mat, file = out_fn)
