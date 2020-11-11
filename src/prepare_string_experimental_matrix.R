# this script creates a gene x gene ppi matrix using only experimental scores 
# (not combined score) from STRING.

library(argparser)
library(ioutil)
library(reshape2)

args <- arg_parser("program");
args <- add_argument(args, "--gene", help="file with gene names", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/genes.txt")
args <- add_argument(args, "--link", help="file with STRING links", default="/Users/ashissaha/github/spice_analysis/results/spice_results/shared_data/string.9606/9606.protein.links.detailed.v11.0.txt.gz")
args <- add_argument(args, "--info", help="file with STRING info", default="/Users/ashissaha/github/spice_analysis/results/spice_results/shared_data/string.9606/9606.protein.info.v11.0.txt.gz")
args <- add_argument(args, "--directed", help="Directed edge?", default=FALSE)
args <- add_argument(args, "--norm", help="Return normalized score?", default=TRUE)
args <- add_argument(args, "--o", help="Output file (rds)", default="results/string_experimantal.rds")

### parse args
argv = parse_args(args)
gene_fn = argv$gene
link_fn = argv$link
info_fn = argv$info
directed = as.logical(argv$directed)
normalize_score = as.logical(argv$norm)
out_fn = argv$o

### get gene names
gene_df = read_df(gene_fn, header = F, row.names = F, check.names = F)
genes = as.character(gene_df[,1])

### read string data
link_df = read_df(link_fn, header = T, row.names = F, sep = " ")
info_df = read_df(info_fn, header = T, row.names = F, sep = "\t")

### subset string data to contain given genes only
info_df = info_df[info_df$preferred_name %in% genes, c("protein_external_id", "preferred_name") , drop = F]
colnames(info_df) = c('protein', 'gene')
proteins = unique(c(info_df$protein))
link_df = link_df[link_df$experimental>0, c("protein1", "protein2", "experimental"), drop = F]
link_df = link_df[(link_df$protein1 %in% proteins) & (link_df$protein2 %in% proteins), , drop = F]
link_df = merge(link_df, info_df, by.x = "protein1", by.y = "protein")
link_df = merge(link_df, info_df, by.x = "protein2", by.y = "protein", suffixes = c("1", "2"))

### covert to a matrix
score_mat = suppressWarnings(
  reshape2::acast(
    data = link_df,
    formula = gene1 ~ gene2,
    value.var = "experimental",
    fun.aggregate = max,
    na.rm = T
  )
)
score_mat[!is.finite(score_mat)] = 0
if(normalize_score){
  score_mat = score_mat/1000
}

full_mat = matrix(0, nrow = length(genes), ncol = length(genes), dimnames = list(genes, genes))
full_mat[rownames(score_mat), colnames(score_mat)] = score_mat
if (!directed) 
  full_mat = pmax(full_mat, t(full_mat), na.rm = T)

### save
saveRDS(full_mat, file = out_fn)
