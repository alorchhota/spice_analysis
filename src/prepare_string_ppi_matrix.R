library(spice)
library(argparser)
library(ioutil)

args <- arg_parser("program");
args <- add_argument(args, "--gene", help="file with gene names", default="/Users/ashissaha/github/spice_analysis/results/spice_results/shared_data/dummy_genes_for_string_download.txt")
args <- add_argument(args, "--version", help="STRING version", default="11")
args <- add_argument(args, "--directed", help="Directed edge?", default=FALSE)
args <- add_argument(args, "--species", help="Species id in STRING", default=9606)
args <- add_argument(args, "--threshold", help="STRING score threshold [0-999]", default=0)
args <- add_argument(args, "--dir", help="Directory to save STRING db, empty string to use temp direactory", default="")
args <- add_argument(args, "--homology", help="Maximum homology score, Inf to get all edges", default=Inf)
args <- add_argument(args, "--norm", help="Return normalized score?", default=TRUE)
args <- add_argument(args, "--o", help="Output file (rds)", default="results/string.rds")

### parse args
argv = parse_args(args)
gene_fn = argv$gene
string_version = argv$version
directed = as.logical(argv$directed)
species = argv$species
score_threshold = argv$threshold
string_dir = argv$dir
max_homology_score = argv$homology
normalize_score = as.logical(argv$norm)
out_fn = argv$o

### get gene names
gene_df = read_df(gene_fn, header = F, row.names = F, check.names = F)
genes = as.character(gene_df[,1])

### get ppi scores
string_ppi_mat = spice::get_string_ppi_matrix(
  genes = genes,
  version = string_version, 
  directed = directed, 
  species = species, 
  score.threshold = score_threshold, 
  string.dir = string_dir, 
  max.homology.bitscore = max_homology_score, 
  normalize.score = normalize_score)

### save
saveRDS(string_ppi_mat, file = out_fn)
