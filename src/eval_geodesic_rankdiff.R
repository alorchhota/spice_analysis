# this script computes difference between mean rank of all geodesic distance pairs.

library(argparser)
library(spice)
library(miscutil)

args <- arg_parser("program");
args <- add_argument(args, "--net", help="network file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/pearson_network.rds")
args <- add_argument(args, "--known", help="known interaction file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/string_ppi.rds")
args <- add_argument(args, "--threshold", help="minimum known score for an interaction to be considered true. range [0,1].", default=0.7)
args <- add_argument(args, "--d", help="known geodesic distances (comma separated)", default="1,2,3,4,5")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/string_geodesic_rankdif.rds")

### parse args
argv = parse_args(args)
net_fn = argv$net
known_fn = argv$known
known_threshold = argv$threshold
dvals_input = argv$d
out_fn = argv$o

dvals = suppressWarnings(as.numeric(parse_delimitted_string(dvals_input, delim = ',', rm.empty = T)))

### check variables
stopifnot(dir.exists(dirname(out_fn)))
stopifnot(length(dvals) > 1)
stopifnot(length(known_threshold) == 1)
stopifnot(all(is.finite(known_threshold)))

### load from inputs
net_mat = readRDS(net_fn)
known_interaction_mat = readRDS(known_fn)

### check network and known matrix
stopifnot(is.matrix(net_mat))
stopifnot(is.matrix(known_interaction_mat))
stopifnot(length(rownames(net_mat))>0)
stopifnot(length(rownames(known_interaction_mat))>0)
if(!all(rownames(known_interaction_mat) %in% rownames(net_mat)))
  stop("known interaction matrix must contain genes only network matrix only.")

### take common genes in net and ppi for evaluation with string
common_genes = intersect(rownames(known_interaction_mat), rownames(net_mat))
net_mat = abs(net_mat[common_genes, common_genes])
known_interaction_mat = known_interaction_mat[common_genes, common_genes]

### compute rank differences
mean_ranks = spice::mean_rank_for_known_geodesic_distance(
  net = net_mat,
  known = known_interaction_mat, d = dvals, known.threshold = known_threshold)

rankdiffs = sapply(mean_ranks, function(mr1){
  sapply(mean_ranks, function(mr2){
    return(mr1 - mr2)
  })
})
rownames(rankdiffs) = colnames(rankdiffs) = as.character(dvals)

### save
saveRDS(rankdiffs, file = out_fn)
