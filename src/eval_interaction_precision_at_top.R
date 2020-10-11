library(argparser)
library(spice)
library(miscutil)

args <- arg_parser("program");
args <- add_argument(args, "--net", help="network file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/pearson_network.rds")
args <- add_argument(args, "--known", help="known interaction file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/string_ppi.rds")
args <- add_argument(args, "--threshold", help="minimum known score for an interaction to be considered true. comma-separated values in [0,1].", default="0.7,0.4,0.15")
args <- add_argument(args, "--top", help="number of top edges. comma-separated numeric values. non-numeric values are be converted to NA", default="1000,5000,10000,NA")
args <- add_argument(args, "--na", help="NA handling ('net' / 'known' / 'any').", default="known")
args <- add_argument(args, "--neg", help="negative value handling", default="error")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/string_precision.rds")

### parse args
argv = parse_args(args)
net_fn = argv$net
known_fn = argv$known
score_thresholds_input = argv$threshold
n_top_edges_input = argv$top
na_ignore = argv$na
neg_treat = argv$neg
out_fn = argv$o

score_thresholds = as.numeric(parse_delimitted_string(score_thresholds_input, delim = ',', rm.empty = T))
n_top_edges = suppressWarnings(as.numeric(parse_delimitted_string(n_top_edges_input, delim = ',', rm.empty = T)))

### check variables
stopifnot(na_ignore %in% c("net", "known", "any"))
stopifnot(neg_treat %in% c("none", "warn", "error"))
stopifnot(dir.exists(dirname(out_fn)))
stopifnot(length(n_top_edges) > 0)
stopifnot(length(score_thresholds) > 0)
stopifnot(all(is.finite(score_thresholds)))

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
net_mat = net_mat[common_genes, common_genes]
known_interaction_mat = known_interaction_mat[common_genes, common_genes]

### compute precision
interaction_precision = spice::coexpression_known_interactions_precision_at_top(
  net = net_mat,
  known = known_interaction_mat,
  score.thresholds = score_thresholds,
  n.top.edges = n_top_edges,
  na.ignore = na_ignore,
  neg.treat = neg_treat)

### save
saveRDS(interaction_precision, file = out_fn)
