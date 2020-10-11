library(argparser)
library(spice)

args <- arg_parser("program");
args <- add_argument(args, "--net", help="network file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/pearson_network.rds")
args <- add_argument(args, "--known", help="known interaction file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/string_ppi.rds")
args <- add_argument(args, "--curve", help="Should the curves be returned?", default=TRUE)
args <- add_argument(args, "--max", help="compute max auc?", default=TRUE)
args <- add_argument(args, "--min", help="compute min auc?", default=TRUE)
args <- add_argument(args, "--rand", help="compute random auc?", default=FALSE)
args <- add_argument(args, "--dg", help="compute auc using Davis and Goadrich algorithm?", default=FALSE)
args <- add_argument(args, "--na", help="NA handling ('net' / 'known' / 'any').", default="known")
args <- add_argument(args, "--neg", help="negative value handling", default="error")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/string_auc.rds")

### parse args
argv = parse_args(args)
net_fn = argv$net
known_fn = argv$known
curve = as.logical(argv$curve)
max_compute = as.logical(argv$max)
min_compute = as.logical(argv$min)
rand_compute = as.logical(argv$rand)
dg_compute = as.logical(argv$dg)
na_ignore = argv$na
neg_treat = argv$neg
out_fn = argv$o

### check variables
stopifnot(na_ignore %in% c("net", "known", "any"))
stopifnot(neg_treat %in% c("none", "warn", "error"))
stopifnot(dir.exists(dirname(out_fn)))

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

### compute auc
interaction_auc = spice::coexpression_known_interactions_auc(
  net = net_mat, 
  known = known_interaction_mat, 
  curve = curve, 
  max.compute = max_compute, 
  min.compute = min_compute, 
  rand.compute = rand_compute,
  dg.compute = dg_compute,
  na.ignore = na_ignore, 
  neg.treat = neg_treat)

### save
saveRDS(interaction_auc, file = out_fn)
