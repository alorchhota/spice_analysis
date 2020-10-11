library(argparser)
library(spice)

args <- arg_parser("program");
args <- add_argument(args, "--net", help="network file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/pearson_network.rds")
args <- add_argument(args, "--known", help="known interaction file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/string_ppi.rds")
args <- add_argument(args, "--method", help="method to compute correlation ('pearson' / 'spearman' / 'kendall')", default="spearman")
args <- add_argument(args, "--na", help="NA handling ('net' / 'known' / 'any').", default="known")
args <- add_argument(args, "--neg", help="negative value handling", default="error")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/string_cor.rds")

### parse args
argv = parse_args(args)
net_fn = argv$net
known_fn = argv$known
method = argv$method
na_ignore = argv$na
neg_treat = argv$neg
out_fn = argv$o

### check variables
stopifnot(method %in% c("spearman", "pearson", "kendall"))
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

### compute cor
interaction_cor = spice::coexpression_known_interactions_cor(
  net = net_mat,
  known = known_interaction_mat,
  method = method,
  na.ignore = na_ignore,
  neg.treat = neg_treat)

### save
saveRDS(interaction_cor, file = out_fn)
