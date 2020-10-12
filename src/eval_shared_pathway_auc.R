library(argparser)
library(spice)

args <- arg_parser("program");
args <- add_argument(args, "--net", help="network file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/pearson_network.rds")
args <- add_argument(args, "--pathway", help="pathway file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/shared_data/msigdb/kegg_genesets.rds")
args <- add_argument(args, "--curve", help="Should the curves be returned?", default=TRUE)
args <- add_argument(args, "--max", help="compute max auc?", default=TRUE)
args <- add_argument(args, "--min", help="compute min auc?", default=TRUE)
args <- add_argument(args, "--rand", help="compute random auc?", default=FALSE)
args <- add_argument(args, "--dg", help="compute auc using Davis and Goadrich algorithm?", default=FALSE)
args <- add_argument(args, "--na_rm", help="remove NA?", default=FALSE)
args <- add_argument(args, "--neg", help="negative value handling", default="error")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/shared_pathway_auc.rds")

### parse args
argv = parse_args(args)
net_fn = argv$net
pathway_fn = argv$pathway
curve = as.logical(argv$curve)
max_compute = as.logical(argv$max)
min_compute = as.logical(argv$min)
rand_compute = as.logical(argv$rand)
dg_compute = as.logical(argv$dg)
na_rm = as.logical(argv$na_rm)
neg_treat = argv$neg
out_fn = argv$o

### check variables
stopifnot(file.exists(net_fn))
stopifnot(file.exists(pathway_fn))
stopifnot(neg_treat %in% c("none", "warn", "error"))
stopifnot(dir.exists(dirname(out_fn)))

### load from inputs
net_mat = readRDS(net_fn)
pathways = readRDS(pathway_fn)

### check network
stopifnot(is.matrix(net_mat))
stopifnot(length(rownames(net_mat))>0)

### compute shared-pathway auc
shared_pathway_auc = spice::coexpression_shared_pathway_auc(
  net = net_mat, 
  pathways = pathways,
  curve = curve,
  max.compute = max_compute,
  min.compute = min_compute,
  rand.compute = rand_compute,
  dg.compute = dg_compute,
  na.rm = na_rm,
  neg.treat = neg_treat)

### save
saveRDS(shared_pathway_auc, file = out_fn)
