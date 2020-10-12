library(argparser)
library(spice)

args <- arg_parser("program");
args <- add_argument(args, "--net", help="network file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500/pearson_network.rds")
args <- add_argument(args, "--pathway", help="pathway file (rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/shared_data/msigdb/kegg_genesets.rds")
args <- add_argument(args, "--min_gene", help="minimum number of genes in a pathway", default=5)
args <- add_argument(args, "--max_gene", help="maximum number of genes in a pathway", default=100)
args <- add_argument(args, "--iter", help="number of iterations", default=10000)
args <- add_argument(args, "--seed", help="numeric random number generation seed. non-numeric to avoid seeding.", default="NA")
args <- add_argument(args, "--na_rm", help="remove NA?", default=FALSE)
args <- add_argument(args, "--neg", help="negative value handling", default="error")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/pathway_enrichment.rds")

### parse args
argv = parse_args(args)
net_fn = argv$net
pathway_fn = argv$pathway
min_gene = as.numeric(argv$min_gene)
max_gene = as.numeric(argv$max_gene)
iter = as.numeric(argv$iter)
seed = suppressWarnings(as.numeric(argv$seed))
na_rm = as.logical(argv$na_rm)
neg_treat = argv$neg
out_fn = argv$o

if(!is.finite(seed)){
  seed = NULL
}

### check variables
stopifnot(file.exists(net_fn))
stopifnot(file.exists(pathway_fn))
stopifnot((min_gene > 0) & (max_gene >= min_gene) )
stopifnot(neg_treat %in% c("none", "warn", "error"))
stopifnot(dir.exists(dirname(out_fn)))

### load from inputs
net_mat = readRDS(net_fn)
pathways = readRDS(pathway_fn)

### check network
stopifnot(is.matrix(net_mat))
stopifnot(length(rownames(net_mat))>0)

### compute shared-pathway auc
pathway_enrich = spice::coexpression_pathway_enrichment(
  net = net_mat,
  pathways = pathways,
  min.gene = min_gene,
  max.gene = max_gene,
  iter = iter,
  seed = seed,
  na.rm = na_rm,
  neg.treat = neg_treat)

### save
saveRDS(pathway_enrich, file = out_fn)
