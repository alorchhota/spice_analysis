library(msigdbr)
library(argparser)

args <- arg_parser("program");
args <- add_argument(args, "--species", help="species for msigdb (see msigdbr::msigdbr_show_species())", default="Homo sapiens")
args <- add_argument(args, "--cat", help="category (e.g., 'H', 'C2', etc.)", default="H")
args <- add_argument(args, "--subcat", help="sub-category (e.g., '', 'CP:KEGG', etc.)", default="")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/msigdb.rds")

### parse args
argv = parse_args(args)
pathway_db = argv$db
species = argv$species
cat = argv$cat
subcat = argv$subcat
out_fn = argv$o

### get pathways
msigdb_df = as.data.frame(msigdbr(species = species, category = cat, subcategory = subcat))
pathways = tapply(msigdb_df$human_gene_symbol, msigdb_df$gs_name, FUN = c)

### save
saveRDS(pathways, file = out_fn)
