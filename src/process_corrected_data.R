# process data to replace ensembl gene id by gene symbol
library(ioutil)
library(miscutil)
library(genomicsutil)
library(argparser)

args <- arg_parser("program");
args <- add_argument(args, "--expr", help="expression data", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_analysis/data/corrected_expression/Muscle_Skeletal.v7.corrected.txt")
args <- add_argument(args, "--count", help="read count data", default="/work-zfs/abattle4/lab_data/GTEx_v7/rna_seq_from_portal/gene_count/Muscle-Skeletal.txt")
args <- add_argument(args, "--annot", help="gene annotation data", default="/work-zfs/abattle4/ashis/progres/misc/cross_mappability/manuscript_data/misc/gencode.v19.annotation.gene.txt")
args <- add_argument(args, "--n", help="number of genes", default=500)
args <- add_argument(args, "--select", help="VAR/CV/AXPVAR/AXPCV for variance/coefficient of variation, AXP stands for protein coding genes in autsomal and X chr", default='CV')
args <- add_argument(args, "--qnorm", help="per-gene quantile normalize?", default=TRUE)
args <- add_argument(args, "--o", help="output file", default="results/processed.txt")

argv = parse_args(args)
expr_fn = argv$expr
count_fn = argv$count
n_genes = argv$n
gene_selection_method = argv$select
do_qnorm = as.logical(argv$qnorm)
annot_fn = argv$annot
out_expr_fn = argv$o

stopifnot(gene_selection_method %in% c('VAR', 'CV', 'AXPVAR', 'AXPCV'))

expr_df = read_df(expr_fn)
count_df = read_df(count_fn)
annot_df = read_df(annot_fn)

### set colnames of count data after parsing
cn = colnames(count_df)
cn = sapply(cn, function(s) paste(strsplit(s, split = '-')[[1]][1:2], sep = '-', collapse = '-') )
colnames(count_df) = cn
stopifnot(all(rownames(expr_df) %in% rownames(count_df)))
stopifnot(all(colnames(expr_df) %in% colnames(count_df)))
count_df = count_df[rownames(expr_df), colnames(expr_df)]

### map ensembl gene ids to gene symbols
genes = rownames(expr_df)
symbols = annot_df[genes, 'gene_name']
symbols_counts = table(symbols)
multi_symbols = names(symbols_counts[symbols_counts>1])
symbols_to_take = !(symbols %in% c(multi_symbols, NA))

new_expr_df = expr_df[symbols_to_take, ,drop = F]
rownames(new_expr_df) = symbols[symbols_to_take]
new_count_df = count_df[symbols_to_take, ,drop = F]
rownames(new_count_df) = symbols[symbols_to_take]

new_genes = genes[symbols_to_take]
new_annot_df = annot_df[new_genes,,drop=F]
rownames(new_annot_df) = new_annot_df$gene_name 

### take top n genes
if(gene_selection_method == 'CV'){
  selected_expr_df = filter_expr_by_coeff_of_variation(expr.df = new_expr_df, raw.df = new_count_df, n = n_genes, min.var = 1e-6, min.mean = -Inf)
} else if(gene_selection_method == 'VAR'){
  selected_expr_df = filter_expr_by_variance(expr.df = new_expr_df, raw.df = new_count_df, n = n_genes, min.var = 1e-6, min.mean = -Inf)
} else if(gene_selection_method == 'AXPVAR'){
  selected_expr_df = filter_expr_by_chr(expr.df = new_expr_df, 
                                        chr.include = as.character(c(1:22,'X')), 
                                        annot.gene = new_annot_df, 
                                        chr.col = "chr")
  selected_expr_df = filter_expr_by_gene_type(expr.df = selected_expr_df, 
                                              annot.gene = new_annot_df, 
                                              type.col = "gene_type", 
                                              type.values = "protein_coding")
  gene.variance = as.numeric(apply(new_count_df[rownames(selected_expr_df),], MARGIN = 1, var))
  selected_expr_df = filter_expr_by_coordinate_overlap(expr.df = selected_expr_df, 
                                                       annot.gene = new_annot_df, 
                                                       priority.gene = gene.variance, 
                                                       chr.col = "chr", 
                                                       pos.start.col = "start_pos", 
                                                       pos.end.col = "end_pos")
  selected_expr_df = filter_expr_by_variance(expr.df = selected_expr_df, 
                                             raw.df = new_count_df[rownames(selected_expr_df), ,drop=F], 
                                             n = n_genes, 
                                             min.var = 1e-6, 
                                             min.mean = -Inf)
} else if(gene_selection_method == 'AXPCV'){
  selected_expr_df = filter_expr_by_chr(expr.df = new_expr_df, 
                                        chr.include = as.character(c(1:22,'X')), 
                                        annot.gene = new_annot_df, 
                                        chr.col = "chr")
  selected_expr_df = filter_expr_by_gene_type(expr.df = selected_expr_df, 
                                              annot.gene = new_annot_df, 
                                              type.col = "gene_type", 
                                              type.values = "protein_coding")
  gene.sd = as.numeric(apply(new_count_df[rownames(selected_expr_df),], MARGIN = 1, sd))
  gene.mean = as.numeric(apply(new_count_df[rownames(selected_expr_df),], MARGIN = 1, mean))
  gene.cv = gene.sd / gene.mean
  gene.cv[is.infinite(gene.cv)] = 0
  selected_expr_df = filter_expr_by_coordinate_overlap(expr.df = selected_expr_df, 
                                                       annot.gene = new_annot_df, 
                                                       priority.gene = gene.cv, 
                                                       chr.col = "chr", 
                                                       pos.start.col = "start_pos", 
                                                       pos.end.col = "end_pos")
  selected_expr_df = filter_expr_by_coeff_of_variation(expr.df = selected_expr_df, 
                                                       raw.df = new_count_df[rownames(selected_expr_df), ,drop=F], 
                                                       n = n_genes, 
                                                       min.var = 1e-6, 
                                                       min.mean = -Inf)
} else {
  stop('gene selection not implemented!')
}

### per-gene quantile normalize
if(do_qnorm == T)
  selected_expr_df = to_inv_normal(as.matrix(selected_expr_df))

write_df(selected_expr_df, file = out_expr_fn)
