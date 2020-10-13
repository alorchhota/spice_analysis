library(feather)
library(argparser)
library(ioutil)
library(genomicsutil)
library(mappabilityutil)
library(miscutil)
library(reshape2)
source('src/coexpression_network.R')
library(spice)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "--expr", help="expression data", default="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/v8_cis_eqtl_expr_corrected_sym/Whole_Blood.v8.corrected.VAR.500.txt")
args <- add_argument(args, "--method", help="method to run", default="pearson")
args <- add_argument(args, "--thread", help="number of threads (cores)", default=1)
args <- add_argument(args, "--annot", help="gene annotation file", default="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/gencode.v26.annotation.gene.txt")
args <- add_argument(args, "--cross", help="crossmappability file", default="/work-zfs/abattle4/lab_data/annotation/mappability_hg38_gencode26/hg38_cross_mappability_strength_symmetric_mean_sorted.txt")
args <- add_argument(args, "--cross_th", help="crossmappability threshold", default=10)
args <- add_argument(args, "--old_glasso", help="if TRUE, use saved graphical lasso networks from file (if any).", default=TRUE)
args <- add_argument(args, "--seed", help="random seed number. non-numeric to avoid setting seed.", default="NA")
args <- add_argument(args, "--o", help="output directory", default="results/outdir")

argv = parse_args(args)
expr_fn = argv$expr
annot_fn = argv$annot
crossmap_fn = argv$cross
crossmap_threshold = argv$cross_th
out_dir = argv$o
method = argv$method
use_old_glasso = argv$old_glasso
seed = suppressWarnings(as.numeric(argv$seed))
if(!is.finite(seed))
  seed = NULL
n.cores = argv$thread

### other settings
data_sep = '\t'

### process arguments
implemented_methods = c('wgcna_pearson_unsigned_0.6', 'wgcna_pearson_unsigned_0.7', 'wgcna_pearson_unsigned_0.8', 'wgcna_pearson_unsigned_0.9', 
                        'wgcna_pearson_signed_0.6', 'wgcna_pearson_signed_0.7', 'wgcna_pearson_signed_0.8', 'wgcna_pearson_signed_0.9',
                        'wgcna_spearman_unsigned_0.6', 'wgcna_spearman_unsigned_0.7', 'wgcna_spearman_unsigned_0.8', 'wgcna_spearman_unsigned_0.9', 
                        'wgcna_spearman_signed_0.6', 'wgcna_spearman_signed_0.7', 'wgcna_spearman_signed_0.8', 'wgcna_spearman_signed_0.9',
                        
                        'wgcna',

                        'pearson', 'spearman', 
                        # 'pearson_rp', 'spearman_rp',
                        
                        'genie3', 'et_genie3',
                        
                        'aracne_spearman_e0_freq', 'aracne_spearman_e.05_freq', 'aracne_spearman_e.1_freq',
                        'aracne_spearman_e0_width', 'aracne_spearman_e.05_width', 'aracne_spearman_e.1_width',
                        
                        'aracne_pearson_e0_freq', 'aracne_pearson_e.05_freq', 'aracne_pearson_e.1_freq',
                        'aracne_pearson_e0_width', 'aracne_pearson_e.05_width', 'aracne_pearson_e.1_width',
                        
                        'aracne_empirical_e0_freq', 'aracne_empirical_e.05_freq', 'aracne_empirical_e.1_freq',
                        'aracne_empirical_e0_width', 'aracne_empirical_e.05_width', 'aracne_empirical_e.1_width',
                        
                        'aracne',
                        
                        'mrnetb_spearman_freq', 'mrnetb_pearson_freq', 'mrnetb_empirical_freq',
                        'mrnetb_spearman_width', 'mrnetb_pearson_width', 'mrnetb_empirical_width',
                        
                        'mrnetb',
                        
                        'mrnet_spearman_freq', 'mrnet_pearson_freq', 'mrnet_empirical_freq',
                        'mrnet_spearman_width', 'mrnet_pearson_width', 'mrnet_empirical_width',
                        
                        'mrnet',
                        
                        'clr_spearman_freq', 'clr_pearson_freq', 'clr_empirical_freq', 
                        'clr_spearman_width', 'clr_pearson_width', 'clr_empirical_width', 
                        
                        'clr',
                        
                        'glasso_likelihood', 
                        'glasso_5k', 'glasso_10k', 'glasso_25k', 'glasso_50k',
                        
                        'pcorr_pearson', 'pcorr_spearman', 
                        
                        'pcorr',
                        
                        'pcit', 'rlowpc', 'random',
                        
                        'spice_Mp_G1_Ra_A0', 'spice_Mp_G1_Rm_A0', 'spice_Mp_G1_RM_A0', 'spice_Mp_G1_Ra_A1', 'spice_Mp_G1_Rm_A1', 'spice_Mp_G1_RM_A1',
                        'spice_Mp_G.9_Ra_A0', 'spice_Mp_G.9_Rm_A0', 'spice_Mp_G.9_RM_A0', 'spice_Mp_G.9_Ra_A1', 'spice_Mp_G.9_Rm_A1', 'spice_Mp_G.9_RM_A1',
                        'spice_Mp_G.8_Ra_A0', 'spice_Mp_G.8_Rm_A0', 'spice_Mp_G.8_RM_A0', 'spice_Mp_G.8_Ra_A1', 'spice_Mp_G.8_Rm_A1', 'spice_Mp_G.8_RM_A1',
                        'spice_Mp_G.7_Ra_A0', 'spice_Mp_G.7_Rm_A0', 'spice_Mp_G.7_RM_A0', 'spice_Mp_G.7_Ra_A1', 'spice_Mp_G.7_Rm_A1', 'spice_Mp_G.7_RM_A1',
                        
                        'spice_Ms_G1_Ra_A0', 'spice_Ms_G1_Rm_A0', 'spice_Ms_G1_RM_A0', 'spice_Ms_G1_Ra_A1', 'spice_Ms_G1_Rm_A1', 'spice_Ms_G1_RM_A1',
                        'spice_Ms_G.9_Ra_A0', 'spice_Ms_G.9_Rm_A0', 'spice_Ms_G.9_RM_A0', 'spice_Ms_G.9_Ra_A1', 'spice_Ms_G.9_Rm_A1', 'spice_Ms_G.9_RM_A1',
                        'spice_Ms_G.8_Ra_A0', 'spice_Ms_G.8_Rm_A0', 'spice_Ms_G.8_RM_A0', 'spice_Ms_G.8_Ra_A1', 'spice_Ms_G.8_Rm_A1', 'spice_Ms_G.8_RM_A1',
                        'spice_Ms_G.7_Ra_A0', 'spice_Ms_G.7_Rm_A0', 'spice_Ms_G.7_RM_A0', 'spice_Ms_G.7_Ra_A1', 'spice_Ms_G.7_Rm_A1', 'spice_Ms_G.7_RM_A1',
                        
                        'spice_ME_G1_Ra_A0', 'spice_ME_G1_Rm_A0', 'spice_ME_G1_RM_A0', 'spice_ME_G1_Ra_A1', 'spice_ME_G1_Rm_A1', 'spice_ME_G1_RM_A1',
                        'spice_ME_G.9_Ra_A0', 'spice_ME_G.9_Rm_A0', 'spice_ME_G.9_RM_A0', 'spice_ME_G.9_Ra_A1', 'spice_ME_G.9_Rm_A1', 'spice_ME_G.9_RM_A1',
                        'spice_ME_G.8_Ra_A0', 'spice_ME_G.8_Rm_A0', 'spice_ME_G.8_RM_A0', 'spice_ME_G.8_Ra_A1', 'spice_ME_G.8_Rm_A1', 'spice_ME_G.8_RM_A1',
                        'spice_ME_G.7_Ra_A0', 'spice_ME_G.7_Rm_A0', 'spice_ME_G.7_RM_A0', 'spice_ME_G.7_Ra_A1', 'spice_ME_G.7_Rm_A1', 'spice_ME_G.7_RM_A1',
                        
                        'spice_MP_G1_Ra_A0', 'spice_MP_G1_Rm_A0', 'spice_MP_G1_RM_A0', 'spice_MP_G1_Ra_A1', 'spice_MP_G1_Rm_A1', 'spice_MP_G1_RM_A1',
                        'spice_MP_G.9_Ra_A0', 'spice_MP_G.9_Rm_A0', 'spice_MP_G.9_RM_A0', 'spice_MP_G.9_Ra_A1', 'spice_MP_G.9_Rm_A1', 'spice_MP_G.9_RM_A1',
                        'spice_MP_G.8_Ra_A0', 'spice_MP_G.8_Rm_A0', 'spice_MP_G.8_RM_A0', 'spice_MP_G.8_Ra_A1', 'spice_MP_G.8_Rm_A1', 'spice_MP_G.8_RM_A1',
                        'spice_MP_G.7_Ra_A0', 'spice_MP_G.7_Rm_A0', 'spice_MP_G.7_RM_A0', 'spice_MP_G.7_Ra_A1', 'spice_MP_G.7_Rm_A1', 'spice_MP_G.7_RM_A1',
                        
                        'spice' )
stopifnot(method %in% implemented_methods)


create_glass_dir <- function(){
  ### glasso optimization settings
  glasso_dir = paste0(out_dir, '/glasso')
  if(!dir.exists(glasso_dir))
    dir.create(glasso_dir)
  return(glasso_dir)
}

create_wgcna_dir <- function(){
  wgcna_dir = paste0(out_dir, '/wgcna')
  if(!dir.exists(wgcna_dir))
    dir.create(wgcna_dir)
  return(wgcna_dir)
}

### log i/o
verbose_print(sprintf("expr_fn: %s", expr_fn))
verbose_print(sprintf("method: %s", method))
verbose_print(sprintf("out_dir: %s", out_dir))

### read and process data
if(endsWith(x = expr_fn, suffix = '.feather')){
  expr_df = read_feather_df(fn = expr_fn, rownames.col = 1)
} else {
  expr_df = read_df(expr_fn, sep = data_sep)
}

if(any(dim(expr_df)<10))
  warning(sprintf('small number of features/samples in expression matrix: %s x %s', nrow(expr_df), ncol(expr_df)))

##### read cross-mappability data #####
get_crossmap_mat <- function(crossmap_fn, annot_fn){
  crossmap_mat = NULL
  if(file.exists(annot_fn) && file.exists(crossmap_fn)){
    annot_df = read_df(annot_fn, header = T, row.names = F)
    annot_df = annot_df[annot_df$gene_name %in% rownames(expr_df), , drop=F]
    
    crossmap_df = read_df(crossmap_fn, header = F, row.names = F)
    colnames(crossmap_df) = c('gene1', 'gene2', 'crossmap')
    crossmap_df = crossmap_df[crossmap_df$crossmap > crossmap_threshold, , drop = F]
    crossmap_df =  filter_crossmap_by_genes(crossmap_df = crossmap_df, incl.genes = ensembl_genes)
    
    crossmap_df = merge(crossmap_df, annot_df[,c('gene_id', 'gene_name')], by.x = 'gene1', by.y = "gene_id", all = F)
    crossmap_df = merge(crossmap_df, annot_df[,c('gene_id', 'gene_name')], by.x = 'gene2', by.y = "gene_id", all = F, suffixes = c('1', '2'))
    crossmap_mat = acast(crossmap_df, gene_name1 ~ gene_name2, fun.aggregate = mean, value.var = 'crossmap', fill = 0)
  }
  
  return(crossmap_mat)
}

######## network reconstruction methods ######
##### wgcna methods #####
# 'wgcna_pearson_unsigned_0.6', 'wgcna_pearson_unsigned_0.7', 'wgcna_pearson_unsigned_0.8', 'wgcna_pearson_unsigned_0.9', 
# 'wgcna_pearson_signed_0.6', 'wgcna_pearson_signed_0.7', 'wgcna_pearson_signed_0.8', 'wgcna_pearson_signed_0.9',
# 'wgcna_spearman_unsigned_0.6', 'wgcna_spearman_unsigned_0.7', 'wgcna_spearman_unsigned_0.8', 'wgcna_spearman_unsigned_0.9', 
# 'wgcna_spearman_signed_0.6', 'wgcna_spearman_signed_0.7', 'wgcna_spearman_signed_0.8', 'wgcna_spearman_signed_0.9',

run_wgcna_generic <- function(){
  wgcna.method = NULL
  if(grepl(pattern = "_spearman_", x = method, fixed = TRUE)){
    wgcna.method = "spearman"
  } else if(grepl(pattern = "_pearson_", x = method, fixed = TRUE)){
    wgcna.method = "pearson"
  } else {
    stop("wgcna method not implemented: check wgcna.method")
  }
  
  r2 = NULL
  if(grepl(pattern = "_0.6", x = method, fixed = TRUE)){
    r2 = 0.6
  } else if(grepl(pattern = "_0.7", x = method, fixed = TRUE)){
    r2 = 0.7
  } else if(grepl(pattern = "_0.8", x = method, fixed = TRUE)){
    r2 = 0.8
  } else if(grepl(pattern = "_0.9", x = method, fixed = TRUE)){
    r2 = 0.9
  } else {
    stop("wgcna method not implemented: check r2")
  }
  
  wgcna.type = NULL
  if(grepl(pattern = "_unsigned_", x = method, fixed = TRUE)){
    wgcna.type = "unsigned"
  } else if(grepl(pattern = "_signed_", x = method, fixed = TRUE)){
    wgcna.type = "signed"
  } else {
    stop("wgcna method not implemented: check wgcna.type")
  }
  
  net <- get_wgcna_net(expr_mat = expr_df, 
                       method = wgcna.method,
                       type = wgcna.type,
                       RsquaredCut = r2, 
                       verbose = T)
  return(net)
}

run_wgcna_pearson_unsigned_0.6 = run_wgcna_generic 
run_wgcna_pearson_unsigned_0.7 = run_wgcna_generic 
run_wgcna_pearson_unsigned_0.8 = run_wgcna_generic 
run_wgcna_pearson_unsigned_0.9 = run_wgcna_generic

run_wgcna_pearson_signed_0.6 = run_wgcna_generic 
run_wgcna_pearson_signed_0.7 = run_wgcna_generic 
run_wgcna_pearson_signed_0.8 = run_wgcna_generic 
run_wgcna_pearson_signed_0.9 = run_wgcna_generic

run_wgcna_spearman_unsigned_0.6 = run_wgcna_generic 
run_wgcna_spearman_unsigned_0.7 = run_wgcna_generic 
run_wgcna_spearman_unsigned_0.8 = run_wgcna_generic 
run_wgcna_spearman_unsigned_0.9 = run_wgcna_generic

run_wgcna_spearman_signed_0.6 = run_wgcna_generic 
run_wgcna_spearman_signed_0.7 = run_wgcna_generic 
run_wgcna_spearman_signed_0.8 = run_wgcna_generic 
run_wgcna_spearman_signed_0.9 = run_wgcna_generic

run_wgcna <- function(){
  wgcna.method = "pearson"
  r2 = 0.6
  wgcna.type = "signed"
  net <- get_wgcna_net(expr_mat = expr_df, 
                       method = wgcna.method,
                       type = wgcna.type,
                       RsquaredCut = r2, 
                       verbose = T)
  return(net)
}

##### correlation methods #####
run_pearson <- function(){
  net <- get_cor_net(expr_mat = expr_df, method = 'pearson')
  return(net)
}

run_spearman <- function(){
  net <- get_cor_net(expr_mat = expr_df, method = 'spearman')
  return(net)
}

# # 'pearson_rp', 'spearman_rp',

##### genie3 methods #####
# 'genie3', 'et_genie3',
run_genie3 <- function(){
  net <- get_genie3_net(expr_df, n.cores = n.cores, verbose = T, seed = seed)
  return(net)
}

run_et_genie3 <- function(){
  net <- get_genie3_net(expr_df, n.cores = n.cores, verbose = T, tree.method = 'ET', seed = seed)
  return(net)
}

##### aracne methods #####
# 'aracne_spearman_e0_freq', 'aracne_spearman_e.05_freq', 'aracne_spearman_e.1_freq',
# 'aracne_spearman_e0_width', 'aracne_spearman_e.05_width', 'aracne_spearman_e.1_width',
# 
# 'aracne_pearson_e0_freq', 'aracne_pearson_e.05_freq', 'aracne_pearson_e.1_freq',
# 'aracne_pearson_e0_width', 'aracne_pearson_e.05_width', 'aracne_pearson_e.1_width',
# 
# 'aracne_empirical_e0_freq', 'aracne_empirical_e.05_freq', 'aracne_empirical_e.1_freq',
# 'aracne_empirical_e0_width', 'aracne_empirical_e.05_width', 'aracne_empirical_e.1_width',

run_aracne_generic <- function(){
  mi.estimator = NULL
  if(grepl(pattern = "_spearman_", x = method, fixed = TRUE)){
    mi.estimator = "spearman"
  } else if(grepl(pattern = "_pearson_", x = method, fixed = TRUE)){
    mi.estimator = "pearson"
  } else if(grepl(pattern = "_empirical_", x = method, fixed = TRUE)){
    mi.estimator = "mi.empirical"
  } else {
    stop("arcne method not implemented: check mi.estimator")
  }
  
  eps = NULL
  if(grepl(pattern = "_e0_", x = method, fixed = TRUE)){
    eps = 0
  } else if(grepl(pattern = "_e.05_", x = method, fixed = TRUE)){
    eps = 0.05
  } else if(grepl(pattern = "_e.1_", x = method, fixed = TRUE)){
    eps = 0.1
  } else {
    stop("arcne method not implemented: check eps")
  }
  
  disc = NULL
  if(grepl(pattern = "_freq", x = method, fixed = TRUE)){
    disc = "equalfreq"
  } else if(grepl(pattern = "_width", x = method, fixed = TRUE)){
    disc = "equalwidth"
  } else {
    stop("arcne method not implemented: check disc")
  }
  
  net <- get_aracne_net(expr_df, mi.estimator = mi.estimator, eps = eps, disc = disc)
  return(net)
  
}

run_aracne_spearman_e0_freq= run_aracne_generic 
run_aracne_spearman_e.05_freq = run_aracne_generic 
run_aracne_spearman_e.1_freq = run_aracne_generic

run_aracne_spearman_e0_width = run_aracne_generic 
run_aracne_spearman_e.05_width = run_aracne_generic 
run_aracne_spearman_e.1_width = run_aracne_generic

run_aracne_pearson_e0_freq = run_aracne_generic 
run_aracne_pearson_e.05_freq = run_aracne_generic 
run_aracne_pearson_e.1_freq = run_aracne_generic

run_aracne_pearson_e0_width = run_aracne_generic 
run_aracne_pearson_e.05_width = run_aracne_generic 
run_aracne_pearson_e.1_width = run_aracne_generic

run_aracne_empirical_e0_freq = run_aracne_generic 
run_aracne_empirical_e.05_freq = run_aracne_generic 
run_aracne_empirical_e.1_freq = run_aracne_generic

run_aracne_empirical_e0_width = run_aracne_generic 
run_aracne_empirical_e.05_width = run_aracne_generic 
run_aracne_empirical_e.1_width = run_aracne_generic

run_aracne <- function(){
  mi.estimator = "pearson"
  eps = 0.1
  disc = "equalfreq"
  net <- get_aracne_net(expr_df, mi.estimator = mi.estimator, eps = eps, disc = disc)
  return(net)
}


##### mrnetb methods ######
# 'mrnetb_spearman_freq', 'mrnetb_pearson_freq', 'mrnetb_empirical_freq',
# 'mrnetb_spearman_width', 'mrnetb_pearson_width', 'mrnetb_empirical_width',

run_mrnetb_generic <- function(){
  mi.estimator = NULL
  if(grepl(pattern = "_spearman_", x = method, fixed = TRUE)){
    mi.estimator = "spearman"
  } else if(grepl(pattern = "_pearson_", x = method, fixed = TRUE)){
    mi.estimator = "pearson"
  } else if(grepl(pattern = "_empirical_", x = method, fixed = TRUE)){
    mi.estimator = "mi.empirical"
  } else {
    stop("mrnetb method not implemented: check mi.estimator")
  }
  
  disc = NULL
  if(grepl(pattern = "_freq", x = method, fixed = TRUE)){
    disc = "equalfreq"
  } else if(grepl(pattern = "_width", x = method, fixed = TRUE)){
    disc = "equalwidth"
  } else {
    stop("mrnetb method not implemented: check disc")
  }
  
  net <- get_mrnetb_net(expr_df, mi.estimator = mi.estimator, disc = disc)
  return(net)
}

run_mrnetb_spearman_freq = run_mrnetb_generic 
run_mrnetb_pearson_freq = run_mrnetb_generic 
run_mrnetb_empirical_freq = run_mrnetb_generic
run_mrnetb_spearman_width = run_mrnetb_generic 
run_mrnetb_pearson_width = run_mrnetb_generic 
run_mrnetb_empirical_width = run_mrnetb_generic

run_mrnetb <- function(){
  mi.estimator = "pearson"
  disc = "equalfreq"
  net <- get_mrnetb_net(expr_df, mi.estimator = mi.estimator, disc = disc)
  return(net)
}

##### mrnet methods #####
# 'mrnet_spearman_freq', 'mrnet_pearson_freq', 'mrnet_empirical_freq',
# 'mrnet_spearman_width', 'mrnet_pearson_width', 'mrnet_empirical_width',

run_mrnet_generic <- function(){
  mi.estimator = NULL
  if(grepl(pattern = "_spearman_", x = method, fixed = TRUE)){
    mi.estimator = "spearman"
  } else if(grepl(pattern = "_pearson_", x = method, fixed = TRUE)){
    mi.estimator = "pearson"
  } else if(grepl(pattern = "_empirical_", x = method, fixed = TRUE)){
    mi.estimator = "mi.empirical"
  } else {
    stop("mrnet method not implemented: check mi.estimator")
  }
  
  disc = NULL
  if(grepl(pattern = "_freq", x = method, fixed = TRUE)){
    disc = "equalfreq"
  } else if(grepl(pattern = "_width", x = method, fixed = TRUE)){
    disc = "equalwidth"
  } else {
    stop("mrnet method not implemented: check disc")
  }
  
  net <- get_mrnet_net(expr_df, mi.estimator = mi.estimator, disc = disc)
  return(net)
}

run_mrnet_spearman_freq = run_mrnet_generic 
run_mrnet_pearson_freq = run_mrnet_generic 
run_mrnet_empirical_freq = run_mrnet_generic
run_mrnet_spearman_width = run_mrnet_generic 
run_mrnet_pearson_width = run_mrnet_generic 
run_mrnet_empirical_width = run_mrnet_generic

run_mrnet <- function(){
  mi.estimator = "pearson"
  disc = "equalfreq"
  net <- get_mrnet_net(expr_df, mi.estimator = mi.estimator, disc = disc)
  return(net)
}

##### clr methods #######
# 'clr_spearman_freq', 'clr_pearson_freq', 'clr_empirical_freq', 
# 'clr_spearman_width', 'clr_pearson_width', 'clr_empirical_width', 
run_clr_generic <- function(){
  mi.estimator = NULL
  if(grepl(pattern = "_spearman_", x = method, fixed = TRUE)){
    mi.estimator = "spearman"
  } else if(grepl(pattern = "_pearson_", x = method, fixed = TRUE)){
    mi.estimator = "pearson"
  } else if(grepl(pattern = "_empirical_", x = method, fixed = TRUE)){
    mi.estimator = "mi.empirical"
  } else {
    stop("mrnet method not implemented: check mi.estimator")
  }
  
  disc = NULL
  if(grepl(pattern = "_freq", x = method, fixed = TRUE)){
    disc = "equalfreq"
  } else if(grepl(pattern = "_width", x = method, fixed = TRUE)){
    disc = "equalwidth"
  } else {
    stop("mrnet method not implemented: check disc")
  }
  
  net <- get_clr_net(expr_df, mi.estimator = mi.estimator, disc = disc)
  return(net)
}

run_clr_spearman_freq = run_clr_generic 
run_clr_pearson_freq = run_clr_generic 
run_clr_empirical_freq = run_clr_generic
run_clr_spearman_width = run_clr_generic 
run_clr_pearson_width = run_clr_generic 
run_clr_empirical_width = run_clr_generic

run_clr <- function(){
  mi.estimator = "pearson"
  disc = "equalfreq"
  net <- get_clr_net(expr_df, mi.estimator = mi.estimator, disc = disc)
  return(net)
}

##### glasso methods #####
# 'glasso_likelihood', 
run_glasso_likelihood <- function(){
  lambdas = seq(0.05, 0.5, 0.05)
  fold = 1
  glasso_dir = create_glass_dir()
  out.pfx = sprintf("%s/glasso", glasso_dir)
  validation.frac = 0.2
  glasso_seed = ifelse(is.numeric(seed), seed, 101) # must be numeric
  net <- get_glasso_net_optimized_by_likelihood(expr_df = expr_df, 
                                                lambdas = lambdas, 
                                                out.pfx = out.pfx, 
                                                fold = fold, 
                                                validation.frac = validation.frac,
                                                save.net = T,
                                                use.old.net = use_old_glasso,
                                                verbose = T, 
                                                seed = glasso_seed)
  return(net)
}

# 'glasso_5k', 'glasso_10k', 'glasso_25k', 'glasso_50k',

run_glasso_net_optimized_by_n_edge_generic <- function(){
  lambda.max = 0.5
  lambda.min = 0.05
  lambda.interval = 0.05
  lambda.fine.interval = 0.01
  
  expected.n.edge = NULL
  if(grepl(pattern = "_5k", x = method, fixed = TRUE)){
    expected.n.edge = 5000
  } else if(grepl(pattern = "_10k", x = method, fixed = TRUE)){
    expected.n.edge = 10000
  } else if(grepl(pattern = "_25k", x = method, fixed = TRUE)){
    expected.n.edge = 25000
  } else if(grepl(pattern = "_50k", x = method, fixed = TRUE)){
    expected.n.edge = 50000
  } else {
    stop("mrnet method not implemented: check mi.estimator")
  }
  
  glasso_dir = create_glass_dir()
  out.pfx = sprintf("%s/glasso", glasso_dir)
  
  net <- get_glasso_net_optimized_by_edge_count(expr_df = expr_df, 
                                                lambda.max =   lambda.max,
                                                lambda.min = lambda.min, 
                                                lambda.interval = lambda.interval, 
                                                lambda.fine.interval = lambda.fine.interval, 
                                                expected.n.edge = expected.n.edge,
                                                out.pfx = out.pfx, 
                                                inverse.normal = T,
                                                save.net = T,
                                                use.old.net = use_old_glasso,
                                                verbose = T)
  return(net)
}

run_glasso_5k = run_glasso_net_optimized_by_n_edge_generic 
run_glasso_10k = run_glasso_net_optimized_by_n_edge_generic 
run_glasso_25k = run_glasso_net_optimized_by_n_edge_generic 
run_glasso_50k = run_glasso_net_optimized_by_n_edge_generic

##### partial correlation ######
# 'pcorr_pearson', 'pcorr_spearman', 
run_pcorr_pearson <- function(){
  net <- get_pcorr_net(expr_df, method = 'pearson')
  return(net)
}

run_pcorr_spearman <- function(){
  net <- get_pcorr_net(expr_df, method = 'spearman')
  return(net)
}

run_pcorr <- function(){
  net <- get_pcorr_net(expr_df, method = 'spearman')
  return(net)
}

##### pcit, rlowpc, random #########
# 'pcit', 'rlowpc', 'random',
run_pcit <- function(){
  net <- get_pcit_net(expr_df, method = "pearson")
  return(net)
}

run_rlowpc <- function(){
  net <- get_rlowpc_net(expr_df)
  return(net)
}

run_random <- function(){
  net <- get_random_net(expr_df, seed = seed)
  return(net)
}

##### spice methods ####################
# 'spice_Mp_G1_Ra_A0', 'spice_Mp_G1_Rm_A0', 'spice_Mp_G1_RM_A0', 'spice_Mp_G1_Ra_A1', 'spice_Mp_G1_Rm_A1', 'spice_Mp_G1_RM_A1',
# 'spice_Mp_G.9_Ra_A0', 'spice_Mp_G.9_Rm_A0', 'spice_Mp_G.9_RM_A0', 'spice_Mp_G.9_Ra_A1', 'spice_Mp_G.9_Rm_A1', 'spice_Mp_G.9_RM_A1',
# 'spice_Mp_G.8_Ra_A0', 'spice_Mp_G.8_Rm_A0', 'spice_Mp_G.8_RM_A0', 'spice_Mp_G.8_Ra_A1', 'spice_Mp_G.8_Rm_A1', 'spice_Mp_G.8_RM_A1',
# 'spice_Mp_G.7_Ra_A0', 'spice_Mp_G.7_Rm_A0', 'spice_Mp_G.7_RM_A0', 'spice_Mp_G.7_Ra_A1', 'spice_Mp_G.7_Rm_A1', 'spice_Mp_G.7_RM_A1',
# 
# 'spice_Ms_G1_Ra_A0', 'spice_Ms_G1_Rm_A0', 'spice_Ms_G1_RM_A0', 'spice_Ms_G1_Ra_A1', 'spice_Ms_G1_Rm_A1', 'spice_Ms_G1_RM_A1',
# 'spice_Ms_G.9_Ra_A0', 'spice_Ms_G.9_Rm_A0', 'spice_Ms_G.9_RM_A0', 'spice_Ms_G.9_Ra_A1', 'spice_Ms_G.9_Rm_A1', 'spice_Ms_G.9_RM_A1',
# 'spice_Ms_G.8_Ra_A0', 'spice_Ms_G.8_Rm_A0', 'spice_Ms_G.8_RM_A0', 'spice_Ms_G.8_Ra_A1', 'spice_Ms_G.8_Rm_A1', 'spice_Ms_G.8_RM_A1',
# 'spice_Ms_G.7_Ra_A0', 'spice_Ms_G.7_Rm_A0', 'spice_Ms_G.7_RM_A0', 'spice_Ms_G.7_Ra_A1', 'spice_Ms_G.7_Rm_A1', 'spice_Ms_G.7_RM_A1',
# 
# 'spice_ME_G1_Ra_A0', 'spice_ME_G1_Rm_A0', 'spice_ME_G1_RM_A0', 'spice_ME_G1_Ra_A1', 'spice_ME_G1_Rm_A1', 'spice_ME_G1_RM_A1',
# 'spice_ME_G.9_Ra_A0', 'spice_ME_G.9_Rm_A0', 'spice_ME_G.9_RM_A0', 'spice_ME_G.9_Ra_A1', 'spice_ME_G.9_Rm_A1', 'spice_ME_G.9_RM_A1',
# 'spice_ME_G.8_Ra_A0', 'spice_ME_G.8_Rm_A0', 'spice_ME_G.8_RM_A0', 'spice_ME_G.8_Ra_A1', 'spice_ME_G.8_Rm_A1', 'spice_ME_G.8_RM_A1',
# 'spice_ME_G.7_Ra_A0', 'spice_ME_G.7_Rm_A0', 'spice_ME_G.7_RM_A0', 'spice_ME_G.7_Ra_A1', 'spice_ME_G.7_Rm_A1', 'spice_ME_G.7_RM_A1',
# 
# 'spice_MP_G1_Ra_A0', 'spice_MP_G1_Rm_A0', 'spice_MP_G1_RM_A0', 'spice_MP_G1_Ra_A1', 'spice_MP_G1_Rm_A1', 'spice_MP_G1_RM_A1',
# 'spice_MP_G.9_Ra_A0', 'spice_MP_G.9_Rm_A0', 'spice_MP_G.9_RM_A0', 'spice_MP_G.9_Ra_A1', 'spice_MP_G.9_Rm_A1', 'spice_MP_G.9_RM_A1',
# 'spice_MP_G.8_Ra_A0', 'spice_MP_G.8_Rm_A0', 'spice_MP_G.8_RM_A0', 'spice_MP_G.8_Ra_A1', 'spice_MP_G.8_Rm_A1', 'spice_MP_G.8_RM_A1',
# 'spice_MP_G.7_Ra_A0', 'spice_MP_G.7_Rm_A0', 'spice_MP_G.7_RM_A0', 'spice_MP_G.7_Ra_A1', 'spice_MP_G.7_Rm_A1', 'spice_MP_G.7_RM_A1',

run_spice_generic <- function(){
  spice_method = NULL
  if(grepl(pattern = "_Mp_", x = method, fixed = TRUE, ignore.case = FALSE)){
    spice_method = "pearson"
  } else if(grepl(pattern = "_Ms_", x = method, fixed = TRUE, ignore.case = FALSE)){
    spice_method = "spearman"
  } else if(grepl(pattern = "_ME_", x = method, fixed = TRUE, ignore.case = FALSE)){
    spice_method = "mi.empirical"
  } else if(grepl(pattern = "_MP_", x = method, fixed = TRUE, ignore.case = FALSE)){
    spice_method = "mi.pearson"
  } else {
    stop("spice method not implemented: check spice_method")
  }
  
  
  frac.gene = NULL
  if(grepl(pattern = "_G1_", x = method, fixed = TRUE, ignore.case = FALSE)){
    frac.gene = 1.0
  } else if(grepl(pattern = "_G.9_", x = method, fixed = TRUE, ignore.case = FALSE)){
    frac.gene = 0.9
  } else if(grepl(pattern = "_G.8_", x = method, fixed = TRUE, ignore.case = FALSE)){
    frac.gene = 0.8
  } else if(grepl(pattern = "_G.7_", x = method, fixed = TRUE, ignore.case = FALSE)){
    frac.gene = 0.7
  } else {
    stop("spice method not implemented: check frac.gene")
  }
  
  rank.ties = NULL
  if(grepl(pattern = "_Ra_", x = method, fixed = TRUE, ignore.case = FALSE)){
    rank.ties = "average"
  } else if(grepl(pattern = "_Rm_", x = method, fixed = TRUE, ignore.case = FALSE)){
    rank.ties = "min"
  } else if(grepl(pattern = "_RM_", x = method, fixed = TRUE, ignore.case = FALSE)){
    rank.ties = "max"
  } else {
    stop("spice method not implemented: check rank.ties")
  }
  
  adjust.weight = NULL
  if(grepl(pattern = "_A0", x = method, fixed = TRUE, ignore.case = FALSE)){
    adjust.weight = F
  } else if(grepl(pattern = "_A1", x = method, fixed = TRUE, ignore.case = FALSE)){
    adjust.weight = T
  } else {
    stop("spice method not implemented: check adjust.weight")
  }
  
  net <- spice(expr = expr_df, 
               method = spice_method, 
               iter = 100, 
               frac.gene = frac.gene, 
               frac.sample = 0.8, 
               n.cores = n.cores, 
               rank.ties = rank.ties, 
               filter.mat = NULL, 
               weight.method = "qnorm", 
               adjust.weight = adjust.weight, 
               adjust.clr = F, 
               verbose = T, 
               seed = seed)
  
  return(net)
}

run_spice_Mp_G1_Ra_A0 = run_spice_generic 
run_spice_Mp_G1_Rm_A0 = run_spice_generic 
run_spice_Mp_G1_RM_A0 = run_spice_generic 
run_spice_Mp_G1_Ra_A1 = run_spice_generic 
run_spice_Mp_G1_Rm_A1 = run_spice_generic 
run_spice_Mp_G1_RM_A1 = run_spice_generic
run_spice_Mp_G.9_Ra_A0 = run_spice_generic 
run_spice_Mp_G.9_Rm_A0 = run_spice_generic 
run_spice_Mp_G.9_RM_A0 = run_spice_generic 
run_spice_Mp_G.9_Ra_A1 = run_spice_generic 
run_spice_Mp_G.9_Rm_A1 = run_spice_generic 
run_spice_Mp_G.9_RM_A1 = run_spice_generic
run_spice_Mp_G.8_Ra_A0 = run_spice_generic 
run_spice_Mp_G.8_Rm_A0 = run_spice_generic 
run_spice_Mp_G.8_RM_A0 = run_spice_generic 
run_spice_Mp_G.8_Ra_A1 = run_spice_generic 
run_spice_Mp_G.8_Rm_A1 = run_spice_generic 
run_spice_Mp_G.8_RM_A1 = run_spice_generic
run_spice_Mp_G.7_Ra_A0 = run_spice_generic 
run_spice_Mp_G.7_Rm_A0 = run_spice_generic 
run_spice_Mp_G.7_RM_A0 = run_spice_generic 
run_spice_Mp_G.7_Ra_A1 = run_spice_generic 
run_spice_Mp_G.7_Rm_A1 = run_spice_generic 
run_spice_Mp_G.7_RM_A1 = run_spice_generic

run_spice_Ms_G1_Ra_A0 = run_spice_generic 
run_spice_Ms_G1_Rm_A0 = run_spice_generic 
run_spice_Ms_G1_RM_A0 = run_spice_generic 
run_spice_Ms_G1_Ra_A1 = run_spice_generic 
run_spice_Ms_G1_Rm_A1 = run_spice_generic 
run_spice_Ms_G1_RM_A1 = run_spice_generic
run_spice_Ms_G.9_Ra_A0 = run_spice_generic 
run_spice_Ms_G.9_Rm_A0 = run_spice_generic 
run_spice_Ms_G.9_RM_A0 = run_spice_generic 
run_spice_Ms_G.9_Ra_A1 = run_spice_generic 
run_spice_Ms_G.9_Rm_A1 = run_spice_generic 
run_spice_Ms_G.9_RM_A1 = run_spice_generic
run_spice_Ms_G.8_Ra_A0 = run_spice_generic 
run_spice_Ms_G.8_Rm_A0 = run_spice_generic 
run_spice_Ms_G.8_RM_A0 = run_spice_generic 
run_spice_Ms_G.8_Ra_A1 = run_spice_generic 
run_spice_Ms_G.8_Rm_A1 = run_spice_generic 
run_spice_Ms_G.8_RM_A1 = run_spice_generic
run_spice_Ms_G.7_Ra_A0 = run_spice_generic 
run_spice_Ms_G.7_Rm_A0 = run_spice_generic 
run_spice_Ms_G.7_RM_A0 = run_spice_generic 
run_spice_Ms_G.7_Ra_A1 = run_spice_generic 
run_spice_Ms_G.7_Rm_A1 = run_spice_generic 
run_spice_Ms_G.7_RM_A1 = run_spice_generic

run_spice_ME_G1_Ra_A0 = run_spice_generic 
run_spice_ME_G1_Rm_A0 = run_spice_generic 
run_spice_ME_G1_RM_A0 = run_spice_generic 
run_spice_ME_G1_Ra_A1 = run_spice_generic 
run_spice_ME_G1_Rm_A1 = run_spice_generic 
run_spice_ME_G1_RM_A1 = run_spice_generic
run_spice_ME_G.9_Ra_A0 = run_spice_generic 
run_spice_ME_G.9_Rm_A0 = run_spice_generic 
run_spice_ME_G.9_RM_A0 = run_spice_generic 
run_spice_ME_G.9_Ra_A1 = run_spice_generic 
run_spice_ME_G.9_Rm_A1 = run_spice_generic 
run_spice_ME_G.9_RM_A1 = run_spice_generic
run_spice_ME_G.8_Ra_A0 = run_spice_generic 
run_spice_ME_G.8_Rm_A0 = run_spice_generic 
run_spice_ME_G.8_RM_A0 = run_spice_generic 
run_spice_ME_G.8_Ra_A1 = run_spice_generic 
run_spice_ME_G.8_Rm_A1 = run_spice_generic 
run_spice_ME_G.8_RM_A1 = run_spice_generic
run_spice_ME_G.7_Ra_A0 = run_spice_generic 
run_spice_ME_G.7_Rm_A0 = run_spice_generic 
run_spice_ME_G.7_RM_A0 = run_spice_generic 
run_spice_ME_G.7_Ra_A1 = run_spice_generic 
run_spice_ME_G.7_Rm_A1 = run_spice_generic 
run_spice_ME_G.7_RM_A1 = run_spice_generic

run_spice_MP_G1_Ra_A0 = run_spice_generic 
run_spice_MP_G1_Rm_A0 = run_spice_generic 
run_spice_MP_G1_RM_A0 = run_spice_generic 
run_spice_MP_G1_Ra_A1 = run_spice_generic 
run_spice_MP_G1_Rm_A1 = run_spice_generic 
run_spice_MP_G1_RM_A1 = run_spice_generic
run_spice_MP_G.9_Ra_A0 = run_spice_generic 
run_spice_MP_G.9_Rm_A0 = run_spice_generic 
run_spice_MP_G.9_RM_A0 = run_spice_generic 
run_spice_MP_G.9_Ra_A1 = run_spice_generic 
run_spice_MP_G.9_Rm_A1 = run_spice_generic 
run_spice_MP_G.9_RM_A1 = run_spice_generic
run_spice_MP_G.8_Ra_A0 = run_spice_generic 
run_spice_MP_G.8_Rm_A0 = run_spice_generic 
run_spice_MP_G.8_RM_A0 = run_spice_generic 
run_spice_MP_G.8_Ra_A1 = run_spice_generic 
run_spice_MP_G.8_Rm_A1 = run_spice_generic 
run_spice_MP_G.8_RM_A1 = run_spice_generic
run_spice_MP_G.7_Ra_A0 = run_spice_generic 
run_spice_MP_G.7_Rm_A0 = run_spice_generic 
run_spice_MP_G.7_RM_A0 = run_spice_generic 
run_spice_MP_G.7_Ra_A1 = run_spice_generic 
run_spice_MP_G.7_Rm_A1 = run_spice_generic 
run_spice_MP_G.7_RM_A1 = run_spice_generic

run_spice_generic <- function(){
  spice_method = "pearson"
  frac.gene = 0.8
  rank.ties = "average"
  adjust.weight = T
  net <- spice(expr = expr_df, 
               method = spice_method, 
               iter = 100, 
               frac.gene = frac.gene, 
               frac.sample = 0.8, 
               n.cores = n.cores, 
               rank.ties = rank.ties, 
               filter.mat = NULL, 
               weight.method = "qnorm", 
               adjust.weight = adjust.weight, 
               adjust.clr = F, 
               verbose = T, 
               seed = seed)
  
  return(net)
}

##### run appropriate method expression data #####
function_name = paste0("run_", method)
coexpression_functions = mget(function_name)
if(length(coexpression_functions) != 1){
  stop(sprintf("could not find an appropriate function for the method: %s.", method))
}
coexpression_function = coexpression_functions[[1]]

net_time = system.time(net <- coexpression_function())

# take absolute values and make symmetric
absnet <- abs(net)  # absolute value
absnet <- pmax(absnet, t(absnet), na.rm = T) # symmetric

##### save network and time #####
net_data_fn = sprintf('%s/%s_network.rds', out_dir, method)
saveRDS(net, file = net_data_fn)
absnet_data_fn = sprintf('%s/%s_absnet.rds', out_dir, method)
saveRDS(absnet, file = absnet_data_fn)
time_data_fn = sprintf('%s/%s_time.rds', out_dir, method)
saveRDS(net_time, file = time_data_fn)
