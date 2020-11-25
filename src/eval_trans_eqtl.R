library(ioutil)
library(miscutil)
library(MatrixEQTL)
library(argparser)
library(mappabilityutil)
library(reshape2)
library(flock)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "--net", help="network prefix (before method)", default="/work-zfs/abattle4/ashis/progres/spice_anlysis/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/5000/pearson_absnet.rds")
args <- add_argument(args, "--top", help="number of top edges to evaluate tran-eqtls", default=10000)
args <- add_argument(args, "--expr", help="expression file", default="/work-zfs/abattle4/ashis/progres/spice_anlysis/gtex_v8/data/normalized/Muscle_Skeletal.txt")
args <- add_argument(args, "--cov", help="covariate file", default="/work-zfs/abattle4/ashis/progres/spice_anlysis/gtex_v8/data/GTEx_Analysis_v8_eQTL_covariates/Muscle_Skeletal.v8.covariates.txt")
args <- add_argument(args, "--geno_pfx", help="genotype file name prefix with path (before chromosome in the name)", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_")
args <- add_argument(args, "--geno_sfx", help="genotype file name suffix (after chromosome in the name)", default="_dosage_MAF_05_not_in_repeat.RData")
args <- add_argument(args, "--gene_snp", help="file with gene-to-snps map", default="/work-zfs/abattle4/ashis/prog/spice_analysis/results/snps_nearby_gene.rds")
args <- add_argument(args, "--annot", help="gene annotation file (txt)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt")
args <- add_argument(args, "--crossmap", help="crossmapping file", default='/work-zfs/abattle4/lab_data/annotation/mappability_hg38_gencode26/hg38_cross_mappability_strength.txt')
args <- add_argument(args, "--d", help="distance threshold for cross-mappability", default=1e6)
args <- add_argument(args, "--repo", help="repository of all tests statistics per tissue (*.rds)", default="results/repo.rds")
args <- add_argument(args, "--max_snp", help="maximum number of snps in a matrix-eqtl call", default=10000)
args <- add_argument(args, "--o", help="output file (*.rds)", default="results/n_egenes.rds")

argv = parse_args(args)
net_fn = argv$net
n_top_edges = argv$top
expr_fn = argv$expr
cov_fn = argv$cov
geno_pfx = argv$geno_pfx
geno_sfx = argv$geno_sfx
gene_to_snps_fn = argv$gene_snp
gene_annot_fn = argv$annot
crossmap_fn = argv$crossmap
snp_gene_dist_for_crossmap = argv$d
repository_fn = argv$repo
max_snp_per_meqtl = argv$max_snp
out_fn = argv$o

repository_lock_fn = sprintf("%s.lock", repository_fn)

### check inputs
stopifnot(file.exists(net_fn))
stopifnot(endsWith(tolower(repository_fn), ".rds"))

### read network
net_mat = readRDS(net_fn)

### read expression
expr_df = read_df(expr_fn, header = T, row.names = T)

### read covariate
cov_df = read_df(cov_fn, header = T, row.names = T)
stopifnot(length(intersect(colnames(expr_df), colnames(cov_df))) == ncol(expr_df))
cov_df = cov_df[,colnames(expr_df), drop = F]
n_uniq = sapply(rownames(cov_df), function(idx)length(unique(as.numeric(cov_df[idx,]))))
if(any(n_uniq<=1)){
  const_covariates = names(n_uniq[n_uniq<=1])
  cov_df = cov_df[-which(rownames(cov_df) %in% const_covariates), , drop=F]
}

### read gene annotation and filter for genes in network and expression matrix
annot_df = read_df(gene_annot_fn, header = T, row.names = F)
rownames(annot_df) = annot_df$gene_id
annot_df = annot_df[annot_df$gene_id %in% rownames(expr_df), , drop = F]
annot_df = annot_df[annot_df$gene_name %in% rownames(net_mat), , drop = F]
sym_2_ensg = new.env(hash = T)
tmp <- mapply(function(g,ensg){
  sym_2_ensg[[g]] <<- ensg
  return(NA)
}, annot_df$gene_name, annot_df$gene_id)

### read gene-to-snp mapping
gene_2_snps = readRDS(gene_to_snps_fn)

### read cross-mappability data
crossmap_df = read_df(crossmap_fn, header = F, row.names = F)

### function to get top edges from a network 
### after filtering edges between cross-mappable genes
### and genes from the same chromosome
### return a data-frame with top n edges (by weights)
get_top_edges <- function(weight_mat, crossmap_df, annot_df, n.max.edges = NULL){
  ### ensure symmetry
  stopifnot(setequal(rownames(weight_mat), colnames(weight_mat)))
  genes = rownames(weight_mat)
  weight_mat = weight_mat[genes, genes]
  weight_mat = pmax(weight_mat, t(weight_mat), na.rm = T)
  weight_mat[!is.finite(weight_mat)] = NA
  
  ### map gene names to ensembl ids
  cur_annot_df = annot_df[annot_df$gene_name %in% genes, , drop = F]
  stopifnot(setequal(cur_annot_df$gene_name, genes) && nrow(cur_annot_df) == length(genes))
  # names(cur_annot_df) = cur_annot_df$gene_id
  ensgids = cur_annot_df$gene_id
  names(ensgids) = cur_annot_df$gene_name
  ensgids = ensgids[genes]
  
  ### remove cross-mappable gene pairs
  cur_crossmap_df = mappabilityutil::filter_crossmap_by_genes(crossmap_df, incl.genes = as.character(ensgids))
  colnames(cur_crossmap_df) = c('gene1', 'gene2', 'crossmap')
  cur_crossmap_df$sym1 = cur_annot_df[cur_crossmap_df[,1], 'gene_name']
  cur_crossmap_df$sym2 = cur_annot_df[cur_crossmap_df[,2], 'gene_name']
  cur_crossmap_mat0 = suppressWarnings(
    acast(
      data = cur_crossmap_df,
      formula = sym1 ~ sym2,
      value.var = "crossmap",
      fun.aggregate = max,
      na.rm = T
    )
  )
  cur_crossmap_mat = matrix(
    NA,
    nrow = nrow(weight_mat),
    ncol = ncol(weight_mat),
    dimnames = list(genes, genes)
  )
  cur_crossmap_mat[rownames(cur_crossmap_mat0), colnames(cur_crossmap_mat0)] = cur_crossmap_mat0
  cur_crossmap_mat = pmax(cur_crossmap_mat, t(cur_crossmap_mat), na.rm = T)
  weight_mat[cur_crossmap_mat>0] = NA
  
  ### remove gene pairs from same chr
  genes_chromosomes = cur_annot_df[as.character(ensgids), "chr"]
  chromosomes = unique(genes_chromosomes)
  for(chr in chromosomes){
    gidx = which(genes_chromosomes == chr)
    weight_mat[gidx, gidx] = NA
  }
  
  ### take top N edges
  weight_mat[lower.tri(weight_mat, diag = T)] = NA
  weights = as.numeric(weight_mat[!is.na(weight_mat)])
  q = ifelse(is.null(n.max.edges), 0, 1-n.max.edges/length(weights))
  edge.threshold = as.numeric(quantile(weights, probs = q))
  weight_mat[weight_mat < edge.threshold] = NA
  
  edges_df <- reshape2::melt(weight_mat, na.rm = TRUE)
  colnames(edges_df) <- c("gene1", "gene2", "weight")
  edges_df$gene1 = as.character(edges_df$gene1)
  edges_df$gene2 = as.character(edges_df$gene2)
  edges_df <- edges_df[order(edges_df$weight, decreasing = TRUE), ,drop = F]
  if (!is.null(n.max.edges)) {
    edges_df <- edges_df[seq_len(min(nrow(edges_df), n.max.edges)), , drop = F]
  }
  rownames(edges_df) <- NULL
  return(edges_df)
}

### function to get snp-gene pairs corresponding to given network edges
### TODO: remove cross-mappable snp-gene pairs
get_trans_snp_gene_pairs_to_test <- function(edges_df){
  # map snps nearby one gene to the other gene (both both genes)
  # remember to send unique snp-gene pairs only
  snp_gene_pair_list = mapply(function(gene1, gene2){
    snps1 = gene_2_snps[[gene1]]
    snps2 = gene_2_snps[[gene2]]
    pairs = NULL
    if(length(snps1)>0 || length(snps2) > 0){
      pairs = data.frame(
        snp = c(snps1, snps2),
        gene = c(rep(gene2, length(snps1)), rep(gene1, length(snps2))),
        stringsAsFactors = F
      )
    }
    return(pairs)
  }, edges_df[,1], edges_df[,2], SIMPLIFY = F)
  snp_gene_pairs = do.call(rbind, snp_gene_pair_list)
  rm(snp_gene_pair_list)
  gc()
  snp_gene_pairs = unique(snp_gene_pairs)
  rownames(snp_gene_pairs) = NULL
  return(snp_gene_pairs)
}

### function to read genotype file
read_genotype_file <- function(geno_fn, na_genotype_value="-"){
  stopifnot(endsWith(geno_fn, ".RData"))
  load(geno_fn)
  colnames(genotype_mat_not_in_repeat) = gsub('\\.', '-', colnames(genotype_mat_not_in_repeat))
  genotype_mat_not_in_repeat[genotype_mat_not_in_repeat==na_genotype_value] = NA
  for(cn in colnames(genotype_mat_not_in_repeat))
    genotype_mat_not_in_repeat[,cn] = as.numeric(genotype_mat_not_in_repeat[,cn])
  return(genotype_mat_not_in_repeat)
}

### function to get rows of x that are not present in y
### similar to setdiff() where each row is an element of the set.
setdiff_df <- function(x, y){
  stopifnot(is.data.frame(x) && is.data.frame(y))
  if (ncol(x) != ncol(y) || any(colnames(x) != colnames(y))) {
    stop("x and y must have set of columns in the same order.")
  }
  z = rbind(y, x)
  is.uniq = !duplicated(z)
  x.uniq = is.uniq[seq_len(nrow(x)) + nrow(y)]
  z = x[x.uniq, , drop = F]
  return(z)
}

### gather snp-gene pairs to test
verbose_print("extracting top edges ...")
top_edges = get_top_edges(
  weight_mat = net_mat,
  crossmap_df = crossmap_df,
  annot_df = annot_df,
  n.max.edges = n_top_edges
)
verbose_print("generating snp-gene pairs ...")
all_snp_gene_pairs = get_trans_snp_gene_pairs_to_test(top_edges)
snp_gene_pairs_to_test = all_snp_gene_pairs
rm(net_mat, top_edges)
gc()

### remove already-tested snp-gene pairs
if(file.exists(repository_fn)){
  fl = flock::lock(repository_lock_fn)
  tested_stats = readRDS(repository_fn)
  flock::unlock(fl)
  
  tested_pairs = tested_stats[, c('snps', 'gene_name'), drop = F]
  names(tested_pairs) = c('snp', 'gene')
  snp_gene_pairs_to_test = setdiff_df(snp_gene_pairs_to_test, tested_pairs)
  
  rm(tested_stats, tested_pairs)
  gc()
}

### add chr, pos and ensembl id
snp_gene_pairs_to_test = unique(snp_gene_pairs_to_test)
snp_parts = strsplit(snp_gene_pairs_to_test$snp, split = "_")
snps_chr = sapply(snp_parts, function(x) x[1])
snps_pos = sapply(snp_parts, function(x) as.integer(x[2]))
snp_gene_pairs_to_test$snp_chr = as.character(snps_chr) # chr of snp
snp_gene_pairs_to_test$snp_pos = as.integer(snps_pos)   # pos of snp
gene_ensgids = sapply(snp_gene_pairs_to_test$gene, function(g) sym_2_ensg[[g]])
snp_gene_pairs_to_test$ensgid = gene_ensgids # ensembl gene id

rm(snp_parts, snps_chr, snps_pos, gene_ensgids)
gc()

# sort pairs by snp locations (chr)
# take a batch of snps and test with target genes
chromosomes = sort(unique(snp_gene_pairs_to_test$snp_chr))
all_tested_eqtls = NULL
for(chr in chromosomes){
  verbose_print(sprintf('calling trans-eqtls for snps in %s',chr))
  chr_snp_gene_pairs = snp_gene_pairs_to_test[snp_gene_pairs_to_test$snp_chr == chr, , drop =F]
  if(nrow(chr_snp_gene_pairs) == 0) {
    next
  }
  
  chr_snps = unique(chr_snp_gene_pairs$snp)
  chr_genes = unique(chr_snp_gene_pairs$ensgid)
  
  ### read genotypes
  geno_fn = sprintf("%s%s%s", geno_pfx, chr, geno_sfx)
  geno_df = read_genotype_file(geno_fn)
  geno_df = geno_df[chr_snps, , drop = F]
  
  ### 
  chr_expr_df = expr_df[chr_genes, , drop = F]
  
  common_samples = Reduce(intersect, list(colnames(chr_expr_df), colnames(geno_df), colnames(cov_df)))
  meqtl_expr_slice = SlicedData$new(as.matrix(chr_expr_df[,common_samples,drop=F]))
  meqtl_cov_slice = SlicedData$new(as.matrix(cov_df[,common_samples,drop=F]))
  
  # split snps and call eqtls
  snp_splits = sort(unique(c(seq(max_snp_per_meqtl, nrow(geno_df), max_snp_per_meqtl), nrow(geno_df))))
  snp_splits = snp_splits[snp_splits<=nrow(geno_df)]
  for(ss in seq_len(length(snp_splits))){
    print(sprintf("calling matrix-eqtl for snps in %s: part %d of %d", chr, ss, length(snp_splits) ))
    start_idx = ifelse(ss==1, 1, snp_splits[ss-1]+1)
    end_idx = snp_splits[ss]
    # meqtl_snp_slice = SlicedData$new(geno_df[start_idx:end_idx,,drop=F])
    meqtl_snp_slice = SlicedData$new(as.matrix(geno_df[start_idx:end_idx,common_samples,drop=F]))
    
    ### run matrix-eQTL
    me = Matrix_eQTL_engine(snps = meqtl_snp_slice,
                            gene = meqtl_expr_slice,
                            cvrt = meqtl_cov_slice,
                            output_file_name = NULL,
                            pvOutputThreshold = 1,
                            useModel = modelLINEAR, 
                            verbose = FALSE,
                            pvalue.hist = FALSE,
                            min.pv.by.genesnp = FALSE,
                            noFDRsaveMemory = FALSE)
    
    me$all$eqtls$snps = as.character(me$all$eqtls$snps)  # convert factor to chracter
    me$all$eqtls$gene = as.character(me$all$eqtls$gene)  # convert factor to chracter
    
    ### filter tests not present in chr_snp_gene_pairs
    split_chr_snp_gene_pairs = chr_snp_gene_pairs[chr_snp_gene_pairs$snp %in% rownames(geno_df)[start_idx:end_idx], , drop = F]
    target_tests_df = merge(
      me$all$eqtls,
      split_chr_snp_gene_pairs,
      by.x = c("snps", "gene"),
      by.y = c("snp", "ensgid"),
      suffixes = c(".x", ".y")
    )
    colnames(target_tests_df) <- gsub(x = colnames(target_tests_df), pattern = sprintf("^gene.y$"), replacement = "gene_name")
    target_tests_df = target_tests_df[,-which(colnames(target_tests_df) %in% c("snp_chr", "snp_pos"))]
    all_tested_eqtls = rbind(all_tested_eqtls, target_tests_df)
    
    rm(me, meqtl_snp_slice, split_chr_snp_gene_pairs, target_tests_df)
    gc()
  }
}

### aggregate all tests
fl = flock::lock(repository_lock_fn)
if(file.exists(repository_fn)){
  repo_stats = readRDS(repository_fn)
  combined_stats = rbind(repo_stats, all_tested_eqtls)
  is.uniq = !duplicated(combined_stats[, c("snps", "gene"), drop = F])
  combined_stats = combined_stats[is.uniq, , drop = F]
} else {
  combined_stats = all_tested_eqtls
}
if(nrow(combined_stats) > 0){
  combined_stats$FDR = NA
}
saveRDS(combined_stats, file = repository_fn)
flock::unlock(fl)

### perform fdr and compute #eGenes
net_eqtl_stats = merge(
  combined_stats,
  all_snp_gene_pairs,
  by.x = c("snps", 'gene_name'),
  by.y = c("snp", "gene")
)

n_sig_egenes = 0
n_sig_esnps = 0
n_sig_eqtls = 0
net_sig_eqtls = NULL
if(nrow(net_eqtl_stats) > 0){
  net_eqtl_stats$FDR = p.adjust(net_eqtl_stats$p, method = "BH")  
  net_sig_eqtls = net_eqtl_stats[net_eqtl_stats$FDR <= 0.05, , drop = F]
  n_sig_egenes = length(net_sig_eqtls$gene)
  n_sig_esnps = length(net_sig_eqtls$snps)
  n_sig_eqtls = nrow(net_sig_eqtls)
}

### save results
res = list(
  n_tested_pairs = nrow(net_eqtl_stats),
  sig_eqtls = net_sig_eqtls,
  n_sig_egenes = n_sig_egenes,
  n_sig_esnps = n_sig_esnps,
  n_sig_eqtls = n_sig_eqtls
)
saveRDS(res, file = out_fn)

