# this script creates a map from each gene to it's nearby snps
library(argparser)
library(ioutil)
library(miscutil)

args <- arg_parser("program");
args <- add_argument(args, "--geno_pfx", help="genotype file name prefix with path (before chromosome in the name)", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_")
args <- add_argument(args, "--geno_chr", help="comma separated chromosomes", default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX")
args <- add_argument(args, "--geno_sfx", help="genotype file name suffix (after chromosome in the name)", default="_dosage_MAF_05_not_in_repeat.RData")
args <- add_argument(args, "--gene_annot", help="gene annotation file", default="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/gencode.v26.annotation.gene.txt")
args <- add_argument(args, "--d", help="max distance between a snp and a gene's TSS", default=1e6)
args <- add_argument(args, "--o", help="Output plot file", default="results/snps_nearby_gene.rds")

### parse args
argv = parse_args(args)
geno_pfx = argv$geno_pfx
geno_chromosomes_input = argv$geno_chr
geno_sfx = argv$geno_sfx
gene_annot_fn = argv$gene_annot
d = argv$d
out_fn = argv$o

### process settings and annotation
geno_chromosomes = parse_delimitted_string(geno_chromosomes_input, delim = ',', rm.empty = T)

### read gene annotations
gencode_df = read_df(gene_annot_fn, header = T, row.names = F)
# add TSS
tss_values =  as.integer(apply(gencode_df, MARGIN = 1, FUN = function(row){
  ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
}))
gencode_df$tss = tss_values
# remove genes with non-unique symbols
gene_counts = table(gencode_df[,"gene_name"])
non_unique_genes = names(gene_counts[gene_counts>1])
gencode_df = gencode_df[!(gencode_df$gene_name %in% non_unique_genes), ,drop = F]

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

### map gene to snps
gene_2_snps = list()
for(chr in geno_chromosomes){
  verbose_print(sprintf("processing %s ...", chr))
  
  ### read genotypes
  geno_fn = sprintf("%s%s%s", geno_pfx, chr, geno_sfx)
  geno_df = read_genotype_file(geno_fn)
  
  ### get snps and their locations
  snp_ids = rownames(geno_df)
  parts = strsplit(snp_ids, split = '_')
  snp_positions = sapply(parts, function(x) as.integer(x[2]))
  
  ### delete genotype matrix
  rm(geno_df)
  gc()
  
  ### genes in chr
  chr_genes = gencode_df[gencode_df$chr == chr, c('gene_name', 'tss')]
  
  ### for each, find snps
  for(gi in seq_len(nrow(chr_genes))){
    g = chr_genes[gi,1]
    g_pos = chr_genes[gi, 2]
    g_snps = snp_ids[(snp_positions >= g_pos-d) & (snp_positions <= g_pos+d)]
    gene_2_snps[[g]] = g_snps
  }
}

### save
saveRDS(gene_2_snps, file = out_fn)
