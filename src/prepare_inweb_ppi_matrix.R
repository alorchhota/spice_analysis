library(spice)
library(argparser)
library(ioutil)
library(reshape2)

args <- arg_parser("program");
args <- add_argument(args, "--gene", help="file with gene names", default="/work-zfs/abattle4/ashis/progres/spice_anlysis/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/5000/genes.txt")
args <- add_argument(args, "--inweb", help="inWeb file", default="/work-zfs/abattle4/lab_data/inWeb_inbiomap/human_interactions_with_evidence.tab")
args <- add_argument(args, "--pathway", help="Pathway gene file (rds). Interactions will be restricted to pathways only.", default="")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/inweb.rds")

### parse args
argv = parse_args(args)
gene_fn = argv$gene
inweb_fn = argv$inweb
pathway_fn = argv$pathway
out_fn = argv$o

### get gene names
gene_df = read_df(gene_fn, header = F, row.names = F, check.names = F)
genes = as.character(gene_df[,1])

### get inweb ppi scores
lines = readLines(inweb_fn)
line_splits = strsplit(x = lines, split = "\t")
inweb_df = sapply(line_splits, function(splits){
  gene1 = strsplit(splits[1], split = "_")[[1]][1]
  gene2 = strsplit(splits[2], split = "_")[[1]][1]
  score = as.numeric(splits[5])
  return(list(gene1 = gene1, gene2 = gene2, score = score))
})
colnames(inweb_df) = NULL
inweb_df = t(inweb_df)

inweb_df = as.data.frame(inweb_df)
inweb_df$gene1 = as.character(inweb_df$gene1)
inweb_df$gene2 = as.character(inweb_df$gene2)
inweb_df$score = as.numeric(inweb_df$score)

inweb_df = inweb_df[inweb_df$gene1 %in% genes, , drop = F ]
inweb_df = inweb_df[inweb_df$gene2 %in% genes, , drop = F ]

inweb_ppi_score_mat = suppressWarnings(acast(data = inweb_df,
                                             formula = gene1 ~ gene2,
                                             value.var = "score",
                                             fun.aggregate = max,
                                             na.rm = T))
inweb_ppi_score_mat[is.infinite(inweb_ppi_score_mat)] = NA

### create a matrix with all genes mapped in string
ppi_score_mat = matrix(NA,
                       nrow = length(genes),
                       ncol = length(genes),
                       dimnames = list(genes, genes))
ppi_score_mat[rownames(inweb_ppi_score_mat), colnames(inweb_ppi_score_mat)] = inweb_ppi_score_mat
ppi_score_mat = pmax(ppi_score_mat, t(ppi_score_mat), na.rm = T)
ppi_score_mat[is.na(ppi_score_mat)] = 0

### restrict to interactions where both genes are in same pathway
if(is.character(pathway_fn) && nchar(trimws(pathway_fn)) > 0){
  stopifnot(file.exists(pathway_fn))
  allowed_ppi_mat = matrix(
    0,
    nrow = nrow(ppi_score_mat),
    ncol = ncol(ppi_score_mat),
    dimnames = list(rownames(ppi_score_mat), colnames(ppi_score_mat))
  )
  
  pathways = readRDS(pathway_fn)
  for(pw in pathways){
    pw = intersect(pw, rownames(ppi_score_mat))
    if(length(pw)>1){
      allowed_ppi_mat[pw, pw] = 1
    }
  }
  
  ppi_score_mat = ppi_score_mat * allowed_ppi_mat
}

### save
saveRDS(ppi_score_mat, file = out_fn)
