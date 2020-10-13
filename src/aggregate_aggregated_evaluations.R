library(argparser)
library(miscutil)

args <- arg_parser("program");
args <- add_argument(args, "--dir", help="directory where tissues are separeted", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results")
args <- add_argument(args, "--tissue", help="comma-separated tissue names", default="Whole_Blood,Muscle_Skeletal")
args <- add_argument(args, "--correction", help="comma-separated correction labels", default="corrected")
args <- add_argument(args, "--gene_selection", help="comma-separated gene selection options", default="AXPVAR")
args <- add_argument(args, "--n_genes", help="comma-separated number of genes options", default="1500")
args <- add_argument(args, "--fn", help="aggregated evaluation file name", default="aggregated_evaluations_validation.rds")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/master_aggregated_evaluations_validation.rds")

### parse args
argv = parse_args(args)
res_dir = argv$dir
tissues = parse_delimitted_string(argv$tissue)
correction_labels = parse_delimitted_string(argv$correction)
gene_selections = parse_delimitted_string(argv$gene_selection)
n_genes = as.numeric(parse_delimitted_string(argv$n_genes))
file_name = parse_delimitted_string(argv$fn)
out_fn = argv$o

### aggregate
aggregated_df = NULL
for(tissue in tissues){
  for(correction_label in correction_labels){
    for(gene_selection in gene_selections){
      for(n_gene in n_genes){
        fn = sprintf("%s/%s/%s/%s/%s/%s", 
                     res_dir, 
                     tissue, 
                     correction_label, 
                     gene_selection, 
                     n_gene, 
                     file_name)
        if(!file.exists(fn)){
          next
        }
        cur_df = readRDS(fn)
        if(is.null(aggregated_df)){
          # initialize aggregated_df
          aggregated_df = data.frame(tissue = tissue, 
                                     correction_label = correction_label, 
                                     gene_selection = gene_selection, 
                                     n_gene = n_gene, 
                                     stringsAsFactors = F)
          aggregated_df = aggregated_df[rep(1, nrow(cur_df)), , drop = F]
          for(new_col in colnames(cur_df)){
            aggregated_df[,new_col] = cur_df[,new_col]
          }
        } else {
          # for each new column, put NA for previous methods
          for(new_col in setdiff(colnames(cur_df), colnames(aggregated_df))){
            aggregated_df[,new_col] = NA
          }
          # add current results with NA for unavailable columns
          cur_df2 = aggregated_df[rep(1, nrow(cur_df)),,drop=F]
          cur_df2[,] = NA
          cur_df2[,"tissue"] = tissue
          cur_df2[,"correction_label"] = correction_label
          cur_df2[,"gene_selection"] = gene_selection
          cur_df2[,"n_gene"] = n_gene
          cur_df2[,colnames(cur_df)] = cur_df
          aggregated_df = rbind(aggregated_df, cur_df2)
          rownames(aggregated_df) = NULL
        }
      }
    }
  }
}

### save
saveRDS(aggregated_df, file = out_fn)
