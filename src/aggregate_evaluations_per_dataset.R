library(argparser)
library(miscutil)

args <- arg_parser("program");
args <- add_argument(args, "--dir", help="results directory for a dataset", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/results/Muscle_Skeletal/corrected/AXPVAR/1500")
args <- add_argument(args, "--methods", help="name of methods. comma-separated values", default="pearson,wgcna_pearson_unsigned_0.8")
args <- add_argument(args, "--o", help="Output file (rds)", default="results/summary.rds")

### parse args
argv = parse_args(args)
res_dir = argv$dir
methods_input = argv$methods
out_fn = argv$o

methods = parse_delimitted_string(methods_input, delim = ",", rm.empty = T)

### check variables
stopifnot(dir.exists(res_dir))
stopifnot(dir.exists(dirname(out_fn)))

### read string ppi auc results
get_string_auc_res <- function(res_fn){
  if(!file.exists(res_fn)){
    return(list())
  }
  res = readRDS(res_fn)
  res_list = list(string_aupr = res$pr$auc.integral,
                  string_auroc = res$roc$auc)
  return(res_list)
}

### read string ppi hub auc results
get_string_hub_auc_res <- function(res_fn){
  if(!file.exists(res_fn)){
    return(list())
  }
  res = readRDS(res_fn)
  res_list = list(string_hub_aupr = res$pr$auc.integral,
                  string_hub_auroc = res$roc$auc)
  return(res_list)
}

### read string ppi spearman cor
get_string_spearman_cor_res <- function(res_fn){
  if(!file.exists(res_fn)){
    return(list())
  }
  res = readRDS(res_fn)
  res_list = list(spearman_r = as.numeric(res$estimate))
  return(res_list)
}

### read string ppi precision
get_string_precision <- function(res_fn){
  if(!file.exists(res_fn)){
    return(list())
  }
  res = readRDS(res_fn)
  var_names_df = merge(rownames(res), colnames(res))
  var_names = matrix(sprintf("string_precision_top%s_th%s", 
                               as.character(var_names_df[,2]),  
                               as.character(var_names_df[,1])), 
                       nrow = nrow(res))
  res_list = lapply(as.numeric(res), function(x)x)
  names(res_list) = var_names
  return(res_list)
}

### read shared pathway auc results
get_shared_pathway_auc_res <- function(res_fn, pathway_name){
  if(!file.exists(res_fn)){
    return(list())
  }
  res = readRDS(res_fn)
  res_list = list(res$pr$auc.integral, 
                  res$roc$auc)
  names(res_list) = c(sprintf("%s_shared_aupr", pathway_name), 
                      sprintf("%s_shared_auroc", pathway_name))
  return(res_list)
}

### read pathway enrichment results
get_pathway_enrichment_res <- function(res_fn, pathway_name){
  if(!file.exists(res_fn)){
    return(list())
  }
  res = readRDS(res_fn)
  res_list = list(sum(res$fdr <= 0.05))
  names(res_list) = sprintf("%s_enriched_pathway", pathway_name)
  return(res_list)
}

### aggregate results per tissue
aggregated_df = NULL
for(method in methods){
  ### aggregate all results
  all_res = unlist(list(
    get_string_auc_res(res_fn = sprintf("%s/%s_string_ppi_auc.rds", res_dir, method)),
    get_string_hub_auc_res(res_fn = sprintf("%s/%s_string_ppi_hub_auc.rds", res_dir, method)),
    get_string_spearman_cor_res(res_fn = sprintf("%s/%s_string_ppi_spearman_cor.rds", res_dir, method)),
    get_string_precision(res_fn = sprintf("%s/%s_string_ppi_precision.rds", res_dir, method)),
    lapply(c("hallmark", "kegg", "biocarta", "reactome", "go"), function(pathway){
      get_shared_pathway_auc_res(
        res_fn = sprintf("%s/%s_shared_pathway_auc_%s.rds", res_dir, method, pathway),
        pathway_name = pathway)
    }),
    lapply(c("hallmark", "kegg", "biocarta", "reactome", "go"), function(pathway){
      get_pathway_enrichment_res(
        res_fn = sprintf("%s/%s_pathway_enrichment_%s.rds", res_dir, method, pathway),
        pathway_name = pathway)
    })
  ))
  
  ### update aggregated_df
  if(is.null(aggregated_df)){
    # initialize aggregated_df
    aggregated_df = data.frame(method = method, stringsAsFactors = F)
    for(new_col in names(all_res)){
      aggregated_df[,new_col] = all_res[new_col]
    }
  } else {
    # for each new column, put NA for previous methods
    for(new_col in setdiff(names(all_res), colnames(aggregated_df))){
      aggregated_df[,new_col] = NA
    }
    # add current results with NA for unavailable columns
    cur_df = aggregated_df[1,,drop=F]
    cur_df[1,] = NA
    cur_df[,"method"] = method
    cur_df[,names(all_res)] = all_res
    aggregated_df = rbind(aggregated_df, cur_df)
  }
}

rownames(aggregated_df) = NULL

### save
saveRDS(aggregated_df, file = out_fn)
