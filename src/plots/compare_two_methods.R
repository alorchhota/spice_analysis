# this script create a plot to visualize relative performance of two methods.

library(argparser)
library(ggplot2)
library(miscutil)

args <- arg_parser("program");
args <- add_argument(args, "--res", help="aggregate results file (*.rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/aggregated/all_evaluations_test.rds")
args <- add_argument(args, "--method1", help="Method 1", default="spice")
args <- add_argument(args, "--method2", help="Method 2", default="wgcna")
args <- add_argument(args, "--metrics", help="list of metrics (comma-separated)", default="string_aupr,string_spearman_r,string_precision_topNA_th0.7,string_hub_aupr,hallmark_shared_aupr,hallmark_enriched_pathway,n_sig_egenes,string_geodesic_rankdiff_2_1,string_geodesic_rankdiff_3_2")
args <- add_argument(args, "--labels", help="labels for metrics (comma-separated)", default="Interaction,Spearman,Precision,Hub,Shared,#Pathways,#eGenes,Rank:2-1,Rank:3-2")
args <- add_argument(args, "--o", help="Output data file (*.rds)", default="results/spice_vs_wgcna.rds")
args <- add_argument(args, "--nsample", help="file with number of samples per tissue", default="")
args <- add_argument(args, "--plt", help="Output plot file (*.pdf)", default="results/spice_vs_wgcna.pdf")

### parse args
argv = parse_args(args)
res_fn = argv$res
method1 = argv$method1
method2 = argv$method2
metrics_input = argv$metrics
labels_input = argv$labels
nsample_fn = argv$nsample
plt_data_fn = argv$o
plt_fn = argv$plt

### process input
metrics = parse_delimitted_string(metrics_input, delim = ',', rm.empty = T)
labels = parse_delimitted_string(labels_input, delim = ',', rm.empty = T)
stopifnot(length(metrics) == length(labels))
names(labels) = metrics
metrics = labels

nsample_df = NULL
if(file.exists(nsample_fn)){
  nsample_df = read.table(nsample_fn, header = F, stringsAsFactors = F, row.names = 1)
  nsample_df = nsample_df[order(nsample_df[,1], decreasing = T), , drop = F]
}

# ### metrics to plot
# metrics = c(string_aupr = "Interaction", 
#             string_spearman_r = "Spearman",
#             string_precision_topNA_th0.7 = "Precision",
#             string_hub_aupr = "Hub",
#             # string_kegg_aupr = "Interaction AUPR (STRING-KEGG)", 
#             # string_kegg_spearman_r = "Spearman r (STRING-KEGG)",
#             # string_kegg_precision_topNA_th0.7 = "Precision (STRING-KEGG)",
#             # string_kegg_hub_aupr = "Hub AUPR (STRING-KEGG)",
#             # string_exp_aupr = "Interaction AUPR (STRING-EXP)", 
#             # string_exp_spearman_r = "Spearman r (STRING-EXP)",
#             # string_exp_precision_topNA_th0.7 = "Precision (STRING-EXP)",
#             # string_exp_hub_aupr = "Hub AUPR (STRING-EXP)",
#             # kegg_interaction_aupr = "Interaction AUPR (KEGG)", 
#             # kegg_interaction_spearman_r = "Spearman r (KEGG)",
#             # kegg_interaction_precision_topNA_th0.7 = "Precision (KEGG)",
#             # kegg_interaction_hub_aupr = "Hub AUPR (KEGG)",
#             # chipseq_interaction_aupr = "Interaction AUPR (Chip-seq)", 
#             # chipseq_interaction_spearman_r = "Spearman r (Chip-seq)",
#             # chipseq_interaction_precision_topNA_th0.7 = "Precision (Chip-seq)",
#             # chipseq_interaction_hub_aupr = "Hub AUPR (Chip-seq)",
#             # inweb_aupr = "Interaction AUPR (InWeb_IM)", 
#             # inweb_spearman_r = "Spearman r (InWeb_IM)",
#             # inweb_precision_topNA_th0.7 = "Precision (InWeb_IM)",
#             # inweb_hub_aupr = "Hub AUPR (InWeb_IM)",
#             # kegg_shared_aupr = "Shared pathway AUPR (KEGG)", 
#             hallmark_shared_aupr = "Shared", 
#             # reactome_shared_aupr = "Shared pathway AUPR (Reactome)", 
#             # go_shared_aupr = "Shared pathway AUPR (GO)",
#             # kegg_enriched_pathway = "#Enriched pathways (KEGG)",
#             hallmark_enriched_pathway = "#Pathways", 
#             # reactome_enriched_pathway = "#Enriched pathways (Reactome)",
#             # go_enriched_pathway = "#Enriched pathways (GO)",
#             n_sig_egenes = "#eGenes")


### read results data
res_df = readRDS(res_fn)
stopifnot(all(c(method1, method2) %in% unique(res_df$method)))
stopifnot(all(names(metrics) %in% colnames(res_df) ))

### prepare plot data
tissues = sort(unique(res_df$tissue))
if(!is.null(nsample_df)){
  stopifnot(all(tissues %in% rownames(nsample_df)))
  tissues = intersect(rownames(nsample_df), tissues)
}
methods = unique(res_df$method)

plt_df = NULL
for(tissue in tissues){
  tissue_res = res_df[res_df$tissue == tissue, ,drop=F]
  for(metric in names(metrics)){
    metric_res = tissue_res[,c("method", metric)]
    method1_value = metric_res[metric_res$method == method1, metric]
    method2_value = metric_res[metric_res$method == method2, metric]
    
    if(method2_value == 0){
      if(method1_value == 0){
        ratio = 1
      } else {
        ratio = Inf
      }
    } else {
      ratio = method1_value / method2_value
    }
    
    df = data.frame(tissue = rep(tissue,2),
                    method = c(method1, method2),
                    metric = rep(metric, 2),
                    value = c(method1_value, method2_value),
                    ratio = c(ratio, 1),
                    stringsAsFactors = F)
    plt_df = rbind(plt_df, df)
  }
}

plt_df$pair = sprintf("%s-%s", plt_df$tissue, plt_df$metric)
plt_df$improvement_percent = (plt_df$ratio - 1) * 100
improved = rep(0, nrow(plt_df))
improved[plt_df$ratio > 1] = 1
improved[plt_df$ratio < 1] = -1
plt_df$improved = improved
plt_df$metric_label = metrics[plt_df$metric]

plt_df$tissue = factor(plt_df$tissue, levels = tissues)
plt_df$method = factor(plt_df$method)
plt_df$metric = factor(plt_df$metric, levels = names(metrics))
plt_df$metric_label = factor(plt_df$metric_label, levels = as.character(metrics))
plt_df$pair = factor(plt_df$pair)
plt_df$improved = factor(plt_df$improved)

### save plot data
saveRDS(plt_df, file = plt_data_fn)

### plot
pdf(plt_fn, width = 8, height = 10)
plt_df[plt_df$improvement_percent >= 200, "improvement_percent"] = 200
ggplot(data = plt_df, aes(tissue, improvement_percent)) +
  theme_bw() +
  geom_point(aes(color=improved), show.legend = F) +
  geom_line(aes(group = pair, color = improved), show.legend = F) +
  facet_grid(metric_label ~ ., scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Tissue") +
  ylab("Improvement (%)")
dev.off()

