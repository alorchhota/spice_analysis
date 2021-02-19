library(argparser)
library(ggplot2)
library(cowplot)
library(data.table)
library(miscutil)

args <- arg_parser("program");
args <- add_argument(args, "--res", help="aggregate results file (*.rds)", default="results/spice_results/gtex_v8/aggregated/all_evaluations_test.rds")
args <- add_argument(args, "--methods", help="list of methods (comma-separated)", default="random,pcorr,glasso_likelihood,clr,aracne,et_genie3,mrnet,mrnetb,wgcna,spice")
args <- add_argument(args, "--method_labels", help="labels for methods (comma-separated)", default="Random,PCor,GLasso,CLR,ARACNE,GENIE3,MRNET,MRNETB,WGCNA,SPICE")
args <- add_argument(args, "--metrics", help="list of metrics (comma-separated)", default="string_aupr,string_spearman_r,string_precision_topNA_th0.7,string_hub_aupr,hallmark_shared_aupr,hallmark_enriched_pathway,n_sig_egenes,string_geodesic_rankdiff_2_1")
args <- add_argument(args, "--metric_labels", help="labels for metrics (comma-separated)", default="Interaction AUPR (STRING),Spearman r (STRING),Precision (STRING),Hub AUPR (STRING),Shared pathway AUPR (Hallmark),#Enriched pathways (Hallmark),#Trans-eGenes (GTEx),Rankdiff (2-1)")
args <- add_argument(args, "--o", help="Output plot data (*.rds)", default="results/method_comparison.rds")
args <- add_argument(args, "--plt", help="Output plot file (pdf)", default="results/method_comparison.pdf")

### parse args
argv = parse_args(args)
res_fn = argv$res
methods_input = argv$methods
method_labels_input = argv$method_labels
metrics_input = argv$metrics
metric_labels_input = argv$metric_labels
plt_data_fn = argv$o
plt_fn = argv$plt

### process input
metrics = parse_delimitted_string(metrics_input, delim = ',', rm.empty = T)
metric_labels = parse_delimitted_string(metric_labels_input, delim = ',', rm.empty = T)
stopifnot(length(metrics) == length(metric_labels))
names(metric_labels) = metrics
metrics = metric_labels

methods = parse_delimitted_string(methods_input, delim = ',', rm.empty = T)
method_labels = parse_delimitted_string(method_labels_input, delim = ',', rm.empty = T)
stopifnot(length(methods) == length(method_labels))
names(method_labels) = methods
methods = method_labels

### get pathways
res_df = readRDS(res_fn)

### prepare plot data
available_metric_names = intersect(names(metrics), colnames(res_df))
metrics = metrics[available_metric_names]
available_method_names = intersect(names(methods), res_df$method)
available_method_labels = methods[available_method_names]
method_df = data.frame(method = available_method_names,
                       Method = factor(available_method_labels, levels = available_method_labels), 
                       stringsAsFactors = F)
metric_df = data.frame(metric = names(metrics),
                       Metric = factor(metrics, levels = metrics), 
                       stringsAsFactors = F)
# res_df0 = res_df
res_df = res_df[res_df$method %in% available_method_names,
                colnames(res_df) %in% c(available_metric_names, 'method'), drop = F]

plt_df = melt(as.data.table(res_df), measure.vars = available_metric_names, value.name = "Value")

plt_df = merge(plt_df, method_df, by = "method", all.x = T, all.y = F)
plt_df = merge(plt_df, metric_df, by.x = "variable", by.y = "metric", all.x = T, all.y = F)

### save 
saveRDS(plt_df, file = plt_data_fn)

### beging plot
pdf(plt_fn, width = 8, height = 6)

### plot
p <- ggplot(plt_df, aes(x = Method, y = Value)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1) +
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "pointrange",
    fatten = 2,
    color = "red"
  ) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  ylab("Metric value")+
  facet_wrap(~ Metric, ncol = 2, scales = "free_y") 

print(p)

### close plots
dev.off()
