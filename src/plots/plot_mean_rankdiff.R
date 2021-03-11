library(argparser)
library(ggplot2)
library(cowplot)
library(data.table)

args <- arg_parser("program");
args <- add_argument(args, "--res", help="aggregate results file (*.rds)", default="results/spice_results/gtex_v8/aggregated/all_evaluations_test.rds")
args <- add_argument(args, "--o", help="Output plot data file (*.rds)", default="results/rankdiff_comparison.rds")
args <- add_argument(args, "--plt", help="Output plot file (*.pdf)", default="results/rankdiff_comparison.pdf")

### parse args
argv = parse_args(args)
res_fn = argv$res
plt_data_fn = argv$o
plt_fn = argv$plt

### methods to plot
methods = c(
  random = "Random",
  pcorr = "PCor",
  glasso_likelihood = "GLasso",
  clr = "CLR",
  aracne = "ARACNE",
  et_genie3 = "GENIE3",
  mrnet = "MRNET",
  mrnetb = "MRNETB",
  wgcna = "WGCNA",
  spice = "SPICE",
  spice_Mp_G.8_Ra_A1_I200="SPICE_I200",
  spice_Mp_G1_Ra_A1="SPICE_G1"
)

### metrics to plot
metrics = c(string_geodesic_rankdiff_2_1 = "2-1", 
            string_geodesic_rankdiff_3_2 = "3-2",
            string_geodesic_rankdiff_4_3 = "4-3")


### get pathways
res_df = readRDS(res_fn)

### plot
available_metric_names = intersect(names(metrics), colnames(res_df))
metrics = metrics[available_metric_names]
available_method_names = intersect(names(methods), res_df$method)
available_method_labels = methods[available_method_names]
method_df = data.frame(
  method = available_method_names,
  Method = factor(available_method_labels, levels = available_method_labels),
  stringsAsFactors = F
)
metric_df = data.frame(
  metric = names(metrics),
  Metric = factor(metrics, levels = metrics),
  stringsAsFactors = F
)
# res_df0 = res_df
res_df = res_df[res_df$method %in% available_method_names,
                colnames(res_df) %in% c(available_metric_names, 'method'), drop = F]

# compute mean and standard error for each metric
plt_df = NULL
for(metric in names(metrics)){
  means = tapply(res_df[,metric], res_df$method, mean)
  sds = tapply(res_df[,metric], res_df$method, sd)
  ns = tapply(res_df[,metric], res_df$method, length)
  
  means = means[names(methods)]
  sds = sds[names(methods)]
  ns = ns[names(methods)]
  stderrs = sds / sqrt(ns)
  
  stat_df = data.frame(metric = as.character(metrics[metric]),
                       method = as.character(methods),
                       rankdiff_mean = as.numeric(means),
                       rankdiff_sd = as.numeric(sds),
                       rankdiff_stderr =  as.numeric(stderrs),
                       stringsAsFactors = F)
  plt_df = rbind(plt_df, stat_df)
}

plt_df$method = factor(plt_df$method, levels = as.character(methods))

### save plot data
saveRDS(plt_df, file = plt_data_fn)

### plot
pdf(plt_fn)

p <- ggplot(plt_df, aes(x = metric, y = rankdiff_mean)) +
  theme_bw() +
  geom_bar(aes(fill = method), stat = "identity",
           position = position_dodge(0.8), width = 0.7)+
  geom_errorbar(
    aes(ymin = rankdiff_mean-rankdiff_sd, ymax = rankdiff_mean+rankdiff_sd, group = method),
    width = 0.2, position = position_dodge(0.8)
  ) +
  ylab("Difference in mean rank") +
  xlab("Known geodesic distances")

print(p)

dev.off()
