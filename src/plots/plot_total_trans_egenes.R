library(argparser)
library(ggplot2)
library(data.table)

args <- arg_parser("program");
#args <- add_argument(args, "--res", help="aggregate results file (*.rds)", default="/work-zfs/abattle4/ashis/progres/spice_anlysis/gtex_v8/aggregated/all_evaluations_test.rds")
args <- add_argument(args, "--res", help="aggregate results file (*.rds)", default="results/spice_results/gtex_v8/aggregated/all_evaluations_test.rds")
args <- add_argument(args, "--o", help="Output plot data file (*.rds)", default="results/trans_egenes_comparison.rds")
args <- add_argument(args, "--plt", help="Output plot file (pdf)", default="results/trans_egenes_comparison.pdf")

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
  spice = "SPICE"
)

### metrics to plot
metrics = c(n_sig_egenes = "#Trans-eGenes")

### read inputs
res_df = readRDS(res_fn)
stopifnot(all(names(metrics) %in% colnames(res_df) ))

### prepare plot data
tissues = sort(unique(res_df$tissue))

# check if trans-eqtl called exactly once for each tissue-method pair.
tissue_count_per_method = tapply(res_df$tissue, res_df$method, length)
stopifnot(all(tissue_count_per_method == length(tissues)))

# add number of egenes
total_n_egenes = tapply(res_df[,names(metrics)], res_df$method, sum, na.rm = T)
total_n_egenes = total_n_egenes[intersect(names(methods), names(total_n_egenes))]

plt_df = data.frame(
  method = names(total_n_egenes),
  method_label = methods[names(total_n_egenes)],
  n_egene = as.numeric(total_n_egenes),
  stringsAsFactors = F
)
plt_df$method = factor(plt_df$method, levels = names(methods))
plt_df$method_label = factor(plt_df$method_label, levels = as.character(methods))

### save plot data
saveRDS(plt_df, file = plt_data_fn)

### start plots
pdf(plt_fn)

### plot
p <- ggplot(plt_df, aes(x = method_label, y = n_egene)) +
  theme_bw() +
  geom_bar(aes(fill = method_label), stat = "identity",
           position = position_dodge(0.8), width = 0.7, show.legend = F) +
  xlab("Method") +
  ylab("Total number of trans-eGenes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(p)

### close plots
dev.off()
