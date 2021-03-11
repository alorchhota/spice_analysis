library(argparser)
library(ggplot2)
library(cowplot)

args <- arg_parser("program");
args <- add_argument(args, "--rankdiff", help="rankdiff plot data file (*.rds)", default="results/rankdiff_plot.rds")
args <- add_argument(args, "--egene", help="number of egene plot data file (*.rds)", default="results/total_trans_egene_plot.rds")
args <- add_argument(args, "--plt", help="Output plot file (*.pdf)", default="results/rankdiff_negene.pdf")


### parse args
argv = parse_args(args)
egene_fn = argv$egene
rankdiff_fn = argv$rankdiff
plt_fn = argv$plt


### rankdiff plot
rankdiff_plt_df = readRDS(rankdiff_fn)
rankdiff_plt <- ggplot(rankdiff_plt_df, aes(x = metric, y = rankdiff_mean)) +
  theme_bw() +
  geom_bar(aes(fill = method), stat = "identity",
           position = position_dodge(0.8), width = 0.7, show.legend = F)+
  geom_errorbar(
    aes(ymin = rankdiff_mean-rankdiff_sd, ymax = rankdiff_mean+rankdiff_sd, group = method),
    width = 0.2, position = position_dodge(0.8)
  ) +
  xlab("Distances between genes") +
  ylab("Rank difference") +
  labs(fill = "Method")

# print(rankdiff_plt)

### egene plot
egene_plt_df = readRDS(egene_fn)
egene_plt <- ggplot(egene_plt_df, aes(x = method_label, y = n_egene)) +
  theme_bw() +
  geom_bar(aes(fill = method_label), stat = "identity",
           position = position_dodge(0.8), width = 0.7, show.legend = T) +
  xlab("Method") +
  ylab("Total number of trans-eGenes") +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(axis.text.x = element_blank()) +
  #theme(legend.title = element_blank()) +
  labs(fill = "Method")

# plot(egene_plt)

### combine plots and save
pdf(plt_fn, width = 8, height = 4)
plot_grid(rankdiff_plt, egene_plt, labels = c('A', 'B'), rel_widths = c(3, 2))
dev.off()
