library(argparser)
library(ggplot2)
library(cowplot)

args <- arg_parser("program");
args <- add_argument(args, "--res", help="aggregate results file (*.rds)", default="/Users/ashissaha/github/spice_analysis/results/spice_results/gtex_v8/aggregated/all_evaluations_test.rds")
args <- add_argument(args, "--o", help="Output plot file (pdf)", default="results/mean_performance_per_method_trans_egenes.pdf")

### parse args
argv = parse_args(args)
res_fn = argv$res
out_fn = argv$o

### metrics to plot
metrics = c(string_aupr = "Interaction AUPR (STRING)", 
            string_spearman_r = "Spearman r (STRING)",
            string_precision_topNA_th0.7 = "Precision (STRING)",
            string_hub_aupr = "Hub AUPR (STRING)",
            string_kegg_aupr = "Interaction AUPR (STRING-KEGG)", 
            string_kegg_spearman_r = "Spearman r (STRING-KEGG)",
            string_kegg_precision_topNA_th0.7 = "Precision (STRING-KEGG)",
            string_kegg_hub_aupr = "Hub AUPR (STRING-KEGG)",
            string_exp_aupr = "Interaction AUPR (STRING-EXP)", 
            string_exp_spearman_r = "Spearman r (STRING-EXP)",
            string_exp_precision_topNA_th0.7 = "Precision (STRING-EXP)",
            string_exp_hub_aupr = "Hub AUPR (STRING-EXP)",
            kegg_interaction_aupr = "Interaction AUPR (KEGG)", 
            kegg_interaction_spearman_r = "Spearman r (KEGG)",
            kegg_interaction_precision_topNA_th0.7 = "Precision (KEGG)",
            kegg_interaction_hub_aupr = "Hub AUPR (KEGG)",
            chipseq_interaction_aupr = "Interaction AUPR (Chip-seq)", 
            chipseq_interaction_spearman_r = "Spearman r (Chip-seq)",
            chipseq_interaction_precision_topNA_th0.7 = "Precision (Chip-seq)",
            chipseq_interaction_hub_aupr = "Hub AUPR (Chip-seq)",
            inweb_aupr = "Interaction AUPR (InWeb_IM)", 
            inweb_spearman_r = "Spearman r (InWeb_IM)",
            inweb_precision_topNA_th0.7 = "Precision (InWeb_IM)",
            inweb_hub_aupr = "Hub AUPR (InWeb_IM)",
            kegg_shared_aupr = "Shared pathway AUPR (KEGG)", 
            hallmark_shared_aupr = "Shared pathway AUPR (Hallmark)", 
            reactome_shared_aupr = "Shared pathway AUPR (Reactome)", 
            go_shared_aupr = "Shared pathway AUPR (GO)",
            kegg_enriched_pathway = "#Enriched pathways (KEGG)",
            hallmark_enriched_pathway = "#Enriched pathways (Hallmark)", 
            reactome_enriched_pathway = "#Enriched pathways (Reactome)",
            go_enriched_pathway = "#Enriched pathways (GO)",
            n_sig_egenes = "#Trans-eGenes")


### get pathways
res_df = readRDS(res_fn)

### plot function
get_comparison_plots <- function(plt_df,
                                 xvar,
                                 metrics = c(
                                   string_aupr = "Interaction AUPR (STRING)",
                                   spearman_r = "Spearman r (STRING)",
                                   string_precision_topNA_th0.7 = "Precision (STRING)",
                                   string_hub_aupr = "Hub AUPR (STRING)"
                                 ),
                                 title = "") {
  plotlist <- lapply(names(metrics), function(metric) {
    p <- ggplot(plt_df, aes_string(x = xvar, y = metric)) +
      geom_violin(trim = FALSE) +
      theme_bw() +
      geom_dotplot(binaxis = "y", stackdir = "center") +
      stat_summary(
        fun.data = "mean_sdl",
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = "red"
      ) +
      ggtitle(title) + 
      ylab(metrics[metric]) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      xlab("")
    p
  })
  return(plotlist)
}

### start plots
pdf(out_fn, width = 10, height = 6)

### plot
available_metric_names = intersect(names(metrics), colnames(res_df))
metrics = metrics[available_metric_names]
for(mi in seq(from = 1, to = length(metrics), by = 4)){
  p <- cowplot::plot_grid(
    plotlist = get_comparison_plots(
      res_df,
      xvar = "method",
      metrics = metrics[mi:min(mi+3, length(metrics))],
      title = ""
    ),
    ncol = 2
  )
  print(p)
}

### close plots
dev.off()

