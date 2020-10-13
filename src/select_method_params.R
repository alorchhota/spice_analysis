library(argparser)
library(ggplot2)
library(cowplot)

args <- arg_parser("program");
args <- add_argument(args, "--res", help="aggregate results file (*.rds)", default="results/validation/all_evaluations_validation.rds")
args <- add_argument(args, "--o", help="Output plot file (pdf)", default="results/best_param_per_method.pdf")

### parse args
argv = parse_args(args)
res_fn = argv$res
out_fn = argv$o

metrics = c("string_aupr", "string_hub_aupr", "spearman_r", "string_precision_topNA_th0.7", "kegg_shared_aupr", "reactome_shared_aupr", "kegg_enriched_pathway", "reactome_enriched_pathway")

### get pathways
res_df = readRDS(res_fn)

### plot function
get_comparison_plots <- function(plt_df,
                                 xvar,
                                 metrics = c(
                                   "string_aupr",
                                   "string_hub_aupr",
                                   "spearman_r",
                                   "string_precision_topNA_th0.7",
                                   "kegg_shared_aupr",
                                   "reactome_shared_aupr",
                                   "kegg_enriched_pathway",
                                   "reactome_enriched_pathway"
                                 ),
                                 title = "") {
  plotlist <- lapply(metrics, function(metric) {
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
      ggtitle(title)
    p
  })
  return(plotlist)
}

### start plots
pdf(out_fn)

############ wgcna parameters############
wgcna_df = res_df[grepl("wgcna", res_df$method, fixed = TRUE), , drop = F]
wgcna_df$cor = sapply(wgcna_df$method, function(s) strsplit(s, split = "_")[[1]][2])
wgcna_df$sign = sapply(wgcna_df$method, function(s) strsplit(s, split = "_")[[1]][3])
wgcna_df$rsquare = sapply(wgcna_df$method, function(s) strsplit(s, split = "_")[[1]][4])

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    wgcna_df,
    xvar = "cor",
    metrics = metrics,
    title = "wgcna:cor"
  ),
  ncol = 3
)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    wgcna_df,
    xvar = "sign",
    metrics = metrics,
    title = "wgcna:sign"
  ),
  ncol = 3
)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    wgcna_df,
    xvar = "rsquare",
    metrics = metrics,
    title = "wgcna:rsquare"
  ),
  ncol = 3
)

### selected param: pearson, signed, r2 = 0.6?
selected_method = "wgcna_pearson_signed_0.6"

mean_performance = aggregate(wgcna_df[, metrics], list(wgcna_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)


############ genie3 parameters############
genie3_df = res_df[grepl("genie3", res_df$method, fixed = TRUE), , drop = F]

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    genie3_df,
    xvar = "method",
    metrics = metrics,
    title = "genie3:method"
  ),
  ncol = 3
)

### selected param: et_genie3
selected_method = "et_genie3"

mean_performance = aggregate(genie3_df[, metrics], list(genie3_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)


############ aracne parameters############
aracne_df = res_df[grepl("aracne", res_df$method, fixed = TRUE), , drop = F]
aracne_df$cor = sapply(aracne_df$method, function(s) strsplit(s, split = "_")[[1]][2])
aracne_df$eps = sapply(aracne_df$method, function(s) strsplit(s, split = "_")[[1]][3])
aracne_df$disc = sapply(aracne_df$method, function(s) strsplit(s, split = "_")[[1]][4])

# pdf(out_fn)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    aracne_df,
    xvar = "cor",
    metrics = metrics,
    title = "aracne:cor"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    aracne_df,
    xvar = "eps",
    metrics = metrics,
    title = "aracne:eps"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    aracne_df,
    xvar = "disc",
    metrics = metrics,
    title = "aracne:disc"
  ),
  ncol = 3
)

### selected param
selected_method = "aracne_pearson_e.1_freq"

mean_performance = aggregate(aracne_df[, metrics], list(aracne_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)


############ mrnetb parameters############
mrnetb_df = res_df[grepl("mrnetb", res_df$method, fixed = TRUE), , drop = F]
mrnetb_df$cor = sapply(mrnetb_df$method, function(s) strsplit(s, split = "_")[[1]][2])
mrnetb_df$disc = sapply(mrnetb_df$method, function(s) strsplit(s, split = "_")[[1]][3])

# pdf(out_fn)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    mrnetb_df,
    xvar = "cor",
    metrics = metrics,
    title = "mrnetb:cor"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    mrnetb_df,
    xvar = "disc",
    metrics = metrics,
    title = "mrnetb:disc"
  ),
  ncol = 3
)

### selected param
selected_method = "mrnetb_pearson_freq"

mean_performance = aggregate(mrnetb_df[, metrics], list(mrnetb_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)

############ mrnet parameters############
mrnet_df = res_df[grepl("mrnet_", res_df$method, fixed = TRUE), , drop = F]
mrnet_df$cor = sapply(mrnet_df$method, function(s) strsplit(s, split = "_")[[1]][2])
mrnet_df$disc = sapply(mrnet_df$method, function(s) strsplit(s, split = "_")[[1]][3])

# pdf(out_fn)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    mrnet_df,
    xvar = "cor",
    metrics = metrics,
    title = "mrnet:cor"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    mrnet_df,
    xvar = "disc",
    metrics = metrics,
    title = "mrnet:disc"
  ),
  ncol = 3
)

### selected param
selected_method = "mrnet_pearson_freq"

mean_performance = aggregate(mrnet_df[, metrics], list(mrnet_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)


############ clr parameters############
clr_df = res_df[grepl("clr_", res_df$method, fixed = TRUE), , drop = F]
clr_df$cor = sapply(clr_df$method, function(s) strsplit(s, split = "_")[[1]][2])
clr_df$disc = sapply(clr_df$method, function(s) strsplit(s, split = "_")[[1]][3])

# pdf(out_fn)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    clr_df,
    xvar = "cor",
    metrics = metrics,
    title = "clr:cor"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    clr_df,
    xvar = "disc",
    metrics = metrics,
    title = "clr:disc"
  ),
  ncol = 3
)

### selected param
selected_method = "clr_pearson_freq"

mean_performance = aggregate(clr_df[, metrics], list(clr_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)

############ pcorr parameters############
pcorr_df = res_df[grepl("pcorr_", res_df$method, fixed = TRUE), , drop = F]
pcorr_df$cor = sapply(pcorr_df$method, function(s) strsplit(s, split = "_")[[1]][2])

# pdf(out_fn)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    pcorr_df,
    xvar = "cor",
    metrics = metrics,
    title = "pcorr:cor"
  ),
  ncol = 3
)

### selected param
selected_method = "pcorr_spearman"

mean_performance = aggregate(pcorr_df[, metrics], list(pcorr_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)


############ spice parameters############
spice_df = res_df[grepl("spice_", res_df$method, fixed = TRUE), , drop = F]
spice_df$cor = sapply(spice_df$method, function(s) strsplit(s, split = "_")[[1]][2])
spice_df$gene_frac = sapply(spice_df$method, function(s) strsplit(s, split = "_")[[1]][3])
spice_df$rank = sapply(spice_df$method, function(s) strsplit(s, split = "_")[[1]][4])
spice_df$adjust = sapply(spice_df$method, function(s) strsplit(s, split = "_")[[1]][5])

# pdf(out_fn)
cowplot::plot_grid(
  plotlist = get_comparison_plots(
    spice_df,
    xvar = "cor",
    metrics = metrics,
    title = "spice:cor"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    spice_df,
    xvar = "gene_frac",
    metrics = metrics,
    title = "spice:gene_frac"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    spice_df,
    xvar = "rank",
    metrics = metrics,
    title = "spice:rank"
  ),
  ncol = 3
)

cowplot::plot_grid(
  plotlist = get_comparison_plots(
    spice_df,
    xvar = "adjust",
    metrics = metrics,
    title = "spice:adjust"
  ),
  ncol = 3
)

### selected param
selected_method = "spice_Mp_G.8_Ra_A1"

mean_performance = aggregate(spice_df[, metrics], list(spice_df$method), mean)
max_vals = apply(mean_performance[, metrics, drop=F], MARGIN = 2, max)
selected_vals = apply(mean_performance[mean_performance$Group.1 == selected_method, metrics,drop=F], MARGIN = 2, max)

print(sprintf("selected_method: %s", selected_method))
print("mean performance per method:")
print(mean_performance)
print("selected - max: ")
print(selected_vals - max_vals)

### close plots
dev.off()
