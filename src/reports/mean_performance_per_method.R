library(miscutil)
library(dplyr)

res_fn = "results/spice_results/gtex_v8/aggregated/all_evaluations_test.rds"
res = readRDS(res_fn)
head_df(res)

means_tbl = res %>%
  group_by(correction_label, gene_selection, n_gene, method) %>% 
  summarise_all(mean)

stopifnot(length(unique(means_tbl$method)) == nrow(means_tbl))

metrics = c(
  'string_aupr',
  'string_spearman_r',
  'string_precision_topNA_th0.7',
  'string_hub_aupr',
  'string_geodesic_rankdiff_2_1',
  'hallmark_shared_aupr',
  'hallmark_enriched_pathway')
   # 'n_sig_genes'

# means_tbl[,c('method', metrics)]

spice_means = means_tbl[means_tbl$method == 'spice', metrics]
wgcna_means = means_tbl[means_tbl$method == 'wgcna', metrics]
improvement_percentage = (spice_means - wgcna_means) / wgcna_means * 100
improvement_percentage
