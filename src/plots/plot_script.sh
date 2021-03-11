#########################################################################
### compare spice and wgcna #############################################
#########################################################################

aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
nsample_fn="results/nsamples_per_tissue.txt"
# nsample_fn="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/nsamples_per_tissue.txt"
method1="spice"
method2="wgcna"
metrics="string_aupr,string_spearman_r,string_precision_topNA_th0.7,string_hub_aupr,hallmark_shared_aupr,hallmark_enriched_pathway,n_sig_egenes,string_geodesic_rankdiff_2_1,string_geodesic_rankdiff_3_2"
labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways,#eGenes,Rank:2-1,Rank:3-2"
outpfx="results/${method1}_vs_${method2}"

Rscript src/plots/compare_two_methods.R \
  --res "$aggregated_res_fn" \
  --method1 "$method1" \
  --method2 "$method2" \
  --metrics "$metrics" \
  --labels "$labels" \
  --nsample "$nsample_fn" \
  --o "$outpfx.rds" \
  --plt "$outpfx.pdf"
  

#########################################################################
### compare spice and wgcna using kegg-related resources  ###############
#########################################################################

aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
nsample_fn="results/nsamples_per_tissue.txt"
# nsample_fn="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/nsamples_per_tissue.txt"
method1="spice"
method2="wgcna"
# metrics="kegg_interaction_aupr,kegg_interaction_spearman_r,kegg_interaction_precision_topNA_th0.7,kegg_interaction_hub_aupr,kegg_shared_aupr,kegg_enriched_pathway"
# labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways"
metrics="kegg_interaction_aupr,kegg_interaction_spearman_r,kegg_interaction_precision_topNA_th0.7,kegg_interaction_hub_aupr,kegg_shared_aupr,kegg_enriched_pathway,kegg_interaction_geodesic_rankdiff_2_1,kegg_interaction_geodesic_rankdiff_3_2"
labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways,Rank:2-1,Rank:3-2"
outpfx="results/${method1}_vs_${method2}_kegg"

Rscript src/plots/compare_two_methods.R \
  --res "$aggregated_res_fn" \
  --method1 "$method1" \
  --method2 "$method2" \
  --metrics "$metrics" \
  --labels "$labels" \
  --nsample "$nsample_fn" \
  --o "$outpfx.rds" \
  --plt "$outpfx.pdf"
  
#########################################################################
### compare spice and genie3 ############################################
#########################################################################

aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
nsample_fn="results/nsamples_per_tissue.txt"
# nsample_fn="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/nsamples_per_tissue.txt"
method1="spice"
method2="et_genie3"
metrics="string_aupr,string_spearman_r,string_precision_topNA_th0.7,string_hub_aupr,hallmark_shared_aupr,hallmark_enriched_pathway,n_sig_egenes,string_geodesic_rankdiff_2_1,string_geodesic_rankdiff_3_2"
labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways,#eGenes,Rank:2-1,Rank:3-2"
outpfx="results/${method1}_vs_${method2}"

Rscript src/plots/compare_two_methods.R \
  --res "$aggregated_res_fn" \
  --method1 "$method1" \
  --method2 "$method2" \
  --metrics "$metrics" \
  --labels "$labels" \
  --nsample "$nsample_fn" \
  --o "$outpfx.rds" \
  --plt "$outpfx.pdf"
  
#########################################################################
### compare spice and glasso ############################################
#########################################################################

aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
nsample_fn="results/nsamples_per_tissue.txt"
# nsample_fn="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/nsamples_per_tissue.txt"
method1="spice"
method2="glasso_likelihood"
metrics="string_aupr,string_spearman_r,string_precision_topNA_th0.7,string_hub_aupr,hallmark_shared_aupr,hallmark_enriched_pathway,n_sig_egenes,string_geodesic_rankdiff_2_1,string_geodesic_rankdiff_3_2"
labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways,#eGenes,Rank:2-1,Rank:3-2"
outpfx="results/${method1}_vs_${method2}"

Rscript src/plots/compare_two_methods.R \
  --res "$aggregated_res_fn" \
  --method1 "$method1" \
  --method2 "$method2" \
  --metrics "$metrics" \
  --labels "$labels" \
  --nsample "$nsample_fn" \
  --o "$outpfx.rds" \
  --plt "$outpfx.pdf"

#########################################################################
### compare spice and glasso ############################################
#########################################################################

aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
nsample_fn="results/nsamples_per_tissue.txt"
# nsample_fn="/work-zfs/abattle4/ashis/progdata/scnet/gtex_v8_net/misc/nsamples_per_tissue.txt"
method1="spice"
method2="spice_Mp_G1_Ra_A1"
# method2="wgcna"
# metrics="string_aupr,string_spearman_r,string_precision_topNA_th0.7,string_hub_aupr,hallmark_shared_aupr,hallmark_enriched_pathway,n_sig_egenes,string_geodesic_rankdiff_2_1,string_geodesic_rankdiff_3_2"
# labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways,#eGenes,Rank:2-1,Rank:3-2"
metrics="string_aupr,string_spearman_r,string_precision_topNA_th0.7,string_hub_aupr,hallmark_shared_aupr,hallmark_enriched_pathway,string_geodesic_rankdiff_2_1,string_geodesic_rankdiff_3_2"
labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways,Rank:2-1,Rank:3-2"
outpfx="results/${method1}_vs_${method2}"

Rscript src/plots/compare_two_methods.R \
  --res "$aggregated_res_fn" \
  --method1 "$method1" \
  --method2 "$method2" \
  --metrics "$metrics" \
  --labels "$labels" \
  --nsample "$nsample_fn" \
  --o "$outpfx.rds" \
  --plt "$outpfx.pdf"


#########################################################################
### combine #trans-egenes and rankdiff plots ############################
#########################################################################

aggregated_res_fn="results/all_evals/all_evaluations_test.rds"

egene_plt_data_fn="results/total_trans_egene_plot.rds"
egene_plt_fn="results/total_trans_egene_plot.pdf"

rankdiff_plt_data_fn="results/rankdiff_plot.rds"
rankdiff_plt_fn="results/rankdiff_plot.pdf"

combined_plt_fn="results/trans_egene_rankdiff.pdf"

# call trans-egene
Rscript src/plots/plot_total_trans_egenes.R   \
  --res "$aggregated_res_fn" \
  --o "$egene_plt_data_fn" \
  --plt "$egene_plt_fn"

# call rankdiff
Rscript src/plots/plot_mean_rankdiff.R   \
  --res "$aggregated_res_fn" \
  --o "$rankdiff_plt_data_fn" \
  --plt "$rankdiff_plt_fn"
  
# combine
Rscript src/plots/combine_eqtl_and_rankdiff_plot.R \
  --egene "$egene_plt_data_fn" \
  --rankdiff "$rankdiff_plt_data_fn" \
  --plt "$combined_plt_fn"
  
 
#########################################################################
### compare multiple methods ############################################
#########################################################################

### R-helper to generate comma-separated metric labels
get_comma_separated_metrics <- function(x) paste(names(x), collapse = ",")
get_comma_separated_metric_labels <- function(x) paste(as.character(x), collapse = ",")
metrics_to_plot = c(
  kegg_interaction_aupr = "Interaction AUPR (KEGG)",
  kegg_interaction_precision_topNA_th0.01 = "Precision (KEGG)",
  kegg_interaction_hub_aupr = "Hub AUPR (KEGG)",
  kegg_interaction_geodesic_rankdiff_2_1 = "Rankdiff 2-1 (KEGG)",
  
  inweb_aupr = "Interaction AUPR (InWeb_IM)",
  inweb_spearman_r = "Spearman r (InWeb_IM)",
  inweb_precision_topNA_th0.01 = "Precision (InWeb_IM)",
  inweb_hub_aupr = "Hub AUPR (InWeb_IM)",
  inweb_geodesic_rankdiff_2_1 = "Rankdiff 2-1 (InWeb_IM)",
  
  string_exp_aupr = "Interaction AUPR (STRING-experiment)",
  string_spearman_r = "Spearman r (STRING-experiment)",
  string_exp_precision_topNA_th0.01 = "Precision (STRING-experiment)",
  string_exp_hub_aupr = "Hub AUPR (STRING-experimental)",
  string_exp_geodesic_rankdiff_2_1 = "Rankdiff 2-1 (STRING-experiment)",
  
  kegg_shared_aupr = "Shared pathway AUPR (KEGG)",
  reactome_shared_aupr = "Shared pathway AUPR (REACTOME)",
  go_shared_aupr = "Shared pathway AUPR (GO)",
  
  kegg_enriched_pathway = "#Enriched pathways (KEGG)",
  reactome_enriched_pathway = "#Enriched pathways (REACTOME)",
  go_enriched_pathway = "#Enriched pathways (GO)"
)

### kegg_interaction
aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
methods="random,pcorr,glasso_likelihood,clr,aracne,et_genie3,mrnet,mrnetb,wgcna,spice"
method_labels="Random,PCor,GLasso,CLR,ARACNE,GENIE3,MRNET,MRNETB,WGCNA,SPICE"
metrics="kegg_interaction_aupr,kegg_interaction_precision_topNA_th0.01,kegg_interaction_hub_aupr,kegg_interaction_geodesic_rankdiff_2_1"
metric_labels="Interaction AUPR (KEGG),Precision (KEGG),Hub AUPR (KEGG),Rankdiff 2-1 (KEGG)"
outpfx="results/method_comparison_kegg_interaction"

Rscript src/plots/compare_multiple_methods.R --res "$aggregated_res_fn" \
                                             --methods "$methods" \
                                             --method_labels "$method_labels" \
                                             --metrics  "$metrics"\
                                             --metric_labels "$metric_labels"  \
                                             --o "$outpfx".rds \
                                             --plt "$outpfx".pdf


### inweb
aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
methods="random,pcorr,glasso_likelihood,clr,aracne,et_genie3,mrnet,mrnetb,wgcna,spice"
method_labels="Random,PCor,GLasso,CLR,ARACNE,GENIE3,MRNET,MRNETB,WGCNA,SPICE"
metrics="inweb_aupr,inweb_spearman_r,inweb_precision_topNA_th0.01,inweb_hub_aupr,inweb_geodesic_rankdiff_2_1"
metric_labels="Interaction AUPR (InWeb_IM),Spearman r (InWeb_IM),Precision (InWeb_IM),Hub AUPR (InWeb_IM),Rankdiff 2-1 (InWeb_IM)"
outpfx="results/method_comparison_inweb"

Rscript src/plots/compare_multiple_methods.R --res "$aggregated_res_fn" \
                                             --methods "$methods" \
                                             --method_labels "$method_labels" \
                                             --metrics  "$metrics"\
                                             --metric_labels "$metric_labels"  \
                                             --o "$outpfx".rds \
                                             --plt "$outpfx".pdf

### string-experimental
aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
methods="random,pcorr,glasso_likelihood,clr,aracne,et_genie3,mrnet,mrnetb,wgcna,spice"
method_labels="Random,PCor,GLasso,CLR,ARACNE,GENIE3,MRNET,MRNETB,WGCNA,SPICE"
metrics="string_exp_aupr,string_spearman_r,string_exp_precision_topNA_th0.01,string_exp_hub_aupr,string_exp_geodesic_rankdiff_2_1"
metric_labels="Interaction AUPR (STRING-experiment),Spearman r (STRING-experiment),Precision (STRING-experiment),Hub AUPR (STRING-experimental),Rankdiff 2-1 (STRING-experiment)"
outpfx="results/method_comparison_string_exp"

Rscript src/plots/compare_multiple_methods.R --res "$aggregated_res_fn" \
                                             --methods "$methods" \
                                             --method_labels "$method_labels" \
                                             --metrics  "$metrics"\
                                             --metric_labels "$metric_labels"  \
                                             --o "$outpfx".rds \
                                             --plt "$outpfx".pdf



### string-exp + inweb + kegg
aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
methods="random,pcorr,glasso_likelihood,clr,aracne,et_genie3,mrnet,mrnetb,wgcna,spice"
method_labels="Random,PCor,GLasso,CLR,ARACNE,GENIE3,MRNET,MRNETB,WGCNA,SPICE"
metrics="string_exp_aupr,inweb_aupr,kegg_interaction_aupr,string_exp_precision_topNA_th0.01,inweb_precision_topNA_th0.01,kegg_interaction_precision_topNA_th0.01,string_exp_hub_aupr,inweb_hub_aupr,kegg_interaction_hub_aupr,string_exp_geodesic_rankdiff_2_1,inweb_geodesic_rankdiff_2_1,kegg_interaction_geodesic_rankdiff_2_1,string_spearman_r,inweb_spearman_r"
metric_labels="Interaction AUPR (Experiment),Interaction AUPR (InWeb_IM),Interaction AUPR (KEGG),Precision (Experiment),Precision (InWeb_IM),Precision (KEGG),Hub AUPR (Experiment),Hub AUPR (InWeb_IM),Hub AUPR (KEGG),Rankdiff 2-1 (Experiment),Rankdiff 2-1 (InWeb_IM),Rankdiff 2-1 (KEGG),Spearman r (Experiment),Spearman r (InWeb_IM)"
outpfx="results/method_comparison_string_exp_inweb_kegg"

Rscript src/plots/compare_multiple_methods.R --res "$aggregated_res_fn" \
                                             --methods "$methods" \
                                             --method_labels "$method_labels" \
                                             --metrics  "$metrics"\
                                             --metric_labels "$metric_labels"  \
                                             --o "$outpfx".rds \
                                             --plt "$outpfx".pdf
                                             

### kegg + reactome + go
aggregated_res_fn="results/all_evals/all_evaluations_test.rds"
methods="random,pcorr,glasso_likelihood,clr,aracne,et_genie3,mrnet,mrnetb,wgcna,spice"
method_labels="Random,PCor,GLasso,CLR,ARACNE,GENIE3,MRNET,MRNETB,WGCNA,SPICE"
metrics="kegg_shared_aupr,reactome_shared_aupr,go_shared_aupr,kegg_enriched_pathway,reactome_enriched_pathway,go_enriched_pathway"
metric_labels="Shared pathway AUPR (KEGG),Shared pathway AUPR (REACTOME),Shared pathway AUPR (GO),#Enriched pathways (KEGG),#Enriched pathways (REACTOME),#Enriched pathways (GO)"
outpfx="results/method_comparison_kegg_reactome_go"

Rscript src/plots/compare_multiple_methods.R --res "$aggregated_res_fn" \
                                             --methods "$methods" \
                                             --method_labels "$method_labels" \
                                             --metrics  "$metrics"\
                                             --metric_labels "$metric_labels"  \
                                             --o "$outpfx".rds \
                                             --plt "$outpfx".pdf
