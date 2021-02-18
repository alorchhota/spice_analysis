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

Rscript src/plots/compare_methods.R \
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
metrics="kegg_interaction_aupr,kegg_interaction_spearman_r,kegg_interaction_precision_topNA_th0.7,kegg_interaction_hub_aupr,kegg_shared_aupr,kegg_enriched_pathway"
labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways"
# metrics="kegg_interaction_aupr,kegg_interaction_spearman_r,kegg_interaction_precision_topNA_th0.7,kegg_interaction_hub_aupr,kegg_shared_aupr,kegg_enriched_pathway,kegg_interaction_geodesic_rankdiff_2_1,kegg_interaction_geodesic_rankdiff_3_2"
# labels="Interaction,Spearman,Precision,Hub,Shared,#Pathways,Rank:2-1,Rank:3-2"
outpfx="results/${method1}_vs_${method2}_kegg"

Rscript src/plots/compare_methods.R \
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

Rscript src/plots/compare_methods.R \
  --res "$aggregated_res_fn" \
  --method1 "$method1" \
  --method2 "$method2" \
  --metrics "$metrics" \
  --labels "$labels" \
  --nsample "$nsample_fn" \
  --o "$outpfx.rds" \
  --plt "$outpfx.pdf"
  