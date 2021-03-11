library(ioutil)
library(miscutil)
library(spice)
library(huge)
library(corrplot)
library(doRNG) # for GENIE3

### settings for simulation
sim_dir = "results/debug/simulation_200samples_1000genes_random_p0.01"
sim_seed = 100
n_repeat = 5
edge_prob = 0.01
sim_samples = rep(200, n_repeat)
sim_genes = rep(1000, n_repeat)
sim_types = rep("random", n_repeat)
sim_labels = paste0("random", seq_len(n_repeat))

stopifnot(length(sim_genes) == length(sim_types))
stopifnot(length(sim_genes) == length(sim_samples))
stopifnot(length(sim_genes) == length(sim_labels))

### settings for running networks
methods = c( "random", "pcorr", "glasso_likelihood", "clr", "aracne",
             "et_genie3",
             "mrnet", "mrnetb", "wgcna", "wgcna_pearson_unsigned_0.8", "spice" )
# #methods = c( "wgcna_pearson_unsigned_0.8", "spice" )
# methods = c("et_genie3")
eval_seed = 101
n_cores = 3

# ### settings for debugging
# methods = methods[1:2]
# method_labels = method_labels[1:2]
# sim_labels = sim_labels[1:2]

### function to simulate data
generate_simulated_data <- function(nsample,
                                    ngene,
                                    edge_prob = 0.01,
                                    graph_type = "scale-free") {
  sim = huge::huge.generator(
    n = nsample,
    d = ngene,
    prob = edge_prob,
    graph = graph_type,
    verbose = F
  )
  
  ### add gene and sample names to the data
  genes = paste0("G", seq_len(ngene))
  samples = paste0("S", seq_len(nsample))
  
  # data should be gene x sample
  sim$data = t(sim$data)
  rownames(sim$data) = genes
  colnames(sim$data) = samples
  
  # adjacency
  sim$theta = as.matrix(sim$theta)
  rownames(sim$theta) = colnames(sim$theta) = genes
  
  # precision
  rownames(sim$omega) = colnames(sim$omega) = genes
  
  return(list(
    data = sim$data,
    adjacency = sim$theta,
    precision = sim$omega
  ))
}

### simulate data and save -- one dataset
simulate_and_save <- function(nsample, ngene, graph_type, dir, dataset_name, edge_prob){
  sim = generate_simulated_data(nsample = nsample,
                                ngene = ngene,
                                graph_type = graph_type,
                                edge_prob = edge_prob)
  # process precision to use as a known matrix
  known = abs(sim$precision)
  known = pmax(known, t(known), na.rm = T)
  diag(known) = 0
  known = known / max(known)
  
  # save
  if(!dir.exists(dir)){
    dir.create(dir)
  }
  
  expr_fn = sprintf("%s/%s_expression.txt", dir, dataset_name)
  adj_fn = sprintf("%s/%s_adjacency.rds", dir, dataset_name)
  prec_fn = sprintf("%s/%s_precision.rds", dir, dataset_name)
  known_fn = sprintf("%s/%s_known.rds", dir, dataset_name)
  
  write_df(sim$data, file = expr_fn)
  saveRDS(sim$adjacency, file = adj_fn)
  saveRDS(sim$precision, file = prec_fn)
  saveRDS(known, file = known_fn)
  
  return(invisible(NULL))
}

# simulate data and save -- all datasets
simulate_and_save_all <- function(sim_labels,
                                  sim_samples,
                                  sim_genes,
                                  sim_types,
                                  sim_seed,
                                  sim_dir, 
                                  edge_prob) {
  set.seed(sim_seed)
  for (sidx in seq_along(sim_labels)) {
    verbose_print(sprintf("simulating %s ...", sim_labels[sidx]))
    simulate_and_save(
      nsample = sim_samples[sidx],
      ngene = sim_genes[sidx],
      graph_type = sim_types[sidx],
      dir = sprintf("%s/%s", sim_dir, sim_labels[sidx]),
      dataset_name = sim_labels[sidx],
      edge_prob = edge_prob
    )
  }
  return(invisible(NULL))
}

simulate_and_save_all(sim_labels = sim_labels,
                      sim_samples = sim_samples,
                      sim_genes = sim_genes,
                      sim_types = sim_types,
                      sim_seed = sim_seed,
                      sim_dir = sim_dir, 
                      edge_prob = edge_prob)


##############################################
############## evaluate networks #############
##############################################

run_network <- function(method, sim_dir, sim_label, eval_seed = 101, n_cores = 1){
  expr_fn = sprintf("%s/%s/%s_expression.txt", sim_dir, sim_label, sim_label)
  out_dir = sprintf("%s/%s", sim_dir, sim_label)
  stopifnot(file.exists(expr_fn))
  
  cmd = sprintf("Rscript src/run_network.R \\
                --expr '%s' --method '%s' --seed %s --o '%s' --thread %s", 
                expr_fn, method, eval_seed, out_dir, n_cores)
  status = system(cmd, intern = F)
  if(status != 0){
    stop(sprintf("error in running network with the following parameters:\n%s.", 
                 paste(c(method, sim_dir, sim_label, eval_seed), collapse = ",\n")) )
  }
  return(invisible(status))
}

run_network_all <- function(methods,
                            sim_labels,
                            sim_dir = sim_dir,
                            eval_seed = eval_seed, 
                            n_cores = n_cores) {
  for(sim_label in sim_labels){
    for(method in methods){
      verbose_print(sprintf("running %s on %s ...", method, sim_label))
      run_network(
        method = method,
        sim_dir = sim_dir,
        sim_label = sim_label,
        eval_seed = eval_seed,
        n_cores = n_cores
      )
    }
  }
  
}

eval_network <- function(method,
                         sim_label,
                         sim_dir = sim_dir) {
  res_dir = sprintf("%s/%s", sim_dir, sim_label)
  net_fn = sprintf("%s/%s_absnet.rds", res_dir, method)
  known_fn = sprintf("%s/%s_known.rds", res_dir, sim_label)
  adj_fn = sprintf("%s/%s_adjacency.rds", res_dir, sim_label)
  
  stopifnot(file.exists(net_fn))
  stopifnot(file.exists(known_fn))
  stopifnot(file.exists(adj_fn))
  
  ### interaction AUPR
  known_auc_fn = sprintf("%s/%s_string_ppi_auc.rds", res_dir, method)
  adj_auc_fn = sprintf("%s/%s_kegg_interaction_auc.rds", res_dir, method)
  
  cmd = sprintf("Rscript src/eval_interaction_auc.R \\
                --net '%s' --known '%s' --o '%s' \\
                --curve TRUE --max TRUE --min TRUE --rand FALSE --dg FALSE \\
                --na 'known' --neg 'error'", 
                net_fn, known_fn, known_auc_fn)
  system(cmd, intern = F)
  
  cmd = sprintf("Rscript src/eval_interaction_auc.R \\
                --net '%s' --known '%s' --o '%s' \\
                --curve TRUE --max TRUE --min TRUE --rand FALSE --dg FALSE \\
                --na 'known' --neg 'error'", 
                net_fn, adj_fn, adj_auc_fn)
  system(cmd, intern = F)
  
  ### spearman r
  known_spearman_fn = sprintf("%s/%s_string_ppi_spearman_cor.rds", res_dir, method)
  adj_spearman_fn = sprintf("%s/%s_kegg_interaction_spearman_cor.rds", res_dir, method)
  
  cmd = sprintf("Rscript src/eval_interaction_cor.R \\
                --net '%s' --known '%s' --o '%s' \\
                --method 'spearman' --na 'known' --neg 'error'", 
                net_fn, known_fn, known_spearman_fn)
  system(cmd, intern = F)
  
  cmd = sprintf("Rscript src/eval_interaction_cor.R \\
                --net '%s' --known '%s' --o '%s' \\
                --method 'spearman' --na 'known' --neg 'error'", 
                net_fn, adj_fn, adj_spearman_fn)
  system(cmd, intern = F)
  
  ### precision
  adj_prec_fn = sprintf("%s/%s_kegg_interaction_precision.rds", res_dir, method)
  
  cmd = sprintf("Rscript src/eval_interaction_precision_at_top.R \\
                --net '%s' --known '%s' --o '%s' \\
                --threshold '0.01' --top 'NA' --na 'known' --neg 'error'", 
                net_fn, adj_fn, adj_prec_fn)
  system(cmd, intern = F)
  
  ### hub aupr
  known_hub_auc_fn = sprintf("%s/%s_string_ppi_hub_auc.rds", res_dir, method)
  adj_hub_auc_fn = sprintf("%s/%s_kegg_interaction_hub_auc.rds", res_dir, method)
  
  cmd = sprintf("Rscript src/eval_interaction_hub_auc.R \\
                --net '%s' --known '%s' --o '%s' \\
                --curve TRUE --max TRUE --min TRUE --rand FALSE --dg FALSE \\
                --na 'known' --neg 'error'", 
                net_fn, known_fn, known_hub_auc_fn)
  system(cmd, intern = F)
  
  cmd = sprintf("Rscript src/eval_interaction_hub_auc.R \\
                --net '%s' --known '%s' --o '%s' \\
                --curve TRUE --max TRUE --min TRUE --rand FALSE --dg FALSE \\
                --na 'known' --neg 'error'", 
                net_fn, adj_fn, adj_hub_auc_fn)
  system(cmd, intern = F)
  
  ### geodesic rankdiff
  adj_geodesic_rankdiff_fn = sprintf("%s/%s_kegg_interaction_geodesic_rankdiff.rds", res_dir, method)
  
  cmd = sprintf("Rscript src/eval_geodesic_rankdiff.R \\
                --net '%s' --known '%s' --o '%s' \\
                --threshold 0.01 -d '1,2,3,4,5'",
                net_fn, adj_fn, adj_geodesic_rankdiff_fn)
  system(cmd, intern = F)
  
  return(invisible(NULL))
}

run_and_eval_network_all <- function(methods,
                                     sim_labels,
                                     sim_dir = sim_dir,
                                     eval_seed = eval_seed, 
                                     n_cores = n_cores) {
  for (sim_label in sim_labels) {
    for (method in methods) {
      verbose_print(sprintf("running %s on %s ...", method, sim_label))
      run_network(method = method,
                  sim_label = sim_label,
                  sim_dir = sim_dir,
                  eval_seed = eval_seed, 
                  n_cores = n_cores)
      
      verbose_print(sprintf("evaluating %s on %s ...", method, sim_label))
      eval_network(method = method,
                   sim_label = sim_label,
                   sim_dir = sim_dir)
    }
  }
}

run_and_eval_network_all(methods = methods,
                         sim_labels = sim_labels,
                         sim_dir = sim_dir,
                         eval_seed = eval_seed, 
                         n_cores = n_cores)

##############################################
############## aggregate evaluations #########
##############################################

aggregate_evaluations <- function(method,
                                  sim_label,
                                  sim_dir = sim_dir) {
  
  res_dir = sprintf("%s/%s", sim_dir, sim_label)
  methods_str = paste(methods, sep = ',', collapse = ',')
  agg_fn = sprintf("%s/aggregated_evaluations_simulation.rds", res_dir)
  stopifnot(dir.exists(res_dir))

  cmd = sprintf("Rscript src/aggregate_evaluations_per_dataset.R \\
                --dir '%s' --methods '%s' --o '%s'",
                res_dir, methods_str, agg_fn)
  system(cmd)
}

aggregate_aggregated_results <- function(sim_labels = sim_labels,
                                         sim_dir = sim_dir) {
  ### aggregate
  aggregated_df = NULL
  for(tissue in sim_labels){
    fn = sprintf("%s/%s/aggregated_evaluations_simulation.rds", 
                 sim_dir, 
                 tissue)
    if(!file.exists(fn)){
      next
    }
    cur_df = readRDS(fn)
    if(is.null(aggregated_df)){
      # initialize aggregated_df
      aggregated_df = data.frame(tissue = tissue, 
                                 stringsAsFactors = F)
      aggregated_df = aggregated_df[rep(1, nrow(cur_df)), , drop = F]
      for(new_col in colnames(cur_df)){
        aggregated_df[,new_col] = cur_df[,new_col]
      }
    } else {
      # for each new column, put NA for previous methods
      for(new_col in setdiff(colnames(cur_df), colnames(aggregated_df))){
        aggregated_df[,new_col] = NA
      }
      # add current results with NA for unavailable columns
      cur_df2 = aggregated_df[rep(1, nrow(cur_df)),,drop=F]
      cur_df2[,] = NA
      cur_df2[,"tissue"] = tissue
      cur_df2[,colnames(cur_df)] = cur_df
      aggregated_df = rbind(aggregated_df, cur_df2)
      rownames(aggregated_df) = NULL
    }
  }
  return(aggregated_df)
}

aggregate_evaluations_all <- function(methods,
                                      sim_labels,
                                      agg_fn,
                                      sim_dir = sim_dir) {
  for (sim_label in sim_labels) {
    for (method in methods) {
      verbose_print(sprintf("aggregating %s on %s ...", method, sim_label))
      aggregate_evaluations(method = method,
                            sim_label = sim_label,
                            sim_dir = sim_dir)
    }
  }
  
  verbose_print("aggregating all evaluation metrics ...")
  all_sim_evals = aggregate_aggregated_results(sim_labels = sim_labels, sim_dir = sim_dir)
  saveRDS(all_sim_evals, agg_fn)
}

all_sim_eval_fn = sprintf("%s/all_evaluations_simulation.rds", sim_dir)
aggregate_evaluations_all(
  methods = methods,
  sim_labels = sim_labels,
  agg_fn = all_sim_eval_fn,
  sim_dir = sim_dir
)

##############################################
############## plots ######### ######### #####
##############################################

plot_comparision_multiple_methods <- function(all_sim_eval_fn, methods, metrics, plt_data_fn, plt_fn){
  ### plot all variables
  methods_str = paste(names(methods), sep = ',', collapse = ',')
  method_labels_str = paste(as.character(methods), sep = ',', collapse = ',')
  metrics_str = paste(names(metrics), sep = ',', collapse = ',')
  metric_labels_str = paste(as.character(metrics), sep = ',', collapse = ',')
  
  cmd = sprintf( "Rscript src/plots/compare_multiple_methods.R \\
               --res '%s' --methods '%s' --method_labels '%s' \\
               --metrics '%s' --metric_labels '%s' \\
               --o '%s' --plt '%s'", 
                 all_sim_eval_fn, methods_str, method_labels_str,
                 metrics_str, metric_labels_str, 
                 plt_data_fn, plt_fn)
  system(cmd)
  
}



### plot all variables
methods_to_plot = c(
  random = "Random",
  pcorr = "PCor",
  glasso_likelihood = "GLasso",
  clr = "CLR",
  aracne = "ARACNE",
  et_genie3 = "GENIE3",
  mrnet = "MRNET",
  mrnetb = "MRNETB",
  #wgcna = "WGCNA",
  wgcna_pearson_unsigned_0.8 = "WGCNA",
  spice = "SPICE"
)

adjacency_based_metrics = c(
  kegg_interaction_aupr = "Interaction AUPR (Simulation)",
  kegg_interaction_precision_topNA_th0.01 = "Precision (Simulation)",
  kegg_interaction_hub_aupr = "Hub AUPR (Simulation)",
  kegg_interaction_geodesic_rankdiff_2_1 = "Rankdiff 2-1 (Simulation)"
)

precision_based_metrics = c(
  string_aupr = "Interaction AUPR (Simulation)",
  string_spearman_r = "Spearman r (Simulation)",
  string_hub_aupr = "Hub AUPR (Simulation)"
)

plt_fn = sprintf("%s/simulation_method_comparison_plot_adjacency_based.pdf", sim_dir)
plt_data_fn = sprintf("%s/simulation_method_comparison_plot_adjacency_based.rds", sim_dir)
plot_comparision_multiple_methods(
  all_sim_eval_fn = all_sim_eval_fn,
  methods = methods_to_plot,
  metrics = adjacency_based_metrics,
  plt_fn = plt_fn,
  plt_data_fn = plt_data_fn
)

plt_fn = sprintf("%s/simulation_method_comparison_plot_precision_based.pdf", sim_dir)
plt_data_fn = sprintf("%s/simulation_method_comparison_plot_precision_based.rds", sim_dir)
plot_comparision_multiple_methods(
  all_sim_eval_fn = all_sim_eval_fn,
  methods = methods_to_plot,
  metrics = precision_based_metrics,
  plt_fn = plt_fn,
  plt_data_fn = plt_data_fn
)

