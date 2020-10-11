get_genie3_net <- function(expr_mat, 
                           regulators = NULL,
                           targets = NULL,
                           n.regulators="sqrt", 
                           n.trees=1000, 
                           tree.method='RF', 
                           seed=NULL, 
                           n.cores=1,
                           verbose = F){
  require('GENIE3')
  
  # expression data
  stopifnot(class(expr_mat) %in% c('data.frame', 'matrix', 'ExpressionSet', 'RangedSummarizedExperiment') )
  if(class(expr_mat) == 'data.frame')
    expr_mat = as.matrix(expr_mat)
  
  # set seed
  if(!is.null(seed) && !is.na(seed) && is.finite(seed))
    set.seed(seed)
  
  # construct weight matrices
  weightMat <- GENIE3(expr_mat, regulators = regulators, targets = targets, treeMethod = tree.method, K = n.regulators, nTrees = n.trees, nCores=n.cores, verbose=verbose)
  
  return(weightMat)
}

get_aracne_net <- function(expr_df, mi.estimator = "spearman", eps=0, disc = "equalwidth", nbins = sqrt(ncol(expr_df))){
  require(minet)
  expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
  mim = build.mim(dataset =  expr_df, estimator = mi.estimator, disc = disc, nbins = nbins)
  aracne_net = aracne(mim, eps=eps)
  return(aracne_net)
}

get_clr_net <- function(expr_df, mi.estimator = "spearman", disc = "equalfreq", nbins = sqrt(ncol(expr_df))){
  require(minet)
  expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
  mim = build.mim(dataset =  expr_df, estimator = mi.estimator, disc = disc, nbins = nbins)
  clr_net = clr(mim)
  return(clr_net)
}

get_mrnet_net <- function(expr_df, mi.estimator = "spearman", disc = "equalfreq", nbins = sqrt(ncol(expr_df))){
  require(minet)
  expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
  mim = build.mim(dataset =  expr_df, estimator = mi.estimator, disc = disc, nbins = nbins)
  mrnet_net = mrnet(mim)
  return(mrnet_net)
}

get_mrnetb_net <- function(expr_df, mi.estimator = "spearman", disc = "equalfreq", nbins = sqrt(ncol(expr_df))){
  require(minet)
  expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
  mim = build.mim(dataset =  expr_df, estimator = mi.estimator, disc = disc, nbins = nbins)
  mrnetb_net = mrnetb(mim)
  return(mrnetb_net)
}

get_pcit_net <- function(expr_df, method = 'pearson'){
  require(PCIT)
  expr_df = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
  m = cor(expr_df, method = method)
  pcit_res = pcit(m)
  pcit_net = matrix(0, nrow = nrow(m), ncol = ncol(m), dimnames = list(rownames(m), colnames(m)))
  pcit_net[idx(pcit_res)] = 1  # keep only significant correlations
  diag(pcit_net) = 1
  return(pcit_net)
}

get_rlowpc_net <- function(expr_df){
  require(RLowPC)
  require(corpcor)
  require(matrixcalc)
  data.exp = as.data.frame(t(expr_df)) # convert to samples x genes matrix, required for build.mim()
  
  ##infer correlation network
  inf.cor<-abs(cor(data.exp))
  diag(inf.cor)<-0
  inf.pcor<-abs(pcor.shrink(data.exp, verbose=F)[1:ncol(data.exp),1:ncol(data.exp)])
  diag(inf.pcor)<-0
  
  # RLowPC
  reduction.sapce<-na.omit(adjmatrix2edgelist(adjmatrix = inf.pcor, directed = F,order = T))
  inf.RLowPC.edge<-RLowPC(data.exp = data.exp,edgelist = reduction.sapce, method = 'pearson', pc.estimator = 'shrink', progressbar = T)
  inf.RLowPC.edge$cor.weight<-abs(inf.RLowPC.edge$cor.weight)
  inf.RLowPC<-edgelist2adjmatrix(inf.RLowPC.edge[,c(1,2,4)],genes = colnames(data.exp), directed = T)
  inf.RLowPC<-abs(inf.RLowPC)
  inf.RLowPC<-pmax(inf.RLowPC,t(inf.RLowPC))
  return(inf.RLowPC)
}

get_random_net <- function(expr_df, seed = NULL, directed = F){
  if(!is.null(seed) && !is.na(seed) && is.finite(seed))
    set.seed(seed)
  
  net = matrix(NA, nrow = nrow(expr_df), ncol = nrow(expr_df), dimnames = list(rownames(expr_df), rownames(expr_df)))
  if(directed){
    net[1:nrow(net), 1:ncol(net)] = abs(rnorm(n = nrow(net) * ncol(net)))
  } else {
    net[lower.tri(net, diag = F)] = abs(rnorm(n = nrow(net) * (nrow(net)-1)/2))
    net = pmax(net, t(net), na.rm = T)
  }
  diag(net) = 0
  net = net / max(net)
  diag(net) = 1
  return(net)
}

get_glasso_net <- function(expr_df, lambda, inverse.normal = T, penalize.diagonal = F, thr=1.0e-4, maxit=1000, approx=FALSE){
  require(glasso)
  require(genomicsutil)
  
  if(inverse.normal == T)
    expr_df = to_inv_normal(as.matrix(expr_df))
  
  s = cov(t(expr_df))
  gl = glasso(s, rho = lambda, penalize.diagonal = penalize.diagonal, maxit = maxit, thr = thr, approx = approx)
  prec = gl$wi
  rownames(prec) = rownames(s)
  colnames(prec) = colnames(s)
  return(prec)
}

get_pcorr_net <- function(expr_df, method = "pearson", cov_df=NULL){
  # expr_df: gene x sample
  # cov_df: sample x cov
  require(ppcor)
  stopifnot(method %in% c("pearson", "kendall", "spearman"))
  expr_mat = t(expr_df) # convert to samples x genes matrix, required for ppcor
  genes = colnames(expr_mat)
  
  if(!is.null(cov_df)){
    stopifnot(all(rownames(expr_mat) %in% rownames(cov_df)))
    cov_df = cov_df[rownames(expr_mat),,drop=F]
    expr_mat = cbind(expr_mat, cov_df)
  }
  
  pcor_mat = ppcor::pcor(expr_mat, method = method)$estimate
  rownames(pcor_mat) = colnames(expr_mat)
  colnames(pcor_mat) = colnames(expr_mat)
  pcor_mat = pcor_mat[genes, genes, drop=F]  # exclude covariates
  return(pcor_mat)
}

get_aggregated_pcorr_net <- function(expr_df, cov_df=NULL, method='pearson', iter=100, frac.row = 1, frac.col= 0.8){
  require('miscutil')
  
  # pairwise_assoc_list = lapply(1:iter, function(it){
  #   verbose_print(sprintf('aggregating pcorr - iteration# %d', it))
  #   sampled_expr = sample_df(x = expr_df, 
  #                            size.row = frac.row * 100, 
  #                            size.col = frac.col * 100, 
  #                            unit = 'percent', 
  #                            replace.row = F, 
  #                            replace.col = F)
  #   sampled_assoc = get_pcorr_net(expr_df = sampled_expr, method = method)
  #   return(sampled_assoc)
  # })
  
  pairwise_assoc_list = list()
  it = 1
  while(it <= iter){
    verbose_print(sprintf('aggregating pcorr - iteration# %d', it))
    sampled_expr = sample_df(x = expr_df, 
                             size.row = frac.row * 100, 
                             size.col = frac.col * 100, 
                             unit = 'percent', 
                             replace.row = F, 
                             replace.col = F)
    sampled_assoc = tryCatch(get_pcorr_net(expr_df = sampled_expr, method = method), error = function(e) {return(NULL)})
    if(!is.null(sampled_assoc)){  # avoid if get_pcorr_net gets an error
      pairwise_assoc_list[[it]] = sampled_assoc
      it = it+1
    }
  }
  
  pairwise_assoc_mean = Reduce('+', pairwise_assoc_list) / length(pairwise_assoc_list)
  return(pairwise_assoc_mean)
}


get_glasso_net_optimized_by_likelihood <- function(expr_df, lambdas, fold, out.pfx, validation.frac = 0.2, seed=101, inverse.normal = T, penalize.diagonal = F, thr=1.0e-4, maxit=1000, approx=FALSE, save.net = T, use.old.net = T, verbose = F){
  require(glasso)
  require(genomicsutil)
  
  stopifnot(is.data.frame(expr_df) || is.matrix(expr_df))
  stopifnot(is.numeric(lambdas) && all(lambdas >= 0))
  stopifnot(is.numeric(fold) && length(fold)==1 && fold >= 1)
  stopifnot(dir.exists(dirname(out.pfx)))
  stopifnot(is.numeric(seed))
  
  lambdas = sort(lambdas, decreasing = T)
  get_theta_file <- function(cur_fold, cur_lambda){
    fn <- sprintf("%s_seed%s_fold%s_lambda%s.rds", out.pfx, seed, cur_fold, cur_lambda)
  }
  
  traceSTheta <- function(S, theta){
    sum(sapply(1:nrow(S), function(i) S[i, ] %*% theta[ ,i])) 
  }
  
  n_samples = ncol(expr_df)
  likelihood_mat <- sapply(1:fold, FUN = function(cur_fold){
    set.seed(seed + cur_fold)
    validation_samples = sort(sample(1:n_samples, size = round(n_samples * validation.frac), replace = F))
    train_samples = setdiff(1:n_samples, validation_samples)
    train_expr_df <- expr_df[ ,train_samples, drop = F]
    validation_expr_df <- expr_df[ ,validation_samples, drop = F]
    if(inverse.normal == T){
      train_expr_df = to_inv_normal(as.matrix(train_expr_df))
      validation_expr_df = to_inv_normal(as.matrix(validation_expr_df))
    }
      
    validation_S = cov(t(validation_expr_df))
    
    sapply(lambdas, function(lambda){
      verbose_print(sprintf('running glasso fold: %s, lambda: %s', cur_fold, lambda), verbose = verbose)
      theta_fn <- get_theta_file(cur_fold = cur_fold, cur_lambda = lambda)
      if(use.old.net && file.exists(theta_fn)){
        theta <- readRDS(theta_fn)
      } else {
        theta <- get_glasso_net(expr_df = train_expr_df, lambda = lambda, inverse.normal = F, penalize.diagonal = penalize.diagonal, thr = thr, maxit = maxit, approx = approx)
        if(save.net)
          saveRDS(theta, file = theta_fn)
      }
      
      validation_likelihood = as.numeric(determinant(theta, logarithm = T)$modulus) - traceSTheta(validation_S, theta)
      return(validation_likelihood)
    })
  })
  
  # select lambda
  avg_likelihood_per_lambda = rowMeans(likelihood_mat)
  max_likelihood = max(avg_likelihood_per_lambda)
  selected_lambda = lambdas[which.max(avg_likelihood_per_lambda)];
  
  # warning if selected lambda is the smallest lambda
  if(selected_lambda == min(lambdas))
    warning("the smallest lambda is the optimum.")
  
  # plot likelihood vs lambda
  plt_fn = sprintf("%s_validation_likelihood.pdf", out.pfx)
  pdf(plt_fn)
  plot(x=lambdas, y=avg_likelihood_per_lambda,  type = 'b', xlab = 'Lambda', ylab = 'Log Likelihood', main=sprintf("Optimum lambda: %s", selected_lambda), pch = 19)
  points(selected_lambda, max_likelihood, type = 'p', pch = 8)
  dev.off()
  
  # run graphical lasso with the selected lambda
  verbose_print(sprintf('reconstructing final glasso net with selected lambda: %s', selected_lambda), verbose = verbose)
  final_net_fn = sprintf("%s_lambda%s.rds", out.pfx, selected_lambda)
  if(use.old.net && file.exists(final_net_fn)){
    final_net = readRDS(final_net_fn)
  } else {
    final_net <- get_glasso_net(expr_df = expr_df, lambda = selected_lambda, inverse.normal = inverse.normal, penalize.diagonal = penalize.diagonal, thr = thr, maxit = maxit, approx = approx)
    if(save.net)
      saveRDS(final_net, file = final_net_fn)
  }
  
  return(final_net)
}

get_glasso_net_optimized_by_edge_count <- function(expr_df, 
                                                   lambda.max, lambda.min, lambda.interval, lambda.fine.interval,
                                                   expected.n.edge, out.pfx, 
                                                   inverse.normal = T, penalize.diagonal = F, thr=1.0e-4, maxit=1000, approx=FALSE, 
                                                   save.net = T, use.old.net = T, verbose = F){
  require(glasso)
  require(genomicsutil)
  require(ioutil)
  
  stopifnot(is.data.frame(expr_df) || is.matrix(expr_df))
  stopifnot(is.numeric(lambda.max) && is.numeric(lambda.min) && is.numeric(lambda.interval) && is.numeric(lambda.fine.interval))
  stopifnot(lambda.max > 0 && lambda.min >= 0 && lambda.interval > 0 && lambda.fine.interval > 0)
  stopifnot(dir.exists(dirname(out.pfx)))
  
  get_theta_file <- function(cur_lambda){
    fn <- sprintf("%s_lambda%s.rds", out.pfx, cur_lambda)
  }
  
  lambdas = c()
  edge_counts = c()
  
  best_net = NULL
  best_edge_count = -Inf
  best_lambda = Inf
  
  lambda = lambda.max
  loop_lambda_min = lambda.min
  loop_lambda_interval = lambda.interval
  
  while(lambda >= loop_lambda_min){
    verbose_print(sprintf('running glasso lambda: %s', lambda), verbose = verbose)
    theta_fn <- get_theta_file(cur_lambda = lambda)
    if(use.old.net && file.exists(theta_fn)){
      theta <- readRDS(theta_fn)
    } else {
      theta <- get_glasso_net(expr_df = expr_df, lambda = lambda, inverse.normal = inverse.normal, penalize.diagonal = penalize.diagonal, thr = thr, maxit = maxit, approx = approx)
      if(save.net)
        saveRDS(theta, file = theta_fn)
    }
    
    n_edge = sum(theta[lower.tri(theta, diag = F)] != 0)
    verbose_print(sprintf('edge count: %s', n_edge), verbose = verbose)
    
    lambdas = append(lambdas, lambda)
    edge_counts = append(edge_counts, n_edge)
    
    # select current best network
    if(abs(n_edge - expected.n.edge) < abs(best_edge_count - expected.n.edge)){
      best_lambda = lambda
      best_edge_count = n_edge
      best_net = theta
    }
    
    # break if more than expected edges found
    if(n_edge == expected.n.edge){
      break
    } else if(n_edge > expected.n.edge){
      if (loop_lambda_interval <= lambda.fine.interval)
        break
      loop_lambda_interval = lambda.fine.interval
      lambda = lambda + lambda.interval - lambda.fine.interval
      if (lambda >= lambda.max)
        break
    } else {
      lambda = lambda - loop_lambda_interval
    }
    
  }
  
  # warning if selected lambda is the smallest or the largest lambda
  if(best_lambda == lambda.min)
    warning("the smallest lambda is the optimum.")
  if(best_lambda == lambda.max)
    warning("the largest lambda is the optimum.")
  
  # plot edge_count vs lambda
  plt_fn = sprintf("%s_edge_count.pdf", out.pfx)
  pdf(plt_fn)
  plot(x=lambdas, y=edge_counts,  type = 'b', xlab = 'Lambda', ylab = 'Edge count', main=sprintf("Optimum lambda: %s", best_lambda), pch = 19)
  points(best_lambda, best_edge_count, type = 'p', pch = 8)
  dev.off()
  
  return(best_net)
}

# simulated_data = simulate_network_data(n.features = 10, n.edges = 10, n.samples = 100)
# expr_df = simulated_data$data
# orig_net = simulated_data$precision
# sim_net = simulated_data$pos_def_precision
# lambdas = c(1e-5, 1e-4, 1e-3, seq(0.01, 0.3, 0.01))
# fold = 3
# out.pfx = "results/glasso_test/sim"
# validation.frac = 0.2
# seed=101
# inverse.normal = T
# penalize.diagonal = F
# thr=1.0e-4
# maxit=1000
# approx=FALSE
# reconstructed_net = get_glasso_net_optimized_by_likelihood(expr_df = expr_df, lambdas = lambdas, out.pfx = out.pfx, fold = 5, use.old.net = F, save.net = T)
# 
# library(gplots)
# pdf(sprintf("%s_sim_reconstructed.pdf", out.pfx))
# heatmap.2(orig_net, Rowv = NULL, Colv = NULL, main = "Original Precision Matrix")
# heatmap.2(sim_net, Rowv = NULL, Colv = NULL, main = "Positive Definite Precision Matrix")
# heatmap.2(reconstructed_net, Rowv = NULL, Colv = NULL, main = "Reconstructed Precision Matrix")
# dev.off()


get_wgcna_net <- function(expr_mat, 
                          out.pfx = NULL,
                          method = 'pearson',
                          type = 'unsigned',
                          RsquaredCut = 0.8,
                          powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
                          removeFirst = FALSE, 
                          nBreaks = 10, 
                          blockSize = NULL, 
                          moreNetworkConcepts = FALSE,
                          verbose = 0){
  require('WGCNA')
  
  # expression data
  stopifnot(is.data.frame(expr_mat) || is.matrix(expr_mat))
  stopifnot(type %in% c('signed', 'unsigned'))
  expr_mat = t(expr_mat)  # sample x gene matrix
  
  ### compute similarity
  if(sum(!is.finite(expr_mat)) > 0){
    similarity_mat = cor(expr_mat, use = 'pairwise.complete.obs', method = method)
  } else {
    similarity_mat = cor(expr_mat, use = 'all.obs', method = method)
  }
  if(sum(is.na(similarity_mat)) > 0)
    similarity_mat[is.na(similarity_mat)] = 0
  if (type == 'unsigned'){
    similarity_mat = abs(similarity_mat)
  } else {
    similarity_mat = 0.5 * (1+similarity_mat)
  }
  
  ### estimate power
  threshold = pickSoftThreshold(data = similarity_mat, 
                                dataIsExpr = F,
                                RsquaredCut = RsquaredCut, 
                                powerVector = powerVector, 
                                removeFirst = removeFirst, 
                                nBreaks = nBreaks, 
                                blockSize = blockSize, 
                                networkType = type,
                                moreNetworkConcepts = moreNetworkConcepts, 
                                verbose = verbose)
  powerEstimate = threshold$powerEstimate
  scaleFreeFit = threshold$fitIndices
  
  if(is.na(powerEstimate)){
    warning(sprintf('could not pick a threshold for R^2 >= %s. picking a threshold with max R^2 = %s', RsquaredCut, max(scaleFreeFit[,2])))
    powerEstimate = scaleFreeFit$Power[which.max(scaleFreeFit$SFT.R.sq)]
  }
  
  net=similarity_mat^powerEstimate
  
  ### compute connectivity
  sf_block_size = ifelse(is.numeric(blockSize), blockSize, 5000)
  if(ncol(expr_mat) <= sf_block_size){
    connectivity=as.vector(apply(net,2,sum, na.rm=T)) - 1  # -1 for diagonals
  } else {
    connectivity=softConnectivity(datExpr = expr_mat, 
                                  corFnc = "cor", 
                                  corOptions = sprintf("use = 'p', method = '%s'", method),
                                  type = type,
                                  power = powerEstimate, 
                                  blockSize = sf_block_size, 
                                  verbose = verbose)
  }

  if(!is.null(out.pfx)){
    ### save data
    data_fn = sprintf("%s_wgcna_scale_free.RData", out.pfx)
    save(powerEstimate, scaleFreeFit, connectivity, file = data_fn)
    
    ### plot
    plt_fn = sprintf("%s_wgcna_scale_free.pdf", out.pfx)
    pdf(plt_fn)
    
    plot(scaleFreeFit[,1], -sign(scaleFreeFit[,3])*scaleFreeFit[,2],xlab="Soft Threshold (power)",ylab="R^2",type="b", pch = 20)
    points(scaleFreeFit[,1], -sign(scaleFreeFit[,3])*scaleFreeFit[,2], type = "p", pch = 19, col="white")
    text(scaleFreeFit[,1], -sign(scaleFreeFit[,3])*scaleFreeFit[,2], labels=as.character(scaleFreeFit[,1]))
    abline(h=RsquaredCut,col="red")
    plot(scaleFreeFit[,1], scaleFreeFit[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="b", pch = 20)
    points(scaleFreeFit[,1], scaleFreeFit[,5], type = "p", pch = 19, col="white")
    text(scaleFreeFit[,1], scaleFreeFit[,5], labels = as.character(scaleFreeFit[,1]))
    
    hist(connectivity, main = sprintf("Connectivity with power=%s", powerEstimate), xlab = "Connectivity", breaks = nBreaks)
    scaleFreePlot(connectivity, main= sprintf("Scale-free topology: power=%s\n", powerEstimate))
    dev.off()  
  }

  return(net)
}


get_spqn_wgcna_net <- function( expr_mat,
                                mean_count_per_gene,
                                out.pfx,
                                ngrp = 10, 
                                size_grp = round(nrow(expr_mat)/ngrp*1.5), 
                                ref_grp = ngrp-1,
                                type = 'unsigned',
                                RsquaredCut = 0.8,
                                powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
                                removeFirst = FALSE, 
                                nBreaks = 10, 
                                blockSize = NULL, 
                                moreNetworkConcepts = FALSE,
                                verbose = 0){
  require('WGCNA')
  require('spqn')
  
  # expression data
  stopifnot(class(expr_mat) %in% c('data.frame', 'matrix') )
  stopifnot(type %in% c('signed', 'unsigned'))
  stopifnot(nrow(expr_mat) == length(mean_count_per_gene) && is.numeric(mean_count_per_gene))
  stopifnot(ref_grp>=1 && ref_grp <= ngrp)
  stopifnot(size_grp <= nrow(expr_mat))
  expr_mat = t(expr_mat)  # sample x gene matrix
  
  ### compute similarity
  if(sum(!is.finite(expr_mat)) > 0){
    similarity_mat = cor(expr_mat, use = 'pairwise.complete.obs', method = 'pearson')
  } else {
    similarity_mat = cor(expr_mat, use = 'all.obs', method = 'pearson')
  }
  
  similarity_mat = normalize_correlation(cor_mat = similarity_mat, 
                                         ave_logrpkm = mean_count_per_gene, 
                                         ngrp = ngrp, 
                                         size_grp = size_grp, 
                                         ref_grp = ref_grp)
  rownames(similarity_mat) = colnames(expr_mat)
  colnames(similarity_mat) = colnames(expr_mat)
  
  if(sum(is.na(similarity_mat)) > 0)
    similarity_mat[is.na(similarity_mat)] = 0
  if (type == 'unsigned')
    similarity_mat = abs(similarity_mat)
  
  ### estimate power
  threshold = pickSoftThreshold(data = similarity_mat, 
                                dataIsExpr = F,
                                RsquaredCut = RsquaredCut, 
                                powerVector = powerVector, 
                                removeFirst = removeFirst, 
                                nBreaks = nBreaks, 
                                blockSize = blockSize, 
                                networkType = type,
                                moreNetworkConcepts = moreNetworkConcepts, 
                                verbose = verbose)
  powerEstimate = threshold$powerEstimate
  scaleFreeFit = threshold$fitIndices
  
  if(is.na(powerEstimate)){
    # save data before returning error
    data_fn = sprintf("%s_spqn_wgcna_scale_free.RData", out.pfx)
    save(powerEstimate, scaleFreeFit, file = data_fn)
    stop(sprintf('could not pick a threshold for R^2 >= %s. max R^2 is %s', RsquaredCut, max(scaleFreeFit[,2])))
  }
  
  ### compute connectivity
  if(type == 'signed'){
    net = (0.5 * (1+similarity_mat) )^powerEstimate
  } else {
    net=similarity_mat^powerEstimate
  }
  
  if(ncol(expr_mat) < 5000){
    connectivity=as.vector(apply(net,2,sum, na.rm=T))  # -1 for diagonals
  } else {
    sf_block_size = ifelse(is.numeric(blockSize), blockSize, 1500)
    connectivity=softConnectivity(datExpr = expr_mat, corFnc = "cor", 
                                  type = type,
                                  power = powerEstimate, 
                                  blockSize = sf_block_size, 
                                  verbose = verbose)
  }
  
  ### save data
  data_fn = sprintf("%s_spqn_wgcna_scale_free.RData", out.pfx)
  save(powerEstimate, scaleFreeFit, connectivity, file = data_fn)
  
  ### plot
  plt_fn = sprintf("%s_spqn_wgcna_scale_free.pdf", out.pfx)
  pdf(plt_fn)
  
  plot(scaleFreeFit[,1], -sign(scaleFreeFit[,3])*scaleFreeFit[,2],xlab="Soft Threshold (power)",ylab="R^2",type="b", pch = 20)
  points(scaleFreeFit[,1], -sign(scaleFreeFit[,3])*scaleFreeFit[,2], type = "p", pch = 19, col="white")
  text(scaleFreeFit[,1], -sign(scaleFreeFit[,3])*scaleFreeFit[,2], labels=as.character(scaleFreeFit[,1]))
  abline(h=RsquaredCut,col="red")
  plot(scaleFreeFit[,1], scaleFreeFit[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="b", pch = 20)
  points(scaleFreeFit[,1], scaleFreeFit[,5], type = "p", pch = 19, col="white")
  text(scaleFreeFit[,1], scaleFreeFit[,5], labels = as.character(scaleFreeFit[,1]))
  
  hist(connectivity, main = sprintf("Connectivity with power=%s", powerEstimate), xlab = "Connectivity", breaks = nBreaks)
  scaleFreePlot(connectivity, main= sprintf("Scale-free topology: power=%s\n", powerEstimate))
  dev.off()
  
  return(net)
}

get_cor_net <- function(expr_mat, 
                        method = 'pearson',
                        verbose = 0){
  # expression data
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  expr_mat = t(expr_mat)  # sample x gene matrix
  
  ### compute similarity
  if(sum(!is.finite(expr_mat)) > 0){
    similarity_mat = cor(expr_mat, use = 'pairwise.complete.obs', method = method)
  } else {
    similarity_mat = cor(expr_mat, use = 'all.obs', method = method)
  }
  if(sum(is.na(similarity_mat)) > 0)
    similarity_mat[is.na(similarity_mat)] = 0
  
  return(similarity_mat)
}
