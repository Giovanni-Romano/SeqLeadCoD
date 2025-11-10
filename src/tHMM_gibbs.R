tHMM_gibbs_adaptive_PT = function(
    Y,
    m,
    # Parameters for the likelihood-specific local parameters
    par_likelihood = list(u = NULL, # First hyperparam HIG
                          v = NULL # Second hyperparam HIG
    ),
    # Parameters for the tRPM process prior
    par_tRPM = list(a_alpha = 1, b_alpha = 1, # Beta prior on alpha
                    urn_type = c("Gnedin", "DP"),
                    eta = 0.5 # Gnedin/DP parameter
    ),
    # Control parameters for MCMC settings
    ctr_mcmc = list(seed = 1234, nburnin = 1000, nchain = 5000, print_step = 100,
                    ncl.init = NULL, verbose = "0",
                    parallel = FALSE, nthreads = 1, cluster = "PSOCK"),
    # Control parameters for swap moves
    ctr_swap  = list(nrep = 10, swap_step = 1, ntry = 1, normalize = FALSE, 
                     start = 50, geom_step = 0.9, exp_schedule = -2,
                     rate_gamma = 0.9, rate_delta = 0.1, rate_kappa = 0.75,
                     adaptive = TRUE, deterministic = FALSE),
    # Control parameters for result storing
    ctr_save = list(save = FALSE, filepath = "",
                    filename = paste("res_", Sys.Date(), ".RDS", sep = "")),
    print_step = 100,
    verbose = c("0", "1", "2"), # 0 = no output, 1 = some info, 2 = detailed info
    seed = 4238,
    save = TRUE,
    path_save = "",
    name_save = paste(paste("res", Sys.Date(), sep = "_"), ".RDS", sep = "")
){
  
  if (ctr_save$save) {
    if (!dir.exists(ctr_save$filepath)) {
      warning("The provided path_save does not exist. Creating it. \n")
      dir.create(ctr_save$filepath, recursive = T)
    }
    cat("##############################\n", "Saving out in RDS at path: ", paste0(ctr_save$filepath, ctr_save$filename, sep = ""), "\n", sep = "")
  }
  
  # Set and check the control parameters
  ctr_mcmc = do.call("set_ctr_mcmc", ctr_mcmc)
  ctr_swap = do.call("set_ctr_swap", ctr_swap)
  ctr_save = do.call("set_ctr_save", ctr_save)
  
  verbose = ctr_mcmc$verbose
  print_step = ctr_mcmc$print_step
  ncl_init = ctr_mcmc$ncl.init
  
  # Set the simulation seed
  set.seed(ctr_mcmc$seed)
  
  # Set the CPU cluster for the parallel execution of local updates
  if(ctr_mcmc$parallel){
    suppressPackageStartupMessages({
      library(parallel)
      library(foreach)
      library(doParallel)
      library(iterators)
    })
    
    clust <- parallel::makeCluster(ctr_mcmc$nthreads, type = ctr_mcmc$cluster)
    doParallel::registerDoParallel(clust)
    cat("Cluster ", attr(clust[[1]], "class")[1], " with ", length(clust), " nodes registered. \n", sep = "")
    funlist = c("tSBM_local_update", "tSBM_chain_swap", 
                "update_gamma", "update_label",
                "update_alpha", "create_suffstat_loglik",
                "update_suffstat",
                "eval_norm_const", "compatibility_check_cppwrapper",
                "sample_gumbel_max")
    
    parallel::clusterExport(clust, varlist = funlist)
  }
  
  # Set the number of iterations and swaps
  niter = ctr_mcmc$nburnin + ctr_mcmc$nchain
  nswap = floor(niter / ctr_swap$swap_step)
  
  # Set the data dimensions
  dims = dim(Y)
  n = dims[1]
  p = dims[2]
  .T = dims[3]
  m.inner = apply(Y, 2, function(x) length(unique(c(x)))) %>% unname()
  if (any(m != m.inner)){stop("Provided m is not consistent with Y!")}
  attrlist = lapply(m, function(x) 1:x)
  
  # Set the urn
  urn_type = par_tRPM$urn_type
  urn = switch(urn_type,
               "Gnedin" = urn_GN,
               "DP" = urn_DP,
               stop("Invalid urn_type"))
  
  # Set Parallel Tempering
  # Initialize the adaptive temperature parameters
  rho.tmp = rep(ctr_swap$exp_schedule, ctr_swap$nrep - 1)
  rho = matrix(NA, nrow = ctr_swap$nrep - 1, ncol = nswap)
  rho[ , 1] = rho.tmp
  prob_move = matrix(NA, nrow = ctr_swap$nrep - 1, ncol = nswap)
  prob_acc = matrix(NA, nrow = floor(ctr_swap$nrep/2 - 1/2), ncol = nswap)
  
  if (!(ctr_swap$adaptive)){
    temperatures.tmp = init_temp_schedule(method = "custom", custom = ctr_swap$custom, normalize = FALSE, decreasing = TRUE)
  } else {
    temperatures.tmp = init_temp_schedule(method = "exp", n_replica = ctr_swap$nrep,
                                          exp = rho.tmp[1], normalize = TRUE, decreasing = TRUE)
  }
  tempering_schedule = matrix(NA, nrow = ctr_swap$nrep, ncol = niter)
  tempering_schedule[ , 1] = temperatures.tmp
  
  chains_tracking = matrix(NA, nrow = ctr_swap$nrep, ncol = niter)
  chains_tracking[ , 1] = 1:ctr_swap$nrep
  
  loglik_tracking = matrix(NA, nrow = ctr_swap$nrep, ncol = niter)
  
  # OBJECTS TO STORE RESULTS
  traces = lapply(1:ctr_swap$nrep, function(j) {
    initialized_objects = tHMM_gibbs_initialize(Y = Y, n = n, .T = .T, niter = niter,
                                                ncl.init = ncl_init, p = p)
    
    initialized_objects[["temperature"]] = temperatures.tmp[j]
    initialized_objects
  })
  
  traces_tmp = lapply(traces,
                      function(x) list(C.tmp = x$C.tmp, mu.tmp = x$mu.tmp, 
                                       alpha.tmp = x$alpha.tmp, gamma.tmp = x$gamma.tmp, 
                                       clS.tmp = x$clS.tmp, matches.tmp = x$matches.tmp,
                                       log_lik.tmp = x$log_lik.tmp,
                                       cached_norm_const = x$cached_norm_const, 
                                       denomin_log_lik = x$denomin_log_lik))
  
  # Elapsed time initialization
  time_start = time_last = Sys.time()
  cat("Start Gibbs of ", niter, " iterations", sep = "")
  
  
  attempted_swap = rep(0, ctr_swap$nrep)
  succesful_swap = rep(0, ctr_swap$nrep)
  
  # Counter for intermediate saves
  counter_save = 0L
  # How many iter between intermediate saves
  iter_to_save = 5000
  
  # Gibbs sampling loop
  for (d in 2:niter){
    
    if (ctr_mcmc$parallel){
      traces_tmp = foreach(tr = iterators::iter(traces_tmp)) %dopar% {
        tHMM_local_update(d = d, 
                          Y = Y, m = m, 
                          trace = tr,
                          par_likelihood = par_likelihood,
                          par_tRPM = par_tRPM,
                          print_step = ctr_mcmc$print_step,
                          verbose = ctr_mcmc$verbose)
      }
    } else {
      
      traces_tmp = lapply(traces_tmp, function(tr)
        tHMM_local_update(d = d, 
                          Y = Y, m = m, 
                          trace = tr,
                          par_likelihood = par_likelihood,
                          par_tRPM = par_tRPM,
                          print_step = ctr_mcmc$print_step,
                          verbose = ctr_mcmc$verbose)
      )
      
    }
    
    
    chains_tracking[ , d] = chains_tracking[ , d-1]
    
    if (d %% ctr_swap$swap_step == 0){
      out_swap = tHMM_chain_swap(traces_tmp, d, 
                                 temperatures = temperatures.tmp, 
                                 deterministic = ctr_swap$deterministic,
                                 n_try = ctr_swap$ntry,
                                 family = family, 
                                 prior = pr_family)
      
      traces_tmp = out_swap$traces
      
      # Store tmp results in the final objects
      for (j in seq_along(traces)){
        # Partitions
        for (t in 1:.T){
          traces[[j]]$C[ , t, d] = mat2vec(traces_tmp[[j]]$C.tmp[[t]])
        }
        # gamma
        traces[[j]]$gamma[ , d] = as.integer(colSums(traces_tmp[[j]]$gamma.tmp[ , 1:.T]))
        # alpha
        traces[[j]]$alpha[ , d] = traces_tmp[[j]]$alpha.tmp
        # loglik
        traces[[j]]$log_lik[d] = traces_tmp[[j]]$log_lik.tmp
      }
      
      attempted_swap = attempted_swap + out_swap$diagnostics$att
      succesful_swap = succesful_swap + out_swap$diagnostics$succ
      
      forward_idx = which(out_swap$diagnostics$succ == 1)
      backward_idx = forward_idx + 1
      which_forward = chains_tracking[, d-1] %in% forward_idx
      which_backward = chains_tracking[, d-1] %in% backward_idx
      
      chains_tracking[which_forward, d] = chains_tracking[which_forward, d-1] + 1
      chains_tracking[which_backward, d] = chains_tracking[which_backward, d-1] - 1
      loglik_tracking[ , d] = out_swap$loglik
      
      # Temperature update
      diff_loglik.tmp = diff(out_swap$loglik)
      diff_temp.tmp = diff(temperatures.tmp)
      prob_move.tmp = pmin(1, exp(- diff_temp.tmp * diff_loglik.tmp))
      prob_acc.tmp = exp(out_swap$logprob)
      
      if (d > ctr_swap$start & ctr_swap$adaptive) {
        rate.tmp = update_rate(iter = (d-ctr_swap$start) / ctr_swap$swap_step, gamma = ctr_swap$rate_gamma, 
                               delta = ctr_swap$rate_delta, kappa = ctr_swap$rate_kappa)
        rho.tmp = update_rho(rho = rho.tmp, alpha = prob_move.tmp, rate = rate.tmp, target = ctr_swap$prob_target)
        temperatures.tmp = update_temp_schedule(rho.tmp, normalize = ctr_swap$normalize, decreasing = TRUE)
        
      }
      
      rho[, d / ctr_swap$swap_step] = rho.tmp
      prob_move[, d / ctr_swap$swap_step] = prob_move.tmp
      # prob_acc[, d / swap_step] = prob_acc.tmp
      
      for (t_idx in seq_along(traces)){
        traces_tmp[[t_idx]]$temperature = temperatures.tmp[t_idx]
      }
    }
    
    tempering_schedule[ , d] = temperatures.tmp
    
    if (ctr_mcmc$verbose > 0) {
      if (d %% ctr_mcmc$print_step == 0){
        now = Sys.time()
        diff_time_start = round(difftime(now, time_start, units = "mins"), 2)
        diff_time_last = round(difftime(now, time_last, units = "mins"), 2)
        clock_time = format(Sys.time(), "%H:%M:%S")
        cat("\n", 
            "iter ", d, " done! \n  ",
            "tempering schedule: ", paste(round(temperatures.tmp, 3), collapse = " "), "\n  ",
            "mean prob acc: ", paste(round(rowMeans(prob_move, na.rm = T), 2), collapse = " "), "\n  ",
            "% succ. swaps: ", paste(round(succesful_swap/attempted_swap, 2), collapse = " "), "\n  ",
            "total exe time: ", diff_time_start, " mins \n  ",
            "last ", ctr_mcmc$print_step, " iter time: ", diff_time_last, " mins \n  ", 
            "print clock time: ", clock_time, "\n", sep = "")
        time_last = now
      }
    }
    
    
    if ( (d %% iter_to_save == 0) & ctr_save$save){
      cat("\n Intermediate save at iter. ", d, sep = " ")
      
      counter_save = counter_save + 1L
      
      exec_time = Sys.time() - time_start
      
      cat("\n")
      cat("ATT:", attempted_swap, "\n")
      cat("SUCC:", succesful_swap, "\n")
      
      traces_to_save = list()
      for (j in seq_along(traces)){
        
        traces_to_save[[j]] = list()
        
        # Partitions
        traces_to_save[[j]]$C = traces[[j]]$C[ , , (d-iter_to_save+1):d]
        
        # alpha
        traces_to_save[[j]]$alpha = traces[[j]]$alpha[ , (d-iter_to_save+1):d]
        
        # gamma
        traces_to_save[[j]]$gamma = traces[[j]]$gamma[ , (d-iter_to_save+1):d]
        
        # loglik
        traces_to_save[[j]]$log_lik = traces[[j]]$log_lik[(d-iter_to_save+1):d]
      }
      
      out <- list("traces" = traces_to_save,
                  "exec_time" = exec_time,
                  par = list(par_likelihood = par_likelihood,
                             par_tRPM = par_tRPM),
                  ctr = list(ctr_mcmc = ctr_mcmc,
                             ctr_save = ctr_save,
                             ctr_swap = ctr_swap),
                  Y = Y, m = m,
                  "diagnostics" = list("attempted_swap" = attempted_swap,
                                       "succesfull_swap" = succesful_swap,
                                       "chains_tracking" = chains_tracking[ , (d-iter_to_save+1):d],
                                       "loglik_tracking" = loglik_tracking[ , (d-iter_to_save+1):d],
                                       "prob_acc" = prob_acc[ , (d-iter_to_save+1):d],
                                       "prob_move" = prob_move[ , (d-iter_to_save+1):d],
                                       "rho" = rho[ , (d-iter_to_save+1):d]))
      
      inter_save_filename = paste0(gsub(".RDS", "", ctr_save$filename), "_part", counter_save, ".RDS")
      
      if (ctr_save$save) {
        if (!dir.exists(ctr_save$filepath)) {
          warning("The provided path_save does not exist. Creating it. \n") 
          dir.create(ctr_save$filepath, recursive = T)}
        cat("Saving out in RDS at path: ", paste0(ctr_save$filepath, inter_save_filename, sep = ""), "\n")
        saveRDS(out, paste(ctr_save$filepath, inter_save_filename, sep = ""))
      }
      
      rm(out); rm(traces_to_save); gc()
    } else if (d == 2 & ctr_save$save){
      exec_time = Sys.time() - time_start
      
      cat("\n")
      cat("ATT:", attempted_swap, "\n")
      cat("SUCC:", succesful_swap, "\n")
      
      out <- list("traces" = traces,
                  "exec_time" = exec_time,
                  par = list(par_likelihood = par_likelihood,
                             par_tRPM = par_tRPM),
                  ctr = list(ctr_mcmc = ctr_mcmc,
                             ctr_save = ctr_save,
                             ctr_swap = ctr_swap),
                  Y = Y, m = m,
                  "diagnostics" = list("attempted_swap" = attempted_swap,
                                       "succesfull_swap" = succesful_swap,
                                       "chains_tracking" = chains_tracking,
                                       "loglik_tracking" = loglik_tracking,
                                       "prob_acc" = prob_acc,
                                       "prob_move" = prob_move,
                                       "rho" = rho))
      
      if (ctr_save$save) {
        if (!dir.exists(ctr_save$filepath)) {
          warning("The provided path_save does not exist. Creating it. \n") 
          dir.create(ctr_save$filepath, recursive = T)}
        cat("Saving out in RDS at path: ", paste0(ctr_save$filepath, ctr_save$filename, sep = ""), "\n")
        saveRDS(out, paste(ctr_save$filepath, ctr_save$filename, sep = ""))
      }
    }
  }
  
  if(ctr_mcmc$parallel){
    parallel::stopCluster(clust)
  }
  
  exec_time = Sys.time() - time_start
  
  cat("\n")
  cat("ATT:", attempted_swap, "\n")
  cat("SUCC:", succesful_swap, "\n")
  
  out = list("traces" = traces_to_save,
             "exec_time" = exec_time,
             par = list(par_likelihood = par_likelihood,
                        par_tRPM = par_tRPM),
             ctr = list(ctr_mcmc = ctr_mcmc,
                        ctr_save = ctr_save,
                        ctr_swap = ctr_swap),
             Y = Y, m = m,
             "diagnostics" = list("attempted_swap" = attempted_swap,
                                  "succesfull_swap" = succesful_swap,
                                  "chains_tracking" = chains_tracking,
                                  "loglik_tracking" = loglik_tracking,
                                  "prob_acc" = prob_acc,
                                  "prob_move" = prob_move,
                                  "rho" = rho))
  
  cat("End Gibbs \n")
  
  if (ctr_save$save) {
    if (!dir.exists(ctr_save$filepath)) {
      warning("The provided path_save does not exist. Creating it. \n") 
      dir.create(ctr_save$filepath, recursive = T)}
    cat("Saving out in RDS at path: ", paste0(ctr_save$filepath, ctr_save$filename, sep = ""), "\n")
    saveRDS(out, paste(ctr_save$filepath, ctr_save$filename, sep = ""))
  }
  
  
  
  return(out)
}


tHMM_local_update = function(
    d,
    Y, m, 
    trace,
    par_likelihood,
    par_tRPM,
    verbose, 
    print_step
) {
  
  # Set the data dimensions
  dims = dim(Y)
  n = dims[1]
  p = dims[2]
  .T = dims[3]
  attrlist = lapply(m, function(x) 1:x)
  
  # Extract hyperparameters
  u = par_likelihood$u
  v = par_likelihood$v
  eta = par_tRPM$eta
  a_alpha = par_tRPM$a_alpha
  b_alpha = par_tRPM$b_alpha
  
  # Set the urn
  urn_type = par_tRPM$urn_type
  urn = switch(urn_type,
               "Gnedin" = urn_GN,
               "DP" = urn_DP,
               stop("Invalid urn_type"))
  
  
  # Current values
  C.local = trace$C.tmp
  mu.local = trace$mu.tmp
  alpha.local = trace$alpha.tmp
  gamma.local = trace$gamma.tmp
  clS.local = trace$clS.tmp
  matches.local = trace$matches.tmp
  log_lik.local = 0
  
  # Pre-compute all possible values of normalization constant
  cached_norm_const = trace$cached_norm_const
  
  # Pre-compute the denominator for the marginal loglik
  denomin_log_lik = trace$denomin_log_lik
  
  # Start loop in t ----
  for (t in 1:.T) {
    
    Y_t = Y[ , , t]
    
    if (verbose > 0) {
      if (d %% print_step == 0) {
        cat("\n")
        cat("iter: ", d, "  ", sep = "")
        cat("time: ", t, "  ", sep = "")
      }
    }
    
    
    # UPDATE GAMMA
    if (verbose == 2) {
      if (d %% print_step == 0) {
        cat("\n", "  update: gammas", sep = "")
      }
    }
    
    # Shuffle order update gammas
    tovisit.gamma = sample(1:n, n, replace = FALSE)
    
    # Start loop in j.gamma ----
    for (j.gamma in 1:n) {
      # Select the next node to visit
      i = tovisit.gamma[1]
      tovisit.gamma = tovisit.gamma[-1]
      
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\r",
              "  update: ", "gammas",
              gettextf("  iter: %3d", j.gamma),
              gettextf("  node: %3d", i),
              sep = "")
        }
      }
      
      if (t == 1) {
        gamma.local[i, t] = 0
      } else {
        gamma.local[i, t] = update_gamma(i = i,
                                         gamma_t = gamma.local[, t],
                                         alpha_t = alpha.local[t],
                                         C_t = C.local[[t]],
                                         C_tm1 = C.local[[t - 1]],
                                         eta = eta)
      }
    } 
    # End loop in j.gamma ----
    
    if (verbose > 0) {
      if (d %% print_step == 0) {
        cat(ifelse(verbose == 2, "\r", "\n"),
            sprintf("  update: gammas  done! (%3d gamma=1)%-20s",
                    sum(gamma.local[, t]),
                    ""),
            sep = "")
      }
    }
    
    
    # UPDATE PARTITIONS
    if (verbose == 2) {
      if (d %% print_step == 0) {
        cat("\n", "  update: labels", sep = "")
      }
    }
    
    # Shuffle order update labels
    # visited.labels = c()
    tovisit.labels = sample(1:n, n, replace = FALSE)
    
    # Start loop in j.lab ----
    for (j.lab in 1:n) {
      # Select the next node to visit
      i = tovisit.labels[1]
      tovisit.labels = tovisit.labels[-1]
      # visited.labels = c(visited.labels, i)
      
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\r",
              gettextf("  update: %s", "labels"),
              gettextf("  iterator: %3d", j.lab),
              gettextf("  node: %3d", i),
              gettextf("  ncl: %3d", ncol(C.local[[t]])),
              sep = "")
        }
      }
      
      
      out_uplab = update_label(i = i,
                               Y_it = Y_t[i, ],
                               C_t = C.local[[t]],
                               C_tp1 = C.local[[t + 1]],
                               gamma_t = gamma.local[, t],
                               gamma_tp1 = gamma.local[, t + 1],
                               eta = eta,
                               mu = mu.local[[t]],
                               clS_t = clS.local[[t]],
                               matches_t = matches.local[[t]],
                               cached_I = cached_norm_const,
                               m = m,
                               u = u,
                               v = v,
                               urn = urn)
      
      C.local[[t]] = out_uplab$lab
      mu.local[[t]] = out_uplab$center
      clS.local[[t]] = out_uplab$clS
      matches.local[[t]] = out_uplab$matches
    }
    # End loop in j.lab ----
    
    if (verbose > 0) {
      if (d %% print_step == 0) {
        cat(ifelse(verbose == 2, "\r", "\n"), 
            sprintf("  update: labels  done! (ncl: %3d, sing: %3d)%-20s",
                    ncol(C.local[[t]]),
                    sum(colSums(C.local[[t]]) == 1),
                    ""),
            sep = "")
      }
    }
    
    # UPDATE CENTERS (MU)
    if (verbose == 2) {
      if (d %% print_step == 0) {
        cat("\n  update: cluster parameters")
      }
    }
    ncluster = ncol(C.local[[t]])
    
    # Start loop in h ----
    for (h in 1:ncluster){
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\r",
              gettextf("  update: %s", "centers"),
              gettextf("  iter: %3d", d),
              gettextf("  cluster: %3d", h),
              sep = "")
        }
      }
      
      idx_h = (C.local[[t]][ , h] == 1)
      Y_th = Y_t[idx_h, , drop = FALSE]
      matches_h = count_matches_update_centers(Y_th, m)
      clS_h = sum(idx_h)
      
      logprob.local = eval_logprob_update_centers(M = Y_th, m = m,
                                                  clS = clS_h, 
                                                  tabs = cached_norm_const)
      
      draw.local = sapply(logprob.local, sample_gumbel_max)
      mu.local[[t]][h,] = unname(draw.local)
    }
    # End loop in h ----  
    
    # After the sampling of the new centers I have to update the matches counts
    tmp = matrix(NA, ncluster, p)
    for (k in 1:ncluster){
      idx_t = (C.local[[t]][, k] == 1)
      tmp[k, ] = colSums(Y[idx_t, , t] == matrix(mu.local[[t]][k, ], 
                                                 nrow = sum(idx_t), 
                                                 ncol = p, byrow = TRUE))
    }
    matches.local[[t]] = tmp
    
    
    if (verbose > 0) {
      if (d %% print_step == 0) {
        cat(ifelse(verbose == 2, "\r", "\n"),
            sprintf(" update: parameters  done!%-20s",
                    "")
        )
      }
    }
    
    # UPDATE ALPHA
    if (verbose == 2) {
      if (d %% print_step == 0) {
        cat("\n  update: alpha")
      }
    }
    if (t > 1) {
      alpha.local[t] = update_alpha(gamma_t = gamma.local[, t],
                                    pr_shape1 = a_alpha,
                                    pr_shape2 = b_alpha)
    }
    if (verbose > 0) {
      if (d %% print_step == 0) {
        cat(ifelse(verbose == 2, "\r", "\n"),
            "  update: alpha  done!")
      }
    }
    
    
    log_lik_t.local = NA
    numer_log_lik = eval_norm_const(clS = clS.local[[t]], matches = matches.local[[t]],
                                    tabs = cached_norm_const)
    log_lik_t.local = sum(colSums(numer_log_lik) - denomin_log_lik)
    
    log_lik.local = log_lik.local + log_lik_t.local
    
  }
  # End loop in t ----
  
  trace$C.tmp = C.local
  trace$mu.tmp = mu.local
  trace$alpha.tmp = alpha.local
  trace$gamma.tmp = gamma.local
  trace$clS.tmp = clS.local
  trace$matches.tmp = matches.local
  trace$log_lik.tmp = log_lik.local
  
  return(trace)
}

tHMM_chain_swap = function(traces, d, temperatures, deterministic = FALSE, n_try, family, prior){
  
  n_tr = length(traces)
  n_try = ifelse(n_try <= floor((n_tr-1)/2), n_try, floor((n_tr-1)/2))
  n_try = max(1, n_try)
  
  loglik = sapply(traces, function(tr) tr$log_lik.tmp)
  
  if (deterministic){
    if (n_tr == 2){
      forward = 1
      backward = 2
    } else {
      if (d %% 2 == 0) {
        forward = seq(2, n_tr - 1, by = 2)
      } else {
        forward = seq(1, n_tr - 1, by = 2)
      }
      backward = forward + 1
    }
  } else {
    even = seq(2, n_tr - 1, by = 2)
    odd = seq(1, n_tr - 1, by = 2)
    
    if (runif(1) < 0.5) {
      forward = sample(x = even, size = n_try)
    } else {
      forward = sample(x = odd, size = n_try)
    }
    backward = forward + 1
  }
  
  diff_loglik = loglik[backward] - loglik[forward]
  diff_temp = temperatures[backward] - temperatures[forward]
  
  logprob = pmin(0, - diff_temp * diff_loglik)
  U = runif(length(forward), 0, 1)
  swap = forward[log(U) < logprob]
  
  for (tr in forward){
    if (tr %in% swap){
      # Temporarily save tr objects
      C.tmp = traces[[tr]]$C.tmp
      mu.tmp = traces[[tr]]$mu.tmp
      alpha.tmp = traces[[tr]]$alpha.tmp
      gamma.tmp = traces[[tr]]$gamma.tmp
      Cls.tmp = traces[[tr]]$Cls.tmp
      matches.tmp = traces[[tr]]$matches.tmp
      
      # Put objects of tr+1 in tr
      traces[[tr]]$C.tmp = traces[[tr+1]]$C.tmp
      traces[[tr]]$mu.tmp = traces[[tr+1]]$mu.tmp
      traces[[tr]]$alpha.tmp = traces[[tr+1]]$alpha.tmp
      traces[[tr]]$gamma.tmp = traces[[tr+1]]$gamma.tmp
      traces[[tr]]$Cls.tmp = traces[[tr+1]]$Cls.tmp
      traces[[tr]]$matches.tmp = traces[[tr+1]]$matches.tmp
      
      # Put objects of tr in tr+1
      traces[[tr+1]]$C.tmp = C.tmp
      traces[[tr+1]]$mu.tmp = mu.tmp
      traces[[tr+1]]$alpha.tmp = alpha.tmp
      traces[[tr+1]]$gamma.tmp = gamma.tmp
      traces[[tr+1]]$Cls.tmp = Cls.tmp
      traces[[tr+1]]$matches.tmp = matches.tmp
    }
  }
  
  attempted_swap = rep(0, n_tr)
  attempted_swap[forward] = 1
  succesful_swap = rep(0, n_tr)
  succesful_swap[swap] = 1
  
  return(list("traces" = traces, 
              "diagnostics" = list("att" = attempted_swap, 
                                   "succ" = succesful_swap),
              "diff_loglik" = diff_loglik,
              "diff_temp" = diff_temp,
              "logprob" = logprob,
              "loglik" = loglik))
}


tHMM_gibbs_parameters = function(
    PPE, Y, m,
    par_likelihood = list(u = NULL, # First hyperparam HIG
                          v = NULL # Second hyperparam HIG
    ),
    ctr_mcmc = list(seed = 1234, nburnin = 1000, nchain = 5000, print_step = 100,
                    verbose = "0"),
    ctr_save = list(save = FALSE, filepath = "",
                    filename = paste("par_", Sys.Date(), ".RDS", sep = ""))
){
  # Set the Y dimensions
  dims = dim(Y)
  n = dims[1]
  p = dims[2]
  .T = dims[3]
  ncl = apply(PPE, 2, function(x) length(unique(c(x)))) %>% unname()
  m.inner = apply(Y, 2, function(x) length(unique(c(x)))) %>% unname()
  if (any(m != m.inner)){stop("Provided m is not consistent with Y!")}
  attrlist = lapply(m, function(x) 1:x)
  
  # Set and check the control parameters
  ctr_mcmc = do.call("set_ctr_mcmc_par", ctr_mcmc)
  ctr_save = do.call("set_ctr_save_par", ctr_save)
  niter = ctr_mcmc$nchain + ctr_mcmc$nburnin
  
  # Extract hyperparameters
  u = par_likelihood$u
  v = par_likelihood$v
  
  # Initialize data storage
  mu = list()
  sigma = list()
  mu.tmp = list()
  sigma.tmp = list()
  
  for (t in 1:.T) {
    mu[[t]] = sigma[[t]] = array(NA, dim = c(ncl[t], p, niter))
    mu.tmp[[t]] = sigma.tmp[[t]] = array(NA, dim = c(ncl[t], p))
    
    for (h in 1:ncl[t]){
      for (j in 1:p){
        mu.tmp[[t]][h, j] = sample(1:m[j], 1)
        sigma.tmp[[t]][h, j] = 1
      }
    }
  }
  
  time_start = time_last = Sys.time()
  
  # Gibbs sampling loop
  for (iter in 1:niter) {
    for (t in 1:.T){
      
      Y_t = Y[ , , t]
      ncluster = ncl[t]
      
      for (h in 1:ncluster){
        
        idx_h = PPE[ , t] == h
        Y_th = Y_t[idx_h, , drop = FALSE]
        
        prob.tmp = Center_prob(data = Y_th,
                               sigma = sigma.tmp[[t]][h,],
                               attrisize = m)
        
        mu.tmp[[t]][h,] = Samp_Center(center_prob = prob.tmp,
                                      attriList = attrlist, p = p)
        
        size_h = sum(idx_h)
        
        for (j in 1:p) {
          
          dd = sum(Y_t[idx_h, j] != mu.tmp[[t]][h, j])
          cc = size_h - dd
          
          sigma.tmp[[t]][h, j] = rhyper_sig2(n = 1,
                                             d = v[j]+dd,
                                             c = u[j]+cc,
                                             m = m[j])
        }
      }
      
      mu[[t]][ , , iter] = mu.tmp[[t]]
      sigma[[t]][ , , iter] = sigma.tmp[[t]]
    }
    
    if (ctr_mcmc$verbose > 0 && iter %% ctr_mcmc$print_step == 0) {
      now = Sys.time()
      diff_time_start = round(difftime(now, time_start, units = "mins"), 2)
      diff_time_last = round(difftime(now, time_last, units = "mins"), 2)
      clock_time = format(Sys.time(), "%H:%M:%S")
      cat("\n\n",
          "##############################\n",
          "Iter ", iter, " done!\n",
          "total exe time: ", diff_time_start, " mins \n  ",
          "last ", ctr_mcmc$print_step, " iter time: ", diff_time_last, " mins \n  ", 
          "print clock time: ", clock_time, "\n", 
          "##############################\n",
          sep = "")
      time_last = now
    }
  }
  
  if (ctr_save$save){
    saveRDS(list(mu = mu, sigma = sigma), 
            paste(ctr_save$filepath, ctr_save$filename, sep = ""))
  }
  
  return(list(mu = mu, sigma = sigma))
}

tHMM_gibbs_initialize = function(Y, n, .T, niter,
                                 ncl.init, p){
  # Initialize data storage
  C = array(NA_integer_, dim = c(n, .T, niter))
  alpha = array(NA_real_, dim = c(.T, niter))
  gamma = array(NA_real_, dim = c(.T, niter))
  log_lik = rep(NA_real_, niter)
  
  # Initial values
  C.initmp = sample(1:ncl.init, n, TRUE)
  C.initmp = match(C.initmp, sort(unique(C.initmp))) 
  C.init = matrix(C.initmp, nrow = n, ncol = .T)
  mu.init = replicate(.T, matrix(NA, nrow = ncl.init, ncol = p), simplify = F)
  for (t in 1:.T){
    for (h in 1:ncl.init){
      for (j in 1:p)
        mu.init[[t]][h, j] = sample(1:m[j], 1)
    }
  }
  alpha.init = rep(0.01, .T)
  alpha.init[1] = NA
  gamma.init = matrix(0, nrow = n, ncol = .T+1)
  
  # Current values
  C.tmp = apply(C.init, 2, function(x) vec2mat(x), simplify = F)
  mu.tmp = mu.init
  alpha.tmp = alpha.init
  gamma.tmp = gamma.init
  clS.tmp = lapply(C.tmp, function(c) colSums(c))
  matches.tmp = list()
  for (t in 1:.T){
    ncl_t = ncol(C.tmp[[t]])
    tmp = matrix(NA, ncl_t, p)
    for (k in 1:ncl_t){
      idx_t = (C.tmp[[t]][, k] == 1)
      tmp[k, ] = colSums(Y[idx_t, , t] == matrix(mu.tmp[[t]][k, ], 
                                                 nrow = sum(idx_t), 
                                                 ncol = p, byrow = TRUE))
    }
    matches.tmp[[t]] = tmp
  }
  log_lik.tmp = NA
  
  # Pre-compute all possible values of normalization constant
  cached_norm_const = build_cached_norm_const(Smax = n, u = u, v = v, m = m)
  
  # Pre-compute the denominator for the marginal loglik
  denomin_log_lik = sapply(cached_norm_const, function(x) x[[1]])
  
  
  out = list(C = C, alpha = alpha, gamma = gamma, log_lik = log_lik,
             C.tmp = C.tmp, mu.tmp = mu.tmp, alpha.tmp = alpha.tmp, gamma.tmp = gamma.tmp, 
             clS.tmp = clS.tmp, matches.tmp = matches.tmp, log_lik.tmp = log_lik.tmp,
             cached_norm_const = cached_norm_const, denomin_log_lik = denomin_log_lik)
}