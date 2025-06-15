tHMM_gibbs = function(
    Y,
    m, 
    # Parameters for the likelihood-specific local parameters
    par_likelihood = list(u = NULL, # First hyperparam HIG
                          v = NULL # Second hyperparam HIG
    ),
    # Parameters for the tRPM process prior
    par_tRPM = list(a_alpha = 1, b_alpha = 1, # Beta prior on alpha
                    urn_type = c("Gnedin", "DP"), 
                    eta = 1 # Gnedin/DP parameter
    ),
    # Control parameters for MCMC settings
    ctr_mcmc = list(seed = 1234, nburnin = 1000, nchain = 5000, print_step = 100, 
                    ncl.init = NULL, verbose = "0"),
    # Control parameters for result storing
    ctr_save = list(save = FALSE, filepath = "",
                    filename = paste("res_", Sys.Date(), ".RDS", sep = "")),
    # Control parameters for managing of alpha parameters
    ctr_alpha = list(fix_alpha.flag = FALSE, fix_alpha.value = NULL)
) {
  
  if (ctr_save$save) {
    if (!dir.exists(ctr_save$filepath)) {
      warning("The provided path_save does not exist. Creating it. \n")
      dir.create(ctr_save$filepath, recursive = T)
    }
    cat("##############################\n", "Saving out in RDS at path: ", paste0(ctr_save$filepath, ctr_save$filename, sep = ""), "\n", sep = "")
  }
  
  # Set the data dimensions
  dims = dim(Y)
  n = dims[1]
  p = dims[2]
  .T = dims[3]
  m.inner = apply(Y, 2, function(x) length(unique(c(x)))) %>% unname()
  if (any(m != m.inner)){stop("Provided m is not consistent with Y!")}
  attrlist = lapply(m, function(x) 1:x)
  
  # Set and check the control parameters
  ctr_mcmc = do.call("set_ctr_mcmc", ctr_mcmc)
  ctr_save = do.call("set_ctr_save", ctr_save)
  ctr_alpha = do.call("set_ctr_alpha", ctr_alpha)
  
  # Extract ctr_alpha
  fix_alpha.flag = ctr_alpha$fix_alpha.flag
  fix_alpha.value = ifelse(ctr_alpha$fix_alpha.value, 
                           ctr_alpha$fix_alpha.value, 
                           NA)
  
  verbose = ctr_mcmc$verbose
  print_step = ctr_mcmc$print_step
  
  # Set the simulation seed
  set.seed(ctr_mcmc$seed)
  
  # Set the number of iterations and swaps
  niter = ctr_mcmc$nburnin + ctr_mcmc$nchain
  
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
  
  # Initialize data storage
  C = array(NA_integer_, dim = c(n, .T, niter))
  alpha = array(NA_real_, dim = c(.T, niter))
  gamma = array(NA_real_, dim = c(.T, niter))
  
  
  # Initial values
  ncl.init = ctr_mcmc$ncl.init
  C.initmp = sample(1:ncl.init, n, TRUE)
  C.initmp = match(C.initmp, sort(unique(C.initmp))) 
  C.init = matrix(C.initmp, nrow = n, ncol = .T)
  # mu.init = apply(Y, 3, function(x) x, simplify = F)
  mu.init = replicate(.T, matrix(NA, nrow = ncl.init, ncol = p), simplify = F)
  for (t in 1:.T){
    for (h in 1:ncl.init){
      for (j in 1:p)
        mu.init[[t]][h, p] = sample(1:m[j], 1)
    }
  }
  sigma.init = replicate(.T, matrix(1, nrow = ncl.init, ncol = p), simplify = FALSE)
  alpha.init = rep(0.01, .T)
  alpha.init[1] = NA
  gamma.init = matrix(0, nrow = n, ncol = .T+1)
  
  # Current values
  C.tmp = apply(C.init, 2, function(x) vec2mat(x), simplify = F)
  mu.tmp = mu.init
  sigma.tmp = sigma.init
  alpha.tmp = alpha.init
  gamma.tmp = gamma.init
  
  # Counter for intermediate saves
  counter_save = 0L
  # How many iter between intermediate saves
  iter_to_save = 2500
  
  
  # Elapsed time initialization
  time_start = time_last = Sys.time()
  cat("Start Gibbs of ", niter, " iterations\n", 
      "##############################\n",
      sep = "")
  
  # Gibbs sampling loop
  
  # Start loop in d ----
  for (d in 1:niter) {
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
      visited.gamma = c()
      tovisit.gamma = sample(1:n, n, replace = FALSE)
      
      # Start loop in j.gamma ----
      for (j.gamma in 1:n) {
        # Select the next node to visit
        i = tovisit.gamma[1]
        tovisit.gamma = tovisit.gamma[-1]
        visited.gamma = c(visited.gamma, i)
        
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
          gamma.tmp[i, t] = 0
        } else {
          gamma.tmp[i, t] = update_gamma(i = i,
                                         gamma_t = gamma.tmp[, t],
                                         alpha_t = alpha.tmp[t],
                                         C_t = C.tmp[[t]],
                                         C_tm1 = C.tmp[[t - 1]],
                                         eta = eta)
        }
      } 
      # End loop in j.gamma ----
      
      if (verbose > 0) {
        if (d %% print_step == 0) {
          cat(ifelse(verbose == 2, "\r", "\n"),
              sprintf("  update: gammas  done! (%3d gamma=1)%-20s",
                      sum(gamma.tmp[, t]),
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
      visited.labels = c()
      tovisit.labels = sample(1:n, n, replace = FALSE)
      
      # Start loop in j.lab ----
      for (j.lab in 1:n) {
        # Select the next node to visit
        i = tovisit.labels[1]
        tovisit.labels = tovisit.labels[-1]
        visited.labels = c(visited.labels, i)
        
        if (verbose == 2) {
          if (d %% print_step == 0) {
            cat("\r",
                gettextf("  update: %s", "labels"),
                gettextf("  iterator: %3d", j.lab),
                gettextf("  node: %3d", i),
                gettextf("  ncl: %3d", ncol(C.tmp[[t]])),
                sep = "")
          }
        }
        
        out_uplab = update_label(i = i,
                                 Y_it = Y_t[i, ],
                                 C_t = C.tmp[[t]],
                                 C_tp1 = C.tmp[[t + 1]],
                                 gamma_t = gamma.tmp[, t],
                                 gamma_tp1 = gamma.tmp[, t + 1],
                                 eta = eta,
                                 mu = mu.tmp[[t]],
                                 sigma = sigma.tmp[[t]], 
                                 m = m,
                                 u = u,
                                 v = v,
                                 urn = urn)
        
        C.tmp[[t]] = out_uplab$lab
        mu.tmp[[t]] = out_uplab$center
        sigma.tmp[[t]] = out_uplab$scale
      }
      # End loop in j.lab ----
      
      if (verbose > 0) {
        if (d %% print_step == 0) {
          cat(ifelse(verbose == 2, "\r", "\n"), 
              sprintf("  update: labels  done! (ncl: %3d, sing: %3d)%-20s",
                      ncol(C.tmp[[t]]),
                      sum(colSums(C.tmp[[t]]) == 1),
                      ""),
              sep = "")
        }
      }
      
      # UPDATE CENTERS (MU) AND SCALE (SIGMA)
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\n  update: cluster parameters")
        }
      }
      ncluster = ncol(C.tmp[[t]])
      
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
        
        idx_h = C.tmp[[t]][ , h] == 1
        Y_th = Y_t[idx_h, , drop = FALSE]
        
        prob.tmp = Center_prob(data = Y_th,
                               sigma = sigma.tmp[[t]][h,],
                               attrisize = m)
        
        mu.tmp[[t]][h,] = Samp_Center(center_prob = prob.tmp,
                                      attriList = attrlist, p=p)
        
        size_h = sum(idx_h)
        
        for (j in 1:p) {
          
          dd = sum(Y_t[idx_h , j] != mu.tmp[[t]][h, j])
          cc = size_h - dd
          
          sigma.tmp[[t]][h, j] = rhyper_sig2(n = 1,
                                             d = v[j]+dd,
                                             c = u[j]+cc,
                                             m = m[j])
        }
      }
      # End loop in h ----  
      
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
        if (!fix_alpha.flag){
          alpha.tmp[t] = update_alpha(gamma_t = gamma.tmp[, t],
                                      pr_shape1 = a_alpha,
                                      pr_shape2 = b_alpha)
        } else {
          alpha.tmp[t] = fix_alpha.value
        }
      }
      if (verbose > 0) {
        if (d %% print_step == 0) {
          cat(ifelse(verbose == 2, "\r", "\n"),
              "  update: alpha  done!")
        }
      }
      
      # Store the results
      C[ , t, d] = mat2vec(C.tmp[[t]])
      alpha[t, d] = alpha.tmp[t]
      gamma[t, d] = mean(gamma.tmp[ , t])
    }
    # End loop in t ----
    
    if (verbose > 0) {
      if (d %% print_step == 0) {
        now = Sys.time()
        diff_time_start = round(difftime(now, time_start, units = "mins"), 2)
        diff_time_last = round(difftime(now, time_last, units = "mins"), 2)
        clock_time = format(Sys.time(), "%H:%M:%S")
        cat("\n\n",
            "##############################\n",
            "total exe time: ", diff_time_start, " mins \n  ",
            "last ", ctr_mcmc$print_step, " iter time: ", diff_time_last, " mins \n  ", 
            "print clock time: ", clock_time, "\n", 
            "##############################\n",
            sep = "")
        time_last = now
      }
    }
    
    if ( (d %% iter_to_save == 0) & ctr_save$save){
      cat("\n Intermediate save at iter. ", d, sep = " ")
      counter_save = counter_save + 1L
      inter_save_filename = paste0(gsub(".RDS", "", ctr_save$filename), "_part", counter_save, ".RDS")
      now = Sys.time()
      diff_time_start = round(difftime(now, time_start, units = "mins"), 2)
      
      out = list(input = list(Y = Y, m = m,
                              par = list(par_likelihood = par_likelihood,
                                         par_tRPM = par_tRPM),
                              ctr = list(ctr_mcmc = ctr_mcmc,
                                         ctr_save = ctr_save,
                                         ctr_alpha = ctr_alpha)),
                 output = list(C = C[ , , 1:d], alpha = alpha[ , 1:d]),
                 execution_time = diff_time_start)
      
      saveRDS(out, paste(ctr_save$filepath, inter_save_filename, sep = ""))
      
      rm(out); gc()
    } else if (d == 2 & ctr_save$save){
      now = Sys.time()
      diff_time_start = round(difftime(now, time_start, units = "mins"), 2)
      out = list(input = list(Y = Y, m = m,
                              par = list(par_likelihood = par_likelihood,
                                         par_tRPM = par_tRPM),
                              ctr = list(ctr_mcmc = ctr_mcmc,
                                         ctr_save = ctr_save,
                                         ctr_alpha = ctr_alpha)),
                 output = list(C = C, alpha = alpha),
                 execution_time = diff_time_start)
      saveRDS(out, paste(ctr_save$filepath, ctr_save$filename, sep = ""))
      rm(out); gc()
    }
  } # End loop in d ----
  
  now = Sys.time()
  diff_time_start = round(difftime(now, time_start, units = "mins"), 2)
  diff_time_last = round(difftime(now, time_last, units = "mins"), 2)
  clock_time = format(Sys.time(), "%H:%M:%S")
  cat("\n\n",
      "##############################\n",
      "END OF GIBBS SAMPLER\n",
      "Total exe time: ", diff_time_start, " mins \n  ",
      "##############################\n",
      sep = "")
  
  out = list(input = list(Y = Y, m = m,
                          par = list(par_likelihood = par_likelihood,
                                     par_tRPM = par_tRPM),
                          ctr = list(ctr_mcmc = ctr_mcmc,
                                     ctr_save = ctr_save,
                                     ctr_alpha = ctr_alpha)),
             output = list(C = C, alpha = alpha),
             execution_time = diff_time_start)
  
  if (ctr_save$save) {
    if (!dir.exists(ctr_save$filepath)) {
      warning("The provided path_save does not exist. Creating it. \n") 
      dir.create(ctr_save$filepath, recursive = T)}
    cat("Saving out in RDS at path: ", paste0(ctr_save$filepath, ctr_save$filename, sep = ""), "\n")
    saveRDS(out, paste(ctr_save$filepath, ctr_save$filename, sep = ""))
  }
  
  return(out)
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
            file.path(ctr_save$filepath, paste0(ctr_save$filename, ".RDS")))
  }
  
  return(list(mu = mu, sigma = sigma))
}