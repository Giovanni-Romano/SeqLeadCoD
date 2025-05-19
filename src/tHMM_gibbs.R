tHMM_gibbs = function(
    Y,
    m, 
    # Parameters for the likelihood-specific local parameters
    par_likelihood = list(u = NULL, # First hyperparam HIG
                          v = NULL # Second hyperparam HIG
    ),
    # Parameters for the tRPM process prior
    par_tRPM = list(a_alpha = 1, b_alpha = 1, # Beta prior on alpha
                    eta = 1 # Gnedin parameter
    ),
    # Control parameters for MCMC settings
    ctr_mcmc = list(seed = 1234, nburnin = 1000, nchain = 5000, print_step = 100, verbose = "0"),
    # Control parameters for result storing
    ctr_save = list(save = FALSE, filepath = "",
                    filename = paste("res_", Sys.Date(), ".RDS", sep = ""))
) {
  
  if (ctr_save$save) {
    if (!dir.exists(ctr_save$filepath)) {
      warning("The provided path_save does not exist. Creating it. \n")
      dir.create(ctr_save$filepath, recursive = T)
    }
    cat("Saving out in RDS at path: ", paste0(ctr_save$filepath, ctr_save$filename, sep = ""), "\n")
  }
  
  
  # Set and check the control parameters
  ctr_mcmc = do.call("set_ctr_mcmc", ctr_mcmc)
  ctr_save = do.call("set_ctr_save", ctr_save)
  
  verbose = ctr_mcmc$verbose
  print_step = ctr_mcmc$print_step
  
  # Set the simulation seed
  set.seed(ctr_mcmc$seed)
  
  # Set the number of iterations and swaps
  niter = ctr_mcmc$nburnin + ctr_mcmc$nchain
  
  # Set the data dimensions
  dims = dim(Y)
  n = dims[1]
  p = dims[2]
  .T = dims[3]
  m.inner = apply(Y, 2, function(x) length(unique(c(x)))) %>% unname()
  if (any(m != m.inner)){stop("Provided m is not consistent with Y!")}
  attrlist = lapply(m, function(x) 1:x)
  
  # Extract hyperparameters
  u = par_likelihood$u
  v = par_likelihood$v
  eta = par_tRPM$eta
  a_alpha = par_tRPM$a_alpha
  b_alpha = par_tRPM$b_alpha
  
  # Initialize data storage
  C = array(NA_integer_, dim = c(n, .T, niter))
  alpha = array(NA_real_, dim = c(.T, niter))
  
  
  # Initial values
  ncl.init = n
  C.init = matrix(1:ncl.init, nrow = ncl.init, ncol = .T)
  # mu.init = replicate(.T, sapply(m, function(attrsize){
  #   sample(1:attrsize, ncl.init, replace = TRUE)
  # }), simplify = FALSE)
  mu.init = apply(Y, 3, function(x) x, simplify = F)
  sigma.init = replicate(.T, matrix(1, nrow = ncl.init, ncol = p), simplify = FALSE)
  alpha.init = rep(0.5, .T)
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
  
  # Gibbs sampling loop
  for (d in 1:niter) {
    for (t in 1:.T) { # Loop in t ----
      
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
      
      # Shuffle order update gammass
      visited.gamma = c()
      tovisit.gamma = sample(1:n, n, replace = FALSE)
      
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
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\r", 
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
                                 v = v)
        
        C.tmp[[t]] = out_uplab$lab
        mu.tmp[[t]] = out_uplab$center
        sigma.tmp[[t]] = out_uplab$scale
      }
      
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\r", 
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
          
          dd = sum(Y_t[h , j] != mu.tmp[[t]][h, j])
          cc = size_h - dd
          
          sigma.tmp[[t]][h, j] = rhyper_sig2(n = 1,
                                             d = v[j]+dd,
                                             c = u[j]+cc,
                                             m = m[j])
        }
      }
      
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\r", 
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
        alpha.tmp[t] = update_alpha(gamma_t = gamma.tmp[, t],
                                    pr_shape1 = a_alpha,
                                    pr_shape2 = b_alpha)
      }
      if (verbose == 2) {
        if (d %% print_step == 0) {
          cat("\r", "  update: alpha  done!")
        }
      }
      # Store the results
      C[ , t, d] = mat2vec(C.tmp[[.T]])
      alpha[t, d] = alpha.tmp[t]
    }
    
    
    
    if ( (d %% iter_to_save == 0) & ctr_save$save){
      cat("\n Intermediate save at iter. ", d, sep = " ")
      
      counter_save = counter_save + 1L
      exec_time = Sys.time() - time_start
      inter_save_filename = paste0(gsub(".RDS", "", ctr_save$filename), "_part", counter_save, ".RDS")
      
      out = list(input = list(Y = Y, m = m,
                              par = list(par_likelihood = par_likelihood,
                                         par_tRPM = par_tRPM),
                              ctr = list(ctr_mcmc = ctr_mcmc,
                                         ctr_save = ctr_save)),
                 output = list(C = C[ , , 1:d], alpha = alpha[ , 1:d]))
      
      saveRDS(out, paste(ctr_save$filepath, inter_save_filename, sep = ""))
      
      rm(out); rm(traces_to_save); gc()
    } else if (d == 2 & ctr_save$save){
      exec_time = Sys.time() - time_start
      out = list(input = list(Y = Y, m = m,
                              par = list(par_likelihood = par_likelihood,
                                         par_tRPM = par_tRPM),
                              ctr = list(ctr_mcmc = ctr_mcmc,
                                         ctr_save = ctr_save)),
                 output = list(C = C, alpha = alpha))
      saveRDS(out, paste(ctr_save$filepath, ctr_save$filename, sep = ""))
    }
  }
  
  out = list(input = list(Y = Y, m = m,
                          par = list(par_likelihood = par_likelihood,
                                     par_tRPM = par_tRPM),
                          ctr = list(ctr_mcmc = ctr_mcmc,
                                     ctr_save = ctr_save)),
             output = list(C = C, alpha = alpha))
  
  if (ctr_save$save) {
    if (!dir.exists(ctr_save$filepath)) {
      warning("The provided path_save does not exist. Creating it. \n") 
      dir.create(ctr_save$filepath, recursive = T)}
    cat("Saving out in RDS at path: ", paste0(ctr_save$filepath, ctr_save$filename, sep = ""), "\n")
    saveRDS(out, paste(ctr_save$filepath, ctr_save$filename, sep = ""))
  }
}