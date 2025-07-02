# Use Gumbel-max trick to sample from logprob vector ----
sample_gumbel_max = function(logp){
  gumbel = -log(-log(runif(length(logp))))
  return(which.max(logp + gumbel))
}

# Compatibility check ----
compatibility_check = function(i, H, k, gamma_tp1, c_t, c_tp1, safe = FALSE){
  
  # Compatibility checks
  if (safe) {
    # SAFER VERSION OF COMPATIBILITY CHECK
    if (all(gamma_tp1 == 0)){
      # We enter here for sure when t = .T because I set \gamma_{i, .T} = 0 for 
      # every i. Nonetheless, it's still correct for every t, in fact if all 
      # observations are free to move, then it's impossible to have compatibility issues.
      comp_checks.safe = rep(1, H + 1)
    } else {
      R_tp1.safe = which(gamma_tp1 == 1)
      adj_tp1.safe = tcrossprod(c_tp1[R_tp1.safe, , drop = FALSE], 
                                c_tp1[R_tp1.safe, , drop = FALSE])
      
      comp_checks.safe = c()
      for (h in 1:(H + 1)){
        c_t.mod = c_t
        c_t.mod[i, k] = 0
        if (h == (H + 1)){
          c_t.mod = cbind(c_t.mod, 0)
        }
        c_t.mod[i, h] = 1
        adj_t.safe = tcrossprod(c_t.mod[R_tp1.safe, , drop = FALSE], 
                                c_t.mod[R_tp1.safe, , drop = FALSE])
        
        comp_checks.safe[h] = all(adj_t.safe == adj_tp1.safe)
      }
    }
    comp_checks = comp_checks.safe
  } else {
    # FASTER VERSION OF COMPATIBILITY CHECK
    # If the i-th observation is free to move from time t to time t+1, then 
    # compatibility I[\rho_t^{h R_{t+1}} = \rho_{t+1}^{h R_{t+1}}] is guaranteed, 
    # because we are moving only obs. i and whatever c_{i, t} we sample in this 
    # update, then obs. i can move again to c_{i, t+1}.
    # If instead i-th observation is constrained to stay in the same cluster from
    # time t to time t+1 (i.e. gamma_tp1[i] == 1), then we actually need to check 
    # the compatibility. To do so, I'd compare the co-clustering matrices at time
    # t and t+1, relative to the sub-partition of the constrained obs. (i.e. 
    # which(gamma_tp1 == 1)). If they are equal, then the proposal is valid.
    # Actually, I think that it's necessary to check the i-th row/column (co-clust 
    # matrix is symm.) of this matrices because the other rows cannot change just
    # updating c_{i, t}.
    if (gamma_tp1[i] == 0){
      comp_checks = rep(1, H + 1)
    } else {
      R_tp1 = which(gamma_tp1 == 1)
      R_tp1mi = R_tp1[R_tp1 != i]
      
      # If obs. i is the only one to be constrained, then the compatibility is 
      # guaranteed for every possible allocations.
      if (length(R_tp1mi) == 0){
        comp_checks = rep(1, H + 1)
      } else {
        # I want to build a matrix with all the possible allocations of the i-th
        # observation. Each row is the membership vector of the i-th obs. for one 
        # of such allocations. Hence it is a (H+1)x(H+1) diagonal matrix with 
        # diagonal equal to 1. Then to check the co-clustering of i-th obs. in each
        # of these H+1 possible allocations I should do 
        # tcrossprod(diag(1, H + 1), t(cbind(c_t[R_tp1, ], 0))), which is simply
        # t(cbind(c_t[R_tp1, ], 0)).
        tmp1 = t(cbind(c_t[R_tp1mi, , drop = FALSE], 0))
        # Compute i-th row of the co-clustering matrix at time t+1.
        tmp2 = tcrossprod(c_tp1[i, , drop = FALSE], c_tp1[R_tp1mi, , drop = FALSE])
        comp_checks = as.integer(apply(tmp1, 1, function(x) all(x == tmp2)))
      }
    }
  }
  
  return(comp_checks)
}
compatibility_check_cppwrapper = function(i, H, k, gamma_tp1, c_t, c_tp1){
  if (gamma_tp1[i] == 0){
    return(rep(1, H + 1))
  } else {
    compatibility_check_cpp(i-1, H, gamma_tp1, c_t, c_tp1)
  }
}

# From vector of labels to binary matrix and viceversa ----
vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- length(unique(clust_lab))
  M <- matrix(0, V, H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

mat2vec <- function(clust_matrix){
  # in: binary nC x "number of clusters" matrix M s.t. M[v,h]=1{node v is in cluster h}
  # out: vector clust_lab of length nC s.t. clust_lab[v]=h if node v is in cluster h
  
  clust_lab = clust_matrix %*% c(1:ncol(clust_matrix)) %>% as.integer()
  return(clust_lab)
  
}

# Setting of control objects ----
set_ctr_mcmc = function(
    seed = 1234, nburnin = 1000, nchain = 5000, print_step = 100, ncl.init = NULL,
    verbose = c("0", "1", "2")
){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(nburnin) & nburnin >= 0)
  stopifnot(is.numeric(nchain) & nchain >= 1)
  stopifnot(is.numeric(print_step) & print_step >= 1)
  stopifnot(is.numeric(ncl.init) & ncl.init >= 1 & ncl.init <= n)
  
  verbose = match.arg(verbose)
  
  list(seed = floor(seed), nburnin = floor(nburnin), nchain = floor(nchain), 
       print_step = floor(print_step), ncl.init = floor(ncl.init),
       verbose = as.integer(verbose))
}

set_ctr_save = function(
    save = FALSE, filepath = "", 
    filename = paste("res_", Sys.Date(), ".RDS", sep = "")
){
  stopifnot(is.logical(save))
  stopifnot(is.character(filepath))
  stopifnot(is.character(filename))
  
  list(save = save, filepath = filepath, filename = filename)
}


set_ctr_mcmc_par = function(
    seed = 1234, nburnin = 1000, nchain = 5000, print_step = 100,
    verbose = c("0", "1")
){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(nburnin) & nburnin >= 0)
  stopifnot(is.numeric(nchain) & nchain >= 1)
  stopifnot(is.numeric(print_step) & print_step >= 1)
  
  verbose = match.arg(verbose)
  
  list(seed = floor(seed), nburnin = floor(nburnin), nchain = floor(nchain), 
       print_step = floor(print_step), verbose = as.integer(verbose))
}

set_ctr_save_par = function(
    save = FALSE, filepath = "", 
    filename = paste("par_", Sys.Date(), ".RDS", sep = "")
){
  stopifnot(is.logical(save))
  stopifnot(is.character(filepath))
  stopifnot(is.character(filename))
  
  list(save = save, filepath = filepath, filename = filename)
}

set_ctr_alpha = function(
    fix_alpha.flag = FALSE, fix_alpha.value = NULL
){
  stopifnot(is.logical(fix_alpha.flag))
  stopifnot(fix_alpha.flag & is.numeric(fix_alpha.value) & (fix_alpha.value > 0 & fix_alpha.value < 1))
  
  list(fix_alpha.flag = fix_alpha.flag, fix_alpha.value = fix_alpha.value)
}


# Distributions ----
dhamming = function(x, mu, sigma, m, logscale = FALSE){
  # x: data (vector of length p, where p is the number of attributes)
  # mu: center (vector of length p)
  # sigma: scale (vector of length p)
  # m: number of attributes (vector of length p)
  # logscale: if TRUE, return log-likelihood
  
  denom = 1 + (m - 1) / exp(1/sigma)
  lognum = - (x != mu) / sigma
  if (logscale) {
    out = sum(lognum) - sum(log(denom))
  } else {
    out = exp(lognum) / denom
  }
  
  return(out)
}

dhamming.multicluster = function(x, mu, sigma, m, logscale = FALSE){
  # x: data (vector of length p, where p is the number of attributes)
  # mu: center for every cluster (matrix of size H x p, where H is nmbr of clusters)
  # sigma: scale for every cluster (matrix of size H x p, where H is nmbr of clusters)
  # m: number of attributes (vector of length p)
  # logscale: if TRUE, return log-likelihood
  
  if (any(dim(mu) != dim(sigma))) {
    stop("mu and sigma must have the same dimensions")
  }
  
  H = nrow(mu)
  p = ncol(mu)
  
  M = matrix(m, nrow = H, ncol = p, byrow = TRUE)
  X = matrix(x, nrow = H, ncol = p, byrow = TRUE)
  
  logdenom = log(1 + (M - 1) / exp(1/sigma))
  lognum = - (X != mu) / sigma
  
  logout = rowSums(lognum) - rowSums(logdenom)
  if (logscale) {
    out = logout
  } else {
    out = exp(logout)
  }
  
  return(out)
}

urn_GN <- function(sizecl, gamma){
  # sizecl: vector of length H (nmbr of existing clusters) with cluster sizes
  # gamma: Gnedin parameter
  
  H = length(sizecl)
  oldcl = (sizecl+1)*(sum(sizecl)-H+gamma)
  newcl = H^2-H*gamma
  out = c(oldcl, newcl)
  
  return(out)
}

urn_DP <- function(sizecl, M){
  return(c(sizecl, M))
}

# Probability of having h (occupied) clusters given n observations assuming Gnedin
HGnedin <- function(n, h, gamma=0.5){
  exp(lchoose(n, h) + lgamma(h-gamma) - lgamma(1-gamma) + log(gamma) + lgamma(n+ gamma - h) - lgamma(n +gamma))
}

HbarGnedin <- function(h, gamma=0.5){
  lpoch = lgamma(1 - gamma + h - 1) - lgamma(1-gamma)
  logout = log(gamma) + lpoch - lfactorial(h) 
  exp(logout)
}

# Expected number of clusters DP
E_Nclusters_DP <- function(n, M){
  # n: number of observations
  # M: concentration parameter
  
  if (M == 0){
    return(0)
  } else {
    return(sum(M / (M + 1:n - 1)))
  }
}


# Simulate data for simstudy ----
simdata = function(seed){
  set.seed(seed)
  
  # Dimensions
  n = 100
  p = 5
  nT = 8
  ncl = 6
  list_values = list(list(c(1, 2), c(2, 3, 5, 7), c(3, 4, 5, 7), c(3, 4, 5, 6), c(4, 6)),
                     list(c(1, 2), c(2, 3, 5, 7), c(3, 4, 5, 7), c(3, 4, 5, 6), c(4, 6)),
                     list(c(1, 2), c(2, 3, 5, 7), c(3, 4, 5, 7), c(3, 4, 5, 6), c(4, 6)),
                     list(c(1, 2), c(2, 3, 5, 7), c(3, 4, 5, 7), c(3, 4, 5, 6), c(4, 6)),
                     list(c(1, 2, 9), c(2, 3, 5, 7, 8, 9), c(3, 4, 5, 7, 8, 9), c(3, 4, 5, 6, 8, 9), c(4, 6, 9)),
                     list(c(1, 2), c(2, 3, 5, 7, 8), c(3, 4, 5, 7, 8), c(3, 4, 5, 6, 8), c(4, 6)),
                     list(c(1, 2), c(2, 3, 5, 7, 8, 10), c(3, 4, 5, 7, 8, 10), c(3, 4, 5, 6, 8, 10), c(4, 6, 10)),
                     list(c(1, 2), c(2, 5, 7, 8, 10), c(3, 4, 5, 8, 10), c(3, 4, 5, 8, 10), c(4, 6, 10)))
  
  # Store objects
  C = array(NA, dim = c(n, nT), dimnames = list("Subject" = 1:n, "Year" = 1:nT))
  Theta = array(NA, dim = c(ncl, 5, 8), dimnames = list("Cluster" = 1:ncl, "Age" = 1:p, "Year" = 1:nT)) 
  Y = array(NA, dim = c(n, p, nT), dimnames = list("Subject" = 1:n, "Age" = 1:p, "Year" = 1:nT))
  to_move = list()
  
  # Time 1
  C[ , 1] = c(rep(1, 30), rep(2, 15), rep(3, 25), rep(4, 10), rep(5, 5), rep(6, 15))
  Theta[ , , 1] = matrix(c(1, 2, 3, 4, 4,
                           2, 2, 3, 3, 4,
                           2, 2, 5, 4, 6,
                           1, 5, 5, 6, 6,
                           1, 5, 3, 3, 6,
                           1, 7, 3, 3, 4), 
                         nrow = ncl, ncol = p, byrow = "T")
  
  # Time 2
  Theta[ , , 2] = Theta[ , , 1]
  to_move[[2]] = list(sample(which(C[ , 1] == 2), 5, replace = FALSE),
                      sample(which(C[ , 1] == 6), 5, replace = FALSE))
  C[ , 2] = C[ , 1]
  C[to_move[[2]][[1]], 2] = 1
  C[to_move[[2]][[2]], 2] = 5
  Y[ , , 2] = Theta[ , , 2][C[ , 2], ]
  
  # Time 3
  Theta[ , , 3] = Theta[ , , 2]
  Theta[2, 4, 3] = Theta[1, 4, 3] 
  Theta[5, , 3] = c(1, 7, 3, 3, 6)
  Theta[6, , 3] = NA
  to_move[[3]] = list(which(C[ , 2] == 6),
                      sample(which(C[ , 2] == 3), 5, replace = FALSE),
                      sample(which(C[ , 2] == 2), 5, replace = FALSE))
  C[ , 3] = C[ , 2]
  C[to_move[[3]][[1]], 3] = 5
  C[to_move[[3]][[2]], 3] = 4
  C[to_move[[3]][[3]], 3] = 1
  
  # Time 4
  Theta[ , , 4] = Theta[ , , 3]
  Theta[3, 2, 4] = Theta[4, 2, 4]
  Theta[2, , 4] = NA
  Theta[5, 1, 4] = 2; Theta[5, 3, 4] = 5
  C[ , 4] = C[ , 3]
  to_move[[4]] = list(which(C[ , 3] == 2),
                      sample(which(C[ , 2] == 3), 5, replace = FALSE))
  C[to_move[[4]][[1]], 4] = 1
  C[to_move[[4]][[2]], 4] = 4
  
  # Time 5
  Theta[ , , 5] = Theta[ , , 4]
  Theta[2, , 5] = c(1, 8, 8, 8, 4)
  Theta[6, , 5] = c(9, 9, 9, 9, 9)
  C[ , 5] = C[ , 4]
  to_move[[5]] = list(sample(which(C[ , 4] == 3), 10, replace = FALSE),
                      # The populations that move from cluster 1 to cluster 2 must
                      # be sampled from those who were in cluster 1 not only at 
                      # time 4, but also at time 3.
                      sample(which(C[ , 4] == 1 & C[ , 3] == 1), 10, replace = FALSE))
  C[to_move[[5]][[1]][1:4], 5] = 6
  C[to_move[[5]][[1]][5:10], 5] = 4
  C[to_move[[5]][[2]], 5] = 2
  
  # Time 6
  Theta[ , , 6] = Theta[ , , 5]
  Theta[3, , 6] = NA
  Theta[6, , 6] = NA
  Theta[4, 1, 6] = 2
  C[ , 6] = C[ , 5]
  to_move[[6]] = list(c(which(C[ , 5] == 3), which(C[ , 5] == 6)),
                      sample(which(C[ , 5] == 1), 10, replace = FALSE),
                      sample(which(C[ , 5] == 4), 12, replace = FALSE))
  C[to_move[[6]][[1]], 6] = 4
  C[to_move[[6]][[2]], 6] = 2
  C[to_move[[6]][[3]], 6] = 5
  
  # Time 7
  Theta[ , , 7] = Theta[ , , 6]
  Theta[1, 5, 7] = 10
  Theta[2, 4:5, 7] = 10
  Theta[4, 4:5, 7] = 10
  Theta[5, 4:5, 7] = 10
  
  C[ , 7] = C[ , 6]
  to_move[[7]] = list(sample(which(C[ , 6] == 4), 12, replace = FALSE),
                      sample(which(C[ , 6] == 1), 10, replace = FALSE))
  C[to_move[[7]][[1]], 7] = 5
  C[to_move[[7]][[2]], 7] = 2
  
  # Time 8
  Theta[ , , 8] = Theta[ , , 7]
  Theta[5, , 8] = NA
  Theta[2, 2:3, 8] = 10
  Theta[1, 4, 8] = 10
  
  C[ , 8] = C[ , 7]
  to_move[[8]] = list(which(C[ , 7] == 5),
                      sample(which(C[ , 7] == 1), 5, replace = FALSE))
  C[to_move[[8]][[1]], 8] = 4
  C[to_move[[8]][[2]], 8] = 2
  
  
  # Y
  prtop = 0.95
  prother = 1-prtop
  for (t in 1:nT){
    
    for (cl in sort(unique(C[ , t]))){
      
      idx = which(C[ , t] == cl)
      
      for (j in 1:p){
        val = list_values[[t]][[j]]
        m = length(val)
        pr = rep(0, max(val))
        top = Theta[cl, j, t]
        others = val[val != top]
        pr[top] = prtop
        if (t == 5 & cl != 6){
          pr[others] = rep(prother/(m-2), m-1)
          pr[9] = 1e-6
        } else if (t== 5 & cl == 6){
          pr[top] = 0.99
          pr[others] = 1e-6
        } else {
          pr[others] = rep(prother/(m-1), m-1)
        }
        Y[idx, j, t] = sample(val, length(idx), prob = pr[pr!=0], replace = T)
      }
      
    }
    
  }
  
  randomorder = sample(1:n, n, replace = FALSE)
  
  out = list(Y = Y[randomorder, , ], C = C[randomorder, ], Theta = Theta, to_move = to_move,
             orig = list(Y = Y, C = C),
             info = list(n = n, p = p, nT = nT, ncl = ncl,
                         list_values = list_values))
  
  return(out)
}

# Plots ----
plot_tHMM_flow = function(ppe, tosort = TRUE, show.text = TRUE,
                          alpha_flow = 0.4, alpha_stratum = 1, add2title = ""){
  library(ggalluvial)
  library(grafify)
  
  if (tosort){
    ppe = ppe %>% seqsort()
  }
  
  # clnms = colnames(ppe)
  colnames(ppe) = 1:ncol(ppe)
  
  max_ncl = max(ppe)
  
  clust_flow = ppe %>% 
    as.data.frame() %>% tibble::rownames_to_column(var = "obs") %>% dplyr::mutate(freq = 1) %>% 
    reshape2::melt(id.vars = c("obs", "freq"), variable.name = "year", value.name = "cluster") %>% 
    dplyr::mutate(cluster = factor(cluster)) 
  
  
  plt_tks = c("1"="2000", "2"="2001", "3"="2002", "4"="2003", "5"="2004",
              "6"="2005", "7"="2006", "8"="2007", "9"="2008", "10"="2009",
              "11"="2010", "12"="2011", "13"="2012", "14"="2013", "15"="2014",
              "16"="2015", "17"="2016", "18"="2017", "19"="2018", "20"="2019",
              "21"="2020", "22"="2021")
  
  plt_thm = theme_grey() + 
    theme(legend.position = "none", 
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank())
  
  plt = clust_flow %>%
    ggplot(map = aes(x = year, y = freq, stratum = cluster, 
                     alluvium = obs, fill = cluster, label = cluster)) +
    geom_flow(alpha = alpha_flow) + 
    geom_stratum(alpha = alpha_stratum, color = "white") +
    scale_x_discrete(expand = c(0.025, 0.025), labels = plt_tks) +
    scale_y_continuous(expand = c(0.025, 0.025)) +
    labs(x = "Age class", y = "Causes of death") +
    ggtitle(paste("Cluster evolution over age classes", add2title)) + 
    plt_thm# + theme(axis.text.x = element_blank(),
            #        axis.ticks.x = element_blank())
  
  if (show.text){
    plt = plt + geom_text(stat = "stratum", size = 3)
  }
  
  return(plt)
}
