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
    seed = 1234, nburnin = 1000, nchain = 5000, print_step = 100, 
    verbose = c("0", "1", "2")
){
  stopifnot(is.numeric(seed))
  stopifnot(is.numeric(nburnin) & nburnin >= 0)
  stopifnot(is.numeric(nchain) & nchain >= 1)
  stopifnot(is.numeric(print_step) & print_step >= 1)
  
  verbose = match.arg(verbose)
  
  list(seed = floor(seed), nburnin = floor(nburnin), nchain = floor(nchain), 
       print_step = floor(print_step), verbose = as.integer(verbose))
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

urn_GN <- function(sizecl, eta){
  # sizecl: vector of length H (nmbr of existing clusters) with cluster sizes
  # eta: Gnedin parameter
  
  H = length(sizecl)
  oldcl = (sizecl+1)*(sum(sizecl)-H+eta)
  newcl = H^2-H*eta
  out = c(oldcl, newcl)
  
  return(out)
}

# Probability of having h (occupied) clusters given n observations assuming Gnedin
HGnedin <- function(n, h, gamma=0.5){
  exp(lchoose(n, h) + lgamma(h-gamma) - lgamma(1-gamma) + log(gamma) + lgamma(n+ gamma - h) - lgamma(n +gamma))
}