# MISCELLANEOUS FUNCTIONS ----
## From vector of labels to binary matrix and viceversa ----
vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0, V, H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

mat2vec <- function(clust_matrix){
  # in: binary nC x "number of clusters" matrix M s.t. M[v,h]=1{node v is in cluster h}
  # out: vector clust_lab of length nC s.t. clust_lab[v]=h if node v is in cluster h
  
  clust_lab = clust_matrix %*% c(1:ncol(clust_matrix))
  return(clust_lab)
  
}


# POSTPROCESSING FUNCTIONS ----
# Re-labeling ----
# Sort the labels of a partition by cluster dimension
dimsort = function(u){
  tab = table(u)
  idx = names(sort(tab, decreasing = TRUE))
  lab = setNames(seq_along(idx), idx)
  out = lab[as.character(u)]
  return(out)
}

# Sort the labels of a partition by matching the labels of a reference partition
refsort.cri = function(u, v){
  # Partition dimensions
  nu = max(u)
  nv = max(v)
  
  # Membership matrices
  U = vec2mat(u)
  V = vec2mat(v)
  
  # Matrix of cross-membership counts
  mat = crossprod(V, U)
  tab = table(u)
  
  # Used and unused labels, and empty ordering vector
  used = c()
  unused = c(1:nu)
  order = rep(NA, length = nu)
  
  # Order initialization
  h = 1
  idx = which.max(mat[h, ])
  used = c(used, idx)
  unused = setdiff(unused, idx)
  order[h] = idx
  
  # Sorting loop
  for(h in 2:nu){
    # Positions of the most frequent co-clustering values in decreasing order
    if (h <= nv) {
      tmp = order(mat[h, ], decreasing = TRUE)
    } else {
      tmp = order(tab, decreasing = TRUE)
    }
    # Exclude the labels that are already assigned
    idx = setdiff(tmp, used)[1]
    # Update used and unused labels
    used = c(used, idx)
    unused = setdiff(unused, idx)
    # Set the h-th entry of the order vector
    order[h] = idx
  }
  
  # Re-order the columns of the membership matrix and
  out = mat2vec(U[, order])
  
  # Return the re-labeled partition
  return(out)
}

# Sort the labels of a partition by matching the labels of a reference partition
refsort = function(u, v){
  # Partition dimensions
  nu = max(u)
  nv = max(v)
  
  # Membership matrices
  U = vec2mat(u)
  V = vec2mat(v)
  
  # Matrix of cross-membership counts
  mat = crossprod(V, U)
  tab = table(u)
  
  # Used and unused labels, and empty ordering vector
  used = c()
  unused = c(1:nu)
  order = rep(NA, length = nu)
  
  # Order initialization
  h = 1
  idx = which.max(mat[h, ])
  used = c(used, idx)
  unused = setdiff(unused, idx)
  order = idx
  
  # Sorting loop
  for(h in 2:nu){
    # Positions of the most frequent co-clustering values in decreasing order
    if (h <= nv) {
      tmp = order(mat[h, ], decreasing = TRUE)
    } else {
      tmp = order(tab, decreasing = TRUE)
    }
    
    idx = tmp[1]
    
    if (idx %in% used){
      next
    } else {
      # Update used and unused labels
      used = c(used, idx)
      unused = setdiff(unused, idx)
      # Set the h-th entry of the order vector
      order = c(order, idx)
    }
  }
  
  order = c(order, unused)
  
  # Re-order the columns of the membership matrix and
  out = mat2vec(U[, order])
  
  # Return the re-labeled partition
  return(out)
}

# Sort the labels of a sequence of partitions by sequential applying refsort 
seqsort = function(partition, reference = TRUE, dimsort = TRUE){
  
  if (dimsort){
    # Sort all the partitions by cluster dimension 
    out = apply(partition, 2, dimsort)
    dimnames(out) = dimnames(partition)
  } else {
    # If not sorting by cluster dimension, just copy the input
    out = partition
  }
  
  # Relabel each partition by sequential reference sorting 
  if(reference){
    for(t in 2:ncol(out)){
      out[,t] = refsort(out[,t], out[,t-1])
    }
  }
  return(out)
}


moda = function(x){
  tx = table(x)
  return(as.numeric(dimnames(tx)$x[which.max(tx)]))
}
