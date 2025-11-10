# # ZOO ----
# m = c(rep(2, 15), 6)
# u = c(rep(6, 15), 3)
# v = c(rep(0.25, 15), 0.25)
# nsig = 10^4
# S  = sapply(1:16, function(x) rhyper_sig2(n = nsig, c = u[x], d = v[x], m[x]), simplify = T)
# 
# Gmax <- 1 - prod(1 / m)
# 
# gini_values = {
#   expS <- (exp(2 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m)))) /
#     (exp(1 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m))))^2
#   expS[S < 0.01] <- 1
#   tmp <- apply(expS, 1, prod)
#   gini_values <- (1 - tmp) / Gmax
# 
#     # expS = t( (t(exp(2/S)) + (m - 1)) / (t(exp(1/S)) + (m - 1))^2 )
#   # expS[S < 0.01] = 1
#   # tmp = apply(expS, 1, prod)
#   # (1 - tmp)/(1 - prod(1/m))
# }
# 
# par(mfrow = c(1, 1))
# hist(gini_values, freq = F, breaks = 21, main = "Normalized Gini", xlab = "Gini values", ylab = "Density", col = "lightblue", border = "black")
# 
# 
# 
# # USPS ----
# m = rep(6, 256)
# u = rep(3, 256)
# v = rep(0.5, 256)
# nsig = 10^3
# S  = sapply(1:256, function(x) rhyper_sig2(n = nsig, c = u[x], d = v[x], m[x]), simplify = T)
# 
# Gmax <- 1 - prod(1 / m)
# 
# gini_values = {
#   expS <- (exp(2 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m)))) /
#     (exp(1 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m))))^2
#   expS[S < 0.01] <- 1
#   tmp <- apply(expS, 1, prod)
#   gini_values <- (1 - tmp) / Gmax
#   
#   # expS = t( (t(exp(2/S)) + (m - 1)) / (t(exp(1/S)) + (m - 1))^2 )
#   # expS[S < 0.01] = 1
#   # tmp = apply(expS, 1, prod)
#   # (1 - tmp)/(1 - prod(1/m))
# }
# 
# par(mfrow = c(1, 1))
# hist(gini_values, freq = F, breaks = 21, main = "Normalized Gini", xlab = "Gini values", ylab = "Density", col = "lightblue", border = "black")
# 
# # SOY BEAN ----
# m = c(as.matrix(read.table("src/ArgientoPaci/data/soybean_mm.txt", quote="\"", comment.char="")))[-36]
# u = c(6, 5, 4.5, 4.25, 3, 3)[m-1]
# v = rep(0.25, 35)
# nsig = 10^3
# S  = sapply(1:35, function(x) rhyper_sig2(n = nsig, c = u[x], d = v[x], m[x]), simplify = T)
# 
# Gmax <- 1 - prod(1 / m)
# 
# gini_values = {
#   expS <- (exp(2 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m)))) /
#     (exp(1 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m))))^2
#   expS[S < 0.01] <- 1
#   tmp <- apply(expS, 1, prod)
#   gini_values <- (1 - tmp) / Gmax
#   
#   # expS = t( (t(exp(2/S)) + (m - 1)) / (t(exp(1/S)) + (m - 1))^2 )
#   # expS[S < 0.01] = 1
#   # tmp = apply(expS, 1, prod)
#   # (1 - tmp)/(1 - prod(1/m))
# }
# 
# par(mfrow = c(1, 1))
# hist(gini_values, freq = F, breaks = 21, main = "Normalized Gini", xlab = "Gini values", ylab = "Density", col = "lightblue", border = "black")


# SOY BEAN 2 ----
m = 2:7
u = c(6.5, 5.5, 4.5, 4.25, 4, 4)
v = rep(0.25, 6)
val = cbind(m, u, v)
nsig = 10^4
S  = apply(val, 1, function(x) rhyper_sig2(n = nsig, c = x[2], d = x[3], m = x[1]), simplify = T)

Gmax <- 1 - 1 / m

gini_values = {
  expS = (exp(2 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m)))) /
    (exp(1 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m))))^2
  expS[S < 0.01] = 1
  gini_values = t( (1 - t(expS)) / Gmax )
}

par(mfrow = c(2, 3))
for (j in 1:ncol(gini_values)){
  hist(gini_values[ , j], freq = F, breaks = 21, main = paste0("N-Gini, v=", u[j], ", w=", v[j]),
       xlab = "Gini values", ylab = "Density", col = "lightblue", border = "black")
}


# Figure S3 supplementary ----
m = 4
u = 5
v = 1
nsig = 10^3
S  = sapply(1, function(x) rhyper_sig2(n = nsig, c = u[x], d = v[x], m[x]), simplify = T)

Gmax <- 1 - prod(1 / m)

gini_values = {
  expS <- (exp(2 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m)))) /
    (exp(1 / S) + (matrix(rep(m - 1, each = nrow(S)), ncol = length(m))))^2
  expS[S < 0.01] <- 1
  tmp <- apply(expS, 1, prod)
  gini_values <- (1 - tmp) / Gmax
  
  # expS = t( (t(exp(2/S)) + (m - 1)) / (t(exp(1/S)) + (m - 1))^2 )
  # expS[S < 0.01] = 1
  # tmp = apply(expS, 1, prod)
  # (1 - tmp)/(1 - prod(1/m))
}

par(mfrow = c(1, 1))
hist(gini_values, freq = F, breaks = 21, main = "Normalized Gini", xlab = "Gini values", ylab = "Density", col = "lightblue", border = "black")
