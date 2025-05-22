#Joint covariance matrix of coefficients across multiple models. Requires same units in all models,
#use NA weights to subset (e.g., weights of 1 for present and NA for absent). Should give
#same results as M-estimation (HC0 vcov). Cluster-robust vcov available.
vcovSUEST <- function(fits, cluster = NULL) {
  
  if (is_null(names(fits)) || anyDuplicated(names(fits)) > 0) {
    names(fits) <- paste0("fit_", seq_along(fits))
  }
  
  rn <- lapply(fits, function(x) rownames(naresid(x$na.action, model.matrix(x))))
  u <- unique(unlist(rn, use.names = FALSE))
  
  if (is_null(cluster)) {
    cluster <- list(u)
  }
  else {
    if (inherits(cluster, "formula")) {
      cluster_tmp <- expand.model.frame(fits[[1L]], cluster, na.expand = TRUE)
      cluster_f <- as.list(model.frame(cluster, cluster_tmp, na.action = na.pass))
    }
    else {
      cluster_f <- as.list(as.data.frame(cluster))
    }
    
    for (i in seq_along(cluster_f)) {
      cluster_f[[i]] <- .fill_vec(u, setNames(cluster_f[[i]], rn[[1L]]))
    }
    
    if (inherits(cluster, "formula")) {
      for (f in seq_along(fits)[-1L]) {
        cluster_tmp <- expand.model.frame(fits[[f]], cluster, na.expand = TRUE)
        cluster_fi <- as.list(model.frame(cluster, cluster_tmp, na.action = na.pass))
        
        
        for (i in seq_along(cluster_f)) {
          cluster_f[[i]] <- .fill_vec(u, setNames(cluster_fi[[i]], rn[[f]]))
        }
      }
    }
    
    cluster <- c(cluster_f, list(u))
  }
  
  p <- length(cluster)
  if (p > 1L) {
    clu <- lapply(seq_len(p), function(i) utils::combn(seq_len(p), i, simplify = FALSE))
    clu <- unlist(clu, recursive = FALSE)
    sign <- vapply(clu, function(i) (-1)^(length(i) + 1), numeric(1L))
    paste_ <- function(...) paste(..., sep = "_")
    for (i in (p + 1L):length(clu)) {
      cluster <- c(cluster, list(Reduce(paste_, unclass(cluster[clu[[i]]]))))
    }
  }
  else {
    clu <- list(1)
    sign <- 1
  }
  
  #Small sample adjustment (setting cadjust = TRUE in vcovCL)
  g <- lapply(seq_along(fits), function(i) {
    vapply(seq_along(clu), function(u) {
      cu <- cluster[[u]]
      length(unique(cu[-fits[[i]]$na.action]))
    }, numeric(1L))
  })
  
  breads <- lapply(fits, function(f) .bread(f) / nobs(f))
  
  ef <- lapply(seq_along(fits), function(i) {
    ef_i <- sandwich::estfun(fits[[i]])
    
    if (anyNA(ef_i)) {
      ef_i[is.na(ef_i)] <- 0
    }
    
    rownames(ef_i) <- rn[[i]]
    
    ef_i <- .fill_mat(u, ef_i)
    
    lapply(seq_along(clu), function(u) {
      rowsum(ef_i, cluster[[u]], reorder = FALSE)
    })
  })
  
  coef_lengths <- vapply(breads, ncol, numeric(1L))
  
  coef_inds <- split(seq_len(sum(coef_lengths)),
                     rep(seq_along(coef_lengths), coef_lengths))
  
  #VCOV matrix to be returned
  V <- matrix(NA_real_, nrow = sum(coef_lengths), ncol = sum(coef_lengths))
  dimnames(V) <- rep.int(list(unlist(lapply(seq_along(fits), function(i) {
    paste(names(fits)[i], colnames(breads[[i]]), sep = "_")
  }))), 2L)
  
  for (i in seq_along(fits)) {
    ind_i <- coef_inds[[i]]
    
    #Usual within-model HC0 vcov
    S <- 0
    for (u in seq_along(clu)) {
      adj <- g[[i]][u] / (g[[i]][u] - 1)
      S <- S + sign[u] * adj * crossprod(ef[[i]][[u]])
    }
    
    V[ind_i, ind_i] <- breads[[i]] %*% S %*% breads[[i]]
    
    for (j in seq_along(fits)[-seq_len(i)]) {
      ind_j <- coef_inds[[j]]
      
      S <- 0
      for (u in seq_along(clu)) {
        adj <- sqrt(g[[i]][u] / (g[[i]][u] - 1)) * sqrt(g[[j]][u] / (g[[j]][u] - 1))
        
        S <- S + sign[u] * adj * crossprod(ef[[i]][[u]], ef[[j]][[u]])
      }
      
      #between-model vcov components
      V[ind_i, ind_j] <- breads[[i]] %*% S %*% breads[[j]]
      V[ind_j, ind_i] <- t(V[ind_i, ind_j])
    }
  }
  
  V
}

#Quickly get bread matrix; for non-lm and non-glm objects, uses sandwich::bread()
.bread <- function(x) {
  if (!class(x)[1L] %in% c("lm", "glm")) {
    return(sandwich::bread(x))
  }
  
  p <- x$rank
  
  if (p == 0) {
    return(matrix(NA_real_, 0L, 0L))
  }
  
  Qr <- x$qr
  p1 <- seq_len(p)
  coef.p <- x$coefficients[Qr$pivot[p1]]
  cov.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  dimnames(cov.unscaled) <- list(names(coef.p), names(coef.p))
  
  out <- cov.unscaled * (p + x$df.residual)
  
  if (class(x)[1L] == "glm" && !substr(x$family$family, 1L, 17L) %in% c("poisson", "binomial", "Negative Binomial")) {
    ww <- weights(x, "working")
    wres <- as.vector(residuals(x, "working")) * ww
    dispersion <- sum(wres^2, na.rm = TRUE) / sum(ww, na.rm = TRUE)
    out <- out * dispersion
  }
  
  out
}