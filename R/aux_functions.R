#Fast (weighted) mean, optionally with subset
.wtd_mean <- function(x, w = NULL, subset = NULL) {
  if (is.null(subset)) {
    if (is.null(w)) {
      sum(x) / length(x)
    }
    else {
      sum(x * w) / sum(w)
    }
  }
  else {
    x <- x[subset]
    if (is.null(w)) {
      sum(x) / length(x)
    }
    else {
      w <- w[subset]
      sum(x * w) / sum(w)
    }
  }
}

.wtd_sd <- function(x, w = NULL, subset = NULL) {
  if (is.null(subset)) {
    if (is.null(w)) {
      sqrt(sum((x - .wtd_mean(x))^2)/(length(x) - 1))
    }
    else {
      sum_w <- sum(w)
      sqrt((sum_w / (sum_w^2 - sum(w^2))) * sum(w * (x - .wtd_mean(x, w))^2))
    }
  }
  else {
    x <- x[subset]
    if (is.null(w)) {
      sqrt(sum((x - .wtd_mean(x))^2)/(length(x) - 1))
    }
    else {
      w <- w[subset]
      sum_w <- sum(w)
      sqrt((sum_w / (sum_w^2 - sum(w^2))) * sum(w * (x - .wtd_mean(x, w))^2))
    }
  }
}

# Create a list indexing position of unlist elements of list x
.list_ind <- function(x) {
  l <- lengths(x)
  lapply(seq_along(x), function(i) {
    sum(l[seq_len(i - 1L)]) + seq_len(l[i])
  })
}

#Checks if a given family specification is okay
.okay_family <- function(family) {
  if (is.character(family)) {
    if (length(family) != 1 || anyNA(family)) return(FALSE)
    if (family %in% c("negbin", "negative.binomial", "Negative Binomial")) return(TRUE)
    family <- get(family, mode = "function", envir = parent.frame(2))
  }
  if (is.function(family)) {
    family <- family()
  }
  
  !is.null(family$family) && is.function(family$variance) &&
    is.function(family$linkinv)
}

.ordinal <- function(x) {
  if (length(x) == 0L || !is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }
  
  if (length(x) > 1L) {
    out <- setNames(vapply(x, .ordinal, character(1L)), names(x))
    return(out)
  }
  
  x0 <- abs(x)
  out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)), 
                           `1` = "st", `2` = "nd", `3` = "rd", "th"))
  if (x < 0) {
    out <- sprintf("-%s", out)
  }
  
  setNames(out, names(x))
}

.firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

.fill_vec <- function(u, ...) {
  out <- setNames(rep(0, length(u)), u)
  
  for (i in seq_len(...length())) {
    x <- ...elt(i)
    out[names(x)] <- x
  }
  
  out
}

.fill_mat <- function(u, x) {
  out <- matrix(0, nrow = length(u), ncol = ncol(x),
                dimnames = list(u, colnames(x)))
  
  out[rownames(x),] <- x
  
  out
}