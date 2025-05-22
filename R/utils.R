#Fast (weighted) mean, optionally with subset
.wtd_mean <- function(x, w = NULL, subset = NULL) {
  if (!is_null(subset)) {
    if (is_null(w)) {
      return(Recall(x[subset]))
    }
    
    return(Recall(x[subset], w = w[subset]))
  }
  
  if (is_null(w)) {
    sum(x) / length(x)
  }
  else {
    sum(x * w) / sum(w)
  }
}

.wtd_sd <- function(x, w = NULL, subset = NULL) {
  if (!is_null(subset)) {
    if (is_null(w)) {
      return(Recall(x[subset]))
    }
    
    return(Recall(x[subset], w = w[subset]))
  }
  
  if (is_null(w)) {
    sqrt(sum((x - .wtd_mean(x))^2) / (length(x) - 1))
  }
  else {
    sum_w <- sum(w)
    sqrt((sum_w / (sum_w^2 - sum(w^2))) * sum(w * (x - .wtd_mean(x, w))^2))
  }
}

.colMax <- function(x, na.rm = TRUE) {
  apply(x, 2L, max, na.rm = na.rm)
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
    if (length(family) != 1L || anyNA(family)) {
      return(FALSE)
    }
    
    if (family %in% c("negbin", "negative.binomial", "Negative Binomial")) {
      return(TRUE)
    }
    
    family <- get(family, mode = "function", envir = parent.frame(2L))
  }
  
  if (is.function(family)) {
    family <- family()
  }
  
  !is_null(family$family) &&
    is.function(family$variance) &&
    is.function(family$linkinv)
}

.ordinal <- function(x) {
  if (is_null(x) || !is.numeric(x)) {
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
  substr(x, 1L, 1L) <- toupper(substr(x, 1L, 1L))
  x
}

.fill_vec <- function(u, ...) {
  out <- setNames(rep.int(0, length(u)), u)
  
  for (i in seq_len(...length())) {
    x <- ...elt(i)
    out[names(x)] <- x
  }
  
  out
}

.fill_mat <- function(u, x) {
  out <- matrix(0, nrow = length(u), ncol = ncol(x),
                dimnames = list(u, colnames(x)))
  
  out[rownames(x), ] <- x
  
  out
}

#More informative and cleaner version of base::match.arg(). Uses chk.
.match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    chk::err("no argument was supplied to `match_arg()`")
  }
  
  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)
  
  if (missing(choices)) {
    sysP <- sys.parent()
    formal.args <- formals(sys.function(sysP))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }
  
  if (is_null(arg)) {
    return(choices[1L])
  }
  
  if (several.ok) {
    chk::chk_character(arg, x_name = .add_quotes(arg.name, "`"))
  }
  else {
    chk::chk_string(arg, x_name = .add_quotes(arg.name, "`"))
    if (identical(arg, choices)) {
      return(arg[1L])
    }
  }
  
  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    chk::err(sprintf("the argument to `%s` should be %s%s",
                     arg.name,
                     ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                     .word_list(choices, and.or = "or", quotes = 2)))
  i <- i[i > 0L]
  
  choices[i]
}

.word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  
  word.list <- setdiff(word.list, c(NA_character_, ""))
  
  if (is_null(word.list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }
  
  word.list <- .add_quotes(word.list, quotes)
  
  L <- length(word.list)
  
  if (L == 1L) {
    out <- word.list
    if (is.are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }
  
  if (is_null(and.or) || isFALSE(and.or)) {
    out <- toString(word.list)
  }
  else {
    and.or <- .match_arg(and.or, c("and", "or"))
    
    if (L == 2L) {
      out <- sprintf("%s %s %s",
                     word.list[1L],
                     and.or,
                     word.list[2L])
    }
    else {
      out <- sprintf("%s, %s %s",
                     toString(word.list[-L]),
                     and.or,
                     word.list[L])
    }
  }
  
  if (is.are) {
    out <- sprintf("%s are", out)
  }
  
  attr(out, "plural") <- TRUE
  
  out
}

.add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }
  
  if (isTRUE(quotes)) {
    quotes <- '"'
  }
  
  if (chk::vld_string(quotes)) {
    quotes_rev <- {
      if (nchar(quotes) == 1L) quotes
      else vapply(lapply(strsplit(quotes_rev, NULL), rev),
                  paste, character(1L), collapse = "")
    }
    
    return(paste0(quotes, x, quotes_rev))
  }
  
  if (!chk::vld_count(quotes) || quotes > 2) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }
  
  if (quotes == 0) {
    return(x)
  }
  
  x <- {
    if (quotes == 1) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }
  
  x
}

.list_modify <- function(old, new) {
  new_names <- names(new)
  new_names <- new_names[nzchar(new_names)]
  
  for (v in new_names) {
    old[[v]] <- new[[v]]
  }
  
  old
}

is_null <- function(x) {
  length(x) == 0L
}