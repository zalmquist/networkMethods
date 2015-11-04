#' Block Fit
#'
#'\code{block.fit} returns an object of class "blockmodel."
#'
#' @param g  a network to be analyzed
#' @param model  a vector or matrix specifying the blockmodel to be fit.  If passed
#'     as a vector, the argument should be given column-wise.  For instance,
#'     two-class blockmodels are fit by passing a vector of the form
#'     c(b11,b21,b12,b22), where bij is the value of the i,j cell in the
#'     blockmodel.  Usually, these values are 1 or 0, although NAs are allowed
#'     (and have the effect of removing the indicated cells from consideration
#'     when evaluating goodness of fit).
#' @param measure  "pearson" for the standard correlation coefficent (may produce
#'     warnings, which are harmless) or "negdist" for the inverse absolute
#'     distance.
#' @param diag  should the diagonal be included?
#' @param   iter  number of annealing iterations to employ.
#' @param   temp.init initial temperature for the annealer.
#' @param   cool  cooling factor for the annealer.
#' @param   hill.climb.refine  should the annealing solution be further refined using
#' @param     hill-climbing?  (Can be slow, but guarantees convergence to a local optimum.)
#'     optimum.
#' @param   seed  an optional vector of initial group memberships.
#' @param   verbose  should the function provide diagnostic output?
#'
#' @return An object of class "blockmodel," with the additional list element block.gof.
#'   (Note that this object can be printed, etc. using the standard methods
#'   for this class.)
#'
#' @author Carter T. Butts, University of California, Irvine
#'   @examples
#'   require(sna)
#'   data(kaptail.ins)
#'   plot(kaptail.ins)
#' # Let's try fitting some blockmodels.  Here are several variants on C/P:
#' kb1<-block.fit(kaptail.ins,c(1,1,1,0))   # Core w/in,out ties
#' kb2<-block.fit(kaptail.ins,c(1,1,0,0))   # Core w/in ties
#' kb3<-block.fit(kaptail.ins,c(1,0,1,0))   # Core w/out ties
#' kb4<-block.fit(kaptail.ins,c(1,0,0,0))   # Isolated core
#' kb5<-block.fit(kaptail.ins,c(1,NA,NA,0)) # Ignore core/periphery interactions
#'
#' # Examine the output directly
#' kb1                                       # These are "blockmodel" objects
#' summary(kb2)
block.fit<-function (g, model, measure = c("pearson", "negdist"), diag = FALSE,
          iter = 3000, temp.init = 10, cool = 0.995, hill.climb.refine = TRUE,
          seed = NULL, verbose = TRUE)
{
  goffun <- switch(match.arg(measure), pearson = gcor, negdist = function(x,
                                                                          y) {
    -hdist(x, y)
  })
  g <- as.sociomatrix.sna(g)
  if (is.list(g))
    return(lapply(g, block.fit, model = model, measure = measure,
                  iter = iter, temp.init = temp.init, cool = cool,
                  seed = seed))
  else if (length(dim(g)) > 2)
    return(apply(g, 1, block.fit, model = model, measure = measure,
                 iter = iter, temp.init = temp.init, cool = cool,
                 seed = seed))
  n <- NROW(g)
  nb <- sqrt(length(model))
  if (floor(nb) != nb)
    stop("Specified model must be in the form of an adjacency structure.")
  model <- matrix(model, nb, nb)
  if (is.null(seed))
    memb.cur <- sample(1:nb, n, replace = TRUE)
  else memb.cur <- seed
  memb.best <- memb.cur
  gof.best <- goffun(g, model[memb.cur, memb.cur])
  gof.cur <- gof.best
  if (verbose)
    cat("Entering annealing loop...\n")
  temp <- temp.init
  for (i in 1:iter) {
    if (verbose && (i%%100 == 0))
      cat("\tIteration ", i, ", current GOF is ", gof.cur,
          ".  (Best GOF=", gof.best, ")\n", sep = "")
    memb.cand <- memb.cur
    tochange <- sample(1:n, min(1 + rgeom(1, 0.5), n))
    memb.cand[tochange] <- sample(1:nb, length(tochange),
                                  replace = TRUE)
    gof.cand <- goffun(g, model[memb.cand, memb.cand])
    if (is.na(gof.cand))
      gof.cand <- -Inf
    if (log(runif(1)) < (gof.cand - gof.cur)/temp) {
      memb.cur <- memb.cand
      gof.cur <- gof.cand
      if (gof.cand > gof.best) {
        memb.best <- memb.cand
        gof.best <- gof.cand
      }
    }
    temp <- temp * cool
  }
  if (verbose)
    cat("Annealing completed.\n")
  if (hill.climb.refine) {
    if (verbose)
      cat("Refining solution via hill-climbing procedure...\n")
    converged <- FALSE
    while (!converged) {
      if (verbose)
        cat("\tRefining; current GOF is ", gof.best,
            "\n")
      converged <- TRUE
      for (i in 1:n) {
        for (j in 1:nb) if (memb.best[i] != j) {
          ostate <- memb.best[i]
          memb.best[i] <- j
          gof.cand <- goffun(g, model[memb.best, memb.best])
          if (is.na(gof.cand))
            gof.cand <- -Inf
          if (gof.cand > gof.best) {
            gof.best <- gof.cand
            converged <- FALSE
          }
          else memb.best[i] <- ostate
        }
      }
    }
  }
  if (verbose)
    cat("Preparing and returning output.\n")
  bm <- list()
  bm$block.membership <- memb.best
  bm$order.vector <- order(memb.best)
  bm$block.content <- "types"
  bm$blocked.data <- g[bm$order.vector, bm$order.vector]
  bm$block.model <- model
  bm$rlabels <- paste("Block", 1:nb)
  bm$cluster.method <- "confirmatory"
  bm$equiv.fun <- "block.fit"
  bm$equiv.metric <- measure
  bm$block.gof <- gof.best
  class(bm) <- "blockmodel"
  bm
}
