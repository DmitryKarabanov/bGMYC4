#' Log-prior density of the GMYC parameters
#'
#' @param params Numeric vector `c(py, pc, t)`.
#' @inheritParams bgmyc.singlephy
#' @details Sum of three uniform log-densities. For the integer threshold this
#'   is a discrete uniform up to a constant (cancels in the MH ratio). Modify
#'   this function to use gamma or lognormal priors; it plugs directly into
#'   [bgmyc.gibbs()].
#' @return Log prior density (scalar).
#' @export
bgmyc.prior <-
function(params, py1, py2, pc1, pc2, t1, t2){
	loglik<-dunif(params[1], min=py1, max=py2, log=TRUE)+
			dunif(params[2], min=pc1, max=pc2, log=TRUE)+
			dunif(params[3], min=t1, max=t2, log=TRUE)
	return(loglik)
	}
