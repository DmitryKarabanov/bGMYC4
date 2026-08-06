bgmyc.singlephy <-
function(phylo, mcmc, burnin, thinning, py1=0, py2=2, pc1=0, pc2=2, t1=2, t2=51, scale=c(20, 10, 5.00), start=c(1.0, 0.5, 50.0), sampler=bgmyc.gibbs, likelihood=bgmyc.lik, prior=bgmyc.prior){
	
  message("You are running bGMYC on a single phylogenetic tree.\n")
  message("This tree contains ", length(phylo$tip.label), " tips.\n")
  message("The Yule process rate change parameter has a uniform prior ranging from ", py1, " to ", py2, ".\n")
  message("The coalescent process rate change parameter has a uniform prior ranging from ", pc1, " to ", pc2, ".\n")
  message("The threshold parameter, which is equal to the number of species, has a uniform prior ranging from ", t1, " to ", t2, ". The upper bound of this prior should not be more than the number of tips in your trees.\n")
  message("The MCMC will start with the Yule parameter set to ", start[1], ".\n")
  message("The MCMC will start with the coalescent parameter set to ", start[2], ".\n")
  message("The MCMC will start with the threshold parameter set to ", start[3], ". If this number is greater than the number of tips in your tree, an error will result.\n")
  message("Given your settings for mcmc, burnin and thinning, your analysis will result in ", (mcmc-burnin)/thinning, " samples being retained.\n")
  
	bgmyc.dataprep(phylo)->data
	sampler(data, m=mcmc, burnin, thinning, py1, py2, pc1, pc2, t1, t2, scale=scale, start=start, likelihood=likelihood, prior=prior)->output
	return(output)
	
	}
