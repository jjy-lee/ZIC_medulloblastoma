#Zhen You
##### Beta Binomial Bayesian Model for allele imbalance ratio estimation;
##### Project: ZIC1 is a context dependent cancer driver in the rhombic lip
##### This is used for estimating true allele imbalance ratio of heterozygous SNPs from ChIP-Seq data or
##### RNA-Seq data, with the matched heterozygous SNPs in WGS

##### Usage: 1. source this R script to load the functions to the current R working environment;
##### Usage: 2. use main_fun to calculate allele imbalance and statistic significance;
##### output: the main function will return list containing p value and allele imbalance




### the main function
# Alt_count: alternative allele counts for ChIP/RNA-Seq;
# total_count: counts from two alleles in ChIP/RNA-Seq;
# ratios: heterozygous SNPs WGS ratios
# variance: set to 0.05
# precision: set to 1000

main_fun <- function(Alt_count, total_count, ratios, variance, precision) {

	rd = seq(0.00001,0.99999,length.out=19999);   # simmulated allele ratios as independent variables

	# the allele imbalace ratio is estimated as the ratio corresponding to the summit of the posterior pdf
	llh_true = rd[which.max(exp(post_pdf(rd, Alt_count, total_count, ratios, variance, precision)))]
	# simulate the probability that allele imbalance ratios are smaller than 0.5
	pdfsum = sum(exp(post_pdf(rd, Alt_count, total_count, ratios, variance, precision)))
	llh_p <- sum(exp(post_pdf(rd[1:10000], Alt_count, total_count, ratios, variance, precision)))/pdfsum

	llh_p = ifelse(llh_true < 0.5, 1 - llh_p, llh_p)
	results = list("Pval" = llh_p, "Ratio" = llh_true)
	return(results)
}




##  Bayesian Model some basic function

logtranfac <- function(count) {
	# write factorial as gamma function a! = Γ(a + 1)
	lgamma(count + 1)   ## natural logarithm of gamma function.
	# this allows the binomial coefficient to be generalized to noninteger arguments
}


## beta function
logtranbeta <- function(alpha,beta) {
	# use beta as gamma function
	lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta)

}

# natural logarithm of binomial coefficient
logtranbincoef <- function(n,k) {
	# generalized by factorial mentioned before
	logtranfac(n) - logtranfac(k) - logtranfac(n - k)

}


# join probability density function of beta binomial
betabinpdf <- function(k, n, alpha, beta){
	# This is written using beta properties
	logtranbincoef(n, k) + logtranbeta(alpha + k, n - k + beta) - logtranbeta(alpha, beta)

}




# likelihood function
likelihoodfunc <- function(Alt_count, total_count, wgs_allele, allele_ratio, precision){

	mu <- wgs_allele		# WGS het. SNP counterpart ratio
	# bridging WGS raitos with ChIP/RNA-Seq ratios
	xi <- (allele_ratio * mu)/(1-allele_ratio-mu + 2*allele_ratio*mu)		## referenced from de Santiago et al. 2017 DOI: 10.1186/s13059-017-1165-7

	param_a <- xi * precision
	param_b <- (1 - xi) * precision
	return(betabinpdf(Alt_count, total_count, param_a, param_b))
}



# Posterior distribution by Bayesian inference
post_pdf <- function(allele_ratio, Alt_count, total_count, ratios, variance, precision) {
	precision <- precision
	pi <- 0.5
	wgs_allele <- ratios
	# likelihood function
	llh = log(pi) + likelihoodfunc(Alt_count, total_count, wgs_allele, allele_ratio, precision)

	# prior model
	mu <- ratios  # reference mapping bias from WGS
	var <- variance			# set to 0.05 as de Santiago et al. suggested
	# Method of moments
	alpha <- ((1-mu)/var- 1/mu) * mu^2
	beta <- alpha*(1/mu - 1)
	prior <- log(dbeta(allele_ratio,alpha,beta))  # prior

	# posterior pdf
	total_pos <- llh + prior + log(wgs_allele/(allele_ratio*wgs_allele+(1-allele_ratio)*(1-wgs_allele)) - allele_ratio*wgs_allele*(2*wgs_allele-1)/(allele_ratio*wgs_allele+(1-allele_ratio)*(1-wgs_allele))^2) # Derivation
	return(total_pos)

}




################################################
################################################  End
################################################
