data {
	int<lower=1> N;																				// number of observations
	int<lower=1> S;																				// number of diplotype states
	int<lower=1> J;																				// number of additive haplotype effects
	int<lower=1> M;																				// number of allelic series
	int<lower=0> K;																				// number of covariate effects
	vector[N] y;																				// vector of observations
	matrix<lower=0, upper=1>[N, S] P;															// matrix of diplotype state probabilities
	matrix<lower=0, upper=1>[S, J] A;															// additive design matrix
	matrix<lower=0>[J, J] L[M];																	// array of allelic series covariance cholesky factors
	vector<upper=0>[M] ln_prior_M;																// vector of log allelic series probabilities
	real<lower=0> kappa;																		// shape for prior precision of error
	real<lower=0> lambda;																		// rate for prior precision of error
	real<lower=0> tau_mu;																		// scaled standard deviation of prior mean of additive effects
	real<lower=0> tau_delta;																	// scaled standard deviation of covariate effects
	real<lower=0> nu;																			// degrees of freedom for scaled standard deviation of additive effects
	vector<lower=0, upper=1>[N] w_inv;															// vector of weights (inverse number of replicates)
	matrix[N,K] Z;																				// design matrix for covariates
	vector<lower=1, upper=1>[J] ones;															// vector of ones for prior mean of additive effects
	vector<lower=0, upper=0>[N] zeros;															// vector of zeros for linear combination of covariate effects when K=0
}
transformed data {
	vector[N] y_std;
	matrix[N, S] log_P;
	vector[N] sqrt_w_inv;
	y_std = (y - mean(y)) / sd(y);																// standardized vector of observations
	log_P = log(P);																				// matrix of log diplotype state probabilities
	sqrt_w_inv = sqrt(w_inv);																	// vector of square root weights
}
parameters {
	vector[J] beta_h;																			// additive effects
	real<lower=0> sigma_sq_inv;																	// precision of error
	real mu;																					// mean of additive effects
	real<lower=0> phi;																			// scaled standard deviation of additive effects
	vector[K] delta;																			// covariate effects
}
transformed parameters {
	real<lower=0> sigma;
	sigma = sigma_sq_inv^(-0.5);																// standard deviation of error
}
model {
	vector[S] beta_d = A * beta_h;																// diplotype state effects
	vector[M] lpa = ln_prior_M;																	// prior for allelic series
	vector[N] Z_delta;																			// linear combination of covariate effects
	if (K>0)
		Z_delta = Z * delta;																	
	else
		Z_delta = zeros;
	target += gamma_lpdf(sigma_sq_inv | 0.5*kappa, 0.5*lambda);									// prior for precision of error
	target += normal_lpdf(mu | 0, sigma*tau_mu);												// prior for mean of additive effects
	target += normal_lpdf(delta | 0, sigma*tau_delta);											// prior for covariate effects
	target += student_t_lpdf(phi | nu, 0, 1);													// prior for scaled standard deviation of additive effects
	target += log(2);                                     										// correction for left-truncation of phi (half-t)
	for (m in 1:M) {
		lpa[m] += multi_normal_cholesky_lpdf(beta_h | mu*ones, sigma*phi*L[m]);					// calculate prior for additive effects and current allelic series
	}
	target += log_sum_exp(lpa);																	// prior for additive effects
	for (n in 1:N) {											 
		row_vector[S] lps = log_P[n];															// prior diplotype state probabilties for current observation
		for (s in 1:S) {										 
			lps[s] += normal_lpdf(y_std[n] | Z_delta[n] + beta_d[s], sqrt_w_inv[n]*sigma);		// calculate likelihood and prior for current observation and diplotype state
		}
		target += log_sum_exp(lps);																// likelihood for current observation 
	}
}
generated quantities {
	vector[M] ln_post_M;
	ln_post_M = ln_prior_M;
	for (m in 1:M) {
		ln_post_M[m] += multi_normal_cholesky_lpdf(beta_h | mu*ones, sigma*phi*L[m]);			// prior for additive effects and current allelic series
	}
	ln_post_M = ln_post_M - log_sum_exp(ln_post_M);												// posterior allelic series probabilities
}
