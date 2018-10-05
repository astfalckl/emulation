
# library(plotly)
# library(MASS)
# library(ggplot2)
# library(optimx) 
# library(coda)
# library(mvtnorm)
# library(reshape)
# library(LaplacesDemon)
# library(akima)
# library(optimx)
# library(wesanderson)

set.seed(12345)

#---------------------------------------------------------------------------------------------------

pow_exp <- function(X1, X2, sigma, l1, gamma, noise){
	#X1 and X2 are vectors of univariate inputs
	Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
	dist <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))

	for (i in 1:nrow(Sigma)){
		for (j in 1:ncol(Sigma)){
			theta <- min(abs(X2[j] - X1[i]), 2*pi - abs(X2[j] - X1[i]))
			dist[i,j] <- theta
			Sigma[i,j] <- sigma^2 * exp(-0.5*(theta/l1)^gamma)
		}
	}

	Sigma <- Sigma + diag(noise, nrow = length(X1), ncol = length(X2))
	return(Sigma)
} 

circ_wend2 <- function(X1, X2, sigma, tau, noise){
	#X1 and X2 are vectors of univariate inputs
	Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
	dist <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))

	if (tau < 4) {return("Tau must be more than 4")}

	for (i in 1:nrow(Sigma)){
		for (j in 1:ncol(Sigma)){
			theta <- min(abs(X2[j] - X1[i]), 2*pi - abs(X2[j] - X1[i]))
			dist[i,j] <- theta
			Sigma[i,j] <- sigma^2 * (1 + tau * theta/pi) * max((1 - theta/pi)^tau, 0)
		}
	}

	Sigma <- Sigma + diag(noise, nrow = length(X1), ncol = length(X2))
	return(Sigma)
} 

circ_wend4 <- function(X1, X2, sigma, tau, noise){
	#X1 and X2 are vectors of univariate inputs
	Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
	dist <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))

	if (tau < 6) {return("Tau must be more than 6")}

	for (i in 1:nrow(Sigma)){
		for (j in 1:ncol(Sigma)){
			theta <- min(abs(X2[j] - X1[i]), 2*pi - abs(X2[j] - X1[i]))
			dist[i,j] <- theta
			Sigma[i,j] <- sigma^2 * (1 + tau * theta/pi + (tau^2 - 1)/3 * (theta^2)/(pi^2)) * max((1 - theta/pi)^tau, 0)
		}
	}

	Sigma <- Sigma + diag(noise, nrow = length(X1), ncol = length(X2))
	return(Sigma)
}

circ_gen_cauchy <- function(X1, X2, sigma, alpha, tau, noise){
	#X1 and X2 are vectors of univariate inputs
	Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
	dist <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))

	for (i in 1:nrow(Sigma)){
		for (j in 1:ncol(Sigma)){
			theta <- min(abs(X2[j] - X1[i]), 2*pi - abs(X2[j] - X1[i]))
			dist[i,j] <- theta
			Sigma[i,j] <- sigma^2 * (1 + (theta/pi)^alpha)^(-tau/alpha)
		}
	}

	Sigma <- Sigma + diag(noise, nrow = length(X1), ncol = length(X2))
	return(Sigma)
}

circ_dagum <- function(X1, X2, sigma, alpha, tau, noise){
	#X1 and X2 are vectors of univariate inputs
	Sigma <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))
	dist <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))

	for (i in 1:nrow(Sigma)){
		for (j in 1:ncol(Sigma)){
			theta <- min(abs(X2[j] - X1[i]), 2*pi - abs(X2[j] - X1[i]))
			dist[i,j] <- theta
			Sigma[i,j] <- sigma^2 * (1 - (((theta/pi)^tau)/((1 + theta/pi)^tau))^(alpha/tau))
		}
	}

	Sigma <- Sigma + diag(noise, nrow = length(X1), ncol = length(X2))
	return(Sigma)
}

prod_exppow_wend4 <- function(X1, X2, sigma = 1, tau = 6, l = 1, gamma = 2, noise = 1e-06) {
	# We are now in multivariate land, where X1 and X2 are multivariate (but separable?)
	# First column of X1 is linear magnitude and second is circular. THIS NEEDS TO CHANGE TO BE MORE GENERAL.

	Sigma <- matrix(rep(0, length(X1[ ,1])*length(X2[ ,1])), nrow = length(X1[ ,1]))
	circ_dist <- matrix(rep(0, length(X1[ ,1])*length(X2[ ,1])), nrow = length(X1[ ,1]))

	for (i in 1:nrow(Sigma)){
		for (j in 1:ncol(Sigma)){
			theta <- min(abs(X2[j,2] - X1[i,2]), 2*pi - abs(X2[j,2] - X1[i,2]))
			# print(theta)
			circ_dist[i,j] <- theta

			# Covariance function is a product kernel of the powered exponential, and the Wendland-4 function.
			# print(sigma)
			# print(tau)
			# print(l)
			# print(gamma)
			Sigma[i,j] <- sigma^2 * (exp(-0.5* (abs(X1[i,1] - X2[j,1])^gamma) / l )) * (1 + tau * theta/pi + (tau^2 - 1)/3 * (theta^2)/(pi^2)) * max((1 - theta/pi)^tau, 0)
		}
	}

	Sigma <- Sigma + diag(noise, nrow = length(X1[ ,1]), ncol = length(X2[ ,1]))
	return(Sigma)
}

k_current_wind <- function(d, sigma = p1, tauwind = p2, taucurr = p3, lwind = p4, lcurr = p5, noise = noise) {

	k1 <- (exp(-0.5* (abs(d[[1]])^2) / lwind ))
	k2 <- (1 + tauwind * d[[2]]/pi + (tauwind^2 - 1)/3 * (d[[2]]^2)/(pi^2)) * pmax((1 - d[[2]]/pi)^tauwind, 0)
	k3 <- (exp(-0.5* (abs(d[[3]])^2) / lcurr ))
	k4 <- (1 + taucurr * d[[4]]/pi + (taucurr^2 - 1)/3 * (d[[4]]^2)/(pi^2)) * pmax((1 - d[[4]]/pi)^taucurr, 0)

	Sigma <- sigma^2 * k1 * k2 * k3 * k4

	Sigma <- Sigma + diag(noise, nrow = length(Sigma[ ,1]), ncol = length(Sigma[1, ]))
	return(Sigma)

	# for (i in 1:nrow(Sigma)){
	# 	for (j in 1:ncol(Sigma)){
	# 		theta1 <- min(abs(X2[j,2] - X1[i,2]), 2*pi - abs(X2[j,2] - X1[i,2]))
	# 		theta2 <- min(abs(X2[j,4] - X1[i,4]), 2*pi - abs(X2[j,4] - X1[i,4]))

	# 		k1 <- (exp(-0.5* (abs(X1[i,1] - X2[j,1])^2) / l1 ))
	# 		k2 <- (1 + tau1 * theta1/pi + (tau1^2 - 1)/3 * (theta1^2)/(pi^2)) * max((1 - theta1/pi)^tau1, 0)
	# 		k3 <- (exp(-0.5* (abs(X1[i,3] - X2[j,3])^2) / l2 ))
	# 		k4 <- (1 + tau2 * theta2/pi + (tau2^2 - 1)/3 * (theta2^2)/(pi^2)) * max((1 - theta2/pi)^tau2, 0)

	# 		Sigma[i,j] <- sigma^2 * k1 * k2 * k3 * k4
	# 	}
	# }

	# Sigma <- Sigma + diag(noise, nrow = length(X1[ ,1]), ncol = length(X2[ ,1]))
	# return(Sigma)

}


k_CWSW <- function(d, sigma = 1, tauwind = 6, taucurr = 6, tauswell = 6, lwind = 1, lcurr = 1, lhs = 1, ltp = 1, noise = 1e-06) {

	#d is a list of distances 

	k1 <- (exp(-0.5* (abs(d[[1]])^2) / lhs ))
	k2 <- (exp(-0.5* (abs(d[[2]])^2) / ltp ))
	k3 <- (1 + tauswell * d[[3]]/pi + (tauswell^2 - 1)/3 * (d[[3]]^2)/(pi^2)) * pmax((1 - d[[3]]/pi)^tauswell, 0)
	k4 <- (exp(-0.5* (abs(d[[4]])^2) / lwind ))
	k5 <- (1 + tauwind * d[[5]]/pi + (tauwind^2 - 1)/3 * (d[[5]]^2)/(pi^2)) * pmax((1 - d[[5]]/pi)^tauwind, 0)
	k6 <- (exp(-0.5* (abs(d[[6]])^2) / lcurr ))
	k7 <- (1 + taucurr * d[[7]]/pi + (taucurr^2 - 1)/3 * (d[[7]]^2)/(pi^2)) * pmax((1 - d[[7]]/pi)^taucurr, 0)

	Sigma <- sigma^2 * k1 * k2 * k3 * k4 * k5 * k6 * k7

	Sigma <- Sigma + diag(noise, nrow = length(Sigma[ ,1]), ncol = length(Sigma[1, ]))
	return(Sigma)

}

k_CWSW_matern <- function(d, sigma = 1, tauwind = 6, taucurr = 6, tauswell = 6, pwind = 1, pcurr = 1, phs = 1, ptp = 1, noise = 1e-06) {

	#d is a list of distances 

	k1 <- (1 + sqrt(3)*abs(d[[1]])/phs) * (exp(-sqrt(3) * abs(d[[1]]) / phs ))
	k2 <- (1 + sqrt(3)*abs(d[[2]])/ptp) * (exp(-sqrt(3) * abs(d[[2]]) / ptp ))
	k3 <- (1 + tauswell * d[[3]]/pi + (tauswell^2 - 1)/3 * (d[[3]]^2)/(pi^2)) * pmax((1 - d[[3]]/pi)^tauswell, 0)
	k4 <- (1 + sqrt(3)*abs(d[[4]])/pwind) * (exp(-sqrt(3) * abs(d[[4]]) / pwind ))
	k5 <- (1 + tauwind * d[[5]]/pi + (tauwind^2 - 1)/3 * (d[[5]]^2)/(pi^2)) * pmax((1 - d[[5]]/pi)^tauwind, 0)
	k6 <- (1 + sqrt(3)*abs(d[[6]])/pcurr) * (exp(-sqrt(3) * abs(d[[6]]) / pcurr ))
	k7 <- (1 + taucurr * d[[7]]/pi + (taucurr^2 - 1)/3 * (d[[7]]^2)/(pi^2)) * pmax((1 - d[[7]]/pi)^taucurr, 0)

	Sigma <- sigma^2 * k1 * k2 * k3 * k4 * k5 * k6 * k7

	Sigma <- Sigma + diag(noise, nrow = length(Sigma[ ,1]), ncol = length(Sigma[1, ]))
	return(Sigma)

}


k_All <- function(d, hp_MAP) {

	#d is a list of distances 

	sigma <- hp_MAP[1]
	tauwind <- hp_MAP[2]
	taucurr <- hp_MAP[3]
	tauswell <- hp_MAP[4]
	tausea <- hp_MAP[5]
	lwind <- hp_MAP[6]
	lcurr <- hp_MAP[7]
	lhsswell <- hp_MAP[8]
	ltpswell <- hp_MAP[9]
	lhssea <- hp_MAP[10]
	ltpsea <- hp_MAP[11]
	noise <- hp_MAP[12]

	k1 <- (exp(-0.5* (abs(d[[1]])^2) / lhsswell ))
	k2 <- (exp(-0.5* (abs(d[[2]])^2) / ltpswell ))
	k3 <- (1 + tauswell * d[[3]]/pi + (tauswell^2 - 1)/3 * (d[[3]]^2)/(pi^2)) * pmax((1 - d[[3]]/pi)^tauswell, 0)
	k4 <- (exp(-0.5* (abs(d[[4]])^2) / lhssea ))
	k5 <- (exp(-0.5* (abs(d[[5]])^2) / ltpsea ))
	k6 <- (1 + tausea * d[[6]]/pi + (tausea^2 - 1)/3 * (d[[6]]^2)/(pi^2)) * pmax((1 - d[[6]]/pi)^tausea, 0)
	k7 <- (exp(-0.5* (abs(d[[7]])^2) / lwind ))
	k8 <- (1 + tauwind * d[[8]]/pi + (tauwind^2 - 1)/3 * (d[[8]]^2)/(pi^2)) * pmax((1 - d[[8]]/pi)^tauwind, 0)
	k9 <- (exp(-0.5* (abs(d[[9]])^2) / lcurr ))
	k10 <- (1 + taucurr * d[[10]]/pi + (taucurr^2 - 1)/3 * (d[[10]]^2)/(pi^2)) * pmax((1 - d[[10]]/pi)^taucurr, 0)

	Sigma <- sigma^2 * k1 * k2 * k3 * k4 * k5 * k6 * k7 * k8 * k9 * k10

	Sigma <- Sigma + diag(noise, nrow = length(Sigma[ ,1]), ncol = length(Sigma[1, ]))
	return(Sigma)

}


k_All_matern <- function(d, sigma = 1, tauwind = 6, taucurr = 6, tauswell = 6, tausea = 6, pwind = 1, pcurr = 1, phsswell = 1, ptpswell = 1, phssea = 1, ptpsea = 1, noise = 1e-06) {

	#d is a list of distances 
	k1 <- (1 + sqrt(3)*abs(d[[1]])/phsswell) * (exp(-sqrt(3) * abs(d[[1]]) / phsswell ))
	k2 <- (1 + sqrt(3)*abs(d[[2]])/ptpswell) * (exp(-sqrt(3) * abs(d[[2]]) / ptpswell ))
	k3 <- (1 + tauswell * d[[3]]/pi + (tauswell^2 - 1)/3 * (d[[3]]^2)/(pi^2)) * pmax((1 - d[[3]]/pi)^tauswell, 0)
	k4 <- (1 + sqrt(3)*abs(d[[4]])/phssea) * (exp(-sqrt(3) * abs(d[[4]]) / phssea ))
	k5 <- (1 + sqrt(3)*abs(d[[5]])/ptpsea) * (exp(-sqrt(3) * abs(d[[5]]) / ptpsea ))
	k6 <- (1 + tausea * d[[6]]/pi + (tausea^2 - 1)/3 * (d[[6]]^2)/(pi^2)) * pmax((1 - d[[6]]/pi)^tausea, 0)
	k7 <- (1 + sqrt(3)*abs(d[[7]])/pwind) * (exp(-sqrt(3) * abs(d[[7]]) / pwind ))
	k8 <- (1 + tauwind * d[[8]]/pi + (tauwind^2 - 1)/3 * (d[[8]]^2)/(pi^2)) * pmax((1 - d[[8]]/pi)^tauwind, 0)
	k9 <- (1 + sqrt(3)*abs(d[[9]])/pcurr) * (exp(-sqrt(3) * abs(d[[9]]) / pcurr ))
	k10 <- (1 + taucurr * d[[10]]/pi + (taucurr^2 - 1)/3 * (d[[10]]^2)/(pi^2)) * pmax((1 - d[[10]]/pi)^taucurr, 0)

	Sigma <- sigma^2 * k1 * k2 * k3 * k4 * k5 * k6 * k7 * k8 * k9 * k10

	Sigma <- Sigma + diag(noise, nrow = length(Sigma[ ,1]), ncol = length(Sigma[1, ]))
	return(Sigma)

}

dist_circ <- function(X1, X2) {

	dist <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))

	for (i in 1:nrow(dist)){
		for (j in 1:ncol(dist)){
			theta <- min(abs(X2[j] - X1[i]), 2*pi - abs(X2[j] - X1[i]))
			dist[i,j] <- theta
		}
	}

	return(dist)

}

dist_lin <- function(X1, X2) {

	dist <- matrix(rep(0, length(X1)*length(X2)), nrow = length(X1))

	for (i in 1:nrow(dist)){
		for (j in 1:ncol(dist)){
			r <- abs(X2[j] - X1[i])
			dist[i,j] <- r
		}
	}

	return(dist)

}

k_cw <- function(d1, d2, d3, d4, sigma = p1, tau1 = p2, tau2 = p3, l1 = p4, l2 = p5, noise = noise) {

	k1 <- (exp(-0.5* (d1^2) / l1 ))
	k2 <- (1 + tau1 * d2/pi + (tau1^2 - 1)/3 * (d2^2)/(pi^2)) * max((1 - d2/pi)^tau1, 0)
	k3 <- (exp(-0.5* (d3^2) / l2 ))
	k4 <- (1 + tau2 * d4/pi + (tau2^2 - 1)/3 * (d4^2)/(pi^2)) * max((1 - d4/pi)^tau1, 0)

	Sigma <- sigma^2 * k1 * k2 * k3 * k4

	return(Sigma)

}

polar_GP_sample <- function(method, points = 100, sigma = 1, l = 1, gamma = 1, tau = 4, alpha = 1, noise = 0.0000001, plot = FALSE) {
	X <- seq(0,2*pi, 2*pi/(points-1))

	if (method == "Wendland-2") {
		Cov_circ <- circ_wend2(X, X, sigma, tau, noise)
	} else if (method == "Wendland-4") {
		Cov_circ <- circ_wend4(X, X, sigma, tau, noise)
	} else if (method == "Exp-1") {
		Cov_circ <- pow_exp(X, X, sigma, l, gamma, noise)
	} else if (method == "GenCauchy") {
		Cov_circ <- circ_gen_cauchy(X, X, sigma, alpha, tau, noise)
	} else if (method == "Dagum") {
		Cov_circ <- circ_dagum(X, X, sigma, alpha, tau, noise)
	} else {
		return("No valid covariance function specified")
	}

	sample_circ <- rmvnorm(1, mean = rep(0,length(Cov_circ[ ,1])), Cov_circ, method = "chol")
	sample_circ <- as.vector(sample_circ)
	df_circ <- data.frame(X = X, Y = sample_circ)
	# df_circ <- data.frame(X = c(X,X[1]), Y = c(sample_circ,sample_circ[1]))

	if (plot) {
	plot_sample <- ggplot(data = df_circ) + 
		geom_line(aes(x = X, y = Y)) +
		scale_x_continuous(breaks = seq(0,2*pi,2*pi/8)) +
		coord_polar() +
		# xlim(0,2*pi) +
		ylim(-3,3) + 
		theme_bw() +
		theme(panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5),
				axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
				axis.title.x = element_blank(), axis.text.x = element_blank(),
				axis.title.y = element_blank(),panel.border = element_blank())

	print(plot_sample)
	}

	return(df_circ)
}

polar_GP_samples <- function(method, points = 100, lines = 3, sigma = 1, l = 1, gamma = 1, tau = 4, alpha = 1, noise = 0.0000001, plot = FALSE){
	samples_df <- NULL
	samples_melt <- NULL

	samples_df <- data.frame(X = seq(0,2*pi, 2*pi/(points-1)))

	for (i in 1:lines){
		sample_temp <- polar_GP_sample(method, points, sigma, l, gamma, tau, alpha, noise, plot = FALSE)
		samples_df[i+1] <- sample_temp[ ,2]
		names(samples_df)[i+1] <- paste("Y",i,sep="")
	}

	samples_melt <- melt(samples_df[ ,-1])
	samples_melt[3] <- rep(seq(0,2*pi, 2*pi/(points-1)),lines)
	names(samples_melt)[3] <- "X"

	circles <- data.frame(X = seq(0,2*pi, 2*pi/(points-1)), Zero = rep(0,points), Minus = rep(-3*sigma,points), Plus = rep(3*sigma,points))

	pal <- wes_palette("Zissou", lines, type = "continuous")

	plot_sample <- ggplot(data = samples_melt) + 
	geom_line(aes(x = X, y = value, colour = variable), alpha = 0.5) +
	scale_color_manual(values = pal) +
	geom_line(data = circles, aes(x = X, y = Zero), color = "#000000", size = 1.5) +
	geom_line(data = circles, aes(x = X, y = Minus), color = "#909090", size = 1.5, linetype = "dashed") +
	geom_line(data = circles, aes(x = X, y = Plus), color = "#909090", size = 1.5, linetype = "dashed") +
	scale_x_continuous(breaks = seq(0,2*pi,2*pi/8)) +
	coord_polar() +
	# xlim(0,2*pi) +
	ylim(-5*sigma,4*sigma) + 
	theme_bw() +
	theme(panel.grid.minor=element_blank(), plot.title = element_text(hjust = 0.5),
			axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
			axis.title.x = element_blank(), axis.text.x = element_blank(),
			axis.title.y = element_blank(),panel.border = element_blank(),legend.position="none")

	return(plot_sample)

}

log_likelihood_chol <- function(X, Y, method, p1 = 0, p2 = 0, p3 = 0, p4 = 0, p5 = 0, p6 = 0, p7 = 0, p8 = 0, p9 = 0, p10 = 0, p11 = 0, noise = 0){
	# n <- length(Y)

	if (method == "Wendland-2") {
		Kxx <- circ_wend2(X, X, sigma = p1, tau = p2, noise = noise)
	} else if (method == "Wendland-4") {
		Kxx <- circ_wend4(X, X, sigma = p1, tau = p2, noise = noise)
	} else if (method == "Exp-powered") {
		Kxx <- pow_exp(X, X, sigma = p1, l = p2, gamma = p3, noise = noise)
	} else if (method == "GenCauchy") {
		Kxx <- circ_gen_cauchy(X, X, sigma = p1, alpha = p2, tau = p3, noise = noise)
	} else if (method == "Dagum") {
		Kxx <- circ_dagum(X, X, sigma = p1, alpha = p2, tau = p3, noise = noise)
	} else if (method == "Exp-pow/Wendland-4") {
		Kxx <- prod_exppow_wend4(X, X, sigma = p1, tau = p2, l = p3, gamma = p4, noise = noise)
	} else if (method == "W4-W4-SE-SE") {
		Kxx <- k_current_wind(d = X, sigma = p1, tauwind = p2, taucurr = p3, lwind = p4, lcurr = p5, noise = noise)
	} else if (method == "CWSW") {
		Kxx <- k_CWSW(d = X, sigma = p1, tauwind = p2, taucurr = p3, tauswell = p4, lwind = p5, lcurr = p6, lhs = p7, ltp = p8, noise = noise)
	} else if (method == "CWSW_matern") {
		Kxx <- k_CWSW_matern(d = X, sigma = p1, tauwind = p2, taucurr = p3, tauswell = p4, pwind = p5, pcurr = p6, phs = p7, ptp = p8, noise = noise)
	} else if (method == "All") {
		Kxx <- k_All(d = X, sigma = p1, tauwind = p2, taucurr = p3, tauswell = p4, tausea = p5, lwind = p6, lcurr = p7, lhsswell = p8, ltpswell = p9, lhssea = p10, ltpsea = p11, noise = noise)
	} else if (method == "All_matern") {
		Kxx <- k_All_matern(d = X, sigma = p1, tauwind = p2, taucurr = p3, tauswell = p4, tausea = p5, pwind = p6, pcurr = p7, phsswell = p8, ptpswell = p9, phssea = p10, ptpsea = p11, noise = noise)
	} else {
		return("No valid covariance function specified")
	}

	U <- chol(Kxx)
	# U_inv <- solve(U)
	# Kxx_logdet <- log(det(U)^2)
	# Y_matrix <- as.matrix(Y)
	# # log_likelihood <- -n/2*log(2*pi)  - (0.5 * t(Y_matrix) %*% U_inv %*% t(U_inv) %*% Y_matrix) - Kxx_logdet/2

	# print(-n/2*log(2*pi)  - (0.5 * t(Y_matrix) %*% U_inv %*% t(U_inv) %*% Y_matrix) - Kxx_logdet/2)

	log_likelihood <- - sum(log(diag(U))) - 0.5 * sum(backsolve(U, Y, transpose = TRUE) ^ 2)

	# print(log_likelihood)

	return(log_likelihood)
}

chol_solve <- function(R, b) {
	backsolve(R, backsolve(R,b, transpose = TRUE))
}

SE_Kxx_2dim <- function(X1, X2, sigma, l1, l2, noise){
	Sigma <- matrix(rep(0, length(X1[ ,1])*length(X2[ ,1])), nrow = length(X1[ ,1]))
	for (i in 1:nrow(Sigma)){
		for (j in 1:ncol(Sigma)){
			#theta <- min(abs(X2[j,2] - X1[i,2]), 2*pi - abs(X2[j,2] - X1[i,2]))
			#Sigma[i,j] <- sigma^2 * exp(-0.5*((abs(X1[i,1] - X2[j,1])^2)/l1)) * exp(-0.5*(theta)/l2)
			Sigma[i,j] <- sigma^2 * exp(-0.5*((abs(X1[i,1] - X2[j,1])^2)/l1 + (abs(X1[i,2] - X2[j,2])^2)/l2))
		}
	}
	Sigma <- Sigma + diag(noise, nrow = length(X1[ ,1]), ncol = length(X2[ ,1]))
	return(Sigma)
} 

plot_gp_realisations <- function(params, data_circ) {
	I <- length(params$y_predict[,1])

	plot(1, type="n", xlab="x", ylab="y", main="Posterior Realisations", xlim = c(0,2*pi), ylim = c(0,10))
	for (i in 1:I){
	       lines(data_circ$x_predict, params$y_predict[i,], col=c("#8F272705"))}
	points(data_circ$x_predict, data_circ$y_predict, col="white", pch=16, cex=0.6)
	points(data_circ$x_predict, data_circ$y_predict, col=c_mid_teal, pch=16, cex=0.4)
	points(data_circ$x, data_circ$y, col="white", pch=16, cex=1.2)
	points(data_circ$x, data_circ$y, col="black", pch=16, cex=0.8)

	p <- recordPlot()

	return(p)
}

coverage <- function(samples, y_actual) {
	coverage <- data.frame(actual = seq(0,0.99,0.01))
	for (k in 1:100) {
		upper <-  apply(samples, 2, quantile, 0.50 + k/200)
		lower <-  apply(samples, 2, quantile, 0.50 - k/200)

		in_quantile <- 0

		for (m in 1:length(samples[1, ])) {
		  if (y_actual[m] < upper[m] & y_actual[m] > lower[m]) {
		    in_quantile <- in_quantile + 1
		  }  
		}

		coverage[k,2] <- in_quantile/length(samples[1, ])
		cat(sprintf('\rIteration %d/%d', k, 100))
	}
	return(coverage)
}

coverage_plot <- function(samples, y_actual) {

	cov <- coverage(samples, y_actual)

	p <- ggplot(data = cov) +
			geom_line(aes_string(x = colnames(cov)[1], y = colnames(cov)[1])) +
			geom_point(aes_string(x = colnames(cov)[1], y = colnames(cov)[2])) +
			theme_bw() +
			xlab("Nominal Coverage") +
			ylab("Empirical Coverage") +
			ggtitle("Coverage Plot") +
			coord_fixed(ratio = 1) +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	return(p)
}

IPE_plot <- function(samples, y_actual) {
	IPE <- data.frame(index = seq(1,length(y_actual),1))
	E <- apply(samples, 2, mean)
	V <- apply(samples, 2, var)
	IPE[ ,2] <- (y_actual - E) / sqrt(V)

	p <- ggplot(data = IPE) +
		geom_point(aes_string(x = colnames(IPE)[1], y = colnames(IPE)[2]), shape = 4) +
		geom_hline(yintercept = -3, linetype = 2) + 
		geom_hline(yintercept = 0) + 
		geom_hline(yintercept = 3, linetype = 2) +
		theme_bw() +
		ylim(-6,6) +
		xlab("Point Index") +
		ylab("Individual Prediction Error") +
		ggtitle("Individual Prediction Errors") +
		# coord_fixed(ratio = 500/12) +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)

	return(p)
}

IPE_qqplot <- function(samples, y_actual) {
	IPE <- data.frame(index = seq(1,length(y_actual),1))
	E <- apply(samples, 2, mean)
	V <- apply(samples, 2, var)
	IPE[ ,2] <- (y_actual - E) / sqrt(V)
	IPE_sort <- IPE[order(IPE[,2]), ]

	p <- ggplot(IPE_sort, aes(sample = V2)) + 
		stat_qq(distribution = qt, dparams = length(y_actual), alpha = 0.15) +
		xlim(-4,4) +
		ylim(-4,4) +
		geom_abline(intercept = 0, slope = 1) +
		theme_bw() +
		xlab("Nominal Quantiles") +
		ylab("Empirical Quantiles") +
		ggtitle("Q-Q Plot of Errors") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)

	return(p)
}

plot_diagnostics <- function(y_draw, y_predict, data_resid, sigma) {

	resid <- data.frame(index = seq(1,length(data_resid[,1]),1))
	resid[ ,2] <- (data_resid$y_actual - data_resid$y_resid)/sigma

	p1 <- var_scatter(y_draw)
	p2 <- coverage_plot(y_draw, y_predict)
	p3 <- IPE_plot(y_draw, y_predict)
	p4 <- IPE_qqplot(y_draw, y_predict)

	grobs <- list(p1, p2, p3, p4)
	lay <- rbind(c(1,2),
				 c(3,4))
	grid.arrange(grobs = grobs, layout_matrix = lay)
}


chol_solve <- function(R, b) {
	backsolve(R, backsolve(R,b, transpose = TRUE))
}

residual_plot <- function(data, sigma){
	resid <- data.frame(index = seq(1,length(data[,1]),1))
	resid[ ,2] <- (data$y_actual - data$y_resid)/sigma
	p1 <- ggplot(data, aes(x_lin, x_circ, colour = y_actual)) + geom_point() + scale_colour_gradient2(low = "red", high = "blue", mid = "white", midpoint = 1.6)
	p2 <- ggplot(data, aes(x_lin, x_circ, colour = y_resid)) + geom_point() + scale_colour_gradient2(low = "red", high = "blue", mid = "white", midpoint = 1.6)
	p3 <- residuals_scatter(resid)
	grobs <- list(p1, p2, p3)
	lay <- rbind(c(1,2,3))
	grid.arrange(grobs = grobs, layout_matrix = lay)
}

residuals_scatter <- function(resid) {
	p <- ggplot(resid, aes(resid[ ,1], resid[ ,2])) + 
	geom_point(shape = 4) +
	geom_hline(yintercept = -2, linetype = 2) + 
	geom_hline(yintercept = 0) + 
	geom_hline(yintercept = 2, linetype = 2) +
	xlab("Point Index") +
	ylab("Residual") +
	ggtitle("Residual Errors") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	return(p)
}

cov_wendland <- function(dist, tau) {
	# i <- length(dist[ ,1])
	# j <- length(dist[1, ])
	# K <- (1 + tau/pi*dist + (2*tau^2 - 3)/(5 * pi^2) * dist^2 + (tau^3 - 4*tau)/(15 * pi^3) * dist^3) * (1 - dist/pi)^tau
	K <- (1 + tau/pi*dist + (tau^2 - 1)/(3 * pi^2) * dist^2) * pmax((1 - dist/pi)^tau, 0)
	return(K)
}

cov_exp <- function(dist, rho) {
	K <- exp(-0.5*dist^2/rho^2)
	return(K)
}

cov_matern32 <- function(dist, rho) {
	K <- (1 + sqrt(3)/(rho)*dist) * exp(-sqrt(3)/(rho)*dist)
	return(K)
}

cov_prods <- function(d, param, types) {
	k <- matrix(1, nrow = nrow(d[[1]]), ncol = ncol(d[[1]]))
	for (i in 1:length(types)) {
		if(types[i] == "circular") {
			k <- k * cov_wendland(d[[i]], param[i])
		} else {
			# k <- k * cov_exp(d[[i]], param[i])
			k <- k * cov_matern32(d[[i]], param[i])
		}
	}
	return(k)
}

var_scatter <- function(y_draw, var_accept) {
	df <- data.frame(index = seq(1,length(y_draw[1,]),1))
	df[ ,2] <- apply(y_draw, 2, var)
	p <- ggplot(df, aes(df[ ,1], df[ ,2])) + 
	geom_point(shape = 4, alpha = 0.7) +
	geom_hline(yintercept = var_accept , linetype = 2) + 
	# geom_hline(yintercept = -2, linetype = 2) + 
	# geom_hline(yintercept = 0) + 
	# geom_hline(yintercept = 2, linetype = 2) +
	xlab("Point Index") +
	ylab("Variance") +
	ggtitle("Variances") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)
	return(p)
}

var_density <- function(y_draw, var_accept) {
	V <- data.frame(var = apply(y_draw, 2, var))
	p <- ggplot(V) + 
		geom_density(aes(x = var)) +
		geom_vline(xintercept = var_accept , linetype = 2) + 
		# xlim(0,0.10) +
		xlab("Prediction Variance") +
		ylab("Density") +
		ggtitle("Prediction Variance Density") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)
	return(p)
}

IPE_density <- function(samples,y_actual) {
	IPE <- data.frame(index = seq(1,length(y_actual),1))
	E <- apply(samples, 2, mean)
	V <- apply(samples, 2, var)
	IPE[ ,2] <- (y_actual - E) / sqrt(V)

	p <- ggplot() + 
		geom_density(data = IPE, aes_string(x = colnames(IPE)[2])) +
		geom_line(data = data.frame(x = seq(-6,6,0.01), den = dnorm(seq(-6,6,0.01))), aes(x,den), linetype = 2, color = "#8b8b8b") +
		# geom_vline(xintercept = .02 , linetype = 2) + 
		xlim(-6,6) +
		xlab("Individual Prediction Error") +
		ylab("Density") +
		ggtitle("Individual Prediction Error Density") +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)
	return(p)
}

cred_intervals <- function(y_draw, y_predict){
	df <- data.frame(index = seq(1,length(y_draw[1,]),1))
	df[ ,2] <- apply(y_draw, 2, var)

	output_percentile <- data.frame(y = y_predict, quantiles = t(apply(y_draw, 2, quantile, probs = c(0.025,0.25,0.50,0.75,0.975))))

	ggplot(output_percentile, aes(x = y)) +   
		geom_boxplot(aes(ymin = quantiles.25., lower = quantiles.25., middle = quantiles.50., upper = quantiles.75., ymax = quantiles.75.), stat = "identity", alpha = 0.15) +
		geom_point(aes(x = y, y = quantiles.50.), size = 0.2) +
		# xlim(0.5,3) + ylim(0.5,3) +
		geom_abline(intercept = 0, slope = 1) +
		xlab("Ariane7 Output") +
		ylab("Emulator Prediction") +
		ggtitle("Emulator Credible Intervals") +
		coord_fixed(ratio = 1) +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


gg_histogram <- function(points, x) {
	ggplot() + 
	geom_histogram(data = data.frame(y = points), aes(x = y, y = ..density..), color = "black", alpha = 0.1, bins = 25) + 
	# geom_histogram(data = data.frame(y = y_predict), aes(x = y, y = ..density..), color = "blue", alpha = 0.5, bins = 25) + 
	xlim(0,3) + ylim(0,1.1) +
	xlab(x) +
	ylab("Histogram Density") +
	# ggtitle("Emulator Credible Intervals") +
	# coord_fixed(ratio = 1) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = 1)
}







