
library(rstan)
library(mvtnorm)
library(circular)
library(ggplot2)
library(gridExtra)
library(reshape)
library(akima)
library(MASS)

source("resub_functions.R")
theme_set(theme_bw())

# Read in synthetic data
data_metocean <- readRDS("data_clean_6000.RDS")

# Read in optimised hyperparameter values
hp_MAP <- readRDS("hp_optim.RDS")

# Label if each variable is linear or circular
types <- c("linear", "linear", "circular", "linear", "linear", "circular", "linear", "circular", "linear", "circular")

df_length <- length(data_metocean[ ,1]) 

N <- 2500
N_predict <- 500

sample_index <- sample(seq(1, df_length, 1), N, replace = FALSE)
prediction_index <- sample(seq(1, df_length - N,1), 500, replace = FALSE)

y_samples <- data_metocean[sample_index,11]
y_predict <- data_metocean[-sample_index, ][prediction_index,11]

X_all <- rbind(data_metocean[sample_index,1:10], data_metocean[-sample_index, ][prediction_index,1:10])

d <- list()
for (i in 1:length(types)) {
	if(types[i] == "circular") {
		options(warn = -1)
		d[[i]] <- as.matrix(dist.circular(X_all[ ,i], method = "geodesic"))
		options(warn = 0)
	} else {
		d[[i]] <- as.matrix(dist(X_all[ ,i]))
	}
}

d_xx <- lapply(d, function(x) x[c(1:N), c(1:N)])
d_xxs <- lapply(d, function(x) x[c(1:N), c((N+1):(N + N_predict))])
d_xsxs <- lapply(d, function(x) x[c((N+1):(N + N_predict)), c((N+1):(N + N_predict))])

k_xx <- k_All(d_xx, hp_MAP)
k_xxs <- k_All(d_xxs, hp_MAP)
k_xsxs <-  k_All(d_xsxs, hp_MAP)

k_xx_inv <- solve(k_xx)
k_xx_chol <- solve(t(chol(k_xx)))

M <- t(k_xxs) %*% k_xx_inv %*% as.matrix(y_samples)
Cov <- k_xsxs - crossprod(k_xx_chol %*% k_xxs)

y_draw <- as.matrix(rmvnorm(1000, mean = M, sigma = Cov))

p1 <- var_scatter(y_draw)
p2 <- coverage_plot(y_draw, y_predict)
p3 <- IPE_plot(y_draw, y_predict)
p4 <- var_density(y_draw)
p5 <- IPE_density(y_draw, y_predict)
p6 <- cred_intervals(y_draw, y_predict, var_accept = 0.1 * var(data_metocean[,11]))

grobs <- list(p1, p2, p3, p4, p5, p6)
lay <- rbind(c(1,4,2),
			 c(3,6,7))
diagnostics <- grid.arrange(grobs = grobs, layout_matrix = lay)








