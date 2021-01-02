

### Bayesian Semiparametric Longitudinal Drift-Diffusion Mixed Models ###

# Codes accompanying "Bayesian Semiparametric Longitudinal Drift-Diffusion Mixed 
# Models for Tone Learning in Adults".
# Codes written by Giorgio Paulon (giorgio.paulon@utexas.edu), last modified 
# on June 10, 2020, in Austin, TX.

# This file implements a novel generic framework for longitudinal
# functional mixed models that allows automated assessment of an associated 
# predictor’s local time-varying influence. We build on this to develop a novel 
# inverse-Gaussian drift-diffusion mixed model for multi-alternative
# decision-making processes in longitudinal settings. Our proposed model and 
# associated computational machinery make use of B-spline mixtures, hidden 
# Markov models (HMM) and factorial hidden Markov models (fHMM), locally 
# informed Hamming ball samplers etc. to address statistical challenges.


### Input ###

# tau <- vector of response times
# ind <- vector of indicators denoting the participants
# time <- vector of training blocks
# trial <- vector of indicators denoting the trial number within each block
# cens <- vector of censoring indicators
# D <- matrix of covariates X = \{s, d\}
# knots <- knots for the spline basis
# xgrid <- grid where to evaluate the functional parameters
# hypers <- hyperparameters of the MCMC
# Niter <- total number of iterations (default = 5000)
# burnin <- burnin of the chain (default = 2000)
# thin <- thinning factor (default = 5)


### Output ###

# Output is saved in a .RData file and comprises the following objects.

# Z <- posterior samples for the cluster membership labels
# post_mean_delta <- posterior samples for the population-level offset parameters
# post_mean_mu <- posterior samples for the population-level drift curves
# post_mean_b <- posterior samples for the population-level boundary curves
# post_ind_delta <- posterior samples for the individual-level offset parameters
# post_ind_mu <- posterior samples for the individual-level drift curves
# post_ind_b <- posterior samples for the individual-level boundary curves
# sigma2_mu_us <- posterior samples for the smoothness of the drift random effects
# sigma2_mu_ua <- posterior samples for the variance of the drift random effects
# sigma2_b_us <- posterior samples for the smoothness of the boundary random effects
# sigma2_b_ua <- posterior samples for the variance of the boundary random effects
# sigma2_1_mu <- posterior samples for the smoothness of the drift main effects
# sigma2_1_b <- posterior samples for the smoothness of the boundary main effects
# pred_ans <- predicted distribution for the population-level category decision
# pred_time <- predicted distribution for the population-level response times
# pred_ans_ind <- predicted distribution for the individual category decisions
# pred_time_ind <- predicted distribution for the individual response times


### Additional Comments ###

# The current file reproduces the figures in Section 2 and Section 5 of the 
# main manuscript. 




# Load relevant libraries and data ----------------------------------------

# Set the working directory
# setwd('~/Desktop/Code/')

# Load the C++ libraries
library(Rcpp) # 1.0.4.6
library(RcppArmadillo) # 0.9.880.1.0
library(RcppProgress) # 0.4.2
library(rgen) # 0.0.1

# Load the plotting libraries
library(ggplot2) # 3.3.1
library(tidyr) # 1.1.0
library(reshape2) # 1.4.4
library(latex2exp) # 0.4.0
library(RColorBrewer) # 1.1-2
library(plyr) # 1.8.6
library(tidyverse) # 1.3.0

# Load other libraries
library(mvtnorm) # 3.8.2
library(gtools) # 1.1-0
library(LaplacesDemon) # 16.1.4

theme_set(theme_bw(base_size = 16))
cols <- brewer.pal(9, "Set1")

source("R/drift_diff_fcts.R")
sourceCpp('src/drift_diff_fcts.cpp')


data <- read.csv(file = 'data/data.csv')






# Descriptive plots (Section 2 in the manuscript) -------------------------

phonemes <- c('High-level: ā', 'Low-rising: á', 'Low-dipping: ǎ', 'High-falling: à')
data$s <- mapvalues(factor(data$s, levels = 1:4), from = 1:4, to = phonemes)


data_aggr <- data %>%
  mutate(d1 = as.numeric(d == 1), d2 = as.numeric(d == 2), 
         d3 = as.numeric(d == 3), d4 = as.numeric(d == 4)) %>%
  group_by(block, s) %>%
  summarise(d1 = mean(d1), d2 = mean(d2), d3 = mean(d3), d4 = mean(d4))

# Generate Figure 1 (left) in the manuscript
data %>%
  mutate(d1 = as.numeric(d == 1), d2 = as.numeric(d == 2), 
         d3 = as.numeric(d == 3), d4 = as.numeric(d == 4)) %>%
  group_by(ind, block, s) %>%
  summarise(d1 = mean(d1), d2 = mean(d2), d3 = mean(d3), d4 = mean(d4)) %>%
  ggplot() + 
  geom_point(aes(x = jitter(block, 0.5), y = d1, group = ind), alpha = 0.4, col = cols[1]) + 
  geom_line(aes(x = jitter(block, 0.5), y = d1, group = ind), alpha = 0.4, col = cols[1]) + 
  geom_point(aes(x = jitter(block, 0.5), y = d2, group = ind), alpha = 0.4, col = cols[2]) + 
  geom_line(aes(x = jitter(block, 0.5), y = d2, group = ind), alpha = 0.4, col = cols[2]) + 
  geom_point(aes(x = jitter(block, 0.5), y = d3, group = ind), alpha = 0.4, col = cols[3]) + 
  geom_line(aes(x = jitter(block, 0.5), y = d3, group = ind), alpha = 0.4, col = cols[3]) + 
  geom_point(aes(x = jitter(block, 0.5), y = d4, group = ind), alpha = 0.4, col = cols[4]) + 
  geom_line(aes(x = jitter(block, 0.5), y = d4, group = ind), alpha = 0.4, col = cols[4]) + 
  geom_line(aes(x = block, y = d1), col = cols[1], size = 1.5, data = data_aggr) + 
  geom_line(aes(x = block, y = d2), col = cols[2], size = 1.5, data = data_aggr) + 
  geom_line(aes(x = block, y = d3), col = cols[3], size = 1.5, data = data_aggr) + 
  geom_line(aes(x = block, y = d4), col = cols[4], size = 1.5, data = data_aggr) + 
  facet_wrap( ~ s, nrow = 2, ncol = 2) +
  labs(x = 'block', y = 'accuracy') +
  scale_color_brewer(name = "", palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks=seq(1, 10, by = 3)) + 
  theme(legend.position = "none")

# Generate Figure 1 (right) in the manuscript
data %>%
  group_by(block, s, d) %>%
  summarise(count = n(), mean_r_time = mean(r_time)) %>%
  ggplot() + 
  geom_line(aes(x = block, y = mean_r_time, col = factor(d)), size = 1.5) +
  facet_wrap( ~ s, nrow = 2, ncol = 2) +
  labs(x = 'block', y = 'response time') +
  scale_color_brewer(name = "", palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks=seq(1, 10, by = 3)) + 
  theme(legend.position = "none")

data$s <- as.numeric(data$s)



# Run the MCMC ------------------------------------------------------------

# Choose the number of knots (default is between the beginning and the end of 
# the study, at every block)
T_min <- min(data$block)
T_max <- max(data$block)
knots <- T_min:T_max
K <- length(knots)
# Choose a fine grid in order to plot the curves resulting from the spline basis
# expansion
xgrid <- seq(T_min, T_max, length.out = 100)
B <- B_basis(xgrid, knots)
n_ind <- length(unique(data$ind))

# Choose the hyperparameters
hypers <- NULL
hypers$P <- P_smooth1(K + 1)
# hypers$a_sigma <- 100; hypers$b_sigma <- 1

Niter <- 5000
burnin <- 2000
thin <- 5
samp_size <- (Niter - burnin) / thin


# The following lines run the MCMC algorithm. This takes 8 hours on a Dell 
# machine with 16 Gb RAM. Alternatively, the user can load the results from 
# a previous run of the algorithm. In order to do so, leave the boolean variable 
# "run_MCMC" to FALSE below



run_MCMC <- FALSE

{if (run_MCMC == T){
  set.seed(123)
  fit <- loc_clust_fHMM_deltai(tau = data$r_time, 
                               ind = data$ind, 
                               time = data$block, 
                               trial = data$trial, 
                               cens = rep(0, nrow(data)), 
                               D = cbind(data$s, data$d), 
                               knots = knots,
                               xgrid = xgrid, 
                               hypers = hypers, 
                               Niter = Niter, 
                               burnin = burnin, 
                               thin = thin)
  save(fit, file = './results.Rdata')
  
}
  else if (run_MCMC == F){
    load(file = './results.Rdata')
  }}





# Posterior estimates (Section 6 in the manuscript) -----------------------


phonemes <- c('High-level: ā', 'Low-rising: á', 'Low-dipping: ǎ', 'High-falling: à')

post_mean <- melt(apply(fit$post_mean_mu, 1:3, mean))
post_quant <- melt(apply(fit$post_mean_mu, 1:3, quantile, probs = 0.05))
post_mean$Var1 <- rep(xgrid, nrow(post_mean)/length(xgrid))
post_mean$low <- post_quant$value
post_quant <- melt(apply(fit$post_mean_mu, 1:3, quantile, probs = 0.95))
post_mean$upp <- post_quant$value
post_mean$Var2 <- mapvalues(factor(post_mean$Var2, levels = 1:4), 
                            from = 1:4, to = phonemes)
post_mean$Var3 <- mapvalues(factor(post_mean$Var3, levels = 1:4), 
                            from = 1:4, to = phonemes)

# Generate Figure 8 (left) in the manuscript
ggplot(post_mean) + 
  geom_line(aes(x = Var1, y = value, col = factor(Var3)), size = 1.5) + 
  geom_ribbon(aes(x = Var1, ymin = low, ymax = upp, fill = factor(Var3)), alpha = 0.4) + 
  facet_wrap( ~ Var2, nrow = 2, ncol = 2) +
  labs(x = 'block', y = TeX("$\\mu_{d,s}(t)$")) +
  scale_color_brewer(name = "", palette = "Set1") + 
  scale_fill_brewer(name = "", palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, 2) + 
  scale_x_continuous(breaks=seq(1, K, by = 3)) + 
  theme(legend.position = "none")

# Generate Figure S.16 in the Supplementary Materials
ggplot(post_mean[post_mean$Var2 == post_mean$Var3,]) + 
  geom_line(aes(x = Var1, y = value, col = factor(Var3)), size = 1.5) + 
  geom_ribbon(aes(x = Var1, ymin = low, ymax = upp, fill = factor(Var3)), alpha = 0.4) + 
  labs(x = 'block', y = TeX("$\\mu_{d,s}(t)$")) +
  scale_color_brewer(name = "", palette = "Set1") + 
  scale_fill_brewer(name = "", palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, 2) + 
  scale_x_continuous(breaks=seq(1, K, by = 3)) + 
  theme(legend.position = "none")


post_mean <- melt(apply(fit$post_mean_b, 1:3, mean))
post_quant <- melt(apply(fit$post_mean_b, 1:3, quantile, probs = 0.05))
post_mean$Var1 <- rep(xgrid, nrow(post_mean)/length(xgrid))
post_mean$low <- post_quant$value
post_quant <- melt(apply(fit$post_mean_b, 1:3, quantile, probs = 0.95))
post_mean$upp <- post_quant$value
post_mean$Var2 <- mapvalues(factor(post_mean$Var2, levels = 1:4), 
                            from = 1:4, to = phonemes)
post_mean$Var3 <- mapvalues(factor(post_mean$Var3, levels = 1:4), 
                            from = 1:4, to = phonemes)

# Generate Figure 8 (right) in the manuscript
ggplot(post_mean) + 
  geom_line(aes(x = Var1, y = value, col = factor(Var3)), size = 1.5) + 
  geom_ribbon(aes(x = Var1, ymin = low, ymax = upp, fill = factor(Var3)), alpha = 0.4) + 
  facet_wrap( ~ Var2, nrow = 2, ncol = 2) +
  labs(x = 'block', y = TeX("$\\b_{d,s}(t)$")) + 
  scale_color_brewer(name = "", palette = "Set1") + 
  scale_fill_brewer(name = "", palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(1.25, 3.25) + 
  scale_x_continuous(breaks=seq(1, K, by = 3)) + 
  theme(legend.position = "none")


# Plot the predictive distributions for the repsonses (results not shown in the 
# manuscript)
melt_time <- cbind(melt(fit$pred_time), melt(fit$pred_ans)[,4])
names(melt_time) <- c('time','s','iteration','r_time','d')
melt_time$s <- mapvalues(factor(melt_time$s, levels = 1:4), 
                         from = 1:4, to = phonemes)


data %>%
  mutate(d1 = as.numeric(d == 1), d2 = as.numeric(d == 2), 
         d3 = as.numeric(d == 3), d4 = as.numeric(d == 4)) %>%
  group_by(block, s) %>%
  summarise(d1 = mean(d1), d2 = mean(d2), d3 = mean(d3), d4 = mean(d4)) %>%
  ggplot() + 
  geom_point(aes(x = block, y = d1), col = cols[1], alpha = 0.75) + 
  geom_line(aes(x = block, y = d1), col = cols[1], alpha = 0.75) + 
  geom_point(aes(x = block, y = d2), col = cols[2], alpha = 0.75) + 
  geom_line(aes(x = block, y = d2), col = cols[2], alpha = 0.75) + 
  geom_point(aes(x = block, y = d3), col = cols[3], alpha = 0.75) + 
  geom_line(aes(x = block, y = d3), col = cols[3], alpha = 0.75) + 
  geom_point(aes(x = block, y = d4), col = cols[4], alpha = 0.75) + 
  geom_line(aes(x = block, y = d4), col = cols[4], alpha = 0.75) + 
  facet_wrap( ~ s, nrow = 2, ncol = 2) +
  labs(x = 'block', y = 'accuracy',
       title = 'Population empirical accuracies') +
  ylim(0, 1) + 
  scale_x_continuous(breaks=seq(1, K, by = 3)) + 
  scale_color_brewer(name = "", palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")

melt_time %>%
  mutate(d1 = as.numeric(d == 1), d2 = as.numeric(d == 2), 
         d3 = as.numeric(d == 3), d4 = as.numeric(d == 4)) %>%
  group_by(time, s) %>%
  summarise(d1 = mean(d1), d2 = mean(d2), d3 = mean(d3), d4 = mean(d4)) %>%
  ggplot() + 
  geom_point(aes(x = time, y = d1), col = cols[1], alpha = 0.75) + 
  geom_line(aes(x = time, y = d1), col = cols[1], alpha = 0.75) + 
  geom_point(aes(x = time, y = d2), col = cols[2], alpha = 0.75) + 
  geom_line(aes(x = time, y = d2), col = cols[2], alpha = 0.75) + 
  geom_point(aes(x = time, y = d3), col = cols[3], alpha = 0.75) + 
  geom_line(aes(x = time, y = d3), col = cols[3], alpha = 0.75) + 
  geom_point(aes(x = time, y = d4), col = cols[4], alpha = 0.75) + 
  geom_line(aes(x = time, y = d4), col = cols[4], alpha = 0.75) + 
  facet_wrap( ~ s, nrow = 2, ncol = 2) +
  labs(x = 'block', y = 'accuracy',
       title = 'Population predictive accuracies') +
  ylim(0, 1) + 
  scale_x_continuous(breaks=seq(1, K, by = 3)) + 
  scale_color_brewer(name = "", palette = "Set1") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")
# The predictive draws seem to match well with the empirical data!


# Plot the coclustering probabilities
P_coclust <- array(NA, dim = c(6, K))
for (k in 1:K){
  temp <- apply(fit$Z[c(1,6),k:(k+1),], 3, function(x) nrow(unique(x)) == 1)
  P_coclust[1,k] <- mean(temp)
  
  temp <- apply(fit$Z[c(1,11),k:(k+1),], 3, function(x) nrow(unique(x)) == 1)
  P_coclust[2,k] <- mean(temp)
  
  temp <- apply(fit$Z[c(1,16),k:(k+1),], 3, function(x) nrow(unique(x)) == 1)
  P_coclust[3,k] <- mean(temp)
  
  temp <- apply(fit$Z[c(6,11),k:(k+1),], 3, function(x) nrow(unique(x)) == 1)
  P_coclust[4,k] <- mean(temp)
  
  temp <- apply(fit$Z[c(6,16),k:(k+1),], 3, function(x) nrow(unique(x)) == 1)
  P_coclust[5,k] <- mean(temp)
  
  temp <- apply(fit$Z[c(11,16),k:(k+1),], 3, function(x) nrow(unique(x)) == 1)
  P_coclust[6,k] <- mean(temp)
  
}
cl_prob <- data.frame("block" = rep(1:K, each = 6), 
                      "Covariate" = rep(5:0, K), 
                      "value" = as.numeric(P_coclust))
cl_prob$value <- round(cl_prob$value, digits = 2)

# Generate Figure 9 in the manuscript
ggplot() + 
  geom_tile(aes(block, Covariate, fill = value), show.legend = FALSE, 
            color = "white", data = cl_prob) +
  geom_text(aes(block, Covariate, label = value), color = "black", size = 3, 
            data = cl_prob) +
  geom_text(aes(x = 0, y = 5, label = "(1,2)"), color = "black", size = 4) +
  geom_text(aes(x = 0, y = 4, label = "(1,3)"), color = "black", size = 4) +
  geom_text(aes(x = 0, y = 3, label = "(1,4)"), color = "black", size = 4) +
  geom_text(aes(x = 0, y = 2, label = "(2,3)"), color = "black", size = 4) +
  geom_text(aes(x = 0, y = 1, label = "(2,4)"), color = "black", size = 4) +
  geom_text(aes(x = 0, y = 0, label = "(3,4)"), color = "black", size = 4) +
  scale_fill_gradient2(low = cols[1], mid = 'white', high = cols[2], 
                       midpoint = 0.5, limit = c(0, 1)) +
  labs(y = "") + 
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_x_continuous(breaks=1:K) + 
  coord_equal()


# This part needs the MCMC algorithm to have been run. 
# The results uploaded in 'results.Rdata' do not contain these 
# parameters because the file would have been too large to upload.

# Let us now plot the individual specific parameters
# Individual drift parameters
ind_idx <- 20

post_mean <- melt(apply(fit$post_ind_mu[,ind_idx,,,], 1:3, mean))
post_quant <- melt(apply(fit$post_ind_mu[,ind_idx,,,], 1:3, quantile, probs = 0.05))
post_mean$Var1 <- rep(xgrid, nrow(post_mean)/length(xgrid))
post_mean$low <- post_quant$value
post_quant <- melt(apply(fit$post_ind_mu[,ind_idx,,,], 1:3, quantile, probs = 0.95))
post_mean$upp <- post_quant$value
post_mean$Var2 <- mapvalues(factor(post_mean$Var2, levels = 1:4), 
                            from = 1:4, to = phonemes)
post_mean$Var3 <- mapvalues(factor(post_mean$Var3, levels = 1:4), 
                            from = 1:4, to = phonemes)

ind_idx <- 10

post_mean_temp <- melt(apply(fit$post_ind_mu[,ind_idx,,,], 1:3, mean))
post_quant <- melt(apply(fit$post_ind_mu[,ind_idx,,,], 1:3, quantile, probs = 0.05))
post_mean_temp$Var1 <- rep(xgrid, nrow(post_mean_temp)/length(xgrid))
post_mean_temp$low <- post_quant$value
post_quant <- melt(apply(fit$post_ind_mu[,ind_idx,,,], 1:3, quantile, probs = 0.95))
post_mean_temp$upp <- post_quant$value
post_mean_temp$Var2 <- mapvalues(factor(post_mean_temp$Var2, levels = 1:4), 
                                 from = 1:4, to = phonemes)
post_mean_temp$Var3 <- mapvalues(factor(post_mean_temp$Var3, levels = 1:4), 
                                 from = 1:4, to = phonemes)

post_mean <- rbind(post_mean, post_mean_temp)
post_mean$idx <- c(rep(1, nrow(post_mean)/2), rep(2, nrow(post_mean)/2))

# Generate Figure 10(a) in the manuscript
post_mean %>%
  filter(Var2 == Var3) %>%
  ggplot() + 
  geom_line(aes(x = Var1, y = value, col = factor(Var3), group = factor(idx), 
                lty = factor(idx)), size = 1.5) + 
  geom_ribbon(aes(x = Var1, ymin = low, ymax = upp, fill = factor(Var3), 
                  group = factor(idx)), alpha = 0.4) + 
  facet_wrap( ~ Var2, nrow = 2, ncol = 2) +
  labs(x = 'block', y = TeX("$\\mu_{d,s}^{(i)}(t)$")) +
  scale_color_brewer(name = "", palette = "Set1") + 
  scale_fill_brewer(name = "", palette = "Set1") + 
  scale_linetype_manual(values=c(3, 5)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks=seq(1, K, by = 3)) + 
  ylim(0, 4) + 
  theme(legend.position = "none")


# Individual boundary parameters
ind_idx <- 20

post_mean <- melt(apply(fit$post_ind_b[,ind_idx,,,], 1:3, mean))
post_quant <- melt(apply(fit$post_ind_b[,ind_idx,,,], 1:3, quantile, probs = 0.05))
post_mean$Var1 <- rep(xgrid, nrow(post_mean)/length(xgrid))
post_mean$low <- post_quant$value
post_quant <- melt(apply(fit$post_ind_b[,ind_idx,,,], 1:3, quantile, probs = 0.95))
post_mean$upp <- post_quant$value
post_mean$Var2 <- mapvalues(factor(post_mean$Var2, levels = 1:4), 
                            from = 1:4, to = phonemes)
post_mean$Var3 <- mapvalues(factor(post_mean$Var3, levels = 1:4), 
                            from = 1:4, to = phonemes)

ind_idx <- 10

post_mean_temp <- melt(apply(fit$post_ind_b[,ind_idx,,,], 1:3, mean))
post_quant <- melt(apply(fit$post_ind_b[,ind_idx,,,], 1:3, quantile, probs = 0.05))
post_mean_temp$Var1 <- rep(xgrid, nrow(post_mean_temp)/length(xgrid))
post_mean_temp$low <- post_quant$value
post_quant <- melt(apply(fit$post_ind_b[,ind_idx,,,], 1:3, quantile, probs = 0.95))
post_mean_temp$upp <- post_quant$value
post_mean_temp$Var2 <- mapvalues(factor(post_mean_temp$Var2, levels = 1:4), 
                                 from = 1:4, to = phonemes)
post_mean_temp$Var3 <- mapvalues(factor(post_mean_temp$Var3, levels = 1:4), 
                                 from = 1:4, to = phonemes)

post_mean <- rbind(post_mean, post_mean_temp)
post_mean$idx <- c(rep(1, nrow(post_mean)/2), rep(2, nrow(post_mean)/2))

# Generate Figure 10(b) in the manuscript
post_mean %>%
  filter(Var2 == Var3) %>%
  ggplot() + 
  geom_line(aes(x = Var1, y = value, col = factor(Var3), group = factor(idx), 
                lty = factor(idx)), size = 1.5) + 
  geom_ribbon(aes(x = Var1, ymin = low, ymax = upp, fill = factor(Var3), 
                  group = factor(idx)), alpha = 0.4) + 
  facet_wrap( ~ Var2, nrow = 2, ncol = 2) +
  labs(x = 'block', y = TeX("$b_{d,s}^{(i)}(t)$")) +
  scale_color_brewer(name = "", palette = "Set1") + 
  scale_fill_brewer(name = "", palette = "Set1") + 
  scale_linetype_manual(values=c(3, 5)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, 4) + 
  scale_x_continuous(breaks=seq(1, K, by = 3)) + 
  theme(legend.position = "none")


# Additional plots --------------------------------------------------------

ind_idx <- 1
traces <- melt(fit$post_ind_mu[seq(1, 100, by = 11),ind_idx,,,], varnames = c('t','s','d','iteration'))
traces$param <- rep('drift', nrow(traces))
traces$cumulative <- rep(NA, dim = nrow(traces))
for (s in 1:4){
  for (d in 1:4){
    for (t in 1:10){
      traces$cumulative[which((traces$t == t) & (traces$s == s) & (traces$d == d))] <- cummean(traces$value[(traces$t == t) & (traces$s == s) & (traces$d == d)])
    }
  }
}
traces_temp <- melt(fit$post_ind_b[seq(1, 100, by = 11),ind_idx,,,], varnames = c('t','s','d','iteration'))
traces_temp$param <- rep('boundary', nrow(traces_temp))
traces_temp$cumulative <- rep(NA, dim = nrow(traces_temp))
for (s in 1:4){
  for (d in 1:4){
    for (t in 1:10){
      traces_temp$cumulative[which((traces_temp$t == t) & (traces_temp$s == s) & (traces_temp$d == d))] <- cummean(traces_temp$value[(traces_temp$t == t) & (traces_temp$s == s) & (traces_temp$d == d)])
    }
  }
}
traces <- rbind(traces, traces_temp)

# Generate Figure S.7 in the Supplementary Materials
traces %>% 
  dplyr::filter(s == 1,
                d == 1) %>%
  ggplot() + 
  geom_line(aes(x = iteration, y =  value)) + 
  geom_line(aes(x = iteration, y =  cumulative), col = 'red') + 
  facet_grid(rows = vars(param), col = vars(t)) + 
  labs(y = '')



ind_idx <- 1
traces <- melt(fit$post_ind_delta[,ind_idx,], varnames = c('s','iteration'))
traces$cumulative <- rep(NA, dim = nrow(traces))
for (s in 1:4){
  traces$cumulative[which(traces$s == s)] <- cummean(traces$value[( (traces$s == s) )])
}

# Generate Figure S.8 in the Supplementary Materials
traces %>% 
  ggplot() + 
  geom_line(aes(x = iteration, y =  value)) + 
  geom_line(aes(x = iteration, y =  cumulative), col = 'red') + 
  facet_wrap( ~ s, nrow = 1, ncol = 4) +
  labs(y = '')
