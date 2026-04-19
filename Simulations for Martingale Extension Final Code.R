# MC based proof of concept for extension to the asymptotic distribution 
#     of group-sequentially computed logrank statistic for non-contiguous
#     alternatives utilizing staggered uniform entry,exponential survival, 
#     administrative censoring, under equal random allocation to 
#     treatment/control arm, under the assumption of proportional hazards.

################################################################################
# Required packages
################################################################################

library(parallel)
library(tidyverse)
library(patchwork)
library(KMsurv)
library(survival)
library(mvtnorm)
library(beepr)

################################################################################
# Initial Values
################################################################################

# Number of patients in the trial
n = 300
# Calendar times of analysis; the last time of analysis also specifies the 
#     end of the trial.
t = c(1.5,2,2.5,3)
# Total number of analyses
K = length(t)
# Calendar time of end of accrual 
ta = 2

# Distribution of the time to failure variable; here we assume exponentially
#     distributed survival 
method = "Weibull"
# Hazard rate for the control arm 
haz_ct = 1
lam = haz_ct
par = c(haz_ct,1)

# Log-Hazard Ratio
theta = -0.1

# Probability of allocation to the treatment arm
p = 1/2

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
# Functions to compute the asymptotic means, variances, and covariances
################################################################################

# Hazard rate for the control arm
# This can be made a function to deal with distributions other than Exp
lam0 = lam
# Hazard rate for the treatment arm
lam1 = lam0 * exp(theta)

# E(Z|Y(x))
pi_func = function(x)
{
  (p * pexp(x,rate=lam1,lower.tail=FALSE)) / 
    ( (p * pexp(x,rate=lam1,lower.tail=FALSE)) + 
        ((1-p) * pexp(x,rate=lam0,lower.tail=FALSE)) )
}

# P(Y(x)|Z=1)
q_t_x_1_func = function(x,tt)
{
  punif(tt-x,0,ta) * pexp(x,rate=lam1,lower.tail=FALSE)
}

# P(Y(x)|Z=0)
q_t_x_0_func = function(x,tt)
{
  punif(tt-x,0,ta) * pexp(x,rate=lam0,lower.tail=FALSE)
}

# E(Y(x))
m_func = function(x,tt)
{
  punif(tt-x,0,ta) *
    ( p * pexp(x,rate=lam1,lower.tail=FALSE) + 
        (1-p) * pexp(x,rate=lam0,lower.tail=FALSE) )
}

# E(Z Y(x))
m_1_func = function(x,tt)
{
  p * punif(tt-x,0,ta) * pexp(x,rate=lam1,lower.tail=FALSE)
}

# E((1-Z) Y(x))
m_0_func = function(x,tt)
{
  (1-p) * punif(tt-x,0,ta) * pexp(x,rate=lam0,lower.tail=FALSE)
}

################################################################################

# E[\bar{D}(t)]
d_bar_mean = 0

################################################################################

# V[\bar{D}(t)]
d_bar_var = function(tt)
{
  integrate(f=function(x){
    p * lam1 * q_t_x_1_func(x,tt) * ((1-pi_func(x))^2) +
      (1-p) * lam0 * q_t_x_0_func(x,tt) * ((pi_func(x))^2)
  },lower=0,upper=tt,subdivisions=1000L)$value
}

################################################################################

# E[\widetilde{D}(t)]
d_tilde_mean = function(tt)
{
  integrate(f=function(x){
    sqrt(n) * (exp(theta)-1) * lam0 *
      (1 - pi_func(x)) * m_1_func(x,tt) *
      (1 - (1-(1-m_func(x,tt))^n)/(n*m_func(x,tt)))
  },lower=0,upper=tt,subdivisions=1000L)$value
}

################################################################################

# V[\widetilde{D}(t)]
d_tilde_var = function(tt)
{
  d_tilde_var_func = function(x,y)
  {
    (lam0^2) * m_func(y,tt) * pi_func(y) * (1-pi_func(y)) *
      (pi_func(x)*pi_func(x)*(1+m_func(x,tt)) - 
         pi_func(x)*(2+m_func(x,tt)) + pi_func(y)*(2*pi_func(x)-1) + 1)
  }
  
  d_tilde_var_inner_func = function(y)
  {
    sapply(y,
           function(yy){
             integrate(function(x){d_tilde_var_func(x,yy)},
                       lower=0,upper=yy,subdivisions=1000L)$value
           })
  }
  
  d_tilde_var_val = 2 * ((exp(theta)-1)^2) * 
    integrate(function(y){d_tilde_var_inner_func(y)},
              lower=0,upper=tt,subdivisions=1000L)$value
  
  d_tilde_var_val
}

################################################################################

# Cov( \bar{D}(t) , \widetilde{D}(t) )
d_bar_tilde_cov = function(tt)
{
  d_bar_tilde_cov_func = function(x,y)
  {
    exp(theta) * p * (lam0^2) * ((1-pi_func(y))^2) * (1-pi_func(x)) * 
      pexp(y,rate=lam1,lower.tail=FALSE) * punif(tt-y,0,ta) +
    (1-p) * (lam0^2) * ((pi_func(y))^2) * (-pi_func(x)) * 
      pexp(y,rate=lam0,lower.tail=FALSE) * punif(tt-y,0,ta)
  }
  
  d_bar_tilde_cov_inner_func = function(y)
  {
    sapply(y,
           function(yy){
             integrate(function(x){
               d_bar_tilde_cov_func(x,yy)
             },lower=0,upper=yy,subdivisions=1000L)$value
           })
  }
  
  d_bar_tilde_cov_val = (1-exp(theta)) *
    integrate(function(y){d_bar_tilde_cov_inner_func(y)},
              lower=0,upper=tt,subdivisions=1000L)$value
  
 d_bar_tilde_cov_val
}

################################################################################

# Cov( D(t_1) , D(t_2) )
d_cross_cov = function(tt1,tt2)
{
  d_cross_cov_1_val = d_bar_var(tt1)
  
  ##############################################################################
  
  d_cross_cov_2_func = function(x,y)
  {
    -p * exp(theta) * (lam0^2) * ((1-pi_func(x))^2) * (1-pi_func(y)) * 
      punif(pmin(tt2-x,tt1-y),0,ta) * pexp(x,rate=lam1,lower.tail=FALSE) +
    (1-p) * (lam0^2) * ((pi_func(x))^2) * pi_func(y) * 
      punif(pmin(tt2-x,tt1-y),0,ta) * pexp(x,rate=lam0,lower.tail=FALSE)
  }
  d_cross_cov_2_inner_func = function(x)
  {
    sapply(x,
           function(xx){
             integrate(function(y){
               d_cross_cov_2_func(xx,y)
             },lower=0,upper=pmin(tt1,xx),subdivisions=1000L)$value
           })
  }
  d_cross_cov_2_val = (exp(theta)-1) * 
    integrate(function(x){d_cross_cov_2_inner_func(x)},
              lower=0,upper=tt2,subdivisions=1000L)$value
  
  ##############################################################################
  
  d_cross_cov_3_func = function(x,y)
  {
    -p * exp(theta) * (lam0^2) * ((1-pi_func(x))^2) * (1-pi_func(y)) *
      punif((tt1-x),0,ta) * pexp(x,rate=lam1,lower.tail=FALSE) +
    (1-p) * (lam0^2) * ((pi_func(x))^2) * (pi_func(y)) *
      punif((tt1-x),0,ta) * pexp(x,rate=lam0,lower.tail=FALSE)
      
  }
  d_cross_cov_3_inner_func = function(x)
  {
    sapply(x,
           function(xx){
             integrate(function(y){
               d_cross_cov_3_func(xx,y)
             },lower=0,upper=xx,subdivisions=1000L)$value
           })
  }
  d_cross_cov_3_val = (exp(theta)-1) * 
    integrate(function(x){d_cross_cov_3_inner_func(x)},
              lower=0,upper=tt1,subdivisions=1000L)$value
  
  ##############################################################################  
    
  d_cross_cov_4_func = function(x,y)
  {
     (lam0^2) *
    ( ( ((1-pi_func(x))^2) * ((1-pi_func(y))^2) *
          ( p * punif(pmin(tt1-x,tt2-y),0,ta) * pexp(pmax(x,y),rate=lam1,lower.tail=FALSE) -
              m_1_func(x,tt1) * m_1_func(y,tt2) ) ) -
    ( ((1-pi_func(x))^2) * ((pi_func(y))^2) * m_1_func(x,tt1) * m_0_func(y,tt2) ) -
    ( (pi_func(x)^2) * ((1-pi_func(y))^2) * m_0_func(x,tt1) * m_1_func(y,tt2) ) + 
    ( (pi_func(x)^2) * (pi_func(y)^2) *
      ( p * punif(pmin(tt1-x,tt2-y),0,ta) * pexp(pmax(x,y),rate=lam0,lower.tail=FALSE) -
         m_0_func(x,tt1) * m_0_func(y,tt2) ) ) )
  }
  d_cross_cov_4_inner_func = function(x)
  {
    sapply(x,
           function(xx){
             integrate(function(y){
               d_cross_cov_4_func(xx,y)
             },lower=0,upper=tt2,subdivisions=1000L)$value
           })
  }
  d_cross_cov_4_val = ((exp(theta)-1)^2) *
    integrate(function(x){d_cross_cov_4_inner_func(x)},
              lower=0,upper=tt1,subdivisions=1000L)$value
  
  ##############################################################################
  
  d_cross_cov_1_val + d_cross_cov_2_val + d_cross_cov_3_val + d_cross_cov_4_val
}

################################################################################

# These function compute the asymptotic mean, variance, and covariance 
#     expressions at fixed calendar times.
# To compute the expressions at multiple calendar times we can use apply(...).

# E(D(t))
d_t_asymp_mean = function(tt)
{
  d_tilde_mean(tt)
}

# V(D(t))
d_t_asymp_var = function(tt)
{
  d_bar_var(tt) + d_tilde_var(tt) + 2 * d_bar_tilde_cov(tt)
}

# Cov(D(t_1),D(t_2))
d_t1_t2_asymp_cov = function(tt1,tt2)
{
  d_cross_cov(tt1,tt2)
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# Empirical Values of the mean, variance, and covariance for 
#     group-sequentially logrank statistic

################################################################################
# Function to simulate the data
################################################################################

# Function to generate the survival time and death indicator data structure
#     given the entry time (S), event-time (V) and times of analyses (t).

f_XD = function(S,V,Z,t)
{
  S=c(S)
  V=c(V)
  t=c(t)
  
  K = length(t)
  n = length(S)
  
  # We create an array of K matrices, each of dimension nx2.
  # Within each individual matrix, the first column represents time of entry,
  #     and the second column represents time to death measured from entry.
  dat = array(c(rbind(outer(S,t,function(x,y){y-x}),matrix(rep(V,K),ncol=K))),
              dim=c(n,2,K))
  
  # Measuring the time on test for each individual as measured at each analysis.
  # Each column corresponds to an analysis. 
  req1 = apply(dat,c(1,3),function(x){max(0,min(x[1],x[2]))})
  # Measuring the indicator of death for each individual at each analysis.
  req2 = apply(dat,c(1,3),function(x){as.numeric(x[2]<=x[1])})
  # Creating an array for the treatment indicators
  req3 = matrix(rep(Z,K),ncol=K)
  
  array(c(rbind(req1,req2,req3)),dim=c(n,3,K))
}

f_surv_data = function(n,t,ta,method="Weibull",par,theta)
{
  # n is the required number of subjects to be registered during the accrual
  # t is the vector of times at which the data is analyzed
  # ta is the end of accrual
  # method is the baseline survival distribution
  # par is the vector of parameters for the survival distribution
  # For "Weibull" the first element of par is the rate, and the second element
  #     is the shape parameter
  # Exponential Survival is a special case of the above with shape = 1
  # theta is the log hazard ratio
  
  # Log HR for which we want to generate the data
  n_thet = length(theta)
  
  # Number of interim analyses
  K = length(t)
  
  # Generating the arrival times using uniform distribution of arrival times 
  #     for a poisson process (conditional ...)
  S = runif(n,0,ta)
  
  # Vector of covariates indicating treatment
  Z = numeric(n)
  Z[sample(1:n,ceiling(n/2))] = 1
  
  # Generating survival times from the Cox PH Model
  if(method=="Weibull")
  {
    # Generating required number of uniform(0,1) variates
    U = runif(n)
    # Using the inverse probability transform to generate the survival times
    # V has n_thet columns, each column corresponding to a different log HR
    V = ( diag(-log(U)) %*% ( 1 / (par[1] * exp(Z %o% theta)) ) ) ^ (1/par[2])
  }
  
  # Creating an array to store the generated data.
  # The first dim is the number of subjects, and the second dim is 2 because
  #     each observation is of the type (X_i(t_k),\Delta_i(t_k)), the third
  #     dim is the total number of analyses, and the fourth dim is the 
  #     number of \theta for which the trial is being simulated
  # dat_arr = array(dim=c(n_t,3,K,n_thet))
  dat_arr_XD = array(c(apply(matrix(V,ncol=n_thet),
                             2,function(x){f_XD(S,x,Z,t)})),dim=c(n,3,K,n_thet))
  
  return(dat_arr_XD)
  
}

################################################################################
################################################################################

B = 100000

master_seed = 20250924L
set.seed(master_seed)
seeds = sample.int(.Machine$integer.max, B)

# Generate the raw data

tt = proc.time()
raw_data_list = mclapply(seq_len(B), 
                         function(r){
                           set.seed(seeds[r])
                           AA1 = f_surv_data(n,t,ta,method="Weibull",par,theta)
                           colnames(AA1) = c("X","D","Z")
                           AA1
                         }, mc.cores=13)
proc.time() - tt
beepr::beep(sound = 2)

# Generating the values of \overline{D}(t) and \widetilde{D}(t)

tt = proc.time()
init_d_list = mclapply(seq_len(B),
                       function(r){
                         AA1 = raw_data_list[[r]]
                         
                         LR_stats = c(apply(AA1,c(3,4),
                                          function(mat){
                                            fit <- survdiff(Surv(mat[,1],
                                                                 mat[,2]) ~ mat[,3])
                                            OE  <- fit$obs - fit$exp
                                            -OE[1]
                                          }))
                         
                         LR_stats
                         
                       },mc.cores=13)
proc.time() - tt
beepr::beep(sound = 2)

init_d_mat = do.call(rbind, init_d_list)

rm(raw_data_list)
rm(init_d_list)

################################################################################
################################################################################

# Expressions for the means and variances assuming contiguous alternatives.

# Variance expressions are obtained using expressions provided in Tsiatis 1981

# Hazard rate for the control arm
# This can be made a function to deal with distributions other than Exp
lam0 = lam
# Hazard rate for the treatment arm
lam1 = lam0 * exp(theta)

# Variance of treatment covariate Z
sig2z = p*(1-p)
# Minimum of calendar times of analysis and end of accrual
t.ta.min = pmin(c(t),ta)
# Asymptotic variances at calendar times of analysis
sig2s = sig2z * ( (1-p) * (1/ta) * 
                    (t.ta.min - exp(-lam0*t)*((exp(t.ta.min*lam0)-1)/lam0)) +
                    (p) * (1/ta) * 
                    (t.ta.min - exp(-lam1*t)*((exp(t.ta.min*lam1)-1)/lam1)))

# Mean vector
mean_vec = sqrt(n) * theta * sig2s

################################################################################
################################################################################

# Comparing means, and variances based on MC averages, asymptotic expressions
#     assuming contiguous alternatives, and asymptotic expressions for more
#     general alternatives

#########################################
# Comparison of mean vector
#########################################

# Empirical Mean
apply(init_d_mat,2,function(vec){mean(vec/sqrt(n))})
# Mean assuming general theta
apply(cbind(t),1,function(tt){d_t_asymp_mean(tt)})
# Mean assuming contiguous theta
mean_vec

# Comparison of maximum absolute error for the two mean approximations from the 
#     empirical mean
# Proposed Approximation
max(abs(apply(init_d_mat,2,function(vec){mean(vec/sqrt(n))}) -
          apply(cbind(t),1,function(tt){d_t_asymp_mean(tt)}) ))
# Existing approximation
max(abs(apply(init_d_mat,2,function(vec){mean(vec/sqrt(n))}) - mean_vec ))

# Standard error of Estimated Mean 
# We consider the component-wise standard errors only
apply(init_d_mat,2,function(vec){sd(vec/sqrt(n))/(sqrt(B))})


max(abs(apply(init_d_mat,2,function(vec){mean(vec/sqrt(n))}) -
          apply(cbind(t),1,function(tt){d_t_asymp_mean(tt)})))
#########################################
# Comparison of variances and covariances
#########################################

# Empirical variances
apply(init_d_mat,2,function(vec){var(vec/sqrt(n))})
# Variances assuming general theta 
apply(cbind(t),1,function(tt){d_t_asymp_var(tt)})
# Variances assuming contiguous
sig2s

# Comparison of maximum absolute error for the two variance approximations from 
#     the empirical variance
# Proposed Approximation
max(abs(apply(init_d_mat,2,function(vec){var(vec/sqrt(n))}) -
          apply(cbind(t),1,function(tt){d_t_asymp_var(tt)}) ))
# Existing Approximation
max(abs(apply(init_d_mat,2,function(vec){var(vec/sqrt(n))}) - sig2s ))

# Covariance matrix obtained using asymptotic expressions for general theta
general_var_vec = c(apply(cbind(t),1,function(tt){d_t_asymp_var(tt)}))
general_cov_mat = diag(c(general_var_vec))
general_upper_cov_mat_ind = which(upper.tri(general_cov_mat),arr.ind = TRUE)
general_cov_vec = apply(general_upper_cov_mat_ind,1,
                function(ind){d_t1_t2_asymp_cov(tt1=t[ind[1]],tt2=t[ind[2]])})
general_cov_mat[upper.tri(general_cov_mat)] = general_cov_vec
general_cov_mat[lower.tri(general_cov_mat)] = (t(general_cov_mat))[lower.tri(general_cov_mat)]

# Covariance matrix obtained using asymptotic expressions for contiguous theta
contig_var_vec = sig2s
contig_cov_mat = diag(contig_var_vec)
contig_cov_mat[lower.tri(contig_cov_mat)] = sig2s[rep(c(1:(K-1)),times=c((K-1):1))]
contig_cov_mat[upper.tri(contig_cov_mat)] = (t(contig_cov_mat))[upper.tri(contig_cov_mat)]

# Empirical Covariance Matrix
emp_var_vec = apply(init_d_mat,2,function(vec){var(vec/sqrt(n))})
emp_cov_mat = diag(emp_var_vec)
emp_upper_cov_mat_ind = which(upper.tri(emp_cov_mat),arr.ind = TRUE)
emp_cov_vec = apply(emp_upper_cov_mat_ind,1,
                    function(ind){cov(init_d_mat[,ind[1]]/sqrt(n),init_d_mat[,ind[2]]/sqrt(n))})
emp_cov_mat[upper.tri(emp_cov_mat)] = emp_cov_vec
emp_cov_mat[lower.tri(emp_cov_mat)] = (t(emp_cov_mat))[lower.tri(emp_cov_mat)]

# Maximum Absolute Entry-wise Error for (General - Empirical) covariance matrix
max(general_cov_mat - emp_cov_mat)

# Maximum Absolute Entry-wise Error for (Contiguous - Empirical) covariance matrix
max(contig_cov_mat - emp_cov_mat)

# Frobenius norm for General versus Empirical
sqrt(sum((general_cov_mat - emp_cov_mat)^2))

# Frobenius norm for Contiguous versus Empirical
sqrt(sum((contig_cov_mat - emp_cov_mat)^2))

################################################################################
################################################################################

# Plots using base R

par(mfrow=c(1,3))

hist(init_d_mat[,1]*(1/sqrt(n)),
     freq=FALSE,breaks=21,main="D(t=1.5)/sqrt(n)",xlab="")
lines(density(init_d_mat[,1]*(1/sqrt(n)), na.rm = TRUE), lwd = 2, lty=1)
curve(dnorm(x,mean=d_t_asymp_mean(1.5),sd=sqrt(d_t_asymp_var(1.5))) , 
      add=TRUE, col=2, lwd=2, lty=2)

hist(init_d_mat[,2]*(1/sqrt(n)),
     freq=FALSE,breaks=21,main="D(t=2)/sqrt(n)",xlab="")
lines(density(init_d_mat[,2]*(1/sqrt(n)), na.rm = TRUE), lwd = 2, lty=1)
curve(dnorm(x,mean=d_t_asymp_mean(2),sd=sqrt(d_t_asymp_var(2))) , 
      add=TRUE, col=2, lwd=2, lty=2)

hist((init_d_mat[,1]+init_d_mat[,2])*(1/sqrt(n)),
     freq=FALSE,breaks=21,main="(D(t=1.5)+D(t=2))/sqrt(n)",
     xlab="")
lines(density((init_d_mat[,1]+init_d_mat[,2])*(1/sqrt(n)), 
              na.rm = TRUE), lwd = 2, lty=1)
curve(dnorm(x,mean=d_t_asymp_mean(1.5)+d_t_asymp_mean(2),
            sd=sqrt(d_t_asymp_var(1.5)+d_t_asymp_var(2)+
                      2*d_t1_t2_asymp_cov(1.5,2))) , 
      add=TRUE, col=2, lwd=2, lty=2)

################################################################################
################################################################################

# Plots using ggplot

p1 = ggplot(data.frame(x = (init_d_mat[, 1] * (1 / sqrt(n)))), 
            aes(x = (init_d_mat[, 1] * (1 / sqrt(n))) )) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 21,
                 color = "black",
                 fill = "grey",
                 na.rm = TRUE) +
  geom_density(linewidth = 1.0, na.rm = TRUE) +
  stat_function(fun = dnorm,
                args = list(mean = d_t_asymp_mean(1.5), 
                            sd = sqrt(d_t_asymp_var(1.5))),
                color = "red",
                linewidth = 1.0,
                linetype = "dashed") +
  labs(title = "D(t=1.5)/sqrt(n)",
       x = "",
       y = "Density") +
  theme_bw()

p2 = ggplot(data.frame(x = (init_d_mat[, 2] * (1 / sqrt(n)))), 
            aes(x = (init_d_mat[, 2] * (1 / sqrt(n))) )) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 21,
                 color = "black",
                 fill = "grey",
                 na.rm = TRUE) +
  geom_density(linewidth = 1.0, na.rm = TRUE) +
  stat_function(fun = dnorm,
                args = list(mean = d_t_asymp_mean(2), 
                            sd = sqrt(d_t_asymp_var(2))),
                color = "red",
                linewidth = 1.0,
                linetype = "dashed") +
  labs(title = "D(t=2)/sqrt(n)",
       x = "",
       y = "Density") +
  theme_bw()

p3 = ggplot(data.frame(x = (init_d_mat[,1]+init_d_mat[,2])*(1/sqrt(n))), 
            aes(x = (init_d_mat[,1]+init_d_mat[,2])*(1/sqrt(n)))) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 21,
                 color = "black",
                 fill = "grey",
                 na.rm = TRUE) +
  geom_density(linewidth = 1.0, na.rm = TRUE) +
  stat_function(fun = dnorm,
                args = list(mean = (d_t_asymp_mean(1.5)+d_t_asymp_mean(2)), 
                            sd = sqrt(d_t_asymp_var(1.5)+d_t_asymp_var(2)+
                                      2*d_t1_t2_asymp_cov(1.5,2))),
                color = "red",
                linewidth = 1.0,
                linetype = "dashed") +
  labs(title = "(D(t=1.5)+D(t=2))/sqrt(n)", x = "", y = "Density") +
  theme_bw()

p1+p2+p3

################################################################################
#
################################################################################

theta_vec = -seq(15,1,-1)/10
gen_abs_diff = c(0.0103,0.0099,0.0095,0.0092,0.0087,0.0082,0.0075,0.0070,0.0063,
                 0.0056,0.0047,0.0040,0.0034,0.0023,0.0014)
contig_abs_diff = c(0.3715,0.3233,0.2787,0.2374,0.1995,0.1648,0.1334,0.1049,0.0800,
                    0.0584,0.0404,0.0259,0.0149,0.0076,0.0032)

abs_diff_comp_dat = data.frame(theta_vec,gen_abs_diff,contig_abs_diff)

ggplot(abs_diff_comp_dat , aes(theta_vec)) +
  geom_point(aes(y=gen_abs_diff,color="red")) +
  geom_point(aes(y=contig_abs_diff,color="blue")) +
  labs(x=expression(theta),y="Maximum Absolute Difference") +
  theme_bw() +
  scale_color_discrete(labels = c("Contiguous","General"))+
  labs(color = "Approximation")
