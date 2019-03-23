# define our function of MCMC method
gibbs_metropolis <- function(observations, M, eta_start, k_v_start, k_u_start, 
                             u_start, alpha_u = 1, alpha_v = 1, beta_u=0.01, beta_v=0.01, r_matrix)
{
  n <- length(observations$Y)
  eta <- matrix(nrow = n, ncol = M)
  eta[,1] <- eta_start
  #####################################
  #We make k_v,k_u,u to be global variables
  k_v <<- rep(NA,M) 
  k_u <<- rep(NA,M) 
  test2 <<- rep(NA,M) 
  u <<- matrix(nrow = n, ncol = M)

  v <<- matrix(nrow = n, ncol = M) #we now create a variable for the noise v
  
  k_v[1] <<- k_v_start
  k_u[1] <<- k_u_start
  u[,1] <<- u_start
  v[,1]<<-rnorm(n,0,k_v_start)
  #####################################
  # start the iteration from the 2nd step
  for(i in 2:M)
  {
    # k_u and k_v are drawed from gamma distributions according to the execise 1c
    k_u[i] <<- rgamma(1, shape = (((n-1)/2)+ alpha_u),       #we assign the value to the global variable
                      scale = 1/(beta_u+0.5*(u[,i-1] %*% r_matrix %*% u[,i-1])))
    k_v[i] <<- rgamma(1,shape = (0.5*n + alpha_u),           #we assign the value to the global variable
                      scale = 1/(beta_v+0.5*((eta[,i-1]-u[,i-1]) %*% (eta[,i-1]-u[,i-1]))))
    # u is drawed from normal distribution in Canonical Form 
    u[,i] <<- rmvnorm.canonical(1, b = k_v[i]*eta[,i-1], 
                               Q = (k_v[i] * diag(1,n) + k_u[i] * r_matrix))
    v[,i]<<-rnorm(n,0,k_v[i])
    # calculate the b and c values according to second order Taylor series expansion in execise 1b
    bi <- observations$Y + observations$E * exp(eta[,i-1]) * (eta[,i-1]-1)
    c <- observations$E * exp(eta[,i-1])
    # According to execise 1c, eta is drawed from normal distribution 
    eta[,i] <- rmvnorm.canonical(1, b = (k_v[i]*u[,i] + bi), Q = (k_v[i] + diag.spam(c)))
      # calculate the log density of eta_i
    lden_now <- -0.5*t(eta[,i])%*%(k_v[i]*diag(1,n))%*%eta[,i]+
      t(eta[,i])%*%(k_v[i]*u[,i])+
      t(eta[,i])%*%observations$Y-
      t(exp(eta[,i])) %*% observations$E
    # calculate the log density of the previous eta
    lden_ex <- -0.5*t(eta[,i-1])%*%(k_v[i]*diag(1,n))%*%eta[,i-1]+
      t(eta[,i-1])%*%(k_v[i]*u[,i])+
      t(eta[,i-1])%*%observations$Y-
      t(exp(eta[,i-1])) %*% observations$E
    # calculate the log value of the transfer probability from the previous eta to eta_i
    lq_now <- dmvnorm.canonical(x = eta[,i], 
                                b = (k_v[i]*u[,i] + bi), Q = (k_v[i] + diag.spam(c)))
    # Calculate the previous b and previous c values
    # according to second order Taylor series expansion in execise 1b
    b_ex <- observations$Y + observations$E * exp(eta[,i]) * (eta[,i]-1)
    c_ex <- observations$E * exp(eta[,i])
    # calculate the log value of the transfer probability from eta_i to previous eta
    lq_ex <- dmvnorm.canonical(x = eta[,i-1], 
                               b = (k_v[i]*u[,i] + b_ex), Q = (k_v[i] + diag.spam(c_ex)))
    # calulate the log value of the accept rate
    laccp <- lden_now + lq_ex - lden_ex - lq_now
        # use reject method to accept or reject the sample of eta
    if (log(runif(1)) > laccp)
    {eta[,i] <- eta[,i-1]}
  }
  
  return(eta)
}
