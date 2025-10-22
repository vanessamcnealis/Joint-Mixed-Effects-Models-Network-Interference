jmm_fit <- function(data, outcome_formula, propensity_formula, A, alpha.vec){
  
  N <- nrow(data)
  
  m <- length(unique(data$component.id))
  
  
  # Initial values for the regression estimator -----------------------------
  
  
  counterfactual <- data %>% uncount(degree+1) %>% 
    dplyr::select(-k_treated) %>%
    mutate(prop_treated = 0) %>% group_by(id) %>% 
    mutate_at(vars(prop_treated), ~. + 0:(first(degree)))  %>%
    rename(k_treated = prop_treated) %>% 
    mutate(prop_treated = ifelse(degree==0, 0, k_treated/degree)) %>%
    dplyr::select(id, prop_treated, k_treated, degree) %>% 
    rename(counterfactual_prop_treated = prop_treated)
  
  
  
  marginalize_b <- function(z, m, N, beta, sigma_b1, sigma_b2, sigma_e, rho, tau, 
                            Sigma_b, delta.dat, outcome_formula){
    
    b.pred <- mvrnorm(n=m, mu=c(0,0), Sigma = Sigma_b)
    colnames(b.pred) <- c("b1", "b2")
    b.dat <- data.frame(component.id = 1:m, b.pred)
    
    Xm <- model.matrix(object = outcome_formula, data = data)
    Xm <- Xm[,colnames(Xm)!="1 | component.idTRUE"]
    
    Xm_df <- data.frame(Xm)
    Xm_df$id <- data$id
    Xm_df$component.id <-  data$component.id
    Xm_df <- left_join(Xm_df, counterfactual, by = "id")
    Xm_df <- left_join(Xm_df, b.dat, by = "component.id")
    Xm_df <- left_join(Xm_df, delta.dat, by = "id")
    
    Xmz <- Xm_df %>% mutate(prop_treated = counterfactual_prop_treated,
                            Z=z) %>%
      mutate(Z.prop_treated = Z*prop_treated) %>%
      dplyr::select(-c(id, component.id, counterfactual_prop_treated, k_treated, degree, b1, b2,
                       delta)) %>%
      as.matrix()
    colnames(Xmz) <- colnames(Xm)
    
    # Generate potential outcomes from the posterior predictive distribution
    
    pred <- Xmz %*% beta + Xm_df$b2 + Xm_df$delta
    
    return(pred)
  }
  
  pred_reg <- function(mcmc_dat, outcome_formula, B){
    
    beta <- mcmc_dat$beta
    delta <- mcmc_dat$delta
    sigma_b1 <- mcmc_dat$sigma_b1
    sigma_b2 <- mcmc_dat$sigma_b2
    rho <- mcmc_dat$rho
    sigma_e <- mcmc_dat$sigma_e
    tau <- mcmc_dat$tau
    
    # Create covariance matrix for b
    cov <- rho*sigma_b1*sigma_b2
    Sigma_b <- matrix(c(sigma_b1^2, cov, cov, sigma_b2^2), nrow=2, byrow=TRUE)
    
    # block_matrices <- list()
    # for (nu in 1:m){
    #   ind <- data$id[data$component.id == nu]
    #   n_nu <- length(ind)
    #   
    #   # extract the relevant adjacency matrix
    #   A_nu <- A[ind, ind]
    #   
    #   W <- solve(diag(n_nu) - tau*A_nu)
    #   
    #   block_matrices[[nu]]<- (sigma_e)^2*(W%*%t(W))
    # }
    # 
    # Sigma_delta <- bdiag(block_matrices)
    # 
    # delta <- mvrnorm(n = 1, mu = rep(0, N), 
    #                  Sigma = Sigma_delta)
    
    delta.dat <- data.frame(id = data$id, delta=delta)
    
    
    temp0 <- replicate(B, marginalize_b(z=0, m = m, N=N, beta=beta, sigma_b1=sigma_b1, 
                                        sigma_b2 = sigma_b2, rho = rho, tau=tau,
                                        sigma_e=sigma_e, Sigma_b = Sigma_b, delta.dat= delta.dat,
                                        outcome_formula = outcome_formula))
    
    
    
    temp1 <- replicate(B, marginalize_b(z=1, m = m, N=N, beta=beta, sigma_b1=sigma_b1, 
                                        sigma_b2 = sigma_b2, rho = rho, tau=tau,
                                        sigma_e=sigma_e, Sigma_b = Sigma_b, delta.dat= delta.dat,
                                        outcome_formula = outcome_formula))
    
    
    pred.0 <- apply(temp0, 1, mean)
    pred.1 <- apply(temp1, 1, mean)
    
    Xm <- model.matrix(object = outcome_formula, data = data)
    Xm <- Xm[,colnames(Xm)!="1 | component.idTRUE"]
    
    Xm_df <- data.frame(Xm)
    Xm_df$id <- data$id
    Xm_df$component.id <-  data$component.id
    Xm_df <- left_join(Xm_df, counterfactual, by = "id")
    
    reg.0.vec <- numeric(length(alpha.vec))
    reg.1.vec <- numeric(length(alpha.vec))
    
    for (j in 1:length(alpha.vec)){
      Xm_df <- Xm_df %>% mutate(weights = dbinom(x = k_treated, size = degree,
                                                 prob = alpha.vec[j]))
      Xm_df <- Xm_df %>% add_column(pred.0, pred.1)
      temp <- Xm_df %>% group_by(id) %>% dplyr::summarise(reg.0.ind = sum(weights*pred.0),
                                                          reg.1.ind = sum(weights*pred.1))
      temp <- temp %>% left_join(., data %>% dplyr::select(id, component.id), by = "id")
      
      
      temp_by_component <- temp %>% group_by(component.id) %>%
        dplyr::summarise(reg.0.component = mean(reg.0.ind),
                         reg.1.component = mean(reg.1.ind))
      
      reg.0.vec[j] <- mean(temp_by_component$reg.0.component)
      reg.1.vec[j] <- mean(temp_by_component$reg.1.component)
    }
    
    return(c(reg.0.vec, reg.1.vec))
  }
  
  
  temp <- data %>% group_by(component.id) %>% summarise(n.id = n())
  n.id <- temp$n.id
  J <- length(unique(data$component.id))
  start <- numeric(J)
  finish <- numeric(J)
  
  for (j in 1:J){
    start[j] = (ifelse(j==1, 1, sum(n.id[1:(j-1)])+1))
    finish[j] = (sum(n.id[1:j]))
  }
  Xm <- model.matrix(object = outcome_formula, data = data)
  Xm <- Xm[,colnames(Xm)!="1 | component.idTRUE"]
  
  Xe <- model.matrix(object = propensity_formula, data = data)
  Xe <- Xe[,colnames(Xe)!="1 | component.idTRUE"]
  
  Imat <- diag(N)
  
  jags.dat <- list("start" = start, "finish" = finish, "J" = J, "Xm" = Xm,
                   "Xe" = Xe, "Z" = data$Z, "Y" = data$Y,
                   "r" = ncol(Xm), "q" = ncol(Xe),
                   "Mean" = c(0,0), "A" = A, "Imat" = Imat)
  
  n.chains = 4
  n.adapt = 2000
  n.burn = 10000
  n.iter = 30000
  thin = 50
  
  model <- "model
    {
      
      # Stage 1: Likelihood
      for ( j in 1:J){
        
        for( i in start[j]:finish[j]){
          Z[i] ~ dbern(p[i])
          logit(p[i]) = Xe[i,] %*% gamma + b[j,1]
        }
        
        Y[start[j]:finish[j]] ~ dmnorm(mu[start[j]:finish[j]], Phi[j, start[j]:finish[j], start[j]:finish[j]])
        mu[start[j]:finish[j]] <- Xm[start[j]:finish[j],] %*% beta + b[j,2]
        delta[start[j]:finish[j]] <- Y[start[j]:finish[j]] - mu[start[j]:finish[j]] 
        Phi[j,start[j]:finish[j], start[j]:finish[j]] <- (1/(sd_epsilon^2))*t(Imat[start[j]:finish[j], start[j]:finish[j]] - tau*A[(start[j]:finish[j]), (start[j]:finish[j])]) %*% (Imat[start[j]:finish[j], start[j]:finish[j]] - tau*A[start[j]:finish[j], start[j]:finish[j]])  
        }
      
      
      # Stage 2: Prior specification
      
      for (j in 1:J) {
        b[j, 1:2] ~ dmnorm(Mean[1:2],Omega[1:2, 1:2])
      }
      
      # Stage 3
      ### Define the hyperpriors
      Omega <- inverse(Sigma)
      Sigma[1,1] <- sd_b1^2
      Sigma[2,2] <- sd_b2^2
      Sigma[1,2] <- rho*sd_b1*sd_b2
      Sigma[2,1] <- rho*sd_b1*sd_b2
      
      sd_b1 ~ dt(0, 1/(5^2), 1)T(0,)
      sd_b2 ~ dt(0, 1/(5^2), 1)T(0,)
      rho ~ dunif(-1,1)
      tau ~ dunif(-0.99,0.99)
      
      sd_epsilon ~ dt(0, 1/(5^2), 1)T(0,)
      
      for(k in 1:r){
        beta[k] ~ dnorm(0, 0.01)
      }
      
      for(l in 1:q){
        gamma[l] ~ dnorm(0, 0.01)
      }
      
    }"
  
  
  jags.params<-c("beta","gamma", "sd_b1", "sd_b2", "rho", "sd_epsilon", "tau",
                 "delta")
  
  
  jagsfit<-jags.model(textConnection(model), data = jags.dat, 
                      n.chains = n.chains, n.adapt = n.adapt)
  
  
  update(jagsfit, n.burn)
  
  
  samps <- coda.samples(jagsfit, jags.params, n.iter = n.iter, thin=thin)
  
  
  
  output_list1 <-apply(samps[[1]],1,
                       function(x){
                         beta=x[(1):(ncol(Xm))]
                         delta=x[(ncol(Xm)+1):(ncol(Xm)+nrow(data))]
                         gamma=x[(ncol(Xm)+nrow(data)+1):(ncol(Xm)+nrow(data)+ncol(Xe))]
                         rho = x[(ncol(Xm)+nrow(data)+ncol(Xe)+1)]
                         sigma_b1 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+2)]
                         sigma_b2 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+3)]
                         sigma_e = x[(ncol(Xm)+nrow(data)+ncol(Xe)+4)]
                         tau = x[(ncol(Xm)+nrow(data)+ncol(Xe)+5)]
                         list(beta=beta, delta=delta, gamma=gamma, rho = rho, tau = tau,
                              sigma_b1 = sigma_b1,
                              sigma_b2 = sigma_b2,
                              sigma_e=sigma_e)
                       })
  
  output_list2 <-apply(samps[[2]],1,
                       function(x){
                         beta=x[(1):(ncol(Xm))]
                         delta=x[(ncol(Xm)+1):(ncol(Xm)+nrow(data))]
                         gamma=x[(ncol(Xm)+nrow(data)+1):(ncol(Xm)+nrow(data)+ncol(Xe))]
                         rho = x[(ncol(Xm)+nrow(data)+ncol(Xe)+1)]
                         sigma_b1 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+2)]
                         sigma_b2 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+3)]
                         sigma_e = x[(ncol(Xm)+nrow(data)+ncol(Xe)+4)]
                         tau = x[(ncol(Xm)+nrow(data)+ncol(Xe)+5)]
                         list(beta=beta, delta=delta, gamma=gamma, rho = rho, tau = tau,
                              sigma_b1 = sigma_b1,
                              sigma_b2 = sigma_b2,
                              sigma_e=sigma_e)
                       })
  
  output_list3 <-apply(samps[[3]],1,
                       function(x){
                         beta=x[(1):(ncol(Xm))]
                         delta=x[(ncol(Xm)+1):(ncol(Xm)+nrow(data))]
                         gamma=x[(ncol(Xm)+nrow(data)+1):(ncol(Xm)+nrow(data)+ncol(Xe))]
                         rho = x[(ncol(Xm)+nrow(data)+ncol(Xe)+1)]
                         sigma_b1 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+2)]
                         sigma_b2 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+3)]
                         sigma_e = x[(ncol(Xm)+nrow(data)+ncol(Xe)+4)]
                         tau = x[(ncol(Xm)+nrow(data)+ncol(Xe)+5)]
                         list(beta=beta, delta=delta, gamma=gamma, rho = rho, tau = tau,
                              sigma_b1 = sigma_b1,
                              sigma_b2 = sigma_b2,
                              sigma_e=sigma_e)
                       })
  
  output_list4 <-apply(samps[[4]],1,
                       function(x){
                         beta=x[(1):(ncol(Xm))]
                         delta=x[(ncol(Xm)+1):(ncol(Xm)+nrow(data))]
                         gamma=x[(ncol(Xm)+nrow(data)+1):(ncol(Xm)+nrow(data)+ncol(Xe))]
                         rho = x[(ncol(Xm)+nrow(data)+ncol(Xe)+1)]
                         sigma_b1 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+2)]
                         sigma_b2 = x[(ncol(Xm)+nrow(data)+ncol(Xe)+3)]
                         sigma_e = x[(ncol(Xm)+nrow(data)+ncol(Xe)+4)]
                         tau = x[(ncol(Xm)+nrow(data)+ncol(Xe)+5)]
                         list(beta=beta, delta=delta, gamma=gamma, rho = rho, tau = tau,
                              sigma_b1 = sigma_b1,
                              sigma_b2 = sigma_b2,
                              sigma_e=sigma_e)
                       })
  
  
  output_list <- append(output_list1, output_list2)
  output_list <- append(output_list, output_list3)
  output_list <- append(output_list, output_list4)
  
  summary_list<-pblapply(output_list, pred_reg, outcome_formula = outcome_formula, B = 50)
  
  temp <- array(unlist(summary_list), dim=c(length(alpha.vec), 2, n.chains*((n.iter/thin))))
  
  temp2 <- base::apply(temp, MARGIN = c(1,2), mean)
  
  # Posterior distribution of direct effect
  temp22 <- base::apply(temp, MARGIN = 3, function(x) x[,2] - x[,1])
  
  temp3 <- base::apply(temp, MARGIN = c(1,2), var)
  temp4 <- base::apply(temp, MARGIN = c(1,2), quantile, probs = c(0.025))
  temp5 <- base::apply(temp, MARGIN = c(1,2), quantile, probs = c(0.975))
  
  
  # Posterior distribution of the marginal potential outcome
  temp23 <- base::apply(temp, MARGIN = 3, function(x) alpha.vec*x[,2] + (1-alpha.vec)*x[,1])
  
  # Posterior distribution of the overall effect
  temp222 <- array(apply(temp23, 2, function(x) outer(x, x, FUN = "-")), dim = c(length(alpha.vec), length(alpha.vec), 4*(n.iter/thin)))
  
  # Posterior distribution of the indirect effect
  reg00 <- temp[,1,]
  temp333 <- array(apply(reg00, 2, function(x) outer(x, x, FUN = "-")), dim = c(length(alpha.vec), length(alpha.vec), 4*(n.iter/thin)))
  
  # Posterior distribution of total effect
  temp444 <- array(apply(temp, 3, function(x) outer(x[,2], x[,1], FUN = "-")), dim = c(length(alpha.vec), length(alpha.vec), 4*(n.iter/thin)))
  
  mu.0.hat = as.numeric(temp2[,1])
  names(mu.0.hat) <- alpha.vec
  
  mu.1.hat = as.numeric(temp2[,2])
  names(mu.1.hat) <- alpha.vec
  
  mu.hat = base::apply(temp23, MARGIN = 1, mean)
  names(mu.hat) <- alpha.vec

  DE.hat = as.numeric(temp2[,2] - temp2[,1])
  names(DE.hat) <- alpha.vec
  
  IE.hat = apply(temp333, MARGIN = c(1,2), mean)
  rownames(IE.hat) <- alpha.vec
  colnames(IE.hat) <- alpha.vec
  
  OE.hat = apply(temp222, MARGIN = c(1,2), mean)
  rownames(OE.hat) <- alpha.vec
  colnames(OE.hat) <- alpha.vec
  
  TE.hat = apply(temp444, MARGIN = c(1,2), mean)
  rownames(TE.hat) <- alpha.vec
  colnames(TE.hat) <- alpha.vec
  
  var.mu.0 = as.numeric(temp3[,1])
  names(var.mu.0) <- alpha.vec
  
  var.mu.1 = as.numeric(temp3[,2])
  names(var.mu.1) <- alpha.vec
  
  var.mu = apply(temp23, MARGIN = 1, var)
  names(var.mu) <- alpha.vec
  
  var.DE = apply(temp22, MARGIN = 1, var)
  names(var.DE) <- alpha.vec
  
  var.IE = apply(temp333, MARGIN = c(1,2), var)
  rownames(var.IE) <- alpha.vec
  colnames(var.IE) <- alpha.vec
  
  var.OE = apply(temp222, MARGIN = c(1,2), var)
  rownames(var.OE) <- alpha.vec
  colnames(var.OE) <- alpha.vec
  
  var.TE = apply(temp444, MARGIN = c(1,2), var)
  rownames(var.TE) <- alpha.vec
  colnames(var.TE) <- alpha.vec
  
  CrI.mu.0 = t(cbind(temp4[,1], temp5[,1]))
  colnames(CrI.mu.0) <- alpha.vec
  rownames(CrI.mu.0) <- c("2.5%", "97.5%")
  
  CrI.mu.1 = t(cbind(temp4[,2], temp5[,2]))
  colnames(CrI.mu.1) <- alpha.vec
  rownames(CrI.mu.1) <- c("2.5%", "97.5%")
  
  CrI.mu = apply(temp23, MARGIN = 1, quantile, probs = c(0.025, 0.975))
  colnames(CrI.mu) <- alpha.vec
  
  CrI.DE = apply(temp22, MARGIN = 1, quantile, probs = c(0.025, 0.975))
  colnames(CrI.DE) <- alpha.vec
  
  CrI.IE.lb = apply(temp333, MARGIN = c(1,2), quantile, probs = c(0.025))
  colnames(CrI.IE.lb) <- alpha.vec
  rownames(CrI.IE.lb) <- alpha.vec
  
  CrI.IE.ub = apply(temp333, MARGIN = c(1,2), quantile, probs = c(0.975))
  colnames(CrI.IE.ub) <- alpha.vec
  rownames(CrI.IE.ub) <- alpha.vec
  
  CrI.OE.lb = apply(temp222, MARGIN = c(1,2), quantile, probs = c(0.025))
  colnames(CrI.OE.lb) <- alpha.vec
  rownames(CrI.OE.lb) <- alpha.vec
  
  CrI.OE.ub = apply(temp222, MARGIN = c(1,2), quantile, probs = c(0.975))
  colnames(CrI.OE.ub) <- alpha.vec
  rownames(CrI.OE.ub) <- alpha.vec
  
  CrI.TE.lb = apply(temp444, MARGIN = c(1,2), quantile, probs = c(0.025))
  colnames(CrI.TE.lb) <- alpha.vec
  rownames(CrI.TE.lb) <- alpha.vec
  
  CrI.TE.ub = apply(temp444, MARGIN = c(1,2), quantile, probs = c(0.975))
  colnames(CrI.TE.ub) <- alpha.vec
  rownames(CrI.TE.ub) <- alpha.vec
  
  return(list(mu.0.hat=mu.0.hat,
              mu.1.hat=mu.1.hat,
              mu.hat=mu.hat,
              DE.hat=DE.hat,
              IE.hat=IE.hat,
              OE.hat=OE.hat,
              TE.hat=TE.hat,
              var.mu.0=var.mu.0,
              var.mu.1=var.mu.1,
              var.mu=var.mu,
              var.DE=var.DE,
              var.IE=var.IE,
              var.OE=var.OE,
              var.TE=var.TE,
              CrI.mu.0=CrI.mu.0,
              CrI.mu.1=CrI.mu.1,
              CrI.mu = CrI.mu,
              CrI.DE=CrI.DE,
              CrI.IE.lb=CrI.IE.lb,
              CrI.IE.ub=CrI.IE.ub,
              CrI.OE.lb=CrI.OE.lb,
              CrI.OE.ub=CrI.OE.ub,
              CrI.TE.lb=CrI.TE.lb,
              CrI.TE.ub=CrI.TE.ub))
}
