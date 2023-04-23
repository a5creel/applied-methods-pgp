getBeta_SE <- function(df){
  # --------------------------------------------------------------
  # read in data 
  # --------------------------------------------------------------
  myData <- df %>%
    mutate(education = as.factor(education)) %>%
    mutate(Y = re78) %>%
    select(Y, education)
  
  # --------------------------------------------------------------
  # p_j, probability of X = x_j
  # --------------------------------------------------------------
  p <- summary(myData$education)/length(myData$education)
  
  rows <- length(p)
  
  
  # --------------------------------------------------------------
  # pi_j, \tau percentile of Y, conditional on the value of X
  # --------------------------------------------------------------
  
  myPi <- as.data.frame(matrix(nrow = rows, ncol = 9))
  colnames(myPi) <- c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9")
  rownames(myPi) <- levels(myData$education)
  
  for (i in 1:length(levels(myData$education))) {
    # condition on education level 
    myWorking<- myData %>%
      filter(education == levels(myData$education)[i])
    
    # calculate percentiles of outcome
    myPi[i,] <- round(unname(quantile(myWorking$Y, probs = c(.1, .2, .3, .4, .5, .6, .7, .8, .9))),0)
  }
  
  
  
  # --------------------------------------------------------------
  # Calculate beta
  # --------------------------------------------------------------
  X <- as.numeric(levels(myData$education))
  W <- diag(p)
  matrix_pi <- unname(as.matrix(myPi))
  
  myBeta <- (t(X) %*% W %*% X)^{-1} %*% t(X) %*% W %*% matrix_pi
  
  
  # clean up beta 
  myBeta_t <- as.data.frame(round(myBeta, 0))
  colnames(myBeta_t) <- c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9")
  rownames(myBeta_t) <- c("Beta_t")
  
  # RETURN myBeta_t
  
  # --------------------------------------------------------------
  # Calculate b_jt and t_jt
  # --------------------------------------------------------------
  
  myTau <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
  myN <- table(myData$education)
  
  b_j_tau <- matrix(nrow = rows, ncol = 9)
  t_j_tau <- matrix(nrow = rows, ncol = 9)
  
  
  for (j in 1:rows) {
    for (t in 1:9) {
      b_j_tau[j, t] <- max(1, round(myTau[t]*myN[j]- 1.96 * (myTau[t]*(1-myTau[t])*myN[j])^(1/2), 
                                     digits = 0))
      t_j_tau[j,t]<- min(myN[j], round(myTau[t]*myN[j]+ 1.96 * (myTau[t]*(1-myTau[t])*myN[j])^(1/2), 
                                     digits = 0))
    }
  
  }
  
  # --------------------------------------------------------------
  # Calculating sigma (a rows x 9 matrix)
  # --------------------------------------------------------------
  educ_j <- levels(myData$education)
  
  sigma <- matrix(nrow = rows, ncol = 9)
  
  #calculating y_j(t_j)
  for (j in 1:rows) {
    for (t in 1:9) {
      # calculating y_j()
      y_j <- myData %>%
        filter(education == educ_j[j]) %>%
        select(Y) %>%
        arrange(Y)
      
      y_j_t <- y_j$Y[t_j_tau[j,t]]
      y_j_b <- y_j$Y[b_j_tau[j,t]]
  
      sigma[j, t] <- myN[j]*((y_j_t - y_j_b)/2*1.96)^2
    
    }
  }
  
  # --------------------------------------------------------------
  # Calculating Sigma: nine rowsxrows matrices
  # --------------------------------------------------------------
  Sigma <- vector("list", 9)
  
  for (t in 1:9) {
    Sigma[[t]] <- matrix(nrow = rows, ncol = rows,
                         data = diag(sigma[,t]))
  }
  
  # --------------------------------------------------------------
  # Calculating V
  # --------------------------------------------------------------
  V <- vector(length = 9)
  
  for (t in 1:9) {
    V[t] <- (t(X) %*% X)^{-1} %*% t(X) %*% W %*% Sigma[[t]] %*% W %*% X%*% (t(X) %*% W %*% X)^{-1}
  }
  
  # --------------------------------------------------------------
  # Calculating Delta
  # --------------------------------------------------------------
  #storing each delta in a list
  Delta <- vector("list", 9)
  
  #each item is a matrix that is rowsxrows with values along diagonal
  Delta_vector <- vector(length = rows)
  
  #populating Delta
  for (i in 1:9) {
    for (j in 1:rows) {
      Delta_vector[j] <- (myPi[j, i] - X[j]*myBeta[i])^2/p[j]
    }
    Delta[[i]] <- diag(Delta_vector)
  } 
   
  # --------------------------------------------------------------
  # Calculating D
  # --------------------------------------------------------------
  D <- vector(length = 9)
  
  for (t in 1:9) {
    D[t] <- (t(X) %*% X)^{-1} %*% t(X) %*% W %*% Delta[[t]] %*% W %*% X%*% (t(X) %*% W %*% X)^{-1}
  }
  
  # --------------------------------------------------------------
  # Calculating and reporting SE
  # --------------------------------------------------------------
  mySE <- vector(mode = "numeric", length = 9)
  
  for (t in 1:9) {
    mySE[t] <- sqrt((V[t]+D[t])/sum(myN))
  }
  
  
  # clean up beta 
  mySE_t <- as.data.frame(round(t(mySE), 3))
  colnames(mySE_t) <- c(".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9")
  rownames(mySE_t) <- c("Stand_Error")
  
  return(list(SE = mySE_t, Beta = myBeta_t))
  
}
  
