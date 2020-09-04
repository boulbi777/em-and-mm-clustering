##################### Fonctions utiles ############################

#' Fonction permettant de simuler les données de l'échantillon
#' dans un espace R^p avec un nombre de classes prédéfini
#'
#' @param n integer. la taille de l'échantillon à simuler 
#' @param p integer. le nombre de covariables à considérer
#' @param K integer. Nombre de classes réelle dans l'échantillon (K < 2^p+1)
#' @param delta double. Disparité des points de l'échantillon
#'
#' @return La liste des échantillons X, Z et R.
simulerEchantillon <- function(n, p=2, K, delta){
  
  z <- sample(1:K, n, replace = TRUE) #classes réelles des données
  centers <- permutations(2, p, c(-delta, delta),repeats=TRUE)
  
  x <- matrix(rnorm(n * p), n, p) #Génération de x
  r <- ifelse(matrix(runif(n*p),nrow = n) < 0.1, 0, 1) #10% des données sont manquantes
  x <- ifelse(r, x, NA) # X avec des NA
  
  for (k in 1:K) x[which(z==k),] <- sweep(x[which(z==k),], 2, centers[k,], "+")
  
  return(list('x' = x, 'z' = z, 'r' = r))
}
# La fonction sweep ici permet de prendre l'échantillon d'un k fixé
# et d'ajouter le centre de ce k à tous l'échantillon sélectionné.
# Ca permet de bien disperser les points dans les clusters donnés.


#' Fonction d'initialisation des paramètres de l'algo EM
#'
#' @param K integer. nombre de classes
#' @param ech list. l'échantillon complet généré
#'
#' @return list. La liste complète des paramètres du modèle
initialisation <- function(K, ech, seed){
  
  p = ncol(ech$x)
  
  #Initialisation des poids
  set.seed(seed)
  u  = runif(K)
  Pi = u/sum(u)
  
  #Initialisation des alphas et betas
  set.seed(seed); alpha = matrix(runif(K*p, min=-1, max = 1), nrow = K)
  set.seed(seed); beta  = matrix(runif(K*p, min=-1, max = 1), nrow = K)
  
  #Initialisation des moyennes (on les sélectionne dans l'échantillon observé)
  set.seed(seed)
  x_obs = na.omit(ech$x)
  mu = x_obs[sample(1:nrow(x_obs), size = K),]
  
  #Sigma sous forme de matrice
  Sigma = matrix(rep(1, K*p), nrow = K)
  
  return(list('Pi'= Pi, 'alpha' = alpha, "beta" = beta, 'mu' = mu, 'Sigma'=Sigma))
}

#' Algorithme EM dans le cas des données manquantes avec beta_kj=0.
#'
#' @param x matrix. L'échantillon utilisé pour le clustering
#' @param K integer. Le nombre de classe à prédire (K < 2^p+1)
#' @param param list. Ensemble des paramètres initialisés de l'algo
#' @param tol double. Critère d'arrêt de l'algorithme
#'
#' @return list. L'ensemble des paramètres
singleEM0 <- function(x, K, param, tol){
  p <- ncol(x)
  n <- nrow(x)
  rmat <- !is.na(x)
  logprobcond <- matrix(0, n, K) #Contient les proba conditionnelles de x_i dans la classe k
  for (k in 1:K){
    logprobcond[,k] <- log(param$Pi[k]) #Initialisation avec les probas à priori
    for (j in 1:p){
      who <- which(rmat[,j]==1)
      
      #Lorsque x_ij est disponible
      logprobcond[who, k] <- logprobcond[who, k] + dnorm(x[who,j], param$mu[k,j], sqrt(param$Sigma[k,j]), log=TRUE) + param$alpha[k,j] - log(1 + exp(param$alpha[k,j]))
      
      #Lorsque x_ij est une donnée manquante
      logprobcond[-who,k] <- logprobcond[-who,k] + log(1/(1+exp(param$alpha[k,j])))
    }
  }
  
  normval <- apply(logprobcond, 1, max)
  logprobcond <- exp(sweep(logprobcond, 1, normval, "-"))
  prec <- -Inf
  loglike <- sum(normval) + sum(log(rowSums(logprobcond))) #Log-vraisemblance actuel
  evolution = NA
  it <- 0 #Itéraions de l'EM
  while ((loglike - prec)>tol){
    print(paste0(c("Iteration ", it, "..."), collapse = ""))
    print(paste0(c("Log-like = ",  loglike), collapse = ""))
    
    it <- it + 1
    
    # Estep
    tik <- logprobcond / rowSums(logprobcond)
    
    # Mstep
    param$Pi <- colSums(tik) / n
    for (k in 1:K){
      for (j in 1:p){
        impute <- x[,j]
        impute[which(rmat[,j]==0)] <- param$mu[k,j]
        param$mu[k,j] <- sum(impute * tik[,k]) / sum(tik[,k])
        param$Sigma[k,j] <- sum( (impute - param$mu[k,j]) * (impute - param$mu[k,j]) * tik[,k]) / sum(tik[,k])
        tmp <- coefficients(glm(rmat[,j] ~ 1, family = "quasibinomial", weights = tik[,k])) #famille 'quasibinomial' pour éviter les warnings
        param$alpha[k,j] <- tmp[1]
      }
    }
    for (k in 1:K){
      logprobcond[,k] <- log(param$Pi[k])
      for (j in 1:p){
        who <- which(rmat[,j]==1)
        logprobcond[who, k] <- logprobcond[who, k] + dnorm(x[who,j], param$mu[k,j], sqrt(param$Sigma[k,j]), log=TRUE) + param$alpha[k,j] - log(1 + exp(param$alpha[k,j]))
        logprobcond[-who,k] <- logprobcond[-who,k] + log(1/(1+exp(param$alpha[k,j])))
      }
    }
    normval <- apply(logprobcond, 1, max)
    logprobcond <- exp(sweep(logprobcond, 1, normval, "-"))
    prec <- loglike
    loglike <- sum(normval) + sum(log(rowSums(logprobcond)))
    evolution = c(evolution, loglike)
  }
  return(list(evolution = evolution, 
              param     = param, 
              loglike   = loglike, 
              prec      = prec, 
              error     = loglike<prec, 
              zhat      =apply(logprobcond, 1, which.max))
  )
}

#' Algorithme EM dans le cas des données manquantes avec beta_kj non nul.
#' Elle est identique à l'algorithme singleEM0 précédant sauf dans les parties commentées.
#'
#' @param x matrix. L'échantillon utilisé pour le clustering
#' @param K integer. Le nombre de classe à prédire (K < 2^p+1)
#' @param param list. Ensemble des paramètres initialisés de l'algo
#' @param tol double. Critère d'arrêt de l'algorithme
#'
#' @return list. L'ensemble des paramètres
singleEM <- function(x, K, param, tol){
  p <- ncol(x)
  n <- nrow(x)
  rmat <- !is.na(x)
  logprobcond <- matrix(0, n, K)
  for (k in 1:K){
    logprobcond[,k] <- log(param$Pi[k])
    for (j in 1:p){
      who <- which(rmat[,j]==1)
      
      # Prise en compte du cas où beta_kj est non nul dans le calcul des probas conditionnelles dans le cas où x_ij est connue
      logprobcond[who, k] <- logprobcond[who, k] + dnorm(x[who,j], param$mu[k,j], sqrt(param$Sigma[k,j]), log=TRUE) + param$alpha[k,j] + param$beta[k,j] * x[who,j] - log(1 + exp(param$alpha[k,j] + param$beta[k,j] * x[who,j]))
      
      # Dans le cas x_ij manquant, un calcul intégral s'impose
      integrand <- function(u) dnorm(u, param$mu[k,j], sqrt(param$Sigma[k,j])) / (1 + exp(param$alpha[k,j] + param$beta[k,j] * u))
      res <- integrate(integrand, lower = -Inf, upper = Inf)
      logprobcond[-who,k] <- logprobcond[-who,k] + log(res$value)
    }
  }
  normval <- apply(logprobcond, 1, max)
  logprobcond <- exp(sweep(logprobcond, 1, normval, "-"))
  prec <- -Inf
  loglike <- sum(normval) + sum(log(rowSums(logprobcond)))
  evolution = NA #evolution du loglike
  it <- 0
  while ((loglike - prec)>tol){
    print(paste0(c("Iteration ", it, "..."), collapse = ""))
    print(paste0(c("Log-like = ",  loglike), collapse = ""))
    
    it <- it + 1
    
    # Estep
    tik <- logprobcond / rowSums(logprobcond)
    
    # Mstep
    param$Pi <- colSums(tik) / n
    for (k in 1:K){
      for (j in 1:p){
        impute <- x[,j]
        impute[which(rmat[,j]==0)] <- param$mu[k,j]
        param$mu[k,j] <- sum(impute * tik[,k]) / sum(tik[,k])
        param$Sigma[k,j] <- sum( (impute - param$mu[k,j]) * (impute - param$mu[k,j]) * tik[,k]) / sum(tik[,k])
        tmp <- coefficients(glm(rmat[,j]~impute, family = "quasibinomial", weights = tik[,k]))
        param$alpha[k,j] <- tmp[1]
        param$beta[k,j] <- tmp[2]
      }
    }
    for (k in 1:K){
      logprobcond[,k] <- log(param$Pi[k])
      for (j in 1:p){
        who <- which(rmat[,j]==1)
        
        #Refaire le même calcul jusqu'à l'arrêt de l'algorithme (cas x_ij connue)
        logprobcond[who, k] <- logprobcond[who, k] + dnorm(x[who,j], param$mu[k,j], sqrt(param$Sigma[k,j]), log=TRUE) + param$alpha[k,j] + param$beta[k,j] * x[who,j] - log(1 + exp(param$alpha[k,j] + param$beta[k,j] * x[who,j]))
        
        # Cas où x_ij est manquant
        integrand <- function(u) dnorm(u, param$mu[k,j], sqrt(param$Sigma[k,j])) / (1 + exp(param$alpha[k,j] + param$beta[k,j] * u))
        res <- integrate(integrand, lower = -Inf, upper = Inf)
        logprobcond[-who,k] <- logprobcond[-who,k] + log(res$value  )
      }
    }
    normval <- apply(logprobcond, 1, max)
    logprobcond <- exp(sweep(logprobcond, 1, normval, "-"))
    prec <- loglike
    loglike <- sum(normval) + sum(log(rowSums(logprobcond)))
    evolution = c(evolution, loglike)
  }
  return(list(evolution = evolution, 
              param     = param, 
              loglike   = loglike, 
              prec      = prec, 
              error     = loglike<prec, 
              zhat      =apply(logprobcond, 1, which.max))
  )
}


################## Main ################## 

#' Fonction principale du package. A lancer pour tout éxécuter automatiquement.
#' Lancer les autres fonctions singleEM ou singleEM0 lorsqu'on possède déjà l'échantillon.
#'
#' @param n integer. La taille de l'échantillon à générer
#' @param p integer. Le nombre de covariables
#' @param delta double Disparité des points dans l'espace
#' @param K integer. Nombre de classes
#' @param epsilon double. Critère d'arrêt de l'algorithme
#' @param seed integer. Graine générateur pour la reproductibilité
#' @param beta bool. False (par défaut) si nous voulons fixer beta à zéro
#' @param useC bool. True si nous voulons utiliser le langage C
#'
#' @return list. Ensemble des paramètres ainsi que l'évolution de la vraisemblance
mainEM <- function(n, p, delta, K, epsilon, seed, beta = FALSE, useC = FALSE){
  
  set.seed(seed) #fixer la graine pour un souci de reproductibilité
  
  # Génération de l'échantillon
  ech = simulerEchantillon(n, p, K, delta)
  
  # Affichage de l'échantillon sur un plan 2D
  plot(ech$x, col = ech$z, xlab = "X1", ylab = "X2")
  
  #Initialisation des paramètres
  print("Initialisation...")
  params = initialisation(K, ech, seed)
  
  if(beta){
    if(useC){
      print("On utilise C dans le cas beta_kj non nul")
      res = singleEM_C(K, params, ech$x, epsilon)
      print(res)
      plot(res$evolution, type = 'o', ylab = "Log-Vraisemblance", xlab = "itération EM", main = "EM avec beta & C")
    } 
    else{
      print("On utilise R dans le cas beta_kj non nul")
      res = singleEM(ech$x, K, params, epsilon)
      print(res)
      plot(res$evolution, type = 'o', ylab = "Log-Vraisemblance", xlab = "itération EM", main = "EM avec beta sur R")
    }
  }
  else{
    if(useC){
      print("On utilise C avec beta_kj=0")
      res = singleEM0_C(K, params, ech$x, epsilon)
      print(res)
      plot(res$evolution, type = 'o', ylab = "Log-Vraisemblance", xlab = "itération EM", main = "EM simple et avec C")
    } 
    else{
      print("On utilise R avec beta_kj=0")
      res = singleEM0(ech$x, K, params, epsilon)
      print(res)
      plot(res$evolution, type = 'o', ylab = "Log-Vraisemblance", xlab = "itération EM", main = "EM simple avec R")
    }
  }
}
