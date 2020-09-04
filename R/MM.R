##################### Fonctions utiles ############################

#' Fonction permettant de simuler les données de l'échantillon
#' dans un espace R^p avec un nombre de classes prédéfini
#'
#' @param n integer. la taille de l'échantillon à simuler 
#' @param d integer. le nombre de covariables à considérer
#' @param K integer. Nombre de classes réelle dans l'échantillon (K < 2^p+1)
#' @param delta double. Disparité des points de l'échantillon
#'
#' @return La liste des échantillons X, Z, R
simulerEchantillon <- function(n, d=2, K, delta){
  
  z <- sample(1:K, n, replace = TRUE) #classes réelles des données
  centers <- permutations(2, d, c(-delta, delta), repeats=TRUE)
  
  x <- matrix(rnorm(n * d), n, d) #Génération de x
  r <- ifelse(matrix(runif(n*d),nrow = n) < 0.1, 0, 1) #10% des données sont manquantes
  x <- ifelse(r, x, NA) # X avec des NA
  
  for (k in 1:K) x[which(z==k),] <- sweep(x[which(z==k),], 2, centers[k,], "+")
  
  return(list('x' = x, 'z' = z, 'r' = r))
}

#' Fonction d'initialisation des paramètres de l'algo MM
#'
#' @param K integer. nombre de classes
#' @param ech list. l'échantillon complet généré
#' @param seed integer. Graine générateur pour la reproductibilité
#'
#' @return list. La liste complète des paramètres du modèle
initialisation <- function(K, ech, seed){
  
  n = nrow(ech$x)
  d = ncol(ech$x)
  
  #Initialisation des poids
  set.seed(seed)
  u  = runif(K)
  Pi = u/sum(u)
  
  #Initialisation de to
  set.seed(seed)
  to = matrix(runif(K*d), nrow = K)
  to = sweep(to, 2, colSums(to), "/")
  
  #Initialisatipn de la matrice de fonctions p par des fonctions de densité de loi normales 
  func = function(u) dnorm(u*n^(1/5))
  l = list()
  nb=1
  for(k in 1:K){
    for(j in 1:d){
      l[[nb]] <- func
      nb = nb+1
    }
  }
  
  # Matrice des fonctions
  p = matrix(l, nrow = K)
  
  return(list('Pi'= Pi, 'p' = p, "to" = to))
}



#' Fonction permettant de calculer des probas conditionnelles d'appartenance 
#' des individus aux classes (est très similaire à l'Etape Expectation de l'EM)
#'
#' @param K integer. nombre de classes
#' @param params list. Liste des paramètres initialisés du modèle
#' @param ech list. Echantillon complet généré
#'
#' @return matrix. La matrice des proba conditionnelles
probCond <- function(K, params, ech){
  n = nrow(ech$x)
  d = ncol(ech$x)
  h = n^(-1/5) #Critère de pouce pour la largeur de la fenêtre
  
  #Calcul des t_ik
  t = matrix(nrow = n, ncol=K)
  
  for(i in 1:n){
    for(k in 1:K){
      # Premier terme dans l'article
      term1 = prod(sapply(1:d, function(j) params$to[k,j]^ech$r[i,j] * (1 - params$to[k,j])^(1-ech$r[i,j])))
      
      res = numeric(d) # Pour stocker les éléments du second terme
      for(j in 1:d){
        f_integral = function(u) (1/h) * dnorm((ech$x[i,j] - u) / h) * params$p[[k,j]](u)
        res[j] = unlist(ifelse(ech$r[i,j]==0, 1, integrate(f_integral, -Inf, Inf)))
      }
      
      term2 = exp(sum(log(res))) #second terme
      
      smoothed_term = term1 * term2 #Terme final qui est le produit des deux termes
      
      t[i,k] = params$Pi[k] * smoothed_term 
    }
    #normalisation
    t[i,] = t[i,]/sum(t[i,])
  }
  
  return(t)
}


#' Fonction de mise à jour des paramètres du modèle
#'
#' @param t matrix. Matrice des probas conditionnelles de l'étape précédente
#' @param params list. Liste des paramètres du modèle de l'étape précédente
#' @param K integer. Nombre de classes 
#' @param ech list. Echantillon complet.
#'
#' @return list. Liste des paramètres du modèle mis à jour
update <- function(t, params, K, ech){
  n = nrow(ech$x)
  d = ncol(ech$x)
  h = n^(-1/5) #Critère de pouce pour la largeur de la fenêtre
  
  # Mise à jour du paramètre Pi (proba à priori des classes)
  params$Pi = colMeans(t) 
  
  #Mise à jour du paramètre to
  params$to = (t(t) %*% ech$r) / colSums(t)
  
  # Mise à jour des fonctions de p
  for(k in 1:K){
    for(j in 1:d){
      params$p[[k,j]] = function(u) sum(sapply(1:n, function(i) ech$r[i,j]*t[i,k]*(1/h)* dnorm((ifelse(is.na(ech$x[i,j]), 0, ech$x[i,j]) - u) / h))) / 
        sum(sapply(1:n, function(i) ech$r[i,j]*t[i,k]))
    }
  }
  
  return(params)
}


#' Fonction qui calcule la log-vraisemblance lissée du modèle
#'
#' @param params list. La liste des paramètres permettant le calcul de la vraisemblance
#' @param K list. Nombre de classes considéré
#' @param ech list. Liste complète de l'échantillon (x, z & r)
#'
#' @return double. Valeur de la valeur de la vraisemblance lissée
log.vraisemblance <- function(params, K, ech){
  n = nrow(ech$x)
  d = ncol(ech$x)
  h = n^(-1/5) #Critère de pouce pour la largeur de la fenêtre
  
  value = 0
  for(i in 1:n){
    inner_sum = 0
    for(k in 1:K){
      #1er terme
      term1 = prod(sapply(1:d, function(j) params$to[k,j]^ech$r[i,j] * (1 - params$to[k,j])^(1-ech$r[i,j])))
      
      res = numeric(d)
      for(j in 1:d){
        f_integral = function(u) (1/h) * dnorm((ech$x[i,j] - u) / h) * params$p[[k,j]](u)
        res[j] = unlist(ifelse(ech$r[i,j]==0, 1, integrate(f_integral, -Inf, Inf)))
      }
      
      #Second terme
      term2 = exp(sum(log(res)))
      
      #cas 2
      smoothed_term = term1 * term2
      
      inner_sum = inner_sum + params$Pi[k] * smoothed_term
    }
    value = value + log(inner_sum)
  }
  
  return(value)
}


#' Fonction globale permettant d'exécuter l'algorithme à partir d'un échantillon donné
#' et d'une initialisation de paramètres.
#'
#' @param ech list. Ensemble de l'échantillon (x & r)
#' @param K integer. Nombre de classes
#' @param param list. Les paramètres initialisé du modèle
#' @param tol double. Critère d'arrêt de l'algorithme
#'
#' @return list. Ensemble des paramètres ainsi que l'évolution de la vraisemblance
MMnonParam <- function(ech, K, param, tol){
  l = log.vraisemblance(param, K, ech) #calul de la vraisemblance à l'étape 0
  evolution = NA
  i = 1
  while(TRUE){ #Itération
    print(paste0(c("Iteration ", i, "..."), collapse = ""))
    t = probCond(K, param, ech)
    param = update(t, param, K, ech)
    l_new = log.vraisemblance(param, K, ech)
    evolution = c(evolution, l_new)
    print(l_new)
    if(abs(l_new - l) < tol) break
    else l = l_new
    i = i+1
  }
  
  plot(evolution, type = 'o', ylab = "Log-Vraisemblance", xlab = "itération EM")
  
  return (list(evolution=evolution, params=params))
}


#' Fonction principale du package. A lancer pour tout éxécuter automatiquement.
#'
#' @param n integer. La taille de l'échantillon à générer
#' @param d integer. Le nombre de covariables
#' @param delta double Disparité des points dans l'espace
#' @param K integer. Nombre de classes
#' @param tol double. Critère d'arrêt de l'algorithme
#' @param seed integer. Graine générateur pour la reproductibilité
#'
#' @return list. Ensemble des paramètres ainsi que l'évolution de la vraisemblance
mainMMnonParam <- function(n, d, K, delta, seed, tol){
  set.seed(seed)
  ech = simulerEchantillon(n, d, K, delta)
  plot(ech$x, col = ech$z, xlab = "X1", ylab = "X2") #Affichage de l'échantillon.
  params = initialisation(K, ech, seed)
  result = MMnonParam(ech, K, params, tol = tol)
  return(result)
}
