# EM and MM clustering

Clustering en présence des données manquantes MAR et MCAR dans une base de données en utilisant les algorithmes Expectation Maximization (EM) et Minorization Maximization (MM).

Ce projet est focalisé sur la comparaison des performances de certains algorithmes declustering  en  présence  de données  manquantes  dans  notre  échantillon  de  travail.  Plusieursétudes ont été faites sur la possibilité de tenir compte des valeurs manquantes dans la basede donnée pour effectuer un meilleur clustering. Notre étude s’est essentiellement baséesur l’écriture et l’implémentation de l’algorithme EM dans le cas des données manquantes.Nous avons adopté une méthode très particulière pour contrôler le processus générateur desdonnées manquantes via des fonctions de répartition bien connues. Les résultats de l’implé-mentation de cet algorithme nous montrent une croissance continue de la log-vraisemblanceau fur et à mesure que le nombre d’itération augmente. Nous sommes arrivés à la conclusionque les écritures théoriques de l’EM en tenant compte des données manquantes MCAR et MARsont satisfaisantes et peuvent être sujet de comparaison par rapport à d’autres algorithmes declustering gérant également des données manquantes.

## Structure du code du projet
Le code du projet a été mis en place selon la structure d'un package. Nous détaillons les dossiers et fichiers importants du projet.


- **man** :  Ce dossier contient l'ensemble des fichiers markdown de description des différentes fonctions utilisées dans le projet.

- **R** : Ce dossier contient tous les fichiers R. Nous y distinguons 3 fichiers :

  - **mainEM.R** : Ce fichier contient le programme principal de lancement de l'algorithme EM. Cette fonction principale propose plusieurs paramètres qui permettent de générer des échantillons, d'exécuter l'algorithme EM sur l'échantillon selon le nombre de classes en choisissant le langage R ou C++.

  - **MMnonParam.R** : Ce fichier contient l'implémentation de la méthode MM.

  - **RcppExports.R** : Ce dernier est un fichier automatiquement généré par "RStudio". Il contient les appels des      fonctions en R des méthodes implémentées en C++.


- **src**: Contient les fichiers C++ et d'autres fichiers annexes nécessaires  à la bonne compilation du package.

  - **EM.cpp** : Ce fichier contient toutes les méthodes C++ de l'algorithme EM, la simulation des données et l'initialisation des paramètres n'en faisant pas partie car sont gérées par le fichier mainEM.R.



- **DESCRIPTION** : Contient une description succincte  du projet.



## Les méthodes R utilisées
Nous détaillons dans cette section l'utilité des différentes méthodes implémentées en R.


**simulerEchantillon <- function(n, p=2, K, delta)** : Cette fonction permettant de simuler les données de l'échantillon dans un espace $`\mathbb{R}^p`$ avec un nombre de classes prédéfinies.

@param n *integer*. la taille de l'échantillon à simuler 

@param p *integer*. le nombre de covariables à considérer (égal 2 par défaut)

@param K *integer*. Nombre de classes réelles dans l'échantillon ($K < 2^p+1$)

@param delta *double*. Disparité des points de l'échantillon

@return *list*.Retourne la liste des échantillons X, Z, R.




**initialisation <- function(K, ech, seed)** : Cette fonction permet d'initialisation des paramètres de l'algo EM.

@param K *integer*. nombre de clusters

@param ech *list*. L'échantillon complet généré (généralement l'output de la fonction **simulerEchantillon}).

@return *list*. retourne La liste complète des paramètres du modèle.



**singleEM0 <- function(x, K, param, tol)** : Algorithme EM dans le cas des données manquantes avec $\beta_{kj}=0$ (cas MCAR). Cette fonction peut être exécutée sans préalablement utiliser les deux précédentes (i.e. Lorsqu'un utilisateur possède déjà ses propres données et une initialisation des paramètres de son modèle).

@param x *matrix*. L'échantillon utilisé pour le clustering

@param K *integer*. Le nombre de classe à prédire ($K < 2^p+1$)

@param param *list*. Ensemble des paramètres initialisés de l'algo

@param tol *double*. Critère d'arrêt de l'algorithme

@return *list*. Retourne l'ensemble des paramètres\\


**singleEM <- function(x, K, param, tol)** : Algorithme EM dans le cas des données manquantes avec $\beta_{kj} \neq 0$.

@param x *matrix*. L'échantillon utilisé pour le clustering

@param K *integer*. Le nombre de classe à prédire ($K < 2^p+1$)

@param param *list*. Ensemble des paramètres initialisés de l'algo

@param tol *double*. Critère d'arrêt de l'algorithme

@return *list*. Retourne l'ensemble des paramètres\\



**mainEM <- function(n, p, delta, K, epsilon, seed, beta = FALSE, useC = FALSE)** : Cette fonction principale du package. A lancer pour tout éxécuter automatiquement. Lancer les autres fonctions singleEM ou singleEM0 lorsqu'on possède déjà l'échantillon. Cette fonction fait automatiquement appel aux fonctions **simulerEchantillon** et **initialisation**. Elle fait également appel aux fonction C++ lorsque **useC = TRUE**.

@param n *integer*. La taille de l'échantillon à générer

@param p *integer*. Le nombre de covariables

@param delta *double*. Disparité des points dans l'espace

@param K *integer*. Nombre de classes

@param epsilon *double*. Critère d'arrêt de l'algorithme

@param seed *integer*. Graine générateur pour la reproductibilité

@param beta *bool*. False (par défaut) si nous voulons fixer beta à zér

@param useC *bool*. True si nous voulons utiliser le langage C

@return *list*. Retourne l'ensemble des paramètres ainsi que l'évolution de la vraisemblance




## C++ et des libraires d'algèbre linéaire
L'utilisation de C++ a un avantage énorme, dû à sa rapidité d'exécution. En effet, le langage C++ est un langage compilé tandis que le langage R est un langage interprété. C'est d'ailleurs la raison pour laquelle une majeure partie des packages R voient leur code source implémenter en C++.

Durant ce projet, la librairie Eigen a beaucoup été utilisée. Il s'agit d'une librairie d'algèbre linéaire qui permet de facilement manipuler les matrices (multiplication, transposition, etc...). Elle est très flexible et bien documentée. Rstudio met déjà à disposition de ces utilisateurs la possibilité de commencer directement un projet avec RcppEigen qui se base sur la librairie standard de Rcpp. Nous avons également utilisé la librairie RcppNumerical spécialement pour les calculs d'intégrales présentes dans l'algorithme EM.

## Résulats des implémentations

Pour vérifier que notre algorithme fonctionne bien, il est d'usage de faire appel au graphe d'évolution de la vraisemblance. Comme énoncé précédemment, nous utilisons plutôt la log-vraisemblance qui n'est juste qu'une transformation monotone de cette dernière.

L'EM est un algorithme itératif, qui est construit de telle manière que l'augmentation de la valeur de la log-vraisemblance à chaque itération est garantie (i.e. $`l(\Theta^{(t)}) \leq l(\Theta^{(t+1)}) ~~ \forall t \geq 0`$).

Les cas montrent de bons résultats.

