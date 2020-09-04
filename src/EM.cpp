// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
#include <RcppNumerical.h>

using namespace Numer;

// [[Rcpp::depends(RcppEigen)]]

// GLM coefs with Newton-Raphson. See the link for more info : https://tomroth.com.au/logistic/
// [[Rcpp::export]]
Rcpp::NumericVector glm_coefs(Rcpp::NumericVector y1, Rcpp::NumericMatrix X1, Rcpp::NumericVector w1, double epsilon){
    //Conversion to Eigen Data Structures...
    Eigen::Map<Eigen::MatrixXd> X = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(X1);
    Eigen::Map<Eigen::VectorXd> y = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(y1);
    Eigen::Map<Eigen::VectorXd> w = Rcpp::as<Eigen::Map<Eigen::VectorXd> >(w1); //Vecteurs des pondérations du glm
    
    // Initialisation de beta 0
    Eigen::VectorXd beta(X.cols());
    
    double diff = 1e5; //fixé très grand au départ
    
    while(diff > epsilon){
        
        // Calcul de la proba sachant X et beta
        Eigen::VectorXd p(X.rows());
        Eigen::VectorXd X_Beta(X*beta);
        
        p = exp(X_Beta.array()) / (exp(X_Beta.array()) + 1); // le vecteur 'p'
        
        //La matrice des vecteurs de probablité
        Eigen::MatrixXd W = (p.array()*(1-p.array())).matrix().asDiagonal();
        
        // Changement de beta à l'itération actuelle
        Eigen::VectorXd beta_change = (X.transpose() * W * w.matrix().asDiagonal() * X).inverse() * X.transpose() * w.matrix().asDiagonal() * (y.array() - p.array()).matrix();
        
        //Mise à jour de beta
        beta += beta_change;
        
        //Mise à jour de la différence
        diff = beta_change.array().pow(2).sum();
    }
    
    //Convertion de beta au bon format
    SEXP s = Rcpp::wrap(beta);
    Rcpp::NumericVector beta_hat(s);
    
    return beta_hat;
}


//////////////////////////// Miao algo /////////////////////////////

//// Premier Cas : beta_kj = 0 /////

// [[Rcpp::export]]
Rcpp::List singleEM0_C(int K, Rcpp::List params, Rcpp::NumericMatrix x1, double epsilon){
    //Conversion to Eigen Data Structures...
    Eigen::Map<Eigen::MatrixXd> x = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x1);
    
    // Dimensions de l'échantillon
    int n = x.rows(), d = x.cols();
    
    //Extract all Data Structures from params
    Rcpp::NumericVector pi    = params["Pi"];
    Rcpp::NumericMatrix alpha = params["alpha"];
    Rcpp::NumericMatrix mu    = params["mu"];
    Rcpp::NumericMatrix sigma = params["Sigma"];
    
    //evolution vraisemblance
    std::list<double> evolution;
    
    //Construction de rmat
    Rcpp::NumericMatrix rmat(n,d);
    for(int i=0; i<n; i++){
        for(int j=0; j<d; j++){
            rmat(i,j) = (Rcpp::NumericVector::is_na(x.coeff(i,j)) ? 0 : 1);
        }
    }
    
    //logprobcond
    Eigen::MatrixXd logprobcond(n,K);
    for(int k=0; k<K; k++){
        // Initialisation 
        for(int i=0; i<n; i++){
            logprobcond(i,k) = log(pi(k));
        }
        //MAJ
        for(int j=0; j<d; j++){
            for(int i=0; i<n; i++){
                if(rmat(i,j)==1){ 
                    logprobcond(i,k) += R::dnorm(x(i,j), mu(k,j), sqrt((double) sigma(k,j)), 1) + alpha(k,j)  - log(1 + exp(alpha(k,j)));
                }
                else{ 
                    logprobcond(i,k) +=  log(1/(1 + exp(alpha(k,j))));
                }
            }
        }
    }
    
    Eigen::VectorXd normval(n);
    double maxi;
    for(int i=0; i<n; i++){
        maxi = logprobcond(i,0);
        for(int k=0; k<K; k++){
            maxi = std::max(logprobcond(i,k), maxi);   
        }
        normval(i) = maxi;
    }
    
    //Sweep function    
    for(int i=0; i<n; i++){
        for(int k=0; k<K; k++){
            logprobcond(i,k) = exp(logprobcond(i,k) - normval(i));
        }        
    }
    
    double prec = R_NegInf;
    double loglike = normval.sum() + logprobcond.rowwise().sum().array().log().sum();
    
    int it = 0;
    while(abs(loglike - prec) > epsilon){
        Rcpp::Rcout << "Iteration " << it << "...\n";
        Rcpp::Rcout << "Log_like with C++ :  " << loglike << "\n\n";
        it++;
        
        // E_step
        Eigen::MatrixXd t_ik(n,K);        
        for(int i=0; i<n; i++){
            t_ik.row(i) = logprobcond.row(i).array() / logprobcond.row(i).sum();
        }
        
        // M_step
        pi = t_ik.colwise().sum().array() / n;
        for(int k=0; k<K; k++){
            for(int j=0; j<d; j++){
                Eigen::VectorXd impute(n);
                for(int i=0; i<n; i++){
                    impute(i) = (rmat(i,j)==0 ? mu(k,j) : x(i,j)); 
                }
                
                //Mu et Sigma
                mu(k,j)    = (impute.array() * t_ik.col(k).array()).sum() / t_ik.col(k).sum();
                sigma(k,j) = (pow(impute.array() - mu(k,j), 2) * t_ik.col(k).array()).sum() / t_ik.col(k).sum();
                
                //alpha et beta (faire matrice de 1 en Eigen, puis convertir)
                Eigen::MatrixXd mat; mat.setOnes(n,1);
                SEXP s1 = Rcpp::wrap(mat);
                Rcpp::NumericMatrix reg(s1);
                SEXP s2 = Rcpp::wrap(t_ik.col(k));
                Rcpp::NumericVector weights(s2);
                
                // Faire la regression logistique
                Rcpp::NumericVector coefs = glm_coefs(rmat(Rcpp::_, j), reg, weights, 1e-2); //regression
                alpha(k,j) = coefs(0);
            }
        }
        
        for(int k=0; k<K; k++){
            // Initialisation 
            for(int i=0; i<n; i++){
                logprobcond(i,k) = log(pi(k));
            }
            
            //MAJ
            for(int j=0; j<d; j++){
                for(int i=0; i<n; i++){
                    if(rmat(i,j)==1){ 
                        logprobcond(i,k) += R::dnorm(x(i,j), mu(k,j), sqrt((double) sigma(k,j)), 1) + alpha(k,j) - log(1 + exp(alpha(k,j)));
                    }
                    else{
                        logprobcond(i,k) +=  log(1/(1 + exp(alpha(k,j))));
                    }
                }
                
            }
            
        }
        
        //function apply
        double maxi;
        for(int i=0; i<n; i++){
            maxi = logprobcond(i,0);
            for(int k=0; k<K; k++){
                maxi = std::max(logprobcond(i,k), maxi);   
            }
            normval(i) = maxi;
        }
        
        //Sweep function    
        for(int i=0; i<n; i++){
            for(int k=0; k<K; k++){
                logprobcond(i,k) = exp(logprobcond(i,k) - normval(i));
            }        
        }
        
        double prec = loglike;
        loglike = normval.sum() + logprobcond.rowwise().sum().array().log().sum();
        evolution.push_back(loglike);
    }
    
    return  Rcpp::List::create(Rcpp::Named("pi")= pi,
                               Rcpp::Named("alpha") = alpha,
                               Rcpp::Named("beta")  = "beta est nul dans ce cas",
                               Rcpp::Named("mu")    = mu,
                               Rcpp::Named("sigma") = sigma,
                               Rcpp::Named("evolution") = evolution);
    
}

///// Second cas : beta_kj non nul //////
// Création de la classe "missInt" (missing Integrate) pour intégrer dans le cas beta_kj non nul
// Cette classe est utilisée par RcppNumerical pour le calcul de l'integrale
class missInt: public Func{
private:
    double alpha_kj;
    double beta_kj;
    double mu_kj;
    double sigma_kj;
public:
    missInt(double alpha_kj_, double beta_kj_, double mu_kj_, double sigma_kj_) : alpha_kj(alpha_kj_), beta_kj(beta_kj_), mu_kj(mu_kj_), sigma_kj(sigma_kj_) {}
    
    double operator()(const double& u) const
    {
        return R::dnorm(u, mu_kj, sigma_kj, 0) / (1 + exp(alpha_kj + beta_kj*u));
    }
};


// [[Rcpp::export]]
Rcpp::List singleEM_C(int K, Rcpp::List params, Rcpp::NumericMatrix x1, double epsilon){
    
    //Conversion to Eigen Data Structures...
    Eigen::Map<Eigen::MatrixXd> x = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(x1);
    
    // Dimensions de l'échantillon
    int n = x.rows(), d = x.cols();
    
    //Extract all Data Structures from params
    Rcpp::NumericVector pi    = params["Pi"];
    Rcpp::NumericMatrix alpha = params["alpha"];
    Rcpp::NumericMatrix beta  = params["beta"];
    Rcpp::NumericMatrix mu    = params["mu"];
    Rcpp::NumericMatrix sigma = params["Sigma"];
    
    //evolution vraisemblance
    std::list<double> evolution;
    
    //Construction de rmat
    Rcpp::NumericMatrix rmat(n,d);
    for(int i=0; i<n; i++){
        for(int j=0; j<d; j++){
            rmat(i,j) = (Rcpp::NumericVector::is_na(x.coeff(i,j)) ? 0 : 1);
        }
    }
    
    //logprobcond
    Eigen::MatrixXd logprobcond(n,K);
    for(int k=0; k<K; k++){
        // Initialisation 
        for(int i=0; i<n; i++){
            logprobcond(i,k) = log(pi(k));
        }
        
        //MAJ
        for(int j=0; j<d; j++){
            for(int i=0; i<n; i++){
                if(rmat(i,j)==1){ 
                    logprobcond(i,k) += R::dnorm(x(i,j), mu(k,j), sqrt((double) sigma(k,j)), 1) + alpha(k,j) + beta(k,j)*x(i,j)  - log(1 + exp(alpha(k,j)+beta(k,j)*x(i,j)));
                }
                else{ 
                    // integration
                    missInt f(alpha(k,j), beta(k,j), mu(k,j), sigma(k,j));
                    double err_est;
                    int err_code;
                    double res = integrate(f, R_NegInf, R_PosInf, err_est, err_code);
                    
                    //stocker le résultat
                    logprobcond(i,k) += log(res); 
                }
            }
        }
    }
    
    Eigen::VectorXd normval(n);
    double maxi;
    for(int i=0; i<n; i++){
        maxi = logprobcond(i,0);
        for(int k=0; k<K; k++){
            maxi = std::max(logprobcond(i,k), maxi);   
        }
        normval(i) = maxi;
    }
    
    //Sweep function    
    for(int i=0; i<n; i++){
        for(int k=0; k<K; k++){
            logprobcond(i,k) = exp(logprobcond(i,k) - normval(i));
        }        
    }
    
    double prec = R_NegInf;
    double loglike = normval.sum() + logprobcond.rowwise().sum().array().log().sum();
    
    int it = 0;
    while(abs(loglike - prec) > epsilon){
        Rcpp::Rcout << "Iteration " << it << "...\n";
        Rcpp::Rcout << "Log_like with C++ :  " << loglike << "\n\n";
        it++;
        
        // E_step
        Eigen::MatrixXd t_ik(n,K);        
        for(int i=0; i<n; i++){
            t_ik.row(i) = logprobcond.row(i).array() / logprobcond.row(i).sum();
        }
        
        // M_step
        pi = t_ik.colwise().sum().array() / n;
        for(int k=0; k<K; k++){
            for(int j=0; j<d; j++){
                Eigen::VectorXd impute(n);
                for(int i=0; i<n; i++){
                    impute(i) = (rmat(i,j)==0 ? mu(k,j) : x(i,j)); 
                }
                
                //Mu et Sigma
                mu(k,j)    = (impute.array() * t_ik.col(k).array()).sum() / t_ik.col(k).sum();
                sigma(k,j) = (pow(impute.array() - mu(k,j), 2) * t_ik.col(k).array()).sum() / t_ik.col(k).sum();
                
                //alpha et beta (faire matrice de 1 en Eigen, puis convertir)
                Eigen::MatrixXd mat; mat.setOnes(n,2); mat.col(1) = impute;
                SEXP s1 = Rcpp::wrap(mat);
                Rcpp::NumericMatrix reg(s1);
                SEXP s2 = Rcpp::wrap(t_ik.col(k));
                Rcpp::NumericVector weights(s2);
                
                // Regression logiste à 2 coefficients
                Rcpp::NumericVector coefs = glm_coefs(rmat(Rcpp::_, j), reg, weights, 1e-2); //regression
                alpha(k,j) = coefs(0);
                beta(k,j)  = coefs(1);
            }
        }
        
        for(int k=0; k<K; k++){
            // Initialisation 
            for(int i=0; i<n; i++){
                logprobcond(i,k) = log(pi(k));
            }
            
            //MAJ
            for(int j=0; j<d; j++){
                for(int i=0; i<n; i++){
                    if(rmat(i,j)==1){ 
                        logprobcond(i,k) += R::dnorm(x(i,j), mu(k,j), sqrt((double) sigma(k,j)), 1) + alpha(k,j) + beta(k,j)*x(i,j)  - log(1 + exp(alpha(k,j)+beta(k,j)*x(i,j)));
                    }
                    else{ 
                        // integration
                        missInt f(alpha(k,j), beta(k,j), mu(k,j), sigma(k,j));
                        double err_est;
                        int err_code;
                        double res = integrate(f, R_NegInf, R_PosInf, err_est, err_code);
                        // Rcpp::Rcout << res;
                        logprobcond(i,k) += log(res);
                    }
                }
                
            }
            
        }
        
        //function apply
        double maxi;
        for(int i=0; i<n; i++){
            maxi = logprobcond(i,0);
            for(int k=0; k<K; k++){
                maxi = std::max(logprobcond(i,k), maxi);   
            }
            normval(i) = maxi;
        }
        
        
        //Sweep function    
        for(int i=0; i<n; i++){
            for(int k=0; k<K; k++){
                logprobcond(i,k) = exp(logprobcond(i,k) - normval(i));
            }        
        }
        
        double prec = loglike;
        loglike = normval.sum() + logprobcond.rowwise().sum().array().log().sum();
        evolution.push_back(loglike);
    }
    
    return  Rcpp::List::create(Rcpp::Named("pi")= pi,
                               Rcpp::Named("alpha") = alpha,
                               Rcpp::Named("beta")  = beta,
                               Rcpp::Named("mu")    = mu,
                               Rcpp::Named("sigma") = sigma,
                               Rcpp::Named("evolution") = evolution);
    
}
