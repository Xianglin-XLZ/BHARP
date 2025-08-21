/*
 * BHARP.cpp
 * 
 * Implementation of the Bayesian Hierarchical Adjustable Random Partition (BHARP) model
 * with reversible jump MCMC for split–merge operations. 
 *
 * 
 * I/O:
 *   - CSVs for beta, Delta, cpre, z, q, tau, w, mu, sigma,
 *     and RJMCMC diagnostics (move type, accept prob, accepted).
 *
 * Dependencies:
 *   - C++14, Rcpp, RcppArmadillo.
 *
 * Exported R function:
 *   - Bharp() : main entry for running MCMC sampling and writing results
 *
 * Author: Xianglin Zhao
 */



// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo)]]


#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <fstream>

using namespace Rcpp;
using namespace std;


////////////////////function definition//////////////////////
// 



// // calculate Omega from z

List getOmega(IntegerVector z,int q) {
  
  List Omega(q);
  for(int t = 1; t <= q; t++) {
    std::vector<int> indices;
    for(int i = 0; i < z.size(); i++) {
      if(z[i] == t) {
        // Adding 1 to make it 1-indexed, similar to R's indexing
        indices.push_back(i + 1);
      }
    }
    Omega[t-1] = indices;
  }
  return Omega;
}

// calculate m from z
IntegerVector getm(IntegerVector z,int q) {
  
  IntegerVector m(q, 0); // Initialize a vector of zeros with length q (for counts)
  
  IntegerVector::iterator p;
  // Count occurrences of each value
  for( p =z.begin(); p < z.end(); p++) {
    // Assuming z contains values that are between 1 and q, inclusive
    if(*p >= 1 && *p <= q) {
      m[*p - 1]++; // Subtract 1 for 0-based indexing in C++
    }
  }
  return m;
}

// rdirichlet 
NumericVector rdirichlet(IntegerVector alpha_m) {
  int size= alpha_m.size();
  NumericVector w (size);
  double sum_ = 0.0;
  for (int t = 0; t < size; t++) {
    double cur = R::rgamma(alpha_m[t], 1.0);
    w[t] = cur;
    sum_ += cur;
  }
  //  normalize
  for (int j = 0; j < size; j++) {
    w[j] /= sum_;
  }
  return(w);
}


//1-dimensional constraint mean
//1T %*% D %*% 1 = sum(diag)
//D %*% 1 %*% 1T %*% theta [x]= diag[x]*sum(theta)
arma::vec cond1mean(const arma::vec& theta, const arma::vec& diag) {
  double sumTheta = arma::sum(theta);
  double sumDiag = arma::sum(diag);
  return theta - (diag * sumTheta) / sumDiag;
}

//1-dimensional constraint var
//1T %*% D %*% 1 = sum(diag)
//D %*% 1 %*% 1T %*% D [x,y]= diag[x]*diag[y]
arma::mat cond1Var(const arma::vec& theta, const arma::vec& diag) {
  //int n = diag.n_elem;
  arma::mat D = arma::diagmat(diag); // Convert diag vector to diagonal matrix D
  //arma::vec oneVec = arma::ones<arma::vec>(n);
  double sumDiag = arma::sum(diag);
  // Calculate variance using matrix operations
  //arma::mat Var = D - (D * oneVec * oneVec.t() * D) / sumDiag;
  arma::mat Var = D - diag*diag.t()/ sumDiag;
  return Var;
}

//sample controlling one-dimensional sum=0
//sample the first n-1 elements and control the last one to be -sum of the previous
NumericVector cond1samp(const NumericVector& theta, const NumericVector& diag){
  int n = theta.length();
  // convert to armadillo
  arma::vec theta_arma = Rcpp::as<arma::vec>(theta);
  arma::vec diag_arma = Rcpp::as<arma::vec>(diag);
  //calculate mean and variance
  arma::vec mean=cond1mean(theta_arma,diag_arma);
  arma::mat var=cond1Var(theta_arma,diag_arma);
  
  // the first n-1 elements
  arma::vec mean_adj = mean.head(n-1);
  arma::mat var_adj = var.submat(0, 0, n-2, n-2);
  
  // Sample n-1 elements from a multivariate normal distribution
  arma::vec sample = mean_adj + arma::chol(var_adj) * arma::randn(n-1);
  
  arma::vec result(n);
  result.head(n-1) = sample;
  result(n-1) = -arma::accu(sample);
  
  return wrap(result);
}


//calculate P(q+1 merge)/P( q split)
double PmoveRatio( const int q, const int s){
  if (q==1 && q+1<s) {return 0.40/1.0;}
  else if (q+1==s && q!=1) {return 1.0/0.60;}
  else if (q+1==s && q==1) {return 1.0;}
  else {return 0.40/0.60;}
}

//calculate power prior for q
NumericVector power_prior_seq(int s) {
  NumericVector result(s);
  double sum = 0.0;
  
  // Compute squares and their sum
  for (int i = 0; i < s; ++i) {
    result[i] = pow(i+1,2);
    sum += result[i];
  }
  
  // Normalize by the sum
  for (int i = 0; i < s; ++i) {
    result[i] /= sum;
  }
  
  return result;
}

void write_csv(const Rcpp::NumericMatrix& mat, const std::string& outdir, const std::string& varname, int datasetid, int chainid) {
  std::string filename = outdir + "/" +varname+"_"+std::to_string(datasetid)+"-" + std::to_string(chainid)+".csv";
  std::ofstream file(filename.c_str());
  if (!file.is_open())
    Rcpp::stop("Failed to open file.");
  
  for (int i = 0; i < mat.nrow(); i++) {
    for (int j = 0; j < mat.ncol(); j++) {
      file << mat(i, j);
      if (j < mat.ncol() - 1) file << ",";
    }
    file << "\n";
  }
  file.close();
}
void write_csv(const Rcpp::IntegerMatrix& mat, const std::string& outdir, const std::string& varname, int datasetid, int chainid) {
  std::string filename=outdir+"/"+varname+"_"+std::to_string(datasetid)+"-" + std::to_string(chainid)+".csv";
  std::ofstream file(filename.c_str());
  if (!file.is_open())
    Rcpp::stop("Failed to open file.");
  
  for (int i = 0; i < mat.nrow(); i++) {
    for (int j = 0; j < mat.ncol(); j++) {
      file << mat(i, j);
      if (j < mat.ncol() - 1) file << ",";
    }
    file << "\n";
  }
  file.close();
}
void write_csv(const Rcpp::NumericVector& vec, const std::string& outdir, const std::string& varname, int datasetid, int chainid) {
  std::string filename=outdir+"/"+varname+"_"+std::to_string(datasetid)+"-" + std::to_string(chainid)+".csv";
  std::ofstream file(filename.c_str());
  if (!file.is_open())
    Rcpp::stop("Failed to open file.");
  
  for (int i = 0; i < vec.size(); i++) {
    file << vec[i] << "\n";
  }
  file.close();
}
void write_csv(const Rcpp::IntegerVector& vec, const std::string& outdir, const std::string& varname, int datasetid, int chainid) {
  std::string filename=outdir+"/"+varname+"_"+std::to_string(datasetid)+"-" + std::to_string(chainid)+".csv";
  std::ofstream file(filename.c_str());
  if (!file.is_open())  Rcpp::stop("Failed to open file.");
  
  for (int i = 0; i < vec.size(); i++) {
    file << vec[i] << "\n";
  }
  file.close();
}



//helper function to flatten Matrix row-major
NumericVector MatToVecByRow(NumericMatrix mat) {
  int n = mat.nrow();
  int m = mat.ncol();
  NumericVector vec(n * m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      vec[i * m + j] = mat(i, j);
    }
  }
  return vec;
}
IntegerVector MatToVecByRow(IntegerMatrix mat) {
  int n = mat.nrow();
  int m = mat.ncol();
  IntegerVector vec(n * m);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      vec[i * m + j] = mat(i, j);
    }
  }
  return vec;
}

//helper function to flatten vector of vector
NumericVector VecOfVecToVec(int unit_length, std::vector<NumericVector> vec_vec) {
  int n_vectors = vec_vec.size();
  NumericVector result( n_vectors*unit_length); // Initialize the result vector with NA
  std::fill(result.begin(), result.end(), NA_REAL);
  
  for (int i = 0; i < n_vectors; i++) {
    NumericVector vec = vec_vec[i];
    int start = i*unit_length;
    for (int k = 0; k < vec.size(); k++){
      result[start+k] = vec[k];
    }
  }
  return result;
}


//function to calculate the probability to allocate to two clusters
//used log-sum-exp trick
NumericVector getLogProbt1t2(double wt1, double wt2, double mut1, double mut2, double sigmat1, double sigmat2,
                             double Delta_ik){
  NumericVector LogProbt1t2(2); //the logprobability of assigning to t1 and t2
  //P(z[k]=t)=w[t]*dnorm(Delta[ik], mu[t], sqrt(sigma[t]))
  
  //P1 = wt1 * sqrt(1.0/sigmat1) * exp( -0.5 * pow( Delta_ik-mut1 , 2.0)/ sigmat1);
  double LogP1 = std::log(wt1) - 0.5*std::log(sigmat1)
    - 0.5* std::pow( Delta_ik-mut1 , 2.0)/sigmat1;
  //P2 = wt2 * sqrt(1.0/sigmat2) * exp( -0.5 * pow( Delta_ik-mut2 , 2.0)/ sigmat2);
  double LogP2 = std::log(wt2) - 0.5*std::log(sigmat2)
    - 0.5* std::pow( Delta_ik-mut2 , 2.0)/sigmat2;
  double M = std::max(LogP1, LogP2);
  
  double sum_exp = std::exp(LogP1 - M) + std::exp(LogP2 - M);
  double log_sum_exp = M + std::log(sum_exp);
  
  LogProbt1t2[0] = LogP1-log_sum_exp;
  LogProbt1t2[1] = LogP2-log_sum_exp;
  return LogProbt1t2;
}



/////////////////////////----class definition----///////////////////////////
class Model{
public:
  int s;                    // number of grps, ncol for Delta
  int m;                    // number of arms, nrow for Delta
  double cpre;              // cell precision (assume homoskedastic)
  NumericVector beta;       //beta1,...,betam, arm average
  NumericMatrix Delta;      // shift from betas, row sum = 0
  
  IntegerVector iq;                  //q:npeak
  std::vector<NumericVector> iw;     //peak weights length q
  IntegerMatrix iz;                  //z: allocation vector. z_k length s
  
  std::vector<List> iOmega;          //Omega: length q. Omega_t: the cluster 
  std::vector<IntegerVector> im;     //m: the sample size of each peak. length q
  NumericVector itau;                //tau:between peak precision
  std::vector<NumericVector> imu;    //peak center length q
  std::vector<NumericVector> isigma; //peak variance length q
  
  IntegerVector iMove;        //type of move
  NumericVector iP_Accept;    //accept rate
  IntegerVector iAccepted;    //decision of accepting the move
  
  
  
  void update( const IntegerVector & Data_i, const IntegerVector& Data_k, const NumericVector& Data_Y,
               NumericVector c_beta, NumericVector p_beta, 
               double a_cell, double b_cell, double a_between, double b_between, double a_within, double b_within);
  void split(int i, double a_within, double b_within, 
             double au1, double au2, double au3);
  void merge(int i, double a_within, double b_within, 
             double au1, double au2, double au3);
  void print() const;
  
  
  Model &operator=(const Model& Mprev){
    if (this!=&Mprev){
      //copy previous model
      this->s = (Mprev.s);
      this->m = (Mprev.m);
      this->cpre = (Mprev.cpre);
      
      this->beta = Rcpp::clone(Mprev.beta);
      this->Delta = Rcpp::clone(Mprev.Delta);
      
      this->iq = Rcpp::clone(Mprev.iq);
      
      // Deep copy std::vector<NumericVector> iw
      this->iw.resize(Mprev.iw.size());
      for (size_t i = 0; i < Mprev.iw.size(); ++i) {
        this->iw[i] = Rcpp::clone(Mprev.iw[i]);
      }
      
      this->iz = Rcpp::clone(Mprev.iz);
      
      // Deep copy std::vector<List> iOmega
      this->iOmega.resize(Mprev.iOmega.size());
      for (size_t i = 0; i < Mprev.iOmega.size(); ++i) {
        this->iOmega[i] = Rcpp::clone(Mprev.iOmega[i]);
      }
      
      // Deep copy std::vector<IntegerVector> im
      this->im.resize(Mprev.im.size());
      for (size_t i = 0; i < Mprev.im.size(); ++i) {
        this->im[i] = Rcpp::clone(Mprev.im[i]);
      }
      
      this->itau = Rcpp::clone(Mprev.itau);
      
      // Deep copy std::vector<NumericVector> imu and isigma
      this->imu.resize(Mprev.imu.size());
      this->isigma.resize(Mprev.isigma.size());
      for (size_t i = 0; i < Mprev.imu.size(); ++i) {
        this->imu[i] = Rcpp::clone(Mprev.imu[i]);
        this->isigma[i] = Rcpp::clone(Mprev.isigma[i]);
      }
      this->iMove =Rcpp::clone(Mprev.iMove);
      this->iP_Accept = Rcpp::clone(Mprev.iP_Accept);
      this->iAccepted = Rcpp::clone(Mprev.iAccepted);

    }
    return*this;
  }
  
  //Default constructor
  Model() = default;
  
  //Constructor initialize
  Model(int ngrp, int narm,
        NumericVector c_beta, NumericVector p_beta,
        double a_cell, double b_cell, double a_between, double b_between, double a_within, double b_within) 
    : s(ngrp), m(narm){
    // norm (mean, sd)
    // gamma (shape, scale) scale=1/rate
    // Vector<INTSXP> sample(int n, int size, bool replace = false, sugar::probs_t probs = R_NilValue, bool one_based = true);
    
    // Initialize vectors with proper sizes
    iq = IntegerVector(m);                        // Size m
    iw = std::vector<NumericVector>(m);           // m vectors, one for each arm
    iz = IntegerMatrix(m, s);                     // Matrix with dimensions m x s
    iOmega = std::vector<List>(m);                // m Lists
    im = std::vector<IntegerVector>(m);           // m IntegerVectors
    itau = NumericVector(m);                      // Size m
    imu = std::vector<NumericVector>(m);          // m NumericVectors
    isigma = std::vector<NumericVector>(m);       // m NumericVectors
    iMove = std::vector<unsigned int>(m, 9);      // Initialize to 0
    iP_Accept = std::vector<double>(m, -9.0);      // Initialize to 0.0
    iAccepted = std::vector<unsigned int>(m, 9);  // Initialize to 0
    
    
    cpre = R::rgamma( a_cell, 1.0/b_cell);            //cell precision
    beta = NumericVector(m);
    Delta = NumericMatrix(m, s);                      //create Delta_{ik}
    
    for (int i=0; i<m;i++){
      beta[i] = R::rnorm(c_beta[i], sqrt(1.0/p_beta[i]));     //create beta1...betam
      iq[i] = sample(s,1,TRUE,power_prior_seq(s))[0];                                 //sample number of component from discrete uniform
      iw[i] = rdirichlet(IntegerVector(iq[i], 1));            //sample weight from Dirichlet
      iz(i,_) = sample(iq[i], s, TRUE, iw[i]);                //sample assignment vector with weight
      iOmega[i] = getOmega(iz(i,_), iq[i]);                   // grps in each cluster
      im[i] = getm(iz(i,_), iq[i]);                           // ngrps in each cluster
      itau[i] = R::rgamma(a_between, 1.0 / b_between);                    // between cluster precision
      imu[i] = rnorm(iq[i], 0, sqrt(1.0 / itau[i]));                      // component centers
      isigma[i] = 1.0 / Rcpp::rgamma(iq[i], a_within, 1.0 / b_within);    // within component variances
      for (int k=0; k<s; k++){                                            // Delta
        int t = iz(i,k)-1;                                          //component id 0-based
        Delta(i,k) = R::rnorm(imu[i][t], sqrt(isigma[i][t]));      //sample Delta from t-th component
      }
      
    }
  }
  
  // Copy constructor
  Model(const Model& Mprev) {
    *this = Mprev;
  }
  
}; //class Model end 

void Model::print() const{
  Rcpp::Rcout << "Model Information:\n";
  Rcpp::Rcout << "Number of subgroups (s): " << s << "\n";
  Rcpp::Rcout << "Number of arms (m): " << m << "\n";
  Rcpp::Rcout << "Cell precision (cpre):"<<cpre<<"\n";
  
  
  Rcpp::Rcout << "\n Arm Mean Beta  (size " << beta.size() << "):\n";
  for (int i = 0; i < beta.size(); ++i) {
    Rcpp::Rcout << beta[i] << " ";
  }
  Rcpp::Rcout << "\n";
  
  Rcpp::Rcout << "\n Subgroup Shift Delta matrix (ncol "<< Delta.ncol()<<", nrow "<<Delta.nrow()<<"):\n";
  for (int i = 0; i < Delta.nrow(); ++i) {
    for (int j = 0; j < Delta.ncol(); ++j) {
      Rcpp::Rcout << Delta(i, j) << " ";
    }
    Rcpp::Rcout << "\n";
  }
  Rcpp::Rcout << "\n";
  
  
  for(int i=0 ; i<m; i++){
    Rcpp::Rcout <<"Partition Delta row"<<i+1<<" \n";
    Rcpp::Rcout <<"  Number of peaks q: "<<iq[i] <<"\n";
    Rcpp::Rcout << "  Weight vector w (size " << iw[i].size() << "): ";
    for (int j = 0; j < iw[i].size(); ++j) {
      Rcpp::Rcout << iw[i][j] << " ";
    }
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "  Allocation vector z : ";
    for (int j = 0; j < iz.ncol(); ++j) {
      Rcpp::Rcout << iz(i, j) << " ";
    }
    Rcpp::Rcout << "\n";
    
    Rcpp::Rcout << "  Omega (size " << iOmega[i].size() << "): ";
    for (int j = 0; j < iOmega[i].size(); ++j) {
      IntegerVector omega_j = iOmega[i][j];
      Rcpp::Rcout << "[";
      for (int k = 0; k < omega_j.size(); ++k) {
        Rcpp::Rcout << omega_j[k] << (k < omega_j.size() - 1 ? ", " : "");
      }
      Rcpp::Rcout << "] ";
    }
    Rcpp::Rcout << "\n";
    
    Rcpp::Rcout << "  Peak size m (size " << im[i].size() << "): ";
    for (int j = 0; j < im[i].size(); ++j) {
      Rcpp::Rcout << im[i][j] << " ";
    }
    Rcpp::Rcout << "\n";
    
    Rcpp::Rcout << "  Between cluster precision tau: " << itau[i] << "\n";
    
    Rcpp::Rcout << "  Cluster centers mu (size " << imu[i].size() << "): ";
    for (int j = 0; j < imu[i].size(); ++j) {
      Rcpp::Rcout << imu[i][j] << " ";
    }
    Rcpp::Rcout << "\n";
    
    Rcpp::Rcout << "  Within cluster variances sigma (size " << isigma[i].size() << "): ";
    for (int j = 0; j < isigma[i].size(); ++j) {
      Rcpp::Rcout << isigma[i][j] << " ";
    }
    Rcpp::Rcout << "\n";
    
    Rcpp::Rcout << "  Move type: " << iMove[i] << "\n";
    Rcpp::Rcout << "  Probability of acceptance: " << iP_Accept[i] << "\n";
    Rcpp::Rcout << "  Move accepted: " << iAccepted[i] << "\n";
    
    
    
  }
  
  Rcpp::Rcout << "\n";
  
};




void Model::update (const IntegerVector & Data_i, const IntegerVector& Data_k, const NumericVector& Data_Y,
                    NumericVector c_beta, NumericVector p_beta,
                    double a_cell, double b_cell, double a_between, double b_between, double a_within, double b_within) {
  
  for (int i=0; i<m;i++){
    
    //update w
    iw[i] = rdirichlet( 1 + im[i]);
    
    //update z
    NumericVector iz_k_prob_(iq[i]);     //length of prob = q[i]
    
    for(int k =0; k<s; k++){
      NumericVector sigma_vec = as<NumericVector>(isigma[i]);
      sigma_vec = pmax(sigma_vec, 1e-8);  
      
      NumericVector log_prob = log(iw[i]) - 0.5*log(sigma_vec)
        - 0.5* pow( Delta(i,k)-imu[i] , 2.0)/sigma_vec;
      
      //iz_k_prob_ = iw[i] * sqrt(1.0/isigma[i]) * exp( -0.5 * pow( Delta(i,k)-imu[i] , 2.0)/ isigma[i]); //calculate P(iz_k=t)
      //for(int t=0; t<iq[i]; t++){
      // iz_k_prob_[t] = iw[i][t] * R::dnorm(Delta(i,k), imu[i][t], sqrt(isigma[i][t]), 0);
      //}
      iz_k_prob_ = exp(log_prob - max(log_prob));
      iz(i,k) = sample(iq[i],1,TRUE,iz_k_prob_)[0];   
    }
    
    //Calculate Omega and m
    iOmega[i] = getOmega(iz(i,_),iq[i]);
    im[i] =getm(iz(i,_),iq[i]);
    
    //update sigma and mu for each cluster
    NumericVector Deltai_ssq_(iq[i]); //sum of squared deviation in each cluster
    NumericVector Deltai_sum_(iq[i]); //sum in each cluster
    for (int t=0; t<iq[i]; t++){
      IntegerVector iOmega_t = as<IntegerVector>(iOmega[i][t])-1; // list of members in cluster t, convert to 0-based
      for(int k : iOmega_t){
        //calculate sum and sum of sq in each cluster of Deltai
        Deltai_ssq_[t]+=pow(Delta(i,k)-imu[i][t] , 2.0);
        Deltai_sum_[t]+=Delta(i,k);
      }
      //update sigma, preventing underflow with std::max
      isigma[i][t] = std::max(1e-8,
                              1.0/R::rgamma( a_within+0.5*im[i][t], 1.0/(b_within + 0.5*Deltai_ssq_[t]) ) );
      //update mu
      double var_imu_t = 1.0/(std::max(1e-8,itau[i] + im[i][t]/isigma[i][t])) ;
      imu[i][t] = R::rnorm( Deltai_sum_[t]/isigma[i][t] * var_imu_t,
                            sqrt(var_imu_t) );
    }
    
    
    //update tau
    itau[i] = R::rgamma(a_between+0.5*iq[i], 1.0/( b_between + 0.5*sum(imu[i]*imu[i]) ) );
    
    
  }//for(int i=0; i<m;i++)
  
  
  
  
  
  //update cell precision cpre
  //summarize data by cell: cell count, cell mean, cell sum of square
  NumericMatrix r(m,s);
  NumericMatrix ybar(m,s);
  NumericMatrix ssq(m,s);
  for (int i=0;i<m;i++){
    for(int k=0;k<s;k++){
      double sumY=0; //container
      for(int l=0; l<Data_Y.size(); l++){
        if(Data_i[l]==i+1 && Data_k[l]==k+1){                 //filter cell ik
          r(i,k)++;                                           //count+1
          ssq(i,k)+= pow(Data_Y[l]-beta[i]-Delta(i,k), 2.0);  //sum of (Y-theta[ik])^2
          sumY+=Data_Y[l];
        }}
      ybar(i,k) = r(i,k) > 0 ? sumY / r(i,k) : 0.0;
      
    }}
  
  //update cell precision
  cpre = R::rgamma(a_cell+0.5*sum(r), 1.0/( b_cell+0.5*sum(ssq) ) );
  
  //recenter ybar for different components
  NumericMatrix ybar_delta(m,s);
  NumericMatrix ybar_beta(m,s);
  
  for (int i=0;i<m;i++){
    ybar_delta(i,_)=ybar(i,_)-beta(i);
    ybar_beta(i,_) = ybar(i,_)-Delta(i,_); 
  }
  
  
  //update Delta
  //calculate unconditional mean and variance for each Delta;
  NumericMatrix V(m,s);     //unconditional variance
  NumericMatrix Theta(m,s); //unconditional mean
  for (int i=0;i<m;i++){
    for(int k=0;k<s;k++){
      //z_k need to change to 0-based
      V(i,k) = 1.0 / ( 1.0/isigma[i][iz(i,k)-1] + r(i,k)*cpre );
      Theta(i,k) = V(i,k) *
        ( imu[i][iz(i,k)-1]/isigma[i][iz(i,k)-1] + r(i,k)*cpre*ybar_delta(i,k) );
    }
    Delta(i,_)=cond1samp(Theta(i,_),V(i,_));
  }
  
  
  //update beta
  
  double var_beta_i;
  double sumrybarbeta_i; //sumk(r_ik*ybeta_ik)
  
  for (int i=0; i<m; i++){
    var_beta_i = 1.0/( p_beta(i) + sum(r(i,_))*cpre );
    sumrybarbeta_i = sum(r(i,_)*ybar_beta(i,_));
    
    beta(i) = R::rnorm(var_beta_i*(c_beta(i)*p_beta(i)+cpre*sumrybarbeta_i),
         sqrt(var_beta_i));
  }
  
  
  
  
  
  
  
  
  
  
  
  
  //set beta parameters :original 2 2 1
  double au1=4.0;
  double au2=2.0;
  double au3=2.0;
  
  
  //rjMCMC step for Deltai
  for (int i=0;i<m;i++){
    //choose Move: 1merge 2split
    if(iq[i]==1){
      iMove[i]=2;
    }else if (iq[i]==s){
      iMove[i]=1;
    }else{
      iMove[i]=R::rbinom(1,0.60)+1;  //adjust the probability of split 
    }
    
    if (iMove[i]==2) {
      this->split(i,a_within,b_within,au1,au2,au3);
    } else if (iMove[i]==1){
      this->merge(i,a_within,b_within,au1,au2,au3);
    }
    
  }
  
  
  
}



void Model::merge(int i, double a_within, double b_within, 
                  double au1, double au2, double au3){
  //local variables according to the variable choosing move
  int q = iq[i]-1;            //q is new peak count for the candidate stage!!
  int t1 = sample(q,1)[0]-1;   //select a peak to merge; transfer to 0-based
  int t2 = iq[i]-1;          //t2 is the last peak (0-based)
  int mt1 = im[i][t1];               int mt2=im[i][t2]; 
  NumericVector Delta_i = this->Delta(i,_);
  
  //auxiliary variables 
  
  double u1=R::rbeta(au1,au1); 
  double u2=R::rbeta(au2,au2)*2.0-1.0;
  double u3=R::rbeta(au3,au3);
  
  double wt1=iw[i][t1];            double wt2=iw[i][t2];
  double mut1=imu[i][t1];          double mut2=imu[i][t2];
  double sigmat1=isigma[i][t1];    double sigmat2=isigma[i][t2];
  
  IntegerVector Omegat1 = as<IntegerVector>(iOmega[i][t1])-1;  //0-based
  IntegerVector Omegat2 = as<IntegerVector>(iOmega[i][t2])-1;  //0-based;
  
  std::vector<int> Omegat0;
  
  IntegerVector z_new=IntegerVector(iz(i,_));
  double tau=itau[i];
  
  
  //calculate candidate values
  double wt0 = wt1+wt2;
  double mut0 = (wt1*mut1+wt2*mut2)/wt0;
  double sigmat0 = wt1/wt0*(mut1*mut1+sigmat1)+wt2/wt0*(mut2*mut2+sigmat2)-mut0*mut0;
  
  
  //need to calculate Palloc
  double log_Palloc=0.0;
  
  
  for(int k : Omegat1){
    NumericVector LogProbt1t2(2); //the log probability of assigning to t1 and t2
    LogProbt1t2=getLogProbt1t2(wt1,wt2,mut1,mut2,sigmat1,sigmat2,Delta_i(k));
    
    //Palloc *= Probt1t2[0];
    log_Palloc += LogProbt1t2[0];
    Omegat0.push_back(k+1); //save 1-based
  }
  for(int k : Omegat2){
    NumericVector LogProbt1t2(2); //the log probability of assigning to t1 and t2
    LogProbt1t2=getLogProbt1t2(wt1,wt2,mut1,mut2,sigmat1,sigmat2,Delta_i(k));
    
    //Palloc *= Probt1t2[1];
    log_Palloc += LogProbt1t2[1];
    
    z_new[k] = t1+1; //move members in the last peak to t1 peak
    Omegat0.push_back(k+1); //save 1-based
  }
  
  //calculate sum of squares
  double sst1=0.0;
  double sst2=0.0;
  double sst1t2=0.0;
  //Omegat1 is 0-based;
  for(int k : Omegat1){
    sst1+=   (Delta_i(k)-mut1)*(Delta_i(k)-mut1)/sigmat1;
    sst1t2+= (Delta_i(k)-mut0)*(Delta_i(k)-mut0)/sigmat0;
  }
  //Omegat2 is 0-based;
  for(int k : Omegat2){
    sst2+=   (Delta_i(k)-mut2)*(Delta_i(k)-mut2)/sigmat2;
    sst1t2+= (Delta_i(k)-mut0)*(Delta_i(k)-mut0)/sigmat0;
  }
  
  
  // //calculate P_Accept
  // double Prob_Accept=
  //   sqrt( pow(sigmat1/sigmat0,mt1)*pow(sigmat2/sigmat0,mt2) )
  //   * exp( 0.5*(sst1+sst2-sst1t2) )
  //   * sqrt(2.0*M_PI/tau) * exp(0.5*tau*(mut1*mut1+mut2*mut2-mut0*mut0))
  //   * R::gammafn(a_within)/pow(b_within,a_within) * pow(sigmat1*sigmat2/sigmat0,a_within+1.0) * exp(b_within/sigmat1 + b_within/sigmat2 - b_within/sigmat0)
  //   * R::gammafn(q)/R::gammafn(q+1) * pow(wt0/wt1,mt1)*pow(wt0/wt2,mt2)
  //   * ( 0.5*R::dbeta(u1,2,2,false) * R::dbeta((u2+1.0)/2.0,2,2,false) * R::dbeta(u3,1,1,false) )
  //   * Palloc/PmoveRatio(q,s)
  //   * 1.0/(wt0*(1.0-u2*u2)) * pow(sigmat0/u1/(1.0-u1),-1.5);
  
  //calculate log_P_Accept
  double log_Prob_Accept=
    0.5*mt1*log(sigmat1/sigmat0)+0.5*mt2*log(sigmat2/sigmat0)
    + 0.5*(sst1+sst2-sst1t2) 
    + 0.5*log(2.0*M_PI/tau) + (0.5*tau*(mut1*mut1+mut2*mut2-mut0*mut0))
    +log(R::gammafn(a_within))-a_within*log(b_within) +(a_within+1.0)*log(sigmat1*sigmat2/sigmat0) + (b_within/sigmat1 + b_within/sigmat2 - b_within/sigmat0)
    +log(R::gammafn(q))-log(R::gammafn(q+1)) + mt1*log(wt0/wt1) +mt2*log(wt0/wt2)
    +log( 0.5*R::dbeta(u1,au1,au1,false) * R::dbeta((u2+1.0)/2.0,au2,au2,false) * R::dbeta(u3,au3,au3,false) )
    +log_Palloc-log(PmoveRatio(q,s)) 
    -log( wt0*(1.0 -u2*u2)*pow(sigmat0/u1/(1.0-u1),1.5)) ;
    
    double Prob_Accept = std::min(exp(log_Prob_Accept),1.0);
    int Accept = R::rbinom(1, Prob_Accept);
    this->iP_Accept[i]=Prob_Accept;
    this->iAccepted[i]=Accept;
    
    
    if (Accept){
      this->iq[i]-=1;
      this->iz(i,_) = z_new;
      this->iw[i][t1] = wt0;         this->iw[i].erase(t2);
      this->imu[i][t1] = mut0;       this->imu[i].erase(t2);
      this->isigma[i][t1] = sigmat0; this->isigma[i].erase(t2);
      this->iOmega[i][t1] = Omegat0; this->iOmega[i].erase(t2);
      this->im[i][t1] = mt1+mt2;     this->im[i].erase(t2);
    }
    return;
}

void Model::split(int i, double a_within, double b_within, 
                  double au1, double au2, double au3){
  //local variables according to the variable choosing move
  int q = iq[i];
  int t0 = sample(q,1)[0]-1;       //select a cluster to split; transfer to 0-based
  NumericVector Delta_i = this->Delta(i,_);
  
  //auxiliary variables
  double u1=R::rbeta(au1,au1);
  double u2=R::rbeta(au2,au2)*2.0 - 1.0;
  double u3=R::rbeta(au3,au3);
  
  //try truncate
  if(FALSE){
    double u2_lim = 0.8;   //abs of mu 
    double u3_low = 0.1, u3_high = 0.9;
    
    int max_attempts = 2;  //try _ times to truncate
    int tries = 0;
    while ( (std::fabs(u2) > u2_lim) && (tries < max_attempts) ) {
      u2 = R::rbeta(au2,au2)*2.0 - 1.0; 
      tries++;
    }
    tries = 0;
    while ( ((u3 < u3_low)||(u3>u3_high)) && (tries<max_attempts)) {
      u3 = R::rbeta(au3,au3);
      tries++;
    }
  }

  
  
  double wt0 = iw[i][t0];
  double mut0 = imu[i][t0];
  double sigmat0 = isigma[i][t0];

  IntegerVector Omegat0 = as<IntegerVector>(iOmega[i][t0])-1;  //0-based
  if (Omegat0.size() < 2) {
    // dont split singleton!
    this->iP_Accept[i]=-1;
    this->iAccepted[i]=-1;
    return;
  }
    IntegerVector z_new = IntegerVector(iz(i,_));
  double tau = itau[i];
  
  
  //calculate candidate values
  double wt1 =     u1*wt0;
  double wt2 = (1.0-u1)*wt0;
  double mut1 = mut0 - u2*sqrt(sigmat0*wt2/wt1);
  double mut2 = mut0 + u2*sqrt(sigmat0*wt1/wt2);
  double sigmat1 =     u3*(1.0-u2*u2)*sigmat0*wt0/wt1;
  double sigmat2 = (1.0-u3)*(1.0-u2*u2)*sigmat0*wt0/wt2;
  
  //decide how to re-allocate each subgrp in component t0; record m and omega
  double log_Palloc=0.0;
  
  int mt1=0;
  int mt2=0;
  std::vector<int> Omegat1;
  std::vector<int> Omegat2;
  for(int k : Omegat0){
    NumericVector LogProbt1t2(2); //the log probability of assigning to t1 and t2
    LogProbt1t2=getLogProbt1t2(wt1,wt2,mut1,mut2,sigmat1,sigmat2,Delta_i(k));
    NumericVector Probt1t2 = exp(LogProbt1t2);
    
    
    //for each subgrp in t0 reallocate to cluster t1 and t2
    IntegerVector t1t2= IntegerVector::create(t0 + 1, q + 1); //1-based
    
    //hard allocation
    double ratio = Probt1t2[0] / (Probt1t2[0] + Probt1t2[1]);
    if(ratio < 0) {   //hard allocation boundaries canbe changed
      z_new[k] = q+1;   // hard allocate to q
      log_Palloc += LogProbt1t2[1];
      mt2++;
      Omegat2.push_back(k+1);
    }
    else if(ratio > 1) {
      z_new[k] = t0+1; // hard allocate to t0
      log_Palloc += LogProbt1t2[0];
      mt1++;
      Omegat1.push_back(k+1);
    }
    else {
      // original sample : soft allocate
      z_new[k] = sample( t1t2, 1, false, Probt1t2 )[0]; // save 1-based in z
      if (z_new[k] == t0+1){
        mt1++;
        Omegat1.push_back(k+1); //record grp 1-based
        //Palloc *= Probt1t2[0];
        log_Palloc += LogProbt1t2[0];
      }
      if (z_new[k] ==  q+1){
        mt2++;
        Omegat2.push_back(k+1); //record grp 1-based
        //Palloc *= Probt1t2[1];
        log_Palloc += LogProbt1t2[1];
      }
    }
    
  }
  if (Omegat1.size() < 1||Omegat2.size()<1 ) {
    // dont split to empty clusters
    this->iP_Accept[i]=-2;
    this->iAccepted[i]=-2;
    return;
  }
  
  
  //calculate sum of squares
  double sst1=0.0;
  double sst2=0.0;
  double sst1t2=0.0;
  //Omegat1 is one based;
  for(int k : Omegat1){
    sst1+=   (Delta_i(k-1)-mut1)*(Delta_i(k-1)-mut1)/sigmat1;
    sst1t2+= (Delta_i(k-1)-mut0)*(Delta_i(k-1)-mut0)/sigmat0;
  }
  //Omegat2 is one based;
  for(int k : Omegat2){
    sst2+=   (Delta_i(k-1)-mut2)*(Delta_i(k-1)-mut2)/sigmat2;
    sst1t2+= (Delta_i(k-1)-mut0)*(Delta_i(k-1)-mut0)/sigmat0;
  }
  
  
  
  //calculate log(Paccept)
  double log_Prob_Accept=
    0.5*mt1*log(sigmat0/sigmat1)+0.5*mt2*log(sigmat0/sigmat2)
    -0.5*(sst1+sst2-sst1t2) 
    +0.5*log(tau/2.0/M_PI) -0.5*tau*(mut1*mut1+mut2*mut2-mut0*mut0)
    +a_within*log(b_within)-log(R::gammafn(a_within)) + (a_within+1.0)*log(sigmat0/sigmat1/sigmat2)-b_within/sigmat1-b_within/sigmat2+b_within/sigmat0
    +log(R::gammafn(q+1))-log(R::gammafn(q)) + mt1*log(wt1/wt0)+mt2*log(wt2/wt0)
    -log( ( 0.5 * R::dbeta(u1,au1,au1,false) * R::dbeta((u2+1.0)/2.0,au2,au2,false) * R::dbeta(u3,au3,au3,false) ))
    + log(PmoveRatio(q,s))-log_Palloc
    +log( wt0*(1.0 -u2*u2)*pow(sigmat0/u1/(1.0-u1),1.5));
    
    double Prob_Accept = std::min(exp(log_Prob_Accept),1.0);
    int Accept = R::rbinom(1, Prob_Accept);
    this->iP_Accept[i]=Prob_Accept;
    this->iAccepted[i]=Accept;
    if(FALSE){
      //log 
      std::ofstream logfile("split_debug_log.txt", std::ios::app); 
      
      if (logfile.is_open()) {
        logfile << "------ Split Proposal ------\n";
        logfile << "--u1: " <<u1<<", u2: "<<u2<<", u3: "<<u3<<"--\n";
        logfile << "Splitting cluster t0: " << t0 << "\n";
        logfile << "wt0: " << wt0 << "\n";
        logfile << "mut0: " << mut0 << "\n";
        logfile << "sigmat0: " << sigmat0 << "\n";
        logfile << "Omegat0: " << Omegat0 << "\n\n";
        logfile << "wt1: " << wt1 << ",    wt2: " << wt2 << "\n";
        logfile << "mut1: " << mut1 << ",    mut2: " << mut2 << "\n";
        logfile << "sigmat1: " << sigmat1 << ",    sigmat2: " << sigmat2 << "\n";
        logfile << "Omegat1: " ; for (int x=0; x<mt1;x++){logfile<<Omegat1[x]<<" ";}
        logfile<< ",    Omegat2: "; for (int x=0; x<mt2;x++){logfile<<Omegat2[x]<<" ";}
        logfile<<  "\n\n";
        logfile << "log_Paccept: " << Prob_Accept << "\n";
        logfile << "Accepted: " << Accept << "\n";
        logfile << "----------------------------\n\n";
      }
    }

    
    if (Accept){
      this->iq[i]+=1;
      this->iz(i,_) = z_new;
      this->iw[i][t0] = wt1;         this->iw[i].push_back(wt2);
      this->imu[i][t0] = mut1;       this->imu[i].push_back(mut2);
      this->isigma[i][t0] = sigmat1; this->isigma[i].push_back(sigmat2);
      this->iOmega[i][t0] = Omegat1; this->iOmega[i].push_back(Omegat2);
      this->im[i][t0] = mt1;         this->im[i].push_back(mt2);
    }
    return; 
    
}


//snapshot of all the parameters in the chain model vector
void collectModel(  const std::string& results_dir,
                    int dataid,  int chainid, const std::vector<Model>& Model_n,
                    int ngrp, int narm){
  
  int n = Model_n.size(); //a vector of models. 
  
  //mean and precision component
  NumericMatrix beta_n(n, narm);
  NumericMatrix Delta_n(n, narm*ngrp);
  NumericVector cpre_n(n);
  
  IntegerMatrix iMove_n(n,narm);
  NumericMatrix iP_Accept_n(n,narm);
  IntegerMatrix iAccepted_n(n,narm);
  
  //partitioning each arm
  IntegerMatrix iq_n(n, narm);       //iq1, ... iqm
  std::fill(iq_n.begin(), iq_n.end(), NA_REAL);
  NumericMatrix itau_n(n,narm);      //itau1,....itaum
  std::fill(itau_n.begin(), itau_n.end(), NA_REAL);
  
  IntegerMatrix iz_n(n,narm*ngrp);   //iz1,...,izm
  std::fill(iz_n.begin(), iz_n.end(), NA_REAL);
  NumericMatrix iw_n(n,narm*ngrp);  //iw1,...,iwm
  std::fill(iw_n.begin(), iw_n.end(), NA_REAL);
  NumericMatrix imu_n(n,narm*ngrp);
  std::fill(imu_n.begin(), imu_n.end(), NA_REAL);
  NumericMatrix isigma_n(n,narm*ngrp);
  std::fill(isigma_n.begin(), isigma_n.end(), NA_REAL);
  
  
  //RjMCMC monitoring
  for (int cyc=0; cyc<n ;cyc++){
    
    const Model& M=Model_n[cyc];
    
    //mean and precision
    beta_n(cyc,_) = (M.beta);
    Delta_n(cyc,_) = MatToVecByRow(M.Delta); 
    cpre_n(cyc) = (M.cpre);
    //partitioning
    iz_n(cyc,_)=MatToVecByRow(M.iz);
    iq_n(cyc,_)=(M.iq);
    itau_n(cyc,_)=(M.itau);
    
    iw_n(cyc,_) = VecOfVecToVec(ngrp, M.iw);
    imu_n(cyc,_) = VecOfVecToVec(ngrp, M.imu);
    isigma_n(cyc,_) = VecOfVecToVec(ngrp, M.isigma);
    
    iMove_n(cyc,_) = M.iMove;
    iP_Accept_n(cyc,_) = M.iP_Accept;
    iAccepted_n(cyc,_) = M.iAccepted;
  }
  
  
  write_csv(beta_n,results_dir,"beta",dataid,chainid);
  write_csv(Delta_n,results_dir,"Delta",dataid, chainid);
  write_csv(cpre_n,results_dir,"cpre",dataid, chainid);
  write_csv(iz_n,results_dir,"iz",dataid,chainid);
  
  write_csv(iq_n,results_dir,"iq",dataid, chainid);
  write_csv(itau_n,results_dir,"itau",dataid, chainid);
  write_csv(iw_n,results_dir,"iw",dataid, chainid);
  write_csv(imu_n,results_dir,"imu",dataid, chainid);
  write_csv(isigma_n,results_dir,"isigma",dataid, chainid);
  
  write_csv(iMove_n,results_dir,"iMoveType",dataid, chainid);
  write_csv(iP_Accept_n,results_dir,"iP_Accept",dataid, chainid);
  write_csv(iAccepted_n,results_dir,"iAccepted",dataid, chainid);
  
  
}

// [[Rcpp::export]]
void Bharp(
    IntegerVector Data_i, IntegerVector Data_k, NumericVector Data_Y,
    int ngrp, int narm,
    NumericVector c_beta, NumericVector p_beta,
    double a_cell, double b_cell, double a_between, double b_between, double a_within, double b_within,
    const std::string& results_dir,
    int nsamp, int thn = 1, int burn = 0, int dataid = 0, int nchain = 4,
    int max_attempts_per_chain = 5   
)  {
  Rcpp::Rcout<< "Dataset "<< dataid <<":\n";
  for (int chainid = 1; chainid <= nchain; chainid++) {
    
    bool chain_success = false;
    int attempt_count = 0;
    int step_count = -1;
    
    while (!chain_success && attempt_count < max_attempts_per_chain) {
      attempt_count++;
      Rcpp::Rcout <<"  Starting chain " << chainid 
                  << ", attempt " << attempt_count << "...\n";
      
      try {
        // 1) initialize
        std::vector<Model> ModelVec;
        Model Mcur(ngrp, narm,
                   c_beta, p_beta,
                   a_cell, b_cell, a_between, b_between, a_within, b_within);
        step_count = 0;
        // 2) burn-in
        for (int i = 0; i < burn; i++) {
          Mcur.update(Data_i, Data_k, Data_Y,
                      c_beta, p_beta,
                      a_cell, b_cell, a_between, b_between, a_within, b_within);
          step_count++;
        }
        
        ModelVec.push_back(Mcur);
        
        // 3) MCMC sample
        for (int i = 1; i < nsamp; i++) {
          for (int j = 0; j < thn; j++) {
            Mcur.update(Data_i, Data_k, Data_Y,
                        c_beta, p_beta,
                        a_cell, b_cell, a_between, b_between, a_within, b_within);
            step_count++;
          }
          ModelVec.push_back(Mcur);
        }
        
        // 4) save
        collectModel(results_dir, dataid, chainid, ModelVec, ngrp, narm);
        
        // 
        chain_success = true;
        Rcpp::Rcout << "  Completed successfully.\n";
        
      } catch (std::exception &ex) {
        
        Rcpp::Rcout << "  Failed on attempt " << attempt_count
                    << " step "<< step_count
                    << " with error:\n      [" << ex.what() << "]\n";
        if (attempt_count < max_attempts_per_chain) {
          Rcpp::Rcout << "  Retrying... \n";
        } else {
          Rcpp::Rcout << "  Failed after " << max_attempts_per_chain 
                      << " attempts. Moving to next chain.\n";
        }
      }
      
    } // end while(!chain_success && attempt_count < max_attempts_per_chain)
    
  } // end for chainid
  
  Rcpp::Rcout << "Dataset "<< dataid << ": All chains finished.\n";
}




