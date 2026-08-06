// -------------------------------------------------------------------------
// BHARP MCMC Sampler
//
// This file implements the MCMC and reversible-jump MCMC algorithms for
// the Bayesian Hierarchical Adjustable Random Partition (BHARP) model.
//
// The exported Bharp() function runs multiple MCMC chains and writes
// posterior samples to CSV files for subsequent analysis in R.
//
// Model notation:
//   I: number of treatment arms
//   K: number of subgroups
//   q: number of mixture components within an arm
//   z: subgroup-to-component allocation
//
// This file is compiled from R using Rcpp::sourceCpp().
// -------------------------------------------------------------------------



// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(RcppArmadillo)]]

#define ARMA_USE_CURRENT
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <fstream>
#include <memory>


using namespace Rcpp;
using namespace std;



// calculate Omega from z
List getOmega(IntegerVector z,int q) {
  
  List Omega(q);
  for(int t = 1; t <= q; t++) {
    std::vector<int> indices;
    for(int i = 0; i < z.size(); i++) {
      if(z[i] == t) {
        indices.push_back(i + 1);
      }
    }
    Omega[t-1] = indices;
  }
  return Omega;
}

// calculate m from z
IntegerVector getm(IntegerVector z,int q) {
  IntegerVector m(q, 0); 
  IntegerVector::iterator p;
  for( p =z.begin(); p < z.end(); p++) {
    if(*p >= 1 && *p <= q) {
      m[*p - 1]++; 
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
  for (int j = 0; j < size; j++) {
    w[j] /= sum_;
  }
  return(w);
}



//calculate P(q+1 merge)/P( q split)
double PmoveRatio( const int q, const int K){
  if (q==1 && q+1<K) {return 0.40/1.0;}
  else if (q+1==K && q!=1) {return 1.0/0.60;}
  else if (q+1==K && q==1) {return 1.0;}
  else {return 0.40/0.60;}
}

//calculate power prior for q
NumericVector power_prior_seq(int K) {
  NumericVector result(K);
  double sum = 0.0;
  for (int i = 0; i < K; ++i) {
    result[i] = pow(i+1,2);
    sum += result[i];
  }
  for (int i = 0; i < K; ++i) {
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
  NumericVector result( n_vectors*unit_length); 
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
                             double theta_ik){
  NumericVector LogProbt1t2(2); 

  double LogP1 = std::log(wt1) - 0.5*std::log(sigmat1)
    - 0.5* std::pow( theta_ik-mut1 , 2.0)/sigmat1;

  double LogP2 = std::log(wt2) - 0.5*std::log(sigmat2)
    - 0.5* std::pow( theta_ik-mut2 , 2.0)/sigmat2;
  double M = std::max(LogP1, LogP2);
  
  double sum_exp = std::exp(LogP1 - M) + std::exp(LogP2 - M);
  double log_sum_exp = M + std::log(sum_exp);
  
  LogProbt1t2[0] = LogP1-log_sum_exp;
  LogProbt1t2[1] = LogP2-log_sum_exp;
  return LogProbt1t2;
}





struct Hyper {
  double a_cell,  b_cell,  a_between,  b_between,  a_within,  b_within;
  // rjMCMC auxiliary 
  double au1 = 4.0, au2 = 2.0, au3 = 2.0;
};


struct Model{
  std::shared_ptr<const Hyper> hp;
  int K;                    
  int I;                    
  double cpre;              
  NumericMatrix theta;      
  
  IntegerVector iq;                 
  std::vector<NumericVector> iw;    
  IntegerMatrix iz;                  
  
  std::vector<List> iOmega;          
  std::vector<IntegerVector> im;     
  NumericVector itau;                
  std::vector<NumericVector> imu;    
  std::vector<NumericVector> isigma; 
  
  IntegerVector iMove;        
  NumericVector iP_Accept;   
  IntegerVector iAccepted;    
  
  
  void update( const IntegerVector & Data_i, const IntegerVector& Data_k, const NumericVector& Data_Y);
  void split(int i, double a_within, double b_within, 
             double au1, double au2, double au3);
  void merge(int i, double a_within, double b_within, 
             double au1, double au2, double au3);

  
  
  Model &operator=(const Model& Mprev){
    if (this!=&Mprev){
      this->hp = Mprev.hp;
      this->K = (Mprev.K);
      this->I = (Mprev.I);
      this->cpre = (Mprev.cpre);
      this->theta = Rcpp::clone(Mprev.theta);
      this->iq = Rcpp::clone(Mprev.iq);
      
      this->iw.resize(Mprev.iw.size());
      for (size_t i = 0; i < Mprev.iw.size(); ++i) {
        this->iw[i] = Rcpp::clone(Mprev.iw[i]);
      }
      
      this->iz = Rcpp::clone(Mprev.iz);
      this->iOmega.resize(Mprev.iOmega.size());
      for (size_t i = 0; i < Mprev.iOmega.size(); ++i) {
        this->iOmega[i] = Rcpp::clone(Mprev.iOmega[i]);
      }
      
      this->im.resize(Mprev.im.size());
      for (size_t i = 0; i < Mprev.im.size(); ++i) {
        this->im[i] = Rcpp::clone(Mprev.im[i]);
      }
      
      this->itau = Rcpp::clone(Mprev.itau);
      
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
  
  Model() = default;
  Model(int ngrp, int narm,
        std::shared_ptr<const Hyper> hyper) 
    : hp(std::move(hyper)), K(ngrp), I(narm){
    
    const Hyper& H = *hp;
    iq = IntegerVector(I);                        
    iw = std::vector<NumericVector>(I);           
    iz = IntegerMatrix(I, K);                     
    iOmega = std::vector<List>(I);                
    im = std::vector<IntegerVector>(I);           
    itau = NumericVector(I);                      
    imu = std::vector<NumericVector>(I);          
    isigma = std::vector<NumericVector>(I);       
    iMove      = IntegerVector(I, 9);
    iP_Accept  = NumericVector(I, -9.0);
    iAccepted  = IntegerVector(I, 9);
    cpre = R::rgamma( H.a_cell, 1.0/H.b_cell);            
    theta = NumericMatrix(I, K);                      
    
    for (int i=0; i<I;i++){
      iq[i] = sample(K,1,TRUE,power_prior_seq(K))[0];                                
      iw[i] = rdirichlet(IntegerVector(iq[i], 1));            
      iz(i,_) = sample(iq[i], K, TRUE, iw[i]);                
      iOmega[i] = getOmega(iz(i,_), iq[i]);                  
      im[i] = getm(iz(i,_), iq[i]);                           
      itau[i] = R::rgamma(H.a_between, 1.0 / H.b_between);                    
      imu[i] = rnorm(iq[i], 0, sqrt(1.0 / itau[i]));                      
      isigma[i] = 1.0 / Rcpp::rgamma(iq[i], H.a_within, 1.0 / H.b_within);    
      for (int k=0; k<K; k++){                                            
        int t = iz(i,k)-1;                                          
        theta(i,k) = R::rnorm(imu[i][t], sqrt(isigma[i][t]));      
      }
    }
  }

  Model(const Model& Mprev) {
    *this = Mprev;
  }
  
}; 



void Model::update (const IntegerVector & Data_i, const IntegerVector& Data_k, const NumericVector& Data_Y) {
  const Hyper& H = *hp;
  for (int i=0; i<I;i++){
    iw[i] = rdirichlet( 1 + im[i]);
    NumericVector iz_k_prob_(iq[i]);     
    for(int k =0; k<K; k++){
      NumericVector sigma_vec = as<NumericVector>(isigma[i]);
      sigma_vec = pmax(sigma_vec, 1e-8);  
      NumericVector log_prob = log(iw[i]) - 0.5*log(sigma_vec)
        - 0.5* pow( theta(i,k)-imu[i] , 2.0)/sigma_vec;
      iz_k_prob_ = exp(log_prob - max(log_prob));
      iz(i,k) = sample(iq[i],1,TRUE,iz_k_prob_)[0];   
    }
    iOmega[i] = getOmega(iz(i,_),iq[i]);
    im[i] =getm(iz(i,_),iq[i]);
    NumericVector thetai_ssq_(iq[i]); 
    NumericVector thetai_sum_(iq[i]); 
    for (int t=0; t<iq[i]; t++){
      IntegerVector iOmega_t = as<IntegerVector>(iOmega[i][t])-1; 
      for(int k : iOmega_t){
        thetai_ssq_[t]+=pow(theta(i,k)-imu[i][t] , 2.0);
        thetai_sum_[t]+=theta(i,k);
      }
      isigma[i][t] = std::max(
        1e-8,
        1.0/R::rgamma( H.a_within+0.5*im[i][t], 1.0/(H.b_within + 0.5*thetai_ssq_[t]) ) 
      );
      double var_imu_t = 1.0/(std::max(1e-8,itau[i] + im[i][t]/isigma[i][t])) ;
      imu[i][t] = R::rnorm( 
        thetai_sum_[t]/isigma[i][t] * var_imu_t,
        sqrt(var_imu_t) 
      );
    }
    itau[i] = R::rgamma(H.a_between+0.5*iq[i], 1.0/( H.b_between + 0.5*sum(imu[i]*imu[i]) ) );
  }
  NumericMatrix r(I,K);
  NumericMatrix ybar(I,K);
  NumericMatrix ssq(I,K);
  for (int i=0;i<I;i++){
    for(int k=0;k<K;k++){
      double sumY=0; 
      for(int l=0; l<Data_Y.size(); l++){
        if(Data_i[l]==i+1 && Data_k[l]==k+1){                 
          r(i,k)++;                                           
          ssq(i,k)+= pow(Data_Y[l]-theta(i,k), 2.0);  
          sumY+=Data_Y[l];
        }
      }
      ybar(i,k) = r(i,k) > 0 ? sumY / r(i,k) : 0.0;
    }
  }
  cpre = R::rgamma(H.a_cell+0.5*sum(r), 1.0/( H.b_cell+0.5*sum(ssq) ) );
  double V;     
  double M; 
  for (int i=0;i<I;i++){
    for(int k=0;k<K;k++){
      V = 1.0 / ( 1.0/isigma[i][iz(i,k)-1] + r(i,k)*cpre );
      M = V *
        ( imu[i][iz(i,k)-1]/isigma[i][iz(i,k)-1] + r(i,k)*cpre*ybar(i,k) );
      theta(i,k) = R::rnorm (M,sqrt(V));
    }
  }

  for (int i=0;i<I;i++){
    //choose Move: 1merge 2split
    if(iq[i]==1){
      iMove[i]=2;
    }else if (iq[i]==K){
      iMove[i]=1;
    }else{
      iMove[i]=R::rbinom(1,0.60)+1;  
    }
    if (iMove[i]==2) {
      this->split(i,H.a_within,H.b_within,H.au1,H.au2,H.au3);
    } else if (iMove[i]==1){
      this->merge(i,H.a_within,H.b_within,H.au1,H.au2,H.au3);
    }
  }
}



void Model::merge(int i, double a_within, double b_within, 
                  double au1, double au2, double au3){
  int q = iq[i]-1;            //q is new peak count for the candidate stage!!
  int t1 = sample(q,1)[0]-1;   //select a peak to merge; transfer to 0-based
  int t2 = iq[i]-1;          //t2 is the last peak (0-based)
  int mt1 = im[i][t1];               int mt2=im[i][t2]; 
  NumericVector theta_i = this->theta(i,_);
  
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
    LogProbt1t2=getLogProbt1t2(wt1,wt2,mut1,mut2,sigmat1,sigmat2,theta_i(k));
    log_Palloc += LogProbt1t2[0];
    Omegat0.push_back(k+1); //save 1-based
  }
  for(int k : Omegat2){
    NumericVector LogProbt1t2(2); //the log probability of assigning to t1 and t2
    LogProbt1t2=getLogProbt1t2(wt1,wt2,mut1,mut2,sigmat1,sigmat2,theta_i(k));
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
    sst1+=   (theta_i(k)-mut1)*(theta_i(k)-mut1)/sigmat1;
    sst1t2+= (theta_i(k)-mut0)*(theta_i(k)-mut0)/sigmat0;
  }
  //Omegat2 is 0-based;
  for(int k : Omegat2){
    sst2+=   (theta_i(k)-mut2)*(theta_i(k)-mut2)/sigmat2;
    sst1t2+= (theta_i(k)-mut0)*(theta_i(k)-mut0)/sigmat0;
  }
  
  double log_Prob_Accept=
    0.5*mt1*log(sigmat1/sigmat0)+0.5*mt2*log(sigmat2/sigmat0)
    + 0.5*(sst1+sst2-sst1t2) 
    + 0.5*log(2.0*M_PI/tau) + (0.5*tau*(mut1*mut1+mut2*mut2-mut0*mut0))
    +log(R::gammafn(a_within))-a_within*log(b_within) +(a_within+1.0)*log(sigmat1*sigmat2/sigmat0) + (b_within/sigmat1 + b_within/sigmat2 - b_within/sigmat0)
    +log(R::gammafn(q))-log(R::gammafn(q+1)) + mt1*log(wt0/wt1) +mt2*log(wt0/wt2)
    +log( 0.5*R::dbeta(u1,au1,au1,false) * R::dbeta((u2+1.0)/2.0,au2,au2,false) * R::dbeta(u3,au3,au3,false) )
    +log_Palloc-log(PmoveRatio(q,K)) 
    -log( wt0*(1.0 -u2*u2)*pow(sigmat0/u1/(1.0-u1),1.5)) ;
    
    double Prob_Accept = std::min(exp(log_Prob_Accept),1.0);
    int Accept = R::rbinom(1, Prob_Accept);
    this->iP_Accept[i]=Prob_Accept;
    this->iAccepted[i]=Accept;
    
    
    if (Accept){
      this->iq[i]-=1;
      this->iz(i,_) = z_new;
      this->iw[i][t1] = wt0;         this->iw[i].erase(this->iw[i].begin() + t2);
      this->imu[i][t1] = mut0;       this->imu[i].erase(this->imu[i].begin() + t2);
      this->isigma[i][t1] = sigmat0; this->isigma[i].erase(this->isigma[i].begin() + t2);
      this->iOmega[i][t1] = Omegat0; this->iOmega[i].erase(this->iOmega[i].begin() + t2);
      this->im[i][t1] = mt1+mt2;     this->im[i].erase(this->im[i].begin() + t2);
    }
    return;
}

void Model::split(int i, double a_within, double b_within, 
                  double au1, double au2, double au3){
  //local variables according to the variable choosing move
  int q = iq[i];
  int t0 = sample(q,1)[0]-1;       //select a cluster to split; transfer to 0-based
  NumericVector theta_i = this->theta(i,_);
  
  //auxiliary variables
  double u1=R::rbeta(au1,au1);
  double u2=R::rbeta(au2,au2)*2.0 - 1.0;
  double u3=R::rbeta(au3,au3);
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
    NumericVector LogProbt1t2(2); 
    LogProbt1t2=getLogProbt1t2(wt1,wt2,mut1,mut2,sigmat1,sigmat2,theta_i(k));
    NumericVector Probt1t2 = exp(LogProbt1t2);
    
    IntegerVector t1t2= IntegerVector::create(t0 + 1, q + 1); //1-based

    double ratio = Probt1t2[0] / (Probt1t2[0] + Probt1t2[1]);
    if(ratio < 0) {   
      z_new[k] = q+1;   
      log_Palloc += LogProbt1t2[1];
      mt2++;
      Omegat2.push_back(k+1);
    }
    else if(ratio > 1) {
      z_new[k] = t0+1; 
      log_Palloc += LogProbt1t2[0];
      mt1++;
      Omegat1.push_back(k+1);
    }
    else {
     
      z_new[k] = sample( t1t2, 1, false, Probt1t2 )[0]; // save 1-based in z
      if (z_new[k] == t0+1){
        mt1++;
        Omegat1.push_back(k+1); //record grp 1-based
        log_Palloc += LogProbt1t2[0];
      }
      if (z_new[k] ==  q+1){
        mt2++;
        Omegat2.push_back(k+1); //record grp 1-based
        log_Palloc += LogProbt1t2[1];
      }
    }
    
  }
  if (Omegat1.size() < 1||Omegat2.size()<1 ) {
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
    sst1+=   (theta_i(k-1)-mut1)*(theta_i(k-1)-mut1)/sigmat1;
    sst1t2+= (theta_i(k-1)-mut0)*(theta_i(k-1)-mut0)/sigmat0;
  }
  //Omegat2 is one based;
  for(int k : Omegat2){
    sst2+=   (theta_i(k-1)-mut2)*(theta_i(k-1)-mut2)/sigmat2;
    sst1t2+= (theta_i(k-1)-mut0)*(theta_i(k-1)-mut0)/sigmat0;
  }
  
  
  
  //calculate log(Paccept)
  double log_Prob_Accept=
    0.5*mt1*log(sigmat0/sigmat1)+0.5*mt2*log(sigmat0/sigmat2)
    -0.5*(sst1+sst2-sst1t2) 
    +0.5*log(tau/2.0/M_PI) -0.5*tau*(mut1*mut1+mut2*mut2-mut0*mut0)
    +a_within*log(b_within)-log(R::gammafn(a_within)) + (a_within+1.0)*log(sigmat0/sigmat1/sigmat2)-b_within/sigmat1-b_within/sigmat2+b_within/sigmat0
    +log(R::gammafn(q+1))-log(R::gammafn(q)) + mt1*log(wt1/wt0)+mt2*log(wt2/wt0)
    -log( ( 0.5 * R::dbeta(u1,au1,au1,false) * R::dbeta((u2+1.0)/2.0,au2,au2,false) * R::dbeta(u3,au3,au3,false) ))
    + log(PmoveRatio(q,K))-log_Palloc
    +log( wt0*(1.0 -u2*u2)*pow(sigmat0/u1/(1.0-u1),1.5));
    
    double Prob_Accept = std::min(exp(log_Prob_Accept),1.0);
    int Accept = R::rbinom(1, Prob_Accept);
    this->iP_Accept[i]=Prob_Accept;
    this->iAccepted[i]=Accept;
    
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
  
  int n = Model_n.size();

  NumericMatrix theta_n(n, narm*ngrp);
  NumericVector cpre_n(n);
  
  IntegerMatrix iMove_n(n,narm);
  NumericMatrix iP_Accept_n(n,narm);
  IntegerMatrix iAccepted_n(n,narm);
  
  //partitioning each arm
  IntegerMatrix iq_n(n, narm);       //iq1, ... iqm
  std::fill(iq_n.begin(), iq_n.end(), NA_INTEGER);
  NumericMatrix itau_n(n,narm);      //itau1,....itaum
  std::fill(itau_n.begin(), itau_n.end(), NA_REAL);
  
  IntegerMatrix iz_n(n,narm*ngrp);   //iz1,...,izm
  std::fill(iz_n.begin(), iz_n.end(), NA_INTEGER);
  NumericMatrix iw_n(n,narm*ngrp);  //iw1,...,iwm
  std::fill(iw_n.begin(), iw_n.end(), NA_REAL);
  NumericMatrix imu_n(n,narm*ngrp);
  std::fill(imu_n.begin(), imu_n.end(), NA_REAL);
  NumericMatrix isigma_n(n,narm*ngrp);
  std::fill(isigma_n.begin(), isigma_n.end(), NA_REAL);

  for (int cyc=0; cyc<n ;cyc++){
    const Model& M=Model_n[cyc];
    theta_n(cyc,_) = MatToVecByRow(M.theta); 
    cpre_n(cyc) = (M.cpre);
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
  

  write_csv(theta_n,results_dir,"theta",dataid, chainid);
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
    double a_cell, double b_cell, double a_between, double b_between, double a_within, double b_within,
    const std::string& results_dir,
    int nsamp, int thn = 1, int burn = 0, int dataid = 0, int nchain = 4,
    int max_attempts_per_chain = 5   
)  {
  Rcpp::Rcout<< "Dataset "<< dataid <<":\n";
  //hyperparameter
  auto hyper = std::make_shared<Hyper>();
  hyper->a_cell    = a_cell;
  hyper->b_cell    = b_cell;
  hyper->a_between = a_between;
  hyper->b_between = b_between;
  hyper->a_within  = a_within;
  hyper->b_within  = b_within;
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
        Model Mcur(ngrp, narm, hyper);
        step_count = 0;
        // 2) burn-in
        for (int i = 0; i < burn; i++) {
          Mcur.update(Data_i, Data_k, Data_Y);
          step_count++;
        }
        
        ModelVec.push_back(Mcur);
        
        // 3) MCMC sample
        for (int i = 1; i < nsamp; i++) {
          for (int j = 0; j < thn; j++) {
            Mcur.update(Data_i, Data_k, Data_Y);
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
      
    }
    
  } 
  
  Rcpp::Rcout << "Dataset "<< dataid << ": Chain processing finished.\n";
}




