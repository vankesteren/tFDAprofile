/*
 Performant functions for calculating tFDA scores
 Erik-Jan van Kesteren
 Thesis - UMC Utrecht Julius Centre
 License: https://opensource.org/licenses/MIT

 Version information:
 platform       x86_64-w64-mingw32          
 arch           x86_64                      
 os             mingw32                     
 system         x86_64, mingw32             
 status                                     
 major          3                           
 minor          3.2                         
 year           2016                        
 month          10                          
 day            31                          
 svn rev        71607                       
 language       R                           
 version.string R version 3.3.2 (2016-10-31)
 nickname       Sincere Pumpkin Patch       

 Package: Rcpp
 Title: Seamless R and C++ Integration
 Version: 0.12.9
 Date: 2017-01-14

 Package: RcppArmadillo
 Title: 'Rcpp' Integration for the 'Armadillo' Templated Linear Algebra Library
 Version: 0.7.700.0.0
 Date: 2017-02-07

 Package: RcppProgress
 Title: An Interruptible Progress Bar with OpenMP Support for C++ in R Packages
 Version: 0.3
 Date: 2016-12-14
*/

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppArmadillo.h>
#include <progress.hpp>

using namespace Rcpp;

// [[Rcpp::export]]
double tScore(arma::vec x1, arma::vec x2) {
  double mdif = arma::mean(x1)-arma::mean(x2);
  double sd1 = arma::var(x1)/x1.size();
  double sd2 = arma::var(x2)/x2.size();
  return mdif/sqrt(sd1+sd2);
}

// [[Rcpp::export]]
double projT(arma::mat X, arma::vec cl) {
  // X needs to be an n*2 matrix (columns are variables, rows are observations)
  // cl needs to be indicator class labels (0 & 1)
  
  // first, subset X by class labels
  arma::mat cl0 = X.rows(find(1-cl));
  arma::mat cl1 = X.rows(find(cl));
  
  // then, generate empirical covariance matrix and mean vector
  arma::mat cov0 = cov(cl0);
  arma::mat cov1 = cov(cl1);
  arma::mat mean0 = arma::mean(cl0,0);
  arma::mat mean1 = arma::mean(cl1,0);
  
  // calculate normal to project original data
  arma::mat normal = (cov0+cov1).i()*(mean1-mean0).t();
  
  
  // calculate tFDA score from projected t and max t
  double projt = fabs(tScore(cl0*normal, cl1*normal));
  
  return projt;
}

// [[Rcpp::export]]
double maxT(arma::mat X, arma::vec cl) {
  // X needs to be an n*2 matrix (columns are variables, rows are observations)
  // cl needs to be indicator class labels (0 & 1)
  
  // first, subset X by class labels
  arma::mat cl0 = X.rows(find(1-cl));
  arma::mat cl1 = X.rows(find(cl));
  
  double t0 = fabs(tScore(cl0.col(0), cl1.col(0)));
  double t1 = fabs(tScore(cl0.col(1), cl1.col(1)));
  double maxt = (t0<t1)?t1:t0; 
  
  return maxt;
}

// [[Rcpp::export]]
double tFDA(arma::mat X, arma::vec cl) {
  // X needs to be an n*2 matrix (columns are variables, rows are observations)
  // cl needs to be indicator class labels (0 & 1)
  
  // first, subset X by class labels
  arma::mat cl0 = X.rows(find(1-cl));
  arma::mat cl1 = X.rows(find(cl));
  
  // then, generate empirical covariance matrix and mean vector
  arma::mat cov0 = cov(cl0);
  arma::mat cov1 = cov(cl1);
  arma::mat mean0 = arma::mean(cl0,0);
  arma::mat mean1 = arma::mean(cl1,0);
  
  // calculate normal to project original data
  arma::mat normal = (cov0+cov1).i()*(mean1-mean0).t();

  
  // calculate tFDA score from projected t and max t
  double projt = fabs(tScore(cl0*normal, cl1*normal));
  double t0 = fabs(tScore(cl0.col(0), cl1.col(0)));
  double t1 = fabs(tScore(cl0.col(1), cl1.col(1)));
  double maxt = (t0<t1)?t1:t0; 
  
  return projt/maxt;
}

// [[Rcpp::export]]
arma::mat rotate(arma::mat X, arma::vec cl) {
  // X needs to be an n*2 matrix (columns are variables, rows are observations)
  // cl needs to be indicator class labels (0 & 1)
  
  // first, subset X by class labels
  arma::mat cl0 = X.rows(find(1-cl));
  arma::mat cl1 = X.rows(find(cl));
  
  // then, generate empirical covariance matrix and mean vector
  arma::mat cov0 = cov(cl0);
  arma::mat cov1 = cov(cl1);
  arma::mat mean0 = arma::mean(cl0,0);
  arma::mat mean1 = arma::mean(cl1,0);
  
  // calculate normal to project original data
  arma::mat normal = (cov0+cov1).i()*(mean1-mean0).t();
  
  return X*normal;
}

// [[Rcpp::export]]
List tFDAprofile(arma::mat X, arma::vec cl, double maxt_max = 1.96, 
               double projt_min = 1.96, bool display_progress = true) {
  // Function to calculate Q-prop and mean tFDA in fourth quadrant
  // X needs to be an n*p matrix (columns are variables, rows are observations)
  // cl needs to be indicator class labels (0 & 1)
  // maxt_max and projt_min define the search quadrant
  // display_progress shows a progress bar
  
  // initialise info
  int ncols = X.n_cols;
  int outlen = (ncols*(ncols+1)/2)-ncols;
  Progress p(outlen, display_progress);

  int i = 0;
  double t = 0;
  
  for (int a = 0; a < ncols; a++) {
    for (int b = 0; b < a; b++) {
      // concatenate column numbers as name
      std::ostringstream ossOut;
      ossOut << a+1 << ":" << b+1;
      std::string name = ossOut.str();
      
      // select columns to pass to tfda function
      arma::uvec ind;
      ind << a << b;
      arma::mat X_sel = X.cols(ind);
      
      // calculate maxt and test for quadrant
      double max_t = maxT(X_sel, cl);
      if (max_t>=maxt_max) {
        p.increment();
        continue;
      }
      
      // calculate projected t and test for quadrant
      double proj_t = projT(X_sel, cl);
      if (proj_t<=projt_min) {
        p.increment();
        
        continue;
      }
      
      // if tests passed increment indicator and add tFDA
      i++;
      t += proj_t/max_t;
      p.increment();
    };
  };
  
  double qprop = (double) i / (double) outlen;
  double tfda = t / (double) i;
  
  return List::create(Named("total") = outlen,
                      Named("qn") = i,
                      Named("qprop") = qprop,
                      Named("tfda") = tfda);
}


// [[Rcpp::export]]
List tFDAquadrant(arma::mat X, arma::vec cl, double maxt_max = 1.96, 
                 double projt_min = 1.96, int top = 100, 
                 bool display_progress = true) {
  // Function to output tFDA scores from fourth quadrant and their column index
  // X needs to be an n*p matrix (columns are variables, rows are observations)
  // cl needs to be indicator class labels (0 & 1)
  // maxt_max and projt_min define the search quadrant
  // top indicates how many topscores from this quadrant should be output
  // display_progress shows a progress bar
  
  // initialise info
  int ncols = X.n_cols;
  int outlen = (ncols*(ncols+1)/2)-ncols;
  Progress p(outlen, display_progress);
  
  int i = 0;
  NumericVector tfda(top);
  NumericVector projt(top);
  NumericVector maxt(top);
  StringVector names(top);
  double meansel = 0;
  
  for (int a = 0; a < ncols; a++) {
    for (int b = 0; b < a; b++) {
      // concatenate column numbers as name
      std::ostringstream ossOut;
      ossOut << a+1 << ":" << b+1;
      std::string name = ossOut.str();
      
      // select columns to pass to tfda function
      arma::uvec ind;
      ind << a << b;
      arma::mat X_sel = X.cols(ind);
      
      // calculate maxt and test for quadrant
      double max_t = maxT(X_sel, cl);
      if (max_t>=maxt_max) {
        p.increment();
        continue;
      }
      
      // calculate projected t and test for quadrant
      double proj_t = projT(X_sel, cl);
      if (proj_t<=projt_min) {
        p.increment();
        continue;
      }
      
      // if tests passed calculate tfdascore
      double score = proj_t/max_t;
      
      // add to top-output if larger than min of top
      if (score > min(tfda)){ 
        int pos = which_min(tfda);
        tfda(pos) = score;
        names(pos) = name;
        projt(pos) = proj_t;
        maxt(pos) = max_t;
        
        // update mean
        i++;
        meansel = meansel + (score-meansel)/i;
      };
      
      p.increment();
    };
  };
  
  tfda.attr("names") = names;
  projt.attr("names") = names;
  maxt.attr("names") = names;
  
  return List::create(Named("tfda") = tfda,
                      Named("projt") = projt,
                      Named("maxt") = maxt,
                      Named("Meansel") = meansel);
}

// [[Rcpp::export]]
List tFDAqprop(arma::mat X, arma::vec cl, double maxt_max = 1.96, 
                  double projt_min = 1.96, bool display_progress = true) {
  // Function to calculate proportion of pairs in fourth quadrant
  // X needs to be an n*p matrix (columns are variables, rows are observations)
  // cl needs to be indicator class labels (0 & 1)
  // maxt_max and projt_min define the search quadrant
  // display_progress shows a progress bar
  
  // initialise info
  int ncols = X.n_cols;
  int outlen = (ncols*(ncols+1)/2)-ncols;
  Progress p(outlen, display_progress);
  
  int i = 0;
  
  for (int a = 0; a < ncols; a++) {
    for (int b = 0; b < a; b++) {
      // concatenate column numbers as name
      std::ostringstream ossOut;
      ossOut << a+1 << ":" << b+1;
      std::string name = ossOut.str();
      
      // select columns to pass to tfda function
      arma::uvec ind;
      ind << a << b;
      arma::mat X_sel = X.cols(ind);
      
      // calculate maxt and test for quadrant
      double max_t = maxT(X_sel, cl);
      if (max_t>=maxt_max) {
        p.increment();
        continue;
      }
      
      // calculate projected t and test for quadrant
      double proj_t = projT(X_sel, cl);
      if (proj_t<=projt_min) {
        p.increment();
        continue;
      }
      
      // if tests passed increment indicator
      i++;
      p.increment();
    };
  };
  
  double qprop = (double) i / (double) outlen;
  
  return List::create(Named("total") = outlen,
                      Named("qn") = i,
                      Named("qprop") = qprop);
}

