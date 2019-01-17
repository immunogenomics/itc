#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//typedef arma::mat MATTYPE;
//typedef arma::vec VECTYPE;
//typedef arma::fmat MATTYPE;
//typedef arma::fvec VECTYPE;



// [[Rcpp::export]]
arma::mat exp_mean(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, int nrow, const arma::uvec& groups, const arma::uvec& group_sizes) {
    int ngroups = group_sizes.n_elem;
    arma::mat res = arma::zeros<arma::mat>(nrow, ngroups);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            res(i[j], groups[c]) += std::expm1(x[j]);
        }
    }
    
    for (int c = 0; c < ngroups; c++) {
        for (int r = 0; r < nrow; r++) {
            res(r, c) /= group_sizes[c];
        }
    }
        
    return(res);
}


// Updated: 01/16/2019
// Written by: Ilya Korsunsky, ilya.korsunsky@gmail.com
// Functions to more efficiently process raw single-cell RNA-seq data
// Used in utils_ilya.R
// (Functions adapted from Seurat package [Butler et al., 2018])

// [[Rcpp::export]]
arma::mat log_vmr(const arma::vec& x, const arma::vec& p, const arma::vec& i, 
                  int ncol, int nrow, const arma::mat& means,
                  const arma::uvec& groups, const arma::uvec& group_sizes) {
    
    int ngroups = group_sizes.n_elem;
    arma::mat res = arma::zeros<arma::mat>(nrow, ngroups);
    arma::mat nnzero = arma::zeros<arma::mat>(nrow, ngroups);
    double tmp;
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            // i[j] gives the row num
            tmp = std::expm1(x[j]) - means(i[j], groups(c));
            res(i[j], groups[c]) += tmp * tmp;
            nnzero(i[j], groups(c))++;
        }
    }
    
    for (int c = 0; c < ngroups; c++) {
        for (int r = 0; r < nrow; r++) {
            res(r, c) += (group_sizes[c] - nnzero(r, c)) * means(r, c) * means(r, c);
            res(r, c) /= (group_sizes[c] - 1);
        }
    }
    
    res = log(res / means);
    res.replace(arma::datum::nan, 0);
    
    return(res);
}

// [[Rcpp::export]]
arma::vec normalizeCLR_dgc(const arma::vec& x, const arma::vec& p, const arma::vec& i, int ncol, int nrow, int margin) {    
    arma::vec res = x;
    if (margin == 1) {
        // first compute scaling factors for each row
        arma::vec geo_mean = arma::zeros<arma::vec>(nrow);
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                // i[j] gives the row num
                geo_mean(i[j]) += std::log1p(x[j]);
            }
        }
        for (int i = 0; i < nrow; i++) {
//            geo_mean(i) = (geo_mean(i) / (1 + ncol));    
            geo_mean(i) = std::exp(geo_mean(i) / ncol);    
        }
        // then  scale data
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                res(j) = std::log1p(res(j) / geo_mean(i[j]));
            }
        }        
    } else {
        // first compute scaling factors for each column
        arma::vec geo_mean = arma::zeros<arma::vec>(ncol);
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                geo_mean(c) += std::log1p(x[j]);
            }
            geo_mean(c) = std::exp(geo_mean(c) / nrow);
        }
        
        // then  scale data
        for (int c = 0; c < ncol; c++) {
            for (int j = p[c]; j < p[c + 1]; j++) {
                res(j) = std::log1p(res(j) / geo_mean(c));
            }
        }        
        
    }
    
    return res;
}
