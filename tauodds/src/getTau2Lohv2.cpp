//Credit: Inspired and amended from https://github.com/HopkinsIDD/IDSpatialStats/commit/2782d6dcc9ee4be9855b5e468ce789425b81d49a by @gilesjohnr @jlessler
//Author: @t-pollington
//Date: 27 Sep 2019
#include <cstdlib>
#include <iostream>
#include <Rcpp.h>
#define NTYPE unsigned short //fine if N<=65535
using namespace std;
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector getTau2Lohv2(const NumericVector x, const NumericVector y, const NumericVector onset, const NumericVector r_low, const NumericVector r, const IntegerVector index, const IntegerVector importT1, const IntegerVector importT2){ //see how multiple vector arguments are parsed in here rather than the single posmat matrix
  NTYPE i,j,k;
  float dist2 = 0;
  float r2 = 0;
  float r2_low = 0;
  unsigned long num_cnt = 0, denom_cnt = 0; //counters for those filling conditions//
  unsigned short r_size = r.size();
  NTYPE N = x.size();
  int *inds = INTEGER(index);
  NumericVector tau2(r_size, NULL);
  float thetaInf = 0;
  bool bstrapconflict = 0;
  bool temprelated = 0;
  bool withindist = 0;
  unsigned short int T1 = importT1[0];
  unsigned short int T2 = importT2[0];
  //printf("T1 = %d\n",T1);
  //printf("T2 = %d\n",T2);
  //printf("N = %d\n",N);
  //calc thetaInf
  for (i=0;i<N;i++) { // i loop only goes through indices chosen
    for (j=0;j<N;j++) {
      //printf("inds[i] = %d, j = %d\n", inds[i], j);
      bstrapconflict = ((inds[i] - 1) == j); //do not compare someone with themself if bootstrapping*/
      temprelated = (abs(onset[j]-onset[inds[i]-1]) <= T2) && (abs(onset[j]-onset[inds[i]-1]) >= T1); //could add temporal restrictions here that the pair be within the start/end dates of the study
      num_cnt = num_cnt + (!bstrapconflict && temprelated);
      denom_cnt = denom_cnt + (!bstrapconflict && !temprelated);
      //printf("num_cnt = %d, denom_cnt = %d\n", num_cnt, denom_cnt);
    }
  }
  thetaInf = (float)num_cnt/denom_cnt;
  //printf("num_cnt = %d, denom_cnt = %d, thetaInf = %f\n", num_cnt, denom_cnt,thetaInf);
  
  //calc pi(r2_low,r2)
  for (k=0;k<r_size;k++) {
    num_cnt = 0;
    denom_cnt = 0;
    r2  = r[k]*r[k]; //transformation to squared distances to make sqrt() redundant in this line and the one below
    r2_low = r_low[k]*r_low[k];
    //printf("k = %d, r2 = %f, r2_low = %f\n", k, r2, r2_low);
    
    for (i=0;i<N;i++) {
      for (j=0;j<N;j++) { 
        //printf("inds[i] = %d, j = %d\n", inds[i], j);
        dist2 = (x[inds[i]-1] - x[j])*(x[inds[i]-1] - x[j]) + (y[inds[i]-1] - y[j])*(y[inds[i]-1] - y[j]); //calculate the distance
        withindist = ((dist2 >= r2_low) && (dist2 < r2));
        if (!withindist) continue;
        bstrapconflict = ((inds[i] - 1) == j); //do not compare someone with themself if bootstrapping
        temprelated = (abs(onset[j]-onset[inds[i]-1]) <= T2) && (abs(onset[j]-onset[inds[i]-1]) >= T1);
        num_cnt = num_cnt + (!bstrapconflict && temprelated);
        denom_cnt = denom_cnt + (!bstrapconflict && !temprelated);
        //printf("within distance and num_cnt = %d, denom_cnt = %d\n", num_cnt, denom_cnt);
      }
    }
    //printf("num_cnt = %d, denom_cnt = %d\n", num_cnt, denom_cnt);
    tau2[k] = (float)num_cnt/denom_cnt; // pi(r.min,r.max)
    tau2[k] = (float)tau2[k]/thetaInf;
    //Rcout << "The value of tau2 : " << tau2[k] << " for k = " << k << "\n";
  }
  return(tau2);
}