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

NumericVector getTau2Loh(const NumericVector x, const NumericVector y, const NumericVector onset, const NumericVector r_low, const NumericVector r, const IntegerVector index, const IntegerVector importT1, const IntegerVector importT2){ //see how multiple vector arguments are parsed in here rather than the single posmat matrix
  NTYPE i,j,k;
  float dist2 = 0;
  float r2 = 0;
  float r2_low = 0;
  unsigned short r_size = r.size();
  NTYPE N = x.size();
  int *inds = INTEGER(index);
  NumericVector tau2(r_size, NULL);
  NumericVector thetaiInf(N, NULL);
  float thetai = 0;
  bool bstrapconflict = 0;
  bool temprelated = 0;
  bool withindist = 0;
  unsigned short int T1 = importT1[0];
  unsigned short int T2 = importT2[0];
  int mik1 = 0;
  int mik0 = 0;
  int mik1d = 0;
  int mik0d = 0;
  float taui = 0;
  float tausum = 0;
  int successcounter = 0; //times when taui is not nan and not Inf
  
  // printf("T1 = %d\n",T1);
  // printf("T2 = %d\n",T2);
  // printf("N = %d\n",N);
  
  //calc thetaInf
  for (i=0;i<N;i++) { // i loop only goes through indices chosen
    //thetaiInf calculation
    mik1 = 0;
    mik0 = 0;
    for (j=0;j<N;j++) {
      //printf("inds[i] = %d, j = %d\n", inds[i], j);
      bstrapconflict = ((inds[i] - 1) == j); //do not compare someone with themself if bootstrapping*/
      temprelated = (abs(onset[j]-onset[inds[i]-1]) <= T2) && (abs(onset[j]-onset[inds[i]-1]) >= T1); //could add temporal restrictions here that the pair be within the start/end dates of the study
      mik1 = mik1 + (!bstrapconflict && temprelated);
      mik0 = mik0 + (!bstrapconflict && !temprelated);
    }
    thetaiInf[i] = (float)mik1/mik0;
    //Rcout << "The value of thetaInf: " << thetaiInf[i] << " for i = " << i << "mik1=" << mik1 << "mik0=" << mik0 << "\n";
  }
  
  //calc thetai(r2_low,r2)
  for (k=0;k<r_size;k++) {
    r2 = r[k]*r[k]; //transformation to squared distances to make sqrt() redundant in this line and the one below
    r2_low = r_low[k]*r_low[k];
    //printf("k = %d, r2 = %f, r2_low = %f\n", k, r2, r2_low);
    tausum = 0;
    successcounter = 0;
    for (i=0;i<N;i++) {
      thetai = 0;
      taui = 0;
      mik1d = 0;
      mik0d = 0;
      for (j=0;j<N;j++) {
        //printf("inds[i] = %d, j = %d\n", inds[i], j);
        dist2 = (x[inds[i]-1] - x[j])*(x[inds[i]-1] - x[j]) + (y[inds[i]-1] - y[j])*(y[inds[i]-1] - y[j]); //calculate the distance
        withindist = ((dist2 >= r2_low) && (dist2 < r2));
        if (!withindist) continue;
        bstrapconflict = ((inds[i] - 1) == j); //do not compare someone with themself if bootstrapping
        temprelated = (abs(onset[j]-onset[inds[i]-1]) <= T2) && (abs(onset[j]-onset[inds[i]-1]) >= T1);
        mik1d = mik1d + (!bstrapconflict && temprelated);
        mik0d = mik0d + (!bstrapconflict && !temprelated);
        //printf("within distance and mik1d = %d, mik0d = %d, indsiminus1=%d, j=%d\n", mik1d, mik0d, (inds[i]-1), j);
      }
      thetai = (float)mik1d/mik0d;
      //printf("mik1d = %d, mik0d = %d, thetai = %f\n", mik1d, mik0d, thetai);
      taui = thetai/thetaiInf[i];
      //printf("thetai = %f, thetaiInf = %f, taui = %f\n", thetai, thetaiInf[i], taui);
      if (!isnan(taui) && !isinf(taui)){
        successcounter = successcounter + 1;
        tausum = tausum + taui;
      }
      if (isnan(taui) || isinf(taui)){
        //printf("Warning: nan or inf occurred at distance band k = %d\n", k);
      }
      //printf("tausum = %f\n", tausum);
    }
    tau2[k] = (float)tausum/successcounter;
    //Rcout << "The value of tau2 : " << tau2[k] << " for k = " << k << "\n";
  }
  return(tau2);
}