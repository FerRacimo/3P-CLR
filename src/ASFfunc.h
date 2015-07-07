//
//  ASFfunc.h
//  
//
//  Created by Fernando Racimo on 2/9/14.
//
//

#ifndef ____ASFfunc__
#define ____ASFfunc__

#include <iostream>
#include <vector>
using namespace std;


double neutpath(double x, double y, double driftxy);
double selpath(double x, double y, double driftxy, double selcoeffxy, double recdistxy, double Nexy);
double sel_normal_traj(double beta);
double neut_normal_traj(double beta);
double sistfreq_given_outfreq(double pA, double pB, double pC, double driftA, double driftB, double driftC, double selcoeff, double recdist, double mode, double branch);
double Binom(double nValue, double nValue2);
double samp_given_freq(double m, double n, double x);
int tri_integrand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
double sistsamp_given_outfreq(double ma, double na, double mb, double nb, double pc, double driftA, double driftB, double driftC, double selcoeff, double recdist, double mode, double branch);
double pseg_tri(double ma, double na, double mb, double nb, double pc, double driftA, double driftB, double driftC, double selcoeff, double recdist, double mode, double branch);
double chensamp_given_outfreq(double ma, double na, double pc, double drift, double selcoeff, double recdist, double mode);
double pseg_chen(double ma, double na, double pc, double drift, double selcoeff, double recdist, double mode);
double *getdrifts(vector<double> &mavec,vector<double> &navec,vector<double> &mbvec,vector<double> &nbvec,vector<double> &mcvec,vector<double> &ncvec,int totalsites);
vector<double> SupLogLikTri(vector<int> indices, vector<double> &posvec, vector<double> &mavec, vector<double> &navec, vector<double> &mbvec, vector<double> &nbvec, vector<double> &mcvec, vector<double> &ncvec, double driftA, double driftB, double driftC, double pos, vector<double> &sPossible, vector<double> &weightsvec, double mode, double branch);
vector<double> SupLogLikChen(vector<int> indices, vector<double> &posvec, vector<double> &mavec, vector<double> &navec, vector<double> &mcvec, vector<double> &ncvec, double drift, double pos, vector<double> &sPossible, vector<double> &weightsvec, double mode);
vector<int> SampleIndices(vector<int> Original, int numberSNPs);


#endif /* defined(____ASFfunc__) */
