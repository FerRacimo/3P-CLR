//
//  ASFfunc.cpp
//  
//
//  Created by Fernando Racimo on 2/9/14.
//
//

#include "ASFfunc.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <string>
#include <cstring>
#include <cmath>
#include <sstream>
#include <array>
#include <iomanip>



// Using cubature package functions
#include "cubature.h"

using namespace std;

// Global definitions
#define PI 3.141592654





// FROM HUA CHEN
#define NR_END 1
#define FREE_ARG char*
void free_vector(double *v, long nl, long nh)
// free a double vector allocated with vector()
{
    free((FREE_ARG) (v+nl-NR_END));
}

// FROM HUA CHEN
double *allocvec(long nl, long nh)
{ double *v;
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if(!v) perror("allocation failure in vector");
    return v-nl+NR_END;
}

// FROM HUA CHEN
#define FUNC(x) ((*func)(x))
double trapzd(double (*func)(double), double a, double b, int n)
{
    double x,tnm,sum,del;
    static double s;
    int it,j;
    
    if (n == 1) {
        return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
    } else {
        for (it=1,j=1;j<n-1;j++) it <<= 1;
        tnm=it;
        del=(b-a)/tnm;
        x=a+0.5*del;
        for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
        s=0.5*(s+(b-a)*sum/tnm);
        return s;
    }
}
#undef FUNC

// FROM HUA CHEN
#define NRANSI
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
    int i,m,ns=1;
    double den,dif,dift,ho,hp,w;
    double *c,*d;
    
    dif=fabs(x-xa[1]);
    c=allocvec(1,n);
    d=allocvec(1,n);
    for (i=1;i<=n;i++) {
        if ( (dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1;m<n;m++) {
        for (i=1;i<=n-m;i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0) perror("Error in routine polint");
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }
    free_vector(d,1,n);
    free_vector(c,1,n);
}
#undef NRANSI

// FROM HUA CHEN
#define EPSQ 1.0e-2
#define JMAX 30
#define JMAXP (JMAX+1)
//#define K 5
#define K 5
double qromb(double (*func)(double), double a, double b)
{
    void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
    double trapzd(double (*func)(double), double a, double b, int n);
    
    double ss,dss;
    double s[JMAXP],h[JMAXP+1];
    int j;
    
    h[1]=1.0;
    for (j=1;j<=JMAX;j++) {
        s[j]=trapzd(func,a,b,j);
        if (j >= K) {
            polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
            if (fabs(dss) <= EPSQ*fabs(ss)) return ss;
        }
        h[j+1]=0.25*h[j];
    }
    perror("Too many steps in routine qromb");
    return 0.0;
}
#undef EPSQ
#undef JMAX
#undef JMAXP
#undef K






// Neutral Normal path
double neutpath(double x, double y, double driftxy){
    
    double variance;
    double density;
    
	variance = driftxy * y * (1.0-y);

    density = 1.0/sqrt(2.0*PI*variance)*exp(-pow((x-y),2)/(2.0*variance));

	return density;
}


// Selection Normal path
double selpath(double x, double y, double driftxy, double selcoeffxy, double recdistxy, double Nexy){
    
	double g = 1.0 - pow((1/(2*Nexy)),(recdistxy/selcoeffxy));
	double variance = driftxy * y * (1.0-y);
	double initial = 0.0;
    double predensity;
    double density;
    
    if (x < 1.0 && x > (1.0-g)){
        predensity = initial + (1.0/sqrt(2.0 * PI * variance)) * ((x+g-1.0)/pow(g,2.0))*exp(-pow(x+g-1.0-g*y,2.0)/(2.0*pow(g,2.0)*variance));
    }
    else {
        predensity = initial;
    }
    
    if (x > 0.0 && x < g) {
        density = predensity + (1.0/sqrt(2.0 * PI * variance)) * ((g-x)/pow(g,2.0))*exp(-pow(x-g*y,2.0)/(2.0*pow(g,2.0)*variance));
    
    }
    else {
        density = predensity;
    }
    
    return density;
    
}



// Integrand for sistfreq_given_outfreq under neutrality
double neut_normal_traj(double beta){

    extern double p1, p2, p3, drift1, drift2, drift3;
    double density;
    
    density = neutpath(p1,beta,drift1) * neutpath(p2,beta,drift2) * neutpath(beta,p3,drift3);
    
    return density;
    
}


/*
int neut_normal_traj(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    vector<double> datapass = *((vector<double> *) fdata);
    
    double pA0 = datapass[0]; double pB0 = datapass[1]; double pC0 = datapass[2];
    double drift01 = datapass[3]; double drift02 = datapass[4]; double drift03 = datapass[5];
    
    fval[0] = neutpath(pA0,x[0],drift01) * neutpath(pB0,x[0],drift02) * neutpath(x[0],pC0,drift03);
    
    return 0;

}
*/





// Integrand for sistfreq_given_outfreq under selection

double sel_normal_traj(double beta){
    
    extern double p1, p2, p3, drift1, drift2, drift3, selcoeff123, recdist123, branch123, Ne;
    double density;
    
    if (branch123 == 0){
        density = neutpath(p1,beta,drift1) * neutpath(p2,beta,drift2) * selpath(beta,p3,drift3,selcoeff123,recdist123, Ne);
    }
    else if (branch123 == 1){
        density = selpath(p1,beta,drift1,selcoeff123,recdist123,Ne) * neutpath(p2,beta,drift2) * neutpath(beta,p3,drift3);
    }
    else if (branch123 == 2){
        density = neutpath(p1,beta,drift1) * selpath(p2,beta,drift2,selcoeff123,recdist123,Ne) * neutpath(beta,p3,drift3);
    }
    else{
        density = 0;
    }
    return density;
    
}

/*
int sel_normal_traj(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    vector<double> datapass = *((vector<double> *) fdata);
    
    double pA0 = datapass[0]; double pB0 = datapass[1]; double pC0 = datapass[2];
    double drift01 = datapass[3]; double drift02 = datapass[4]; double drift03 = datapass[5];
    double selcoeff0 = datapass[6]; double recdist0 = datapass[7];
    
    extern double Ne;
    
    fval[0] = neutpath(pA0,x[0],drift01) * neutpath(pB0,x[0],drift02) * selpath(x[0],pC0,drift03,selcoeff0,recdist0,Ne);
    
    return 0;
    
}
*/



// Normal sister-group frequencies given outgrup frequency

double sistfreq_given_outfreq(double pA, double pB, double pC, double driftA, double driftB, double driftC, double selcoeff, double recdist, double mode, double branch){

    extern double p1, p2, p3, drift1, drift2, drift3, selcoeff123, recdist123, branch123;
    double density;

    p1 = pA; p2 = pB; p3 = pC;
    drift1 = driftA; drift2 = driftB; drift3 = driftC;
    selcoeff123 = selcoeff; recdist123 = recdist;
    branch123 = branch;
    
    // neutrality
    if (mode == 0.0){
        density = qromb(neut_normal_traj,0.001,0.999);
        return density;
    }
    
    // selection
    else if (mode == 1.0){
        density = qromb(sel_normal_traj,0.001,0.999);
        return density;
    }
    
    else {
        return 0;
    }

}


/*
double sistfreq_given_outfreq(double pA, double pB, double pC, double driftA, double driftB, double driftC, double selcoeff, double recdist, double mode){

    vector<double> datapass(8);
    datapass = {pA,pB,pC,driftA,driftB,driftC,selcoeff,recdist};

    // neutrality
    if (mode == 0.0){
        double xmin[2] = {0.0001}, xmax[2] = {0.9999}, val, err;
        hcubature(1,neut_normal_traj, &datapass, 1, xmin, xmax, 0, 0, 1e-3, ERROR_INDIVIDUAL, &val, &err);
        //cout << val << " " << err << endl;
        //cout << pA << " " << pB << " " << pC << endl;
        return val;
    }
    
    // selection
    else {
        double xmin[2] = {0.0001}, xmax[2] = {0.9999}, val, err;
        pcubature(1,sel_normal_traj, &datapass, 1, xmin, xmax, 0, 0, 1e-3, ERROR_INDIVIDUAL, &val, &err);
        //cout << val << " " << err << endl;
        //cout << pA << " " << " " << pC << endl;
        //cout << "HERE" << endl;
        return val;
    }

}
*/

    
    

// Obtained from Rosetta Code
double Factorial(double nValue)
{
    double result = nValue;
    double result_next;
    double pc = nValue;
    do
    {
        result_next = result*(pc-1);
        result = result_next;
        pc--;
    }while(pc>2);
    nValue = result;
    return nValue;
}

// Obtained from Rosetta Code
double Binom(double nValue, double nValue2)
{
    double result;
    if(nValue2 == 1)return nValue;
    if(nValue2 == (nValue-1))return nValue;
    if(nValue2 == 0)return 1;
    if(nValue2 == nValue)return 1;
    result = (Factorial(nValue))/(Factorial(nValue2)*Factorial((nValue - nValue2)));
    nValue2 = result;
    return nValue2;
}



// Sample counts given frequency
double samp_given_freq(double m, double n, double x){
    
    double result;
    result = Binom(n,m) * pow(x,m) * pow((1.0 - x),(n-m));
    if( result >= 0) return result;
    else return 0;

}




// Integrand for sistsamp_given_outfreq
int tri_integrand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    vector<double> datapass = *((vector<double> *) fdata); // we can pass σ via fdata argument

    double ma0 = datapass[0]; double na0 = datapass[1]; double mb0 = datapass[2]; double nb0 = datapass[3]; double pc0 = datapass[4];
    double drift01 = datapass[5]; double drift02 = datapass[6]; double drift03 = datapass[7];
    double selcoeff0 = datapass[8]; double recdist0 = datapass[9]; double mode0 = datapass[10];
    double branch0 = datapass[11];
    
    fval[0] = samp_given_freq(ma0, na0, x[0]) * samp_given_freq(mb0, nb0, x[1]) * sistfreq_given_outfreq(x[0], x[1], pc0, drift01, drift02, drift03, selcoeff0, recdist0, mode0,branch0);
    
    return 0; // success
}



// Normal sister-group counts given outgroup frequency
double sistsamp_given_outfreq(double ma, double na, double mb, double nb, double pc, double driftA, double driftB, double driftC, double selcoeff, double recdist, double mode, double branch){

    vector<double> datapass(11);
    datapass = {ma,na,mb,nb,pc,driftA,driftB,driftC,selcoeff,recdist,mode,branch};
    
    double xmin[2] = {0.001,0.001}, xmax[2] = {0.999,0.999}, val, err;
    hcubature(1, tri_integrand, &datapass, 2, xmin, xmax, 0, 0, 1e-3, ERROR_INDIVIDUAL, &val, &err);
    
    return val;
    
}


// Integrand for pseg_tri
int pseg_tri_integrand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    vector<double> datapass = *((vector<double> *) fdata); // we can pass σ via fdata argument
    
    double pc0 = datapass[0];
    double drift01 = datapass[1]; double drift02 = datapass[2]; double drift03 = datapass[3];
    double selcoeff0 = datapass[4]; double recdist0 = datapass[5]; double mode0 = datapass[6];
    double branch0 = datapass[7];
    
    fval[0] = sistfreq_given_outfreq(x[0], x[1], pc0, drift01, drift02, drift03, selcoeff0, recdist0, mode0,branch0);
    
    return 0; // success
}


// Probabilty of segregation for tri-sample model
double pseg_tri(double ma, double na, double mb, double nb, double pc, double driftA, double driftB, double driftC, double selcoeff, double recdist, double mode, double branch){
    
    vector<double> datapass(7);
    datapass = {pc,driftA,driftB,driftC,selcoeff,recdist,mode,branch};
    
    double xmin[2] = {0.001,0.001}, xmax[2] = {0.999,0.999}, val, err;
    hcubature(1, pseg_tri_integrand, &datapass, 2, xmin, xmax, 0, 0, 1e-3, ERROR_INDIVIDUAL, &val, &err);
    
    return val;

}


// Integrand for chensamp_given_outfreq
int chen_integrand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    vector<double> datapass = *((vector<double> *) fdata);
    double ma0 = datapass[0]; double na0 = datapass[1]; double pc0 = datapass[2];
    double drift01 = datapass[3];
    double selcoeff0 = datapass[4]; double recdist0 = datapass[5]; double mode0 = datapass[6];
    extern double Ne;
    
    if(mode0 == 0){
        fval[0] = samp_given_freq(ma0,na0,x[0])*neutpath(x[0],pc0,drift01);
    }
    else{
        fval[0] = samp_given_freq(ma0,na0,x[0])*selpath(x[0],pc0,drift01,selcoeff0,recdist0,Ne);
    };


    return 0;
}

// Normal Chen counts given outgroup frequency
double chensamp_given_outfreq(double ma, double na, double pc, double drift, double selcoeff, double recdist, double mode){

    vector<double> datapass(7);
    datapass = {ma,na,pc,drift,selcoeff,recdist,mode};
    
    double xmin[1] = {0.001}, xmax[1] = {0.999}, val, err;
    hcubature(1, chen_integrand, &datapass, 1, xmin, xmax, 0, 0, 1e-3, ERROR_INDIVIDUAL, &val, &err);
    
    return val;

}

// Integrand for pseg_chen
int pseg_chen_integrand(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval){
    vector<double> datapass = *((vector<double> *) fdata);
    double pc0 = datapass[0];
    double drift01 = datapass[1];
    double selcoeff0 = datapass[2]; double recdist0 = datapass[3]; double mode0 = datapass[4];
    
    extern double Ne;
    
    if(mode0 == 0){
        fval[0] = neutpath(x[0],pc0,drift01);
    }
    else{
        fval[0] = selpath(x[0],pc0,drift01,selcoeff0,recdist0,Ne);
    };
    
    return 0; // success

    
};

// Probability of segregation under two-pop model
double pseg_chen(double ma, double na, double pc, double drift, double selcoeff, double recdist, double mode){

    vector<double> datapass(5);
    datapass = {pc,drift,selcoeff,recdist,mode};

    double xmin[1] = {0.001}, xmax[1] = {0.999}, val, err;
    hcubature(1, pseg_chen_integrand, &datapass, 1, xmin, xmax, 0, 0, 1e-3, ERROR_INDIVIDUAL, &val, &err);

    return val;
};




// Get drifts
double *getdrifts(vector<double> &mavec,vector<double> &navec,vector<double> &mbvec,vector<double> &nbvec,vector<double> &mcvec,vector<double> &ncvec, int totalsites){
    
    extern double drifts[3];
    double totalSNPs = totalsites;
    double cuma = 0;
    double cumb = 0;
    double cumc = 0;
    double drifta;
    double driftb;
    double driftc;
    
    int i;
    for (i=0; i<totalsites; i++) {

        double freqa = mavec[i]/navec[i];
        double freqb = mbvec[i]/nbvec[i];
        double freqc = mcvec[i]/ncvec[i];
        
        // filter for SNPs larger than 5% and smaller than 95%
        if (freqa > 0.05 && freqa < 0.95 && freqb > 0.05 && freqb < 0.95 && freqc > 0.05 && freqc < 0.95){
            cuma = cuma + (freqa - freqc)*(freqa - freqb) / (freqa*(1 - freqa));
            cumb = cumb + (freqb - freqa)*(freqb - freqc) / (freqb*(1 - freqb));
            cumc = cumc + (freqc - freqa)*(freqc - freqb) / (freqc*(1 - freqc));
        }
        // don't include non-qualified SNPs in the total SNP count
        else {totalSNPs = totalSNPs - 1;}
        
    }
    drifta = cuma / totalSNPs;
    driftb = cumb / totalSNPs;
    driftc = cumc / totalSNPs;
    
    drifts[0] = drifta;
    drifts[1] = driftb;
    drifts[2] = driftc;
    
    return 0;
    
}



// Window computation for neutral tri-model
vector<double> SupLogLikTri(vector<int> indices, vector<double> &posvec, vector<double> &mavec, vector<double> &navec, vector<double> &mbvec, vector<double> &nbvec, vector<double> &mcvec, vector<double> &ncvec, double driftA, double driftB, double driftC, double pos, vector<double> &sPossible, vector<double> &weightsvec, double mode, double branch){
    
    double mindist = 0.000001;
    double pos0 = pos+0.00001;
    //double maxLogLik = -100000;
    double FinalLogLik = 0;
    double BestS = 0;
    double CurrLik = 0;
    
    if (mode == 0){
        for (vector<int>::iterator it = indices.begin() ; it != indices.end(); ++it){
            double freqc = mcvec[*it]/ncvec[*it];

            double temp1 = sistsamp_given_outfreq(mavec[*it],navec[*it],mbvec[*it],nbvec[*it],freqc,driftA,driftB,driftC,0,0,0,branch);
            if(temp1 == 0 || temp1 != temp1 || fabs(temp1) > 100000000000) {  temp1 = 0.000000000000001; }

            double temp2 = pseg_tri(mavec[*it],navec[*it],mbvec[*it],nbvec[*it],freqc,driftA,driftB,driftC,0,0,0,branch);
            if(temp2 == 0 || temp2 != temp2 || fabs(temp2) > 100000000000) {  temp2 = 0.0000001; }
            
            double temp = log(temp1) - log(temp2);

            double loglik = weightsvec[*it]*temp;
            if(loglik != loglik){
                cout << "1 : Nan is here" << endl;
            }
            else if(loglik < -1000000000000){
                cout << "-Inf is here" << endl;
                cout << temp1 << " " << temp2 << " " << " " << temp << " " << loglik << endl;
            }
            else{
                FinalLogLik += loglik;
            }
        }
    }
    else {
        double BestLik = -10000000;
        for (vector<double>::iterator scurr = sPossible.begin() ; scurr != sPossible.end(); ++scurr){
            CurrLik = 0;
            
            int t = 0;
            
            for (vector<int>::iterator it = indices.begin() ; it != indices.end(); ++it){
                //cout << *scurr << endl;
                // find distance to central SNP
                double distrec = abs(posvec[*it]-pos0);
                //cout << distrec << endl;
                // round number
                distrec = floor(distrec*10000000000+0.5)/10000000000;
                distrec = max(distrec,mindist);
                double freqc = mcvec[*it]/ncvec[*it];

                double temp1 = sistsamp_given_outfreq(mavec[*it],navec[*it],mbvec[*it],nbvec[*it],freqc,driftA,driftB,driftC,*scurr,distrec,1,branch);
                if(temp1 == 0 || temp1 != temp1 || fabs(temp1) > 100000000000) {
                    //cout << mavec[*it] << " " << mbvec[*it] << " " << freqc << endl;
                    //cout << setprecision(20) << " " << *scurr << " " << distrec << " " << temp1 << endl;
                    temp1 = 0.0000000000001; t=t+1;
                    //cout << " " << endl;
                }
                //cout << temp1 << endl;

                double temp2 = pseg_tri(mavec[*it],navec[*it],mbvec[*it],nbvec[*it],freqc,driftA,driftB,driftC,*scurr,distrec,1,branch);
                if(temp2 == 0 || temp2 != temp2 || fabs(temp2) > 100000000000) { temp2 = 0.0000001; }
                //cout << temp2 << endl;

                double temp = log(temp1) - log(temp2);
                //cout << temp << endl;
                double loglik = weightsvec[*it]*temp;
                //cout << mavec[*it] << " " << navec[*it] << " " << mbvec[*it] << " " << nbvec[*it] << " " << mcvec[*it] << " " << ncvec[*it] << endl;
                //cout << loglik << endl;
                //cout << " " << endl;
                if(loglik != loglik){
                    cout << "2 : Nan is here" << endl;
                    cout << temp1 << " " << temp2 << " " << " " << temp << " " << loglik << endl;
                    cout << mavec[*it] << ", " << navec[*it] << ", " << mbvec[*it] << ", " << nbvec[*it] << ", " << freqc << ", " << driftA << ", " << driftB << ", " << driftC << ", " << *scurr << ", " << distrec << ", " << 1 << endl;
                }
                else if(loglik < -1000000000000){
                    cout << "-Inf is here" << endl;
                    cout << temp1 << " " << temp2 << " " << " " << temp << " " << loglik << endl;
                }
                else {
                    CurrLik += loglik;
                }
            }
            
            //cout << *scurr << " " << t << endl;
            
            if( CurrLik > BestLik){
                BestS = *scurr;
                BestLik = CurrLik + 0;
            }
        }
        FinalLogLik = BestLik;
    };
    
    vector<double> result = {FinalLogLik,BestS};
    return result;
    
}





// Window computation for neutral tri-model
vector<double> SupLogLikChen(vector<int> indices, vector<double> &posvec, vector<double> &mavec, vector<double> &navec, vector<double> &mcvec, vector<double> &ncvec, double drift, double pos, vector<double> &sPossible, vector<double> &weightsvec, double mode){
    
    double mindist = 0.000001;
    double pos0 = pos+0.00001;
    //double maxLogLik = -100000;
    double FinalLogLik = 0;
    double BestS = 0;
    double CurrLik = 0;
    
    if (mode == 0){
        for (vector<int>::iterator it = indices.begin() ; it != indices.end(); ++it){
            double freqc = mcvec[*it]/ncvec[*it];
            double temp1 = chensamp_given_outfreq(mavec[*it],navec[*it],freqc,drift,0,0,0);
            if(temp1 == 0 || temp1 != temp1 || fabs(temp1) > 100000000000) { temp1 = 0.0000000000001; }
            double temp2 = pseg_chen(mavec[*it],navec[*it],freqc,drift,0,0,0);
            if(temp2 == 0 || temp2 != temp2 || fabs(temp2) > 100000000000) { temp2 = 0.0000001; }
            if(temp2 != temp2) { temp2 = 0.000001; }
            double temp = log(temp1) - log(temp2);
            double loglik = weightsvec[*it]*temp;
            if(loglik != loglik){
                cout << "3 : Nan is here" << endl;
                cout << temp1 << " " << temp2 << " " << " " << temp << " " << loglik << endl;
            }
            else if(fabs(loglik) > 1000000000000){
                cout << "Inf is here" << endl;
                cout << temp1 << " " << temp2 << " " << " " << temp << " " << loglik << endl;
            }
            else{
                FinalLogLik += loglik;
            }
        }
    }
    else {
        double BestLik = -10000000;
        for (vector<double>::iterator scurr = sPossible.begin() ; scurr != sPossible.end(); ++scurr){
            CurrLik = 0;
            
            for (vector<int>::iterator it = indices.begin() ; it != indices.end(); ++it){
                // find distance to central SNP
                double distrec = abs(posvec[*it]-pos0);
                // round number
                distrec = floor(distrec*10000000000+0.5)/10000000000;
                distrec = max(distrec,mindist);
                double freqc = mcvec[*it]/ncvec[*it];
  
                double temp1 = chensamp_given_outfreq(mavec[*it],navec[*it],freqc,drift,*scurr,distrec,1);
                if(temp1 == 0 || temp1 != temp1 || fabs(temp1) > 100000000000) { temp1 = 0.0000000000001;}

                double temp2 = pseg_chen(mavec[*it],navec[*it],freqc,drift,*scurr,distrec,1);
                if(temp2 == 0 || temp2 != temp2 || fabs(temp2) > 100000000000) { temp2 = 0.0000001; }

                double temp = log(temp1) - log(temp2);
                double loglik = weightsvec[*it]*temp;
                if(loglik != loglik){
                    cout << "4 : Nan is here" << endl;
                    cout << temp1 << " " << temp2 << " " << " " << temp << " " << loglik << endl;
                }
                else if(fabs(loglik) > 1000000000000){
                    cout << "Inf is here" << endl;
                    cout << temp1 << " " << temp2 << " " << " " << temp << " " << loglik << endl;
                    cout << mavec[*it] << ", " << navec[*it] << ", " << freqc << ", " << drift << ", " << *scurr << ", " << distrec << ", " << 1 << endl;
                }
                else {
                    CurrLik += loglik;
                }
            }
            
            if( CurrLik > BestLik){
                BestS = *scurr;
                BestLik = CurrLik + 0;
            }
        }
        FinalLogLik = BestLik;
    };
    
    vector<double> result = {FinalLogLik,BestS};
    return result;
    
}



// Randomly sample a certain number of elements from a vector of indices
vector<int> SampleIndices(vector<int> Original, int numberSNPs){
    
    // Randomly shuffle elements in window
    vector<int> reshuffled = Original;
    random_shuffle ( reshuffled.begin(), reshuffled.end() );
    
    // Get a sample from the indices
    int samplesize = fmin(Original.size(),numberSNPs);
    vector<int> sampleidx = {};
    int j=0;
    for(j=0; j < samplesize; j++){
        sampleidx.push_back(reshuffled[j]);
    }
    
    return sampleidx;

}


