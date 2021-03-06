//
//  main.cpp
//  AncSelFin
//
//  Created by Fernando Racimo on 2/9/14.
//  Copyright (c) 2014 Fernando Racimo. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ASFfunc.h"
#include <iomanip>

// Using cubature package functions
#include "cubature.h"

// Using multiple precision library
//#include "mpir.h"

using namespace std;


// Can't parallelize because of these global variables needed by qromb...
double p1, p2, p3, drift1, drift2, drift3, selcoeff123, recdist123, branch123;


//size_t size = 1;

// data vectors
vector<string> chrvec;
vector<double> physvec;
vector<double> posvec;
vector<double> mavec;
vector<double> navec;
vector<double> mbvec;
vector<double> nbvec;
vector<double> mcvec;
vector<double> ncvec;

// Define frequency vectors
vector<double> freqA;
vector<double> freqB;
vector<double> freqC;

// Define sum vectors
vector<double> mabvec;
vector<double> nabvec;

// drift values
vector<double> drifts;

// Effective pop. size
double Ne = 10000;

// total number of sites
int totalsites;

// Vector of possible selection coefficients
//vector<double> sPossible = {0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1};
vector<double> sPossible = {0.00000001,0.0000001,0.000001,0.00001,0.0001,0.001,0.01,0.1};


int main(int argc, const char * argv[])
{
    
    // GET PARAMETERS FROM COMMAND LINE
    // Infile name
    string infilename = argv[1];
    cout << infilename << endl;
    // Outfile name
    string outfilename = argv[2];
    cout << outfilename << endl;
    // Space between analyzed SNPs
    int SNPspace = stoi(argv[3]);
    cout << SNPspace << endl;
    // Number of SNPs per window
    int SNPsperWin = stoi(argv[4]);
    cout << SNPsperWin << endl;
    // Window size in Morgans
    double winlen = stod(argv[5]);
    cout << winlen << endl;
    // Drifts
    const string drifts_unsplit = argv[6];
    char delim = ',';
    stringstream ss(drifts_unsplit);
    string item;
    while (getline(ss, item, delim)) {
        drifts.push_back(stod(item));
    }
    // Weights file
    string weightsfilename = argv[7];
    cout << weightsfilename << endl;
    
    //double test = sistsamp_given_outfreq(98,100,98,100,0.5,0.2,0.2,2.4,10,0.00001,1,0);
    //double test = sistfreq_given_outfreq(0.99,0.99,0.5,0.2,0.2,2.4,0.1,0.001,1,0);
    //double test2 = sistsamp_given_outfreq(99,100,99,100,0.5,0.2,0.2,2.4,0.1,0.001,1,0);
    //double test3 = sistsamp_given_outfreq(100,100,100,100,0.5,0.2,0.2,2.4,0.1,0.001,1,0);
    //double test4 = Binom(100,100);
    //cout << test << endl;
    //cout << test2 << endl;
    //cout << test3 << endl;
    //cout << test4 << endl;
    

    
    // Open input file and read table
    string line;
    ifstream infile (infilename);
    if (infile.is_open())
    {
        // get header
        getline (infile,line);
        
        // read values
        int i = 0;
        while ( getline (infile,line) )
        {
            //store values in vectors;
            stringstream ssin(line);
            int j = 0;
            while (ssin.good() && j < 10){
                
                physvec.push_back(1);
                posvec.push_back (1);
                mavec.push_back (1);
                navec.push_back (1);
                mbvec.push_back (1);
                nbvec.push_back (1);
                mcvec.push_back (1);
                ncvec.push_back (1);
                
                string test;
                if (j == 0) { ssin >> test; chrvec.push_back(test);}
                else if (j == 1) { ssin >> physvec[i];}
                else if (j == 2) { ssin >> posvec[i];}
                else if (j == 3){ ssin >> mavec[i];}
                else if (j == 4){ ssin >> navec[i];}
                else if (j == 5){ ssin >> mbvec[i];}
                else if (j == 6){ ssin >> nbvec[i];}
                else if (j == 7){ ssin >> mcvec[i];}
                else if (j == 8){ ssin >> ncvec[i];}
                ++j;
            }
            freqA.push_back(mavec[i]/navec[i]);
            freqB.push_back(mbvec[i]/nbvec[i]);
            freqC.push_back(mcvec[i]/ncvec[i]);
            mabvec.push_back(mavec[i]+mbvec[i]);
            nabvec.push_back(navec[i]+nbvec[i]);
            
            ++i;
        }
        totalsites = i;
        infile.close();
    }
    else cout << "Unable to open input file\n";
 
    // Open weights file and fill up weights vector
    vector<double> weightsvec(totalsites);
    if(weightsfilename == "NA"){
        fill (weightsvec.begin(),weightsvec.begin()+totalsites,1);
    }
    else{
        ifstream weightsfile (weightsfilename);
    
        // Open weights file
        if (weightsfile.is_open())
        {
            // Get header
            getline (weightsfile,line);
        
            // read values
            int i = 0;
            while ( getline (weightsfile,line) )
            {
                //store values in vectors;
                stringstream ssin(line);
                int j = 0;
                while (ssin.good() && j < 4){
                    
                    weightsvec.push_back(1);
                    string wchr;
                    string wpos;
                    if (j == 0) { ssin >> wchr; }
                    else if (j == 1) { ssin >> wpos;}
                    else if (j == 2) { ssin >> weightsvec[i];}
                    ++j;
                }
                ++i;
            }
        }
        else cout << "Unable to open weights file\n";
    }
    
    
    // Open output file
    ofstream outfile;
    outfile.open (outfilename);

    
    // Define windows of genome
    vector<double> winstartvec;
    vector<double> winendvec;

    // SNP window 2D vectors
    vector< vector<int> > allSNPwindows;
    vector< vector<int> > Csegwindows;
    vector< vector<int> > Bsegwindows;
    vector< vector<int> > Asegwindows;
    allSNPwindows.resize(totalsites);
    Csegwindows.resize(totalsites);
    Bsegwindows.resize(totalsites);
    Asegwindows.resize(totalsites);
    
    // Iterate over all possible central SNPs
    int i=0;
    for(i = 0; i < totalsites; i++){
    
        // Get ideal left and right boundaries
        double idealleft = fmax(posvec[i] - (winlen / 2.0),0.0);
        double idealright = fmin(posvec[i] + (winlen/2.0),posvec[(totalsites-1)]);
        
        // Compute distance of all SNPs to site (negative = left, positive = right)
        int j=0;
        double best_disttoleft = 100;
        double best_leftidx = 0;
        double best_disttoright = 100;
        double best_rightidx = totalsites-1;
        
        // Get actual left and right boundaries
        for(j = 0; j < totalsites; j++){
            double curr_disttoleft = fabs(posvec[j] - idealleft);
            if(curr_disttoleft < best_disttoleft){
                best_disttoleft = curr_disttoleft;
                best_leftidx = j;
            }
            double curr_disttoright = fabs(posvec[j] - idealright);
            if(curr_disttoright < best_disttoright){
                best_disttoright = curr_disttoright;
                best_rightidx = j;
            }
            
        }
        winstartvec.push_back(best_leftidx);
        winendvec.push_back(best_rightidx);
        //cout << best_leftidx << " " << i << " " << best_rightidx << endl;
        
        // Get all window SNPs
        int h;
        for(h = best_leftidx; h <= best_rightidx; h++){
            //cout << i << ", " << h << endl;
            allSNPwindows[i].push_back(h);
            if(freqC[i] > 0.0 && freqC[i] < 1.0){
                Csegwindows[i].push_back(h);
                //cout << i << ", " << h << ", " << Csegwindows[i].size() << endl;
            }
            if(freqB[i] > 0.0 && freqB[i] < 1.0){
                Bsegwindows[i].push_back(h);
            }
            if(freqA[i] > 0.0 && freqA[i] < 1.0){
                Asegwindows[i].push_back(h);
            }
            
        }
    
    };
    
    cout << "Windows determined. Beginning statistic computation..." << endl;
    
    // Print header
    outfile << "Chr" << " \t" << "PhysPos" << "\t" << "GenPos" << "\t" << "PhysStart" << "\t" << "PhysEnd" << "\t" << "GenStart" << "\t" << "GenEnd" << "\t" << "3PCLR.Anc" << "\t" << "3PCLR.Anc.S" << "\t" << "3PCLR.A" << "\t" << "3PCLR.A.S" << "\t" << "3PCLR.B" << "\t" << "3PCLR.B.S" << "\t" << "XPCLRAC" << "\t" << "XPCLRAC.S" << "\t" << "XPCLRBC" << "\t" << "XPCLRBC.S" << endl;
    
    // Iterate over SNP every X SNPs
    i = 0;
    vector<int> Csampleidx;
    vector<int> Bsampleidx;
    vector<int> Asampleidx;
    double THREEPCLR_ratio;
    double THREEPCLR_bestS;
    double THREEPCLRA_ratio;
    double THREEPCLRA_bestS;
    double THREEPCLRB_ratio;
    double THREEPCLRB_bestS;
    double XPCLRAC_ratio;
    double XPCLRAC_bestS;
    double XPCLRBC_ratio;
    double XPCLRBC_bestS;
    //double XPCLRavg_ratio;
    //double XPCLRavg_bestS;
    //double XPCLRAB_ratio;
    //double XPCLRAB_bestS;
    //double XPCLRBA_ratio;
    //double XPCLRBA_bestS;
    vector<double> selection;
    vector<double> neutrality;
    for(i=0; i < totalsites; i+=SNPspace){
        
        // clock_t start = clock();
        
        // Sample indices
        Csampleidx = SampleIndices(Csegwindows[i],SNPsperWin);
        Bsampleidx = SampleIndices(Bsegwindows[i],SNPsperWin);
        Asampleidx = SampleIndices(Asegwindows[i],SNPsperWin);
        
        
        // 3P-CLR
        selection = SupLogLikTri(Csampleidx, posvec, mavec, navec, mbvec, nbvec, mcvec, ncvec, drifts[0], drifts[1], drifts[2], posvec[i], sPossible, weightsvec, 1,0);
        neutrality = SupLogLikTri(Csampleidx, posvec, mavec, navec, mbvec, nbvec, mcvec, ncvec, drifts[0], drifts[1], drifts[2], posvec[i], sPossible, weightsvec, 0,0);
        THREEPCLR_ratio = 2 * (selection[0] - neutrality[0]);
        THREEPCLR_bestS = selection[1];
        
        
        // 3P-CLR - A BRANCH
        selection = SupLogLikTri(Csampleidx, posvec, mavec, navec, mbvec, nbvec, mcvec, ncvec, drifts[0], drifts[1], drifts[2], posvec[i], sPossible, weightsvec, 1,1);
        neutrality = SupLogLikTri(Csampleidx, posvec, mavec, navec, mbvec, nbvec, mcvec, ncvec, drifts[0], drifts[1], drifts[2], posvec[i], sPossible, weightsvec, 0,1);
        THREEPCLRA_ratio = 2 * (selection[0] - neutrality[0]);
        THREEPCLRA_bestS = selection[1];
        //cout << selection[0] << endl;
        //cout << neutrality[0] << endl;
        //cout << THREEPCLRA_ratio << endl;
        //cout << selection[1] << endl;
        //cout << " " << endl;
        

        // 3P-CLR - B BRANCH
        selection = SupLogLikTri(Csampleidx, posvec, mavec, navec, mbvec, nbvec, mcvec, ncvec, drifts[0], drifts[1], drifts[2], posvec[i], sPossible, weightsvec, 1,2);
        neutrality = SupLogLikTri(Csampleidx, posvec, mavec, navec, mbvec, nbvec, mcvec, ncvec, drifts[0], drifts[1], drifts[2], posvec[i], sPossible, weightsvec, 0,2);
        THREEPCLRB_ratio = 2 * (selection[0] - neutrality[0]);
        THREEPCLRB_bestS = selection[1];
        

        // XP-CLR A vs. C
        selection = SupLogLikChen(Csampleidx, posvec, mavec, navec, mcvec, ncvec, drifts[0]+drifts[2], posvec[i], sPossible, weightsvec, 1);
        neutrality = SupLogLikChen(Csampleidx, posvec, mavec, navec, mcvec, ncvec, drifts[0]+drifts[2], posvec[i], sPossible, weightsvec, 0);
        XPCLRAC_ratio = 2 * (selection[0] - neutrality[0]);
        XPCLRAC_bestS = selection[1];


        // XP-CLR B vs. C
        selection = SupLogLikChen(Csampleidx, posvec, mbvec, nbvec, mcvec, ncvec, drifts[1]+drifts[2], posvec[i], sPossible, weightsvec, 1);
        neutrality = SupLogLikChen(Csampleidx, posvec, mbvec, nbvec, mcvec, ncvec, drifts[1]+drifts[2], posvec[i], sPossible, weightsvec, 0);
        XPCLRBC_ratio = 2 * (selection[0] - neutrality[0]);
        XPCLRBC_bestS = selection[1];
        if(XPCLRBC_ratio != XPCLRBC_ratio){
            cout << selection[0] << " " << neutrality [0] << " " << XPCLRBC_ratio << endl;
        }
        

        // XP-CLR AB avg vs C
        //selection = SupLogLikChen(Csampleidx, posvec, mabvec, nabvec, mcvec, ncvec, ((drifts[0]+drifts[1])/2)+drifts[2], posvec[i], sPossible, weightsvec, 1);
        //neutrality = SupLogLikChen(Csampleidx, posvec, mabvec, nabvec, mcvec, ncvec, ((drifts[0]+drifts[1])/2)+drifts[2], posvec[i], sPossible, weightsvec, 0);
        //XPCLRavg_ratio = 2 * (selection[0] - neutrality[0]);
        //XPCLRavg_bestS = selection[1];
        
        // XP-CLR A vs. B
        //selection = SupLogLikChen(Bsampleidx, posvec, mavec, navec, mbvec, nbvec, drifts[0]+drifts[1], posvec[i], sPossible, weightsvec, 1);
        //neutrality = SupLogLikChen(Bsampleidx, posvec, mavec, navec, mbvec, nbvec, drifts[0]+drifts[1], posvec[i], sPossible, weightsvec, 0);
        //XPCLRAB_ratio = 2 * (selection[0] - neutrality[0]);
        //XPCLRAB_bestS = selection[1];

        // XP-CLR B vs. A
        //selection = SupLogLikChen(Asampleidx, posvec, mbvec, nbvec, mavec, navec, drifts[0]+drifts[1], posvec[i], sPossible, weightsvec, 1);
        //neutrality = SupLogLikChen(Asampleidx, posvec, mbvec, nbvec, mavec, navec, drifts[0]+drifts[1], posvec[i], sPossible, weightsvec, 0);
        //XPCLRBA_ratio = 2 * (selection[0] - neutrality[0]);
        //XPCLRBA_bestS = selection[1];
        
        
        
        //cout << THREEPCLR_ratio << "\t" << THREEPCLRA_ratio << "\t" << THREEPCLRB_ratio << "\t" << XPCLRAC_ratio << "\t" << XPCLRBC_ratio << endl;

        
        cout << chrvec[i] << " \t" << physvec[i] << "\t" << posvec[i] << "\t" << physvec[winstartvec[i]] << "\t" << physvec[winendvec[i]] << "\t" << posvec[winstartvec[i]] << "\t" << posvec[winendvec[i]] << "\t" << THREEPCLR_ratio << "\t" << THREEPCLR_bestS << "\t" << THREEPCLRA_ratio << "\t" << THREEPCLRA_bestS << "\t" << THREEPCLRB_ratio << "\t" << THREEPCLRB_bestS << "\t" << XPCLRAC_ratio << "\t" << XPCLRAC_bestS << "\t" << XPCLRBC_ratio << "\t" << XPCLRBC_bestS << endl;
        
        outfile << chrvec[i] << " \t" << physvec[i] << "\t" << posvec[i] << "\t" << physvec[winstartvec[i]] << "\t" << physvec[winendvec[i]] << "\t" << posvec[winstartvec[i]] << "\t" << posvec[winendvec[i]] << "\t" << THREEPCLR_ratio << "\t" << THREEPCLR_bestS << "\t" << THREEPCLRA_ratio << "\t" << THREEPCLRA_bestS << "\t" << THREEPCLRB_ratio << "\t" << THREEPCLRB_bestS << "\t" << XPCLRAC_ratio << "\t" << XPCLRAC_bestS << "\t" << XPCLRBC_ratio << "\t" << XPCLRBC_bestS << endl;
        
        
        
        //cout << chrvec[i] << " \t" << physvec[i] << "\t" << posvec[i] << "\t" << physvec[winstartvec[i]] << "\t" << physvec[winendvec[i]] << "\t" << posvec[winstartvec[i]] << "\t" << posvec[winendvec[i]] << "\t" << THREEPCLR_ratio << "\t" << THREEPCLR_bestS << "\t" << XPCLRAC_ratio << "\t" << XPCLRAC_bestS << "\t" << XPCLRBC_ratio << "\t" << XPCLRBC_bestS << "\t" << XPCLRavg_ratio << "\t" << XPCLRavg_bestS << "\t" << XPCLRAB_ratio << "\t" << XPCLRAB_bestS << "\t" << XPCLRBA_ratio << "\t" << XPCLRBA_bestS << endl;
        
        //clock_t end = clock();
        //double time = (double) (end-start) / CLOCKS_PER_SEC * 1000.0;
        //cout << time << endl;
        //cout << " " << endl;

        
    }
    
    
    outfile.close();

    return 0;
 
    
}

