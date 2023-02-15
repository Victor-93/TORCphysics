// SIDD_1R.cpp: implementation of the SIDD_1R class.
// This class is designed for one run derived from SIDD_Base
//
// program to implement algorithm developed by Craig Benham
//
// author: Chengpeng Bi
// modifiers: Dina Zhabinskaya, Sally Madden, Ian Korf
// compiler: g++
//
// This is a test version. The author has no responsibility for any outcome
// that incurs when you make any trial run.
//
// UC Davis Genome Center
//////////////////////////////////////////////////////////////////////////////

#ifndef SIDD_1R_CPP
#define SIDD_1R_CPP

#include "SIDD_1R.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SIDD_1R::SIDD_1R():SIDD_Base() {
	for (int t=0; t<=Ns; t++) 
		min_RE[t] = min_E;  //minimun energy from the continuous part 
        flag_minE1 = false;
        StatesInOneRun = 0;
}
SIDD_1R::~SIDD_1R()
{
}

void SIDD_1R::gen_OpenBaseEnergy()
{
    double e;
    //t=0 for melting, t=1 for Z-DNA, and t=2 for cruciforms
	for (int t=0; t<=Ns; t++) { 
		for(int i = MinWindowSize; i <= MaxWindowSize; i++){ // window size: 1 - MaxWindowSize
			for(int j = 0; j < length_seq; j++){ // start position to length_seq - 1
                e = a[t] + calc_OPenBasesEnergy(j, i, t);  
                stat_1R  s1r(j, e, RT);	 //input j=starting position, e=energy, RT=constant			
				Lst_OBE[i][t].push_back(s1r);
				if(e < min_RE[t]) min_RE[t] = e; 
			}
		}
	}
	
	// sorting each list of OBE
	for (int t=0; t<=Ns; t++) {
		for(int k = MinWindowSize; k <= MaxWindowSize; k++){	
			Lst_OBE[k][t].sort();
		}
	}
}


bool SIDD_1R::Search_Low1RE()
{
	long count = 0;
	flag_minE1 = false;
	sum_1RG = sum_1RB = 0.0;
	for (int t=0; t<=Ns; t++) {
		exp_one[t] = 0; 
		runs[t]=0; 
        Prob[t]=0;
    }
	LstState::iterator j;
	//collecting one run states withing the energy threshold
	for (int t=0; t<=Ns; t++) {
		for(int i = MinWindowSize; i <= MaxWindowSize; i++){ // window size: 1 - MaxWindowSize
			for(j = Lst_OBE[i][t].begin(); j != Lst_OBE[i][t].end(); ++j){
				double e = j->get_energy(); 
               if (e >= 10000)
                   break;
				int pos =j->get_pos1();
				e += calc_Gres(i*delta_fnc(t,0),i*delta_fnc(t,1),i*delta_fnc(t,2),t*delta_fnc(t,1));   
				if(e >= max_E) break; //outside of threshold
               				count++; //counting the number of one run states
				if(e < min_E){ 
					min_E = e; //update min_E
					max_E = e + theta;
					min_WS = i; // the window size corresponding to the minimum energy
					flag_minE1 = true;
				}
				if(write_profile){  
                    double exponent= exp(-e/RT);
                    exp_one[t] += exponent;
					runs[t]=exp_one[t];
					update_promatrix(pos, i, t, e, exponent);
					sum_1RB += exponent;  
					sum_1RG = sum_1RG + e*exponent; 
				}                  
			}
		}
	}
	
	StatesInOneRun = count;
	// adding the zero open bases term to the partition function
	double closed = exp(-alpha*alpha*K/2/RT);
	ZsumB = sum_1RB + closed;    
	ZsumG = sum_1RG+alpha*alpha*K/2*closed;  
    if(write_profile && results) {
		cout << "Number of one-run states = " << count << endl;        
        for (int t=0; t<=Ns; t++) {
            Prob[t] += exp_one[t];
        }
    }
    return flag_minE1;
}

#endif

