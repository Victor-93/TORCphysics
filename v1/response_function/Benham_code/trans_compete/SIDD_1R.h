// SIDD_1R.h: interface for the SIDD_1R class.
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

#ifndef SIDD_1R_H_
#define SIDD_1R_H_

#include "SIDD_Base.h"
#include "stat_1R.h"

using namespace std;

typedef	list<stat_1R> LstState;  

class SIDD_1R : public SIDD_Base  
{
protected:
	LstState Lst_OBE[MaxInitialWindowSize+1][3]; // lists of open base energy
	bool flag_minE1; // flagging if new min_E found in one run
	double min_RE[3]; // store minimum run energy (a + NI)
    long StatesInOneRun;
	
public:
	SIDD_1R();
	virtual ~SIDD_1R();
	void gen_OpenBaseEnergy();
	bool Search_Low1RE();
	
	void Update_Low1RE();
	void Show_Low1RE();
	double runs[3];
    double Prob[3];
    double exp_one[3];
    double exp_two[3][3];
    double exp_three[3][3][3];
    double exp_four[3][3][3][3];

	double exp_type_B;

	
};

#endif // !defined(AFX_SIDD_1R_H__2688B961_9B45_44D5_8F6F_F45636E04460__INCLUDED_)
