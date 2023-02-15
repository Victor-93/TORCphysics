// SIDD_3R.h: interface for the SIDD_3R class.
//
// author: chengpeng bi
// last modified: August 1, 2002
//
// This is a C++ implementation of SIDD algorithm designed by Dr. Benham
//
// This class is designed for three runs derived from the class of two runs
//
// Detailed documentation to come
//
// UC Davis Genome Center
///////////////////////////////////////////////////////////////////////////////////////

#ifndef SIDD_3R_H_
#define SIDD_3R_H_

#include "SIDD_2R.h"

class SIDD_3R : public SIDD_2R  
{
protected:
	bool flag_minE3;
	long StatesInThreeRuns;
	
public:
	SIDD_3R();
	virtual ~SIDD_3R();
	bool Search_Low3RE();
    
};

#endif // !defined(AFX_SIDD_3R_H__B5CD0095_10DE_4C22_A817_256ABBAD9425__INCLUDED_)
