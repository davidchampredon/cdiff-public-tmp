//
//  utils.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "utils.hpp"
#include "globalvar.hpp"


void stopif(bool condition,
			string error_msg,
			int error_code){
	if (condition)
	{
		cerr << endl << endl;
		cerr <<	endl <<	" *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* ";
		cerr << endl << " *=*=*=*=*=*=*  MODEL ERROR  *=*=*=*=*=*=* ";
		cerr <<	endl <<	" *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* ";
		cerr << endl << endl;
		cerr << error_msg <<endl;
		cerr <<	endl <<	" *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*";
		cerr <<	endl <<	" *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*" << endl;
		throw error_code;
	}
}

float beta_distribution(float a, float b){
	std::gamma_distribution<float> dgam1(a,1);
	std::gamma_distribution<float> dgam2(b,1);
	float x = dgam1(RANDOM_GENERATOR);
	float y = dgam2(RANDOM_GENERATOR);
	float res = x / (x+y);
	return res;
}
