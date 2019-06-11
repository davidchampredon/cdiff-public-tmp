//
//  Disease.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-07-07.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef Disease_hpp
#define Disease_hpp

#include <stdio.h>
#include <iostream>

using namespace std;

class Disease{
	
	string	_name;
	
	float	_prd_incubation_mean;
	float	_prd_incubation_var;
	
	float	_prd_symptom_mean;
	float	_prd_symptom_var;
	
	// Clearance:
	float	_proba_clear_prm;
	float	_prd_clear_asymptom_mean; // length of time before clearance for asymptomatic.
	float	_prd_clear_asymptom_var;
	float	_prd_clear_symptom_mean;  // length of time before clearance for symptomatic infections.
	float	_prd_clear_symptom_var;
	
	// Level of contamination index upon acquisition:
	float	_contamIdx_max_asymptom;
	float	_contamIdx_ratio_symptomatic;
	
	float	_hazard_relapse;	// hazard to relapse
	float	_shape_relapse;		// shape the distribution of the number of relapses
	
	float	_infectIdx_init; 	// initial probability for C. diff infection (i.e. symptomatic) given colonization.
	
public:
	
	Disease(){}
	
	Disease(float prd_incubation_mean,
			float prd_incubation_var,
			float prd_symptom_mean,
			float prd_symptom_var,
			float proba_clear_prm,
			float prd_clear_asymptom_mean,
			float prd_clear_asymptom_var,
			float prd_clear_symptom_mean,
			float prd_clear_symptom_var,
			float contamIdx_asymptomatic,
			float contamIdx_symptomatic,
			float hazard_relapse,
			float shape_relapse,
			float infectIdx_init);
	
	
	float	get_prd_incubation_mean()	{return _prd_incubation_mean;}
	float	get_prd_incubation_var()	{return _prd_incubation_var;}
	
	float	get_prd_symptom_mean()		{return _prd_symptom_mean;}
	float	get_prd_symptom_var()		{return _prd_symptom_var;}
	
	float	get_proba_clear_prm()			{return _proba_clear_prm;}
	float	get_prd_clear_asymptom_mean()	{return _prd_clear_asymptom_mean;}
	float	get_prd_clear_asymptom_var()	{return _prd_clear_asymptom_var;}
	float	get_prd_clear_symptom_mean()	{return _prd_clear_symptom_mean;}
	float	get_prd_clear_symptom_var()		{return _prd_clear_symptom_var;}
	
	float	get_contamIdx_max_asymptom()		{return _contamIdx_max_asymptom;}
	float	get_contamIdx_ratio_symptomatic()	{return _contamIdx_ratio_symptomatic;}
	
	float	get_hazard_relapse()	{return _hazard_relapse;}
	float	get_shape_relapse()		{return _shape_relapse;}
	float	get_infectIdx_init()	{return _infectIdx_init;}
	
	void	show();
};



#endif /* Disease_hpp */



