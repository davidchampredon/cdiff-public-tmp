//
//  Vaccine.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2018-06-15.
//  Copyright © 2018 David CHAMPREDON. All rights reserved.
//

#ifndef Vaccine_hpp
#define Vaccine_hpp

#include <stdio.h>
#include <iostream>

using namespace std;

class Vaccine{
	
	string 	_name;
	
	// Range for vaccine efficacy:
	float	_efficacy_lo;
	float	_efficacy_hi;
	
	// proba to set susceptibility index to 0:
	float	_proba_prevent_colonization;
	
	float	_halfLife;

public:
	
	Vaccine(){}
	
	Vaccine(string name,
			float efficacy_lo,
			float efficacy_hi,
			float proba_prevent_colonization,
			float halfLife);
	
	
	// ---- Get functions ----
	
	string	get_name() const {return _name;}
	float	get_efficacy_lo() const {return _efficacy_lo;}
	float	get_efficacy_hi() const {return _efficacy_hi;}
	float	get_proba_prevent_colonization() const {return _proba_prevent_colonization;}
	float	get_halfLife() const {return _halfLife;}
	
	// ---- Miscellenaous ----
	float	draw_efficacy() const;
};

#endif /* Vaccine_hpp */
