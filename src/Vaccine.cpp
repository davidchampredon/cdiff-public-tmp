//
//  Vaccine.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2018-06-15.
//  Copyright © 2018 David CHAMPREDON. All rights reserved.
//

#include <random>
#include "Vaccine.hpp"
#include "globalvar.hpp"

Vaccine::Vaccine(string name,
				 float efficacy_lo,
				 float efficacy_hi,
				 float proba_prevent_colonization,
				 float halfLife){
	_name = name;
	_efficacy_lo = efficacy_lo;
	_efficacy_hi = efficacy_hi;
	_proba_prevent_colonization = proba_prevent_colonization;
	_halfLife = halfLife;
}

float Vaccine::draw_efficacy() const {
	uniform_real_distribution<float> runif(_efficacy_lo, _efficacy_hi);
	float x = runif(RANDOM_GENERATOR);
	return x;
}
