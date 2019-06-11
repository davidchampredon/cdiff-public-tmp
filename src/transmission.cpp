//
//  transmission.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-28.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include <random>
#include "transmission.hpp"
#include "globalvar.hpp"

void contact_HCW_patient(HCW* hcw, Patient* p, Room* r,
						 float contact_duration,
						 float contact_indiv_halftime,
						 float contact_xchg_contamIdx,
						 float visit_xchg_contamIdx,
						 float visit_room_halftime,
						 Disease D,
						 float curr_time,
						 float timestep){
	
	// HCW shed and acquire spores to/from patient's room:
	shed_in_room(hcw, r, contact_duration, curr_time, contact_indiv_halftime);
	spore_acq_from_room(hcw, r, contact_duration, curr_time, visit_room_halftime, visit_xchg_contamIdx);
	
	// Spores exchange between HCW and patient:
	spore_exchange(hcw, p, contact_duration,
				   contact_indiv_halftime,
				   contact_xchg_contamIdx);
	
	// Potential colonization event:
	
	// From patient to HCW:
	string acqType = p->get_symptomatic()?"patient_symptomatic":"patient_asymptomatic";
	colonization_attempt(hcw, D, curr_time, timestep,acqType);
	
	// From HCW to patient:
	// (HCW is _always_ asymptomatic)
	colonization_attempt(p  , D, curr_time, timestep,"hcw_asymptomatic");
}



double	draw_contact_HCW_patient_duration(double mean_in_minutes,
										  double var_in_minutes){
	double mean = mean_in_minutes / 1440.0;
	double var  = var_in_minutes / 1440.0;
	lognormal_distribution<float> dlnorm(log(mean), var);
	return(dlnorm(RANDOM_GENERATOR));
}


void shed_in_room(Individual* A, Room* r,
				  float visit_duration,
				  float curr_time,
				  float half_time_in_days){
	
	float c = A->get_contamIdx();
	uniform_real_distribution<double> U(0.0, 1.0);
	
	if(c > 0){
		// Draw the amount that is shed in this room from individual 'A':
		float a = log(2.0)/half_time_in_days;
		double shedding_amount = c * U(RANDOM_GENERATOR) * (1-exp(-a * visit_duration));
		
		// Spores are deposited in room:
		r->incr_contamIdx(shedding_amount, curr_time);
		
		// Note: if the contamIdx of the contaminating
		// individual is decreased, the proportion of
		// HAI infections is very (too) small.
		// Hence, the line below is commented, but kept for my reference:
		// A->decr_contamIdx(shedding_amount);
	}
}


void spore_acq_from_room(Individual* indiv,
						 Room* r,
						 float visit_duration,
						 float curr_time,
						 float visit_room_halftime,
						 float visit_xchg_contamIdx){
	
	float Cr = r->get_contamIdx();
	
	if(Cr > 0){
		// The amount of spores that will increase
		// the individual's contamination index
		// is proportional to the visit duration in that room:
		float a = log(2.0)/visit_room_halftime;
		float x =  Cr * visit_xchg_contamIdx * (1-exp(-a * visit_duration));
		
		// Spores are taken from room
		// and deposited on individual:
		r->decr_contamIdx(x);
		indiv->incr_contamIdx(x);
	}
}


float spore_acq_amount_env(float contamIdx, float scaling){
	// TO DO: think of a better way to choose distribution...
	uniform_real_distribution<float> dunif(0.0, contamIdx * scaling);
	return dunif(RANDOM_GENERATOR);
}


void spore_acq_from_other(Individual* A, float envContam_other, float scaling){
	
	if(envContam_other>0){
		float x = spore_acq_amount_env(envContam_other, scaling);
		A->incr_contamIdx(x);
		// DEBUG:
		cout << "Individual #"<< A->get_uid() <<" increased its contamIdx from OTHER environment";
		cout << " by "<< x << endl;
	}
}


void spore_exchange(Individual* x, Individual* y,
					float contact_duration,
					float half_time_in_days,
					float contact_xchg_contamIdx){

	float Cx = x->get_contamIdx();
	float Cy = y->get_contamIdx();
	
	if(Cx>0 || Cy>0){
		// DEBUG:
//		cout << "DEBUG: individuals #"<<x->get_uid()<< " & #" <<y->get_uid();
//		cout << " are exchanging spores. Before exchange:";
//		x->show(); y->show();
		
//		uniform_real_distribution<float> dunif(0.0, 1.0);
//		float u = dunif(RANDOM_GENERATOR);
//		float v = dunif(RANDOM_GENERATOR);
		float pxchg = contact_xchg_contamIdx;
		
		float a = log(2)/half_time_in_days;
		
		float Tx2y = pxchg * Cx * (1-exp(-a * contact_duration));
		float Ty2x = pxchg * Cy * (1-exp(-a * contact_duration));
		
		y->incr_contamIdx(Tx2y);
		// if x is a vector of transmission,
		// there is no creation of new spores for x:
		if(! x->get_isColonized())
			x->decr_contamIdx(Tx2y);
		
		x->incr_contamIdx(Ty2x);
		// if y is a vector of transmission,
		// there is no creation of new spores for y:
		if(! y->get_isColonized())
			y->decr_contamIdx(Ty2x);
		
		// DEBUG:
//		cout << "After exchange:";
//		x->show(); y->show();
	}
}



void colonization_attempt(Individual* A,
						  Disease D,
						  float current_time,
						  float timestep,
						  string acquisitionType){
	
	float contamIdx = A->get_contamIdx();
	
	bool not_colonized 		= !A->get_isColonized();
	bool contamIdx_positive = (contamIdx > 0);
	bool not_csympt			= !A->get_currentlySymptomatic();
	
	uint nr = A->get_nRelapses();
	
	if(not_colonized &&
	   contamIdx_positive &&
	   not_csympt &&
	   nr < MAX_N_RELAPSES)
	{
		float	suscept = A->get_susceptIdx();
		float	p = suscept * contamIdx * timestep;
		bool 	colonization_success = false;
		
		if(p>1) {
			colonization_success = true;
		}
		else{
			bernoulli_distribution dbern(p);
			colonization_success = dbern(RANDOM_GENERATOR);
		}
		if(colonization_success) {
			A->colonize(current_time, D);
			A->set_acquisitionType(acquisitionType);
		}
	}
}






