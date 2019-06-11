//
//  transmission.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-28.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef transmission_hpp
#define transmission_hpp

#include <stdio.h>

#include "Individual.hpp"
#include "Room.hpp"


void contact_HCW_patient(HCW* hcw, Patient* p, Room* r,
						 float contact_duration,
						 float contact_indiv_halftime,
						 float contact_xchg_contamIdx,
						 float visit_xchg_contamIdx,
						 float visit_room_halftime,
						 Disease D,
						 float curr_time,
						 float timestep);

/// Return a random duration for the contact between a HCW and a patient
double	draw_contact_HCW_patient_duration(double mean_in_minutes, double var_in_minutes);

/// Spore shedding from colonized individual.
void	shed_in_room(Individual* X, Room* r,
					 float visit_duration,
					 float curr_time,
					 float halft_time_in_days);

/// Amount of spores acquired from _any_ environment.
float	spore_acq_amount_env(float contamIdx, float scaling);

/// Spores acquisition by an individual from a contaminated room.
void	spore_acq_from_room(Individual* X, Room* r,
							float visit_duration, // unit = days
							float curr_time,
							float visit_room_halftime, // unit = days
							float visit_xchg_contamIdx);

/// Spores acquisition from contaminated environment, other than patient room.
void	spore_acq_from_other(Individual* X, float envContam_other, float scaling);

/** Increase contamination potential upon
 contact between x and y. NOTE: this function takes care
 of both ways x->y and y->x. */
void	spore_exchange(Individual* x, Individual* y,
					   float contact_duration,
					   float half_time_in_days,
					   float contact_xchg_contamIdx);

/// Draw the colonization event.
void	colonization_attempt(Individual* A,
							 Disease D,
							 float curr_time,
							 float timestep,
							 string acquisitionType);


#endif /* transmission_hpp */
