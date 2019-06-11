//
//  censusDF.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-20.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "censusDF.hpp"


void censusDF::get_census_from_HCS(HealthCareSetting HCS){
	_uid_individual		= HCS.census_uid_individual;
	_isAlive			= HCS.census_isAlive;
	_willDie			= HCS.census_willDie;
	_colectomy			= HCS.census_colectomy;
	_fmt				= HCS.census_fmt;
	
	_type_individual	= HCS.census_type_individual;
	_subtype_individual	= HCS.census_subtype_individual;
	
	_sIdx				= HCS.census_sIdx;
	_infectIdx			= HCS.census_infectIdx;
	
	_DoS				= HCS.census_DoS;
	_current_room_uid	= HCS.census_current_room_uid;
	_current_room_contamIdx = HCS.census_current_room_contamIdx;
	_time_admission		= HCS.census_timeAdmission;
	_wardAssigned		= HCS.census_wardAssigned;
	
	_isColonized		= HCS.census_isColonized;
	_wasColonized		= HCS.census_wasColonized;
	_symptomatic		= HCS.census_symptomatic;
	_currentlySymptomatic = HCS.census_currentlySymptomatic;
	_prdIncubation		= HCS.census_prdIncubation;
	_prdSymptom			= HCS.census_prdSymptom;
	
	_contamIdx			= HCS.census_contamIdx;
	_contamIdx_other	= HCS.census_contamIdx_other;
	
	_onDuty				= HCS.census_onDuty;
	_hcs_name			= HCS.census_hcs_name;
	
	_isVax				= HCS.census_isVax;
	
	_timeAcquisition	= HCS.census_timeAcquisition;
	_timeOnset			= HCS.census_timeOnset;
	_timeIsolationStart = HCS.census_timeIsolationStart;
	_timeIsolationEnd 	= HCS.census_timeIsolationEnd;
	_timeClearance	 	= HCS.census_timeClearance;
	_timeRelapse		= HCS.census_timeRelapse;
	
	_isoDurFirst		= HCS.census_isoDurFirst;
	_isoDurRelapses		= HCS.census_isoDurRelapses;
	
	_totDurSymptoms		= HCS.census_totDurSymptoms;
	
	_acquisitionType	= HCS.census_acquisitionType;
	
	_nColonizations		= HCS.census_nColonizations;
	_nInfections		= HCS.census_nInfections;
	_nRelapses			= HCS.census_nRelapses;
}


// TO DO: test transfer b/w HCS,
// in particular, only one HCS appear in the census = problem --> Fix it!
censusDF::censusDF(HealthCareSetting HCS,
				   float census_time){
	
	_census_time = census_time;
	HCS.census_individuals();
	get_census_from_HCS(HCS);
}


censusDF::censusDF(HealthCareSetting HCS,
				   float census_time,
				   bool dummy){
	_census_time = census_time;
	HCS.census_patients_light();
	get_census_from_HCS(HCS);
}
