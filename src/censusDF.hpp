//
//  censusDF.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-20.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef censusDF_hpp
#define censusDF_hpp

#include <stdio.h>

#include "HealthCareSetting.hpp"


class censusDF{
	
public:
	
	float				_census_time;
	
	vector<string>		_hcs_name;
	
	vector<uint>		_uid_individual;
	vector<string>		_type_individual;
	vector<string>		_subtype_individual;
	vector<bool>		_isAlive;
	vector<bool>		_willDie;
	vector<bool>		_colectomy;
	vector<bool>		_fmt;
	vector<float>		_sIdx;
	vector<float>		_infectIdx;
	vector<float>		_time_admission;
	vector<float>		_DoS;
	vector<uint>		_current_room_uid;
	vector<float>		_current_room_contamIdx;
	vector<uint>		_home_room_uid;
	vector<uint>		_wardAssigned;
	
	vector<bool>		_isColonized;
	vector<bool>		_wasColonized;
	vector<bool>		_symptomatic;
	vector<bool>		_currentlySymptomatic;
	vector<float>		_prdIncubation;
	vector<float>		_prdSymptom;
	
	vector<float>		_contamIdx;
	vector<float>		_contamIdx_other;
	
	vector<bool>		_onDuty;
	
	vector<bool>		_isVax;
	
	vector<float>		_timeAcquisition;
	vector<float>		_timeOnset;
	vector<float>		_timeIsolationStart;
	vector<float>		_timeIsolationEnd;
	vector<float>		_timeClearance;
	vector<float>		_timeRelapse;
	
	vector<float>		_isoDurFirst;
	vector<float>		_isoDurRelapses;
	
	vector<float>		_totDurSymptoms;
	
	vector<string>		_acquisitionType;
	
	vector<uint>		_nColonizations;
	vector<uint>		_nInfections;
	vector<uint>		_nRelapses;
	
	censusDF(){}
	censusDF(HealthCareSetting HCS, float census_time);
	censusDF(HealthCareSetting HCS, float census_time, bool dummy);
	
	void		get_census_from_HCS(HealthCareSetting HCS);
	
	
	
};


#endif /* censusDF_hpp */
