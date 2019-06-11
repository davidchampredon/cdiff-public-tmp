//
//  Simulator.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-15.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef Simulator_hpp
#define Simulator_hpp

#include <stdio.h>
#include <vector>
#include <algorithm>

#include "HealthCareSetting.hpp"
#include "censusDF.hpp"


class Simulator{
	
private:
	
	float	_horizon;
	float	_timeStep;
	
	vector<HealthCareSetting*> _HCS;
	vector<censusDF>           _census;
	vector<censusDF>           _census_light;
	
	bool	_recordContacts;  // Switch to record contact durations. If "on", then memory intensive.
	
	// Square matrix defining the movements between HCS.
	// Rows are 'from', columns are 'to';
	// Element value is the intensity of a Poisson process.
	// _mvt_HCS[i,j] = intensity from _HCS[i] to _HCS[j].
	// Diagonal values not used. 
	vector<vector<float> > _mvt_HCS;
	
public:
	
	Simulator(){}
	
	
	Simulator(vector<HealthCareSetting*> hcs,
			  float horizon,
			  float timeStep,
			  bool recordContacts){
		_HCS 			= hcs;
		_horizon		= horizon;
		_timeStep		= timeStep;
		_recordContacts	= recordContacts;
		_census.clear();
		_census_light.clear();
	}
	
	
	// ==== SET FUNCTIONS ====
	
	void set_mvt_HCS(vector<vector<float> > x){_mvt_HCS = x;}
	
	
	// ==== GET FUNCTIONS ====
	
	vector<HealthCareSetting*> 	get_HCS(){return _HCS;}
	vector<censusDF> 			get_census(){return _census;}
	vector<censusDF> 			get_census_light(){return _census_light;}
	
	
	// ==== MISCELLANEOUS ====
	
	/// Calculates the period of the day, in _timeStep units.
	uint period_of_day(float t);
	
	/// Calculate the total number of periods in one day.
	uint total_periods_in_day();
	
	/// Randomly assign periods worked in a day to _all_ HCWs, for the ith HealthCareSetting.
	void assign_HCW_periodsWorked(uint i, float hoursWorkedInOneDay);
	
	/// Number of HCS in this simulator.
	size_t n_HCS() {return _HCS.size();}
	
	
	void set_current_time(float t);
	
	// ==== MOVEMENTS ====

	/** Draw the number of transfers between HCSs.
	 *  Returns a square matrix of the same dimensions as 
	 *  '_mvt_HCS' with elements representing the number of patients transfered.
	*/
	vector<vector<uint> > draw_n_transfer();
	
	
	/** Transfer one patient from HCS[i], ward, room, patient's UID
	 *  to an available room in HCS[j]
	 */
	void transfer(uint i,
				  uint j,
				  uint ward,
				  uint room,
				  uint patient_uid);
	
	/// Transfer all patients between all facilities.
	void transfer_all_patients(float t, int trace_level);
	
	
	// ==== SIMULATION ====
	
	/// Run all events during a time step for _one_ HCS.
	void events_one_HCS(uint i, float t, int trace_level);
	
	/// Run all events during a time step for _all_ HCS.
	void events_all_HCS(float t, int trace_level);
	
	/// Run simulation on all HealthCareSetting.
	void run_simulation(uint seed = 0);

	

	
};


#endif /* Simulator_hpp */
