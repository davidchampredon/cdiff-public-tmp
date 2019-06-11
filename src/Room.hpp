//
//  Room.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef Room_hpp
#define Room_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

#include "globalvar.hpp"
#include "Individual.hpp"

using namespace std;

class Room{
	
protected:
	
	uint	_uid;		// unique ID of this room
	string	_type;		// "patient_room", "other"  ("ICU"?)
	string	_name;
	
	uint	_ward_uid;	// UID of the ward this room belongs to. '0' if doesn't belong to any ward.
	
	// C. difficile overall environmental contamiation
	// in this room. Probability to acquire
	// C. diff from environment only:
	float			_contamIdx;
	
	// Amount and time of environmental contamination.
	// Record needed in order to apply time decay.
	vector<float>	_contamIdx_amount;
	vector<float>	_contamIdx_time;
	
	vector<Patient*>	_patient;
	vector<HCW*>		_hcw;
	
	uint _max_patient; // maximum capacity for patients
	
	// Half time for environemental contamination:
	float	_decay_room_halftime;
	
	// Time when extra cleaning should take place,
	// following symptoms onset of a patient in that room:
	float	_time_extra_cleaning;
	
public:
	
	// === CONSTRUCTORS ===
	
	Room(){
		_uid = NA_UINT;
		_name = "undefined";
		_contamIdx = 0;
		_max_patient = NA_UINT;
		_time_extra_cleaning = -NA_FLOAT;
	}
	
	
	// === SET FUNCTIONS ===
	
	void set_uid(uint x){_uid = x;}
	void set_type(string x){_type = x;}
	void set_ward_uid(uint x){_ward_uid = x;}
	void set_name(string x){_name = x;}
	void set_max_patient(uint x){_max_patient = x;}
	void set_decay_room_halftime(float x){_decay_room_halftime = x;}
	void set_time_extra_cleaning(float x){_time_extra_cleaning = x;
		/*DEBUG*/
		//cout<<"Room #"<<_uid<<" should be extracleaned at "<<_time_extra_cleaning<<endl;
	}
	
	
	// === GET FUNCTIONS ===
	
	uint		get_uid(){return _uid;}
	uint		get_ward_uid(){return _ward_uid;}
	string		get_type(){return _type;}
	uint		get_max_patient(){return _max_patient;}
	float		get_contamIdx(){return _contamIdx;}
	float		get_decay_room_halftime(){return _decay_room_halftime;}
	string		get_name(){return _name;}
	
	vector<Patient*>	get_patient(){return _patient;}
	Patient*			get_patient(uint i){return _patient[i];}
	Patient*			get_patient_by_uid(uint uid);
	uint				get_patient_uid(uint i){return _patient[i]->get_uid();}
	vector<HCW*>		get_hcw(){return _hcw;}
	HCW*				get_hcw(uint i){return _hcw[i];}
	
	float		get_time_extra_cleaning(){return _time_extra_cleaning;}
	
	
	// === CENSUS ===
	
	uint count_patients();
	uint count_hcw();
	bool is_full();
	
	/// Find the position of a patient in '_patient', given its uid.
	uint find_patient_position(uint uid);
	
	// === MOVEMENTS ===
	
	void add_patient(Patient* indiv);
	void add_hcw(HCW* indiv);
	
	void remove_patient(size_t i);
	void remove_hcw(size_t i);
	
	void remove_hcw(HCW* x);
	
	
	// === EPIDEMIOLOGY ===
	
	/// Update the value of _contamIndex, by summing all of its components _contamIdx_amount.
	void update_contamIdx();
	
	/// Increase the environmental contamination index by 'x'. Check <1. Record amount and time.
	void incr_contamIdx(float x, float time);
	
	/// Decrease the environmental contamination index by 'x'. Check >0. 
	void decr_contamIdx(float x) ;
	
	/// Natural decay of contamination.
	void decay_contamIdx(float curr_time);
	
	// === MISCELLENAOUS ===
	
	/// Perform additional cleaning (typically following symptom onset of a patient in that room).
	void extra_cleaning(float clean_room_efficacy);
	
	void show(bool header=true);
	
};





#endif /* Room_hpp */
