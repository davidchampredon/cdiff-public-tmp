//
//  Room.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include <math.h>

#include "Room.hpp"
#include "utils.hpp"


Patient* Room::get_patient_by_uid(uint uid){
	
	Patient* res = nullptr;
	
	uint n = (uint)(_patient.size());
	stopif(n==0, "Cannot get patient: Room is empty!");
	
	bool found = false;
	
	for(uint i=0; i<n; i++){
		if(_patient[i]->get_uid() == uid){
			found = true;
			res = _patient[i];
			break;
		}
	}
	stopif(!found, "Patient UID_" + to_string(uid) + "not found in this room.");
	return res;
}


void Room::add_patient(Patient *indiv){
	string msg =  "Cannot add patient to room #" + to_string(_uid) + ": room full.";
	stopif(_patient.size() == _max_patient,msg);
	_patient.push_back(indiv);
	indiv->set_current_room_uid(_uid);
	indiv->set_current_ward_uid(_ward_uid);
}


void Room::remove_patient(size_t i){
	
	string msg =  "Cannot remove patient from room #" + to_string(_uid) + ": room is empty.";
	stopif(_patient.size() == 0,msg);
	
	msg =  "Cannot remove patient from room #" + to_string(_uid) + ": patient position ("+to_string(i)+") is larger than number of patients in this room ("+to_string(_patient.size())+").";
	stopif(i >= _patient.size(),msg);
	
	_patient.erase(_patient.begin() + i);
}

void Room::remove_hcw(size_t i){
	string msg =  "Cannot remove HCW from room #" + to_string(_uid) + ": room is empty.";
	stopif(_hcw.size() == 0,msg);
	msg =  "Cannot remove HCW from room #" + to_string(_uid) + ": HCW position is larger than number of HCW in this room.";
	stopif(i >= _hcw.size(),msg);
	_hcw.erase(_hcw.begin() + i);
}

void Room::remove_hcw(HCW *x){
	uint uid = x->get_uid();
	
	uint i=0;
	for(i=0; i<_hcw.size(); i++){
		if(_hcw[i]->get_uid() == uid){
			remove_hcw(i);
			break;
		}
	}
	stopif(i==_hcw.size(), "HCW to remove not found!");
}


void Room::show(bool header){
	
	if(header) cout << endl << "-- Room #" << _uid << ":" << endl;
	
	cout << " Type: " << _type << endl;
	cout << " Name: " << _name << endl;
	cout << " Ward: " << _ward_uid << endl;
	cout << " contamIdx: " << _contamIdx << endl;
	cout << " patient size: "<< _patient.size() << " / "<<_max_patient <<endl;
	
	if(_patient.size()>0){
		for(int i=0; i<_patient.size(); i++){
			_patient[i]->show();
		}
	}
	if(_hcw.size()>0){
		for(int i=0; i<_hcw.size(); i++){
			_hcw[i]->show();
		}
	}
	cout << endl << "-- [end info room #"<<_uid<<"]"<<endl;
	
}

uint Room::count_patients(){
	return (uint)_patient.size();
}

uint Room::count_hcw(){
	return (uint)_hcw.size();
}


void Room::update_contamIdx(){
//	_contamIdx = sum(_contamIdx_amount);
	_contamIdx = std::accumulate(_contamIdx_amount.begin(), _contamIdx_amount.end(), 0.00f);
	
	// If the contamination index is tiny,
	// simply set it to 0.0 for speed performance concerns:
	if(_contamIdx < 1e-7){
		_contamIdx = 0.0;
		_contamIdx_amount.clear();
		_contamIdx_time.clear();
	}
}

void Room::incr_contamIdx(float x, float time) {
	// Make sure the contamination index is <= 1.0:
	if(_contamIdx + x > 1) x = 1.0 - _contamIdx;
	// Update amount added for calculating dacay (in other function):
	_contamIdx_amount.push_back(x);
	_contamIdx_time.push_back(time);
	update_contamIdx();
}


void Room::decr_contamIdx(float x){
	float old_contamIndex= _contamIdx;
	
	// Make sure the contamination index is >= 0:
	if(_contamIdx <= x ) {
		_contamIdx = 0.0;
		_contamIdx_amount.clear();
		_contamIdx_time.clear();
	}
	
	if(_contamIdx > x ){
		size_t n = _contamIdx_amount.size();
		
		// NOTE:implicitly, the oldest amounts are
		// removed first. Should not be important (?).
		for(size_t i=0; i<n; i++){
			
			if(_contamIdx_amount[i] <= x){
				x -= _contamIdx_amount[i];
				_contamIdx_amount[i] = 0;
			}
			else{
				_contamIdx_amount[i] -= x;
				x = 0.0;
				i=n; // <-- to force loop exit
			}
		}
		// Reset value to the sum of its components:
		update_contamIdx();
	}

	stopif(old_contamIndex < _contamIdx	,"Problem in function [decr_contamIdx]");
}

void Room::decay_contamIdx(float curr_time){
	
	float decay_rate = log(2.0)/_decay_room_halftime;

	// Calculate new load for each
	// environmental contamination events:
	for(size_t i=0; i<_contamIdx_time.size(); i++){
		float dt = curr_time - _contamIdx_time[i];
		_contamIdx_amount[i] = _contamIdx_amount[i] * exp(-decay_rate * dt);
	}
	// update overall contamination load in this room:
	update_contamIdx();
}



uint Room::find_patient_position(uint uid){
	
	bool found = false;
	uint pos = 0;
	
	for(uint i=0; i<_patient.size(); i++){
		if(_patient[i]->get_uid() == uid){
			found = true;
			pos = i;
		}
	}
	
	string msg = "Patient uid " + to_string(uid) + " not found in room uid " + to_string(_uid);
	stopif(!found,msg);
	
	return pos;
}


bool Room::is_full(){
	bool res = true;
	if(_max_patient > 0 &&
	   count_patients() < _max_patient ) 
		res = false;
	return res;
}


void Room::extra_cleaning(float clean_room_efficacy){
	float oldidx = _contamIdx;
	
	// Draw the efficacy for _this_ cleaning:
	// (formula that parametrizes the Beta distribution with mean)
	float a = 10.0; // For simplicity, shape is hard-coded here.
	float b =  a*(1-clean_room_efficacy)/clean_room_efficacy;
	float cre = beta_distribution(a, b);
	
	// Decrease the contamination index of the room by the drawn mount:
	float decrease_amount = cre * _contamIdx;
	decr_contamIdx(decrease_amount);
	stopif(oldidx < _contamIdx, "Problem in function Room::extra_cleaning()");
}


