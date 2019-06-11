//
//  Registry.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-15.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "Registry.hpp"



void Registry::delete_record(uint uid_patient){
	
	uint i_found = (uint)_uid_patient.size() + 1;
	for(uint i=0; i<_uid_patient.size(); i++){
		if(_uid_patient[i] == uid_patient){
			i_found = i;
			break;
		}
	}
	stopif(i_found >= _uid_patient.size(),
		   "uid_patient not found in Registry: cannot delete.");
	
	_uid_patient.erase(_uid_patient.begin() + i_found);
	_uid_room.erase(_uid_room.begin() + i_found);
	_uid_ward.erase(_uid_ward.begin() + i_found);
}



void Registry::show(){
	cout << "Patient UID <---> Room UID" <<endl;
	for(uint i=0; i<_uid_patient.size(); i++){
		cout << _uid_patient[i] << " <---> "<<_uid_room[i] << endl;
	}
}
