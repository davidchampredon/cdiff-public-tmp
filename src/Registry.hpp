//
//  Registry.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-15.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef Registry_hpp
#define Registry_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

#include "utils.hpp"

class Registry{
private:
	std::vector<uint> _uid_patient;
	std::vector<uint> _uid_room;
	std::vector<uint> _uid_ward;
public:
	Registry(){
		_uid_room.clear();
		_uid_patient.clear();
		_uid_ward.clear();
	}
	
	void add_record(uint uid_patient, uint uid_room, uint uid_ward){
		_uid_patient.push_back(uid_patient);
		_uid_room.push_back(uid_room);
		_uid_ward.push_back(uid_ward);
	}
	
	void delete_record(uint uid_patient);
	
	
	uint get_size(){return (uint)_uid_patient.size();}
	
	uint get_uid_patient(uint i){return _uid_patient[i];}
	uint get_uid_room(uint i){return _uid_room[i];}
	uint get_uid_ward(uint i){return _uid_ward[i];}
	
	
	void show();
	
};


#endif /* Registry_hpp */
