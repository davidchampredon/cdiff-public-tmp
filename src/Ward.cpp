//
//  Ward.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-27.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "Ward.hpp"



bool Ward::is_full(){
	
	bool res = true;
	
	for(uint i=0; i<_room.size(); i++){
		if( ! _room[i]->is_full() ){
			res = false;
			break;
		}
	}
	return res;
}


Room* Ward::find_available_room(){
	for(uint i=0; i<_room.size(); i++){
		if(! _room[i]->is_full() ) return _room[i];
	}
	return nullptr;
}
