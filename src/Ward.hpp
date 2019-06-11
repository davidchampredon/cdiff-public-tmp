//
//  Ward.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-27.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef Ward_hpp
#define Ward_hpp

#include <stdio.h>
#include "Room.hpp"
#include "Individual.hpp"


class Ward{
	
	vector<Room*> _room;
	
	vector<HCW*>  _hcw;
	
public:
	
	Ward(){_room.clear();};
	
	Ward(vector<Room*> r, vector<HCW*> h){
		for(uint i=0; i<r.size(); i++) _room.push_back(r[i]);
		for(uint i=0; i<h.size(); i++) _hcw.push_back(h[i]);
	}
	
	vector<Room*>	get_room(){return _room;}
	
	/// Return the room in position 'i' in the '_room' vector for this ward.
	Room*			get_room(uint i){return _room[i];}
	
	vector<HCW*>	get_hcw(){return _hcw;}
	HCW*			get_hcw(uint i){return _hcw[i];}
	

	void add_room(Room* r){
		_room.push_back(r);
	}
	
	void add_hcw(HCW* h){
		_hcw.push_back(h);
	}
	
	/// Returns true if all rooms have the maximum number of patients (all beds occupied).
	bool is_full();
	
	/// Find the first room that has an available bed.
	Room* find_available_room();
	
	
};



#endif /* Ward_hpp */
