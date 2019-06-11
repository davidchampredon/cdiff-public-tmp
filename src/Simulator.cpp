//
//  Simulator.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-15.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "Simulator.hpp"



void Simulator::transfer(uint i, uint j, uint ward, uint room, uint patient_uid){
	
	// Retrieve the patient that will be transfered:
	Patient* p = _HCS[i]->get_room_by_pos(room, ward)->get_patient_by_uid(patient_uid);
	
	// Find available room at destination:
	vector<vector<uint> > idx = _HCS[j]->find_available_room();
	string warn_msg = "WARNING: cannot transfer patient #" + to_string(p->get_uid())+": no beds available at destination." ;
	if(idx.size()==0){
		//cout << warn_msg <<endl;
		return;
	}
	// TO DO: here take the first available, but better to draw randomly!
	uint dest_ward = 0;
	while (idx[dest_ward].size()==0 & dest_ward<idx.size()) {
		dest_ward++;
	}
	if(dest_ward == idx.size()){
//		cout << warn_msg << endl;
		return;
	}
	
	uint dest_room = idx[dest_ward][0];
	Room* r = _HCS[j]->get_room_by_pos(dest_room, dest_ward);
	
	// Send the patient in the destination HCS.
	// WARNING: assume susceptibility does not change during transfer.
	_HCS[j]->receive_transfered_patient(p,
										_HCS[j]->draw_DoS(),
										_HCS[j]->get_contamIdx_suscept_prm(),
										_HCS[j]->get_decay_contamIdx_halftime(),
										r,
										_HCS[j]->get_current_time());
	
	// Now that the patient is added in
	// the destination HCS, we can safely remove it
	// from the origin HCS:
	_HCS[i]->discharge_patient_by_uid(ward, room, patient_uid);
	
	// DEBUG:
	//cout << "Patient #"<< p->get_uid() << " has been transfered from HCS#"<< i;
	//cout << " ward #"<<ward<<" room UID#"<< _HCS[i]->get_room_by_pos(room, ward)->get_uid() ;
	//cout << " to HCS#"<<j<< " ward #"<<dest_ward<< " room UID#"<< r->get_uid() << endl;
}


void Simulator::transfer_all_patients(float t, int trace_level){
	
	// How many transfer between HCS:
	vector<vector<uint> > n_tr   = draw_n_transfer();
	
	//DEBUG:
	if(trace_level>0){
		cout<< " number of transfers drawn:";
		displayVector(n_tr);
	}
	
	for(uint i=0; i<_HCS.size(); i++)
	{
		vector<vector<uint> > tr_pat = _HCS[i]->find_transferable_patients();
		
		// if there is at least one transferable patient,
		// loop all possible destinations
		if(tr_pat.size()>0){
			for(uint m=0; m < n_tr[i].size(); m++)
			{
				if(m != i){ // don't look at the same HCS
					
					//DEBUG:
					//cout << "trying to transfer "<< n_tr[i][m] <<" patients from HCS"<<i<< " to HCS"<<m <<endl;
					
					uint k   = 0;
					uint cnt = 0;
					while( (k < tr_pat.size()) &&
						  (cnt < n_tr[i][m]) ){
						
						// Retrieve the room and position
						// of one transferable patient:
						uint ward = tr_pat[k][0];
						uint room = tr_pat[k][1];
						uint puid = tr_pat[k][2];
						
						// transfer the identified patient:
						transfer(i,m,ward, room, puid);
						k++;
					}
				}
			}
		}
	}
}

void Simulator::events_one_HCS(uint i, float t, int trace_level){
	
	// Movements in & out:
	_HCS[i]->discharge_all();
	_HCS[i]->admission_new_all();
	
	// Update the infectivity index of all
	// vaccinated individuals. This reflects
	// vaccine waning.
	_HCS[i]->update_all_infectIdx_vax();
	
	// Clean rooms that had symptomatic patients:
	_HCS[i]->extra_cleaning();
	
	// Interaction between HCW & patients
	_HCS[i]->draw_all_contacts_patient_HCW(_recordContacts);
	
	// Patients shedding and/or acquiring from their own room:
	_HCS[i]->all_patients_envContam_own_room();
	
	// Acquisition from 'other' environment (not in patient room)
	_HCS[i]->draw_acq_envContam_other();
	
	// Colectomy events:
	_HCS[i]->draw_colectomy();
	
	// Fecal microbial transplant events:
	_HCS[i]->draw_fmt();
	
	// Updates
	_HCS[i]->adjust_DoS();
	_HCS[i]->update_all_clinical();
	_HCS[i]->update_decay_contamIdx_patient_rooms();
	_HCS[i]->update_envContam_other();
	_HCS[i]->update_isoDur();
	_HCS[i]->update_totDurSymptoms();
	
	// Isolation of patients starting symptoms:
	_HCS[i]->draw_onset_iso_lag();
	_HCS[i]->start_isolation();
	_HCS[i]->end_isolation();
	
	// Relapses:
	_HCS[i]->relapses();
	
	// Reset the contamination index level to
	// a small value for all HCWs:
	// TO DO: that may be wrong when HCW schedule overlaps 2 days!!!
	if(abs(t - (int)(t)) < 0.01)
		_HCS[i]->reset_contamIdx_all_HCW();
	
	// Census snapshot:
	censusDF df(*_HCS[i], t);
	censusDF df_light(*_HCS[i], t, true);
	_census.push_back(df);
	_census_light.push_back(df_light);
}


void Simulator::events_all_HCS(float t, int trace_level){
	for(uint i=0; i<_HCS.size(); i++){
		events_one_HCS(i, t, trace_level);
	}
}


void Simulator::set_current_time(float t){
	for(size_t i=0; i<_HCS.size(); i++){
		// Update simulation time:
		_HCS[i]->set_current_time(t);
	}
}

void Simulator::run_simulation(uint seed){
	
	
	cout << "Entering run_simulation(seed="<<seed<<")" <<endl;
	
	int trace_level = 0;
	
	float hoursWorkedInOneDay = 8.0; // TO DO: change that to input

	// Integrity checks:
	// TO DO : put in function
	stopif(_HCS.size() != _mvt_HCS.size(),
		   "Number of HCS and movement definition inconsistent!" );
	
	// Set timestep:
	for(uint i=0; i<_HCS.size(); i++){
		_HCS[i]->set_time_step(_timeStep);
		assign_HCW_periodsWorked(i, hoursWorkedInOneDay);
	}
	
	// Main time loop
	float t;
	for(t = 0.0 ; t <_horizon; t += _timeStep){
		
		set_current_time(t);
		
		if((int)(t)%100 == 0 && (t-(int)(t)<0.0001) ){
			cout<<"seed_"<< seed << " simulation time:"<<t<< "/"<< _horizon <<endl;
		}
		
		if(trace_level>1){
			cout << " ; period: " << period_of_day(t);
			cout << "/"<< total_periods_in_day() - 1 ;
			cout << endl <<endl;
		}
		
		// Transfers:
		if(_HCS.size()>1) transfer_all_patients(t, trace_level);
		
		// Events within each HCS
		events_all_HCS(t, trace_level);
		
	} // end for 't'
	cout << endl << " ---> Simulation ended with no error <---" << endl;
}


uint Simulator::total_periods_in_day(){
	uint res = round(1/_timeStep);
	stopif(res<1,"Simulation time step is too large. Must be < 1.");
	return res;
}


uint Simulator::period_of_day(float t){
	float d	= t - (uint)(t);
	uint n	= round(d / _timeStep);
	return n;
}


void Simulator::assign_HCW_periodsWorked(uint i, float hoursWorkedInOneDay){
	
	// How many periods to make the fraction of day worked:
	float frac	= hoursWorkedInOneDay / 24.0;
	uint np		= round(frac / _timeStep);
	
	// All periods of a day:
	uint N = total_periods_in_day();
	vector<uint> x;
	for(uint i=0; i<N; i++) x.push_back(i);

	uint n_ward		= (uint)_HCS[i]->get_ward().size();
	uint n_freehcw	= (uint)_HCS[i]->get_hcw().size();
	
	// --- Free HCWs:
	
	for(uint k=0; k < n_freehcw; k++){
		
		// Randomly select (without replacement) 'np' periods in a day:
		std::shuffle(std::begin(x), std::end(x), RANDOM_GENERATOR);
		vector<uint> y(np);
		for(uint i=0; i<np; i++) y[i] = x[i];
		
		_HCS[i]->get_hcw()[k]->set_periodsWorked(y);
	}
	
	// --- HCWs assigned to a ward:
	for(uint w=0; w < n_ward; w++)
	{
		Ward* WARD = _HCS[i]->get_ward()[w];
		for(uint k=0; k < WARD->get_hcw().size(); k++)
		{
			// Randomly select (without replacement) 'np' periods in a day:
			std::shuffle(std::begin(x), std::end(x), RANDOM_GENERATOR);
			vector<uint> y(np);
			for(uint i=0; i<np; i++) y[i] = x[i];
			WARD->get_hcw(k)->set_periodsWorked(y);
		}
	}
}



vector<vector<uint> > Simulator::draw_n_transfer(){
	
	vector<vector<uint> > res(_mvt_HCS.size());
	
	for(uint i=0; i<_mvt_HCS.size(); i++){
		res[i].resize(_mvt_HCS[i].size());
		for(uint j=0; j<_mvt_HCS[i].size(); j++){
			poisson_distribution<uint> rpois(_mvt_HCS[i][j]);
			res[i][j] = rpois(RANDOM_GENERATOR);
		}
	}
	return res;
}







