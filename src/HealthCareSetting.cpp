//
//  HealthCareSetting.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "HealthCareSetting.hpp"
#include "utils.hpp"
#include "globalvar.hpp"
#include "transmission.hpp"


void HealthCareSetting::show(bool details){
	
	cout << endl << endl << " ===== HealthCareSetting Info =====" << endl <<endl;
	
	int g = 20;
	tabcout("Name", _name, g);
	tabcout("Type", _type, g);
	tabcout("Number of wards", _ward.size(),g);
	tabcout("Number of patients", count_patients(),g);
	tabcout("Patients capacity", sum(max_patient_per_ward()),g);
	tabcout("Number of HCW w/ ward:", count_hcw_ward(),g);
	tabcout("Number of free HCW:", count_hcw_free(),g);
	
	if(details){
		_registry.show();
		
		cout <<endl << "---- Free HCWs ----" << endl;
		for(uint i=0; i<_hcw.size(); i++){
			cout <<endl << "---- HCW "<<i << endl;
			_hcw[i]->show();
		}
		
		for(uint w=0; w<_ward.size();w++){
			cout << endl << " === WARD "<< w << " === " <<endl;
			for(uint i=0; i<_ward[w]->get_room().size(); i++){
				cout << endl << "---- ROOM "<< i << endl;
				_ward[w]->get_room(i)->show();
			}
			for(uint i=0; i<_ward[w]->get_hcw().size(); i++){
				cout << endl << "---- HCW "<<i << endl;
				_ward[w]->get_hcw()[i]->show();
			}
		}
	}
	cout << endl << endl << " =====[end HealthCareSetting Info]====" << endl <<endl;
}


vector<uint> HealthCareSetting::max_patient_per_ward(){
	
	uint n = (uint)_ward.size();
	
	// Count the max number of patients for each ward:
	vector<uint> res(n, 0);
	
	for(uint w=0; w<n; w++){
		for(uint i=0; i<_ward[w]->get_room().size(); i++){
			uint mxp  = _ward[w]->get_room(i)->get_max_patient();
			res[w] += mxp;
		}
	}
	return(res);
}


void HealthCareSetting::setup(string hcs_name,
							  string hcs_type,
							  uint n_iso_rooms,
							  uint n_ward_mean,
							  uint n_room_per_ward_mean,
							  uint n_beds_per_room_mean,							  
//							  vector<string> room_type,
//							  vector<uint>   ward,
//							  vector<string> room_name,
//							  vector<uint>   max_patient,
							  uint  nHCW_ward,
							  uint  nHCW_free,
							  uint  uid_start,
							  float occupancy ,
							  float mean_sIdx,
							  float sd_sIdx,
							  float mean_sIdx_HCW,
							  float sd_sIdx_HCW,
							  float mDoS,
							  float sdDoS,
							  float admission_rate,
							  vector<float> lambda_prm,
							  float decay_room_halftime,
							  float contamIdx_suscept_prm,
							  float decay_contamIdx_halftime,
							  float ratio_contam_other_room,
							  float prev_new_admission,
							  float proba_patient_visit_other,
							  float reset_contamIdx,
							  float mean_contact_HCW_patient_duration_minutes,
							  float var_contact_HCW_patient_duration_minutes,
							  float contact_indiv_halftime,
							  float contact_xchg_contamIdx,
							  float visit_xchg_contamIdx,
							  float visit_room_halftime,
							  float clean_room_efficacy,
							  float onset_isolation_lag_mean,
							  Vaccine vaxObj,
							  string vaxStrategy,
							  float proba_vax_newAdmission,
							  float sIdxMin_vax_frailty,
							  float DoSMin_vax_frailty,
							  float vax_start_time,
							  float proba_colectomy,
							  float proba_fmt){
	
	// Create the infrastructure (e.g., rooms):

	//	build_HCS(hcs_name,  hcs_type,
//			  room_type, ward,  room_name,
//			  max_patient,
//			  decay_room_halftime);
	
	build_HCS_random(hcs_name,
					 hcs_type,
					 n_iso_rooms,
					 n_ward_mean,
					 n_room_per_ward_mean,
					 n_beds_per_room_mean,
					 decay_room_halftime);
	
	cout << "build_HCS() done." <<endl;
	show();
	
	// Set parameters values:
	_contamIdx_suscept_prm	= contamIdx_suscept_prm;
	_decay_contamIdx_halftime	= decay_contamIdx_halftime;
	_ratio_contam_other_room = ratio_contam_other_room;
	_prev_new_admission = prev_new_admission;
	_proba_patient_visit_other = proba_patient_visit_other;
	_reset_contamIdx = reset_contamIdx;
	_mean_contact_HCW_patient_duration_minutes = mean_contact_HCW_patient_duration_minutes;
	_var_contact_HCW_patient_duration_minutes = var_contact_HCW_patient_duration_minutes;
	_contact_indiv_halftime = contact_indiv_halftime;
	_contact_xchg_contamIdx = contact_xchg_contamIdx;
	_visit_xchg_contamIdx = visit_xchg_contamIdx;
	_visit_room_halftime = visit_room_halftime;
	_clean_room_efficacy = clean_room_efficacy;
	_onset_isolation_lag_mean = onset_isolation_lag_mean;
	
	_proba_colectomy = proba_colectomy;
	_proba_fmt = proba_fmt;
	
	cout << "setup()-params." <<endl;
	
	
	// Vaccination:
	set_vaccine(vaxObj);
	set_vaxStrategy(vaxStrategy);
	set_proba_vax_newAdmission(proba_vax_newAdmission);
	set_sIdxMin_vax_frailty(sIdxMin_vax_frailty);
	set_DoSMin_vax_frailty(DoSMin_vax_frailty);
	set_vax_start_time(vax_start_time);
	
	cout << "setup()-vaccine." <<endl;
	
	// Populate with patients:
	set_max_uid_indiv(uid_start); // TO DO: do a check to prevent same uid in 2 different HCS
	set_sIdx_newPatient_prm({mean_sIdx, sd_sIdx});
	set_DoS_newPatient_prm({mDoS, sdDoS});
	populate_patient(occupancy);
	cout << "setup()-individuals" <<endl;
	
	// Populate with HCWs:
	set_sIdx_HCW_prm({mean_sIdx_HCW, sd_sIdx_HCW});
	populate_hcw(nHCW_ward, nHCW_free);
	
	cout << "setup()-populate." <<endl;
	// Movements:
	set_admission_rate(admission_rate);

	// Contacts:
	set_lambda_prm(lambda_prm);
	
	cout << "setup()-contacts." <<endl;
	
	seed_initial_infection(prev_new_admission);
	
	cout << "setup()-seedInitialInfections." <<endl;
}


void HealthCareSetting::build_HCS(string hcs_name,
								  string hcs_type,
								  vector<string> room_type,
								  vector<uint> ward,
								  vector<string> room_name,
								  vector<uint> max_patient,
								  float decay_room_halftime){
	
	_name = hcs_name;
	_type = hcs_type;
	
	_current_time = -NA_FLOAT;
	_time_step    = 0.0;
	
	size_t n = room_type.size();
	
	// Integrity checks:
	string msg = "build_HCS(): All vectors must have the same size, but at least one does not.";
	stopif(n != ward.size() ||
		   n != room_name.size() ||
		   n != max_patient.size(),
		   msg);
	
	// Total number of wards:
	uint mxw = 0;
	for(uint i=0; i<n; i++)
		if(mxw < ward[i]) mxw = ward[i];
	
	// Create Ward objects:
	_ward.clear();
	_ward.resize(mxw+1);
	for(uint i=0; i<mxw+1; i++){
		Ward* tmpw = new Ward;
		_ward[i] = tmpw;
	}
	
	// Populate wards with rooms:
	for(uint i=0; i<n; i++){
		Room *tmp = new Room;
		
		tmp->set_type(room_type[i]);
		tmp->set_uid(i);
		tmp->set_ward_uid(ward[i]);
		tmp->set_name(room_name[i]);
		tmp->set_max_patient(max_patient[i]);
		tmp->set_decay_room_halftime(decay_room_halftime);
		
		_ward[ward[i]]->add_room(tmp);
	}
	
	// Start with a clean environment:
	_envContam_other = 0.0;
	
	_contact_HCW_patient_dur_rec.clear();
}

void HealthCareSetting::build_HCS_random(string hcs_name,
										 string hcs_type,
										 uint n_iso_rooms,
										 uint n_ward_mean,
										 uint n_room_per_ward_mean,
										 uint n_beds_per_room_mean,
										 float decay_room_halftime){
	
	_name = hcs_name;
	_type = hcs_type;
	
	_current_time = -NA_FLOAT;
	_time_step    = 0.0;
	// Start with a clean environment:
	_envContam_other = 0.0;
	_contact_HCW_patient_dur_rec.clear();
	
	// setup distributions to draw random numbers:
	poisson_distribution<uint> rpois_wards(n_ward_mean);
	poisson_distribution<uint> rpois_rooms(n_room_per_ward_mean);
	poisson_distribution<uint> rpois_beds(n_beds_per_room_mean);
	
	// Create Ward objects:
	uint n_wards = rpois_wards(RANDOM_GENERATOR);
	_ward.clear();
	_ward.resize(n_wards+1);
	for(uint i=0; i<n_wards+1; i++){
		Ward* tmpw = new Ward;
		_ward[i] = tmpw;
	}
	
	// Ward #0 is for isolation rooms:
	for(uint i=0; i<n_iso_rooms; i++){
		Room *tmp = new Room;
		
		tmp->set_type("isolation_room");
		tmp->set_uid(i);
		tmp->set_ward_uid(0);
		tmp->set_name("ISO"+to_string(i));
		tmp->set_max_patient(1);
		tmp->set_decay_room_halftime(decay_room_halftime);
		_ward[0]->add_room(tmp);
	}
	
	// Add rooms and beds to each ward:
	for(uint w=1; w <= n_wards; w++){
		uint n_rooms_w = rpois_rooms(RANDOM_GENERATOR);
		for(uint r=0; r<n_rooms_w; r++){
			Room *tmp = new Room;
			tmp->set_type("patient_room");
			tmp->set_uid(r);
			tmp->set_ward_uid(w);
			tmp->set_name("W" + to_string(w) + "_R" + to_string(r));
			// draw number of beds in this room:
			uint n_beds_w_r = rpois_beds(RANDOM_GENERATOR);
			tmp->set_max_patient(n_beds_w_r);
			
			tmp->set_decay_room_halftime(decay_room_halftime);
			_ward[w]->add_room(tmp);
		}
	}
}

void HealthCareSetting::stats_architecture(){
	
	size_t n_wards = _ward.size()-1;
	size_t n_iso = _ward[0]->get_room().size();
	
	uint s_rooms = 0;
	uint s_beds  = 0;
	
	for(uint w=1; w <= n_wards; w++){
		size_t nr = _ward[w]->get_room().size();
		s_rooms += nr;
		for(uint r =0; r<nr ; r++){
			s_beds += _ward[w]->get_room(r)->get_max_patient();
		}
	}
	
	cout << endl << " === Stats Architecture ========================== " << endl;
	tabcout("Number of wards (non-iso)", n_wards);
	tabcout("Number isolation rooms", n_iso);
	tabcout("Number patient rooms", s_rooms);
	tabcout("Number beds", s_beds);
	tabcout("Mean rooms per ward", (float)(s_rooms) / n_wards);
	tabcout("Mean beds per room", (float)(s_beds) / s_rooms);
	cout << " ================================================= " << endl;
}


void HealthCareSetting::admission_new_patient(float sIdx,
											  float DoS,
											  float contamIdx_suscept_prm,
											  float decay_contamIdx_halftime,
											  Room *r,
											  float time){
	Patient* p = new Patient;
	p->set_uid(generate_new_uid_indiv());
	p->set_susceptIdx(sIdx);
	p->set_susceptIdxBaseline(sIdx);
	
	// Assume the proba of infection is
	// the same as colonization proba:
	float infectIdx = _disease.get_infectIdx_init();
	p->set_infectIdx(infectIdx);
	p->set_infectIdx_baseline(infectIdx);
	
	p->set_timeAdmission(time);
	p->set_DoS(DoS);
	p->set_contamIdx_suscept_prm(contamIdx_suscept_prm);
	p->set_decay_contamIdx_halftime(decay_contamIdx_halftime);
	
	// Draw if this new patient has been
	// vaccinated prior admission:
	vaccinate_individual(p);
	
	// Draw if this new patient is colonized
	// at admission or not (conditional on being vaccinated):
	bernoulli_distribution dbern(_prev_new_admission);
	bool isCol = dbern(RANDOM_GENERATOR);
	if(isCol){
		float past_time = time - 1.0; // Gross assumption... TO DO: change this to random times
		if(past_time<0) past_time = 0.0;
		p->colonize(past_time, _disease);
		p->set_acquisitionType("community");
	}
	
	r->add_patient(p);
	
	// DEBUG:
	//cout <<"New patient admitted uid:"<< p->get_uid()<<" ; in room "<< r->get_uid();
	//cout << " ; colonized@admin: "<< isCol << " ; contamIdx = " << p->get_contamIdx() << endl;
	
	// Record in the registry as admission:
	_registry.add_record(p->get_uid(), r->get_uid(), r->get_ward_uid() );
}


float HealthCareSetting::draw_DoS(){
	
	lognormal_distribution<float> dLnorm(_DoS_newPatient_prm[0],
										 _DoS_newPatient_prm[1]);
	float dos = dLnorm(RANDOM_GENERATOR);
	
	// Put a hard cap on the value.
	// It may have unrealistic values by chance.
	if(dos>1000) dos = 1000.0;
	
	return dos;
}


float HealthCareSetting::draw_sIdx(string indiv_type){
	
	double m = 0, s=0;
	bool type_found = false;
	
	if(indiv_type == "patient") {
		m = _sIdx_newPatient_prm[0];
		s = _sIdx_newPatient_prm[1];
		type_found = true;
	}
	else if (indiv_type == "HCW"){
		m = _sIdx_HCW_prm[0];
		s = _sIdx_HCW_prm[1];
		type_found = true;
	}
	
	stopif(!type_found,
		   "Unknown individual type for draw_sIdx(string indiv_type).");
	
	normal_distribution<float> dNorm(m,s);
	float sIdx = dNorm(RANDOM_GENERATOR);
	sIdx = (sIdx<1e-4 ? 1e-4 : sIdx);
	sIdx = (sIdx>1 ? 1 : sIdx);
	
	// Bimodal distribution for patients:
	if(indiv_type == "patient") {
		// Multiplicative factor for mean
		// susceptibility index of low risk group
		// when compared to high risk group:
		float mult_low_risk_group = 1.0 / 4.0;
		
		// Draw normal distribution:
		normal_distribution<float> dNorm2(m * mult_low_risk_group, s);
		float sIdx2 = dNorm2(RANDOM_GENERATOR);
		sIdx2 = (sIdx2<0 ? 0 : sIdx2);
		sIdx2 = (sIdx2>1 ? 1 : sIdx2);
		
		bernoulli_distribution dBern(0.5);
		int pp = dBern(RANDOM_GENERATOR);
		if(pp) sIdx=sIdx2;
	}
	
	return sIdx;
}


void HealthCareSetting::populate_patient(float occupancy){
	
	float t0 = 0.0;
	
	size_t n_w =_ward.size();
	
	for(uint w=0; w<n_w; w++){
		size_t n_r =_ward[w]->get_room().size();
		for(uint i=0; i<n_r; i++)
		{
			string	type	= _ward[w]->get_room(i)->get_type();
			uint	mxp		= _ward[w]->get_room(i)->get_max_patient();
			
			if( type == "patient_room" && mxp >0){
				
				// Draw the number of patient in this room:
				binomial_distribution<uint> dBinom(mxp, occupancy);
				uint n = dBinom(RANDOM_GENERATOR);
				
				// For each patient, formally admit it:
				for(uint k = 0; k<n; k++){
					float sIdx = draw_sIdx();
					float DoS  = draw_DoS();
					admission_new_patient(sIdx, DoS,
										  _contamIdx_suscept_prm,
										  _decay_contamIdx_halftime,
										  _ward[w]->get_room(i), t0);
				}
			}
		}
	}
}


void HealthCareSetting::populate_hcw(uint nHCW_ward,
									 uint nHCW_free){
	
	stopif(nHCW_ward < _ward.size(), "MORE WARDS THAN HCW SUPPLIED!");
	
	// The allocation of HCWs to wards
	// is proportional to the size of the
	// ward, defined as the maximum number
	// patients in the ward.
	
	vector<uint> mxp = max_patient_per_ward();
	
	// Maximum number of patients in the whole HCS:
	uint tot = sum(mxp);
	
	vector<uint> hcw_per_ward;
	for(uint i=0; i<mxp.size(); i++){
		float prop_patients_i = (float)mxp[i]/(float)(tot) ;
		uint npw = uint( prop_patients_i * (float)nHCW_ward );
		hcw_per_ward.push_back(npw);
	}
	
	for(uint w=0; w<_ward.size(); w++){
		uint count = 0;
		while(count <= hcw_per_ward[w]){
			HCW* hcw = new HCW;
			hcw->set_uid(generate_new_uid_indiv());
			hcw->set_type("HCW");
			hcw->set_subtype("allocWard");
			hcw->set_susceptIdx(draw_sIdx("HCW"));
			hcw->set_wardAssigned(w);
			hcw->set_contamIdx_suscept_prm(_contamIdx_suscept_prm);
			hcw->set_contamIdx(0.0);
			hcw->set_infectIdx(0.0); // assume HCW cannot be infected.
			hcw->set_infectIdx_baseline(0.0);
			
			_ward[w]->add_hcw(hcw);
			count++;
		}
	}
	
	// "Free" HCW, not allocated a ward
	
	for(uint i=0; i < nHCW_free; i++){
		HCW* hcw = new HCW;
		hcw->set_uid(generate_new_uid_indiv());
		hcw->set_type("HCW");
		hcw->set_subtype("free");
		hcw->set_susceptIdx(draw_sIdx("HCW"));
		hcw->set_contamIdx_suscept_prm(_contamIdx_suscept_prm);
		hcw->set_contamIdx(0.0);
		hcw->set_infectIdx(0.0); // assume HCW cannot be infected.
		hcw->set_infectIdx_baseline(0.0);
		
		add_hcw(hcw);
	}
}

uint HealthCareSetting::count_patients(){
	uint s = 0;
	for(uint w=0; w<_ward.size();w++){
		for(uint i=0; i<_ward[w]->get_room().size(); i++){
			s += _ward[w]->get_room(i)->count_patients();
		}
	}
	return s;
}

uint HealthCareSetting::count_hcw_ward(){
	uint s = 0;
	for(uint w=0; w<_ward.size();w++)
		s += _ward[w]->get_hcw().size();
	
	return s;
}

uint HealthCareSetting::count_hcw_free(){
	return (uint)_hcw.size();
}


void HealthCareSetting::census_one_individual(Individual *indiv){
	
	census_uid_individual.push_back(indiv->get_uid());
	census_sIdx.push_back(indiv->get_susceptIdx());
	census_infectIdx.push_back(indiv->get_infectIdx());
	census_type_individual.push_back(indiv->get_type());
	census_subtype_individual.push_back(indiv->get_subtype());
	census_isColonized.push_back(indiv->get_isColonized());
	census_wasColonized.push_back(indiv->get_wasColonized());
	census_contamIdx.push_back(indiv->get_contamIdx());
	census_contamIdx_other.push_back(_envContam_other);
	census_current_room_uid.push_back(indiv->get_current_room_uid());
	census_hcs_name.push_back(_name);
	
	census_timeAcquisition.push_back(indiv->get_timeAcquisition());
	census_timeClearance.push_back(indiv->time_clearance());
	
	census_nColonizations.push_back(indiv->get_nColonizations());
	
	census_acquisitionType.push_back(indiv->get_acquisitionType());
}


void HealthCareSetting::census_one_patient(Patient *p){
	
	census_timeAdmission.push_back(p->get_timeAdmission());
	census_DoS.push_back(p->get_DoS());
	
	census_isAlive.push_back(p->get_isAlive());
	census_willDie.push_back(p->get_willDie());
	
	census_colectomy.push_back(p->get_colectomy());
	census_fmt.push_back(p->get_fmt());
	
	census_wardAssigned.push_back(p->get_current_ward_uid());
	
	census_symptomatic.push_back(p->get_symptomatic());
	census_currentlySymptomatic.push_back(p->get_currentlySymptomatic());
	census_prdIncubation.push_back(p->get_prdIncubation());
	census_prdSymptom.push_back(p->get_prdSymptom());
	
	census_nInfections.push_back(p->get_nInfections());
	census_nRelapses.push_back(p->get_nRelapses());
	
	census_onDuty.push_back(false); // Patients are never "on duty", only HCWs can.
	
	census_isVax.push_back(p->get_isVax());
	
	census_timeOnset.push_back(p->get_timeOnset());
	census_timeRelapse.push_back(p->get_timeRelapse());
	census_timeIsolationStart.push_back(p->get_timeIsolationStart());
	census_timeIsolationEnd.push_back(p->get_timeIsolationEnd());
	
	census_isoDurFirst.push_back(p->get_isoDurFirst());
	census_isoDurRelapses.push_back(p->get_isoDurRelapses());
	
	census_totDurSymptoms.push_back(p->get_totDurSymptoms());
}


void HealthCareSetting::census_one_hcw(HCW *hcw){
	
	census_timeAdmission.push_back(0.0);
	census_DoS.push_back(NA_FLOAT);
	
	census_isAlive.push_back(true);
	census_willDie.push_back(false);
	census_colectomy.push_back(false);
	census_fmt.push_back(false);
	
	census_wardAssigned.push_back(hcw->get_wardAssigned());
	
	// assumption: HCW never symptomatic:
	census_symptomatic.push_back(false);
	census_currentlySymptomatic.push_back(false);
	census_prdIncubation.push_back(0.0);
	census_prdSymptom.push_back(0.0);
	
	census_onDuty.push_back(hcw->onDuty(_current_time, _time_step));
	
	census_isVax.push_back(false);
	
	census_timeOnset.push_back(NA_FLOAT);
	census_timeIsolationStart.push_back(-NA_FLOAT);
	census_timeIsolationEnd.push_back(NA_FLOAT);
	census_timeRelapse.push_back(NA_FLOAT);
	
	census_isoDurFirst.push_back(NA_FLOAT);
	census_isoDurRelapses.push_back(NA_FLOAT);
	census_totDurSymptoms.push_back(NA_FLOAT);
}

void HealthCareSetting::census_clear_all(){
	census_uid_individual.clear();
	census_isAlive.clear();
	census_willDie.clear();
	census_colectomy.clear();
	census_fmt.clear();
	census_sIdx.clear();
	census_infectIdx.clear();
	census_type_individual.clear();
	census_subtype_individual.clear();
	census_isColonized.clear();
	census_wasColonized.clear();
	census_contamIdx.clear();
	census_contamIdx_other.clear();
	census_current_room_uid.clear();
	census_timeAdmission.clear();
	census_DoS.clear();
	census_wardAssigned.clear();
	census_symptomatic.clear();
	census_currentlySymptomatic.clear();
	census_prdIncubation.clear();
	census_prdSymptom.clear();
	census_onDuty.clear();
	census_isVax.clear();
	
	census_timeAcquisition.clear();
	census_timeOnset.clear();
	census_timeIsolationStart.clear();
	census_timeIsolationEnd.clear();
	census_timeClearance.clear();
	census_timeRelapse.clear();
	
	census_isoDurFirst.clear();
	census_isoDurRelapses.clear();
	census_totDurSymptoms.clear();
	
	census_acquisitionType.clear();
	
	census_nColonizations.clear();
	census_nInfections.clear();
	census_nRelapses.clear();
}


void HealthCareSetting::census_patients_light(){
	
	// Reset everything:
	census_clear_all();
	
	// Scan all rooms in this HCS
	// and retrieve information:
	for(uint w=0; w<_ward.size(); w++)
	{
		// Patients only:
		for(uint i=0; i<_ward[w]->get_room().size(); i++){
			for(uint j=0; j< _ward[w]->get_room(i)->count_patients(); j++){
				Patient* p = _ward[w]->get_room(i)->get_patient(j);
				
				// Record only if admission or discharge time:
//				float time_adm = p->get_timeAdmission();
				float time_dis = p->get_timeDischarge();
				
//				bool cond_adm = is_in_current_timestep(time_adm);
				bool cond_dis = is_in_current_timestep(time_dis);

//				bool cond_time = (cond_adm || cond_dis);
				bool cond_time = cond_dis;
				
				if( cond_time ){
					census_one_individual(p);
					census_one_patient(p);
					// room contamination
					// (info is redundant [several times for each patient]
					//  because a tall-format dataframe is built)
					float contamidx = _ward[w]->get_room(i)->get_contamIdx();
					census_current_room_contamIdx.push_back(contamidx);
				}
			}
		}
	}
}


void HealthCareSetting::census_individuals(){
	
	// Reset everything:
	census_clear_all();
	
	// Scan all rooms in this HCS
	// and retrieve information:
	for(uint w=0; w<_ward.size(); w++)
	{
		// Patients:
		for(uint i=0; i<_ward[w]->get_room().size(); i++){
			for(uint j=0; j< _ward[w]->get_room(i)->count_patients(); j++){
				Patient* p = _ward[w]->get_room(i)->get_patient(j);
				census_one_individual(p);
				census_one_patient(p);
				// room contamination
				// (info is redundant [several times for each patient]
				//  because a tall-format dataframe is built)
				census_current_room_contamIdx.push_back(_ward[w]->get_room(i)->get_contamIdx());
			}
		}
		// HCW assigned to a ward:
		for(uint j=0; j<_ward[w]->get_hcw().size(); j++){
			HCW* hcw = _ward[w]->get_hcw(j);
			census_one_individual(hcw);
			census_one_hcw(hcw);
			census_current_room_contamIdx.push_back(-NA_FLOAT); // no room is allocated to HCW, but need a dummy value for tall-format dataframe
		}
	}
	// Free HCWs:
	for(uint i=0; i<_hcw.size(); i++){
		census_one_individual(_hcw[i]);
		census_one_hcw(_hcw[i]);
		census_current_room_contamIdx.push_back(-NA_FLOAT); // no room is allocated to HCW, but need a dummy value for tall-format dataframe
	}
}


vector< vector<uint> > HealthCareSetting::find_available_room(){
	
	vector< vector<uint> > res(_ward.size());
	
	for(uint w=0; w<_ward.size(); w++)
	{
		if( w != ISO_WARD){   // <-- Isolation ward is excluded
			vector<uint> tmp;
			for(uint i=0; i<_ward[w]->get_room().size(); i++){
				if(!_ward[w]->get_room(i)->is_full())
					tmp.push_back(i);
			}
			res[w] = tmp;
		}
	}
	return res;
}

void HealthCareSetting::event_new_admission(float time){
	
	vector<vector<uint> > w = find_available_room();
	
	// Check if not full:
	vector<size_t> ws;
	for(uint i=0; i<w.size();i++) ws.push_back(w[i].size());
	
	bool notfull = (sum(ws)>0);
	
	// Check if the HCS is full
	// (i.e. all roms of all wards are full).
	// Not critical if this happens sometimes,
	// but may be if it remains full for a long time.
	if(!notfull & false){
		cout << "DEBUG (time "<<time<<") - patient refused (all rooms full)";
		cout << " Total # patients:"<< count_patients() <<endl;
	}
	
	if(notfull){
		// TO DO: change below using:
		// select_random_non_empty_ward()
		
		// Identify available wards:
		vector<uint> idx_ward_available;
		for(uint k=0; k<w.size(); k++){
			if(w[k].size()>0)
				idx_ward_available.push_back(k);
		}
		
		// Choose one ward randomly among all available ones:
		uniform_int_distribution<uint> unif(0, (uint)idx_ward_available.size()-1);
		uint a = unif(RANDOM_GENERATOR);
		uint ward_slct = idx_ward_available[a];
		
		// Choose randomly the available room within this ward:
		uint n = (uint) w[ward_slct].size();
		uniform_int_distribution<uint> unif2(0, n-1);
		uint b = unif2(RANDOM_GENERATOR);
		uint room_slct = w[ward_slct][b];
		
		// Patient features:
		float DoS	= draw_DoS();
		float sIdx	= draw_sIdx();
		
		admission_new_patient(sIdx, DoS,
							  _contamIdx_suscept_prm,
							  _decay_contamIdx_halftime,
							  _ward[ward_slct]->get_room(room_slct),
							  time);
	}
}

void HealthCareSetting::discharge_patient(uint i_ward,
										  uint i_room,
										  uint i_position){
	_ward[i_ward]->get_room(i_room)->remove_patient(i_position);
}

void HealthCareSetting::discharge_patient_by_uid(uint i_ward,
												 uint i_room,
												 uint patient_uid){
	Room* r = _ward[i_ward]->get_room(i_room);
	
	bool found = false;
	uint i_position = NA_UINT;
	for(uint i=0; i < r->get_patient().size(); i++){
		if(r->get_patient(i)->get_uid() == patient_uid){
			found = true;
			i_position = i;
			break;
		}
	}
	stopif(!found,"Patient " + to_string(patient_uid) + " to discharged not found!");
	_ward[i_ward]->get_room(i_room)->remove_patient(i_position);
}


void HealthCareSetting::discharge_all(){
	int trace_level = 0;
	
	// Loop on all wards and rooms
	// to reach all patients:
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			for(uint p=0; p<_ward[w]->get_room(r)->count_patients(); p++)
			{
				float ta  = _ward[w]->get_room(r)->get_patient(p)->get_timeAdmission();
				float DoS = _ward[w]->get_room(r)->get_patient(p)->get_DoS();
				
				if( ta + DoS < _current_time ) // :discharge time reached.
				{
					discharge_patient(w, r, p);
					// DEBUG:
					if(trace_level>0){
						if(_ward[w]->get_room(r)->count_patients()==0)
							cout << "Ward #"<<w<<" room #"<<_ward[w]->get_room(r)->get_uid()<<" is empty."<<endl;
					}
				}
			}
		}
	}
}


void HealthCareSetting::admission_new_all(){
	poisson_distribution<uint> dpois(_admission_rate * _time_step);
	uint n = dpois(RANDOM_GENERATOR);
	for(uint i=0; i<n; i++) event_new_admission(_current_time);
}


vector<Patient*> HealthCareSetting::draw_patients_contacted_by_HCW(HCW* hcw, uint n){
	
	vector<Patient*> res;
	
	// Check if a ward is allocated or not:
	uint w_first = 0;
	uint w_last = (uint)_ward.size()-1;
	
	if(hcw->get_subtype() == "allocWard"){
		w_first = hcw->get_wardAssigned();
		w_last  = w_first;
	}
	//	cout<<"HCW contacting:"<<endl;
	//	hcw->show();
	
	// Retrieve all patients eligible to be drawn
	vector<Patient*> ep;
	for(uint w=w_first; w<= w_last; w++){
		//		cout << " DEBUG: drawing in ward #"<<w <<endl;
		for(uint i=0; i<_ward[w]->get_room().size(); i++){
			for(uint k=0; k<_ward[w]->get_room(i)->get_patient().size(); k++){
				ep.push_back(_ward[w]->get_room(i)->get_patient(k));
			}
		}
	}
	
	// Select randomly 'n' out of all the eligible patients found:
	
	// if not enough, then return all of them!
	if(ep.size() <= n) res = ep;
	
	if(ep.size() > n){
		uniform_int_distribution<uint> dunif(0, (uint)ep.size()-1);
		for(uint i=0; i<n; i++){
			uint idx = dunif(RANDOM_GENERATOR);
			res.push_back(ep[idx]);
		}
	}
	return res;
}

uint HealthCareSetting::draw_n_visits_HCW_patient(float dt){
	poisson_distribution<uint> dpois(dt * _lambda_prm[0]);
	return(dpois(RANDOM_GENERATOR));
}


void HealthCareSetting::contacts_one_HCW(HCW* hcw, bool recordContacts){
	
	// First, check if this HCW is on duty:
	bool on_duty = hcw->onDuty(_current_time, _time_step);
	
	// The HCW is on duty, so draw contacts with patients.
	// Draws the number of patients contacted as well as
	// the duration of each contact.
	if(on_duty){
		uint n = draw_n_visits_HCW_patient(_time_step);
		vector<Patient*> P = draw_patients_contacted_by_HCW(hcw, n);
		for(uint k=0; k<P.size(); k++)
		{
			uint w    = P[k]->get_current_ward_uid();
			uint ruid = P[k]->get_current_room_uid();
			Room* r = get_room_by_uid(ruid, w);
			
			float tau = draw_contact_HCW_patient_duration(_mean_contact_HCW_patient_duration_minutes,
														  _var_contact_HCW_patient_duration_minutes);
			contact_HCW_patient(hcw, P[k], r, tau,
								_contact_indiv_halftime,
								_contact_xchg_contamIdx,
								_visit_xchg_contamIdx,
								_visit_room_halftime,
								_disease,
								_current_time,
								_time_step);
			
			if(recordContacts) _contact_HCW_patient_dur_rec.push_back(tau);
		}
	}
}


void HealthCareSetting::draw_all_contacts_patient_HCW(bool recordContacts){
	
	// For HCW without any ward assigned:
	for(uint i=0; i<_hcw.size(); i++){
		contacts_one_HCW(_hcw[i], recordContacts);
	}
	
	// For HCW with a ward assigned:
	for(uint w=0; w<_ward.size(); w++){
		for(uint j=0; j<_ward[w]->get_hcw().size(); j++){
			contacts_one_HCW(_ward[w]->get_hcw(j), recordContacts);
		}
	}
}


Room* HealthCareSetting::get_room_by_pos(uint room_pos, uint ward_uid){
	return _ward[ward_uid]->get_room(room_pos);
}

Room* HealthCareSetting::get_room_by_uid(uint room_uid, uint ward_uid){
	
	bool found = false;
	uint ii = 0;
	
	for(uint i=0; i<_ward[ward_uid]->get_room().size(); i++){
		if(_ward[ward_uid]->get_room(i)->get_uid() == room_uid){
			found = true;
			ii = i;
			break;
		}
	}
	string errmsg = "Room uid:" + to_string(room_uid) + " not found in ward uid: " + to_string(ward_uid);
	stopif(!found, errmsg);
	return _ward[ward_uid]->get_room(ii);
}


vector<Patient*> HealthCareSetting::retrieve_all_patients(){
	vector<Patient*> res;
	
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			for(uint p=0; p<_ward[w]->get_room(r)->get_patient().size(); p++){
				res.push_back(_ward[w]->get_room(r)->get_patient(p));
			}
		}
	}
	return res;
}

vector<Patient*> HealthCareSetting::draw_patients(uint n){
	
	vector<Patient*> allp = retrieve_all_patients();
	
	stopif(allp.size()<n, "Asking to draw too many patients!");
	
	//random_shuffle(allp.begin(), allp.end()); // <-- use another default NG and mess up R interface
	std::shuffle(std::begin(allp), std::end(allp), RANDOM_GENERATOR);
	
	vector<Patient*> res(n);
	for(uint i=0; i<n; i++) res[i] = allp[i];
	
	return res;
}


void HealthCareSetting::seed_initial_infection(float prevalence_percent){
	
	vector<Patient*> allp = retrieve_all_patients();
	uint N = (uint)allp.size();
	uint p = (uint)( prevalence_percent * N);
	
	if(p==0){
		cout << "* * WARNING * * Cannot seed infection: no one drawn (initial prevalence too low)."<<endl;
		return;
	}
	vector<Patient*> pinf = draw_patients(p);
	cout << "Seeding initial infection:" << endl;
	float t0 = 0.0;
	for(uint i=0; i<pinf.size(); i++){
		pinf[i]->colonize(t0, _disease);
		pinf[i]->set_acquisitionType("initial");
		cout << " Patient #" << pinf[i]->get_uid() << " initial colonization." <<endl;
	}
}


void HealthCareSetting::update_all_clinical(){
	
	// Free HCW without assigned ward
	for(uint i=0; i<_hcw.size(); i++){
		_hcw[i]->update_clinical(_current_time, _disease);
	}
	
	// HCW with ward assigned
	for(uint w=0; w<_ward.size(); w++){
		for(uint i=0; i<_ward[w]->get_hcw().size(); i++){
			_ward[w]->get_hcw(i)->update_clinical(_current_time, _disease);
		}
	}
	
	// Patients
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			for(uint p=0; p<_ward[w]->get_room(r)->get_patient().size(); p++)
			{
				Patient* Pat = _ward[w]->get_room(r)->get_patient(p);
				Pat->update_clinical_patient(_current_time, _disease);
				
			}
		}
	}
}


void HealthCareSetting::update_decay_contamIdx_patient_rooms(){
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			_ward[w]->get_room(r)->decay_contamIdx(_current_time);
		}
	}
}


void HealthCareSetting::update_envContam_other(){
	// Calculate the mean contamination level
	// across all patient rooms:
	float m = 0.0;
	uint cnt = 0;
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			m += _ward[w]->get_room(r)->get_contamIdx();
			//			cout << "DEBUG:: w="<<w<<" r="<<r<<" m="<<_ward[w]->get_room(r)->get_contamIdx() ;
			//			cout << " uid="<< _ward[w]->get_room(r)->get_uid() <<endl;
			cnt++;
		}
	}
	m = m / cnt;
	
	// Draw the contamination level with a mean
	// at 'm' and a std dev at 0.5*m
	// (such that CV = 0.5)
	normal_distribution<float> dnorm(m, 0.5*m);
	_envContam_other = dnorm(RANDOM_GENERATOR) * _ratio_contam_other_room;
	if(_envContam_other<0) _envContam_other = 0.0;
}



void HealthCareSetting::update_all_infectIdx_vax(){
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			for(uint p=0; p<_ward[w]->get_room(r)->get_patient().size(); p++)
			{
				_ward[w]->get_room(r)->get_patient(p)->update_infectIdx_vax(_current_time, _vaccine);
			}
		}
	}
}


void HealthCareSetting::update_isoDur(){
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			
			vector<Patient* > Pat = _ward[w]->get_room(r)->get_patient();

			for(uint i=0; i<Pat.size(); i++){
				if(Pat[i]->get_isIsolated())
				{
					uint nr = Pat[i]->get_nRelapses();
					if(nr==0)
						Pat[i]->incr_isoDurFirst(_time_step);
					else
						Pat[i]->incr_isoDurRelapses(_time_step);
				}
			}
		}
	}
}



void HealthCareSetting::update_totDurSymptoms(){
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			
			vector<Patient* > Pat = _ward[w]->get_room(r)->get_patient();
			
			for(uint i=0; i<Pat.size(); i++){
				if(Pat[i]->get_currentlySymptomatic() )
				{
					Pat[i]->incr_totDurSymptoms(_time_step);
				}
			}
		}
	}
}

void HealthCareSetting::draw_acq_envContam_other(){
	
	bernoulli_distribution dbern_visit(_proba_patient_visit_other * _time_step);
	bernoulli_distribution dbern_exposed(_envContam_other);
	
	// Acquisition by Patients:
	
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			for(uint i=0; i<_ward[w]->get_room(r)->get_patient().size(); i++){
				// Patient visited other places,
				// and are exposed to pathogen:
				if(dbern_visit(RANDOM_GENERATOR) &&
				   dbern_exposed(RANDOM_GENERATOR))
				{
					colonization_attempt(_ward[w]->get_room(r)->get_patient(i),
										 _disease,
										 _current_time,
										 _time_step,
										 "environment");
				}
			}
		}
	}
	
	// All HCWs are supposed to visit communal spaces all the time,
	// hence no condition as with patients.
	
	// HCWs assigned to ward:
	for(uint w=0; w<_ward.size(); w++){
		for(uint i=0; i<_ward[w]->get_hcw().size() ; i++){
			if(_ward[w]->get_hcw(i)->onDuty(_current_time, _time_step) &&
			   dbern_exposed(RANDOM_GENERATOR))
			{
				colonization_attempt(_ward[w]->get_hcw(i),
									 _disease,
									 _current_time,
									 _time_step,
									 "environment");
			}
		}
	}
	
	// Free HCWs not linked to a ward:
	for(uint i=0; i<_hcw.size() ; i++){
		if( _hcw[i]->onDuty(_current_time, _time_step) &&
		   dbern_exposed(RANDOM_GENERATOR))
		{
			colonization_attempt(_hcw[i],
								 _disease,
								 _current_time,
								 _time_step,
								 "environment");
		}
	}
}


void HealthCareSetting::reset_contamIdx_one_HCW(HCW* hcw, float reset_level_max){
	
	uniform_real_distribution<float> U(0.0, reset_level_max);
	float c = hcw->get_contamIdx();
	if(c > reset_level_max){
		float new_ci = U(RANDOM_GENERATOR);
		hcw->set_contamIdx(new_ci);
	}
}

void HealthCareSetting::reset_contamIdx_all_HCW(){
	// Free HCWs not linked to a ward:
	for(uint i=0; i<_hcw.size() ; i++){
		reset_contamIdx_one_HCW(_hcw[i], _reset_contamIdx);
	}
	// HCWs assigned to ward:
	for(uint w=0; w<_ward.size(); w++){
		for(uint i=0; i<_ward[w]->get_hcw().size() ; i++){
			reset_contamIdx_one_HCW(_ward[w]->get_hcw(i), _reset_contamIdx);
		}
	}
}


void HealthCareSetting::adjust_DoS(){
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			for(uint p=0; p<_ward[w]->get_room(r)->get_patient().size(); p++){
				_ward[w]->get_room(r)->get_patient(p)->adjust_DoS(_current_time);
			}
		}
	}
}

void HealthCareSetting::one_patient_envContam_own_room(Patient* p,
													   float curr_time){
	uint room_uid = p->get_current_room_uid();
	uint ward_uid = p->get_current_ward_uid();
	Room* r = get_room_by_uid(room_uid, ward_uid);
	
	// Patient shed and acquire spores to/from patient's room:
	shed_in_room(p, r, _time_step, curr_time, _contact_indiv_halftime);
	spore_acq_from_room(p, r, _time_step, curr_time,
						_visit_room_halftime,
						_visit_xchg_contamIdx);
	
	// Potential colonization event:
	colonization_attempt(p, _disease, curr_time, _time_step, "environment");
}


void HealthCareSetting::all_patients_envContam_own_room(){
	vector<Patient*> p = retrieve_all_patients();
	for(uint i=0; i<p.size(); i++){
		one_patient_envContam_own_room(p[i], _current_time);
	}
}


void HealthCareSetting::print_ward_room(){
	for(uint w=0; w<_ward.size(); w++){
		cout << endl << "------- Ward #"<< w << endl;
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			cout << " > room uid: " << room[r]->get_uid();
			cout << " ; name: " << room[r]->get_name();
			cout << " ; max beds: " << room[r]->get_max_patient() ;
			cout << " ; patients: " << room[r]->count_patients() <<endl;
		}
	}
}

vector<Patient*> HealthCareSetting::retrieve_onset_patients(float curr_time,
															float timestep){
	vector<Patient*> res;
	
	for(uint w=0; w<_ward.size(); w++){
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			vector<Patient*> P = room[r]->get_patient();
			for(uint k=0; k<P.size(); k++){
				float t_onset = P[k]->get_timeOnset();
				
				// Special case for patients who relapse:
				if(P[k]->get_nRelapses()>0){
					t_onset = P[k]->get_timeRelapse();
				}
				// Retrieve patient if onset occurs
				// during this simulation timestep:
				if(is_in_current_timestep(t_onset)){
					res.push_back(P[k]);
				}
			}
		}
	}
	return res;
}


Room* HealthCareSetting::patient_room(Patient *P){
	Room* r = get_room_by_uid(P->get_current_room_uid(),
							  P->get_current_ward_uid());
	return r;
}


void HealthCareSetting::move_patient(Patient* P, Room* r){
	Room* old_room = patient_room(P);
	uint idx = old_room->find_patient_position(P->get_uid());
	// Add before removing (=deleting) it:
	r->add_patient(P);
	old_room ->remove_patient(idx);
}


void HealthCareSetting::isolate(vector<Patient *> P){

	Ward* iso_ward = _ward[ISO_WARD];
	
	for(uint i=0; i<P.size(); i++){
		
		Room* old_room = patient_room(P[i]);
		
		// Warns that no isolation room is available.
		// Does not stop execution, because the program
		// will keep on trying to isolate symptomatic
		// patients.
		if(iso_ward->is_full()){
			cout << " WARNING: Cannot isolate patient (uid="<<P[i]->get_uid()<<")";
			cout << " because all isolation rooms are full. ";
			cout << "Will keep on trying until one iso room is available.";
			cout << "(current time: "<<_current_time <<")" <<endl;
		}
		else{
			Room* r = iso_ward->find_available_room();
			move_patient(P[i], r);
			P[i]->set_timeIsolationStart(_current_time);
			P[i]->set_isIsolated(true);
			
//			cout<<"DEBUG: uid:"<<P[i]->get_uid();
//			cout<<" ; isolation started at t = "<<P[i]->get_timeIsolationStart()<<endl;
		}
		// Set the time for extra cleaning of patient's room
		// at the next simulation time step:
		old_room->set_time_extra_cleaning(_current_time + _time_step);
	}
}


void HealthCareSetting::draw_onset_iso_lag(){
	// Retrieve all patients that have onset of symptoms:
	vector<Patient*> P = retrieve_onset_patients(_current_time, _time_step);
	
	// If any, draw the onset-to-isolation lag:
	unsigned long n = P.size();
	if(n>0) {
		for(unsigned long i=0; i<n; i++){
			exponential_distribution<float> dExp(1.0 / _onset_isolation_lag_mean);
			float t_iso_start = dExp(RANDOM_GENERATOR);
			P[i]->set_timeIsolationStart(_current_time + t_iso_start);
		}
	}
}


vector<Patient *> HealthCareSetting::retrieve_patients_to_isolate(float curr_time,
																  float timestep){
	vector<Patient*> res;
	
	for(uint w=0; w<_ward.size(); w++){
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			vector<Patient*> P = room[r]->get_patient();
			for(uint k=0; k<P.size(); k++)
			{
				float t_onset = P[k]->get_timeIsolationStart();
				
				// Retrieve patient if isolation start time
				// occurs during this simulation timestep:
				if(is_in_current_timestep(t_onset)){
					res.push_back(P[k]);
				}
			}
		}
	}
	return res;
}

void HealthCareSetting::start_isolation(){
	vector<Patient*> P = retrieve_patients_to_isolate(_current_time, _time_step);
	if(P.size()>0) {
		isolate(P);
	}
}


// WARNING: if d>0, doesn't record isolation end time properly.
vector<Patient*> HealthCareSetting::retrieve_resolution_isopatients(float curr_time,
																	float timestep,
																	float d){
	vector<Patient*> v;
	
	// Look only in isolation ward.
	Ward* W = _ward[ISO_WARD];
	
	for(uint i=0; i < W->get_room().size(); i++){
		for(uint k=0; k < W->get_room(i)->get_patient().size(); k++)
		{
			Patient* P = W->get_room(i)->get_patient(k);
			float t_res = P->time_resolution();
			
			if(t_res + d < curr_time &&
			   curr_time <= t_res + d + timestep)
				v.push_back(P);
		}
	}
	return v;
}


vector<uint> HealthCareSetting::find_non_empty_ward(){
	vector<uint> res ;
	for(uint i=0; i< _ward.size(); i++){
		if(! (i == ISO_WARD) && !_ward[i]->is_full() )
			res.push_back(i);
	}
	return(res);
}

Ward* HealthCareSetting::select_random_non_empty_ward(){
	Ward* W = nullptr;
	
	vector<uint> idx = find_non_empty_ward();
	
	if(idx.size()>0){
		uniform_int_distribution<uint> U(0,(uint)idx.size()-1);
		uint idx_rand = U(RANDOM_GENERATOR);
		W = _ward[idx[idx_rand]];
	}
	return W;
}

// WARNING: if d>0, doesn't record isolation end time properly.
void HealthCareSetting::end_isolation(){
	int trace_level = 0;
	
	// Patients end their isolation as soon as the symptoms resolve:
	vector<Patient*> Pat0 = retrieve_resolution_isopatients(_current_time,
															_time_step, 0);
	
	// Patients who did not have the chance to get
	// a room available when their symptom resolve:
	vector<Patient*> Pat  = retrieve_end_isolation_ASAP();
	
	// Concatenate Pat and Pat0 into Pat, that is
	// patients that had just their resolution and
	// patients who did not find a room at the time
	// of resolution are candidates for ending isolation:
	Pat.insert( Pat.end(), Pat0.begin(), Pat0.end() );
	
	for(uint k=0; k<Pat.size(); k++){
		
		Ward* W =select_random_non_empty_ward();
		
		if(W == nullptr){
			if(trace_level>0){
				cout << "Warning: Patient #"<< Pat[k]->get_uid();
				cout <<" cannot leave isolation room because all other rooms are full." << endl;
			}
			Pat[k]->set_endIsolationASAP(true);
		}
		
		if(W != nullptr){
			Room* r = W->find_available_room();
			
			if(trace_level>0){
				cout << "Patient #"<<Pat[k]->get_uid() ;
				cout << " ends its isolation, go to room #"<<r->get_uid();
				cout << " at time "<<_current_time<< endl;
				if(trace_level>1) r->show();
			}
			Pat[k]->set_timeIsolationEnd(_current_time);
			Pat[k]->set_endIsolationASAP(false);
			Pat[k]->set_isIsolated(false);
			
			float ts = Pat[k]->get_timeIsolationStart();
			float te = Pat[k]->get_timeIsolationEnd();
			stopif(te<ts, "Isolation ends _before_ it started!");
			
			move_patient(Pat[k],r);
		}
	}
}


vector<Patient* > HealthCareSetting::retrieve_end_isolation_ASAP(){
	
	vector<Patient* > res;
	
	for(uint r=0; r<_ward[ISO_WARD]->get_room().size(); r++){
		for(uint p=0; p<_ward[ISO_WARD]->get_room(r)->count_patients(); p++)
		{
			Patient* pat = _ward[ISO_WARD]->get_room(r)->get_patient(p);
			bool isol_asap  = pat->get_endIsolationASAP();
			if(isol_asap){
				res.push_back(pat);
			}
		}
	}
	
	return res;
}


bool HealthCareSetting::is_in_current_timestep(float t){
	bool res = false;
	if(_current_time <= t && t < _current_time + _time_step) res = true;
	return res;
}


void HealthCareSetting::extra_cleaning(){

	// Loop through all rooms in all wards:
	for(uint w=0; w<_ward.size(); w++){
		for(uint i=0; i<_ward[w]->get_room().size(); i++)
		{
			Room* r = _ward[w]->get_room(i);
			
			// Isolation room are always cleaned:
			if(w==ISO_WARD){
				r->extra_cleaning(_clean_room_efficacy);
			}
			
			// If it's time to clean this room:
			else{
				float tec = r->get_time_extra_cleaning();
				if(is_in_current_timestep(tec)) {
					float oldidx = r->get_contamIdx();
					r->extra_cleaning(_clean_room_efficacy);
					stopif(oldidx < r->get_contamIdx(),
						   "Problem in function: HealthCareSetting::extra_cleaning()");
				}
			}
		}
	}
}


vector<vector<uint> > HealthCareSetting::find_transferable_patients(){
	
	vector<vector<uint> > res;
	int trace_level = 0;
	
	for(uint w=0; w<_ward.size(); w++){
		for(uint r=0; r<_ward[w]->get_room().size(); r++){
			for(uint p=0; p<_ward[w]->get_room(r)->count_patients(); p++)
			{
				float ta  = _ward[w]->get_room(r)->get_patient(p)->get_timeAdmission();
				float DoS = _ward[w]->get_room(r)->get_patient(p)->get_DoS();
				if( ta + DoS < _current_time )
				{
					// DEBUG ---
					//					cout << current_time << " DEBUG: discharge patient ";
					//					cout <<_ward[w]->get_room(r)->get_patient(p)->get_uid()<< " in room "<<r;
					//					cout << " DoS=" << DoS << " < "<< current_time - ta << endl;
					// ---
					uint uid_patient =  _ward[w]->get_room(r)->get_patient(p)->get_uid();
					res.push_back({w,r,uid_patient});
					
					//_ward[w]->get_room(r)->remove_patient(p);
					
					// DEBUG:
					if(trace_level>0){
						if(_ward[w]->get_room(r)->count_patients()==0)
							cout << "Ward #"<<w<<" room #"<<_ward[w]->get_room(r)->get_uid()<<" is empty."<<endl;
					}
				}
			}
		}
	}
	return res;
}


void HealthCareSetting::receive_transfered_patient(Patient* p,
												   float DoS,
												   float contamIdx_suscept_prm,
												   float decay_contamIdx_halftime,
												   Room* r,
												   float time){
	// Reset what is necessary for the new HCS:
	p->set_timeAdmission(time);
	p->set_DoS(DoS);
	p->set_contamIdx_suscept_prm(contamIdx_suscept_prm);
	p->set_decay_contamIdx_halftime(decay_contamIdx_halftime);
	
	r->add_patient(p);
	
	//cout <<"Transfered patient  uid:"<< p->get_uid()<<" ; in room "<< r->get_uid();
	//cout << " ; colonized@admin: "<< p->get_isColonized() << " ; contamIdx = " << p->get_contamIdx() << endl;
	
	// Record in the registry as admission:
	_registry.add_record(p->get_uid(), r->get_uid(), r->get_ward_uid() );
}

// DELETE WHEN SURE
void HealthCareSetting::vaccinate(){
//	bool vaxstrategyFound = false;
//	// Switch over all vaccination
//	// strategies implemented:
//	if(_vaxStrategy == "none"){
//		// do nothing!
//		vaxstrategyFound = true;
//	}
//	if(_vaxStrategy == "pre-emptive"){
//		vax_newAdm_random(_proba_vax_newAdmission);
//		vaxstrategyFound = true;
//	}
//	if(_vaxStrategy == "frailty"){
//		vax_newAdm_frailty(_sIdxMin_vax_frailty,
//						   _DoSMin_vax_frailty,
//						   _proba_vax_newAdmission);
//		vaxstrategyFound = true;
//	}
//	// Check if the vaccination strategy was found:
//	string errmsg = "ERROR --> Vaccination strategy name unknown:" + _vaxStrategy;
//	stopif(!vaxstrategyFound, errmsg);
}


void HealthCareSetting::vaccinate_individual(Patient* p){
	
	bool vaxstrategyFound = false;
	
	// Switch over all vaccination
	// strategies implemented:
	
	if(_vaxStrategy == "none")
	{   // do nothing!
		vaxstrategyFound = true;
	}
	
	if(_vaxStrategy == "frailty"){
		float sIdx  = p->get_susceptIdx();
		float dos 	= p->get_DoS();
		float t_adm = p->get_timeAdmission();
		
		// Filter only the patients admitted after
		// start of vaccination campaign and
		// who have a high susceptibility index:
		if(t_adm >= _vax_start_time &&
		   sIdx  >= _sIdxMin_vax_frailty &&
		   dos   >= _DoSMin_vax_frailty)
		{
			// Draws a vaccination event
			// for this new admission:
			bernoulli_distribution dbern(_proba_vax_newAdmission);
			bool doVax = dbern(RANDOM_GENERATOR);
			if(doVax){
				p->vaccinate(_current_time, _vaccine);
			}
		}
		vaxstrategyFound = true;
	}
	
	// Check if the vaccination strategy was found:
	string errmsg = "ERROR --> Vaccination strategy name unknown:" + _vaxStrategy;
	stopif(!vaxstrategyFound, errmsg);
}


void HealthCareSetting::vax_newAdm_random(float proba){
	
	// loop through all patients, keep the
	// new admissions only, and attempt to vaccinate them:
	for(uint w=0; w<_ward.size(); w++){
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			vector<Patient*> P = room[r]->get_patient();
			for(uint k=0; k<P.size(); k++){
				float t_adm = P[k]->get_timeAdmission();
				if( is_in_current_timestep(t_adm) ){
					// Draws a vaccination event
					// for this new admission:
					bernoulli_distribution dbern(proba);
					bool doVax = dbern(RANDOM_GENERATOR);
					if(doVax){
						P[k]->vaccinate(_current_time, _vaccine);
					}
				}
			}
		} // end for r
	} // end for w
}// end function


void HealthCareSetting::vax_newAdm_frailty(float sIdxmin,
										   float DoSmin,
										   float proba){
	// loop through all patients, keep the
	// new admissions only, and attempt to vaccinate them:
	for(uint w=0; w<_ward.size(); w++){
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			vector<Patient*> P = room[r]->get_patient();
			for(uint k=0; k<P.size(); k++)
			{
				float t_adm = P[k]->get_timeAdmission();
				float sIdx  = P[k]->get_susceptIdx();
				float dos 	= P[k]->get_DoS();
				
				// Filter only the patients admitted now
				// who have a high susceptibility index:
				if( is_in_current_timestep(t_adm) &&
				   sIdx >= sIdxmin &&
				   dos >= DoSmin){
					// Draws a vaccination event
					// for this new admission:
					bernoulli_distribution dbern(proba);
					bool doVax = dbern(RANDOM_GENERATOR);
					if(doVax){
						P[k]->vaccinate(_current_time, _vaccine);
					}
				}
			}
		} // end for r
	} // end for w
}


vector<Patient*> HealthCareSetting::draw_patients_relapse(float proba,
														  float shape){
	vector<Patient*> res;
	
	for(uint w=0; w<_ward.size(); w++){
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			vector<Patient*> P = room[r]->get_patient();
			for(uint k=0; k<P.size(); k++)
			{
				uint nr	= P[k]->get_nRelapses();
				
				if( P[k]->get_isColonized() &&
				   (P[k]->get_nInfections()>0) &&
				    nr < MAX_N_RELAPSES )
				{
					float t_res = P[k]->time_resolution();
					float t_clr = P[k]->time_clearance();
					float t_ons = P[k]->get_timeOnset();
					bool cs 	= P[k]->get_currentlySymptomatic();
					
					// Timing condition:
//					bool cond_timing = (_current_time > t_res +1) && (_current_time < t_clr -1);
					bool cond_timing =
					(_current_time > t_res +0) &&
					(_current_time < t_clr) &&
					(t_ons < _current_time) &&
					(!cs);
					
					// First relapse:
					bool cond_0 = (nr==0 && cond_timing);
					
					// Already had at least one relapse:
					bool cond_1 = (nr>0 && cond_timing);
					
					if(cond_0 || cond_1)
					{
						// Higher susceptibility (frailty) makes patient
						// more likely to relapse:
						float sidx = P[k]->get_susceptIdx();
						float sidx2 = sidx<shape? 0 : 1.0;
						float ps = proba * sidx2;
						ps = (ps>1.0 ? 1.0 : ps);
						
						// Draw relapse event:
						bernoulli_distribution dbern(ps);
						bool relapse = dbern(RANDOM_GENERATOR);
						if(relapse)	res.push_back(P[k]);
						
						// DEBUG:
//						cout <<" DEBUG>> ";
//						cout <<" time:"<<_current_time;
//						cout <<" uid:"<< P[k]->get_uid();
//						cout <<" nRelapse="<<nr;
//						cout <<" t_onset="<<t_ons <<" ; t_res="<<t_res;
//						cout <<" ps="<<ps;
//						cout <<" relapse="<<relapse;
//						cout<< endl;
						
					}
				} // endif isColonized
			}
		} // end for r
	} // end for w
	return res;
}


void HealthCareSetting::relapses(){
	
	float hr = _disease.get_hazard_relapse();
	float sr = _disease.get_shape_relapse();
	vector<Patient*> pr = draw_patients_relapse(hr, sr);
	size_t n = pr.size();
	
	if(n>0){
		for(uint i=0; i<n; i++){
			pr[i]->relapse_triggered(_current_time, _disease);
			pr[i]->adjust_DoS_relapse(_current_time);
			
//			cout<<"DEBUG:Relapse triggered uid:"<< pr[i]->get_uid();
//			cout<<" ; relapse #"<<pr[i]->get_nRelapses();
//			cout<<" ; time:"<<_current_time << endl;
		}
		// Relapsing patients are isolated again:
		isolate(pr);
	}
}

void HealthCareSetting::draw_colectomy(){
	
	// loop through all patients
	for(uint w=0; w<_ward.size(); w++){
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			vector<Patient*> P = room[r]->get_patient();
			for(uint k=0; k<P.size(); k++){
				
				float t_acq = P[k]->get_timeAcquisition();
				float iidx  = P[k]->get_infectIdx();
				
				// Filter only the patients who have
				// just acquired infection and
				// who have a high infection index:
				if( is_in_current_timestep(t_acq) &&
				   iidx >= 0.75){ // TO DO: do not hard code the threshold.
					
					// Draws a colectomy event
					// for this new infection:
					bernoulli_distribution dbern(_proba_colectomy);
					bool colectomy = dbern(RANDOM_GENERATOR);
					if(colectomy){
						P[k]->set_colectomy(true);
					}
				}
			}
		} // end for r
	} // end for w
}


void HealthCareSetting::draw_fmt(){
	
	for(uint w=0; w<_ward.size(); w++){
		vector<Room*> room = _ward[w]->get_room();
		for(uint r=0; r<room.size(); r++){
			vector<Patient*> P = room[r]->get_patient();
			for(uint k=0; k<P.size(); k++){
				uint n = P[k]->get_nRelapses();
				// Filter only the patients
				// who had at least one relapse:
				if( n >0 ){
					// Draws a vaccination event
					// for this new admission:
					bernoulli_distribution dbern(_proba_fmt);
					bool fmt = dbern(RANDOM_GENERATOR);
					if(fmt){
						P[k]->set_fmt(true);
					}
				}
			}
		} // end for r
	} // end for w
	
}

