//
//  run.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-12.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "run.hpp"
#include "globalvar.hpp"
#include "Simulator.hpp"


/// Run test on a single HCS.
void run_test_single(){
	
	int trace_level = 0;
	// Matrix of movements. Construction by rows.
	vector<vector<float> > mvt_HCS;
	mvt_HCS.push_back({0.0});
	
	float prd_incubation_mean = 1;
	float prd_incubation_var  = 0.3;
	float prd_symptom_mean    = 3.5;
	float prd_symptom_var	  = 2;
	
	float proba_clear_prm = 2;
	float prd_clear_asymptom_mean = 7;
	float prd_clear_asymptom_var  = 3;
	float prd_clear_symptom_mean  = 20;
	float prd_clear_symptom_var   = 4;
	
	float contamIdx_max_asymptomatic  = 0.3;
	float contamIdx_ratio_symptomatic = 2.0;
	
	float hazardRelapse = 0.5;//0.01;
	float shape_relapse = 0.09;
	
	float infectIdx_init = 0.066;
	
	Disease D(prd_incubation_mean,
			  prd_incubation_var,
			  prd_symptom_mean,
			  prd_symptom_var,
			  proba_clear_prm,
			  prd_clear_asymptom_mean,
			  prd_clear_asymptom_var,
			  prd_clear_symptom_mean,
			  prd_clear_symptom_var,
			  contamIdx_max_asymptomatic,
			  contamIdx_ratio_symptomatic,
			  hazardRelapse,
			  shape_relapse,
			  infectIdx_init);
	
	D.show();
	
	HealthCareSetting* Q_1 = new HealthCareSetting;
	Q_1->set_disease(D);
	
	vector<string> room_type;
	vector<uint>   ward;
	vector<string> room_name;
	vector<uint>   max_patient;
	
	float occupancy 	= 0.7;
	float mean_sIdx 	= 0.5;
	float sd_sIdx   	= 0.2;
	float mean_sIdx_HCW = 0.01;
	float sd_sIdx_HCW   = 0.1;
	
	float mDoS 	= log(4.0);
	float sdDoS = 0.2;
	
	// All these vector must have the same size:
	room_type = {"isolation_room","isolation_room",
		"isolation_room","isolation_room",
		"patient_room","patient_room","patient_room", "patient_room"};
	ward = {0,0,0,0,1,1,2,3};
	room_name = {"ISO1","ISO2",
		"ISO3","ISO4",
		"P1","P2","P3", "P4"};
	max_patient = {1,1,1,1,2,4,6,3};
	
	
	uint n_iso_rooms = 10;
	uint n_ward_mean = 3;
	uint n_room_per_ward_mean = 4;
	uint n_beds_per_room_mean = 3;
	
	uint nHCW_ward = 17;
	uint nHCW_free = 8;
	
	string hcs_name = "Sunnybrook";
	string hcs_type = "hospital";
	
	float admission_rate = 3.1;
	
	vector<float> lambda_prm = {4.123};
	
	float decay_room_halftime = 3.456;
	float contamIdx_suscept_prm = 0.01;
	float decay_contamIdx_halftime = 5.678;
	float ratio_contam_other_room = 0.1;
	float prev_new_admission = 0.10;
	float proba_patient_visit_other = 0.2;
	float reset_contamIdx = 0.05;
	float mean_contact_HCW_patient_duration_minutes = 15.0;
	float var_contact_HCW_patient_duration_minutes = 10.0 ;
	float contact_indiv_halftime = 0.03;
	float contact_xchg_contamIdx = 0.7;
	float visit_xchg_contamIdx = 0.8;
	float visit_room_halftime    = 0.04;
	
	float clean_room_efficacy = 0.85;
	float onset_isolation_lag_mean = 1.3;
	
	float vax_efficacy_lo = 0.60;
	float vax_efficacy_hi = 0.90;
	
	float vax_proba_prev_col = 0.1;
	float vax_halftime       = 222;
	Vaccine vaxObj("the_vax",
				   vax_efficacy_lo,
				   vax_efficacy_hi,
				   vax_proba_prev_col,
				   vax_halftime);
	string vaxStrategy = "frailty";  // "pre-emptive";
	
	float proba_vax_newAdmission = 0.5;
	float sIdxMin_vax_frailty    = 0.4;
	float DoSMin_vax_frailty     = 2.0;
	float vax_start_time         = 300;
	
	float proba_colectomy = 0.01;
	float proba_fmt       = 0.12;
	
	uint uid_start_1 = 0;
	
	Q_1->setup(hcs_name,
			   hcs_type,
			   n_iso_rooms,
			   n_ward_mean,
			   n_room_per_ward_mean,
			   n_beds_per_room_mean,
			   nHCW_ward,
			   nHCW_free,
			   uid_start_1,
			   occupancy,
			   mean_sIdx,  sd_sIdx,
			   mean_sIdx_HCW,  sd_sIdx_HCW,
			   mDoS,  sdDoS,
			   admission_rate,
			   lambda_prm,
			   decay_room_halftime,
			   contamIdx_suscept_prm,
			   decay_contamIdx_halftime,
			   ratio_contam_other_room,
			   prev_new_admission,
			   proba_patient_visit_other,
			   reset_contamIdx,
			   mean_contact_HCW_patient_duration_minutes,
			   var_contact_HCW_patient_duration_minutes,
			   contact_indiv_halftime,
			   contact_xchg_contamIdx,
			   visit_xchg_contamIdx,
			   visit_room_halftime,
			   clean_room_efficacy,
			   onset_isolation_lag_mean,
			   vaxObj,
			   vaxStrategy,
			   proba_vax_newAdmission,
			   sIdxMin_vax_frailty,
			   DoSMin_vax_frailty,
			   vax_start_time,
			   proba_colectomy,
			   proba_fmt);
	
	Q_1->stats_architecture();
	
	if(trace_level>0){
		Q_1->print_ward_room();
		Q_1->show(true);
	}
	float horizon = 500; //1825;
	float timestep = 0.5;
	bool recordContacts = false;
	
	vector<HealthCareSetting*> HCS = {Q_1};
	
	Simulator S(HCS, horizon, timestep, recordContacts);
	S.set_mvt_HCS(mvt_HCS);
	
	S.run_simulation();
	
	vector<censusDF> df = S.get_census();
	vector<censusDF> df_light = S.get_census_light();
	
//	unsigned long n = df.size();
	censusDF a = df_light[df_light.size()-1];
	
	vector<float> x = a._contamIdx;
	cout << " DEBUG:" <<endl;
	displayVector(x);
	
	cout << " DONE run_test_single()" << endl;
}




