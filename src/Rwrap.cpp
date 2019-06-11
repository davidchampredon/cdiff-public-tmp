//
//  Rwrap.cpp
//  naf
//
//  Created by David CHAMPREDON on 2016-07-15.
//  Copyright (c) 2016 David CHAMPREDON. All rights reserved.
//

#include <Rcpp.h>
#include <iostream>
#include <typeinfo>

#include "HealthCareSetting.hpp"
#include "globalvar.hpp"
#include "Simulator.hpp"


using namespace std;
using namespace Rcpp;



// ====================================================================================
// ====================================================================================
// = = = =  M A I N    E X P O R T  = = = =
// ====================================================================================
// ====================================================================================


void print_list(List params, uint seed){
	
	string seedinfo = "(print_list -> seed="+to_string(seed)+") ";
	
	string hcs_name = params["hcs_name"]; cout<<seedinfo<<"hcs_name = "<<hcs_name<<endl;
	string hcs_type = params["hcs_type"];cout<<seedinfo<<"hcs_type = "<<hcs_type<<endl;
	
	float occupancy	= params["occupancy"];cout<<seedinfo<<"occupancy = "<<occupancy<<endl;
	float mean_sIdx	= params["sIdx_newPatient_prm1"];cout<<seedinfo<<"mean_sIdx = "<<mean_sIdx<<endl;
	float sd_sIdx	= params["sIdx_newPatient_prm2"];cout<<seedinfo<<"sd_sIdx = "<<sd_sIdx<<endl;
	float mDoS		= params["DoS_newPatient_prm1"];cout<<seedinfo<<"mDoS = "<<mDoS<<endl;
	float sdDoS		= params["DoS_newPatient_prm2"];cout<<seedinfo<<"sdDoS = "<<sdDoS<<endl;
	
	float mean_sIdx_HCW	= params["sIdx_HCW_prm1"];cout<<seedinfo<<"mean_sIdx_HCW = "<<mean_sIdx_HCW<<endl;
	float sd_sIdx_HCW	= params["sIdx_HCW_prm2"];cout<<seedinfo<<"sd_sIdx_HCW = "<<sd_sIdx_HCW<<endl;
	
	uint nNurses   = params["nNurses"];cout<<seedinfo<<"nNurses = "<<nNurses<<endl;
	uint nDoctors  = params["nDoctors"];cout<<seedinfo<<"nDoctors = "<<nDoctors<<endl;
	uint uid_start = params["uid_start"];cout<<seedinfo<<"uid_start = "<<uid_start<<endl;
	
	float admission_rate = params["admission_rate"];cout<<seedinfo<<"admission_rate = "<<admission_rate<<endl;
	
	float lambda_prm_1	 = params["lambda_prm_1"];cout<<seedinfo<<"lambda_prm_1 = "<<lambda_prm_1<<endl;
	float lambda_prm_2	 = params["lambda_prm_2"];cout<<seedinfo<<"lambda_prm_2 = "<<lambda_prm_2<<endl;
//	vector<float> lambda_prm = {lambda_prm_1, lambda_prm_2};
	
	float decay_room_halftime      = params["decay_room_halftime"];cout<<seedinfo<<"decay_room_halftime = "<<decay_room_halftime<<endl;
	float contamIdx_suscept_prm    = params["contamIdx_suscept_prm"];cout<<seedinfo<<"contamIdx_suscept_prm = "<<contamIdx_suscept_prm<<endl;
	float decay_contamIdx_halftime = params["decay_contamIdx_halftime"];cout<<seedinfo<<"decay_contamIdx_halftime = "<<decay_contamIdx_halftime<<endl;
	float ratio_contam_other_room  = params["ratio_contam_other_room"];cout<<seedinfo<<"ratio_contam_other_room = "<<ratio_contam_other_room<<endl;
	
	float prev_new_admission        = params["prev_new_admission"];cout<<seedinfo<<"prev_new_admission = "<<prev_new_admission<<endl;
	float proba_patient_visit_other = params["proba_patient_visit_other"];cout<<seedinfo<<"proba_patient_visit_other = "<<proba_patient_visit_other<<endl;
	float reset_contamIdx           = params["reset_contamIdx"];cout<<seedinfo<<"reset_contamIdx = "<<reset_contamIdx<<endl;
	
	float mean_contact_HCW_patient_duration_minutes = params["mean_contact_HCW_patient_duration_minutes"];
	float var_contact_HCW_patient_duration_minutes  = params["var_contact_HCW_patient_duration_minutes"];
	cout<<seedinfo<<"mean_contact_HCW_patient_duration_minutes = "<<mean_contact_HCW_patient_duration_minutes<<endl;
	cout<<seedinfo<<"var_contact_HCW_patient_duration_minutes = "<<var_contact_HCW_patient_duration_minutes<<endl;
	
	
	float contact_indiv_halftime = params["contact_indiv_halftime"];cout<<seedinfo<<"contact_indiv_halftime = "<<contact_indiv_halftime<<endl;
	float visit_room_halftime    = params["visit_room_halftime"];cout<<seedinfo<<"visit_room_halftime = "<<visit_room_halftime<<endl;
	
	float clean_room_efficacy	 = params["clean_room_efficacy"];cout<<seedinfo<<"clean_room_efficacy = "<<clean_room_efficacy<<endl;
	
	float onset_isolation_lag_mean	 = params["onset_isolation_lag_mean"];cout<<seedinfo<<"onset_isolation_lag_mean = "<<onset_isolation_lag_mean<<endl;
	
	
	float vax_efficacy_lo	= params["vax_efficacy_lo"]; cout<<seedinfo<<"vax_efficacy_lo = "<<vax_efficacy_lo<<endl;
	float vax_efficacy_hi	= params["vax_efficacy_hi"]; cout<<seedinfo<<"vax_efficacy_hi = "<<vax_efficacy_hi<<endl;
	float vax_ppc	= params["vax_proba_prevent_colonization"];
	cout<<seedinfo<<"vax_proba_prevent_colonization = "<<vax_ppc <<endl;
	
	float vax_halflife	= params["vax_halflife"];cout<<seedinfo<<"vax_halflife = "<<vax_halflife<<endl;
	string vax_name 	= params["vax_name"];cout<<seedinfo<<"vax_name = "<<vax_name<<endl;
	
	string vaxStrategy			 = params["vaxStrategy"];cout<<seedinfo<<"vaxStrategy = "<<vaxStrategy<<endl;
	float proba_vax_newAdmission = params["proba_vax_newAdmission"];cout<<seedinfo<<"proba_vax_newAdmission = "<<proba_vax_newAdmission<<endl;
	float sIdxMin_vax_frailty	 = params["sIdxMin_vax_frailty"];cout<<seedinfo<<"sIdxMin_vax_frailty = "<<sIdxMin_vax_frailty<<endl;
	float DoSMin_vax_frailty	 = params["DoSMin_vax_frailty"];cout<<seedinfo<<"DoSMin_vax_frailty = "<<DoSMin_vax_frailty<<endl;
	float vax_start_time		 = params["vax_start_time"];cout<<seedinfo<<"vax_start_time = "<<vax_start_time<<endl;
	
	float proba_colectomy		 = params["proba_colectomy"];cout<<seedinfo<<"proba_colectomy = "<<proba_colectomy<<endl;
	float proba_fmt				 = params["proba_fmt"];cout<<seedinfo<<"proba_fmt = "<<proba_fmt<<endl;
}

/// Setup a HealthCareSetting from a list which was populated in R.
HealthCareSetting* setup_from_list(List params, Disease D, uint seed){
	
	cout << "seed="<<seed<< ": Unpacking parameters from R list... " << endl;
	
	string seedinfo = "(setup_from_list-> seed = "+to_string(seed)+") ";
	
	print_list(params, seed);
	
	// unpack all parameters:
	string hcs_name = params["hcs_name"];
	string hcs_type = params["hcs_type"];
	
//	vector<string> room_type	= params["room_type"];
//	vector<uint> ward			= params["ward"];
//	vector<string> room_name	= params["room_name"];
//	vector<uint> max_patient	= params["max_patient"];
			
	uint n_iso_rooms			= params["n.iso"];
	uint n_ward_mean			= params["n.wards"];
	uint n_room_per_ward_mean	= params["rooms.per.ward.mean"];
	uint n_beds_per_room_mean	= params["beds.per.room.mean"];
	
//	cout << "seed="<<seed<< ": room_type.size()=" << room_type.size() << endl;
	
	float occupancy	= params["occupancy"];
	float mean_sIdx	= params["sIdx_newPatient_prm1"];
	float sd_sIdx	= params["sIdx_newPatient_prm2"];
	float mDoS		= params["DoS_newPatient_prm1"];
	float sdDoS		= params["DoS_newPatient_prm2"];
	
	float mean_sIdx_HCW	= params["sIdx_HCW_prm1"];
	float sd_sIdx_HCW	= params["sIdx_HCW_prm2"];
	
	uint nNurses   = params["nNurses"];
	uint nDoctors  = params["nDoctors"];
	uint uid_start = params["uid_start"];
	
	float admission_rate = params["admission_rate"];
	
	float lambda_prm_1	 = params["lambda_prm_1"];
	float lambda_prm_2	 = params["lambda_prm_2"];
	vector<float> lambda_prm = {lambda_prm_1, lambda_prm_2};
	
	float decay_room_halftime      = params["decay_room_halftime"];
	float contamIdx_suscept_prm    = params["contamIdx_suscept_prm"];
	float decay_contamIdx_halftime = params["decay_contamIdx_halftime"];
	float ratio_contam_other_room  = params["ratio_contam_other_room"];
	
	float prev_new_admission        = params["prev_new_admission"];
	float proba_patient_visit_other = params["proba_patient_visit_other"];
	float reset_contamIdx           = params["reset_contamIdx"];
	
	float mean_contact_HCW_patient_duration_minutes = params["mean_contact_HCW_patient_duration_minutes"];
	float var_contact_HCW_patient_duration_minutes  = params["var_contact_HCW_patient_duration_minutes"];
	
	float contact_indiv_halftime = params["contact_indiv_halftime"];
	float contact_xchg_contamIdx = params["contact_xchg_contamIdx"];
	float visit_xchg_contamIdx	 = params["visit_xchg_contamIdx"];
	float visit_room_halftime    = params["visit_room_halftime"];
	
	float clean_room_efficacy	 = params["clean_room_efficacy"];
	float onset_isolation_lag_mean	 = params["onset_isolation_lag_mean"];
	
	float vax_efficacy_lo	= params["vax_efficacy_lo"];
	float vax_efficacy_hi	= params["vax_efficacy_hi"];
	float vax_proba_prevent_colonization = params["vax_proba_prevent_colonization"];
	float vax_halflife	= params["vax_halflife"];
	string vax_name 	= params["vax_name"];
	
	Vaccine vaxObj(vax_name,
				   vax_efficacy_lo,
				   vax_efficacy_hi,
				   vax_proba_prevent_colonization,
				   vax_halflife);
	
	string vaxStrategy			 = params["vaxStrategy"];
	float proba_vax_newAdmission = params["proba_vax_newAdmission"];
	float sIdxMin_vax_frailty	 = params["sIdxMin_vax_frailty"];
	float DoSMin_vax_frailty	 = params["DoSMin_vax_frailty"];
	float vax_start_time		 = params["vax_start_time"];
	
	float proba_colectomy		 = params["proba_colectomy"];
	float proba_fmt				 = params["proba_fmt"];
	
	cout << "seed="<<seed <<" Unpacking parameters done."<<endl;
	
	// Create HealthCareSetting object:
	cout << "seed="<<seed<< ": Setting up HCS "<< endl ;
	HealthCareSetting* Q = new HealthCareSetting;
	
	// Warning: disease must be set before setup()
	Q->set_disease(D);
	
	Q->setup(hcs_name,
			 hcs_type,
			 n_iso_rooms,
			 n_ward_mean,
			 n_room_per_ward_mean,
			 n_beds_per_room_mean,
			 nNurses,  nDoctors,
			 uid_start,
			 occupancy,
			 mean_sIdx,  sd_sIdx,
			 mean_sIdx_HCW,  sd_sIdx_HCW,
			 mDoS, sdDoS,
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
	
	cout << "seed="<<seed << "HCS "<< Q->get_name() <<" built & populated. " << endl;
	
	return Q;
}

// [[Rcpp::export]]
List abmcdiff_one_simul(List hcsParams,
						List simParams,
						List diseaseParams,
						List mvtParams){
	
	try{
		int trace_level = 0;
		
		cout << " Entering Rwrap... " ;
		cout << endl;
		
		// === Unpack parameters from R ===
		
		// Disease:
		float prd_incubation_mean	= diseaseParams["prd_incubation_mean"];
		float prd_incubation_var	= diseaseParams["prd_incubation_var"];
		float prd_symptom_mean		= diseaseParams["prd_symptom_mean"];
		float prd_symptom_var		= diseaseParams["prd_symptom_var"];
		
		float proba_clear_prm			= diseaseParams["proba_clear_prm"];
		float prd_clear_asymptom_mean	= diseaseParams["prd_clear_asymptom_mean"];
		float prd_clear_asymptom_var	= diseaseParams["prd_clear_asymptom_var"];
		float prd_clear_symptom_mean	= diseaseParams["prd_clear_symptom_mean"];
		float prd_clear_symptom_var		= diseaseParams["prd_clear_symptom_var"];
		
		float contamIdx_max_asympt			= diseaseParams["contamIdx_max_asympt"];
		float contamIdx_ratio_symptomatic	= diseaseParams["contamIdx_ratio_symptomatic"];
		
		float hazard_relapse		= diseaseParams["hazard_relapse"];
		float shape_relapse			= diseaseParams["shape_relapse"];
		float infectIdx_init		= diseaseParams["infectIdx_init"];
		
		Disease D(prd_incubation_mean,
				  prd_incubation_var,
				  prd_symptom_mean,
				  prd_symptom_var,
				  proba_clear_prm	,
				  prd_clear_asymptom_mean,
				  prd_clear_asymptom_var,
				  prd_clear_symptom_mean,
				  prd_clear_symptom_var,
				  contamIdx_max_asympt,
				  contamIdx_ratio_symptomatic,
				  hazard_relapse,
				  shape_relapse,
				  infectIdx_init);
		
		cout<< "  Disease show in RWarp:" <<endl;
		D.show();
		
		// Simulation:
		
		float horizon	 	= simParams["horizon"];
		float timestep 		= simParams["timestep"];
		uint seedRNG   		= simParams["seedRNG"];
		bool lightCensus 	= simParams["lightCensus"];
		bool recordContacts = simParams["recordContacts"];
		
		cout << " ... parameters unpacked. " << endl;
		
		// =========================
		// === Call C++ function ===
		// =========================
		
		RANDOM_GENERATOR.seed(seedRNG);
		cout << " ... RNG seed set (value: "<< seedRNG <<")"<< endl;
		
		size_t n_hcs = hcsParams.size();
		
		// Setup healthcare settings:
		vector<HealthCareSetting*> Q;
		
		for(size_t q=0; q<n_hcs; q++){
			Q.push_back(setup_from_list(hcsParams[q], D, seedRNG));
			if(trace_level>0) Q[q]->print_ward_room();
			Q[q]->stats_architecture();
		}
		cout << " ... HCS built & populated (seed="<< seedRNG << ")." << endl;
		
		// Setup simulator:
		Simulator S(Q, horizon, timestep, recordContacts);
		cout << " ... Simulator setup (seed="<< seedRNG << ")." << endl;
		
		// Setup movements:
		vector<vector<float> > mvt_HCS;
		for(int i=0; i < mvtParams.size(); i++)
			mvt_HCS.push_back(mvtParams[i]);
		
		S.set_mvt_HCS(mvt_HCS);
		cout << " ... movements set (seed="<< seedRNG << ")." << endl;
		if(trace_level>0){
			cout <<" Movements b/w HCS set to:";
			displayVector(mvt_HCS);
		}
		
		// Run the simulation:
		S.run_simulation(seedRNG);
		
		// Retrieve results through census.
		// Each element of the vector is a time step:
		vector<censusDF> df;
		df = lightCensus? S.get_census_light():S.get_census();
		
		// === Outputs ===
		
		Rcpp::List rcpplist;
		vector<string> list_elem_name;
		
		// Population clinical database:
		
		for(uint i=0; i<df.size(); i++){
			
			vector<uint> uids = df[i]._uid_individual;
			size_t n = uids.size();
			vector<float> times(n, df[i]._census_time);
			
			// List is split because it cannot take more than 20 elements!!!
			// see: http://lists.r-forge.r-project.org/pipermail/rcpp-devel/2013-April/005700.html
			//
			Rcpp::List L = List::create(Named("time")	= times,
										Named("uid")	= uids,
										//Named("hcs_name")= df[i]._hcs_name,
										Named("type")	= df[i]._type_individual,
										//Named("subtype")= df[i]._subtype_individual,
										//										Named("isAlive")= df[i]._isAlive,
										//										Named("willDie")= df[i]._willDie,
										//Named("onDuty") 			= df[i]._onDuty,
										Named("sIdx")	   			= df[i]._sIdx,
										Named("infectIdx")	        = df[i]._infectIdx,
										Named("timeAdmission")		= df[i]._time_admission,
										Named("DoS")				= df[i]._DoS,
										Named("current_room_uid")	= df[i]._current_room_uid,
										Named("current_room_contamIdx") = df[i]._current_room_contamIdx,
										Named("wardAssigned")		= df[i]._wardAssigned,
										Named("isColonized")		= df[i]._isColonized,
										Named("wasColonized")		= df[i]._wasColonized,
										Named("symptomatic")		= df[i]._symptomatic,
										Named("currentlySymptomatic") = df[i]._currentlySymptomatic,
										Named("contamIdx")			= df[i]._contamIdx,
										Named("contamIdx_other")	= df[i]._contamIdx_other
										);
			
			Rcpp::List L1 = List::create(Named("prdIncubation") 	= df[i]._prdIncubation,
										 Named("prdSymptom")		= df[i]._prdSymptom,
										 Named("isVax")				= df[i]._isVax,
										 Named("timeAcquisition")	= df[i]._timeAcquisition,
										 Named("timeOnset")			= df[i]._timeOnset,
										 Named("timeRelapse")		= df[i]._timeRelapse,
										 Named("timeIsolationStart")= df[i]._timeIsolationStart,
										 Named("timeIsolationEnd")	= df[i]._timeIsolationEnd,
										 Named("timeClearance")		= df[i]._timeClearance,
										 Named("nColonizations")	= df[i]._nColonizations,
										 Named("nInfections")		= df[i]._nInfections,
										 Named("nRelapses")			= df[i]._nRelapses,
										 Named("isoDurFirst")		= df[i]._isoDurFirst,
										 Named("isoDurRelapses")	= df[i]._isoDurRelapses,
										 Named("totDurSymptoms")	= df[i]._totDurSymptoms,
										 Named("colectomy")			= df[i]._colectomy,
										 Named("fmt")				= df[i]._fmt,
										 Named("acquisitionType")	= df[i]._acquisitionType
										 );
			
			L = Language("c",L,L1).eval();
			rcpplist.push_back(L);
			string tmp = "pop_" + to_string(i);
			list_elem_name.push_back(tmp);
		}
		
		// Contact durations:
		if(recordContacts){
			for(uint i=0; i<S.n_HCS(); i++){
				rcpplist.push_back(S.get_HCS()[i]->get_contact_HCW_patient_dur_rec());
				string hcs_name = S.get_HCS()[i]->get_name();
				list_elem_name.push_back("contact_duration_" + hcs_name);
			}
		}
		
		// set the name of all elements of the list:
		rcpplist.attr("names") = list_elem_name;
		
		return rcpplist;
	}
	
	catch (...){
		::Rf_error(">>>> my C++ exception (unknown reason) <<<<");
		return NULL;
	}
}









