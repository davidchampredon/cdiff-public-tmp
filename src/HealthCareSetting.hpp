//
//  HealthCareSetting.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef HealthCareSetting_hpp
#define HealthCareSetting_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include "Registry.hpp"
#include "Ward.hpp"
#include "globalvar.hpp"
#include "Disease.hpp"
#include "Vaccine.hpp"

using namespace std;

class HealthCareSetting{
	
protected:
	
	string		_name;
	string		_type;		// LTCF, HOSPITAL
	
	vector<Ward*> _ward;
	
	// HCW not assigned to
	// a specific ward:
	vector<HCW*>  _hcw;
	
	// Simulation time is recoded (allows convenient access):
	float	_current_time;
	float	_time_step;
	
	
	Registry	_registry; // TO DO : DELETE IF NOT USED
	
	Disease		_disease;
	
	// maximum UID among all individual.
	// (used to know the starting value
	// when adding new individuals)
	uint		_max_uid_indiv;
	
	vector<float>	_sIdx_newPatient_prm;
	vector<float>	_sIdx_HCW_prm;
	
	vector<float>	_DoS_newPatient_prm;
	
	float			_admission_rate;  // number of admission per day
	
	
	// Parameters for the distribution of the
	// number of daily patient visits by HCW:
	vector<float>	_lambda_prm;
	
	// Relationship between susceptibility index and shedding index:
	float	_contamIdx_suscept_prm;
	
	// Daily reset level of contamIdx for HCWs:
	float	_reset_contamIdx;
	
	// Decay shedding after symptoms end:
	float	_decay_contamIdx_halftime;
	
	// Environmental contamination, other than patient rooms:
	float	_envContam_other;
	// Ratio of environmental contamination
	// other / patient room:
	float	_ratio_contam_other_room;
	
	// C. difficile prevalence for new admissions from the  community:
	float	_prev_new_admission;
	
	// Probability patients visiting "other" places per day:
	float	_proba_patient_visit_other;
	
	// Parameters for the distribution of the
	// contact duration between HCW and patients:
	float	_mean_contact_HCW_patient_duration_minutes;
	float	_var_contact_HCW_patient_duration_minutes;
	// Record of all contact durations:
	vector<float> _contact_HCW_patient_dur_rec;
	
	// Time to transmit half of one's contamination index
	// during a contact b/w individuals:
	float	_contact_indiv_halftime;
	
	// Proportion of contamIdx transmitted from
	// one individual to another during contact:
	float	_contact_xchg_contamIdx;
	
	// Proportion of contamIdx transmitted from
	// one room to an individual during visit:
	float	_visit_xchg_contamIdx;
	
	
	// Time to transmit half of the contamination index
	// of one romm during a visit by individual:
	float	_visit_room_halftime;
	
	// Room cleaning:
	float	_clean_room_efficacy;
	
	// Symptoms onset to isolation lag (mean)
	float	_onset_isolation_lag_mean;

	// Severe outcomes:
	float	_proba_colectomy;
	float	_proba_fmt;
	
	// Vaccination
	Vaccine	_vaccine;
	string	_vaxStrategy;
	float	_vax_start_time;
	float	_proba_vax_newAdmission;
	float	_sIdxMin_vax_frailty;  // minimum level of susceptibility index  for a patient to be vaccinated under the strategy "frailty"
	float	_DoSMin_vax_frailty;   // minimum duration of stay for a patient to be vaccinated under the strategy "frailty"
	
public:
	
	// === CONSTRUCTORS ===
	
	HealthCareSetting(){
		_name = "undefined";
	}
	
	uint generate_new_uid_indiv(){
		uint tmp = _max_uid_indiv;
		_max_uid_indiv++;
		return tmp;
	}
	
	void setup(string hcs_name,
			   string hcs_type,
			   uint n_iso_rooms,
			   uint n_ward_mean,
			   uint n_room_per_ward_mean,
			   uint n_beds_per_room_mean,
			   uint  nNurses,
			   uint  nDoctors,
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
			   float proba_fmt
			   );
	
	// === GET FUNCTIONS ===
	
	string			get_name(){return _name;}
	
	vector<Ward*>	get_ward(){return _ward;}
	vector<HCW*>	get_hcw(){return _hcw;}
	
	vector<float>	get_contact_HCW_patient_dur_rec(){return _contact_HCW_patient_dur_rec;}
	float			get_contamIdx_suscept_prm(){return _contamIdx_suscept_prm;}
	float			get_decay_contamIdx_halftime(){return _decay_contamIdx_halftime;}
	
	float			get_current_time(){return _current_time;}
	
	Vaccine			get_vaccine(){return _vaccine;}
	float			get_proba_vax_newAdmission(){return _proba_vax_newAdmission;}
	float			get_sIdxMin_vax_frailty(){return _sIdxMin_vax_frailty;}
	float			get_DoSMin_vax_frailty(){return _DoSMin_vax_frailty;}
	
	float			get_onset_isolation_lag_mean(){return _onset_isolation_lag_mean;}
	
	string			get_vaxStrategy(){return _vaxStrategy;}
	float			get_vax_start_time() {return _vax_start_time;}
	
	// === SET FUNCTIONS ===
	
	void set_current_time(float x){_current_time = x;}
	void set_time_step(float x)   {_time_step = x;}
	
	void set_disease(Disease D){_disease = D;}
	
	void set_sIdx_newPatient_prm(vector<float> x){
		for(uint i=0;i<x.size(); i++)_sIdx_newPatient_prm.push_back(x[i]);
	}
	void set_sIdx_HCW_prm(vector<float> x){
		for(uint i=0;i<x.size(); i++)_sIdx_HCW_prm.push_back(x[i]);
	}
	void set_DoS_newPatient_prm(vector<float> x){
		for(uint i=0;i<x.size(); i++)_DoS_newPatient_prm.push_back(x[i]);
	}
	void set_admission_rate(float x){_admission_rate = x;}
	void set_lambda_prm(vector<float> x){_lambda_prm = x;}
	void set_envContam_other(float x) {_envContam_other = x;}
	void set_max_uid_indiv(uint x) {_max_uid_indiv = x;}
	void set_clean_room_efficacy(float x) {_clean_room_efficacy = x;}
	
	void set_vaccine(Vaccine x) {_vaccine = x;}
	void set_vaxStrategy(string x){_vaxStrategy = x;}
	void set_proba_vax_newAdmission(float x){_proba_vax_newAdmission = x;}
	void set_sIdxMin_vax_frailty(float x){_sIdxMin_vax_frailty = x;}
	void set_DoSMin_vax_frailty(float x){_DoSMin_vax_frailty = x;}
	
	void set_vax_start_time(float x){_vax_start_time = x;}
	
	void set_onset_isolation_lag_mean(float x){_onset_isolation_lag_mean = x;}
	
	void set_proba_colectomy(float x){_proba_colectomy = x;}
	void set_proba_fmt(float x){_proba_fmt = x;}
	
	// === ARCHITECTURE ===
	
	/** Build empty HCS from description vectors
	 *  (meant to be used with R wrap)*/
	void build_HCS(string hcs_name,
					 string hcs_type,
					 vector<string> room_type,
					 vector<uint> ward,
					 vector<string> room_name,
					 vector<uint> max_patient,
				   float decay_room_halftime);
	
	void build_HCS_random(string hcs_name,
						  string hcs_type,
						  uint n_iso_rooms,
						  uint n_ward_mean,
						  uint n_room_per_ward_mean,
						  uint n_beds_per_room_mean,
						  float decay_room_halftime);
	
	void stats_architecture();
	
	/// Retrieve the total number of wards.
	uint total_number_of_wards();
	
	/** Retrieve the maximum number of patients per ward. 
	 * (Used to assign proportional number of nurses to wards) */
	vector<uint> max_patient_per_ward();

	
	/** Retrieve the total maximum of number of patients this
	 * HCS can have. */
	uint total_max_patients();
	
	
	/// Get room given by its POSITION in the vector '_room' in a given ward.
	Room* get_room_by_pos(uint room_pos, uint ward_uid);
	
	/// Get room of a given uid, in a given ward.
	Room* get_room_by_uid(uint room_uid, uint ward_uid);
	
	/// Add a free HCW (not assigned to a specific ward).
	void add_hcw(HCW* hcw) {_hcw.push_back(hcw);}
	
	/// Print out the ward/room structure.
	void print_ward_room();
	
	
	// === CLINICAL ===
	
	float draw_DoS();
	float draw_sIdx(string indiv_type="patient");
	
	/// Draw a colectomy event among all symptomatic patients.
	void draw_colectomy();
	
	/// Draw a FMT event from all patients that had at least one relapse.
	void draw_fmt();
	
	// === MOVEMENTS ===
	
	/// Receive a patient transfered from another HCS.
	void receive_transfered_patient(Patient* p,
									float DoS,
									float contamIdx_suscept_prm,
									float decay_contamIdx_halftime,
									Room* r,
									float time);
	
	/// Admission of ONE new patient only. Deals with vaccination at admission.
	void admission_new_patient(float sIdx,
							   float DoS,
							   float contamIdx_suscept_prm,
							   float decay_contamIdx_halftime,
							   Room* r,
							   float time);
	
	/// Discharge the patient at position 'i_position' in ward 'i_ward', room index 'i_room'.
	void discharge_patient(uint i_ward, uint i_room, uint i_position);
	
	/// Discharge the patient at position 'i_position' in ward 'i_ward', room index 'i_room'.
	void discharge_patient_by_uid(uint i_ward, uint i_room, uint patient_uid);
	
	/** Discharge all patients who have reached
	 * the end of their (pre-specified) duration of stay. */
	void discharge_all();
	
	/** Admission of all new patients during that time step. */
	void admission_new_all();
	
	/** Populate randomly patients across all 'patient_rooms'.
	 Typically used when HCS is newly built. */
	void populate_patient(float occupancy);
	
	/** Populate HCW  of type "hcw_type" 
	 * in rooms of type "roomType" of the HCS. */
	void populate_hcw_roomType(string roomType, string hcwType, uint n);
	
	/** Populate HCW in the HCS. */
	void populate_hcw(uint nHCW_ward, uint nHCW_free);
	
	/** Return the _indexes_ position of the available rooms.
	 *  Returned object is a vector of vector, where columns
	 *  represent the ward and rows the rooms, that is:
	 *  result[ward][room].
	 */
	vector< vector<uint> > find_available_room();
	
	/// Adjust duration of stay for all patients according to their symptomatic period.
	void adjust_DoS();
	
	/// Move patient 'P' to new room 'r'.
	void move_patient(Patient* P, Room* r);
	
	/// Move given patients in isolation room. (NOTE: The isolation ward is a global variable)
	void isolate(vector<Patient*> P);
	
	/// Ends isolation of all patients whose symptoms ended 'd' days ago. Move those patient in a normal room.
	void end_isolation();
	
	/// Retireve patients that must end isolation ASAP.
	vector<Patient* > retrieve_end_isolation_ASAP();
	
	/** Returns the 3 _indexes_ of the ward and room and patient's UID for patients that are transferables.
	 *  Each row of the returned matrix has such a triplet of indexes:
	 *  [w1, r1, puid1]  <-- first transferable patient
	 *  [w2, r2, puid2]  <-- 2nd   transferable patient
	 *  [etc.]
	 * 
	 * For now, criteria is end of stay reached. It's like discharge,
	 * SO MUST BE RUN BEFORE DISCHARGE!!!!
	 */
	vector<vector<uint> > find_transferable_patients();
	
	/// Find the index of non-empty wards (exclude isolation ward).
	vector<uint> find_non_empty_ward();
	
	/// Select randomly a non-empty ward (exclude isolation ward).
	Ward* select_random_non_empty_ward();

	
	// === VACCINATION ====
	
	/// Vaccinate individuals according to vaccination strategy.
	void vaccinate();
	
	/// Vaccinate one individual.
	void vaccinate_individual(Patient* p);
	
	
	/// Vaccinate new admissions randomly with probability p.
	void vax_newAdm_random(float proba);
	
	/// Vaccinate new admission base on their
	/// suscpetibility index (i.e., frailty) and
	/// duration of stay, with probability p.
	void vax_newAdm_frailty(float sIdxmin, float DoSmin, float p);
	
	/// Update the infection index of all vaccinated patients.
	void update_all_infectIdx_vax();
	
	// === EPIDEMIOLOGY ===
	
	/// Draw the number of patient visits a single HCW will do during 'dt'.
	uint draw_n_visits_HCW_patient(float dt);
	
	/// Draw the 'n' patients that will be contacted by a HCW.
	vector<Patient*> draw_patients_contacted_by_HCW(HCW* hcw,
													uint n);
	
	/// Contacts made by one HCW.
	void contacts_one_HCW(HCW* hcw, bool recordContacts);
	
	/** Draw all contacts of all HCW during 'dt'.
	 * These contacts trigger:
	 * 1) transmission between patient and HCW,
	 * 2) shedding in room,
	 * 3) acquisition from room.
	 */
	void draw_all_contacts_patient_HCW(bool recordContacts);
	
	/// Seed the initial patient population with infection.
	void seed_initial_infection(float prevalence_percent);
	
	/// Update clinical variables for ALL individuals.
	void update_all_clinical();
	
	/// Update the contamination decay in all patient rooms.
	void update_decay_contamIdx_patient_rooms();
	
	/// Update the level of environmental contamination in other places than patient rooms.
	void update_envContam_other();
	
	/// Update the duration of isolation for all patients
	void update_isoDur();
	
	/// Update the total duration of symptoms for all patients
	void update_totDurSymptoms();
	
	
	/// Draw acquisition from environmental contamination other than patient room, during dt.
	void draw_acq_envContam_other();
	
	/// Reset contamIdx to a very small level for one HCW only.
	void reset_contamIdx_one_HCW(HCW* hcw,
								 float reset_level);
	
	/// Reset contamIdx to a very small level for all HCW. (meant to reset daily).
	void reset_contamIdx_all_HCW();
	
	/// Environmental contamination (shedding & acquisition) from _one single_ patient in its own room.
	void one_patient_envContam_own_room(Patient* p,
										float curr_time);
	
	/** Environmental contamination (shedding & acquisition) 
	 from _all_ patients in their own room.*/
	void all_patients_envContam_own_room();
	
	
	/// Retrieve all patients that have just started symptoms onset.
	vector<Patient*> retrieve_onset_patients(float curr_time,
											 float timestep);
	
	/// Draw the lag between onset time and isolation time
	/// for all patients starting onset of disease.
	void draw_onset_iso_lag();
	
	/// Retrieve patients that have a drawn
	/// isolation time in the current time step.
	vector<Patient *> retrieve_patients_to_isolate(float curr_time,
												   float timestep);
	
	/// Start isolation of patients that have a drawn
	/// isolation time in the current time step.
	void start_isolation();
	
	/// Retrieve patients, from an isolation ward only, that have had their symptoms resolved 'd' days ago.
	vector<Patient*> retrieve_resolution_isopatients(float curr_time,
												  float timestep,
												  float d);
	
	
	/// Return patients that are candidates for a relapse.
	vector<Patient*> draw_patients_relapse(float proba, float shape);
	
	/// Initiate relapses among patients.
	void relapses();
	
	
	// === EVENTS ===
	
	/// Check if 't' is in the current simulation time step.
	bool is_in_current_timestep(float t);
	
	/** Only ONE new patient is admitted (if room available) */
	void event_new_admission(float time);
	
	/// Scan all rooms for extra cleaning.
	void extra_cleaning();
	
	// === CENSUS ===
	
	vector<uint>	census_uid_individual;
	vector<float>	census_sIdx;
	vector<float>	census_infectIdx;
	vector<string>	census_type_individual;
	vector<string>	census_subtype_individual;
	
	vector<bool>	census_isAlive;
	vector<bool>	census_willDie;
	
	vector<bool>	census_colectomy;
	vector<bool>	census_fmt;
	
	vector<float>	census_DoS;
	vector<float>	census_timeAdmission;
	vector<uint>	census_current_room_uid;
	vector<float>	census_current_room_contamIdx;
	vector<uint>	census_wardAssigned;
	
	vector<bool>	census_isColonized;
	vector<bool>	census_wasColonized;
	vector<bool>	census_symptomatic;
	vector<bool>	census_currentlySymptomatic;
	vector<float>	census_prdIncubation;
	vector<float>	census_prdSymptom;
	
	vector<float>	census_contamIdx;
	vector<float>	census_contamIdx_other;
	
	vector<bool>	census_onDuty;
	vector<string>	census_hcs_name;
	
	vector<bool>	census_isVax;
	
	vector<float>	census_timeAcquisition;
	vector<float>	census_timeOnset;
	vector<float>	census_timeIsolationStart;
	vector<float>	census_timeIsolationEnd;
	vector<float>	census_timeClearance;
	vector<float>	census_timeRelapse;
	
	vector<float>	census_isoDurFirst;
	vector<float>	census_isoDurRelapses;
	
	vector<float>	census_totDurSymptoms;
	
	vector<string>  census_acquisitionType;
	
	vector<uint>	census_nColonizations;
	vector<uint>	census_nInfections;
	vector<uint>	census_nRelapses;
	
	/// Reset and clear all census vectors.
	void 	census_clear_all();
	
	void	census_one_individual(Individual* indiv);
	void	census_one_patient(Patient* p);
	void	census_one_hcw(HCW* hcw);
	
	/// Exhaustive census, at all times.
	void	census_individuals();
	
	/// Census for patient only, at admission and discharge times.
	void	census_patients_light();
	
	
	
	// === MISCELLENAOUS ===
	
	uint capacity_patients();
	
	uint count_patients();
	uint count_hcw_ward();
	uint count_hcw_free();

	
	/// Retrieve ALL patients (from all wards).
	vector<Patient*> retrieve_all_patients();
	
	/// Draw randomly 'n' patients.
	vector<Patient*> draw_patients(uint n);
	
	/// Return patient's current room.
	Room* patient_room(Patient* P);
	
	
	void show(bool details = false);
};




#endif /* HealthCareSetting_hpp */



