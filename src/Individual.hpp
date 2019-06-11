//
//  Individual.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef Individual_hpp
#define Individual_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>

#include "globalvar.hpp"
#include "Disease.hpp"
#include "utils.hpp"
#include "Vaccine.hpp"

using namespace std;


class Individual{
	
protected:
	
    uint    _uid;          // Unique identifier
    
    string  _type;              // type of individual e.g., patient, HCW
    string  _subtype;           // subtype for HCW e.g., nurse, doctor
    
    uint    _current_room_uid;  // UID of the room where this Individual currently is.
	uint	_current_ward_id;
    
	// === CLINICAL VARIABLES ===
	
    float   _susceptIdx;        // susceptibility index
	float   _susceptIdx_baseline;  // baseline susceptibility index (used when vaccinated)
	
	float	_infectIdx;			// infection index, i.e., proba to be infected conditional on being colonized
	float	_infectIdx_baseline;
	
    float   _prdIncubation;     // Incubation period (drawn at acquisition time)
    float   _prdSymptom;        // Symptomatic duration (if symptomatic infection)
    
    // period before clearance.
    // Warning: from acquisition for asymptomatic, but
    // from end of symptoms for symptomatic infections.
    float   _prdClearance;
    
    // Skin (and clothes and instruments) contamination index.
    // Represent probability to infect someone during one contact.
    float   _contamIdx;
    
    // Parameter defining the relationship b/w
    // shedding and susceptibility:
    float   _contamIdx_suscept_prm;
    
    // Decay of shedding after symptoms end:
    float   _decay_contamIdx_halftime;
    
    bool    _isColonized;
	bool    _wasColonized;
    bool    _symptomatic;
    bool    _currentlySymptomatic;
    bool    _recovering;
    
    
    // === EPIDEMIOLOGICAL VARIABLES ===
    
    float   _timeAcquisition;   // Time when C. diff acquired
	float	_timeOnset;			// Time of symptoms onset
	
	string	_acquisitionType;	// type of acquisition: "NA", "human" or "environment"
	
    bool    _isVax;             // Vaccinated or not
	uint	_nColonizations;	// Number of colonizations
	uint	_nInfections;		// Number of infections
	uint	_nRelapses;			// Number of relapses
	
	bool	_endIsolationASAP;	// flag when end isolation not possible bc other rooms full.
	
	
    
public:
    
    // === CONSTRUCTORS ===
    
    Individual() {baseConstructor();}
    
    void baseConstructor();
    
    
    // === GET FUNCTIONS ===
    
    string  get_type(){return _type;}
    string  get_subtype(){return _subtype;}
    uint    get_uid(){return _uid;}
	
	float   get_susceptIdx(){return _susceptIdx;}
	float   get_susceptIdxBaseline(){return _susceptIdx_baseline;}
	float	get_infectIdx(){return _infectIdx;}
	float	get_infectIdx_baseline(){return _infectIdx_baseline;}
	
    uint    get_current_room_uid(){return _current_room_uid;}
	uint    get_current_ward_uid(){return _current_ward_id;}
    
    bool    get_isColonized(){return _isColonized;}
	bool    get_wasColonized(){return _wasColonized;}
    bool    get_symptomatic(){return _symptomatic;}
    bool    get_currentlySymptomatic(){return _currentlySymptomatic;}
    
    float   get_contamIdx(){return _contamIdx;}
    float   get_contamIdx_suscept_prm(){return _contamIdx_suscept_prm;}
    
    float   get_prdIncubation(){return _prdIncubation;}
    float   get_prdSymptom(){return _prdSymptom;}
    float   get_timeAcquisition(){return _timeAcquisition;}
	float   get_timeOnset(){return _timeOnset;}
    
	uint	get_nColonizations(){return _nColonizations;}
	uint	get_nInfections(){return _nInfections;}
	uint	get_nRelapses(){return _nRelapses;}
	
	bool	get_endIsolationASAP(){return _endIsolationASAP;}
	
	string	get_acquisitionType(){return _acquisitionType;}
	
    
    // === SET FUNCTIONS ===
    
    void set_uid(uint uid){_uid = uid;}
    void set_type(string x) {stopif(x!="patient" && x!="HCW","unknown type"); _type = x;}
    void set_subtype(string x) {_subtype = x;}
    void set_susceptIdx(float x){_susceptIdx = x;}
	void set_susceptIdxBaseline(float x){_susceptIdx_baseline = x;}
	void set_infectIdx(float x){_infectIdx = x;}
	void set_infectIdx_baseline(float x){_infectIdx_baseline = x;}
    void set_current_room_uid(uint x){_current_room_uid = x;}
	void set_current_ward_uid(uint x){_current_ward_id = x;}
    void set_contamIdx(float x){_contamIdx = x;}
    void set_contamIdx_suscept_prm(float x){_contamIdx_suscept_prm = x;}
    void set_decay_contamIdx_halftime(float x){_decay_contamIdx_halftime = x;}
	void set_endIsolationASAP(bool x){_endIsolationASAP = x;}
	void set_acquisitionType(string x){_acquisitionType = x;}
	void set_currentlySymptomatic(bool x){_currentlySymptomatic = x;}
	
	void incr_nInfections(){_nInfections++;}
    
    // === EPIDEMIOLOGY ===
	
	/// Return time of symptoms onset.
	void calc_timeOnset();
	
    /// Colonize this individual with C. difficile.
    void colonize(float time, Disease D);
	
	/// Draw incubation period value.
	void draw_prdIncubation(Disease D);
	
	/// Draw the symptomatic period value.
	void draw_symptomatic_period(Disease D);
	
	/// Draw the clearance period.
	void draw_prdClearance(Disease D);
	
	/// Symptomatic infection triggered.
	void symptomatic_infection_triggered(Disease D);

	
    /// Clearance of colonization:
    void clearance();
	
	
	void update_clinical_sympt(float curr_time,
							   float t_sympt_start,
							   float t_sympt_end,
							   Disease D);
	
	void update_clinical_clear(float curr_time);
	
    /// Update all clinical variables
    void update_clinical(float curr_time, Disease D);
	
	/// Increase _contamIdx by 'x'. Check <1.
	void incr_contamIdx(float x);
	
	/// Decrease _contamIdx by 'x'. Check >0.
	void decr_contamIdx(float x);
	
	/// Return time of symptoms resolution.
	float time_resolution();
	
	/// Return time of pathogen clearance.
	float time_clearance();
	
    // === MISCELLEANOUS ===
    
    void    show();
	
};


class Patient: public Individual{
    
protected:
    
    // === CLINICAL VARIABLES ===
    
    bool    _isAlive;
	bool 	_colectomy;		// does this patient have colectomy resulting from CDI?
	bool 	_fmt;			// does this patient have Fecal Microbial Transplant resulting from CDI?
	
	
	// === VACCINATION VARIABLES ===
	
	bool	_isVax;		// is this patient vaccinated?
	uint	_nVaxDoses;	// how many vaccine doses patient received over her/his life.
	float	_timeVax;	// time when patient received her/his _last_ vaccine dose.

    // === HOSPITALIZATION VARIABLES ===
    
    float   _timeAdmission;
    float   _DoS;               // duration of stay
    bool    _willDie;           // Will patient die at end of hospital stay (drawn at admission)?

	float 	_timeIsolationStart;
	float 	_timeIsolationEnd;
	float	_timeLastResolution; // time of last resolution (when nRelapse>0)
	float 	_timeRelapse;

	float	_isoDurFirst;		// Isolation duration for first CDI episode
	float	_isoDurRelapses;	// Isolation duration for all CDI relapses after initial episode
	
	bool	_isIsolated;
	
	float	_totDurSymptoms;		// Total duration of symptoms
	
	
public:
    
    Patient(): Individual(){
        
        _type = "patient";
        
        _timeAdmission = -NA_FLOAT;
        _DoS        = 0;
        _isAlive    = true;
        _willDie    = false;
		
		_colectomy 	= false;
		_fmt		= false;
		
		_isVax		= false;
		_nVaxDoses	= 0;
		_timeVax	= -NA_FLOAT;
		
		_timeIsolationStart = -NA_FLOAT;
		_timeIsolationEnd	= +NA_FLOAT;
		_timeLastResolution = -NA_FLOAT;
		_timeRelapse		= NA_FLOAT;
		
		_isIsolated = false;
		_isoDurFirst = 0.0;
		_isoDurRelapses = 0.0;
		
		_totDurSymptoms = 0.0;
    }
    
    void 	set_timeAdmission(float x){_timeAdmission = x;}
	void	set_colectomy(bool x) {_colectomy = x;}
	void 	set_fmt(bool x){_fmt = x;}
    void 	set_DoS(float x){_DoS = x;}
	void 	set_isVax(bool x){_isVax = x;}
	void 	set_nVaxDoses(uint x){_nVaxDoses = x;}
	void 	set_timeVax(float x){_timeVax = x;}
	void 	set_timeIsolationStart(float x){_timeIsolationStart = x;}
	void 	set_timeIsolationEnd(float x){_timeIsolationEnd = x;}
	void	set_isIsolated(bool x){_isIsolated = x;}
	
	void	incr_isoDurFirst(float x){_isoDurFirst += x;}
	void	incr_isoDurRelapses(float x){_isoDurRelapses += x;}
	
	void	incr_totDurSymptoms(float x){_totDurSymptoms += x;}
	
	bool	get_isAlive(){return _isAlive;}
	bool	get_willDie(){return _willDie;}
	float	get_timeAdmission(){return _timeAdmission;}
    float	get_DoS(){return _DoS;}
	bool	get_isVax(){return _isVax;}
	uint	get_nVaxDoses(){return _nVaxDoses;}
	float	get_timeVax(){return _timeVax;}
	float   get_timeIsolationStart(){return _timeIsolationStart;}
	float   get_timeIsolationEnd(){return _timeIsolationEnd;}
	float	get_timeLastResolution(){return _timeLastResolution;}
	float	get_timeRelapse(){return _timeRelapse;}
	bool	get_colectomy(){return _colectomy;}
	bool	get_fmt(){return _fmt;}
	
	bool	get_isIsolated(){return _isIsolated;}
	float	get_isoDurFirst(){return _isoDurFirst;}
	float	get_isoDurRelapses(){return _isoDurRelapses;}
	float	get_totDurSymptoms(){return _totDurSymptoms;}
	
	float	get_timeDischarge(){return _timeAdmission + _DoS;}
	
	/// Triggers relapse (i.e., set new symptomatic period)
	void relapse_triggered(float curr_time, Disease D);
	
	/// Update all clinical variables.
	void update_clinical_patient(float curr_time, Disease D);
	
	/// Update infection index if this individual is vaccinated.
	void 	update_infectIdx_vax(float curr_time, const Vaccine &V);
	
	/// Vaccinate one individual.
	void	vaccinate(float t, const Vaccine &V);
	
    /// Adjust DoS according to symptomatic C. difficile infection
    void	adjust_DoS(float curr_time);
	
	/// Adjust DoS after relapse
	void	adjust_DoS_relapse(float curr_time);
	
	/// Duration of isolation
	float	isolationDuration();
	
    void show();
};


class HCW: public Individual{
    
protected:
    
    int     _wardAssigned;      // Ward assigned for nurses
	
	vector<uint> _periodsWorked; // Periods of the day that are worked by this HCW.
    
public:
    
    // === CONSTRUCTOR ===
    
    HCW():Individual(){
		_type = "HCW";
		_wardAssigned = NA_UINT;
	}
	
    HCW(string subtype):Individual(){
        _type = "HCW";
        _subtype = subtype;
        _wardAssigned = NA_UINT;
    }
    
    // === SET FUNCTIONS ===
    
    void set_wardAssigned(uint x){_wardAssigned = x;}
	void set_periodsWorked(vector<uint> x);
    
    
    // === GET FUNCTIONS ===
    
    uint get_wardAssigned(){return _wardAssigned;}
    
	
	// === Miscellenaous ===
	
	/// Whether this HCW is on duty or not.
	bool onDuty(float curr_time, float timestep);
	
    void show();
};


#endif /* Individual_hpp */


