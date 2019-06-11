//
//  Individual.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include <random>

#include "Individual.hpp"
//#include "globalvar.hpp"


void Individual::baseConstructor(){
    
    _uid        = 0;
    
    _type       = "undefined";
    _subtype    = "undefined";
    
    _current_room_uid = NA_UINT;
    
    _susceptIdx = -NA_FLOAT;
	_susceptIdx_baseline = -NA_FLOAT;
	
	_infectIdx = -NA_FLOAT;
	_infectIdx_baseline = -NA_FLOAT;
    
    _prdSymptom     = 0.0;
    _prdIncubation  = 0.0;
    _prdClearance   = NA_FLOAT;
    
    _isColonized    = false;
	_wasColonized   = false;
    _symptomatic    = false;
    _currentlySymptomatic = false;
    _recovering     = false;
    
    _contamIdx    = 0.0;
    
    // === EPIDEMIOLOGICAL VARIABLES ===
    
	_timeAcquisition= -NA_FLOAT; // TO DO: change to +NA_FLOAT and check it's ok...
	_timeOnset		= NA_FLOAT;
    _isVax  		= false;
	_nColonizations	= 0;
	_nInfections	= 0;
	_nRelapses		= 0;
	
	_endIsolationASAP = false;
	
	_acquisitionType = "NA";
	
}

void Individual::show(){
    cout << endl << "Individual #" << _uid << ":" << endl;
    cout << " Type: " << _type << endl;
    cout << " Subtype: " << _subtype << endl;
    cout << " Current room UID: " << _current_room_uid << endl;
    cout << " susceptIdx: " << _susceptIdx << endl;
	cout << " infectIdx: " << _infectIdx << endl;
    cout << " contamIdx: "<< _contamIdx << endl;
}

void Patient::show(){
    Individual::show();
    cout << " Alive: " << _isAlive << endl;
    cout << " DoS: "  << _DoS << endl;
}

void HCW::show(){
    Individual::show();
    cout <<" Ward assigned: "<< _wardAssigned << endl;
}

void Individual::draw_prdIncubation(Disease D){
	
	string distrib_type = "gamma";
	
	if(distrib_type == "lnorm"){
		float mu  = log(D.get_prd_incubation_mean());
		float sig = sqrt(D.get_prd_incubation_var());
		lognormal_distribution<float> dlogn(mu, sig);
		_prdIncubation = dlogn(RANDOM_GENERATOR);
	}
	
	if(distrib_type == "gamma"){
		float m = D.get_prd_incubation_mean();
		float v = D.get_prd_incubation_var();
		float b = v / m;
		float a = m / b;
		
		stopif(abs(m)<1e-2, "draw_prdIncubation() mean~0 ! ");
		stopif(abs(b)<1e-9, "draw_prdIncubation() b~0 ! ");
		stopif( (m<=0 || v <=0), "draw_prdIncubation() m or v <0 ! ");
		
		gamma_distribution<float> dgam(a, b);
		_prdIncubation = dgam(RANDOM_GENERATOR);
	}
	// put a hard cap on the value:(TO DO: CHANGE THAT)
	if(_prdIncubation>365.0) _prdIncubation=365.0f;
}


void Individual::draw_symptomatic_period(Disease D){
	
	string distrib_type = "exp"; // <--- TO DO: do not hard code, put in prm CSV file
	
	if(distrib_type == "lnorm"){
		float mu  = log(D.get_prd_symptom_mean());
		float sig = sqrt(D.get_prd_symptom_var());
		lognormal_distribution<float> dlogn(mu, sig);
		_prdSymptom = dlogn(RANDOM_GENERATOR);
	}
	
	if(distrib_type == "gamma"){
		float m = D.get_prd_symptom_mean();
		float v = D.get_prd_symptom_var();
		float b = v / m;
		float a = m / b;
		gamma_distribution<float> dgam(a, b);
		_prdSymptom = dgam(RANDOM_GENERATOR);
	}
	
	if(distrib_type == "exp"){
		float m = D.get_prd_symptom_mean();
		exponential_distribution<float> dexp(1/m);
		_prdSymptom = dexp(RANDOM_GENERATOR);
	}
	
	
	// Hard coded cap. TO DO: CHANGE THAT
	// (very large number were generated: investigate!)
	if(_prdSymptom<1.1)	_prdSymptom = 1.1;
	if(_prdSymptom>60)	_prdSymptom = 60;
}


void Individual::symptomatic_infection_triggered(Disease D){
	
//	cout<<"DEBUG: uid:"<<_uid<<" is infected";
//	cout <<" ; time = "<<_cu
	
	_nInfections++;
	draw_prdIncubation(D);
	calc_timeOnset();
	draw_symptomatic_period(D);
}

void Patient::relapse_triggered(float curr_time, Disease D){
	_nRelapses ++;
	incr_nInfections();
	_timeRelapse = curr_time;
	_prdIncubation = 0;
	_symptomatic = true;
	set_currentlySymptomatic(true);
	_timeLastResolution = _timeIsolationEnd; // previous episode
	_timeOnset = curr_time;
	draw_symptomatic_period(D);
	draw_prdClearance(D);
}

void Individual::draw_prdClearance(Disease D){
	float m   = D.get_prd_clear_asymptom_mean();
	float var = D.get_prd_clear_asymptom_var();
	if(_symptomatic){
		m   = D.get_prd_clear_symptom_mean();
		var = D.get_prd_clear_symptom_var();
	}
	lognormal_distribution<float> dlnorm(log(m),var);
	_prdClearance = dlnorm(RANDOM_GENERATOR);
	if(_prdClearance<0) _prdClearance = 0.01;
}


void Individual::colonize(float time, Disease D){
    
    _isColonized = true;
	_wasColonized = true;
    _timeAcquisition = time;
	_nColonizations++;
    // Assume only patient can be symptomatic.
    if(_type=="patient")
    {
        // Symptomatic infection or not:
        bernoulli_distribution dbern(_infectIdx);
        _symptomatic = dbern(RANDOM_GENERATOR);
        // Incubation, symptomatic periods and number of episodes:
		if(_symptomatic) {
			symptomatic_infection_triggered(D);
//			cout<<"DEBUG: time:"<<time<<" ; uid:"<<_uid;
//			cout <<" is colonized and onset will be: ";
//			cout << time + _prdIncubation ;
//			cout << " and prdSympt = "<< _prdSymptom;
//			cout << " ; time resolution: " << time_resolution() << endl;
		}
    }
	
    // Clearance:
    float proba_clear = pow( 1.0 - _susceptIdx, D.get_proba_clear_prm() );
    bernoulli_distribution dbern_clear(proba_clear);
    bool will_clear = dbern_clear(RANDOM_GENERATOR);
	
    if(will_clear) draw_prdClearance(D);
	
    // Contamination index:
    uniform_real_distribution<float> U(0.0, D.get_contamIdx_max_asymptom() );
	_contamIdx = U(RANDOM_GENERATOR);
//    cout << "DEBUG : Individual #"<<_uid<< " has just been colonized, _contamIdx="<<_contamIdx << endl;
}


void Individual::clearance(){
    _isColonized    = false;
    _contamIdx      = 0.0;
	_recovering 	= false;
}

float Individual::time_resolution(){
	float x = NA_FLOAT;
	uint nr = get_nRelapses();
	
	if(nr==0) x = _timeAcquisition + _prdIncubation + _prdSymptom;
	if(nr>0)  x = _timeOnset + _prdSymptom;
	
	if(x >= NA_FLOAT){
		string msg1 = "time acquisition =" + to_string(_timeAcquisition);
		string msg2 = " prd incub =" + to_string(_prdIncubation);
		string msg3 = " prd sympt =" + to_string(_prdSymptom);
		string msg = msg1 + msg2 + msg3;
		stopif(true, msg);
	}
	return x;
}


void Individual::update_clinical_sympt(float curr_time,
									   float t_sympt_start,
									   float t_sympt_end,
									   Disease D){
	
	// During symptomatic period
	// (code should pass only once inside
	// this "if" condition):
	if(curr_time >= t_sympt_start &&
	   _symptomatic &&
	   !_currentlySymptomatic &&
	   !_recovering)
	{
		_currentlySymptomatic = true;
		
		// Contamination index (i.e., shedding) is
		// increased when symptomatic infection:
		float ratio = D.get_contamIdx_ratio_symptomatic();
		_contamIdx = min(1.0f, _contamIdx * ratio );
	}
	
	// After end of symptoms:
	if( t_sympt_end <= curr_time && _symptomatic )
	{
		_currentlySymptomatic = false;
		_recovering = true;
		float rate = -log(2)/_decay_contamIdx_halftime;
		float dt   = curr_time - t_sympt_end;
		_contamIdx = _contamIdx * exp( rate * dt );
	}
}


void Individual::update_clinical_clear(float curr_time){
	// Check if clearance occurs:
	float t_clear = time_clearance();
	if(t_clear <= curr_time) clearance();
}

void Individual::update_clinical(float curr_time, Disease D){
    if(_isColonized){
        // times symptoms start and end:
		float t_sympt_start = get_timeOnset();
		float t_sympt_end   = time_resolution();
		
		update_clinical_sympt(curr_time, t_sympt_start, t_sympt_end, D);
		update_clinical_clear(curr_time);
    }
}

void Patient::update_clinical_patient(float curr_time, Disease D){
	if(_isColonized){
		// times symptoms start and end:
		float t_sympt_start = get_timeOnset();
		if(_nRelapses>0) t_sympt_start = _timeRelapse;
		// TO DO: check this is always OK:
		float t_sympt_end   = time_resolution();
		
		update_clinical_sympt(curr_time, t_sympt_start, t_sympt_end, D);
		update_clinical_clear(curr_time);
	}
}


void Patient::vaccinate(float t, const Vaccine &V){
	//		cout << " DEBUG: UID" << _uid << "is vax as new admission"<<endl;
	_isVax		= true;
	_timeVax 	= t;
	_nVaxDoses++;
	
	// Vaccine directly affects infectivity index:
	float eff = V.draw_efficacy();
	set_infectIdx(_infectIdx * (1.0 - eff) );
	
	// Colonization prevention probability.
	// Setting the susceptibility index to 0
	// guarantees no colonization is possible (proba=0):
	float p = V.get_proba_prevent_colonization();
	if(p>0.0){
		bernoulli_distribution dbern(p);
		bool prevents_colonization = dbern(RANDOM_GENERATOR);
		if(prevents_colonization) set_susceptIdx(0.0);
	}
}


void Patient::update_infectIdx_vax(float curr_time, const Vaccine &V){
	if(_isVax){
		float dt	= curr_time - _timeVax;
		float rate	= -log(2) / V.get_halfLife();
		float eff	= V.draw_efficacy();
		float tmp1 	= 1.0 - exp(rate * dt);
		float idxb	= get_infectIdx_baseline();
		float tmp2 	= ( tmp1 * eff + (1.0 - eff) ) * idxb ;
		set_infectIdx( tmp2 );
	}
}




void Patient::adjust_DoS(float curr_time){
    // If DoS is shorter than symptomatic period,
    // then extend it:
    float tse = time_resolution(); // time symptoms end
    if(_timeAdmission + _DoS < tse){
        set_DoS(tse - _timeAdmission + 1.0); // <-- '+1' to make sure there is no discharge right at the moment of resolution.
    }
}

void Patient::adjust_DoS_relapse(float curr_time){
	if(get_nRelapses()>0){
		float t_end = curr_time + get_prdSymptom();
		if(_timeAdmission + _DoS < t_end){
			set_DoS(t_end - _timeAdmission + 1.0); // <-- '+1' to make sure there is no discharge right at the moment of resolution.
		}
	}
}


void Individual::incr_contamIdx(float x){
	float tmp = _contamIdx + x;
	if(tmp > 1.0) tmp = 1.0;
	_contamIdx = tmp;
}

void Individual::decr_contamIdx(float x){
	float tmp = _contamIdx - x;
	if(tmp < 0) tmp = 0.0;
	_contamIdx = tmp;
}


float Individual::time_clearance(){
	float t_clear;
	if(_symptomatic)
		t_clear = time_resolution() + _prdClearance;
	else
		t_clear = _timeAcquisition + _prdClearance;
	
	return t_clear;
}

void HCW::set_periodsWorked(vector<uint> x){
	int trace_level = 0;
	_periodsWorked = x;
	// DEBUG:
	if(trace_level>0){
		cout << " Periods worked assigned for HCW #" << get_uid();
		cout << " : ";
		displayVector(_periodsWorked);
	}
}


bool HCW::onDuty(float curr_time, float timestep){
	float dec = curr_time - (uint)curr_time; // Decimal part
	uint curr_period = uint(dec / timestep);
	bool onduty = (std::find(_periodsWorked.begin(), _periodsWorked.end(), curr_period) != _periodsWorked.end() ) ;
	return onduty;
}


float Patient::isolationDuration(){
	float res = 0;
	if (_timeIsolationStart >= 0 )
		res = _timeIsolationEnd - _timeIsolationStart;
	return res;
}


void Individual::calc_timeOnset(){
	_timeOnset = _timeAcquisition + _prdIncubation;
}



