//
//  Disease.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-07-07.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include "Disease.hpp"


Disease::Disease(float prd_incubation_mean,
				 float prd_incubation_var,
				 float prd_symptom_mean,
				 float prd_symptom_var,
				 float proba_clear_prm,
				 float prd_clear_asymptom_mean,
				 float prd_clear_asymptom_var,
				 float prd_clear_symptom_mean,
				 float prd_clear_symptom_var,
				 float contamIdx_asymptomatic,
				 float contamIdx_symptomatic,
				 float hazard_relapse,
				 float shape_relapse,
				 float infectIdx_init){
	
	_prd_incubation_mean = prd_incubation_mean;
	_prd_incubation_var = prd_incubation_var;
	_prd_symptom_mean = prd_symptom_mean;
	_prd_symptom_var = prd_symptom_var;
	
	_proba_clear_prm = proba_clear_prm;
	_prd_clear_asymptom_mean = prd_clear_asymptom_mean;
	_prd_clear_asymptom_var = prd_clear_asymptom_var;
	_prd_clear_symptom_mean = prd_clear_symptom_mean;
	_prd_clear_symptom_var = prd_clear_symptom_var;

	_contamIdx_max_asymptom = contamIdx_asymptomatic;
	_contamIdx_ratio_symptomatic = contamIdx_symptomatic;
	
	_hazard_relapse = hazard_relapse;
	_shape_relapse = shape_relapse;
	
	_infectIdx_init = infectIdx_init;
}

void Disease::show(){
	
	cout << " ==== DISEASE show() ====" <<endl;
	
	cout <<"name: "<<	_name <<endl;
	
	cout <<"_prd_incubation_mean: "<<	_prd_incubation_mean <<endl;
	cout <<"_prd_incubation_var: "<<	_prd_incubation_var <<endl;
	cout <<"_prd_symptom_mean: "<<	_prd_symptom_mean <<endl;
	cout <<"_prd_symptom_var: "<<	_prd_symptom_var <<endl;
	cout <<"_proba_clear_prm: "<<	_proba_clear_prm <<endl;
	cout <<"_prd_clear_asymptom_mean: "<<	_prd_clear_asymptom_mean <<endl;
	cout <<"_prd_clear_asymptom_var: "<<	_prd_clear_asymptom_var <<endl;
	cout <<"_prd_clear_symptom_mean: "<<	_prd_clear_symptom_mean <<endl;
	cout <<"_prd_clear_symptom_var: "<<	_prd_clear_symptom_var <<endl;
	cout <<"_contamIdx_max_asymptom: "<<	_contamIdx_max_asymptom <<endl;
	cout <<"_contamIdx_ratio_symptomatic: "<<	_contamIdx_ratio_symptomatic <<endl;
	cout <<"_hazard_relapse: " <<_hazard_relapse <<endl;
	cout <<"_shape_relapse: "  <<_shape_relapse <<endl;
	cout <<"_infectIdx_init: " <<_infectIdx_init <<endl;
	
	cout << " ==== DISEASE END show() ====" <<endl;
}
