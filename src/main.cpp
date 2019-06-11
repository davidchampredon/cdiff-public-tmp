//
//  main.cpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#include <stdio.h>
#include <iostream>


#include "Individual.hpp"
#include "HealthCareSetting.hpp"
#include "run.hpp"

using namespace std;




int main(int argc, const char * argv[]) {
   
    //run_test_multi();
    run_test_single();
    cout << endl << " --- END --- " << endl;
	return 0;
}

void dummy_test(){
    float m = 13.0;
    float v = 2.0;
    float b = v / m;
    float a = m / b;
    gamma_distribution<float> dgam(a, b);
    int n = 1e4;
    float s= 0.0;
    float s2 = 0.0;
    for(int i=0; i<n; i++){
        float tmp =dgam(RANDOM_GENERATOR) ;
        s +=  tmp;
        s2 +=tmp*tmp;
    }
    cout <<"mean= "<< s/n<<endl;
    cout <<"var = "<< s2/n - (s/n)*(s/n) << endl;
}
