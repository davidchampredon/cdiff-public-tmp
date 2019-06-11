//
//  utils.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-09.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;


void stopif(bool condition,
			string error_msg,
			int error_code = 1);

template <class T> void tabcout(string label, T value, uint left_length = 30){
	// Tabulated cout
	unsigned long x = label.length();
	string space_bck = "";
	
	while (x+space_bck.length() < left_length) {
		space_bck += " ";
	}
	cout << label << space_bck << ": " << value <<endl;
}

template <class T> void displayVector(vector<T> v){
	cout << endl<< "(size="<<v.size()<<")"<<endl;
	
	if(v.size()>0){
		cout << "[";
		for (int i=0; i<v.size()-1; i++)
		{
			cout << v[i] << "; ";
			if ((i+1)%10==0) cout<<endl;
		}
		cout << v[v.size()-1];
		cout<< "]" << endl;
	}
}

template <class T> void displayVector(vector< vector<T> > v){
	cout << endl<< "(size="<<v.size()<<")"<<endl<<"[";
	for (int i=0; i<v.size(); i++){
		displayVector(v[i]);
	}
}

template <class T> T sum(vector<T> v){
	T res = 0;
	for(size_t i=0; i<v.size(); i++) res += v[i];
	return res;
}


float beta_distribution(float a, float b);

#endif /* utils_hpp */
