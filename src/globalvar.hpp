//
//  globalvar.hpp
//  abm-cdiff
//
//  Created by David CHAMPREDON on 2017-06-19.
//  Copyright © 2017 David CHAMPREDON. All rights reserved.
//

#ifndef globalvar_hpp
#define globalvar_hpp

#include <stdio.h>
#include <random>

// ==== Random seed ====

// The random number generator
// must be declared as a global variable
// in order to get the random seed
// initialization right for Rcpp.
// (I don't know why this has to work like that...)

extern std::mt19937_64	RANDOM_GENERATOR;


// Isolation ward:
extern uint	ISO_WARD;

// Dummy number coding NA
extern uint HUGENUMBER;
extern uint		NA_UINT;
extern float	NA_FLOAT;

extern uint	MAX_N_RELAPSES;

#endif /* globalvar_hpp */
