//
// Created by Martin Rusilowicz on 30/08/2017.
//

#include <iostream>
#include "test_cases.h"
#include "binomial_test.h"


void test_cases()
{
    double g  = Binomial::test( EHypothesis::GREATER, 148, 16184, 0.0717471 );
    double ge = 1;
    double l  = Binomial::test( EHypothesis::LESS, 148, 16184, 0.0717471 );
    double le = 4.940656e-324;
    double t  = Binomial::test( EHypothesis::TWOTAILED, 148, 16184, 0.0717471 );
    double te = 9.881313e-324;

    std::cerr << "GREATER " << g << " AND " << ge << std::endl;
    std::cerr << "LESS " << l << " AND " << le << std::endl;
    std::cerr << "TWOTAILED " << t << " AND " << te << std::endl;
    
    exit(0);
}