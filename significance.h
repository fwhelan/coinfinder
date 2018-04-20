//
// Created by Martin Rusilowicz on 22/08/2017.
//

#ifndef COINFINDER_SIGNIFICANCE_H
#define COINFINDER_SIGNIFICANCE_H


#include "constants.h"


class Significance
{        
    public:
        static double correct( double sig_level, ECorrection correction, int num_tests ); 

};


#endif //COINFINDER_SIGNIFICANCE_H
