//
// Created by Martin Rusilowicz on 22/08/2017.
//

#include <iostream>
#include <cmath>
#include "significance.h"


double Significance::correct( double sig_level, ECorrection correction, int num_tests )
{
    switch(correction)
    {
        case ECorrection::NONE:
        {
            std::cerr << "No significance correction, the significance level is " << sig_level << "." << std::endl;
            return sig_level;
        }
            
        case ECorrection::BONFERRONI:
        {
            double corrected = sig_level / num_tests;
            std::cerr << "Bonferroni significance correction, given " << num_tests << " tests, the significance level has been reduced from " << sig_level << " to " << corrected << "." << std::endl;
            
            if (std::isinf(corrected))
            {
                throw std::logic_error("Infinite significance level. This is probably due to an earlier error.");
            }
            
            return corrected;
        }
        
        case ECorrection::FRACTION:
        {
            std::cerr << "No significance level, using fractional cutoff of " << sig_level << "." << std::endl;
            return sig_level;
        }
        
        default:
        {
            throw std::logic_error("Invalid switch on Significance::correct(correction).");
        }
    }
}
