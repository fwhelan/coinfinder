/**
 * ADAPTED FROM: 
 *      https://github.com/chrchang/stats
 * 
 * LICENSE:
 *      GNU GENERAL PUBLIC LICENSE
 *      Version 3, 29 June 2007
 */

#ifndef COINFINDER_BINOMIAL_H
#define COINFINDER_BINOMIAL_H


#include "constants.h"


class Binomial
{
    public:
        static constexpr double INVALID_P = 1000;
        
    public:
        static double test( EHypothesis side, int successes, int observations, double rate); 
        static double two_sided(int successes, int observations, double rate);
        static double two_sided(int successes, int observations, double rate, int midp);
        static double one_sided_greater(int successes, int observations, double rate);
        static double one_sided_less(int successes, int observations, double rate);
};


#endif //COINFINDER_BINOMIAL_H
