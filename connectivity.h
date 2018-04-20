//
// Created by Martin Rusilowicz on 21/08/2017.
//

#ifndef COINFINDER_CONNECTIVITY_H
#define COINFINDER_CONNECTIVITY_H


#include <string>
#include "dataset.h"
#include "binomial_test.h"
#include "significance.h"

class Alpha;
class Beta;


class Connectivity
{
    public:
        static void run( const DataSet& dataset );
        
    private:
        static void _write_header(const TParameters& parameters);
        static bool _test_connectivity( const DataSet& dataset, const Alpha& alpha, const Beta& beta, int alpha_beta_gamma_count, double cor_sig_level );
};


#endif //COINFINDER_CONNECTIVITY_H
