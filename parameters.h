//
// Created by Martin Rusilowicz on 22/08/2017.
//

#ifndef COINFINDER_PARAMETERS_H
#define COINFINDER_PARAMETERS_H


#include <string>
#include "constants.h"


struct TParameters
{
    ECorrection correction;
    EHypothesis alt_hypothesis;
    EMethod method;
    double sig_level;
    ESetMode coin_set_mode; 
    EMaxMode coin_max_mode;
    bool verbose;
    bool permit_filter;
    bool Rmsgs;
    std::string alpha_name;
    std::string beta_name;
    std::string gamma_name;
    std::string alpha_file_name;
    std::string beta_file_name;
    std::string combined_file_name;
    std::string prefix;
    std::string phylogeny;
    int num_cores;
    double upper_filt_thres;
    double filt_thres;
    std::string deep_query_alpha;
    bool output_all;
    
    TParameters();
    void print_and_assert() const;
    static TParameters parse(int arg_count, const char** arg_vals);
};

#endif //COINFINDER_PARAMETERS_H
