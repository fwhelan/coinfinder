//
// Created by Martin Rusilowicz on 21/08/2017.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include "connectivity.h"
#include "binomial_test.h"
#include "significance.h"
#include "elements.h"
#include "constants.h"


/**
 * Runs the connectivity tests, dumping p-values to STD.OUT. as appropriate.
 * 
 * @param dataset       Dataset object 
 * @param test          Test mode 
 * @param sig_level     Significance level (before correction) 
 * @param correction    Significance correction mode 
 */
void Connectivity::run( const DataSet& dataset )
{
    double cor_sig_level = Significance::correct( dataset.get_options().sig_level, dataset.get_options().correction, dataset.get_num_edges());

    Connectivity::_write_header( dataset.get_options());

    int count       = 0;
    int significant = 0;

    // For alpha
    for (const auto& it : dataset.get_alphas().get_table())
    {
        const Alpha& alpha = *it.second;

        for (const auto& kvp : alpha.get_edges())
        {
            if (Connectivity::_test_connectivity( dataset, alpha, *kvp.first, kvp.second, cor_sig_level ))
            {
                ++significant;
            }

            ++count;
        }
    }

    if (count != dataset.get_num_edges())
    {
        std::stringstream ss;
        ss << "An internal error occurred. Please submit a bug report. Details: Bonferroni correction calculated using a different number of edges (" << dataset.get_num_edges() << ") to what there actually were (" << count << ").";
        throw std::logic_error( ss.str());
    }

    std::cerr << "Iterated over " << count << " edges, " << significant << " of which were significant." << std::endl;

}


/**
 * Calculates the p-value for the connectivity (overlap) between `alpha` and `beta`.
 * 
 * @param dataset                   Dataset object 
 * @param alpha                     Alpha group
 * @param beta                      Beta group 
 * @param alpha_beta_gamma_count    Number of overlapping gammas between `alpha` and `beta`. This is for optimisation only as it is already available in the `AlphaId` class.
 */
bool Connectivity::_test_connectivity( const DataSet& dataset, const Alpha& alpha, const Beta& beta, const int alpha_beta_gamma_count, double cor_sig_level )
{
    int total_num_gamma = dataset.get_gammas().size();

    int    actual_overlap     = alpha_beta_gamma_count; // overlap
    double actual_rate        = static_cast<double>(actual_overlap) / total_num_gamma;
    int    alpha_observations = alpha.get_num_gammas(); // alpha col sums
    int    beta_observations  = beta.get_num_gammas(); // alpha col sums
    double beta_gamma_rate    = static_cast<double>(beta_observations) / total_num_gamma; // beta col sums
    double alpha_gamma_rate   = static_cast<double>(alpha_observations) / total_num_gamma;
    double alpha_overlap_rate = static_cast<double>(actual_overlap) / alpha_observations; // beta col sums
    double beta_overlap_rate  = static_cast<double>(actual_overlap) / beta_observations;
    double expected_rate      = alpha_gamma_rate * beta_gamma_rate;
    double expected_overlap   = expected_rate * total_num_gamma;
    double p_value;

    if (alpha_observations == 0 || beta_gamma_rate == 0)
    {
        p_value = Binomial::INVALID_P;
    }
    else
    {
        p_value = Binomial::test( dataset.get_options().alt_hypothesis, actual_overlap, total_num_gamma, expected_rate );
    }

    bool is_significant;

    switch (dataset.get_options().correction)
    {
        case ECorrection::FRACTION:
            is_significant = alpha_overlap_rate >= cor_sig_level;
            break;

        default:
            is_significant = p_value <= cor_sig_level;
            break;
    }

    bool filter = dataset.get_options().output_all || is_significant;

    if (!dataset.get_options().deep_query_alpha_file_name.empty())
    {
	filter = false;
    	std::cerr << "Reading deep query alpha file..." << std::endl;
	std::ifstream file_in;
	file_in.open( dataset.get_options().deep_query_alpha_file_name );
	std::string cell;
        if (!file_in)
        {
	        std::stringstream ss;
	        ss << "Failed to open file: " << dataset.get_options().deep_query_alpha_file_name;
	        throw std::logic_error(ss.str());
	}
	int n=0;
	while (getline( file_in, cell))
	{
		if (alpha.get_name() == cell) {
			filter = true;
		}
	}
    }

    if (dataset.get_options().verbose)
    {
        std::cerr << "ALPHA '" << alpha.get_name() << "' IS " << actual_overlap << " OF " << alpha_observations << " " << beta.get_name() << " (" << int( 100 * alpha_overlap_rate ) << "%)" << std::endl;
    }

    if (filter)
    {
        std::cout << alpha.get_name() // alpha::name
                  << "\t" << beta.get_name() // beta::name
                  << "\t" << actual_overlap // overlap::actual
                  << "\t" << actual_rate // overlap::actual::rate
                  << "\t" << alpha_observations // alpha::gamma
                  << "\t" << alpha_gamma_rate // alpha::gamma::rate
                  << "\t" << beta_observations // beta::gamma
                  << "\t" << beta_gamma_rate // beta::gamma::rate
                  << "\t" << expected_overlap // overlap::expected
                  << "\t" << expected_rate // overlap::expected::rate
                  << "\t" << total_num_gamma // gamma
                  << "\t" << alpha_overlap_rate // alpha::overlap::rate
                  << "\t" << beta_overlap_rate // beta::overlap::rate
                  << "\t" << p_value // overlap::p
                  << "\t" << cor_sig_level // overlap::p
                  << "\t" << ( is_significant ? 1 : 0 ) // overlap::is_significant
                  << std::endl;

        return true;
    }

    return false;
}


/**
 * Writes the output header.
 */
void Connectivity::_write_header( const TParameters& parameters )
{
    const std::string& alpha = parameters.alpha_name;
    const std::string& beta  = parameters.beta_name;
    const std::string& gamma = parameters.gamma_name;

    std::cout << "$ " << alpha
              << "\t" << "$ " << beta
              << "\t" << "# Actual overlap"
              << "\t" << "% Actual overlap rate"
              << "\t" << "# Size of " << alpha
              << "\t" << "% Rate of " << alpha
              << "\t" << "# Size of " << beta
              << "\t" << "% Rate of " << beta
              << "\t" << "£ Expected overlap"
              << "\t" << "% Expected overlap rate"
              << "\t" << "# Number of " << gamma
              << "\t" << "% " << alpha << " overlap rate"
              << "\t" << "% " << beta << " overlap rate"
              << "\t" << "£ p-value"
              << "\t" << "£ Significance level"
              << "\t" << "! Is significant"
              << std::endl;
}
