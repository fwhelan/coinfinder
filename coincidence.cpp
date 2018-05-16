//
// Created by Martin Rusilowicz on 18/08/2017.
//

#include <iostream>
#include <fstream>
#include <vector>

#include "coincidence.h"
#include "dataset.h"
#include "elements.h"
#include "parameters.h"
#include "significance.h"
#include "binomial_test.h"

#include <Python.h>

/**
 * Runs coincidence analysis
 */
void Coincidence::run( const DataSet& dataset, /**< Dataset */
		       const std::string& phylogeny )
{
    Coincidence::_write_header();

    std::cerr << "Iterating matrix..." << std::endl;
    const id_lookup<Alpha>& alpha_table = dataset.get_alphas();
    const id_lookup<Edge>& edge_table = dataset.get_edges();

    const double cor_sig = Significance::correct( dataset.get_options().sig_level, dataset.get_options().correction, dataset.get_num_edges());

    //
    // Create matrix
    //
    int num_alphas = static_cast<int>(alpha_table.get_table().size());
    int count      = 0;

    //
    // Determine the maximum phylogenetic distance in given tree and save distance information into hash
    //
    std::cerr << "Calculating phylogenetic distance information..." << std::endl;
    std::map<std::pair<std::string, std::string>, double> phylo_dists = Coincidence::calc_phylogenetic_distances( phylogeny );
    std::cerr << "Calculating the maximum phylogenetic distance..." << std::endl;
    double max_phylo_dist = 0;
    for(const auto &p : phylo_dists) {
	if((p.second) > max_phylo_dist) { max_phylo_dist = (p.second); }
    }
    //
    // Iterate first alpha
    //
    std::cerr << "Running analyses..." << std::endl;
    for (const auto& kvp_yain : alpha_table.get_table())
    {
        const Alpha& alpha_yain = *kvp_yain.second;
        count++;
	if (( count % 1000 ) == 0)
        {
            std::cerr << "- Row " << count << " of " << num_alphas << std::endl;
        }
        
	const std::map<const Beta*, int>& edges_yain = alpha_yain.get_edges();
        int num_edges_yain = static_cast<int>(edges_yain.size());
	
	//
        // Iterate second alpha
        //
        for (const auto& kvp_tain : alpha_table.get_table())
        {
            const Alpha& alpha_tain = *kvp_tain.second;
            
	    if (alpha_tain.get_name().compare( alpha_yain.get_name()) <= 0)
            {
                continue;
            }

            const std::map<const Beta*, int>& edges_tain = alpha_tain.get_edges();
            int num_edges_tain = static_cast<int>(edges_tain.size());

	    //
            // Count overlaps
	    //
            int overlaps = 0;
	    std::vector<std::string> edges_ovlp;

            for (const auto& kvp_edges : edges_yain)
            {
                const Beta& beta = *kvp_edges.first;
                auto it = edges_tain.find( &beta );
                if (it != edges_tain.end())
                {
                    overlaps += 1;
		    edges_ovlp.push_back(beta.get_name());
                }
            }

	    // Count total range
            int total_range = num_edges_yain + num_edges_tain - overlaps;
            Coincidence::_coincidence_to_p( dataset, alpha_yain, alpha_tain, cor_sig, overlaps, total_range, num_edges_yain, num_edges_tain, phylo_dists, edge_table, max_phylo_dist, edges_ovlp );
        }
    }
}

/**
 * Calculate the phylogenetic distances between all nodes in the given phylogeny up front
*/
std::map<std::pair<std::string,std::string>, double> Coincidence::calc_phylogenetic_distances( const std::string& phylogeny )
{
	/*Definitions*/
	std::map<std::pair<std::string,std::string>, double> phylogenetic_distances;

	/*Embedded Python*/
        Py_Initialize();
	PyObject* pValue;

        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append(\".\")");

        PyObject* pName = PyUnicode_DecodeFSDefault("phylomax");
        PyObject* pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if (pModule != NULL) {
                PyObject* pFunc = PyObject_GetAttrString(pModule, "calc");
                if (pFunc && PyCallable_Check(pFunc)) {
                        PyObject* pArgs = PyTuple_New(1);
                        //PyTuple_SetItem(pArgs, 0, PyUnicode_FromString("core.nex.con.tre.newick")); //send tree as first element
			PyTuple_SetItem(pArgs, 0, PyUnicode_FromString(phylogeny.c_str()));
			pValue = PyObject_CallObject(pFunc, pArgs);
                        Py_DECREF(pArgs);
                        if (pValue == NULL) {
                                std::cerr << "pValue is null" << std::endl;
                                PyErr_Print();
                        } else {
                                if (PyUnicode_Check(pValue) == 1) { //return value is a string; there was an error
                                        std::cerr << "There was an error in Python code" << std::endl;
                                        PyErr_Print();
                                } else if (PyFloat_Check(pValue) == 1) {
					std::cerr << "There was an error in Python code" << std::endl;
					PyErr_Print();
				} else if (PyList_Check(pValue) == 1) { //return value is a list
                                        /*Save returned list into a map with beta names as paired keys and dist as value*/
					for(int a=0; a<PyList_Size(pValue); a=a+3) {
                                                PyObject *value1 = PyList_GetItem(pValue, a);
						PyObject *value2 = PyList_GetItem(pValue, a+1);
						PyObject *value3 = PyList_GetItem(pValue, a+2);
						phylogenetic_distances[std::make_pair(PyUnicode_AsUTF8(value1),PyUnicode_AsUTF8(value2))] = PyFloat_AsDouble(value3);
                                        }

                                } else {
                                        std::cerr << "Error: a list isn't being returned" << std::endl; //catch this properly
                                }
				Py_DECREF(pValue);
                        }
                }
        } else {
		PyErr_Print();
	}
        Py_Finalize();
	return(phylogenetic_distances);
}

/**
 * Calculate the maximum observed phylogenetic distances and the average synthetic distance between all overlapping edge pairs.
*/
std::pair<double, double> Coincidence::calc_secondaries(
						const std::map<std::pair<std::string, std::string>, double>& phylo_dists,
						const id_lookup<Edge>& edge_table,
						const Alpha& alpha_yain,
						const Alpha& alpha_tain,
						const std::vector<std::string>& edges_ovlp )
{
	double syn_sums = 0;
	double phy_max = 0;
	//Cycle through edges_yain, find values also in edges_tain
	for(const auto &p : edges_ovlp) {
		//Synthetic distance
		Edge edge1 = edge_table.find_id(alpha_yain.get_name()+"-"+(p.c_str()));
		Edge edge2 = edge_table.find_id(alpha_tain.get_name()+"-"+(p.c_str()));
		syn_sums += abs(edge1.get_weight()-edge2.get_weight());
                for(const auto &q : edges_ovlp) {
			//Phylogenetic distance
                        if (phylo_dists.find(std::make_pair(p.c_str(), q.c_str()))->second > phy_max) {
                                phy_max = phylo_dists.find(std::make_pair(p.c_str(), q.c_str()))->second;
                        }
                }
        }
	double syn_avg = syn_sums/edges_ovlp.size();
	return(std::make_pair(phy_max, syn_avg));
}


/** Calculates a p-value of coincicidence.                                                               */
void Coincidence::_coincidence_to_p( const DataSet& dataset,        /**< Dataset                         */
                                     const Alpha& alpha_yain,       /**< First alpha                     */
                                     const Alpha& alpha_tain,       /**< Second alpha                    */
                                     double cor_sig_level,    /**< Corrected significance level    */
                                     int both_of,             /**< "both of" count (yain AND tain) */
                                     int one_of,              /**< "one of" count (yain OR tain)   */
                                     int any_yain,            /**< yain count (yain)               */
                                     int any_tain,             /**< tain count (tain)               */
				     const std::map<std::pair<std::string, std::string>, double>& phylo_dists,
				     const id_lookup<Edge>& edge_table,
				     const double max_act_phylodist,
				     const std::vector<std::string>& edges_ovlp )
{
    int       num_observations;
    const int max_coincidence = dataset.get_betas().size();

    const TParameters& options = dataset.get_options();

    switch (options.coin_set_mode)
    {
        case ESetMode::INTERSECTION:
            num_observations = one_of;

            if (num_observations == 0)
            {
                if (options.verbose)
                {
                    std::cerr << "Rejected (" << alpha_yain.get_name() << ", " << alpha_tain.get_name() << ") because there are no observations." << std::endl;
                }
                return;
            }
            break;

        case ESetMode::FULL:
            num_observations = max_coincidence;
            break;

        default:
            throw std::logic_error( "Invalid options around SET_MODE_MASK." );
    }

    // From this we can work out the chance of alpha/beta occurring at the same time
    double chance_i = static_cast<double>(any_yain) / static_cast<double>(num_observations);
    double chance_j = static_cast<double>(any_tain) / static_cast<double>(num_observations);

    double not_cross_1_chance = static_cast<double>( num_observations - any_yain ) / static_cast<double>(num_observations);
    double not_cross_2_chance = static_cast<double>( num_observations - any_tain ) / static_cast<double>(num_observations);

    double rate;
    int    successes;

    switch (options.coin_max_mode)
    {
        case EMaxMode::AVOID:
        {
            successes = one_of - both_of;   // note that the upper triangle (i,j) is 1&2 whilst the lower triangle (j,i) is 1|2
            rate      = ( chance_i * not_cross_2_chance ) + ( chance_j * not_cross_1_chance );
            break;
        }
        case EMaxMode::ACCOMPANY:
        {
            successes = both_of;
            rate      = chance_i * chance_j;
            break;
        }
        default:
        {
            throw std::logic_error( "Invalid options around MAX_MODE_MASK." );
        }
    }

    // This causes problems, get rid of it
    if (rate == 0 || rate == 1.0)
    {
        if (options.verbose)
        {
            std::cerr << "Rejected (" << alpha_yain.get_name() << ", " << alpha_tain.get_name() << ") because the rate is " << rate << "." << std::endl;
        }
        return;
    }

    // Binomial test p-value
    double p_value = Binomial::test( options.alt_hypothesis, successes, num_observations, rate );

    if (options.verbose)
    {
        std::cerr << "*******************************" << std::endl;
        std::cerr << "* yain                " << alpha_yain.get_name() << "." << std::endl;
        std::cerr << "* tain                " << alpha_tain.get_name() << "." << std::endl;
        std::cerr << "*------------------------------" << std::endl;
        std::cerr << "* any_yain            " << any_yain << "." << std::endl;
        std::cerr << "* any_tain            " << any_tain << "." << std::endl;
        std::cerr << "* both_of             " << both_of << "." << std::endl;
        std::cerr << "* one_of              " << one_of << "." << std::endl;
        std::cerr << "* max_coincidence     " << max_coincidence << "." << std::endl;
        std::cerr << "*------------------------------" << std::endl;

        std::cerr << "*------------------------------" << std::endl;
        std::cerr << "* chance_i            " << chance_i << "." << std::endl;
        std::cerr << "* chance_j            " << chance_j << "." << std::endl;
        std::cerr << "* not_cross_1_chance  " << not_cross_1_chance << "." << std::endl;
        std::cerr << "* not_cross_2_chance  " << not_cross_2_chance << "." << std::endl;
        std::cerr << "*------------------------------" << std::endl;
        std::cerr << "* rate                " << rate << "." << std::endl;
        std::cerr << "* successes           " << successes << "." << std::endl;
        std::cerr << "* num_observations    " << num_observations << "." << std::endl;
        std::cerr << "*------------------------------" << std::endl;
        std::cerr << "* p_value LESS        " << Binomial::one_sided_less( successes, num_observations, rate ) << "." << std::endl;
        std::cerr << "* p_value GREATER     " << Binomial::one_sided_greater( successes, num_observations, rate ) << "." << std::endl;
        std::cerr << "* p_value TWOTAILED   " << Binomial::two_sided( static_cast<uint32_t>(successes), static_cast<uint32_t>(num_observations), rate ) << "." << std::endl;
        std::cerr << "*******************************" << std::endl;
    }

    if (p_value > cor_sig_level)
    {
        if (options.verbose)
        {
            std::cerr << "Rejected (" << alpha_yain.get_name() << ", " << alpha_tain.get_name() << ") because it isn't significant with p = " << p_value << "." << std::endl;
        }

        if (!options.output_all)
        {
            return;
        }
    }
    else
    {
        if (options.verbose)
        {
            std::cerr << "Accepted (" << alpha_yain.get_name() << ", " << alpha_tain.get_name() << ") because it is significant with p = " << p_value << "." << std::endl;
        }
    }

    //Result is significant: calculate maximum observed phylogenetic distance and average synthetic distance to output to file
    std::pair<double, double> secondaries = calc_secondaries(phylo_dists, edge_table, alpha_yain, alpha_tain, edges_ovlp);
    double max_obs_phylodist = secondaries.first;
    //double max_obs_phylodist = 0;
    double phylo_output = max_act_phylodist - max_obs_phylodist;
    double avg_syndist = secondaries.second;
    //double avg_syndist = 0;
    
    std::cout << alpha_yain.get_name()
              << "\t" << alpha_tain.get_name()
              << "\t" << p_value
	      << "\t" << phylo_output
	      << "\t" << avg_syndist
              << "\t" << successes
              << "\t" << num_observations
              << "\t" << rate
              << "\t" << static_cast<int>(rate * num_observations + 0.5)
              << "\t" << any_yain
              << "\t" << any_tain
              << "\t" << chance_i
              << "\t" << chance_j
              << std::endl;
}

/**
 * Writes the header
 */
void Coincidence::_write_header()
{
    std::cout << "Source"
              << "\t" << "Target"
              << "\t" << "p"
	      << "\t" << "Max phylogenetic distance"
	      << "\t" << "Avg synthetic distance"
              << "\t" << "successes"
              << "\t" << "observations"
              << "\t" << "rate"
              << "\t" << "expected"
              << "\t" << "total source"
              << "\t" << "total target"
              << "\t" << "fraction source"
              << "\t" << "fraction target"
              << std::endl;
}


