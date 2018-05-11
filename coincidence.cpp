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
void Coincidence::run( const DataSet& dataset /**< Dataset */ )
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
    // Determine the maximum phylogenetic distance in given tree
    // and save distance information into hash
    std::cerr << "Calculating the maximum phylogenetic distance..." << std::endl;
    double max_phenoDist = Coincidence::maximum_phylogenetic_distance( );
    std::cerr << "Maximum phylogenetic distance is " << max_phenoDist << std::endl;

    //
    // Iterate first alpha
    //
    for (const auto& kvp_yain : alpha_table.get_table())
    {
        const Alpha& alpha_yain = *kvp_yain.second;
	std::cerr << "alpha is: " << alpha_yain.get_name() << std::endl; /*Print alpha value for sanity check*/
        count++;
        if (( count % 1000 ) == 0)
        {
            std::cout << "- Row " << count << " of " << num_alphas << std::endl;
        }

        const std::map<const Beta*, int>& edges_yain = alpha_yain.get_edges();
        int num_edges_yain = static_cast<int>(edges_yain.size());
	
	/*Print beta values for sanity check*/
	/*for (std::pair<const Beta*, int> p : edges_yain) {
	    std::cerr << "Betas: " << (p.first)->get_name() << ' ' << p.second << '\n';*/
	    /*Get weight for the associated edge*/
	    /*Edge edge = edge_table.find_id(alpha_yain.get_name()+"-"+(p.first)->get_name());
	    std::cerr << "With edge weights: " << edge.get_weight() << std::endl;
	}*/
	std::cerr << "I have " << num_edges_yain << " edges" << std::endl;
	    
	//
        // Iterate second alpha
        //
        for (const auto& kvp_tain : alpha_table.get_table())
        {
            const Alpha& alpha_tain = *kvp_tain.second;
	    std::cerr << "alpha is: " << alpha_tain.get_name() << std::endl; /*Print alpha value for sanity check*/
            if (alpha_tain.get_name().compare( alpha_yain.get_name()) <= 0)
            {
                continue;
            }

            const std::map<const Beta*, int>& edges_tain = alpha_tain.get_edges();
            int num_edges_tain = static_cast<int>(edges_tain.size());

	    /*Print beta values for sanity check*/
	    /*for (std::pair<const Beta*, int> p : edges_tain) {
            	std::cerr << "Betas: " << (p.first)->get_name() << ' ' << p.second << '\n';*/
		/*Get weight for the associated edge*/
            	/*Edge edge = edge_table.find_id(alpha_tain.get_name()+"-"+(p.first)->get_name());
            	std::cerr << "With edge weights: " << edge.get_weight() << std::endl;
            }*/

	    std::cerr << "I have " << num_edges_tain << " edges" << std::endl;

	    //
            // Count overlaps
	    //
	    std::vector<std::string> edges_ovlp;
            int overlaps = 0;

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

	    std::cerr << "I have " << overlaps << " number of overlapping edges" << std::endl;

	    // Count total range
            int total_range = num_edges_yain + num_edges_tain - overlaps;

	    // Calculate phylogenetic distance via Python Interpreter
	    //std::vector<double> phenoList = Coincidence::calculate_phylogenetic_distance( edges_yain, edges_tain );
	    std::vector<double> phenoList = Coincidence::calculate_phylogenetic_distance( edges_ovlp );
	    /*Temp print out vector to decide on next step*/
	    for (auto v : phenoList) { std::cerr << v << " "; }
	    std::cerr << std::endl;
	    /*Calculate the average of the elements in phenoList*/
	    double avg_pheno = 0;
	    double max_pheno = 0;
	    for (auto v : phenoList) {
		avg_pheno += v;
		if (v > max_pheno) { max_pheno = v; }
	    }
	    avg_pheno = avg_pheno/phenoList.size();
	    std::cerr << "The maximum phenotypic distance is " << max_phenoDist << std::endl;
	    std::cerr << "The maximum observed phenotypic distance is " << max_pheno << std::endl;
	    /*Print length of phenoList*/
	    std::cerr << "phenoList is " << phenoList.size() << " long" << std::endl;
	    // Calculate syntentic distance
	    std::vector<int> synList = Coincidence::calculate_syntentic_distance( edges_yain, edges_tain, edge_table, alpha_yain, alpha_tain );
	    /*Temp print out vector to decide on next step*/
	    for (auto v : synList) { std::cerr << v << " "; }
	    std::cerr << std::endl;
	    /*Calculate the average of the elements in synList*/
	    double avg_syn = 0;
	    for (auto v : synList) {
		avg_syn += v;
	    }
	    avg_syn = avg_syn/synList.size();
	    std::cerr << "The average synthetic distance is " << avg_syn << std::endl;
	    /*Print length of synList*/
	    std::cerr << "synList is " << synList.size() << " long" << std::endl;
	    return;

            Coincidence::_coincidence_to_p( dataset, alpha_yain, alpha_tain, cor_sig, overlaps, total_range, num_edges_yain, num_edges_tain );
        }
    }
}

/**
 * Calculate the maximum phylogenetic distance in the given phylogeny
*/
double Coincidence::maximum_phylogenetic_distance( )
{
        Py_Initialize();
        PyRun_SimpleString("import sys");
        PyRun_SimpleString("sys.path.append(\".\")");

        PyObject* pName = PyUnicode_DecodeFSDefault("phylomax");
        PyObject* pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if (pModule != NULL) {
                PyObject* pFunc = PyObject_GetAttrString(pModule, "calc");
                if (pFunc && PyCallable_Check(pFunc)) {
                        /*Convert vector to array for send*/
                        //PyObject* pArgs = PyTuple_New(2+edgesB1.size()+edgesB2.size());
                        PyObject* pArgs = PyTuple_New(1);
                        PyObject* pValue;
                        PyTuple_SetItem(pArgs, 0, PyUnicode_FromString("core.nex.con.tre.newick")); //send tree as first element
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
					return(PyFloat_AsDouble(pValue));
				} else if (PyList_Check(pValue) == 1) { //return value is a list
                                        Py_ssize_t sizey = PyList_Size(pValue);
                                        for(int a=0; a<PyList_Size(pValue); a=a+1) {
                                                PyObject *value = PyList_GetItem(pValue, a);
                                                //dataList.push_back(PyFloat_AsDouble(value));
                                                //std::cerr << "value: " << PyFloat_AsDouble(value) << std::endl;
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
	return(-1);
}

/**
 * Calculate the phylogenetic distances between all combinations of nodes in edges_yain and edges_tain
*/
//std::vector<double> Coincidence::calculate_phylogenetic_distance( 	const std::map<const Beta*, int>& edgeList_B1,
//							const std::map<const Beta*, int>& edgeList_B2)
std::vector<double> Coincidence::calculate_phylogenetic_distance(	std::vector<std::string> edges_ovlp)
{
	/*Definitions*/
	std::vector<double> dataList;
	/*Save map values into vector for push to Python*/
	/*std::vector<std::string> edgesB1;
	for(std::pair<const Beta*, int> p : edgeList_B1) {
		edgesB1.push_back((p.first)->get_name());
	}
	std::vector<std::string> edgesB2;
	for(std::pair<const Beta*, int> q : edgeList_B2) {
		edgesB2.push_back((q.first)->get_name());
	}*/
	
	Py_Initialize();
	PyRun_SimpleString("import sys");
	PyRun_SimpleString("sys.path.append(\".\")");

	PyObject* pName = PyUnicode_DecodeFSDefault("phylodist");
	PyObject* pModule = PyImport_Import(pName);
	Py_DECREF(pName);
	if (pModule != NULL) {
		PyObject* pFunc = PyObject_GetAttrString(pModule, "calc");
		if (pFunc && PyCallable_Check(pFunc)) {
			/*Convert vector to array for send*/
			//PyObject* pArgs = PyTuple_New(2+edgesB1.size()+edgesB2.size());
			PyObject* pArgs = PyTuple_New(1+edges_ovlp.size());
			PyObject* pValue;
			PyTuple_SetItem(pArgs, 0, PyUnicode_FromString("core.nex.con.tre.newick")); //send tree as first element
			//PyTuple_SetItem(pArgs, 1, PyLong_FromLong(edgesB1.size())); //set a break counter to know where to split the input into the 2 edge lists
			/*for (size_t i=0; i < edgesB1.size(); i++) {
				pValue = PyUnicode_FromString(edgesB1[i].c_str());
				if (!pValue) {
					//catch
				}
				PyTuple_SetItem(pArgs, (i+2), pValue);
			}
			for (size_t i=0; i < edgesB2.size(); i++) {
				pValue = PyUnicode_FromString(edgesB2[i].c_str());
				if (!pValue) {
					//catch
				}
				PyTuple_SetItem(pArgs, (i+2+edgesB1.size()), pValue);
			}*/
			for (size_t i=0; i < edges_ovlp.size(); i++) {
				pValue = PyUnicode_FromString(edges_ovlp[i].c_str());
				if (!pValue) {
					//catch
				}
				PyTuple_SetItem(pArgs, (i+1), pValue);
			}
			pValue = PyObject_CallObject(pFunc, pArgs);
			Py_DECREF(pArgs);

			if (pValue == NULL) {
				std::cerr << "pValue is null" << std::endl;
				PyErr_Print();
			} else {
				if (PyUnicode_Check(pValue) == 1) { //return value is a string; there was an error
					std::cerr << "Node " << PyUnicode_AsUTF8(pValue) << " was not found in phylogeny" << std::endl;
					PyErr_Print();
				}
				if (PyList_Check(pValue) == 1) { //return value is a list
					Py_ssize_t sizey = PyList_Size(pValue);
					for(int a=0; a<PyList_Size(pValue); a=a+1) {
						PyObject *value = PyList_GetItem(pValue, a);
						dataList.push_back(PyFloat_AsDouble(value));
						//std::cerr << "value: " << PyFloat_AsDouble(value) << std::endl;
					}

				} else {
					//std::cerr << "Error: a list isn't being returned" << std::endl; //catch this properly
				}
				Py_DECREF(pValue);
			}
		}
	}
	Py_Finalize();

	/*Cycle through the first few elements of distList to see what kind of state its in*/
	//for(int i=0; i<dataList.size(); ++i) {
	//	std::cerr << "dataList: " << dataList[i] << std::endl;
	//}
	return(dataList);
}

std::vector<int> Coincidence::calculate_syntentic_distance( const std::map<const Beta*, int>& edgeList_B1,
						const std::map<const Beta*, int>& edgeList_B2,
						const id_lookup<Edge>& edge_table,
						const Alpha& alpha_yain,
						const Alpha& alpha_tain)
{
	std::vector<int> synList;
	/* For all intersecting betas between the 2 maps (i.e. genomes which contain both alphas), calculate the distance
	 * between the 2 alphas by retrieving the edge weights and looking at the absolute difference betwene them*/
	for (std::pair<const Beta*, int> p : edgeList_B1) {
		for (std::pair<const Beta*, int> q : edgeList_B2) { //TODO- this is so costly!
			if ((p.first)->get_name() == (q.first)->get_name()) {
            			/*std::cerr << "Betas: " << (p.first)->get_name() << " and " << (q.first)->get_name() << " match"  << '\n';
				std::cerr << "Weights are: " << edge1.get_weight() << " and " << edge2.get_weight() << std::endl;
				std::cerr << "Diff is: " << abs(edge1.get_weight()-edge2.get_weight()) << std::endl;*/
				Edge edge1 = edge_table.find_id(alpha_yain.get_name()+"-"+(p.first)->get_name());
				Edge edge2 = edge_table.find_id(alpha_tain.get_name()+"-"+(q.first)->get_name());
				synList.push_back(abs(edge1.get_weight()-edge2.get_weight()));
			}
		}
	}
	return(synList);
}


/** Calculates a p-value of coincicidence.                                                               */
void Coincidence::_coincidence_to_p( const DataSet& dataset,        /**< Dataset                         */
                                     const Alpha& alpha_yain,       /**< First alpha                     */
                                     const Alpha& alpha_tain,       /**< Second alpha                    */
                                     const double cor_sig_level,    /**< Corrected significance level    */
                                     const int both_of,             /**< "both of" count (yain AND tain) */
                                     const int one_of,              /**< "one of" count (yain OR tain)   */
                                     const int any_yain,            /**< yain count (yain)               */
                                     const int any_tain             /**< tain count (tain)               */ )
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

    std::cout << alpha_yain.get_name()
              << "\t" << alpha_tain.get_name()
              << "\t" << p_value
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
    std::cout << "$ Source"
              << "\t" << "$ Target"
              << "\t" << "£ p"
              << "\t" << "# successes"
              << "\t" << "# observations"
              << "\t" << "% rate"
              << "\t" << "£ expected"
              << "\t" << "# total source"
              << "\t" << "# total target"
              << "\t" << "% fraction source"
              << "\t" << "% fraction target"
              << std::endl;
}


