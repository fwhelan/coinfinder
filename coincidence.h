//
// Created by Martin Rusilowicz on 18/08/2017.
//

#ifndef COINFINDER_CALCULATIONS_H
#define COINFINDER_CALCULATIONS_H

class DataSet;

#include <string>
#include <vector>
#include "dataset.h"
#include "elements.h"
#include "constants.h"

struct TParameters;


class Coincidence
{
    public:
        static void run(DataSet& dataset, const std::string& phylogeny, const std::string path);
        
    private:
        //static std::map<double, std::pair<std::string, std::string>> calc_phylogenetic_distances( const std::string& phylogeny, const std::string path);
	//static double pull_phylogenetic_maximum(  const std::map<double, std::pair<std::string, std::string>>& phylo_dists, const std::map<const Beta*, int>& edges_yain, const std::map<const Beta*, int>& edges_tain);
	//static std::vector<double> calculate_phylogenetic_distance(std::vector<std::string> edges_ovlp);
	//static std::vector<int> calculate_syntentic_distance(const std::map<const Beta*, int>& edgeList_B1, const std::map<const Beta*, int>& edgeList_B2, const id_lookup<Edge>& edge_table, const Alpha& alpha_yain, const Alpha& alpha_tain);
	//static std::pair<double, double> calc_secondaries(const std::map<double, std::pair<std::string, std::string>>& phylo_dists, const id_lookup<Edge>& edge_table, const Alpha& alpha_yain, const Alpha& alpha_tain, const std::vector<std::string>& edges_ovlp);
	static double calc_secondaries(const id_lookup<Edge>& edge_table, const Alpha& alpha_yain, const Alpha& alpha_tain, const std::vector<std::string>& edges_ovlp);
	//static std::string calc_common_ancestor(const DataSet& dataset, const std::string& phylogeny, const std::vector<std::string>& edges_union, std::string path);
	//static void _coincidence_to_p( const DataSet& dataset, const std::string& phylogeny, const Alpha& alpha_yain, const Alpha& alpha_tain, double cor_sig_level, int both_of, int one_of, int any_yain, int any_tain, const std::map<double, std::pair<std::string, std::string>>& phylo_dists, const id_lookup<Edge>& edge_table, const double max_act_phylodist, const std::vector<std::string>& edges_ovlp, const std::vector<std::string>& edges_union, std::string path);
	static double  _coincidence_to_p( const DataSet& dataset, const Alpha& alpha_yain, const Alpha& alpha_tain, double cor_sig_level, int both_of, int one_of, int any_yain, int any_tain, const id_lookup<Edge>& edge_table, const std::vector<std::string>& edges_ovlp, const std::vector<std::string>& edges_union, std::string path);
        static void _write_header(const DataSet& dataset);
};

#endif //COINFINDER_CALCULATIONS_H
