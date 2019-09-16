//
// Created by Martin Rusilowicz on 17/08/2017.
//

#ifndef COINFINDER_FAMILYSET_H
#define COINFINDER_FAMILYSET_H

#include <string>
#include <algorithm>
#include "id_lookup.h"
#include "parameters.h"


class Alpha;
class Beta;
class Edge;
class Gamma;


class DataSet
{
    private:
        const TParameters& _options;
        id_lookup<Alpha> _alphas; // alpha or beta
        id_lookup<Beta> _betas;
        id_lookup<Gamma> _gammas;
        int _num_edges;
	int _num_coincident_edges;
	id_lookup<Edge> _edges;
        
    public:
        //DataSet();
	explicit DataSet(const TParameters& options);
	//DataSet(const TParameters& options);
        ~DataSet();
        
        const id_lookup<Alpha>& get_alphas() const;
        const id_lookup<Beta>& get_betas() const;
        const id_lookup<Gamma>& get_gammas() const;
	const id_lookup<Edge>& get_edges() const;
	const int get_num_betas() const;
	const int get_num_alphas() const;

	void read_files( const std::string& alpha_file_name, const std::string& beta_file_name, const std::string& combined_file_name, const std::string& phylogeny_file_name, const double filt_thres, const double upper_filt_thres );
        
        int get_num_edges() const;
	int get_num_coincident_edges() const;
        const TParameters& get_options() const;

	void _generate_coincident_edge( Alpha& alpha1, Alpha& alpha2, double p_value );

	//static bool isForbidden( char c );
	//void _read_alpha_file( const std::string& file_name );

        
    private:
	//std::string RemoveIllegalChars(std::string cell);
	//bool isForbidden( char c );
        void _read_alpha_file( const std::string& file_name );
        void _read_beta_file( const std::string& file_name );
        void _read_combined_file( const std::string& file_name );

        void _drop_empty();
	void _phylo_check( const std::string& phylogeny_file_name );
	void _drop_saturated(const double upper_filt_thres);
	void _drop_rare(const double filt_thres);
        int _drop_empty_alphas();
	int _drop_saturated_alphas(const double upper_filt_thres);
	int _drop_rare_alphas(const double filt_thres);
        int _drop_empty_betas();
        int _drop_empty_gammas();
        
        void _dump_sizes() const;
        void _warn_superfluous_error( const Alpha& alpha );
        void _warn_superfluous_error( const Beta& beta );
        void _warn_superfluous_error( const Gamma& gamma );
};


#endif //COINFINDER_FAMILYSET_H
