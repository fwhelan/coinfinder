//
// Created by Martin Rusilowicz on 17/08/2017.
//

#ifndef COINFINDER_FAMILYSET_H
#define COINFINDER_FAMILYSET_H

#include <string>
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
	id_lookup<Edge> _edges;
        
    public:
        explicit DataSet(const TParameters& options);
        ~DataSet();
        
        const id_lookup<Alpha>& get_alphas() const;
        const id_lookup<Beta>& get_betas() const;
        const id_lookup<Gamma>& get_gammas() const;
	const id_lookup<Edge>& get_edges() const;

        void read_files( const std::string& alpha_file_name, const std::string& beta_file_name, const std::string& combined_file_name );
        
        int get_num_edges() const;
        const TParameters& get_options() const;

        
    private:
        void _read_alpha_file( const std::string& file_name );
        void _read_beta_file( const std::string& file_name );
        void _read_combined_file( const std::string& file_name );
        
        void _drop_empty();
        int _drop_empty_alphas();
        int _drop_empty_betas();
        int _drop_empty_gammas();
        
        void _dump_sizes() const;
        void _warn_superfluous_error( const Alpha& alpha );
        void _warn_superfluous_error( const Beta& beta );
        void _warn_superfluous_error( const Gamma& gamma );
};


#endif //COINFINDER_FAMILYSET_H
