//
// Created by Martin Rusilowicz on 22/08/2017.
//

#ifndef COINFINDER_ELEMENTS_H
#define COINFINDER_ELEMENTS_H

#include <string>
#include <map>   
#include <unordered_set>
#include "constants.h"

class Alpha;


class Beta;


class Edge;


class Gamma;


/**
 * The ALPHA grouping.
 * e.g. Gene families, experimental groups
 */
class Alpha
{
    protected:
	static int nextIndex;
    private:
        const std::string _name;
        int               _index;
	double		  D;
#ifndef NDEBUG
        std::unordered_set<const Gamma*> _gammas;
#else
        int _num_gammas;
	int _num_edges;
	int _num_coincident_edges;
#endif
        std::map<const Beta*, int> _edges; // beta, to the count of gammas in the edge (i.e. the weight)
	std::map<const Alpha*, double> _coincident_edges;

    public:
        explicit Alpha( const std::string& name );

        const std::string& get_name() const;
	const int get_index() const;

        bool register_edge( const Gamma* gamma, const Beta& beta );
	bool register_coincident_edge( const Alpha& alpha, double p_value );
	void register_D( double D );
        void register_gamma( const Gamma& gamma );
        void unregister_gamma( const Gamma& gamma );

        int get_num_gammas() const;
	int get_num_edges() const;
	int get_num_coincident_edges() const;
	double get_D() const;
	Alpha* get_alpha( std::string );
        const std::map<const Beta*, int>& get_edges() const;
	const std::map<const Alpha*, double>& get_coincident_edges() const;
};
//int Alpha::nextIndex = 0;


/**
 * The BETA grouping.
 * e.g. Taxa, Go-terms
 */
class Beta
{
    private:
        const std::string _name;
#ifndef NDEBUG
        std::unordered_set<const Gamma*> _gammas;
#else
        int _num_gammas;    // Number of gammas this beta contains
#endif
        bool _has_edges;
	CoinStatus _coin_status;

    public:
        explicit Beta( const std::string& name );
        int get_num_gammas() const;
        void register_gamma( const Gamma& gamma );
        void unregister_gamma( const Gamma& gamma );
        const std::string& get_name() const;
        void register_edge( const Gamma* gamma, const Alpha& alpha );
        bool has_edges() const;
};

/**
 * EDGE class.
 * e.g. connections between alphas and betas, weights of the edges.
**/
class Edge
{
        private:
#ifndef NDEBUG
        std::string _name;
#endif
        std::unordered_set<const Alpha*> _alphas;
	std::unordered_set<const Beta*> _betas;
	int num_gammas;
	int weight;

        public:
        explicit Edge( const std::string& name ); //, int num_gammas, int weight);
        void register_nodes( const Alpha& alpha, const Beta& beta );
	const std::string get_name() const;
        int get_num_gammas() const;
        void set_weight( const int weight );
	int get_weight() const;
};

/**
 * The GAMMA, which are the elements of both the `Alpha` and `Beta` groups.
 * e.g. Genes or sequences
 */
class Gamma
{
    private:
#ifndef NDEBUG
        std::string                 _name;
        #endif
        bool                       _has_edges;
        std::unordered_set<Alpha*> _alphas; // the alphas that this gamma is a member of. only needed if we are removing stuff.
        std::unordered_set<Beta*>  _betas; // the betas that this gamma is a member of

    public:
        explicit Gamma( const std::string& name );
        void register_edge( const Alpha& alpha, const Beta& beta );
        bool has_edges() const;

        void register_alpha( Alpha& alpha );
        void register_beta( Beta& beta );
        const std::unordered_set<Alpha*>& get_alphas() const;
        const std::unordered_set<Beta*>& get_betas() const;
        std::unordered_set<Beta*>& get_betas();
        const std::string get_name() const;


};


#endif //COINFINDER_ELEMENTS_H
