//
// Created by Martin Rusilowicz on 22/08/2017.
//

#include <cassert>
#include <sstream>
#include "elements.h"


/**
 * Constructor
 * @param name 
 */
Alpha::Alpha( const std::string& name )
        : _name( name )
          , _index( -1 )
	  ,D()
#ifndef  NDEBUG
        , _gammas()
#else
        , _num_gammas( 0 )
#endif
        , _edges()
	, _coincident_edges()
{
    // pass
}


/**
 * Registers an edge to `beta`
 * @param beta 
 * @return Returns if this is the FIRST edge from `this` `Gamma` to the specified `beta`. 
 */
bool Alpha::register_edge( const Gamma* gamma, const Beta& beta )
{
    auto it = this->_edges.find( &beta );

    if (it == _edges.end())
    {
        _edges[ &beta ] = 1;
        return true;
    }
    else
    {
        _edges[ &beta ] = it->second + 1;
        return false;
    }

    ++this->_num_edges;
}

bool Alpha::register_coincident_edge( const Alpha& alpha )
{
	auto it = this->_coincident_edges.find( &alpha );

	if (it == _coincident_edges.end()) {
		_coincident_edges[ &alpha ] = 1;
		++this->_num_coincident_edges;
		return true;
	} else {
		_coincident_edges[ &alpha ] = it->second + 1;
		return false;
	}
}

void Alpha::register_D( double D ) {
	this->D = D;
}

/**
 * Indicates the number of gammas in `this` group.
 * @return 
 */
int Alpha::get_num_gammas() const
{
#ifndef NDEBUG
    return static_cast<int>(this->_gammas.size());
#else
    return this->_num_gammas;
#endif
}

int Alpha::get_num_edges() const {
	return (this->_edges).size();
}

int Alpha::get_num_coincident_edges() const {
	return (this->_coincident_edges).size();
}


/**
 * Returns the map representing the edges (see `_edges`). 
 * @return 
 */
const std::map<const Beta*, int>& Alpha::get_edges() const
{
    return _edges;
}

const std::map<const Alpha*, int>& Alpha::get_coincident_edges() const
{
   return _coincident_edges;
}


/**
 * Obtains the name of this alpha group.
 * @return 
 */
const std::string& Alpha::get_name() const
{
    return _name;
}


/**
 * Adds a gamma to my collection.
 * 
 * ASSUMPTIONS: The Gamma has not already been added to my collection.
 *              No error is returned if this is not true, but it will break the calculation.
 * 
 * @param gamma 
 */
void Alpha::register_gamma( const Gamma& gamma )
{
#ifndef NDEBUG
    auto it = this->_gammas.find( &gamma );

    if (it != this->_gammas.end())
    {
        std::stringstream ss;
        ss << "The alpha group '" << this->get_name() << "' already has the gamma element '" << gamma.get_name() << "' amongst its " << this->_gammas.size() << " elements.";
        throw std::logic_error( ss.str());
    }

    this->_gammas.insert( &gamma );
#else
    ++this->_num_gammas;
#endif
}


/**
 * Removes a gamma from my collection.
 * 
 * ASSUMPTIONS: The Gamma has been added to my collection.
 *              No error is returned if this is not true, but it will break the calculation.
 * 
 * @param gamma 
 */
void Alpha::unregister_gamma( const Gamma& gamma )
{
#ifndef NDEBUG
    auto it = this->_gammas.find( &gamma );
    assert( it != this->_gammas.end());
    this->_gammas.erase( it );
#else
    --this->_num_gammas;
#endif
}


/**
 * Costructor
 * @param name 
 */
Beta::Beta( const std::string& name )
        : _name( name )
#ifndef NDEBUG
        , _gammas()
#else
        , _num_gammas( 0 )
#endif
        , _has_edges( false )
	, _coin_status( CoinStatus::NEITHER )
{
    // pass
}


/**
 * Obtains the number of edges this `Beta` group has. 
 * @return 
 */
int Beta::get_num_gammas() const
{
#ifndef NDEBUG
    return static_cast<int>(this->_gammas.size());
#else
    return this->_num_gammas;
#endif
}


/**
 * Obtains the name of this `Beta` group.
 * @return 
 */
const std::string& Beta::get_name() const
{
    return _name;
}


/**
 * Registers the gamma as a element of this beta.
 * @param gamma 
 */
void Beta::register_gamma( const Gamma& gamma )
{
#ifndef NDEBUG
    auto it = this->_gammas.find( &gamma );

    if (it != this->_gammas.end())
    {
        std::stringstream ss;
        ss << "The beta group '" << this->get_name() << "' already has the gamma element '" << gamma.get_name() << "' amongst its " << this->_gammas.size() << " elements.";
        throw std::logic_error( ss.str());
    }

    this->_gammas.insert( &gamma );
#else
    ++_num_gammas; // For efficiency, only the count is retained, the API is for consistency.
#endif
}


/**
 * Registers the gamma as a element of this beta.
 * @param gamma 
 */
void Beta::unregister_gamma( const Gamma& gamma )
{
#ifndef NDEBUG
    auto it = this->_gammas.find( &gamma );
    assert ( it != this->_gammas.end());
    this->_gammas.erase( it );
#else
    --_num_gammas; // For efficiency, only the count is retained, the API is for consistency.
#endif
}


/**
 * Registers the existence of an edge on this `Beta` group.
 * @param gamma 
 * @param alpha 
 */
void Beta::register_edge( const Gamma* gamma, const Alpha& alpha )
{
    _has_edges = true; // For efficiency, only the predicate is retained, the API is for consistency.
}


/**
 * Indicates if this `Beta` group possesses any edges. 
 * @return 
 */
bool Beta::has_edges() const
{
    return _has_edges;
}

/**
 * Constructor
**/
Edge::Edge( const std::string& name ) //, int num_gammas, int weight )
	: _alphas()
	 , _betas()
	 , num_gammas()
	 , weight()
#ifndef NBEBUG
	//, _name() //( name )
#endif
{
	// pass
}

void Edge::register_nodes( const Alpha& alpha, const Beta& beta )
{
	this->_alphas.insert( &alpha );
	this->_betas.insert( &beta );
}

const std::string Edge::get_name() const
{
#ifndef NDEBUG
    	return this->_name;
#else
    	std::stringstream ss; 
    	ss << this;
    	return ss.str();
#endif
}

int Edge::get_num_gammas() const
{
	return num_gammas;
}

void Edge::set_weight( const int weight )
{
	this->weight = weight;
}

int Edge::get_weight() const
{
	return weight;
}

/**
 * Constructor
 */
Gamma::Gamma( const std::string& name ) // Note that the name is not stored (to save memory) but we need to accept the parameter to because Gamma shares a common interface with Alpha and Beta
        : _has_edges( false )
          , _betas()
          , _alphas()
#ifndef NDEBUG
        , _name( name )
#endif
{
    // pass
}


/**
 * Registers the beta as an owner of this gamma.
 * @param beta 
 */
void Gamma::register_beta( Beta& beta )
{
#ifndef NDEBUG
    auto it = this->_betas.find( &beta );

    if (it != this->_betas.end())
    {
        std::stringstream ss;
        ss << "The gamma element '" << this->get_name() << "' is already registered as a member of the beta group '" << beta.get_name() << "'.";
        throw std::logic_error( ss.str());
    };
#endif

    this->_betas.insert( &beta );
}


/**
 * Registers the alpha as an owner of this gamma.
 * @param beta 
 */
void Gamma::register_alpha( Alpha& alpha )
{
    this->_alphas.insert( &alpha );
}


/**
 * Obtains the `Beta` groups `this` `Gamma` is a member of.
 * @return 
 */
const std::unordered_set<Beta*>& Gamma::get_betas() const
{
    return _betas;
}


/**
 * Obtains the `Alpha` groups `this` `Gamma` is a member of.
 * @return 
 */
const std::unordered_set<Alpha*>& Gamma::get_alphas() const
{
    return _alphas;
}


/**
 * Obtains the `Beta` groups `this` `Gamma` is a member of.
 * @return 
 */
std::unordered_set<Beta*>& Gamma::get_betas()
{
    return _betas;
}


/**
 * Registers the existence of an edge spanning `this` `Gamma`.
 * @param alpha 
 * @param beta 
 */
void Gamma::register_edge( const Alpha& alpha, const Beta& beta )
{
    _has_edges = true;

}


/**
 * Indicates if any edges span `this` `Gamma`. 
 * @return 
 */
bool Gamma::has_edges() const
{
    return _has_edges;
}


const std::string Gamma::get_name() const
{
#ifndef NDEBUG
    return this->_name;
#else
    std::stringstream ss;
    ss << this;
    return ss.str();
#endif
}
