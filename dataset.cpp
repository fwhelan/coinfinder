//
// Created by Martin Rusilowicz on 17/08/2017.
//

#include "dataset.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include "elements.h"
#include "parameters.h"


/**
 * Constructor
 */
DataSet::DataSet( const TParameters& options )
        : _options( options )
          , _alphas()
          , _betas()
          , _gammas()
          , _num_edges( 0 )
	  , _edges()
{
    // pass
}


/**
 * Destructor.
 */
DataSet::~DataSet()
{
    // pass
}


/**
 * Reads the BETA file.
 * This is the first file to be loaded.
 * @param file_name 
 */
void DataSet::_read_beta_file( const std::string& file_name )
{
    std::cerr << "Reading beta groupings..." << std::endl;

    std::ifstream file_in;
    file_in.open( file_name );
    std::string cell;
    bool        left = true;
    Gamma* gamma = nullptr;
    
    if (!file_in)
    {
        std::stringstream ss;
        ss << "Failed to open file: " << file_name; 
        throw std::logic_error(ss.str());
    }

    while (getline( file_in, cell, left ? static_cast<char>('\t') : static_cast<char>('\n')))
    {
        if (left)
        {
            gamma = &this->_gammas.find_id( cell );
        }
        else
        {
            // BETA--[CONTAINS]--GAMMA 
            Beta& beta = this->_betas.find_id( cell );
            gamma->register_beta( beta ); // Γ -> B
            beta.register_gamma( *gamma ); // B -> Γ
            
            if (_options.verbose)
            {
                std::cerr << "Beta " << beta.get_name() << " has gamma " << gamma->get_name() << std::endl;
            }
        }

        left = !left;
    }

    this->_dump_sizes();
}

/**
 * Reads the ALPHA file.
 * This is always done AFTER the beta file.
 * @param file_name 
 */
void DataSet::_read_alpha_file( const std::string& file_name )
{
    std::cerr << "Reading alpha groupings..." << std::endl;

    std::ifstream file_in;
    file_in.open( file_name );
    std::string cell;
    bool        left = true;
    Gamma* gamma = nullptr;
    
    if (!file_in)
    {
        std::stringstream ss;
        ss << "Failed to open file: " << file_name; 
        throw std::logic_error(ss.str());
    }

    while (getline( file_in, cell, left ? static_cast<char>('\t') : static_cast<char>('\n')))
    {
        if (left)
        {
            gamma = &this->_gammas.find_id( cell );
        }
        else
        {
            // ALPHA--[CONTAINS]->GAMMA
            Alpha& alpha = this->_alphas.find_id( cell );

            if (this->get_options().permit_filter)
            {
                gamma->register_alpha( alpha ); // Γ -> A
            }

            alpha.register_gamma( *gamma ); // A -> Γ
            
            if (_options.verbose)
            {
                std::cerr << "Alpha " << alpha.get_name() << " has gamma " << gamma->get_name() << std::endl;
            }

            for (Beta* beta : gamma->get_betas())
            {
                // ALPHA--[CONTAINS]->GAMMA<-[CONTAINS]--BETA
                if (alpha.register_edge( gamma, *beta )) // A -> E
                {
                    ++this->_num_edges;
                }

                beta->register_edge( gamma, alpha ); // B -> E
                gamma->register_edge( alpha, *beta ); // Γ -> E
                
                if (_options.verbose)
                {
                    std::cerr << "Edge formed between " << alpha.get_name() << " and " << beta->get_name() << " using gamma " << gamma->get_name() << std::endl;
                }
            }
        }

        left = !left;
    }

    this->_dump_sizes();
}

/**
 * Reads the ALPHA-BETA file.
 * This is always done in lieu of reading the alpha AND beta files.
 * @param file_name 
 */
void DataSet::_read_combined_file( const std::string& file_name )
{
	//add a bit to this about reading in a separate bit of information, a middle column, that contains the edge information
    std::cerr << "Reading alpha-beta edges..." << std::endl;

    std::ifstream file_in;
    file_in.open( file_name );
    std::string cell;
    bool        left = true;
    bool	middle = false;
    bool	right = false;
    Alpha* alpha = nullptr; 
    Beta* beta = nullptr;;
    Edge* edge = nullptr;
    
    if (!file_in)
    {
        std::stringstream ss;
        ss << "Failed to open file: " << file_name; 
        throw std::logic_error(ss.str());
    }
    
    int n = 0;

    while (getline( file_in, cell, left ? static_cast<char>('\t') : ( middle ? static_cast<char>('\t') : static_cast<char>('\n'))))
    {
        if (left)
        {
            alpha = &this->_alphas.find_id( cell );
		left = false;
		middle = true;
        }
	else if (middle)
        {
            // ALPHA--[WITH]->BETA
            beta = &this->_betas.find_id( cell );
            
            if (alpha->register_edge( nullptr, *beta )) //include that extra bit of information here about the edge
            {
                ++this->_num_edges;
            }
            
            beta->register_edge( nullptr, *alpha ); // B -> E
            
            if (_options.verbose)
            {
                std::cerr << "Edge formed between " << alpha->get_name() << " and " << beta->get_name() << std::endl;
            }
		middle = false;
		right = true;
        }
	else
	{
		edge = &this->_edges.find_id( alpha->get_name()+"-"+beta->get_name() );

		edge->register_nodes( *alpha, *beta );

		edge->set_weight( std::stoi(cell) );

		if (_options.verbose)
		{
			std::cerr << "Weight added " << cell << std::endl;
			std::cerr << "Weight is " << edge->get_weight() << std::endl;
		}
		right = false;
		left = true;
	}

        n+=1;
    }

    this->_dump_sizes();
}


/**
 * Obtains the gamma collection
 * @return 
 */
const id_lookup<Gamma>& DataSet::get_gammas() const
{
    return _gammas;
}


/**
 * Obtains the beta collection
 * @return 
 */
const id_lookup<Beta>& DataSet::get_betas() const
{
    return _betas;
}

const id_lookup<Edge>& DataSet::get_edges() const
{
	return _edges;
}


/**
 * Drops all alphas that reference no betas (through gammas). 
 * @return Number of items dropped. 
 */
int DataSet::_drop_empty_alphas()
{
    std::vector<std::string> to_drop = std::vector<std::string>();

    std::map<std::string, Alpha*>& table = this->_alphas.get_table();

    for (const auto& kvp : table)
    {
        const Alpha& alpha = *kvp.second;

        if (alpha.get_edges().empty())
        {
            const std::string& name = kvp.first;
            to_drop.push_back( name );
        }
    }

    for (const std::string& name : to_drop)
    {
        if (_options.verbose)
        {
            std::cerr << "Deleting superfluous alpha group '" << name << "'." << std::endl;
        }

        auto it = table.find( name );

        if (it == table.end())
        {
            std::stringstream ss;
            ss << "An internal error occurred, please submit a bug report. Error details: Failed to recall an alpha group for deletion. Name = '" << name << "'.";
            throw std::logic_error( ss.str());
        }

        Alpha* alpha = it->second;

        if (alpha->get_num_gammas())
        {
            std::stringstream ss;
            ss << "An internal error occurred, please submit a bug report. Error details: Alpha's gammas weren't removed earlier. Name = '" << name << "'. Count = " << alpha->get_num_gammas() << ".";
            throw std::logic_error( ss.str());
        }

        delete alpha;

        table.erase( it );
    }

    return static_cast<int>(to_drop.size());
}


/**
 * Drops all betas that reference no alphas (through gammas)
 * @return Number of items dropped
 */
int DataSet::_drop_empty_betas()
{
    std::vector<std::string> to_drop = std::vector<std::string>();
    std::map<std::string, Beta*>& table = this->_betas.get_table();

    for (auto kvp : table)
    {
        const Beta& beta = *kvp.second;

        if (!beta.has_edges())
        {
            const std::string& name = kvp.first;
            to_drop.push_back( name );
        }
    }

    for (const std::string& name : to_drop)
    {
        if (_options.verbose)
        {
            std::cerr << "Deleting superfluous beta group '" << name << "'." << std::endl;
        }

        auto it = table.find( name );

        if (it == table.end())
        {
            std::stringstream ss;
            ss << "An internal error occurred, please submit a bug report. Error details: Failed to recall a beta group for deletion. Name = '" << name << "'.";
            throw std::logic_error( ss.str());
        }

        Beta* beta = it->second;

        if (beta->get_num_gammas())
        {
            std::stringstream ss;
            ss << "An internal error occurred, please submit a bug report. Error details: Beta's gammas weren't removed earlier. Name = '" << name << "'. Count = " << beta->get_num_gammas() << ".";
            throw std::logic_error( ss.str());
        }

        delete beta;
        table.erase( it );
    }

    return static_cast<int>(to_drop.size());
}


/**
 * Drops all gammas that don't form a link between alphas and betas.
 * @return Number of items dropped
 */
int DataSet::_drop_empty_gammas()
{
    std::vector<std::string> to_drop = std::vector<std::string>();
    std::map<std::string, Gamma*>& table = this->_gammas.get_table();

    for (auto kvp : table)
    {
        const Gamma& gamma = *kvp.second;

        if (!gamma.has_edges())
        {
            to_drop.push_back( kvp.first );
        }
    }

    for (const std::string& name : to_drop)
    {
        auto it = table.find( name );

        if (it == table.end())
        {
            std::stringstream ss;
            ss << "An internal error occurred, please submit a bug report. Error details: Failed to recall a beta group for deletion. Name = '" << name << "'.";
            throw std::logic_error( ss.str());
        }
        
        if (_options.verbose)
        {
            std::cerr << "Deleting superfluous gamma element '" << name << "'." << std::endl;
        }

        const Gamma& gamma = *it->second;

        DataSet::_warn_superfluous_error( gamma );

        // Remove from counts for my betas/alphas
        for (Alpha* alpha : gamma.get_alphas())
        {
            alpha->unregister_gamma( gamma );
        }

        for (Beta* beta : gamma.get_betas())
        {
            beta->unregister_gamma( gamma );
        }

        // Notice
        // ------
        // We know that this gamma doesn't form an edge.
        // We just go ahead and drop it from the dataset.
        // Why do we even keep a record of the Gammas?
        table.erase( it );
    }

    return static_cast<int>(to_drop.size());
}


/**
 * Gets the alpha collection
 * @return 
 */
const id_lookup<Alpha>& DataSet::get_alphas() const
{
    return this->_alphas;
}


/**
 * Prints the sizes of the collections to STD.OUT.
 */

void DataSet::_dump_sizes() const
{
    std::cerr << "- n.ALPHA = " << this->_alphas.size() << std::endl;
    std::cerr << "- n.BETA  = " << this->_betas.size() << std::endl;
    std::cerr << "- n.GAMMA = " << this->_gammas.size() << std::endl;
    std::cerr << "- n.EDGES = " << this->_num_edges << std::endl;
}


/**
 * Obtains the number of links between alphas and betas (regardless of their weight)
 * @return 
 */
int DataSet::get_num_edges() const
{
    return _num_edges;
}


/**
 * Drops elements that don't form edges from the collection.
 */
void DataSet::_drop_empty()
{
    std::cerr << "Dropping empty sets..." << std::endl;
    int gamma_dropped = this->_drop_empty_gammas();
    int alpha_dropped = this->_drop_empty_alphas();
    int beta_dropped  = this->_drop_empty_betas();

    if (alpha_dropped || beta_dropped || gamma_dropped)
    {
        std::cerr << "Warning: Superfluous data has been dropped!" << std::endl;
        std::cerr << "- d.ALPHA = " << -alpha_dropped << std::endl;
        std::cerr << "- d.BETA  = " << -beta_dropped << std::endl;
        std::cerr << "- d.GAMMA = " << -gamma_dropped << std::endl;
        this->_dump_sizes();
    }
    else
    {
        std::cerr << "Nothing dropped, your data is good to go." << std::endl;
    }
}


/**
 * Reads in the input files
 * @param alpha_file_name 
 * @param beta_file_name 
 */
void DataSet::read_files( const std::string& alpha_file_name, const std::string& beta_file_name, const std::string& combined_file_name )
{
    if (combined_file_name.empty())
    {
        this->_read_beta_file( beta_file_name );
        this->_read_alpha_file( alpha_file_name );
    }
    else
    {
        this->_read_combined_file( combined_file_name );
    }
    
    this->_drop_empty();
}


/**
 * Returns the options used to create this dataset.
 * @return 
 */
const TParameters& DataSet::get_options() const
{
    return _options;
}


/**
 * Throws a superfluous data error if the hadn't flagged that they expected this.
 * @param alpha 
 */
void DataSet::_warn_superfluous_error( const Alpha& alpha )
{

    if (!_options.permit_filter)
    {
        std::stringstream ss;
        ss << "The alpha group called '" << alpha.get_name() << "' forms no links with any beta group, but it does have " << alpha.get_num_gammas() << " gamma member(s) (I can't tell you which because I didn't store that information earlier to save space). This may be indicative that you are analysing unnecessary data, or have specified data-sets which are unrelated. If you are sure you are doing the right thing, please use the 'p' flag to permit superfluous data removal.";
        throw std::logic_error( ss.str().c_str());
    }
}


/**
 * Throws a superfluous data error if the hadn't flagged that they expected this.
 * @param alpha 
 */
void DataSet::_warn_superfluous_error( const Beta& beta )
{
    if (!_options.permit_filter)
    {
        std::stringstream ss;
        ss << "The beta group called '" << beta.get_name() << "' forms no links with any alpha group, but it does have " << beta.get_num_gammas() << " gamma member(s) (I can't tell you which because I didn't store that information earlier to save space). This may be indicative that you are analysing unnecessary data, or have specified data-sets which are unrelated. If you are sure you are doing the right thing, please use the 'p' flag to permit superfluous data removal.";
        throw std::logic_error( ss.str().c_str());
    }
}


/**
 * Throws a superfluous data error if the hadn't flagged that they expected this.
 * @param alpha 
 */
void DataSet::_warn_superfluous_error( const Gamma& gamma )
{
    if (!_options.permit_filter)
    {
        std::stringstream ss;
        ss << "A gamma member, the name of which I have not stored, does not form links between any alpha and beta groups. This may be indicative that you are analysing unnecessary data, or have specified data-sets which are unrelated. If you are sure you are doing the right thing, please use the 'p' flag to permit superfluous data removal.";
        throw std::logic_error( ss.str().c_str());
    }
}
