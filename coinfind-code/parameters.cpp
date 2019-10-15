//
// Created by Martin Rusilowicz on 22/08/2017.
//

#include <iostream>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>
#include <functional>
#include "parameters.h"
#include "test_cases.h"


#pragma clang diagnostic push
#pragma ide diagnostic ignored "ImplicitIntegerAndEnumConversion"


TParameters::TParameters()
        : sig_level( 0.05 )
          , correction( ECorrection::BONFERRONI ) //ECorrection::_INVALID )
          , alt_hypothesis( EHypothesis::GREATER ) //EHypothesis::_INVALID )
          , method( EMethod::COINCIDENCE ) //EMethod::_INVALID )
          , coin_set_mode( ESetMode::FULL ) //ESetMode::_INVALID )
          , coin_max_mode( EMaxMode::_INVALID )
          , verbose( false )
          , permit_filter( false )
          , alpha_name( "Alpha" )
          , beta_name( "Beta" )
          , gamma_name( "Gamma" )
          , alpha_file_name( "" )
          , beta_file_name( "" )
          , output_all( false )
	  , upper_filt_thres( 1.0 )
	  , filt_thres( 0.05 )
          , combined_file_name( "" )
	  , num_cores ( 2 )
	  , Rmsgs ( false )
	  , prefix ( "coincident" )
	  , roary( false )
{
    // pass
}


struct CommandArg
{
    std::string           _long_name;
    char                  _short_name;
    std::function<void()> _action;


    CommandArg( std::string long_name, char short_name, std::function<void()> action )
            : _long_name( std::move( long_name ))
              , _short_name( short_name )
              , _action( std::move( action ))
    {

    }
};


template<typename T>
void local_assert( T& target, T value, const char* name )
{
    if (target != T::_INVALID)
    {
        std::stringstream ss;
        ss << "More than one value for the '" << name << "' setting has been specified." << std::endl;
        throw std::logic_error( ss.str());
    }

    target = value;
}


void local_assert( bool& target, bool value, const char* name )
{
    if (target)
    {
        std::stringstream ss;
        ss << "More than one value for the '" << name << "' setting has been specified." << std::endl;
        throw std::logic_error( ss.str());
    }

    target = value;
}


#define SET_ONCE( target, value ) local_assert(result.target, value, #target)


TParameters TParameters::parse( int arg_count, const char** arg_vals )
{
    TParameters result = TParameters();

    enum class ECommand
    {
            FLAGS,
            ALPHA_FN,
            BETA_FN,
            COMBINED_FN,
            ALPHA_N,
            BETA_N,
            GAMMA_N,
            SIG_LEVEL,
	    PREFIX,
	    PHYLO,
	    CORES,
	    UPFILT,
	    FILT,
            QUERY,
    };

    ECommand next_command = ECommand::FLAGS;

    std::vector<CommandArg> commands = std::vector<CommandArg>();

    commands.emplace_back( "--nocorrection", 'n', [ &result ]() { SET_ONCE( correction, ECorrection::NONE ); } );
    commands.emplace_back( "--bonferroni", 'm', [ &result ]() { SET_ONCE( correction, ECorrection::BONFERRONI ); } );
    commands.emplace_back( "--fraction", 'c', [ &result ]() { SET_ONCE( correction, ECorrection::FRACTION ); } );
    //commands.emplace_back( "--coincidence", 'c', [ &result ]() { SET_ONCE( method, EMethod::COINCIDENCE ); } );
    //commands.emplace_back( "--connectivity", 'e', [ &result ]() { SET_ONCE( method, EMethod::CONNECTIVITY ); } );
    commands.emplace_back( "--greater", 'g', [ &result ]() { SET_ONCE( alt_hypothesis, EHypothesis::GREATER ); } );
    commands.emplace_back( "--less", 'l', [ &result ]() { SET_ONCE( alt_hypothesis, EHypothesis::LESS ); } );
    commands.emplace_back( "--twotailed", 't', [ &result ]() { SET_ONCE( alt_hypothesis, EHypothesis::TWOTAILED ); } );
    //commands.emplace_back( "--full", 'f', [ &result ]() { SET_ONCE( coin_set_mode, ESetMode::FULL ); } );
    //commands.emplace_back( "--union", 'u', [ &result ]() { SET_ONCE( coin_set_mode, ESetMode::INTERSECTION ); } );
    commands.emplace_back( "--associate", 'a', [ &result ]() { SET_ONCE( coin_max_mode, EMaxMode::ACCOMPANY ); } );
    commands.emplace_back( "--dissociate", 'd', [ &result ]() { SET_ONCE( coin_max_mode, EMaxMode::AVOID ); } );
    commands.emplace_back( "--verbose", 'v', [ &result ]() { SET_ONCE( verbose, true ); } );
    commands.emplace_back( "--filter", 'r', [ &result ]() { SET_ONCE( permit_filter, true ); } );
    commands.emplace_back( "--Rmsgs", 'R', [ &result ]() { SET_ONCE( Rmsgs, true); } );
    commands.emplace_back( "--upfilthreshold", 'U', [ &next_command ]() { next_command = ECommand::UPFILT; } );
    commands.emplace_back( "--filthreshold", 'F', [ &next_command ]() { next_command = ECommand::FILT; } );
    commands.emplace_back( "--level", 'L', [ &next_command ]() { next_command = ECommand::SIG_LEVEL; } );
    //commands.emplace_back( "--a", 'a', [ &next_command ]() { next_command = ECommand::ALPHA_FN; } );
    //commands.emplace_back( "--b", 'b', [ &next_command ]() { next_command = ECommand::BETA_FN; } );
    commands.emplace_back( "--input", 'i', [ &next_command ]() { next_command = ECommand::COMBINED_FN; } );
    commands.emplace_back( "--inputroary", 'I', [ &result ]() { SET_ONCE( roary, true); } );
    //commands.emplace_back( "--alpha", 'A', [ &next_command ]() { next_command = ECommand::ALPHA_N; } );
    //commands.emplace_back( "--beta", 'B', [ &next_command ]() { next_command = ECommand::BETA_N; } );
    //commands.emplace_back( "--gamma", 'C', [ &next_command ]() { next_command = ECommand::GAMMA_N; } );
    commands.emplace_back( "--output", 'o', [ &next_command ]() { next_command = ECommand::PREFIX; } );
    commands.emplace_back( "--phylogeny", 'p', [ &next_command ]() { next_command = ECommand::PHYLO; } );
    commands.emplace_back( "--num_cores", 'x', [ &next_command ]() { next_command = ECommand::CORES; } );
    commands.emplace_back( "--query", 'q', [ &next_command ]() { next_command = ECommand::QUERY; } );
    commands.emplace_back( "--all", 'E', [ &result ]() { SET_ONCE( output_all, true ); } );
    commands.emplace_back( "--test", 'T', &test_cases );


    for (int n = 1; n < arg_count; ++n)
    {
        std::string arg = arg_vals[ n ];

        switch (next_command)
        {
            case ECommand::FLAGS:
            {
#ifndef NDEBUG
                std::cerr << "command-parser: Flags '" << arg << "' specified..." << std::endl;
#endif

                bool okay = false;

                for (const CommandArg& cmd : commands)
                {
                    if (arg == cmd._long_name)
                    {
#ifndef NDEBUG
                        std::cerr << "command-parser: Long name of '" << cmd._long_name << "' specified." << std::endl;
#endif
                        cmd._action();
                        okay = true;
                        break;
                    }
                }

                if (!okay)
                {
                    for (const char& f : arg)
                    {
                        for (const CommandArg& cmd : commands)
                        {
                            if (f == cmd._short_name)
                            {
#ifndef NDEBUG
                                std::cerr << "command-parser: Short name of '" << cmd._long_name << "' specified." << std::endl;
#endif

                                cmd._action();
                                okay = true;
                                break;
                            }
                        }
                    }
                }

                if (!okay)
                {
                    std::stringstream ss;
                    ss << "Unrecognised command-like argument '" << arg << "'." << std::endl;
                    throw std::logic_error( ss.str());
                }

                break;
            }

            case ECommand::ALPHA_FN:
#ifndef NDEBUG
                std::cerr << "command-parser: Alpha filename '" << arg << "' specified." << std::endl;
#endif

                result.alpha_file_name = arg;
                next_command = ECommand::FLAGS;
                break;

            case ECommand::BETA_FN:
#ifndef NDEBUG
                std::cerr << "command-parser: Beta filename '" << arg << "' specified." << std::endl;
#endif

                result.beta_file_name = arg;
                next_command = ECommand::FLAGS;
                break;

            case ECommand::COMBINED_FN:
#ifndef NDEBUG
                std::cerr << "command-parser: Beta filename '" << arg << "' specified." << std::endl;
#endif

                result.combined_file_name = arg;
                next_command = ECommand::FLAGS;
                break;

            case ECommand::ALPHA_N:
#ifndef NDEBUG
                std::cerr << "command-parser: Alpha name '" << arg << "' specified." << std::endl;
#endif

                result.alpha_name = arg;
                next_command = ECommand::FLAGS;
                break;

            case ECommand::BETA_N:
#ifndef NDEBUG
                std::cerr << "command-parser: Beta name '" << arg << "' specified." << std::endl;
#endif

                result.beta_name = arg;
                next_command = ECommand::FLAGS;
                break;

            case ECommand::GAMMA_N:
#ifndef NDEBUG
                std::cerr << "command-parser: Gamma name '" << arg << "' specified." << std::endl;
#endif

                result.gamma_name = arg;
                next_command = ECommand::FLAGS;
                break;

            case ECommand::SIG_LEVEL:
#ifndef NDEBUG
                std::cerr << "command-parser: Significance level '" << arg << "' specified." << std::endl;
#endif

                result.sig_level = std::stof( arg );
                next_command = ECommand::FLAGS;
                break;

	    case ECommand::PREFIX:
#ifndef NDEBUG
		std::cerr << "command-parser: Output prefix '" << arg << "' specified." << std::endl;
#endif
		result.prefix = arg;
		next_command = ECommand::FLAGS;
		break;

	    case ECommand::PHYLO:
#ifndef NDEBUG
		std::cerr << "command=parser: Phylogeny '" << arg << "' specified." << std::endl;
#endif
		result.phylogeny = arg;
		next_command = ECommand::FLAGS;
		break;

	    case ECommand::CORES:
#ifndef NDEBUG
		std::cerr << "command=parser: Number of cores '" << arg << "' specified." << std::endl;
#endif
		result.num_cores = stoi(arg);
		next_command = ECommand::FLAGS;
		break;

	    case ECommand::UPFILT:
#ifndef NDEBUG
		std::cerr << "command-parser: Upper filter threshold '" << arg << "' specified." << std::endl;
#endif
		result.upper_filt_thres = std::stod(arg);
		next_command = ECommand::FLAGS;
		break;

	    case ECommand::FILT:
#ifndef NDEBUG
		std::cerr << "command-parser: Filter threshold '" << arg << "' specified." << std::endl;
#endif
		result.filt_thres = std::stod(arg);
		next_command = ECommand::FLAGS;
		break;

            case ECommand::QUERY:
#ifndef NDEBUG
                std::cerr << "command-parser: Query ID '" << arg << "' specified." << std::endl;
#endif

                result.deep_query_alpha = arg;
                next_command = ECommand::FLAGS;
                break;
        }
    }

    return result;
}


void TParameters::print_and_assert() const
{
    switch (this->correction)
    {
        case ECorrection::BONFERRONI:
            std::cerr << "> CORRECTION ······· = BONFERRONI" << std::endl;
            break;

        case ECorrection::NONE:
            std::cerr << "> CORRECTION ······· = NONE" << std::endl;
            break;

        case ECorrection::FRACTION:
            std::cerr << "> CORRECTION ······· = FRACTION (special)" << std::endl;
            break;

        default:
            throw std::logic_error( "Correction method not specified." );
    }

    switch (this->method)
    {
        case EMethod::COINCIDENCE:
            std::cerr << "> METHOD ··········· = COINCIDENCE" << std::endl;
            break;

        case EMethod::CONNECTIVITY:
            std::cerr << "> METHOD ··········· = CONNECTIVITY" << std::endl;
            break;

        default:
            throw std::logic_error( "Analysis method not specified." );
    }

    switch (this->alt_hypothesis)
    {
        case EHypothesis::GREATER:
            std::cerr << "> ALT_HYPOTHESIS ··· = GREATER" << std::endl;
            break;

        case EHypothesis::LESS:
            std::cerr << "> ALT_HYPOTHESIS ··· = LESS" << std::endl;
            break;

        case EHypothesis::TWOTAILED:
            std::cerr << "> ALT_HYPOTHESIS ··· = TWOTAILED" << std::endl;
            break;

        default:
            throw std::logic_error( "Significance method not specified." );
    }

    switch (this->coin_max_mode)
    {
        case EMaxMode::ACCOMPANY:
            if (this->method != EMethod::COINCIDENCE)
            {
                throw std::logic_error( "Max-mode specified as ACCOMPANY when method is not COINCIDENCE." );
            }
            std::cerr << "> MAX_MODE ········· = ACCOMPANY" << std::endl;
            break;

        case EMaxMode::AVOID:
            if (this->method != EMethod::COINCIDENCE)
            {
                throw std::logic_error( "Max-mode specified as AVOID when method is not COINCIDENCE." );
            }

            std::cerr << "> MAX_MODE ········· = ACCOMPANY" << std::endl;
            break;

        default:
            if (this->method == EMethod::COINCIDENCE)
            {
                throw std::logic_error( "Max-mode not specified when method is COINCIDENCE." );
            }
            break;
    }

    switch (this->coin_set_mode)
    {
        case ESetMode::FULL:
            if (this->method != EMethod::COINCIDENCE)
            {
                throw std::logic_error( "Set-mode specified as FULL when method is not COINCIDENCE." );
            }

            std::cerr << "> SET_MODE ········· = FULL" << std::endl;
            break;

        case ESetMode::INTERSECTION:
            if (this->method != EMethod::COINCIDENCE)
            {
                throw std::logic_error( "Set-mode specified as INTERSECTION when method is not COINCIDENCE." );
            }

            std::cerr << "> SET_MODE ········· = INTERSECTION" << std::endl;
            break;

        default:
            if (this->method == EMethod::COINCIDENCE)
            {
                throw std::logic_error( "Set-mode not specified when method is COINCIDENCE." );
            }
            break;
    }

    if (this->permit_filter)
    {
        std::cerr << "> PERMIT_FILTER ···· = YES" << std::endl;
    }
    else
    {
        std::cerr << "> PERMIT_FILTER ···· = NO" << std::endl;
    }

    if (this->verbose)
    {
        std::cerr << "> VERBOSE ·········· = YES" << std::endl;
    }
    else
    {
        std::cerr << "> VERBOSE ·········· = NO" << std::endl;
    }

    if (this->output_all)
    {
        std::cerr << "> OUTPUT_ALL ······· = YES" << std::endl;
    }
    else
    {
        std::cerr << "> OUTPUT_ALL ······· = NO" << std::endl;
    }

    if (this->sig_level < 0)
    {
        throw std::logic_error( "Significance level is invalid." );
    }

    switch (this->correction)
    {
        case ECorrection::FRACTION:
            if (this->method == EMethod::COINCIDENCE)
            {
                throw std::logic_error( "Cutoff mode specified when method != CONNECTIVITY." );
            }

            std::cerr << "> FRACTION CUTOFF ·· = " << this->sig_level << std::endl;
            std::cerr << "> SIGNIFICANCE_LEVEL = N/A" << std::endl;
            break;

        default:
            std::cerr << "> FRACTION CUTOFF ·· = N/A" << std::endl;
            std::cerr << "> SIGNIFICANCE_LEVEL = " << this->sig_level << std::endl;
            break;
    }
    //if (!this->phylogeny.empty()) {
	//std::cerr << "> PHYLOGENY ........ = " << this->phylogeny << std::endl;
    //} else {
	//throw std::logic_error( "Phylogeny file is missing." ) ;
    //}
    
    if (this->combined_file_name.empty())
    {
        if (this->alpha_file_name.empty())
        {
            throw std::logic_error( "Alpha file is missing." );
        }
        
        if (this->beta_file_name.empty())
        {
            throw std::logic_error( "Beta file is missing." );
        }
        
        std::cerr << "> ALPHA_FILE ······· = " << this->alpha_file_name << std::endl;
        std::cerr << "> BETA_FILE ········ = " << this->beta_file_name << std::endl;
    }
    else
    {
        if (this->method == EMethod::CONNECTIVITY)
        {
            throw std::logic_error( "Cannot use a combined file with connectivity analysis." );
        }
        
        if (!this->alpha_file_name.empty())
        {
            throw std::logic_error( "Cannot specify both the combined and alpha filenames." );
        }
        
        if (!this->beta_file_name.empty())
        {
            throw std::logic_error( "Cannot specify both the combined and alpha filenames." );
        }
        
        std::cerr << "> COMBINED_FILE ···· = " << this->combined_file_name << std::endl;
    }

    
    std::cerr << "> ALPHA_NAME ······· = " << this->alpha_name << std::endl;
    std::cerr << "> BETA_NAME ········ = " << this->beta_name << std::endl;
    
    if (this->combined_file_name.empty())
    {
        std::cerr << "> GAMMA_NAME ······· = " << this->gamma_name << std::endl;
    }
    
}


#pragma clang diagnostic pop
