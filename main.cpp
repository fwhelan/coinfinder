#include <iostream>  
#include "dataset.h"
#include "coincidence.h"
#include "constants.h"
#include "parameters.h"
#include "connectivity.h"
#include "test_cases.h"


/**
 * Entry point, boom!
 * 
 * @param argc  System 
 * @param argv  System
 * @return      System 
 */
int main( int argc, const char** argv )
{
    const char* version = "1.0.0.7.";
    
    //
    // Validate arguments
    //
    if (argc == 1)
    {
        std::cerr << "CoinFinder " << version << " " << std::endl;
#ifndef NDEBUG
        std::cerr << "DEBUG BUILD @ ";
#else
        std::cerr << "RELEASE BUILD @ ";
#endif
        std::cerr << __DATE__ << " " << __TIME__ << "." << std::endl;
        std::cerr << "Please see readme.md for details." << std::endl;

        return 1;
    }

    //
    // Read arguments
    //
    std::cerr << "Reading arguments..." << std::endl;
    TParameters options = TParameters::parse( argc, argv );
    options.print_and_assert();

    //
    // Load in relations
    //
    DataSet dataset = DataSet( options );
    dataset.read_files( options.alpha_file_name, options.beta_file_name, options.combined_file_name, options.phylogeny );

    //
    // Do what it is that needs to be done
    //
    switch (options.method)
    {
        case EMethod::CONNECTIVITY:
            Connectivity::run( dataset );
            break;

        case EMethod::COINCIDENCE:
            Coincidence::run( dataset, options.phylogeny );
            break;

        default:
            throw std::logic_error( "Invalid `method` specified." );
    }



    //
    // Clean up
    //
    std::cerr << "All done!" << std::endl;
    return 0;
}


