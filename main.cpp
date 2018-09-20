#include <iostream>  
#include "dataset.h"
#include "coincidence.h"
#include "lineage.h"
#include "gexf.h"
#include "network.h"
#include "constants.h"
#include "parameters.h"
#include "connectivity.h"
#include "test_cases.h"
#include <stdio.h>
#include <cstdlib>
#include <cstdio>
#include <array>

#include "bugfix.h"
/**
 * Entry point, boom!
 * 
 * @param argc  System 
 * @param argv  System
 * @return      System 
 */
int main( int argc, const char** argv )
{

	//Bug test
	//BugFix::test1();
	//std::cerr << "Pass: returned from test1. Exiting now..." << std::endl;
	//return 0;


    const char* version = "1.1";
    
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
    // Check for ask for help
    //
    std::string help = "-h";
    for (int n = 1; n < argc; ++n) {
        std::string arg = argv[ n ];
	if (arg == help) {
		std::cerr << "./confinder [OPTIONS]" << std::endl;
		std::cerr << "File input- specify either: " << std::endl;
		std::cerr << "    -a or --a              The path to the Gamma-to-Alpha family file with (gamma)(TAB)(alpha)" << std::endl;
		std::cerr << "    -b or --b              The path to the Gamma-to-Beta family file with (gamma)(TAB)(beta)" << std::endl;

		std::cerr << "  or" << std::endl;
		std::cerr << "    -d or --d              The path of the Alpha-to-Beta file with (alpha)(TAB)(beta)" << std::endl;
		std::cerr << "Significance- specify: " << std::endl;
		std::cerr << "    -L or --level          Specify the significnace level cutoff" << std::endl;
		std::cerr << "Naming (optional): " << std::endl;
		std::cerr << "    -A or --alpha          Specify the name of the Alpha group" << std::endl;
		std::cerr << "    -B or --beta           Specify the name of the Beta group" << std::endl;
		std::cerr << "    -C or --gamma          Specify the name of the Gamma group" << std::endl;
		std::cerr << "Mode- specify: " << std::endl;
		std::cerr << "    -e or --connectivity   Connectivity analysis" << std::endl;
		std::cerr << "    -c or --coincidence    Coincidence analysis" << std::endl;
		std::cerr << "Significance correction- specify: " << std::endl;
		std::cerr << "    -n or --nocorrection   No correction, use value as-is" << std::endl;
		std::cerr << "    -m or --bonferroni     Bonferroni correction multiple correction" << std::endl;
		std::cerr << "    -r or --fraction       (Connectivity analysis only) Use fraction rather than p-value" << std::endl;
		std::cerr << "Alternative hypothesis- specify: " << std::endl;
		std::cerr << "    -g or --greater        Greater (recommended)" << std::endl;
		std::cerr << "    -l or --less           Less" << std::endl;
		std::cerr << "    -t or --twotailed      Two-tailed" << std::endl;
		std::cerr << "Set mode (mandatory for coincidence analysis):" << std::endl;
		std::cerr << "    -f or --full           Full; consider all data." << std::endl;
		std::cerr << "    -u or --union          Union; only consider the edges which pertain to the groups being analysed." << std::endl;
		std::cerr << "Max mode (mandatory for coincidence analysis):" << std::endl;
		std::cerr << "    -o or --accompany      Overlap; identify groups that tend to coincide." << std::endl;
		std::cerr << "    -s or --avoid          Separation; identify groups that tend to avoid." << std::endl;
		std::cerr << "Miscellaneous:" << std::endl;
		std::cerr << "    -p or --phylogeny      Phylogeny of Betas in Newick format (mandatory)" << std::endl;
		std::cerr << "    -v or --verbose        Verbose output." << std::endl;
		std::cerr << "    -i or --filter         Permit filtering of saturated and low-abundance data." << std::endl;
		std::cerr << "    -U or --upfilthreshold Upper filter threshold for high-abundance data filtering (default: 1.0 i.e. any alpha in >=100/% of betas." << std::endl;
		std::cerr << "    -F or --filthreshold   Threshold for low-abundance data filtering (default: 0.05 i.e. any alpha in <=5\% of betas." << std::endl;
		std::cerr << "    -q or --query          Query a specific Alpha family." << std::endl;
		std::cerr << "    -T or --test           Runs the test cases and exits." << std::endl;
		std::cerr << "    -E or --all            Outputs all results, regardless of significance." << std::endl;
		return 1;
	}
    }

    //
    // Read arguments
    //
    std::cerr << "Reading arguments..." << std::endl;
    TParameters options = TParameters::parse( argc, argv );
    options.print_and_assert();

    //
    // Save coinfinder location to pass to Python
    //
    std::string command("which coinfinder");
    std::array<char, 128> buffer;
    std::string result;
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
	std::cerr << "Couldn't start command." << std::endl;
	return 1;
    }
    while (fgets(buffer.data(), 128, pipe) != NULL) {
 	result += buffer.data();
    }
    std::cerr << "Result of which coinfinder: " << result << std::endl;
    //Strip off the trailing "/coinfinder"
    result = result.substr(0,result.length()-11);
    //Deal with coinfinder not being set
    if (result == "") {
	result = ".";
    }

    //
    // Load in relations
    //
    DataSet dataset = DataSet( options );
    dataset.read_files( options.alpha_file_name, options.beta_file_name, options.combined_file_name, options.phylogeny, options.filt_thres, options.upper_filt_thres );

    //
    // Do what it is that needs to be done
    //
    switch (options.method)
    {
        case EMethod::CONNECTIVITY:
            Connectivity::run( dataset );
            break;

        case EMethod::COINCIDENCE:
	{
            Coincidence::run( dataset, options.phylogeny, result );
	    Lineage::run( dataset );
	    Gexf::run( dataset );
	    Network::run( dataset );
            break;
	}

        default:
            throw std::logic_error( "Invalid `method` specified." );
    }



    //
    // Clean up
    //
    std::cerr << "All done!" << std::endl;
    return 0;
}
