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

#include "boost/filesystem.hpp"
#include "boost/algorithm/string/replace.hpp"

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
		//std::cerr << "    -a or --a              The path to the Gamma-to-Alpha family file with (gamma)(TAB)(alpha)" << std::endl;
		//std::cerr << "    -b or --b              The path to the Gamma-to-Beta family file with (gamma)(TAB)(beta)" << std::endl;

		//std::cerr << "  or" << std::endl;
		//std::cerr << "    -d or --d              The path of the Alpha-to-Beta file with (alpha)(TAB)(beta)" << std::endl;
		std::cerr << "    -i or --input          The path to the gene_presence_absence.csv output from Roary" << std::endl;
		std::cerr << "                           -or-" << std::endl;
		std::cerr << "                           The path of the Alpha-to-Beta file with (alpha)(TAB)(beta)" << std::endl;
		std::cerr << "    -I or --inputroary     Set if -i is in the gene_presence_absence.csv format from Roary" << std::endl;
                std::cerr << "    -p or --phylogeny      Phylogeny of Betas in Newick format (required)" << std::endl;
                std::cerr << "Max mode (mandatory for coincidence analysis):" << std::endl;
                std::cerr << "    -a or --associate      Overlap; identify groups that tend to associate/co-occur." << std::endl;
                std::cerr << "    -d or --dissociate     Separation; identify groups that tend to dissociate/avoid." << std::endl;
		std::cerr << "Significance- specify: " << std::endl;
		std::cerr << "    -L or --level          Specify the significnace level cutoff (default: 0.05)" << std::endl;
		//std::cerr << "Mode- specify: " << std::endl;
		//std::cerr << "    -e or --connectivity   Connectivity analysis" << std::endl;
		//std::cerr << "    -c or --coincidence    Coincidence analysis" << std::endl;
		std::cerr << "Significance correction- specify: " << std::endl;
                std::cerr << "    -m or --bonferroni     Bonferroni correction multiple correction (recommeneded)" << std::endl;
		std::cerr << "    -n or --nocorrection   No correction, use value as-is" << std::endl;
		std::cerr << "    -c or --fraction       (Connectivity analysis only) Use fraction rather than p-value" << std::endl;
		std::cerr << "Alternative hypothesis- specify: " << std::endl;
		std::cerr << "    -g or --greater        Greater (recommended)" << std::endl;
		std::cerr << "    -l or --less           Less" << std::endl;
		std::cerr << "    -t or --twotailed      Two-tailed" << std::endl;
		//std::cerr << "Set mode (mandatory for coincidence analysis):" << std::endl;
		//std::cerr << "    -f or --full           Full; consider all data." << std::endl;
		//std::cerr << "    -u or --union          Union; only consider the edges which pertain to the groups being analysed." << std::endl;
                //std::cerr << "Naming (optional): " << std::endl;
                //std::cerr << "    -A or --alpha          Specify the name of the Alpha group" << std::endl;
                //std::cerr << "    -B or --beta           Specify the name of the Beta group" << std::endl;
                //std::cerr << "    -C or --gamma          Specify the name of the Gamma group" << std::endl;
		std::cerr << "Miscellaneous:" << std::endl;
		std::cerr << "    -x or --num_cores      The number of cores to use (default: 2)" << std::endl;
		std::cerr << "    -v or --verbose        Verbose output." << std::endl;
		std::cerr << "    -r or --filter         Permit filtering of saturated and low-abundance data." << std::endl;
		std::cerr << "    -U or --upfilthreshold Upper filter threshold for high-abundance data filtering (default: 1.0 i.e. any alpha in >=100/% of betas." << std::endl;
		std::cerr << "    -F or --filthreshold   Threshold for low-abundance data filtering (default: 0.05 i.e. any alpha in <=5\% of betas." << std::endl;
		std::cerr << "    -q or --query          Query a specific gene." << std::endl;
		std::cerr << "    -T or --test           Runs the test cases and exits." << std::endl;
		std::cerr << "    -E or --all            Outputs all results, regardless of significance." << std::endl;
		std::cerr << "Output:" << std::endl;
		std::cerr << "    -o or --output         The prefix of all output files (default: coincident)." << std::endl;
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
    std::string source_path;
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe) {
	std::cerr << "Couldn't start command." << std::endl;
	return 1;
    }
    while (fgets(buffer.data(), 128, pipe) != NULL) {
 	source_path += buffer.data();
    }
    //Strip off the trailing "/coinfinder"
    source_path = source_path.substr(0,source_path.length()-11);
    //Deal with coinfinder not being set
    if (source_path == "") {
	source_path = ".";
    }

    //
    // Save directory of coinfinder call to pass to R
    //
    boost::filesystem::path ll_path = boost::filesystem::initial_path();
    std::string call_path = ll_path.string();
    //Replace any possible spaces with backslashes
    boost::replace_all(call_path, " ", "\\ ");

    // 
    // Check for existence of input files
    //
    if ( !boost::filesystem::exists(options.combined_file_name) ) {
	std::cerr << "Input file " << options.combined_file_name << " does not exist." << std::endl;
	std::cerr << "Exiting..." << std::endl;
	return(-1);
    }

    //
    // If no phylogeny is provided, output a warning to the user
    //
    if (options.phylogeny.empty()) {
	std::cerr << std::endl;
	std::cerr << "You have not provided a phylogeny. This is okay, coinfinder works just fine without a phylogeny. The steps involving the detection of lineage dependent elements will be skipped." << std::endl;
	std::cerr << std::endl;
	std::cerr << "**Please be aware that coinfinder will generate a mock phylogeny for the purposes of data presentation in the output heatmaps.**" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Press any key to accept and continue or Ctrl+C to quit." << std::endl;
	getchar();
    }
	//std::cerr << "Input phylogeny " << options.phylogeny << " does not exist." << std::endl;
	//std::cerr << "Exiting..." << std::endl;
	//return(-1);
    //}

    //
    // If phylogeny is provided, ensure that there are no zero length branches
    //
    if (!options.phylogeny.empty()) {
	std::string syscall = "Rscript " + source_path + "/coinfind-code/check_zeroes.R -a " + call_path + " -t " + options.phylogeny;
        if (options.Rmsgs) {
        	system(syscall.c_str());
        }
        std::string out = Lineage::systemSTDOUT(syscall);
        if (out.find("TRUE") != std::string::npos) {
               	std::cerr << "Error: Input phylogeny contains branch lengths which equal zero." << std::endl;
                std::cerr << "The R function caper::phylo.d, which coinfinder uses to define lineage dependence, cannot handle zero length branches." << std::endl;
                std::cerr << "Please adjust your phylogeny before continuing." << std::endl;
                std::cerr << "Exiting..." << std::endl;
                return(-1);
        }
        if ((out.find("Error") != std::string::npos) || (out.find("error") != std::string::npos)) {
               	std::cerr << "ERROR MESSAGE FROM R: " << std::endl;
                if (out.find("Error") != std::string::npos) {
                       	std::cerr << out.substr(out.find("Error")) << std::endl;
                }
                if (out.find("error") != std::string::npos) {
                       	std::cerr << out.substr(out.find("error")) << std::endl;
                }
                return(-1);
        }
        if ((out.find("Killed") != std::string::npos)) {
               	std::cerr << "ERROR MESSAGE FROM R: " << std::endl;
                std::cerr << out.substr(out.find("Killed")) << std::endl;
                return(-1);
        }
    }
    
    //
    // Format gene_p_a to classic coinfinder input file, if necessary
    //
    if (options.roary) {
    	std::cerr << "Formating Roary output for input into coinfinder..." << std::endl;
    	std::string sysStr = "python3 " + source_path + "/coinfind-code/format_roary.py -i " + options.combined_file_name;
    	std::string out = Lineage::systemSTDOUT(sysStr);
    	if ((out.find("Error") != std::string::npos) || (out.find("error") != std::string::npos)) {
    		std::cerr << "ERROR MESSAGE FROM Python: " << std::endl;
		std::cerr << out << std::endl;
        	return(-1);
	}
    }

    //
    // Create gene_p_a for later input into R, if necessary
    //
    if (!options.roary) {
	std::cerr << "Formating input into gene_p_a for input into coinfinder..." << std::endl;
	std::string sysStr = "python3 " + source_path + "/coinfind-code/create_roary.py -i " + options.combined_file_name;
	std::string out = Lineage::systemSTDOUT(sysStr);
	if ((out.find("Error") != std::string::npos) || (out.find("error") != std::string::npos)) {
		std::cerr << "ERROR MESSAGE FROM Python: " << std::endl;
		std::cerr << out << std::endl;
		return(-1);
	}
    }

    //
    // Load in relations
    //
    DataSet dataset = DataSet( options );
    //dataset.read_files( options.alpha_file_name, options.beta_file_name, options.combined_file_name, options.phylogeny, options.filt_thres, options.upper_filt_thres );
    if (options.roary) {
    	dataset.read_files( options.alpha_file_name, options.beta_file_name, "coincident-input-edges.csv", options.phylogeny, options.filt_thres, options.upper_filt_thres );
    } else {
	dataset.read_files( options.alpha_file_name, options.beta_file_name, options.combined_file_name, options.phylogeny, options.filt_thres, options.upper_filt_thres );
    }

    //
    // Check to be sure neither alphas or betas are empty before beginning work.
    //
    if (dataset.get_num_alphas() <= 0 || dataset.get_num_betas() <=0) {
	std::cerr << "Not enough alphas or betas after file load in and dataset trimming to proceed." << std::endl;
	std::cerr << "Exiting..." << std::endl;
	return(-1);
    }

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
	    int retval = 0;
	    std::string gene_pa = "";
	    if (options.roary) {
		   gene_pa = options.combined_file_name;
	    } else {
		   gene_pa = "gene_presence_absence.csv";
	    } 
	    retval = Coincidence::run( dataset, options.phylogeny, source_path, options.num_cores, options.prefix );
	    if (retval != 0) {
		std::cerr << "Coinfinder did not find any significant coinciding pairs in the input." << std::endl;
		std::cerr << "Exiting..." << std::endl;
		return(-1);
	    }
	    retval = Lineage::run( dataset, source_path, call_path, options.phylogeny, gene_pa, options.num_cores, options.Rmsgs, options.prefix );
	    if(retval != 0) {
		return(-1);
	    }
	    retval = Network::run( dataset, source_path, call_path, options.phylogeny, gene_pa, options.Rmsgs, options.prefix );
	    if(retval != 0) {
		return(-1);
	    }
	    Gexf::run( dataset, options.prefix );
	    //Remove intermediate files
	    //This file is only created if Roary input is used; don't error if it can't be removed.
	    //if (boost::filesystem::exists("coincident-input-edges.csv")) {
	    //	const int result = remove( "coincident-input-edges.csv" );
	    	//if (result != 0) {
    	    	//	    printf( "%s\n", strerror( errno ) ); // No such file or directory
	    	//}
	    //}
	    //const char* nodename = (options.prefix + "_nodes_in.csv").c_str();
	    //if (boost::filesystem::exists( std::string(nodename) )) {
	    //const int result2 = remove( nodename );
	    //if (result2 != 0) {
	    //	    printf( "%s\n", strerror( errno ) ); // No such file or directory
	    //}
            break;
	    //}
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
