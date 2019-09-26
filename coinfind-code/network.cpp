#include <iostream>
#include <fstream>
#include <vector>

#include "network.h"
#include "elements.h"
#include "parameters.h"
/**
 * Prunes the dataset based on D-value and draws a network and heatmap of the pruned dataset.
 */

int Network::run( DataSet& dataset, const std::string& source_path, const std::string& call_path, const std::string& phylogeny, const std::string& gene_pa, bool Rmsgs, const std::string& prefix )
{
	//Call R to draw heatmap and network iff phylogeny provided
	if(!phylogeny.empty()) {
		//std::cerr << "Generate network and coincidence heatmap..." << std::endl;
		std::cerr << "Generate coincidence heatmap..." << std::endl;
		std::string syscall = "Rscript " + source_path + "/coinfind-code/network.R -p " + call_path + " -t " + phylogeny + " -g " + gene_pa + " -o " + prefix; 
		if (Rmsgs) {
			system(syscall.c_str());
		} else {
			std::string out = systemSTDOUT(syscall);
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
		}
		return(0);
	//Call R to draw a heatmap and network without a phylogeny
	} else {
		std::cerr << "Generate coincidence heatmap (without phylogeny input)..." << std::endl;
		std::string syscall = "Rscript " + source_path + "/network_nophylogeny.R -p " + call_path + " -g " + gene_pa + " -o " + prefix;
		if (Rmsgs) {
                        system(syscall.c_str());
                } else {
                        std::string out = systemSTDOUT(syscall);
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
                }
                return(0);
	}
}

/**
 * Execute a command and get the result.
 *
 * @param   cmd - The system command to run.
 * @return  The string command line output of the command.
 */
std::string Network::systemSTDOUT(std::string cmd) {

    std::string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    cmd.append(" 2>&1"); // Do we want STDERR?

    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    return data;
}
