#include <iostream>
#include <fstream>
#include <vector>

#include "network.h"
#include "elements.h"
#include "parameters.h"
/**
 * Prunes the dataset based on D-value and draws a network and heatmap of the pruned dataset.
 */

int Network::run( DataSet& dataset )
{
	//Call R to draw heatmap and network
	std::cerr << "Generate network and coincidence heatmap..." << std::endl;
	std::string out = systemSTDOUT("Rscript network.R");
	if (out.find("Error") != std::string::npos) {
                std::cerr << "ERROR MESSAGE FROM R: " << std::endl;
                std::cerr << out.substr(out.find("Error")) << std::endl;
		return(-1);
        }
	return(0);
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
