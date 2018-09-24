#include <iostream>
#include <fstream>
#include <vector>

#include "network.h"
#include "elements.h"
#include "parameters.h"
/**
 * Prunes the dataset based on D-value and draws a network and heatmap of the pruned dataset.
 */

void Network::run( DataSet& dataset )
{
	//Call R to draw heatmap and network
	std::cerr << "Call network.R..." << std::endl;
	system("Rscript network.R");
	std::cerr << "Return from network.R..." << std::endl;
}
