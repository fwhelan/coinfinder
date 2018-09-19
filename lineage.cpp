#include <iostream>
#include <fstream>
#include <vector>

#include "lineage.h"
#include "elements.h"
#include "parameters.h"
/**
 * Calculates the level of lineage-dependence per gene in the output gene_list
 */

void Lineage::run( std::map<std::string, int> gene_list )
{
	std::cerr << "hello, world." << std::endl;
	std::cerr << "gene_list has " << gene_list.size() << " elements." << std::endl;

	std::ofstream nodefile;
	nodefile.open("coincident_nodes_in.csv");

	for(auto& x : gene_list) {
    		nodefile << x.first << ",";
	}
	nodefile << std::endl;
	nodefile.close();

	std::cerr << "Call lineage.R..." << std::endl;
	//system("Rscript lineage.R");
	std::cerr << "Return from lineage.R..." << std::endl;
	
	std::ifstream file_in;
    	file_in.open( "coincident_nodes.csv" );
	std::string line;
	if (file_in.is_open()) {
		while( getline (file_in, line) ) {
			//TODO save D to alpha as an attribute
			std::cerr << line << std::endl;
		}
		file_in.close();
	} else {
		std::cerr << "Error: unable to open coincident_nodes.csv" << std::endl;
	}
}
