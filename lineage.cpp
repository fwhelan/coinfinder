#include <iostream>
#include <fstream>
#include <vector>

#include "lineage.h"
#include "elements.h"
#include "parameters.h"
/**
 * Calculates the level of lineage-dependence per gene in the dataset
 */

void Lineage::run( DataSet&  dataset)
{
	std::cerr << "hello, world --lineage." << std::endl;
	std::cerr << "dataset has " << dataset.get_num_coincident_edges() << " coincident edges" << std::endl;

	//Identify alphas with coincident edges & output to file
	std::ofstream nodefile;
        nodefile.open("coincident_nodes_in.csv");
	const id_lookup<Alpha>& alpha_table = dataset.get_alphas();
	for (const auto& alpha_list : alpha_table.get_table()) {
        	Alpha& alpha = *alpha_list.second;
		if (alpha.get_num_coincident_edges() > 0) {
			nodefile << alpha.get_name() << ",";
			std::cerr << alpha.get_name() << ",";
		}
	}
	nodefile << std::endl;
	nodefile.close();
	std::cerr << std::endl;

	//Call R to calculate D
	std::cerr << "Call lineage.R..." << std::endl;
	//system("Rscript lineage.R");
	std::cerr << "Return from lineage.R..." << std::endl;
	
	//Save D to alpha as an attribute
	std::ifstream file_in;
    	file_in.open( "coincident_nodes.csv" );
	std::string cell;
	bool left = true;
	std::string name;
	double D;
	std::map<std::string, double> Dvalues;
	if (file_in.is_open()) {
		while (getline( file_in, cell, left ? static_cast<char>('\t') : static_cast<char>('\n'))) {   
        		if (left) {
       				name = cell;
			} else {
				if (cell != "NA") {
					D = std::stod(cell);
					Dvalues[ name ] = D;
				}
			}
			left = !left;
		}
		file_in.close();
	} else {
		std::cerr << "Error: unable to open coincident_nodes.csv" << std::endl;
	}
	int count = 0;
	for (const auto& alpha_list : alpha_table.get_table()) {
                Alpha& alpha = *alpha_list.second;
                if (Dvalues.count(alpha.get_name()) > 0) {
                        alpha.register_D(Dvalues[alpha.get_name()]);
			count++;
			std::cerr << alpha.get_name() << ",";
                }
        }
	std::cerr << std::endl;
	std::cerr << "There are " << count << " alphas with coincident edges." << std::endl;
}
