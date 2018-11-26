#include <iostream>
#include <fstream>
#include <vector>

#include "lineage.h"
#include "elements.h"
#include "parameters.h"
/**
 * Calculates the level of lineage-dependence per gene in the dataset
 */

int Lineage::run( DataSet&  dataset, const std::string& source_path, const std::string& call_path, const std::string& phylogeny, const std::string& gene_pa, int num_cores, bool Rmsgs, const std::string& prefix)
{
	//Identify alphas with coincident edges & output to file
	//Write an edges file while I'm at it for use in Network
	std::ofstream nodefile;
	std::string nodename = prefix + "_nodes_in.csv";
	std::ofstream edgefile;
	std::string edgename = prefix + "_edges.csv";
        nodefile.open(nodename);
	edgefile.open(edgename);
	//Write header to edge table for input into gephi
	edgefile << "Source,Target,weight" << std::endl;
	const id_lookup<Alpha>& alpha_table = dataset.get_alphas();
	for (const auto& alpha_list : alpha_table.get_table()) {
        	Alpha& alpha = *alpha_list.second;
		if (alpha.get_num_coincident_edges() > 0) {
			nodefile << alpha.get_name() << ",";
			const std::map<const Alpha*, double>& coincident_edges = alpha.get_coincident_edges();
			for (const auto& edge_list : coincident_edges) {
				edgefile << alpha.get_name() << "," << (edge_list.first)->get_name() << "," << edge_list.second << std::endl;
			}
		}
	}
	nodefile << std::endl;
	nodefile.close();
	edgefile.close();

	//Call R to calculate D
	std::cerr << "Calculating lineage dependence..." << std::endl;
	std::string syscall = "Rscript " + source_path + "/lineage.R -p " + call_path + " -t " + phylogeny + " -g " + gene_pa + " -c " + std::to_string(num_cores) + " -o " + prefix;
	if (Rmsgs) {
		system(syscall.c_str());
	} else {
		std::string out = systemSTDOUT(syscall);
		if (out.find("Phylogeny contains pairs of tips on zero branch lengths, cannot currently simulate") != std::string::npos) {
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
	
	//Save D to alpha as an attribute
	std::ifstream file_in;
	std::string filename = prefix + "_nodes.csv";
    	file_in.open(filename);
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
				if ((cell != "NA") & (cell != "Result")) {
					D = std::stod(cell);
					Dvalues[ name ] = D;
				}
			}
			left = !left;
		}
		file_in.close();
	} else {
		std::cerr << "Error: unable to open " << prefix << "_nodes.csv" << std::endl;
	}
	int count = 0;
	for (const auto& alpha_list : alpha_table.get_table()) {
                Alpha& alpha = *alpha_list.second;
                if (Dvalues.count(alpha.get_name()) > 0) {
                        alpha.register_D(Dvalues[alpha.get_name()]);
			count++;
                }
        }
	return(0);
}

/**
 * Execute a command and get the result.
 *
 * @param   cmd - The system command to run.
 * @return  The string command line output of the command.
 */
std::string Lineage::systemSTDOUT(std::string cmd) {

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
