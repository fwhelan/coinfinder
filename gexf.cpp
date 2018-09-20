#include <iostream>
#include <fstream>
#include <vector>

#include "gexf.h"
#include "elements.h"
#include "parameters.h"
/**
 * Calculates the level of lineage-dependence per gene in the output gene_list
 */

void Gexf::run( DataSet& dataset )
{
	std::cerr << "hello, world --gexf." << std::endl;
	std::ofstream gexf;
        gexf.open("coincident_network.gexf");
	std::string node_attr_xml = "";
	std::string node_xml = "";
	std::string edge_attr_xml = "";
	std::string edge_xml = "";
	//Cycle through nodes and edges to output to gexf format
	std::string alpha1_name;
	double alpha1_D;
	int alpha1_col;
	double alpha1_size;
	int edge_counter = 0;
	const id_lookup<Alpha>& alpha_table = dataset.get_alphas();
        for (const auto& alpha_list : alpha_table.get_table()) {
                Alpha& alpha = *alpha_list.second;
		if (alpha.get_num_coincident_edges() > 0) {
			alpha1_name = alpha.get_name();
			alpha1_D = alpha.get_D();
			alpha1_col = 289*alpha1_D;
			if (alpha1_col > 255) { alpha1_col = 255; }
			if (alpha1_col < 0)   { alpha1_col = 0; } //TODO- make this more comprehensive
			alpha1_size = 10+(alpha1_D*2);
			//push node to node array for gexf
			//node_xml += "<node id=\"" + alpha1_name + "\" label=\"" + alpha1_name + "\"/>";
			node_attr_xml += "<node id=\"" + alpha1_name + "\" label=\"" + alpha1_name + "\">\n" +
						" <attvalues>\n" +
						"  <attvalue for=\"D-value\" value=\"" + std::to_string(alpha1_D) + "\"/>\n" +
						" </attvalues>\n" +
						"<viz:color r=\"" + std::to_string(alpha1_col) + "\" g=\"173\" b=\"66\"/>" +
						"<viz:size value=\"" + std::to_string(alpha1_size) + "\" />" + 
						"</node>\n";
			std::string alpha2_name;
			double p_value;
			double edge_weight;
			const std::map<const Alpha*, double>& edges = alpha.get_coincident_edges();
			for(const auto& edge_list : edges) {
				alpha2_name = (edge_list.first)->get_name();
				p_value = edge_list.second;
				edge_weight = ((1-p_value)*2)+1; //TODO- make this more comprehensive
				//push edge to edge array for gexf
				edge_attr_xml += "<edge id=\"" + std::to_string(edge_counter) + "\" label=\"" + std::to_string(p_value)
					+ "\" source=\"" + alpha1_name + "\" target=\"" + alpha2_name + "\" weight=\"" + std::to_string(edge_weight) + "\">" +
						" <attvalues>\n" +
						"  <attvalue for=\"p-value\" value=\"" + std::to_string(p_value) + "\"/>\n" +
						" </attvalues>\n" +
						"<viz:color r=\"146\" g=\"142\" b=\"142\"/>" + 
						"</edge>\n";
				edge_counter++;
								
			}
		}
	}
	//Output to gexf file
	gexf << "<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\" xmlns:viz=\"http://www.gexf.net/1.1draft/viz\">" << std::endl;
        gexf << "<graph mode=\"static\" defaultedgetype=\"undirected\">" << std::endl;
	gexf << "<attributes class=\"node\">" << std::endl;
        gexf << "<attribute id=\"D-value\" title=\"D-value\" type=\"double\"/>" << std::endl;
        gexf << "</attributes>" << std::endl;
	gexf << "<attributes class=\"edge\">" << std::endl;
	gexf << "<attribute id=\"p-value\" title=\"p-value\" type=\"double\"/>" << std::endl;
	gexf << "</attributes>" << std::endl;
        gexf << "<nodes>" << std::endl;
    	gexf << node_attr_xml;
        gexf << "</nodes>" << std::endl;
        gexf << "<edges>" << std::endl;
    	gexf << edge_attr_xml << std::endl;
        gexf << "</edges>" << std::endl;
        gexf << "</graph>" << std::endl;
        gexf << "</gexf>" << std::endl;

	gexf.close();
}
