#include <iostream>
#include <fstream>
#include <vector>

#include "gexf.h"
#include "elements.h"
#include "parameters.h"
/**
 * Calculates the level of lineage-dependence per gene in the output gene_list
 */

void Gexf::run( DataSet& dataset, const std::string& prefix )
{
	std::ofstream gexf;
	std::string gexfname = prefix + "_network.gexf";
        gexf.open(gexfname);
	std::string node_attr_xml = "";
	std::string edge_attr_xml = "";
	//Cycle through nodes and edges to output to gexf format
	std::string alpha1_name;
	double alpha1_D;
	std::string alpha1_col; //int
	double alpha1_size;
	int edge_counter = 0;
	const id_lookup<Alpha>& alpha_table = dataset.get_alphas();
        for (const auto& alpha_list : alpha_table.get_table()) {
                Alpha& alpha = *alpha_list.second;
		if (alpha.get_num_coincident_edges() > 0) {
			alpha1_name = alpha.get_name();
			alpha1_D = alpha.get_D();
			//alpha1_col = 255*alpha1_D; //Most meaningful values of D are between 0 and 1
			//if (alpha1_col > 255) { alpha1_col = 255; }
			//if (alpha1_col < 0)   { alpha1_col = 0; }
			alpha1_size = 20+(alpha1_D*2);
			//Colour node by component number
			std::string compname = prefix + "_components.csv";
        		std::string syscall = "grep \"" + alpha1_name + "\" " + compname + " | cut -f 1";
                	std::string ret = systemSTDOUT(syscall);
			alpha1_col = componentLookup(stoi(ret));
			//push node to node array for gexf "<viz:color r=\"" + std::to_string(alpha1_col) + "\" g=\"173\" b=\"66\"/>" +
			node_attr_xml += "<node id=\"" + alpha1_name + "\" label=\"" + alpha1_name + "\">\n" +
						" <attvalues>\n" +
						"  <attvalue for=\"D-value\" value=\"" + std::to_string(alpha1_D) + "\"/>\n" +
						" </attvalues>\n" +
						"<viz:color hex=\"" + alpha1_col + "\"/>" +
						"<viz:size value=\"" + std::to_string(alpha1_size) + "\" />" + 
						"</node>\n";
			std::string alpha2_name;
			double p_value;
			double edge_weight;
			const std::map<const Alpha*, double>& edges = alpha.get_coincident_edges();
			for(const auto& edge_list : edges) {
				alpha2_name = (edge_list.first)->get_name();
				p_value = edge_list.second;
				edge_weight = ((1-p_value)*2);
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

/**
 * Execute a command and get the result.
 *
 * @param   cmd - The system command to run.
 * @return  The string command line output of the command.
 */
std::string Gexf::systemSTDOUT(std::string cmd) {

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

/**
 * Lookup component number and return corresponding colour code
 **/
std::string Gexf::componentLookup(int ret) {
	switch(ret % 138) {
		case 1  : return("#7C0051");
		case 2  : return("#76C3FF");
		case 3  : return("#001106");
		case 4  : return("#01ca10");
		case 5  : return("#e56800");
		case 6  : return("#0005a7");
		case 7  : return("#0a3700");
		case 8  : return("#0245ec");
		case 9  : return("#e3db00");
		case 10 : return("#671a00");
		case 11 : return("#00bf72");
		case 12 : return("#ff9dd6");
		case 13 : return("#f386ff");
		case 14 : return("#d57600");
		case 15 : return("#20000a");
		case 16 : return("#ace967");
		case 17 : return("#840019");
		case 18 : return("#0094f7");
		case 19 : return("#006695");
		case 20 : return("#ff5bf7");
		case 21 : return("#b9007e");
		case 22 : return("#ff89c7");
		case 23 : return("#360093");
		case 24 : return("#001733");
		case 25 : return("#005a69");
		case 26 : return("#d5dca2");
		case 27 : return("#61cb00");
		case 28 : return("#00338b");
		case 29 : return("#26a600");
		case 30 : return("#fdd268");
		case 31 : return("#ffbdca");
		case 32 : return("#b9e590");
		case 33 : return("#ff711b");
		case 34 : return("#8ce9c7");
		case 35 : return("#331300");
		case 36 : return("#a40053");
		case 37 : return("#02d6cf");
		case 38 : return("#640015");
		case 39 : return("#c90039");
		case 40 : return("#fecdb7");
		case 41 : return("#230025");
		case 42 : return("#008c51");
		case 43 : return("#ea21dd");
		case 44 : return("#7a63ff");
		case 45 : return("#b43b00");
		case 46 : return("#0068a7");
		case 47 : return("#01cf41");
		case 48 : return("#f6d633");
		case 49 : return("#93a200");
		case 50 : return("#61f0b9");
		case 51 : return("#ffbf91");
		case 52 : return("#990015");
		case 53 : return("#c4a8ff");
		case 54 : return("#02d9b8");
		case 55 : return("#f2007e");
		case 56 : return("#ac8b00");
		case 57 : return("#d0e072");
		case 58 : return("#ffa368");
		case 59 : return("#02b096");
		case 60 : return("#02c1cc");
		case 61 : return("#ff8de1");
		case 62 : return("#fc0072");
		case 63 : return("#dea3ff");
		case 64 : return("#c37d00");
		case 65 : return("#03c4f4");
		case 66 : return("#ff5dad");
		case 67 : return("#8990ff");
		case 68 : return("#d2d9d9");
		case 69 : return("#b744f6");
		case 70 : return("#380006");
		case 71 : return("#5a3700");
		case 72 : return("#002a23");
		case 73 : return("#392ad1");
		case 74 : return("#238100");
		case 75 : return("#ff673c");
		case 76 : return("#9bddff");
		case 77 : return("#b2e3c0");
		case 78 : return("#89ef68");
		case 79 : return("#00726c");
		case 80 : return("#004376");
		case 81 : return("#5d9fff");
		case 82 : return("#210052");
		case 83 : return("#009e40");
		case 84 : return("#ffac48");
		case 85 : return("#9d1fd7");
		case 86 : return("#d7d5ef");
		case 87 : return("#7cf22c");
		case 88 : return("#b7ddef");
		case 89 : return("#a8a7ff");
		case 90 : return("#008634");
		case 91 : return("#8c3300");
		case 92 : return("#8400a7");
		case 93 : return("#783900");
		case 94 : return("#0099b3");
		case 95 : return("#ff3628");
		case 96 : return("#386700");
		case 97 : return("#003c45");
		case 98 : return("#7f7100");
		case 99 : return("#001f5c");
		case 100: return("#b4a200");
		case 101: return("#272200");
		case 102: return("#ff9ea3");
		case 103: return("#ff325c");
		case 104: return("#e19100");
		case 105: return("#c80024");
		case 106: return("#ff6a61");
		case 107: return("#3c5300");
		case 108: return("#b9d400");
		case 109: return("#fc2109");
		case 110: return("#3a3b00");
		case 111: return("#024fd6");
		case 112: return("#5d7b00");
		case 113: return("#a90009");
		case 114: return("#006027");
		case 115: return("#6ab200");
		case 116: return("#016989");
		case 117: return("#4f005e");
		case 118: return("#016dbf");
		case 119: return("#29000b");
		case 120: return("#002700");
		case 121: return("#659700");
		case 122: return("#ffc545");
		case 123: return("#d7d0ff");
		case 124: return("#910085");
		case 125: return("#0064db");
		case 126: return("#ffa9c8");
		case 127: return("#dc0034");
		case 128: return("#620028");
		case 129: return("#ff5664");
		case 130: return("#006842");
		case 131: return("#441500");
		case 132: return("#ff7588");
		case 133: return("#9ce5da");
		case 134: return("#ff55cd");
		case 135: return("#c27fff");
		case 136: return("#004f31");
		case 137: return("#e2da91");
		case 0: return("#7bb4ff");
		default : return("black");
	}
	return("black");
}

