#include <string>
#include <vector>
#include "elements.h"
#include "constants.h"
#include "dataset.h"

struct TParameters;


class Lineage
{
    public:
        static int run(DataSet& dataset, const std::string& source_path, const std::string& call_path, const std::string& phylogeny, const std::string& gene_pa, int num_cores, bool Rmsgs, const std::string& prefix );
	static std::string systemSTDOUT(std::string cmd);
};
