#include <string>
#include <vector>
#include "elements.h"
#include "constants.h"
#include "dataset.h"

struct TParameters;


class Lineage
{
    public:
        static int run(DataSet& dataset, std::string source_path, std::string call_path, const std::string& phylogeny, int num_cores, bool Rmsgs);
        
    private:
	static std::string systemSTDOUT(std::string cmd);
};
