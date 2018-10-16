#include <string>
#include <vector>
#include "elements.h"
#include "constants.h"
#include "dataset.h"

struct TParameters;


class Network
{
    public:
        static int run( DataSet& dataset, std::string source_path, std::string call_path, const std::string& phylogeny, bool Rmsgs );
        
    private:
	static std::string systemSTDOUT(std::string cmd);
};
