#include <string>
#include <vector>
#include "elements.h"
#include "constants.h"
#include "dataset.h"

struct TParameters;


class Gexf
{
    public:
        static void run( DataSet& dataset, const std::string& prefix);
        
    private:
        static std::string systemSTDOUT(std::string cmd);
	static std::string componentLookup(int ret);
};
