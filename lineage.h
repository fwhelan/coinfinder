#include <string>
#include <vector>
#include "elements.h"
#include "constants.h"
#include "dataset.h"

struct TParameters;


class Lineage
{
    public:
        static int run(DataSet& dataset);
        
    private:
	static std::string systemSTDOUT(std::string cmd);
};
