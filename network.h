#include <string>
#include <vector>
#include "elements.h"
#include "constants.h"
#include "dataset.h"

struct TParameters;


class Network
{
    public:
        static int run( DataSet& dataset);
        
    private:
	static std::string systemSTDOUT(std::string cmd);
};
