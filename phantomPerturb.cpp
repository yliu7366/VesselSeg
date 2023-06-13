
#include <iostream>
#include <fstream>
#include <string>


#include "PhantomGenerator.h"

typedef double CoordinateType;
const unsigned int Dimension = 3;

typedef PhantomGenerator<CoordinateType, Dimension> PhantomType;


int main(int argc, char ** argv)
{

    std::cout << "Phantom Perturber" << std::endl;
    if(argc != 4)
    {
        std::cout << "Wrong Number of Parameters" << std::endl;
        return 0;
    }

    std::string inputName = std::string(argv[1]);
    std::string outputName = std::string(argv[2]);
    CoordinateType randomness = atof(argv[3]);
   
    
    PhantomType::GraphType phantomGraph;
    
    std::ifstream inputFile;
    inputFile.open(inputName.c_str(), std::ios::in);
    inputFile >> phantomGraph;
    inputFile.close();

    PhantomType phantom;
    
    PhantomType::GraphType phantomPertGraph = phantom.PerturbGraph(phantomGraph , randomness);
    
    std::ofstream outputFile;
    std::string graphName = outputName;
    graphName.append(".vessel");
    outputFile.open( graphName.c_str(), std::ios::out);
    outputFile << phantomPertGraph;
    outputFile.close();
    
    std::ofstream pythonFile;
    std::string pythonName = outputName;
    pythonName.append(".py");
    pythonFile.open( pythonName.c_str(), std::ios::out);
    phantomPertGraph.ReportMetricsAsPython(pythonFile);
    pythonFile.close();

    return 0;
}
