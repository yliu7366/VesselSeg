
#include <iostream>
#include <fstream>
#include <string>

#include "VesselGraph.h"

typedef double CoordinateType;
const unsigned int Dimension = 3; 
typedef VesselGraph<Dimension, CoordinateType> VesselGraphType; 

int main(int argc, char ** argv)
{
    std::string inputName(argv[1]);
    std::ifstream inputfile;
    inputfile.open(inputName.c_str(), std::ios::in);
    VesselGraphType vesselGraph;
    inputfile >> vesselGraph;
    inputfile.close();    
    
    std::string graphGMLName(argv[2]);
    std::ofstream outputGMLFile;
    outputGMLFile.open( graphGMLName.c_str(), std::ios::out);
    vesselGraph.WriteGML(outputGMLFile);
    outputGMLFile.close();
    
    return 0;
}
