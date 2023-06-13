
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <vector>
#include <cmath>


#include "PhantomGenerator.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"


typedef double CoordinateType;
const unsigned int Dimension = 3;

typedef PhantomGenerator<CoordinateType, Dimension> PhantomType;
typedef itk::Image<CoordinateType, Dimension> ImageType;


int main(int argc, char ** argv)
{

    std::cout << "Phantom Generator" << std::endl;
    if(argc != 19)
    {
        std::cout << "Wrong Number of Parameters" << std::endl;
        return 0;
    }

    std::string outputName = std::string(argv[1]);
    CoordinateType rootRadius = atof(argv[2]);
    
    CoordinateType leftRadialScaling = atof(argv[3]);
    CoordinateType rightRadialScaling = atof(argv[4]);
    CoordinateType straightRadialScaling = atof(argv[5]);
    CoordinateType leftLengthScaling = atof(argv[6]);
    CoordinateType rightLengthScaling =  atof(argv[7]);
    CoordinateType straightLengthScaling = atof(argv[8]);
    CoordinateType angleLeftIn = atof(argv[9]);
    CoordinateType angleLeftOut = atof(argv[10]);
    CoordinateType angleRightIn = atof(argv[11]);
    CoordinateType angleRightOut = atof(argv[12]);
    CoordinateType angleStraightIn = atof(argv[13]);
    CoordinateType angleStraightOut = atof(argv[14]);
    CoordinateType randomnessFactor = atof(argv[15]);
    
    unsigned int numScales = atoi(argv[16]);
    unsigned int maxSegmentLength = atoi(argv[17]);
    unsigned int minSegmentLength = atoi(argv[18]);

    PhantomType phantom;
    
    PhantomType::GraphType::NodeType root;
    PhantomType::GraphType::NodeType::PositionType rootPosition;
    rootPosition.Fill(0.0);
    PhantomType::GraphType::NodeType::VectorType rootTangent;
    rootTangent.Fill(0.0);
    rootTangent[0] = 1.0;
    PhantomType::GraphType::NodeType::VectorType rootNormal;
    rootNormal.Fill(0.0);
    rootNormal[1] = 1.0;
    PhantomType::GraphType::NodeType::VectorType rootBiNormal;
    rootBiNormal.Fill(0.0);
    rootBiNormal[2] = 1.0;
    
    root.SetID(1);
    root.SetPosition(rootPosition);
    root.SetTangent(rootTangent);
    root.SetNormal(rootNormal);
    root.SetBiNormal(rootBiNormal);
    root.SetRadius(rootRadius);
    root.SetMedialness(0.0);
    
    phantom.SetParameters(leftRadialScaling, rightRadialScaling, straightRadialScaling, leftLengthScaling, rightLengthScaling, straightLengthScaling, angleLeftIn, angleLeftOut, angleRightIn, angleRightOut, angleStraightIn, angleStraightOut, numScales, maxSegmentLength, minSegmentLength, root, randomnessFactor);
    
    PhantomType::GraphType phantomGraph = phantom.GenerateGraph();
    

    
    std::ofstream outputFile;
    std::string graphName = outputName;
    graphName.append(".vessel");
    outputFile.open( graphName.c_str(), std::ios::out);
    outputFile << phantomGraph;
    outputFile.close();
    
    std::ofstream outputGMLFile;
    std::string graphGMLName = outputName;
    graphGMLName.append(".gml");
    outputGMLFile.open( graphGMLName.c_str(), std::ios::out);
    phantomGraph.WriteGML(outputGMLFile);
    outputGMLFile.close();
    
/*    
    std::ofstream pythonFile;
    std::string pythonName = outputName;
    pythonName.append(".py");
    pythonFile.open( pythonName.c_str(), std::ios::out);
    phantomGraph.ReportMetricsAsPython(pythonFile);
    pythonFile.close();
*/  
/*    
    //to remove
    PhantomType::GraphType phantomGraphInterp = phantom.InterpolateGraph(phantomGraph, 3);
    std::string graphInterpName = outputName;
    graphInterpName.append("Interp.vessel");
    outputFile.open( graphInterpName.c_str(), std::ios::out);
    outputFile << phantomGraphInterp;
    outputFile.close();
    
    PhantomType::GraphType phantomPertGraph = phantom.PerturbGraph(phantomGraph, 0.4);
    std::string graphPertName = outputName;
    graphPertName.append("Pert.vessel");
    outputFile.open( graphPertName.c_str(), std::ios::out);
    outputFile << phantomPertGraph;
    outputFile.close();
    
    PhantomType::GraphType phantomPertGraphInterp = phantom.InterpolateGraph(phantomPertGraph, 3);
    std::string graphPertInterpName = outputName;
    graphPertInterpName.append("PertInterp.vessel");
    outputFile.open( graphPertInterpName.c_str(), std::ios::out);
    outputFile << phantomPertGraphInterp;
    outputFile.close();
*/    
    /*
    PhantomType::GraphType phantomGraph2;
    std::ifstream inputFile;
    inputFile.open(outputName.c_str(), std::ios::in);
    inputFile >> phantomGraph2;
    inputFile.close();
    

    for(PhantomType::GraphType::ConstIterator graphIter = phantomGraph.Begin(); graphIter != phantomGraph.End(); graphIter++)
    {
        PhantomType::GraphType::SegmentType::IDType currentSegmentID = graphIter->first;
        PhantomType::GraphType::SegmentType currentSegment = graphIter->second;
        
        std::cout << "Chord Length 1 " << currentSegment.GetChordLength() << std::endl;
        std::cout << "Chord Length 2 " << (phantomGraph2.GetSegment(currentSegmentID)).GetChordLength() << std::endl;
    }
    */
    
    PhantomType::VectorType spacing;
    spacing.Fill(1.0);

    ImageType::Pointer phantomImageLowRes = phantom.GraphToImage(phantomGraph, 2*rootRadius, spacing);

//    std::ofstream pythonFile;
//    std::string pythonFileName("phantomMetrics.py");
//    pythonFile.open( pythonFileName.c_str(), std::ios::out);
//    phantomGraph.ReportMetricsAsPython(pythonFile);
//    pythonFile.close();
    
    typedef itk::ImageFileWriter< ImageType >  ImageWriterType;
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    
    std::string imageName = outputName;
    imageName.append("Image.mhd");
    imageWriter->SetFileName( imageName );
    imageWriter->SetInput( phantomImageLowRes);
    imageWriter->Update();
    
    return 0;
}
