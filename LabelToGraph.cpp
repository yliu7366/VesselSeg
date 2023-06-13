#include <time.h>
#include <iostream>
#include <string>
#include <limits>
#include <cmath>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkIndex.h"
#include "itkImageRegionIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"

//not in ITK 
#include "itkBinaryThinningImageFilter3D.h"

#include "itkMedialVesselFilter.h"
#include "itkMedialness3DFunction.h"
#include "itkHessianVesselFilter.h"
#include "itkLineTransform3DFunction.h"

#include "CenterlineGraphConverter.h"

#include "VesselGraph.h"
#include "VesselSegment.h"
#include "VesselNode.h"
#include "FixedArray.h"

const unsigned int Dimension = 3;
typedef unsigned int LabelType;
typedef double RealType;
typedef itk::Image< LabelType, Dimension > LabelImageType;
typedef CenterlineGraphConverter<LabelType, RealType, Dimension> ConverterType;

int main(int argc, char ** argv)
{
  if ( argc < 3 )
  {
    std::cerr << "Missing Parameters: "
    << argv[0]
    << " Input_File"
    << " Output_File"
    << std::endl;
    return 0;
  }  

  std::string inputName = std::string(argv[1]);
  std::string outputName = std::string(argv[2]);

  //Set up reader and read input
  typedef itk::ImageFileReader< LabelImageType > ImageReaderType;
  ImageReaderType::Pointer   reader = ImageReaderType::New();
  reader->SetFileName( inputName ); 
  reader->Update();
  LabelImageType::Pointer inputImage = reader->GetOutput();

  long t1, t2;
  std::cout << "Build Vessel Network" << std::endl;
  std::cout.flush();
  t1 = clock();
  ConverterType converter;
  ConverterType::VesselGraphType vesselGraph = converter.ImageToGraph(inputImage);
  vesselGraph = vesselGraph.RemoveRedundentJunctions();
  t2 = clock(); 
  std::cout << "Build Vessel Network time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;

  std::ofstream outputFile;
  outputFile.open( outputName.c_str(), std::ios::out);
  outputFile << vesselGraph;
  outputFile.close();
  
  return 0;
}
