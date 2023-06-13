#include <time.h>
#include <iostream>
#include <string>
#include <limits>
#include <cmath>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkBoxMeanImageFilter.h"

//not in ITK
#include "itkBinaryThinningImageFilter3D.h"


// our stuff
#include "CenterlineGraphConverter.h"
#include "VesselGraph.h"
#include "VesselSegment.h"
#include "VesselNode.h"
#include "FixedArray.h"

const unsigned int Dimension = 3;
typedef unsigned int LabelType;
typedef double RealType;
typedef double InputPixelType;
typedef itk::Image< InputPixelType, Dimension > InputImageType;
typedef double OutputPixelType;
typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
typedef itk::Image< LabelType, Dimension > LabelImageType;
typedef VesselGraph<Dimension, RealType> VesselGraphType; 
     
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< InputImageType,  OutputImageType >  GradientFilterType;     
typedef itk::SigmoidImageFilter<InputImageType, OutputImageType >  SigmoidFilterType; 
typedef itk::FastMarchingImageFilter< InputImageType, OutputImageType >  FastMarchingFilterType;
typedef itk::BinaryThresholdImageFilter< OutputImageType, LabelImageType > ThresholdFilterType;
typedef itk::VotingBinaryIterativeHoleFillingImageFilter< LabelImageType > HoleFilterType;
typedef CenterlineGraphConverter<LabelType, RealType, Dimension> ConverterType;
typedef itk::BinaryThinningImageFilter3D< LabelImageType, LabelImageType> ThinningFilterType;
typedef itk::SignedMaurerDistanceMapImageFilter< LabelImageType, OutputImageType > DistanceFilterType;
typedef itk::BoxMeanImageFilter< InputImageType, OutputImageType > MeanFilterType;

int main(int argc, char ** argv)
{
    std::cout << "Vessel Segmenter using Fast Marching with Thresholding" << std::endl;
    if(argc != 9)
    {
        std::cout << "Wrong Number of Parameters" << std::endl;
        return 0;
    }

    std::string inputName = std::string(argv[1]);
    std::string graphName = std::string(argv[2]);
    std::string outputName = std::string(argv[3]);


    RealType alpha = atof(argv[4]);
    RealType beta = atof(argv[5]);
    RealType stoppingTime = atof(argv[6]);
    RealType timeThreshold = atof(argv[7]);
    int unconnectedSegments = atoi(argv[8]);
    
    VesselGraphType inputGraph;
    std::ifstream inputGraphFile;
    inputGraphFile.open( graphName.c_str(), std::ios::in);
    inputGraphFile >> inputGraph;
    inputGraphFile.close();



    //Set up reader and read input
    typedef itk::ImageFileReader< InputImageType > ImageReaderType;
    ImageReaderType::Pointer   reader = ImageReaderType::New();
    reader->SetFileName( inputName ); 
    reader->Update();
    InputImageType::Pointer inputImage = reader->GetOutput();

    //Set up readers and writers
    typedef itk::ImageFileWriter< OutputImageType >  OutputImageWriterType;
    OutputImageWriterType::Pointer outputWriter = OutputImageWriterType::New();
       
    typedef itk::ImageFileWriter< LabelImageType >  LabelImageWriterType;
    LabelImageWriterType::Pointer labelWriter = LabelImageWriterType::New();
    
    


    //timing parameters
    clock_t t1;
    clock_t t2;
    
    
    SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
    sigmoid->SetOutputMinimum(  0.0  );
    sigmoid->SetOutputMaximum(  1.0  );
    sigmoid->SetAlpha( alpha );
    sigmoid->SetBeta(  beta  );
    sigmoid->SetInput(inputImage);
    sigmoid->Update();
    OutputImageType::Pointer sigmoidImage = sigmoid->GetOutput();
    
    std::string sigmoidName = std::string(outputName);
    sigmoidName.insert(sigmoidName.find_first_of("."),"Sigmoid");
    outputWriter->SetFileName( sigmoidName.c_str() );
    outputWriter->SetInput( sigmoidImage );
    outputWriter->Update();

    FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
    FastMarchingFilterType::NodeContainer::Pointer seeds = FastMarchingFilterType::NodeContainer::New();
    seeds->Initialize();
    
    unsigned int i = 0;
    for(VesselGraphType::ConstIterator segmentIter = inputGraph.Begin(); segmentIter != inputGraph.End(); ++segmentIter)
    {
        VesselGraphType::SegmentType currentSegment = segmentIter->second;
        VesselGraphType::SegmentType::Iterator nodeIter = currentSegment.Begin();
        for(; nodeIter != currentSegment.End(); ++nodeIter)
        {
            VesselGraphType::NodeIDType currentID = *(nodeIter);
            VesselGraphType::NodeType currentNode = inputGraph.GetNode(currentID);
            VesselGraphType::NodeType::PositionType currentPosition = currentNode.GetPosition();
            OutputImageType::PointType seedPosition;
            seedPosition[0] = currentPosition[0];
            seedPosition[1] = currentPosition[1];
            seedPosition[2] = currentPosition[2];
            OutputImageType::IndexType  seedIndex;
            bool inImage = inputImage->TransformPhysicalPointToIndex(seedPosition, seedIndex);
            if(inImage == true)
            {
                FastMarchingFilterType::NodeType  node;
                const double seedValue = 0.0; 
                node.SetValue( seedValue );
                node.SetIndex( seedIndex );
                seeds->InsertElement( i, node );
                i++;
            }
            
            
        }
    }
    fastMarching->SetTrialPoints(  seeds  );
    fastMarching->SetOutputSize( inputImage->GetBufferedRegion().GetSize() );
    fastMarching->SetStoppingValue(  stoppingTime  );
    fastMarching->SetInput(sigmoidImage);
    std::cout << "Fast Marching computation "<< std::endl;
    t1 = clock(); 
    fastMarching->Update();
    t2 = clock(); 
    std::cout << "Fast Marching time "<< t2-t1 << std::endl;
    OutputImageType::Pointer fastMarchingImage = fastMarching->GetOutput();
    
    std::string fastMarchingName = std::string(outputName);
    fastMarchingName.insert(fastMarchingName.find_first_of("."),"FastMarching");
    outputWriter->SetFileName( fastMarchingName.c_str() );
    outputWriter->SetInput( fastMarchingImage );
    outputWriter->Update();
  


    ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
    thresholdFilter->SetLowerThreshold(0.0);
    thresholdFilter->SetUpperThreshold(timeThreshold);
    thresholdFilter->SetInsideValue(1);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetInput(fastMarchingImage);
    thresholdFilter->Update();
    LabelImageType::Pointer thresholdedImage = thresholdFilter->GetOutput();
    
     
    HoleFilterType::Pointer holeFilter = HoleFilterType::New();
    holeFilter->SetBackgroundValue(   0 );
    holeFilter->SetForegroundValue( 1 );
    holeFilter->SetMaximumNumberOfIterations( 20 );
    holeFilter->SetMajorityThreshold( 0 );
    
    //Do hole fill
    holeFilter->SetInput(thresholdedImage);
    holeFilter->Update();
    LabelImageType::Pointer holeFilledImage = holeFilter->GetOutput();
    std::string holeFilledName = std::string(outputName);
    holeFilledName.insert(holeFilledName.find_first_of("."),"HoleFilled");
    labelWriter->SetFileName( holeFilledName.c_str() );
    labelWriter->SetInput( holeFilledImage );
    labelWriter->Update();
    
    ThinningFilterType::Pointer thinningFilter = ThinningFilterType::New();
    thinningFilter->SetInput(holeFilledImage);
    thinningFilter->Update();
    LabelImageType::Pointer centerlineImage = thinningFilter->GetOutput();
    std::string centerlineName = std::string(outputName);
    centerlineName.insert(centerlineName.find_first_of("."),"Centerline");
    labelWriter->SetFileName( centerlineName.c_str() );
    labelWriter->SetInput( centerlineImage );
    labelWriter->Update();
    
    
    DistanceFilterType::Pointer distanceFilter = DistanceFilterType::New();
    distanceFilter->SetInput(holeFilledImage);
    distanceFilter->SetSquaredDistance(true);
    distanceFilter->SetUseImageSpacing(true);
    distanceFilter->Update();
    OutputImageType::Pointer radiusImage = distanceFilter->GetOutput();
    std::string radiusName = std::string(outputName);
    radiusName.insert(radiusName.find_first_of("."),"Radius");
    outputWriter->SetFileName( radiusName.c_str() );
    outputWriter->SetInput( radiusImage );
    outputWriter->Update();
    


    std::cout << "Build Vessel Network" << std::endl;
    std::cout.flush();
    t1 = clock();
    ConverterType converter;
    ConverterType::VesselGraphType vesselGraph = converter.ImageToGraph(holeFilledImage);
    vesselGraph.RemoveUnconnectedSegments(unconnectedSegments);  
    vesselGraph = vesselGraph.RemoveRedundentJunctions();
    t2 = clock(); 
    std::cout << "Build Vessel Network time "<< t2-t1 << std::endl;
 
    std::ofstream outputFile;
    std::string outputFileName(outputName.c_str(), outputName.find_first_of("."));
    outputFileName.append("Refine.vessel");
    outputFile.open( outputFileName.c_str(), std::ios::out);
    outputFile << vesselGraph;
    outputFile.close();
    
    std::ofstream outputGMLFile;
    std::string graphGMLName(outputName.c_str(), outputName.find_first_of("."));
    graphGMLName.append("Refine.gml");
    outputGMLFile.open( graphGMLName.c_str(), std::ios::out);
    vesselGraph.WriteGML(outputGMLFile);
    outputGMLFile.close();
    
    std::ofstream pythonFile;
    std::string pythonFileName(outputName.c_str(), outputName.find_first_of("."));
    pythonFileName.append("Metrics.py");
    pythonFile.open( pythonFileName.c_str(), std::ios::out);
    vesselGraph.ReportMetricsAsPython(pythonFile);
    pythonFile.close();
 
   
    std::cout << "Build Interpolated Graph" << std::endl;
    std::cout.flush();
    t1 = clock();
    ConverterType::VesselGraphType interpGraph = vesselGraph.Interpolate(4);
    t2 = clock(); 
    std::cout << "Build Interpolated Network time "<< t2-t1 << std::endl;
 
    std::string outputFileSubGraphName(outputName.c_str(), outputName.find_first_of("."));
    outputFileSubGraphName.append("Interp.vessel");
    outputFile.open( outputFileSubGraphName.c_str(), std::ios::out);
    outputFile << interpGraph;
    outputFile.close();
    
    
    
    std::cout << "Build Label Map from Graph" << std::endl;
    std::cout.flush();
    t1 = clock();
    LabelImageType::Pointer upSampledImage = LabelImageType::New();
    upSampledImage->CopyInformation(inputImage);
    LabelImageType::RegionType upSampledRegion = inputImage->GetLargestPossibleRegion();
    LabelImageType::SpacingType upSampledSpacing = inputImage->GetSpacing();
    for(unsigned int i = 0; i < 3; i++)
    {
        upSampledRegion.SetSize(i, 1*upSampledRegion.GetSize(i));
        upSampledSpacing[i] = upSampledSpacing[i]/1.0;
    }
    
    upSampledImage->SetSpacing(upSampledSpacing);
    upSampledImage->SetRegions( upSampledRegion );
    upSampledImage->SetRequestedRegionToLargestPossibleRegion();
    upSampledImage->Allocate(); 
    upSampledImage->FillBuffer(0);
    InputImageType::Pointer fromGraphImage = converter.GraphToImage(upSampledImage, 1, interpGraph);
    t2 = clock(); 
    std::cout << "Build Label Map from Graph time "<< t2-t1 << std::endl;
       
    std::string fromGraphName = std::string(outputName);
    fromGraphName.insert(fromGraphName.find_first_of("."),"FromGraphRefine");
    outputWriter->SetFileName( fromGraphName.c_str() );
    outputWriter->SetInput( fromGraphImage );
    outputWriter->Update();

   
    return 0;
}
