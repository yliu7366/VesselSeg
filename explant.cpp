#include <time.h>
#include <iostream>
#include <string>


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryThinningImageFilter3D.h"

#include "itkLine2DFilter.h"
#include "itkLineTransform2DFunction.h"

#include "CenterlineGraphConverter.h"

#include "VesselGraph.h"


const unsigned int TwoDimension = 2;
const unsigned int ThreeDimension = 3;

typedef double RealType;
typedef unsigned char InputPixelType;
typedef itk::Image< InputPixelType, TwoDimension > InputSliceType;
typedef double OutputPixelType;
typedef itk::Image< OutputPixelType, TwoDimension > OutputSliceType;
typedef itk::Vector< RealType, TwoDimension > TwoVectorType;
typedef itk::Image< TwoVectorType, TwoDimension > TwoVectorSliceType; 
typedef itk::MinimumMaximumImageCalculator< InputSliceType > MinMaxType;

typedef itk::Image< InputPixelType, ThreeDimension > InputVolumeType;
typedef itk::Image< OutputPixelType, ThreeDimension > OutputVolumeType;
typedef itk::Vector< RealType, ThreeDimension > ThreeVectorType;
typedef itk::Image< ThreeVectorType, ThreeDimension > ThreeVectorVolumeType; 

typedef CenterlineGraphConverter<InputPixelType, RealType, ThreeDimension> ConverterType;

typedef itk::Line2DFilter< OutputSliceType, OutputSliceType > LineFilterType;
typedef itk::LineTransform2DFunction<OutputSliceType, RealType> LineFunctionType;
typedef itk::VectorLineTransform2DFunction<TwoVectorSliceType, RealType> VectorLineFunctionType;

typedef itk::BinaryThresholdImageFilter<OutputSliceType, InputSliceType> ThresholdFilterType;
typedef itk::VotingBinaryIterativeHoleFillingImageFilter< InputSliceType > HoleFilterType;

typedef itk::BinaryThinningImageFilter< InputSliceType, InputSliceType > Thinning2FilterType;
typedef itk::BinaryThinningImageFilter3D< InputVolumeType, InputVolumeType > Thinning3FilterType;


//argv[1] input name
//argv[2] output name 
//argv[3] intensityMin
//argv[4] radiusFactor
//argv[5] max angles
//argv[6] downsamplingFactor
//argv[7] threshold

int main(int argc, char ** argv)
{
    if(argc != 8)
    {
        std::cout << "Wrong Number of Parameters" << std::endl;
        return 0;
    }

    std::string inputName = std::string(argv[1]);
    std::string outputName = std::string(argv[2]);
    RealType intensityMin = atof(argv[3]); //15; 
    int radiusFactor = atoi(argv[4]); //25;
    int maxNumAngles = atoi(argv[5]); //25; 
    int downsamplingFactor = atoi(argv[6]); //10;
    RealType threshold = atof(argv[7]);//0.65;
    
    //Set up reader and read input
    typedef itk::ImageFileReader< InputSliceType > ImageReaderType;
    ImageReaderType::Pointer   reader = ImageReaderType::New();
    reader->SetFileName( inputName ); 
    reader->Update();
    InputSliceType::Pointer inputImage = reader->GetOutput();

    //Set up writers
    typedef itk::ImageFileWriter< InputSliceType >  InputSliceWriterType;
    InputSliceWriterType::Pointer inputSliceWriter = InputSliceWriterType::New();
    
    typedef itk::ImageFileWriter< InputVolumeType >  InputVolumeWriterType;
    InputVolumeWriterType::Pointer inputVolumeWriter = InputVolumeWriterType::New();
    
    typedef itk::ImageFileWriter< OutputSliceType >  OutputImageWriterType;
    OutputImageWriterType::Pointer outputWriter = OutputImageWriterType::New();

    
    std::string inputName2 = std::string(argv[2]);
    inputName2.insert(inputName2.find_first_of("."),"Input");
    inputSliceWriter->SetFileName( inputName2.c_str() );
    inputSliceWriter->SetInput( inputImage );
    inputSliceWriter->Update();
    
    //timing parameters
    clock_t t1;
    clock_t t2;
    
    InputSliceType::SpacingType pixelSpacing = inputImage->GetSpacing();
    typedef std::vector<RealType> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(pixelSpacing[0]);
    spacingList.push_back(pixelSpacing[1]);
    
    RealType minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<RealType>() );
    RealType maxSpacing = *std::max_element(spacingList.begin(), spacingList.end(), std::less_equal<RealType>() );
    
   
    
    MinMaxType::Pointer minMaxCalc = MinMaxType::New();
    minMaxCalc->SetImage(inputImage);
    minMaxCalc->Compute();
    //std::cout << "Image Max " << minMaxCalc->GetMaximum() << " Min " << minMaxCalc->GetMinimum() << std::endl;
    RealType intensityMax = minMaxCalc->GetMaximum();
    //RealType intensityMin = minMaxCalc->GetMinimum();
 
    RealType lineFactor = radiusFactor*maxSpacing;
   
    OutputSliceType::Pointer downSampledInputImage = OutputSliceType::New();
    downSampledInputImage->CopyInformation(inputImage);
    OutputSliceType::RegionType downSampledRegion = downSampledInputImage->GetLargestPossibleRegion();
    OutputSliceType::SpacingType downSampledSpacing = downSampledInputImage->GetSpacing();
    
    for(unsigned int i = 0; i < TwoDimension; i++)
    {
        downSampledRegion.SetSize(i, int(ceil(RealType(downSampledRegion.GetSize(i))/RealType(downsamplingFactor))));
        downSampledSpacing[i] = downsamplingFactor*downSampledSpacing[i];
    }
    
    downSampledInputImage->SetSpacing( downSampledSpacing);
    downSampledInputImage->SetRegions( downSampledRegion );
    downSampledInputImage->SetRequestedRegionToLargestPossibleRegion();
    downSampledInputImage->Allocate(); 
    downSampledInputImage->FillBuffer(0);
    
    itk::ImageRegionIterator<OutputSliceType> downIter( downSampledInputImage, downSampledInputImage->GetLargestPossibleRegion() );   
    itk::ConstNeighborhoodIterator< InputSliceType >::RadiusType neighborhoodRadius;
    unsigned int numNeighbors = (2*downsamplingFactor+1) * (2*downsamplingFactor+1);
    neighborhoodRadius.Fill(downsamplingFactor);
    itk::ConstNeighborhoodIterator< InputSliceType > inputNeighIter(neighborhoodRadius, inputImage, 
        inputImage->GetRequestedRegion() );
    
    for(downIter.GoToBegin(); !downIter.IsAtEnd(); ++downIter)
    {
        OutputSliceType::IndexType currentIndex = downIter.GetIndex();
        OutputSliceType::PointType currentPoint;
        downSampledInputImage->TransformIndexToPhysicalPoint(currentIndex, currentPoint);
        OutputSliceType::IndexType currentUpIndex;
        bool isInImage = inputImage->TransformPhysicalPointToIndex(currentPoint, currentUpIndex);
        inputNeighIter.SetLocation(currentUpIndex);
        itk::Neighborhood<InputPixelType, TwoDimension> inHood = inputNeighIter.GetNeighborhood();
        RealType currentValue = 0.0;  
        int numAbove = 0; 
        //std::cout << "Neighborhood " << std::endl; 
        for (unsigned int i = 0; i < numNeighbors; i++)
        {
            RealType currentNeigh = RealType(inHood.GetElement(i));
            //std::cout << currentNeigh  << " ";

            if(/*currentNeigh >= currentValue &&*/ currentNeigh >= intensityMin)
            {
                currentValue += currentNeigh;
                numAbove++;
            }

        }
        downIter.Set(RealType(currentValue)/RealType(numNeighbors));
        //downSampledInputImage->SetPixel(currentIndex, currentValue);
    }
    
    std::string downName = std::string(argv[2]);
    downName.insert(downName.find_first_of("."),"DownSampled");
    outputWriter->SetFileName(downName.c_str() );
    outputWriter->SetInput( downSampledInputImage );
    outputWriter->Update();
   
    //Line   
    LineFunctionType::Pointer lineFunction = LineFunctionType::New();
    lineFunction->SetInputImage(downSampledInputImage);
    lineFunction->SetParameters(lineFactor, maxNumAngles);
    
    VectorLineFunctionType::Pointer vectorLineFunction = VectorLineFunctionType::New();
    vectorLineFunction->SetParameters(lineFactor, maxNumAngles);
    
    LineFilterType::Pointer lineFilter = LineFilterType::New();
    lineFilter->SetInput(downSampledInputImage);
    lineFilter->SetIntensityMaximum( intensityMax );
    lineFilter->SetIntensityMinimum( intensityMin );
    lineFilter->SetLineFunction(&(*lineFunction));
    lineFilter->SetVectorLineFunction(&(*vectorLineFunction));
    lineFilter->Allocate();

    //std::cout << "Line computation "<< std::endl;
    t1 = clock(); 
    lineFilter->ComputeLineFilter();
    t2 = clock(); 
    //std::cout << "Line time "<< t2-t1 << std::endl;
    OutputSliceType::Pointer lineImage = lineFilter->GetLineImage();
    
    std::string lineName = std::string(argv[2]);
    lineName.insert(lineName.find_first_of("."),"Line");
    outputWriter->SetFileName(lineName.c_str() );
    outputWriter->SetInput( lineImage );
    outputWriter->Update();
    
    ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
    thresholdFilter->SetInput(lineImage);
    thresholdFilter->SetLowerThreshold(threshold);
    thresholdFilter->SetUpperThreshold(1.0);
    thresholdFilter->SetInsideValue(1);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();
    InputSliceType::Pointer thresholdImage = thresholdFilter->GetOutput();
    
    std::string thresholdName = std::string(argv[2]);
    thresholdName.insert(thresholdName.find_first_of("."),"Thresholded");
    inputSliceWriter->SetFileName( thresholdName.c_str() );
    inputSliceWriter->SetInput( thresholdImage );
    //inputSliceWriter->Update();
    
    HoleFilterType::Pointer holeFilter = HoleFilterType::New();
    holeFilter->SetBackgroundValue(   0 );
    holeFilter->SetForegroundValue( 1 );
    holeFilter->SetMaximumNumberOfIterations( 2 );
    holeFilter->SetMajorityThreshold( 2 );
    holeFilter->SetInput(thresholdImage);
    holeFilter->Update();
    InputSliceType::Pointer holeFilledImage = holeFilter->GetOutput();
    
    std::string holeName = std::string(argv[2]);
    holeName.insert(holeName.find_first_of("."),"HoleFilled");
    inputSliceWriter->SetFileName( holeName.c_str() );
    inputSliceWriter->SetInput( holeFilledImage );
    inputSliceWriter->Update();
    
    Thinning2FilterType::Pointer thinning2Filter = Thinning2FilterType::New();
    thinning2Filter->SetInput(holeFilledImage);
    thinning2Filter->Update();
    InputSliceType::Pointer thinned2Image = thinning2Filter->GetOutput();
    
    std::string thinName = std::string(argv[2]);
    thinName.insert(thinName.find_first_of("."),"Thinned");
    inputSliceWriter->SetFileName( thinName.c_str() );
    inputSliceWriter->SetInput( thinned2Image );
    //inputSliceWriter->Update();
    
    //Convert to 3D
    ConverterType converter;
    
    InputVolumeType::Pointer volume = InputVolumeType::New();
    InputVolumeType::SpacingType volumeSpacing;
    volumeSpacing[0] = downSampledInputImage->GetSpacing()[0];
    volumeSpacing[1] = downSampledInputImage->GetSpacing()[1];
    volumeSpacing[2] = 1.0;
    
    InputVolumeType::PointType volumeOrigin;
    volumeOrigin[0] = downSampledInputImage->GetOrigin()[0];
    volumeOrigin[1] = downSampledInputImage->GetOrigin()[1];
    volumeOrigin[2] = 0.0;
    
    InputVolumeType::RegionType volumeRegion;
    InputVolumeType::RegionType::SizeType volumeSize;
    volumeSize[0] = downSampledInputImage->GetLargestPossibleRegion().GetSize()[0];
    volumeSize[1] = downSampledInputImage->GetLargestPossibleRegion().GetSize()[1];
    volumeSize[2] = 1;
    volumeRegion.SetSize(volumeSize);
    
    InputVolumeType::IndexType volumeStart;
    volumeStart[0] = downSampledInputImage->GetLargestPossibleRegion().GetIndex()[0];
    volumeStart[1] = downSampledInputImage->GetLargestPossibleRegion().GetIndex()[1];
    volumeStart[2] = 0; 
    volumeRegion.SetIndex(volumeStart);
    
    volume->SetRegions(volumeRegion);
    volume->SetRequestedRegionToLargestPossibleRegion ();
    volume->Allocate();
    volume->FillBuffer(0);
    
    itk::ImageRegionIterator<InputSliceType> twoDIter( holeFilledImage, holeFilledImage->GetLargestPossibleRegion() );
    for(twoDIter.GoToBegin(); !twoDIter.IsAtEnd(); ++twoDIter)
    {
        InputPixelType currentValue = twoDIter.Get();
        InputSliceType::IndexType current2Index = twoDIter.GetIndex();
        
        InputVolumeType::IndexType current3Index;
        current3Index[0] = current2Index[0];
        current3Index[1] = current2Index[1];
        current3Index[2] = 0;
        
        volume->SetPixel(current3Index, currentValue);
    }
   
    //Do 3d thinning 
    Thinning3FilterType::Pointer thinning3Filter = Thinning3FilterType::New();
    thinning3Filter->SetInput(volume);
    thinning3Filter->Update();
    InputVolumeType::Pointer centerlineVolume = thinning3Filter->GetOutput();
    
    InputSliceType::Pointer centerlineSlice = InputSliceType::New();
    centerlineSlice->CopyInformation(inputImage);
    centerlineSlice->SetSpacing( downSampledSpacing);
    centerlineSlice->SetRegions( downSampledRegion );
    centerlineSlice->SetRequestedRegionToLargestPossibleRegion();
    centerlineSlice->Allocate(); 
    centerlineSlice->FillBuffer(0);
    itk::ImageRegionIterator<InputSliceType> centerIter( centerlineSlice, centerlineSlice->GetLargestPossibleRegion() );
    for(centerIter.GoToBegin(); !centerIter.IsAtEnd(); ++centerIter)
    {
        InputSliceType::IndexType current2Index = centerIter.GetIndex();
        
        InputVolumeType::IndexType current3Index;
        current3Index[0] = current2Index[0];
        current3Index[1] = current2Index[1];
        current3Index[2] = 0;
        
        centerlineSlice->SetPixel(current2Index, centerlineVolume->GetPixel(current3Index));
    }

    std::string centerSliceName = std::string(argv[2]);
    centerSliceName.insert(centerSliceName.find_first_of("."),"CenterlineSlice");
    inputSliceWriter->SetFileName( centerSliceName.c_str() );
    inputSliceWriter->SetInput( centerlineSlice );
    inputSliceWriter->Update();
    
    //Prepare other volumes for graph
    std::string centerVolumeName = std::string(argv[2]);
    //centerVolumeName.insert(centerVolumeName.find_first_of("."),"CenterlineVolume");
    inputVolumeWriter->SetFileName( centerVolumeName.c_str() );
    inputVolumeWriter->SetInput( centerlineVolume );
    inputVolumeWriter->Update();
    
    OutputVolumeType::Pointer radiusVolume = OutputVolumeType::New();
    radiusVolume->CopyInformation(centerlineVolume);
    radiusVolume->SetRegions(centerlineVolume->GetLargestPossibleRegion() );
    radiusVolume->SetRequestedRegionToLargestPossibleRegion();
    radiusVolume->Allocate();
    radiusVolume->FillBuffer(1.0); 
    
    OutputVolumeType::Pointer probabilityVolume = OutputVolumeType::New();
    probabilityVolume->CopyInformation(centerlineVolume);
    probabilityVolume->SetRegions(centerlineVolume->GetLargestPossibleRegion() );
    probabilityVolume->SetRequestedRegionToLargestPossibleRegion();
    probabilityVolume->Allocate();
    probabilityVolume->FillBuffer(1.0);
    
    ThreeVectorVolumeType::Pointer tangentVolume = ThreeVectorVolumeType::New();
    tangentVolume->CopyInformation(centerlineVolume);
    tangentVolume->SetRegions(centerlineVolume->GetLargestPossibleRegion() );
    tangentVolume->SetRequestedRegionToLargestPossibleRegion();
    tangentVolume->Allocate();
    ThreeVectorType zeroVector;
    zeroVector.Fill(0);
    tangentVolume->FillBuffer(zeroVector);
    
    //Convert volume to graph
    ConverterType::VesselGraphType vesselGraph = converter.ImageToGraph(centerlineVolume, radiusVolume, probabilityVolume, tangentVolume);
 
    std::ofstream outputFile;
    std::string outputFileName(outputName.c_str(), outputName.find_first_of("."));
    outputFileName.append("Graph.vessel");
    outputFile.open( outputFileName.c_str(), std::ios::out);
    //outputFile << vesselGraph;
    outputFile.close();
    
    //Pretty graph and extract metrics
    ConverterType::VesselGraphType vesselGraphRemove = vesselGraph.RemoveRedundentJunctions();

    std::ofstream outputFileRemove;
    std::string outputFileNameRemove(outputName.c_str(), outputName.find_first_of("."));
    outputFileNameRemove.append("GraphRemove.vessel");
    outputFileRemove.open( outputFileNameRemove.c_str(), std::ios::out);
    outputFileRemove << vesselGraphRemove;
    outputFileRemove.close();
    
    
    std::ofstream outputFileRemoveGML;
    std::string outputFileNameRemoveGML(outputName.c_str(), outputName.find_first_of("."));
    outputFileNameRemoveGML.append("GraphRemove.gml");
    outputFileRemoveGML.open( outputFileNameRemoveGML.c_str(), std::ios::out);
    vesselGraphRemove.WriteGMLJunctionOnly(outputFileRemoveGML);
    outputFileRemoveGML.close();

    
    std::ofstream pythonFile;
    std::string pythonFileName(outputName.c_str(), outputName.find_first_of("."));
    pythonFileName.append("Metrics.py");
    pythonFile.open( pythonFileName.c_str(), std::ios::out);
    vesselGraphRemove.ReportMetricsAsPython(pythonFile);
    pythonFile.close();
    
    return 0;
}
