
#include <iostream>
#include <string>

#include "VesselGraph.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkDiscreteGaussianImageFilter.h"

const unsigned int Dimension = 3;
typedef double RealType;
typedef double InputPixelType;
typedef itk::Image< InputPixelType, Dimension > InputImageType;
typedef double OutputPixelType;
typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

typedef VesselGraph<Dimension, RealType> GraphType; 

typedef itk::DiscreteGaussianImageFilter<InputImageType, OutputImageType> GaussianFilterType;

RealType BoxMuller(RealType mean, RealType sdev)
{
    RealType randNum1 = RealType(rand()) / RealType(RAND_MAX);
    RealType randNum2 = RealType(rand()) / RealType(RAND_MAX);
    RealType term1 = sqrt(-2.0*log(randNum1));
    RealType term2 = 2.0*M_PI*randNum2;
    
    RealType gaussRand = term1*sin(term2)*sdev + mean;
    return(gaussRand);
}

int main(int argc, char ** argv)
{
    std::string inputName = std::string(argv[1]);
    std::string graphName = std::string(argv[2]);
    std::string outputName = std::string(argv[3]);
    RealType sigma = atof(argv[4]);
    RealType backgroundValue = atof(argv[5]);
    RealType minValue = atof(argv[6]);
    RealType maxValue = atof(argv[7]);
    RealType noiseMean = atof(argv[8]);
    RealType noiseSTD = atof(argv[9]);

    //Set up reader and read input
    typedef itk::ImageFileReader< InputImageType > ImageReaderType;
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName( inputName ); 
    reader->Update();
    InputImageType::Pointer inputImage = reader->GetOutput();
    
    GaussianFilterType::Pointer gaussFilter = GaussianFilterType::New();
    gaussFilter->SetInput(inputImage);
    gaussFilter->SetVariance(sqrt(sigma));
    gaussFilter->Update();
    OutputImageType::Pointer processedImage = gaussFilter->GetOutput();
    
    std::ifstream inputfile;
    inputfile.open(graphName.c_str(), std::ios::in);
    GraphType graph;
    inputfile >> graph;
    inputfile.close();
    
    RealType maxCenterlineIntensity = 0.0;
    RealType minCenterlineIntensity = 1.0;
    GraphType::ConstIterator segmentIter = graph.Begin();
    for(; segmentIter != graph.End(); ++segmentIter)
    {
        GraphType::SegmentType currentSegment = segmentIter->second;
        GraphType::SegmentType::ConstIterator nodeIter = currentSegment.Begin();
        for(; nodeIter != currentSegment.End(); ++nodeIter)
        {
            GraphType::NodeIDType currentNodeID = *nodeIter;
            GraphType::NodeType currentNode = graph.GetNode(currentNodeID);
            GraphType::NodeType::PositionType currentPosition = currentNode.GetPosition();
            
            InputImageType::PointType currentPoint;
            for(unsigned int i = 0; i < Dimension; i++)
            {
                currentPoint[i] = currentPosition[i];
            }
            InputImageType::IndexType currentIndex;
            processedImage->TransformPhysicalPointToIndex(currentPoint, currentIndex);
            RealType currentIntensity = processedImage->GetPixel(currentIndex);
            if(currentIntensity > maxCenterlineIntensity)
            {
                maxCenterlineIntensity = currentIntensity;
            }
            if(currentIntensity < minCenterlineIntensity && currentIntensity != 0)
            {
                minCenterlineIntensity = currentIntensity;
            }
            
        }
    }
    
    std::cout << "Max Intensity " << maxCenterlineIntensity << std::endl;
    std::cout << "Min Intensity " << minCenterlineIntensity << std::endl;
    
    itk::ImageRegionIterator< OutputImageType > iter(processedImage, processedImage->GetRequestedRegion() );
    for(iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
    {
        RealType currentValue = iter.Get();
        RealType newValue = 0;
        if(currentValue <= minCenterlineIntensity)
        {
            newValue = backgroundValue;
        }
        else if(currentValue > minCenterlineIntensity)
        {
            newValue = ((maxValue - minValue) / (maxCenterlineIntensity - minCenterlineIntensity)) * (currentValue-minCenterlineIntensity) + minValue;
        }
        RealType currentNoise = BoxMuller(noiseMean, noiseSTD);
        iter.Set(newValue + currentNoise);
    }
    

    //Set up writers
    typedef itk::ImageFileWriter< OutputImageType >  ImageWriterType;
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    
    std::string imageName = outputName;
    imageWriter->SetFileName( imageName );
    imageWriter->SetInput( processedImage);
    imageWriter->Update();
    
    return(0);
}
