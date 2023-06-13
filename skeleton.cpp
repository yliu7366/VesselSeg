#include <time.h>
#include <iostream>
#include <string>


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThinningImageFilter.h"
#include "itkBinaryThinningImageFilter3D.h"
#include "itkBinaryThresholdImageFilter.h"


const unsigned int TwoDimension = 2;
const unsigned int ThreeDimension = 3;

typedef unsigned char InputPixelType;
typedef itk::Image< InputPixelType, TwoDimension > InputSliceType;

typedef itk::Image< InputPixelType, ThreeDimension > InputVolumeType;

typedef itk::ShiftScaleImageFilter< InputSliceType, InputSliceType > ShiftScaleFilterType;
typedef itk::RescaleIntensityImageFilter< InputSliceType, InputSliceType > RescaleFilterType;
typedef itk::BinaryThresholdImageFilter< InputSliceType, InputSliceType > ThresholdFilterType;
typedef itk::BinaryThinningImageFilter< InputSliceType, InputSliceType > Thinning2FilterType;
typedef itk::BinaryThinningImageFilter3D< InputVolumeType, InputVolumeType > Thinning3FilterType;


int main(int argc, char ** argv)
{
    if(argc != 3)
    {
        std::cout << "Wrong Number of Parameters" << std::endl;
        return 0;
    }

    std::string inputName = std::string(argv[1]);
    std::string outputName = std::string(argv[2]);

    clock_t t1;
    clock_t t2;
    
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
    
    RescaleFilterType::Pointer rescaledFilter = RescaleFilterType::New();
    rescaledFilter->SetOutputMinimum(0);
    rescaledFilter->SetOutputMaximum(255);
    
    /*
    ShiftScaleFilterType::Pointer shiftScaleFilter = ShiftScaleFilterType::New();
    shiftScaleFilter->SetInput(inputImage);
    shiftScaleFilter->SetShift(-255.0);
    shiftScaleFilter->SetScale(-255.0);
    shiftScaleFilter->Update();
    InputSliceType::Pointer shiftScaleImage = shiftScaleFilter->GetOutput();
    
    inputImage = shiftScaleImage;
    */
    
    ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
    thresholdFilter->SetInsideValue(1);
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->SetLowerThreshold(128);
    thresholdFilter->SetUpperThreshold(255);
    thresholdFilter->SetInput(inputImage);
    thresholdFilter->Update();
    InputSliceType::Pointer thresholdImage = thresholdFilter->GetOutput();

    inputImage = thresholdImage;
    
    std::string shiftScaleName = std::string(argv[2]);
    shiftScaleName.insert(shiftScaleName.find_first_of("."),"Rescaled");
    inputSliceWriter->SetFileName( shiftScaleName.c_str() );
    rescaledFilter->SetInput(inputImage);
    rescaledFilter->Update();
    inputSliceWriter->SetInput( rescaledFilter->GetOutput() );
    inputSliceWriter->Update();
    
    t1 = clock(); 
    Thinning2FilterType::Pointer thinning2Filter = Thinning2FilterType::New();
    thinning2Filter->SetInput(inputImage);
    thinning2Filter->Update();
    InputSliceType::Pointer thinned2Image = thinning2Filter->GetOutput();
    rescaledFilter->SetInput(thinned2Image);
    rescaledFilter->Update();
    t2 = clock(); 
    std::cout << "Thinning 2d time "<< double(t2-t1)/CLOCKS_PER_SEC << std::endl;
    
    std::string thinName = std::string(argv[2]);
    thinName.insert(thinName.find_first_of("."),"Thinned2d");
    inputSliceWriter->SetFileName( thinName.c_str() );
    inputSliceWriter->SetInput( rescaledFilter->GetOutput() );
    inputSliceWriter->Update();
    
    t1 = clock();
    InputVolumeType::Pointer volume = InputVolumeType::New();
    InputVolumeType::SpacingType volumeSpacing;
    volumeSpacing[0] = inputImage->GetSpacing()[0];
    volumeSpacing[1] = inputImage->GetSpacing()[1];
    volumeSpacing[2] = 1.0;
    
    InputVolumeType::PointType volumeOrigin;
    volumeOrigin[0] = inputImage->GetOrigin()[0];
    volumeOrigin[1] = inputImage->GetOrigin()[1];
    volumeOrigin[2] = 0.0;
    
    InputVolumeType::RegionType volumeRegion;
    InputVolumeType::RegionType::SizeType volumeSize;
    volumeSize[0] = inputImage->GetLargestPossibleRegion().GetSize()[0];
    volumeSize[1] = inputImage->GetLargestPossibleRegion().GetSize()[1];
    volumeSize[2] = 1;
    volumeRegion.SetSize(volumeSize);
    
    InputVolumeType::IndexType volumeStart;
    volumeStart[0] = inputImage->GetLargestPossibleRegion().GetIndex()[0];
    volumeStart[1] = inputImage->GetLargestPossibleRegion().GetIndex()[1];
    volumeStart[2] = 0; 
    volumeRegion.SetIndex(volumeStart);
    
    volume->SetRegions(volumeRegion);
    volume->SetRequestedRegionToLargestPossibleRegion ();
    volume->Allocate();
    volume->FillBuffer(0);
    
    itk::ImageRegionIterator<InputSliceType> twoDIter( inputImage, inputImage->GetLargestPossibleRegion() );
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
    centerlineSlice->SetSpacing( inputImage->GetSpacing() );
    centerlineSlice->SetRegions( inputImage->GetLargestPossibleRegion() );
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
    
    rescaledFilter->SetInput(centerlineSlice);
    rescaledFilter->Update();
    t2 = clock(); 
    std::cout << "Thinning 3d time "<< double(t2-t1)/CLOCKS_PER_SEC << std::endl;
    
    std::string centerSliceName = std::string(argv[2]);
    centerSliceName.insert(centerSliceName.find_first_of("."),"Thinned3d");
    inputSliceWriter->SetFileName( centerSliceName.c_str() );
    inputSliceWriter->SetInput( rescaledFilter->GetOutput() );
    inputSliceWriter->Update();
    
    
    
    return 0;
}
