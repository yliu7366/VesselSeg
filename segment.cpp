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
typedef double InputPixelType;
typedef itk::Image< InputPixelType, Dimension > InputImageType;
typedef double OutputPixelType;
typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
typedef itk::Image< LabelType, Dimension > LabelImageType;

typedef itk::MinimumMaximumImageCalculator< InputImageType > MinMaxType;

typedef itk::MedialVesselFilter< InputImageType, OutputImageType > VesselFilterType;
        
typedef itk::HessianVesselFilter< InputImageType, OutputImageType > HessianVesselFilterType;
        
typedef itk::BinaryThinningImageFilter3D< LabelImageType, LabelImageType > ThinningFilterType;

typedef itk::Vector< RealType, Dimension > VectorType;
typedef itk::Image< VectorType, Dimension > ImageVectorVolumeType;    
typedef itk::Medialness3DFunction< ImageVectorVolumeType, RealType > MedialFunctionType; 
typedef itk::MedialnessMin3DFunction< ImageVectorVolumeType, RealType > MedialMinFunctionType; 
typedef itk::MedialnessMin3DFunction< ImageVectorVolumeType, RealType > MedialMinFunctionType; 
typedef itk::MedialnessVarRadius3DFunction< ImageVectorVolumeType, RealType > MedialVarRadFunctionType; 

typedef itk::GradientRecursiveGaussianImageFilter< InputImageType, ImageVectorVolumeType > GradientFilterType;

typedef itk::VotingBinaryIterativeHoleFillingImageFilter< VesselFilterType::LabelImageType > HoleFilterType;

typedef CenterlineGraphConverter<LabelType, RealType, Dimension> ConverterType;

typedef itk::LineTransform3DFunction<InputImageType, RealType> LineFunctionType;
typedef itk::VectorLineTransform3DFunction<ImageVectorVolumeType, RealType> VectorLineFunctionType;
typedef itk::CorrelateVectorsOnLine3DFunction<ImageVectorVolumeType, RealType> CorrelateVectorsOnLineType;


int main(int argc, char ** argv)
{
    std::cout << "Vessel Segmenter" << std::endl;
    if(argc != 20)
    {
        std::cout << "Wrong Number of Parameters" << std::endl;
        return 0;
    }

    std::string inputName = std::string(argv[1]);
    std::string outputName = std::string(argv[2]);
    RealType startRadius = atof(argv[3]);
    RealType stopRadius =  atof(argv[4]);
    int numRadii = atoi(argv[5]);
    RealType intensityMin = atof(argv[6]);
    RealType intensityMax = atof(argv[7]);
    RealType intensityDifference = atof(argv[8]);
    RealType gradientSigma = atof(argv[9]);
    RealType edgeRatio = atof(argv[10]);
    RealType vesselnessAlpha = atof(argv[11]);
    RealType vesselnessBeta = atof(argv[12]);
    RealType vesselnessGamma = atof(argv[13]);
    RealType vesselnessThreshold = atof(argv[14]);
    RealType medialnessThreshold = atof(argv[15]);
    RealType medialnessPercent = atof(argv[16]);
    int maxNumPredictions = atoi(argv[17]);
    RealType dotProductThreshold = atof(argv[18]);
    int unconnectedSegments = atoi(argv[19]);



    //Set up reader and read input
    typedef itk::ImageFileReader< InputImageType > ImageReaderType;
    ImageReaderType::Pointer   reader = ImageReaderType::New();
    reader->SetFileName( inputName ); 
    reader->Update();
    InputImageType::Pointer inputImage = reader->GetOutput();

    //Set up writers
    typedef itk::ImageFileWriter< OutputImageType >  OutputImageWriterType;
    OutputImageWriterType::Pointer outputWriter = OutputImageWriterType::New();
       
    typedef itk::ImageFileWriter< LabelImageType >  LabelImageWriterType;
    LabelImageWriterType::Pointer labelWriter = LabelImageWriterType::New();
    
    typedef itk::ImageFileWriter< ImageVectorVolumeType >  VectorImageWriterType;
    VectorImageWriterType::Pointer vectorWriter = VectorImageWriterType::New();
    
    std::string inputName2 = std::string(argv[2]);
    inputName2.insert(inputName2.find_first_of("."),"Input");
    outputWriter->SetFileName( inputName2.c_str() );
    outputWriter->SetInput( inputImage );
    outputWriter->Update();
    
    InputImageType::SpacingType pixelSpacing = inputImage->GetSpacing();
    typedef std::vector<RealType> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(pixelSpacing[0]);
    spacingList.push_back(pixelSpacing[1]);
    spacingList.push_back(pixelSpacing[2]);
    
    RealType minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<RealType>() );
    RealType maxSpacing = *std::max_element(spacingList.begin(), spacingList.end(), std::less_equal<RealType>() );


    //timing parameters
    clock_t t1;
    clock_t t2;
    
    //Create filter

    std::vector<RealType> radiiList;
    std::vector<RealType> scaleList;
    std::cout << "Scales " << std::endl;
    for(unsigned int i = 0; i <= numRadii; i++)
    {
        RealType currentRadius = startRadius + RealType(i)/RealType(numRadii)*(stopRadius - startRadius);
        RealType currentScale = exp(log(startRadius) + (RealType(i)/RealType(numRadii)) * 
            (log(stopRadius) - log(startRadius)))*1/*/sqrt(3.0)*/;
        radiiList.push_back( currentRadius );
        scaleList.push_back(minSpacing+currentScale);
        std::cout << currentScale << " ";
    }
    
    std::cout << std::endl;
    
    ThinningFilterType::Pointer thinningFilter = ThinningFilterType::New();
    
    HoleFilterType::Pointer holeFilter = HoleFilterType::New();
    holeFilter->SetBackgroundValue(   0 );
    holeFilter->SetForegroundValue( 1 );
    holeFilter->SetMaximumNumberOfIterations( 20 );
    holeFilter->SetMajorityThreshold( 0 );
    
    MinMaxType::Pointer minMaxCalc = MinMaxType::New();
    minMaxCalc->SetImage(inputImage);
    minMaxCalc->Compute();
    std::cout << "Image Max " << minMaxCalc->GetMaximum() << " Min " << minMaxCalc->GetMinimum() << std::endl;

   
    std::cout << "Medialness Filter " << std::endl;
    
    //RealType gradientSigma = minSpacing/2.0;
    //Medialness Function
    GradientFilterType::Pointer gradientFilter = GradientFilterType::New();
    gradientFilter->SetInput( inputImage );
    gradientFilter->SetNormalizeAcrossScale(false);
    gradientFilter->SetSigma( gradientSigma );
    gradientFilter->SetUseImageDirection(false);
      
    std::cout << "Gradient computation" << std::endl;
    std::cout.flush();
    t1 = clock(); 
    gradientFilter->Update();
    ImageVectorVolumeType::Pointer gradientVolume = gradientFilter->GetOutput();
    t2 = clock(); 
    std::cout << "Gradient time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
    
/*    ImageVectorVolumeType::Pointer gradientImage = gradientFilter->GetOutput();
    std::string gradientName = std::string(argv[2]);
    gradientName.insert(gradientName.find_first_of("."),"Gradient");
    vectorWriter->SetFileName( gradientName.c_str() );
    vectorWriter->SetInput( gradientImage );
    vectorWriter->Update();
*/    
    typedef itk::VectorLinearInterpolateImageFunction<ImageVectorVolumeType> ImageVectorInterpolaterType;
    ImageVectorInterpolaterType::Pointer gradientInterpolator = ImageVectorInterpolaterType::New();
    gradientInterpolator->SetInputImage(gradientVolume);
    
    MedialFunctionType::Pointer medialFunction = MedialFunctionType::New();
    medialFunction->SetInputImage(gradientVolume); 
    medialFunction->SetInterpolater(&(*gradientInterpolator));
    medialFunction->SetRadii(startRadius, stopRadius, numRadii, 100);
    
    MedialVarRadFunctionType::Pointer medialVarRadFunction = MedialVarRadFunctionType::New();
    medialVarRadFunction->SetInputImage(gradientVolume); 
    medialVarRadFunction->SetInterpolater(&(*gradientInterpolator));
    medialVarRadFunction->SetRadii(startRadius, stopRadius, numRadii, 10);
    
    //Hessian
    HessianVesselFilterType::Pointer hessianVesselFilter = HessianVesselFilterType::New();
    hessianVesselFilter->SetInput(inputImage);
    hessianVesselFilter->SetVesselnessAlpha(vesselnessAlpha);
    hessianVesselFilter->SetVesselnessBeta(vesselnessBeta);
    hessianVesselFilter->SetVesselnessGamma(vesselnessGamma);
    hessianVesselFilter->SetScales(scaleList);
    hessianVesselFilter->SetIntensityMaximum( intensityMax );
    hessianVesselFilter->SetIntensityMinimum( intensityMin );
    hessianVesselFilter->SetIntensityDifference( intensityDifference );
    hessianVesselFilter->SetEdgeRatio( edgeRatio );
    hessianVesselFilter->SetMedialFunction(&*(medialFunction));
    //hessianVesselFilter->SetMedialFunction(&*(medialVarRadFunction));
    hessianVesselFilter->Allocate();

    std::cout << "Vesselness computation "<< std::endl;
    t1 = clock(); 
    hessianVesselFilter->ComputeVesselness();
    hessianVesselFilter->ComputeMedialVesselness();
    t2 = clock(); 
    std::cout << "Vesselness time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
    OutputImageType::Pointer vesselnessImage = hessianVesselFilter->GetVesselnessImage();
    std::string vesselnessName = std::string(argv[2]);
    vesselnessName.insert(vesselnessName.find_first_of("."),"Vesselness");
    outputWriter->SetFileName( vesselnessName.c_str() );
    outputWriter->SetInput( vesselnessImage );
    outputWriter->Update();
    
    OutputImageType::Pointer vesselMedialImage = hessianVesselFilter->GetVesselMedialImage();
    std::string vesselMedialName = std::string(argv[2]);
    vesselMedialName.insert(vesselMedialName.find_first_of("."),"VesselMedial");
    outputWriter->SetFileName( vesselMedialName.c_str() );
    outputWriter->SetInput( vesselMedialImage );
    outputWriter->Update();
 
/*    
    OutputImageType::Pointer scaleImage = hessianVesselFilter->GetScaleImage();
    std::string scaleName = std::string(argv[2]);
    scaleName.insert(scaleName.find_first_of("."),"Scale");
    outputWriter->SetFileName( scaleName.c_str() );
    outputWriter->SetInput( scaleImage );
    outputWriter->Update();
*/    
    
    
    //Vessel Segmenter
    
    MedialMinFunctionType::Pointer medialMinFunction = MedialMinFunctionType::New();
    medialMinFunction->SetInputImage(gradientVolume); 
    medialMinFunction->SetInterpolater(&(*gradientInterpolator));
    medialMinFunction->SetRadii(startRadius, stopRadius, numRadii, 100);
    

    
    
     
    VesselFilterType::Pointer vesselFilter = VesselFilterType::New();
    vesselFilter->SetInput(inputImage);
    hessianVesselFilter->Register();
    vesselFilter->SetVesselFilter(&(*hessianVesselFilter));
    medialFunction->Register();
    vesselFilter->SetMedialFunction(&*(medialMinFunction));
    //vesselFilter->SetMedialFunction(&*(medialVarRadFunction));
    vesselFilter->Allocate();
    vesselFilter->SetVesselnessThreshold(vesselnessThreshold);
    vesselFilter->SetMedialnessThreshold(medialnessThreshold);
    vesselFilter->SetMaxNumMedialPrediction(maxNumPredictions);
    vesselFilter->SetMedialnessPercent(medialnessPercent);
    vesselFilter->SetDotProductThreshold(dotProductThreshold);
    vesselFilter->SetIntensityMaximum( intensityMax );
    vesselFilter->SetIntensityMinimum( intensityMin );
    vesselFilter->SetIntensityDifference( intensityDifference );
    vesselFilter->SetEdgeRatio( edgeRatio );

    
    std::cout << "Medialness computation" << std::endl;
    std::cout.flush();
    t1 = clock(); 
    vesselFilter->ComputeMedialness();
    t2 = clock(); 
    std::cout << "Medialness time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
    
    OutputImageType::Pointer medialnessImage = vesselFilter->GetMedialnessImage();
    std::string medialnessName = std::string(argv[2]);
    medialnessName.insert(medialnessName.find_first_of("."),"Medialness");
    outputWriter->SetFileName( medialnessName.c_str() );
    outputWriter->SetInput( medialnessImage );
    outputWriter->Update();
    
//    return 0;
    
    std::cout << "Initial Centerline computation" << std::endl;
    std::cout.flush();
    t1 = clock(); 
    vesselFilter->ComputeCenterlineFromMedialness();
    t2 = clock(); 
    std::cout << "Initial Centerline time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
    VesselFilterType::LabelImageType::Pointer initialCenterlineImage = vesselFilter->GetCenterlineImage();
    
    std::string centerlineInitialName = std::string(argv[2]);
    centerlineInitialName.insert(centerlineInitialName.find_first_of("."),"InitialCenterline");
    labelWriter->SetFileName( centerlineInitialName.c_str() );
    labelWriter->SetInput( initialCenterlineImage );
    labelWriter->Update();
      
    std::cout << "Centerline Prediction computation" << std::endl;
    std::cout.flush();
    t1 = clock(); 
    vesselFilter->PredictCenterlineFromMedialness();
    t2 = clock(); 
    std::cout << "Centerline Prediction time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
    VesselFilterType::LabelImageType::Pointer predictCenterlineImage = vesselFilter->GetCenterlineImage();
/*    
    std::string centerlinePredicitedName = std::string(argv[2]);
    centerlinePredicitedName.insert(centerlinePredicitedName.find_first_of("."),"PredictedCenterline");
    labelWriter->SetFileName( centerlinePredicitedName.c_str() );
    labelWriter->SetInput( predictCenterlineImage );
    labelWriter->Update();
*/    
    std::cout << "Link Planes" << std::endl;
    std::cout.flush();
    t1 = clock(); 
    vesselFilter->LinkMedialnessPlanes();
    t2 = clock(); 
    std::cout << "Link Planes time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
    LabelImageType::Pointer linkPlanesImage = vesselFilter->GetCenterlineImage();
    LabelImageType::Pointer centerlineBinaryImage = vesselFilter->GetCenterlineAsBinary();
    LabelImageType::Pointer centerlineSphereBinaryImage = vesselFilter->GetCenterlineSphereImage();

    std::string centerlineLabelsName = std::string(argv[2]);
    centerlineLabelsName.insert(centerlineLabelsName.find_first_of("."),"CenterlineBinaryBeforeThinning");
    labelWriter->SetFileName( centerlineLabelsName.c_str() );
    labelWriter->SetInput( centerlineBinaryImage );
    labelWriter->Update();
    
    std::string centerlineSphereLabelsName = std::string(argv[2]);
    centerlineSphereLabelsName.insert(centerlineSphereLabelsName.find_first_of("."),"CenterlineSphereBeforeThinning");
    labelWriter->SetFileName( centerlineSphereLabelsName.c_str() );
    labelWriter->SetInput( centerlineSphereBinaryImage );
    labelWriter->Update();
    
    OutputImageType::Pointer radiusImage = vesselFilter->GetRadiusImage();
    std::string radiusName = std::string(argv[2]);
    radiusName.insert(radiusName.find_first_of("."),"Radius");
    outputWriter->SetFileName( radiusName.c_str() );
    outputWriter->SetInput( radiusImage );
    outputWriter->Update();
    
    OutputImageType::Pointer lineCenterImage = vesselFilter->GetLineCenterImage();
    std::string lineCenterName = std::string(argv[2]);
    lineCenterName.insert(lineCenterName.find_first_of("."),"LineCenter");
    outputWriter->SetFileName( lineCenterName.c_str() );
    outputWriter->SetInput( lineCenterImage );
    outputWriter->Update();
    
    OutputImageType::Pointer linePredictImage = vesselFilter->GetLinePredictImage();
    std::string linePredictName = std::string(argv[2]);
    linePredictName.insert(linePredictName.find_first_of("."),"LinePredict");
    outputWriter->SetFileName( linePredictName.c_str() );
    outputWriter->SetInput( linePredictImage );
    outputWriter->Update();
    
    
    std::cout << "Thinning Algorithm On Labeled Centerline" << std::endl;
    std::cout.flush();
    t1 = clock();
    LabelImageType::Pointer centerlineLabeledImage = vesselFilter->LabelSphereOnCenterline(centerlineBinaryImage, 1, radiusImage, 1);
    
    holeFilter->SetInput(centerlineLabeledImage);
    holeFilter->Update();
    LabelImageType::Pointer centerlineLabeledHoleImage = holeFilter->GetOutput();

    thinningFilter->SetInput(centerlineLabeledHoleImage);
    thinningFilter->Update();
    LabelImageType::Pointer thinnedCenterlineLabeledImage = thinningFilter->GetOutput();
    
    t2 = clock(); 
    std::cout << "Thinning Algorithm On Labeled Centerline time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
   
    
    std::string centerlineLabeledName = std::string(argv[2]);
    centerlineLabeledName.insert(centerlineLabeledName.find_first_of("."),"CenterlineSphereLabeled");
    labelWriter->SetFileName( centerlineLabeledName.c_str() );
    labelWriter->SetInput( centerlineLabeledImage );
    labelWriter->Update();
    
    std::string centerlineLabeledHoleName = std::string(argv[2]);
    centerlineLabeledHoleName.insert(centerlineLabeledHoleName.find_first_of("."),"CenterlineLabeledHoleFilled");
    labelWriter->SetFileName( centerlineLabeledHoleName.c_str() );
    labelWriter->SetInput( centerlineLabeledHoleImage );
    labelWriter->Update();
    
    std::string thinCenterlineLabeledName = std::string(argv[2]);
    thinCenterlineLabeledName.insert(thinCenterlineLabeledName.find_first_of("."),"CenterlineLabeledThinned");
    labelWriter->SetFileName( thinCenterlineLabeledName.c_str() );
    labelWriter->SetInput( thinnedCenterlineLabeledImage );
    labelWriter->Update();


    std::cout << "Build Vessel Network" << std::endl;
    std::cout.flush();
    t1 = clock();
    ConverterType converter;
    ConverterType::VesselGraphType vesselGraph = converter.ImageToGraph(thinnedCenterlineLabeledImage, radiusImage, medialnessImage, hessianVesselFilter->GetVesselVectorImage());
    vesselGraph.FixPoorRadii(medialnessThreshold);
    vesselGraph.RemoveUnconnectedSegments(unconnectedSegments);  
    vesselGraph = vesselGraph.RemoveRedundentJunctions();
    t2 = clock(); 
    std::cout << "Build Vessel Network time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
 
    std::ofstream outputFile;
    std::string outputFileName(outputName.c_str(), outputName.find_first_of("."));
    outputFileName.append("MedialCenterline.vessel");
    outputFile.open( outputFileName.c_str(), std::ios::out);
    outputFile << vesselGraph;
    outputFile.close();
    
    std::ofstream outputGMLFile;
    std::string graphGMLName(outputName.c_str(), outputName.find_first_of("."));
    graphGMLName.append(".gml");
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
    std::cout << "Build Interpolated Network time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
 
    std::string outputFileSubGraphName(outputName.c_str(), outputName.find_first_of("."));
    outputFileSubGraphName.append("InterpGraph.vessel");
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
    std::cout << "Build Label Map from Graph time "<< (t2-t1)/CLOCKS_PER_SEC << std::endl;
       
    std::string fromGraphName = std::string(argv[2]);
    fromGraphName.insert(fromGraphName.find_first_of("."),"FromGraph");
    outputWriter->SetFileName( fromGraphName.c_str() );
    outputWriter->SetInput( fromGraphImage );
    outputWriter->Update();

   
    return 0;
}
