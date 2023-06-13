#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>

#include "itkMedialVesselFilter.h"

#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkSphereSpatialFunction.h"

//foo
#include "itkLineTransform3DFunction.h"


namespace itk
{

template <class TInputImage, class TOutputImage>
const typename MedialVesselFilter<TInputImage,TOutputImage>::LabelType 
MedialVesselFilter<TInputImage,TOutputImage>::s_MedialCenterlineUnlinked;

template <class TInputImage, class TOutputImage>
const typename MedialVesselFilter<TInputImage,TOutputImage>::LabelType 
MedialVesselFilter<TInputImage,TOutputImage>::s_MedialSphereUnlinked;

template <class TInputImage, class TOutputImage>
const typename MedialVesselFilter<TInputImage,TOutputImage>::LabelType 
MedialVesselFilter<TInputImage,TOutputImage>::s_MedialCenterlineLinked;

template <class TInputImage, class TOutputImage>
const typename MedialVesselFilter<TInputImage,TOutputImage>::LabelType 
MedialVesselFilter<TInputImage,TOutputImage>::s_MedialSphereLinked;

template <class TInputImage, class TOutputImage>
const typename MedialVesselFilter<TInputImage,TOutputImage>::LabelType 
MedialVesselFilter<TInputImage,TOutputImage>::s_MedialCenterlinePredicted;

template <class TInputImage, class TOutputImage>
const typename MedialVesselFilter<TInputImage,TOutputImage>::LabelType 
MedialVesselFilter<TInputImage,TOutputImage>::s_MedialSpherePredicted;


//Constructor
template <class TInputImage, class TOutputImage>
MedialVesselFilter<TInputImage,TOutputImage>::MedialVesselFilter()
{
    
}

//Destructor
template <class TInputImage, class TOutputImage>
MedialVesselFilter<TInputImage,TOutputImage>::~MedialVesselFilter()
{
}

//GenerateData
template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage, TOutputImage>::GenerateData()
{
    this->Allocate();
    this->ComputeMedialness();
    this->ComputeCenterlineFromMedialness();
    this->PredictCenterlineFromMedialness();
//    this->GraftOutput(?);
//    this->GraftNthOutput(0, ?);
}



template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage, TOutputImage>::Allocate()
{
    typename TInputImage::ConstPointer inputImage = this->GetInput();

    //Scalar images
    
    this->m_MedialnessImage = RealScalarImageType::New();
    this->m_MedialnessImage->CopyInformation(inputImage);
    this->m_MedialnessImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_MedialnessImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_MedialnessImage->Allocate(); 
    this->m_MedialnessImage->FillBuffer(0.0);
      
    this->m_RadiusImage = RealScalarImageType::New();
    this->m_RadiusImage->CopyInformation(inputImage);
    this->m_RadiusImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_RadiusImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_RadiusImage->Allocate(); 
    this->m_RadiusImage->FillBuffer(0.0);
    
    //label images
    
    this->m_CenterlineImage = LabelImageType::New();
    this->m_CenterlineImage->CopyInformation(inputImage);
    this->m_CenterlineImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_CenterlineImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_CenterlineImage->Allocate(); 
    this->m_CenterlineImage->FillBuffer(0);
    
    this->m_CenterlineSphereImage = LabelImageType::New();
    this->m_CenterlineSphereImage->CopyInformation(inputImage);
    this->m_CenterlineSphereImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_CenterlineSphereImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_CenterlineSphereImage->Allocate(); 
    this->m_CenterlineSphereImage->FillBuffer(0);
}


template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage, TOutputImage>::
LabelSphereOnCenterlineAtIndex(const IndexType& centerIndex,
    typename LabelImageType::Pointer labelImage, LabelType centerlineLabel, 
    RealType radius, LabelType radiusLabel) const
{ 
    typename TInputImage::SpacingType imageSpacing = this->m_CenterlineImage->GetSpacing();
    
    typename TInputImage::PointType centerPoint;
    this->m_CenterlineImage->TransformIndexToPhysicalPoint(centerIndex, centerPoint);

    int oversamplingFactor = 3;

    int xRadiusOffset = oversamplingFactor * (int(ceil(radius/imageSpacing[0])) + 1);
    int yRadiusOffset = oversamplingFactor * (int(ceil(radius/imageSpacing[1])) + 1);
    int zRadiusOffset = oversamplingFactor * (int(ceil(radius/imageSpacing[2])) + 1);
    
    typedef std::vector<typename TInputImage::SpacingValueType> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(imageSpacing[0]);
    spacingList.push_back(imageSpacing[1]);
    spacingList.push_back(imageSpacing[2]);
    
    typename TInputImage::SpacingValueType minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<typename TInputImage::SpacingValueType>() );
    
    typename SphereSpatialFunction<3, typename TInputImage::PointType>::Pointer currentSphere = SphereSpatialFunction<3, typename TInputImage::PointType>::New();
    currentSphere->SetCenter(centerPoint);
    
    if(radius > 0)
    {
        currentSphere->SetRadius(radius);
    }
    else
    {
        currentSphere->SetRadius(0.25*minSpacing);
    }

    labelImage->SetPixel(centerIndex, centerlineLabel);

    for(int i = -xRadiusOffset; i <= xRadiusOffset; i++)
    {
        for(int j = -yRadiusOffset; j <= yRadiusOffset; j++)
        {
            for(int k = -zRadiusOffset; k <= zRadiusOffset; k++)
            {
                typename TInputImage::PointType movingPoint;
                movingPoint[0] = centerPoint[0] + (RealType(i)/RealType(oversamplingFactor))*imageSpacing[0];
                movingPoint[1] = centerPoint[1] + (RealType(j)/RealType(oversamplingFactor))*imageSpacing[1];
                movingPoint[2] = centerPoint[2] + (RealType(k)/RealType(oversamplingFactor))*imageSpacing[2];
            
                bool insideSphere = currentSphere->Evaluate(movingPoint);
                IndexType movingIndex;
                bool insideImage = this->m_CenterlineImage->TransformPhysicalPointToIndex(movingPoint, movingIndex);
                
                if(insideSphere == true && insideImage == true
                && labelImage->GetPixel(movingIndex) != centerlineLabel)
                {
                    labelImage->SetPixel(movingIndex, radiusLabel);
                }
            }
        }
    }
    
}    

template <class TInputImage, class TOutputImage>
typename MedialVesselFilter<TInputImage, TOutputImage>::LabelImageType::Pointer 
MedialVesselFilter<TInputImage, TOutputImage>::LabelSphereOnCenterline(
    typename LabelImageType::Pointer centerlineImage, LabelType centerlineLabel, 
    typename RealScalarImageType::Pointer radiusImage, LabelType radiusLabel) const
{
    typename LabelImageType::Pointer sphereBinaryImage = LabelImageType::New();
    sphereBinaryImage->SetSpacing( centerlineImage->GetSpacing() );
    sphereBinaryImage->SetOrigin( centerlineImage->GetOrigin() );
    sphereBinaryImage->SetRegions( centerlineImage->GetLargestPossibleRegion() );
    sphereBinaryImage->Allocate(); 
    sphereBinaryImage->FillBuffer(0);
    
    itk::ImageRegionConstIterator<LabelImageType> centerlineIter(centerlineImage, 
        centerlineImage->GetRequestedRegion() );  

    for(centerlineIter.GoToBegin(); !centerlineIter.IsAtEnd(); ++centerlineIter)
    {
        LabelType isCenterline = centerlineIter.Get();
        IndexType currentIndex = centerlineIter.GetIndex();
        if(isCenterline == centerlineLabel )
        {
            RealType currentRadius = radiusImage->GetPixel(currentIndex);
            this->LabelSphereOnCenterlineAtIndex(currentIndex, sphereBinaryImage, centerlineLabel, currentRadius, radiusLabel);
        }
    }

    return(sphereBinaryImage);
}



template <class TInputImage, class TOutputImage>
typename MedialVesselFilter<TInputImage, TOutputImage>::LabelImageType::Pointer 
MedialVesselFilter<TInputImage, TOutputImage>::GetCenterlineAsBinary() const
{
    typename LabelImageType::Pointer centerlineBinaryImage = LabelImageType::New();
    centerlineBinaryImage->SetSpacing( this->m_CenterlineImage->GetSpacing() );
    centerlineBinaryImage->SetOrigin( this->m_CenterlineImage->GetOrigin() );
    centerlineBinaryImage->SetRegions( this->m_CenterlineImage->GetLargestPossibleRegion() );
    centerlineBinaryImage->Allocate(); 
    centerlineBinaryImage->FillBuffer(0);
    itk::ImageRegionIterator<LabelImageType> centerlineIter( this->m_CenterlineImage, this->m_CenterlineImage->GetLargestPossibleRegion() );
    
    for(centerlineIter.GoToBegin(); !centerlineIter.IsAtEnd(); ++centerlineIter)
    {
        LabelType currentLabel = centerlineIter.Get();
        IndexType currentIndex = centerlineIter.GetIndex();
        if(currentLabel != 0 )
        {
            centerlineBinaryImage->SetPixel(currentIndex, 1);
        }
    }
    
    return(centerlineBinaryImage);
}


template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage, TOutputImage>::LinkMedialnessPlanes()
{
    typename itk::ConstNeighborhoodIterator< LabelImageType >::RadiusType 
        neighborhoodRadius;
    unsigned int numNeighborsEachDirection = 1;
    unsigned int numNeighbors = (2*numNeighborsEachDirection+1)*
        (2*numNeighborsEachDirection+1)*(2*numNeighborsEachDirection+1);
    neighborhoodRadius.Fill(numNeighborsEachDirection);
    typename itk::ConstNeighborhoodIterator< LabelImageType > 
        centerlineNeighIter(neighborhoodRadius, this->m_CenterlineImage, 
        this->m_CenterlineImage->GetRequestedRegion() );
    
    for(centerlineNeighIter.GoToBegin(); !centerlineNeighIter.IsAtEnd(); ++centerlineNeighIter)
    {
        LabelType currentCenterline = centerlineNeighIter.GetCenterPixel();
        typename LabelImageType::IndexType currentIndex = 
            centerlineNeighIter.GetIndex();

        //Get neighborhood around index        
        typename itk::Neighborhood<LabelType, 
            itkGetStaticConstMacro(ImageDimension)> centerlineHood = 
            centerlineNeighIter.GetNeighborhood();
            
        //count number of neighbors that are centerlines  
        std::vector<IndexType> neighboorCenterIndexVectorType;
        for (unsigned int i = 0; i < numNeighbors; i++)
        {
            //foo fix with iterator
            LabelType neighCenterline = centerlineHood.GetElement(i);
            typename LabelImageType::IndexType hoodIndex = centerlineNeighIter.GetIndex(i);
            if(neighCenterline != 0)
            {
                neighboorCenterIndexVectorType.push_back( hoodIndex );          
            }
        }
        
        //see if neighbors are algined
        if(neighboorCenterIndexVectorType.size() >= 3 && currentCenterline != 0)
        {
            RealVectorType currentVesselVector = this->hessianVesselFilter->GetVesselVectorImage()->GetPixel(currentIndex);
            unsigned int alignedCounter = 0;
            for(unsigned int j = 0; j < neighboorCenterIndexVectorType.size(); j++)
            {
                IndexType neighIndex = neighboorCenterIndexVectorType[j];
                if(currentIndex != neighIndex)
                {
                    RealVectorType neighVesselVector = this->hessianVesselFilter->GetVesselVectorImage()->GetPixel(neighIndex);
                    RealType dot = fabs(dot_product(currentVesselVector.GetVnlVector(), neighVesselVector.GetVnlVector()));
                    if(dot > this->m_DotProductThreshold)
                    {
                        ++alignedCounter;
                    }
                }
            }
            //paint circle
            if(alignedCounter >= 2)
            {
                this->m_CenterlineImage->SetPixel(currentIndex, s_MedialCenterlineLinked);
            }
        
        }//end neighborhood condition
        
    }//end centerline loop
}


template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage, TOutputImage>::
PredictCenterlineFromMedialness()
{

    typedef LineTransform3DFunction<TInputImage, RealType> LineFunctionType;
    typedef LineTransformMinVariance3DFunction<TInputImage, RealType> LineMinVarFunctionType;
    typedef VectorLineTransform3DFunction<RealVectorImageType, RealType> VectorLineFunctionType;
    typedef CorrelateVectorsOnLine3DFunction<RealVectorImageType, RealType> VectorLineCorrelateFunctionType;

   
    typename RealScalarImageType::Pointer startImage = this->hessianVesselFilter->GetVesselnessImage();
    typename RealScalarImageType::Pointer startRadiusImage = this->hessianVesselFilter->GetScaleImage();

    typename LineFunctionType::Pointer lineFunction = LineFunctionType::New();
    lineFunction->SetInputImage(startImage);

    RealVectorType zeroVector;
    zeroVector.Fill(0.0); 
    typename RealVectorImageType::Pointer lineVectorImage = RealVectorImageType::New();
    lineVectorImage->CopyInformation(this->GetInput());
    lineVectorImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
    lineVectorImage->SetRequestedRegionToLargestPossibleRegion();
    lineVectorImage->Allocate();
    lineVectorImage->FillBuffer(zeroVector); 
    
    
    typename RealScalarImageType::Pointer lineImage = RealScalarImageType::New();
    lineImage->CopyInformation(this->GetInput());
    lineImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
    lineImage->SetRequestedRegionToLargestPossibleRegion();
    lineImage->Allocate();
    lineImage->FillBuffer(0.0); 
    
    this->m_LinePredictImage = RealScalarImageType::New();
    this->m_LinePredictImage->CopyInformation(this->GetInput());
    this->m_LinePredictImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
    this->m_LinePredictImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_LinePredictImage->Allocate();
    this->m_LinePredictImage->FillBuffer(0.0); 
  
    typename itk::ConstNeighborhoodIterator< RealScalarImageType >::RadiusType 
        neighborhoodRadius;
    unsigned int numNeighborsEachDirection = 1;
    unsigned int numNeighbors = (2*numNeighborsEachDirection+1)*
        (2*numNeighborsEachDirection+1)*(2*numNeighborsEachDirection+1);
    neighborhoodRadius.Fill(numNeighborsEachDirection);
    typename itk::ConstNeighborhoodIterator< RealScalarImageType > 
        startIter(neighborhoodRadius, startImage, startImage->GetLargestPossibleRegion() );
        
    for (startIter.GoToBegin(); !startIter.IsAtEnd();++startIter)
    {   
        RealType currentValue = startIter.GetCenterPixel();
        typename RealScalarImageType::IndexType currentIndex = 
            startIter.GetIndex();
            
        typename itk::Neighborhood<RealType, 
            itkGetStaticConstMacro(ImageDimension)> hood = 
            startIter.GetNeighborhood();
        
        //count number of neighbors with vesselness above threshold    
        unsigned int counter = 0;
        for (unsigned int i = 0; i < numNeighbors; i++)
        {
            RealType neighValue = hood.GetElement(i);
            if(neighValue > this->m_VesselnessThreshold * this->m_MedialnessPercent)
            {
                ++counter;              
            }
        }
        
        if(counter >= 1)
        {
            RealType currentRadius = startRadiusImage->GetPixel(currentIndex);
            
            lineFunction->SetParameters(currentRadius/2.0, 10);
            typename LineFunctionType::TOutput lineOutput = lineFunction->EvaluateAtIndex( currentIndex );
            
            lineImage->SetPixel( currentIndex, lineOutput.GetNorm());
            //this->m_TestImage->SetPixel( currentIndex, lineOutput.GetNorm() );
            
            if(lineOutput.GetNorm() > 0.0)
            {
                lineOutput.Normalize();
            }
            
            lineVectorImage->SetPixel( currentIndex, lineOutput);
        }
    }

    typename VectorLineFunctionType::Pointer vectorLineFunction = VectorLineFunctionType::New();
    vectorLineFunction->SetInputImage(lineVectorImage);
 
    itk::ImageRegionConstIterator<RealScalarImageType>  
        lineIter(lineImage, lineImage->GetRequestedRegion() );

    for (lineIter.GoToBegin(); !lineIter.IsAtEnd(); ++lineIter)
    { 
        if(lineIter.Get() >  this->m_VesselnessThreshold * this->m_MedialnessPercent)
        {
            typename RealScalarImageType::IndexType currentIndex = lineIter.GetIndex();
            RealType currentRadius = startRadiusImage->GetPixel(currentIndex);
            
            vectorLineFunction->SetParameters(currentRadius/2.0, 10);
            typename VectorLineFunctionType::TOutput vectorLineOutput = vectorLineFunction->EvaluateAtIndex( currentIndex );
            RealVectorType lineVector = vectorLineOutput.second;
            this->m_LinePredictImage->SetPixel( currentIndex, vectorLineOutput.first * lineImage->GetPixel(currentIndex));
            
        }
    }

 
    typename RealScalarImageType::Pointer centerInputImage = this->m_LinePredictImage;
    
    typename RealVectorImageType::Pointer normalVectorImage = this->hessianVesselFilter->GetVesselVectorImage();

    itk::ImageRegionConstIterator<LabelImageType>  
        centerlineIter(this->m_CenterlineImage, 
        this->m_CenterlineImage->GetRequestedRegion() );
    
        
    for (centerlineIter.GoToBegin(); !centerlineIter.IsAtEnd(); ++centerlineIter)
    {  
        //start from centerline    
        if(centerlineIter.Get() == s_MedialCenterlineUnlinked )
        {
            //Go forward and backward     
            for(unsigned int oneOrTwo = 1; oneOrTwo <= 2; oneOrTwo++)
            {
                typename LabelImageType::IndexType currentIndex = 
                    centerlineIter.GetIndex();                    
                               
                RealVectorType currentVectorVessel = normalVectorImage->GetPixel(currentIndex);
                RealType currentRadius = this->m_RadiusImage->GetPixel(currentIndex);
                
                RealType lastRadius = currentRadius;
                RealVectorType lastVectorVessel = currentVectorVessel;
                                    
                RealType plusOrMinus = pow(-1.0, RealType(oneOrTwo));
                currentVectorVessel = plusOrMinus * currentVectorVessel;
                
                bool keepPredicting = true;
                IndexListType indexList;
                indexList.push_back(currentIndex);
                
                RealType accumulatedValue =  centerInputImage->GetPixel(currentIndex);
                
                //follow vessel vectors    
                for(unsigned int k = 0; k < this->m_MaxNumMedialPrediction && keepPredicting==true; k++)
                {   
                    //to be proven otherwise later
                    keepPredicting = false;       
                    
                    OffsetType windowOffsetTwo;
                    windowOffsetTwo.Fill(-2);
                    typename RealScalarImageType::SizeType windowSizeFive;
                    windowSizeFive.Fill(5);
                    typename RealScalarImageType::RegionType currentRegion(currentIndex + windowOffsetTwo, windowSizeFive);
                    
                    bool isInside = this->m_CenterlineImage->GetLargestPossibleRegion().IsInside(currentRegion);
                    
                    //need to make sure predicted index and its neighborhood are in image
                    if(isInside == true)
                    {
                        /*
                        RealType currentMedial = this->hessianVesselFilter->GetVesselnessImage()->GetPixel(currentIndex);//this->m_MedialnessImage->GetPixel(currentIndex);
                        if(currentMedial == 0.0 )
                        {
                            RealVectorType currentNormal = this->hessianVesselFilter->GetVesselVectorImage()->GetPixel(currentIndex);
                            this->medialFunction->SetNormal(currentNormal);
                            typename MedialFunctionType::TOutput currentMedial = this->medialFunction->EvaluateAtIndex(currentIndex);
                            this->m_MedialnessImage->SetPixel(currentIndex, currentMedial.first);
                            this->m_RadiusImage->SetPixel(currentIndex, currentMedial.second);
                        }
                        */
                        RealType currentValue = centerInputImage->GetPixel(currentIndex);
                    
                        IndexListType candidateIndexList;
                        if(currentValue >= this->m_MedialnessThreshold)
                        {
                            candidateIndexList = this->GetInDirectionIndices(currentIndex, currentVectorVessel, M_PI/4.0);
                        }
                        else
                        {
                            candidateIndexList = this->GetInDirectionIndices(currentIndex, currentVectorVessel, M_PI/3.0);
                        }
                        
                        IndexType predictedIndex;
                        predictedIndex.Fill(0);
                        RealType predictedValue = 0.0;
                        RealType predictedIntensity = 0.0;
                        
                        typename IndexListType::const_iterator candidateIter = candidateIndexList.begin();
                        for(; candidateIter != candidateIndexList.end(); ++candidateIter)
                        {
                            IndexType currentCandidateIndex = *(candidateIter);
                            
                            /*
                            RealType currentCandidateMedial = this->hessianVesselFilter->GetVesselnessImage()->GetPixel(currentCandidateIndex);//this->m_MedialnessImage->GetPixel(currentCandidateIndex);
                            if(currentCandidateMedial == 0.0 )
                            {
                                RealVectorType currentNormal = this->hessianVesselFilter->GetVesselVectorImage()->GetPixel(currentCandidateIndex);
                                this->medialFunction->SetNormal(currentNormal);
                                typename MedialFunctionType::TOutput currentMedial = this->medialFunction->EvaluateAtIndex(currentCandidateIndex);
                                this->m_MedialnessImage->SetPixel(currentCandidateIndex, currentMedial.first);
                                this->m_RadiusImage->SetPixel(currentCandidateIndex, currentMedial.second);
                            }
                            */
                            RealType currentCandidateValue = centerInputImage->GetPixel(currentCandidateIndex);
                            RealType currentCandidateIntensity = this->GetInput()->GetPixel(currentCandidateIndex);
                            if(currentCandidateValue*currentCandidateIntensity > predictedValue*predictedIntensity )
                            {
                                predictedValue = currentCandidateValue;
                                predictedIndex = currentCandidateIndex;
                                predictedIntensity = currentCandidateIntensity;
                            }
                            
                        }
                        
                                                
                        indexList.push_back(predictedIndex);
                        accumulatedValue += predictedValue;
                        LabelType localMaxCenterline = this->m_CenterlineImage->GetPixel(predictedIndex);
                        
                        //intermediate point in predicted segemnt ... keep going
                        if(localMaxCenterline != s_MedialCenterlineUnlinked  && this->m_CenterlineSphereImage->GetPixel(predictedIndex) != s_MedialSphereUnlinked && predictedIndex != currentIndex)
                        {   
                            keepPredicting = true;

                            IndexType previousIndex = currentIndex;
                            currentIndex = predictedIndex;
                            RealVectorType newVectorVessel;
                            
                            if(this->m_MedialnessImage->GetPixel(currentIndex) > this->m_MedialnessThreshold)
                            {
                                newVectorVessel = normalVectorImage->GetPixel(currentIndex);
                                lastRadius = this->m_RadiusImage->GetPixel(currentIndex);
                                lastVectorVessel = newVectorVessel;
                            }
                            else
                            {
                                typename TInputImage::PointType previousPoint;
                                centerInputImage->TransformIndexToPhysicalPoint(previousIndex, previousPoint );
                                
                                typename TInputImage::PointType currentPoint;
                                centerInputImage->TransformIndexToPhysicalPoint(currentIndex, currentPoint );
                                
                                RealVectorType theVector( (currentPoint - previousPoint).GetDataPointer());
                                theVector.Normalize();
                                newVectorVessel = theVector;
                            }
                            
                            
                            RealType dotProduct = dot_product(
                                currentVectorVessel.GetVnlVector(), 
                                newVectorVessel.GetVnlVector());
                            
                            if(dotProduct < 0.0)
                            {
                                currentVectorVessel = -newVectorVessel;
                            }
                            else if(dotProduct > 0.0)
                            {
                                currentVectorVessel = newVectorVessel;
                            }
                            //else new vector is zero and we use the old one

                        }

                    }//end inside image condition
                    
                }//end prediction loop    
                //add )
                //add addition of segments
                if(indexList.size() >=2 &&
                   accumulatedValue/RealType(indexList.size()) >= this->m_MedialnessThreshold * this->m_MedialnessPercent &&
                    this->m_CenterlineImage->GetPixel(indexList.front()) == s_MedialCenterlineUnlinked &&
                    (this->m_CenterlineImage->GetPixel(indexList.back()) == s_MedialCenterlineUnlinked || this->m_CenterlineSphereImage->GetPixel(indexList.back()) == s_MedialSphereUnlinked ) )
                {
                    RealType averageGoodRadius = 0.0;
                    int numAboveAverage = 0;
                    typename IndexListType::const_iterator indexIter = indexList.begin();
                    for(; indexIter != indexList.end(); ++indexIter)
                    {
                        IndexType prospectiveIndex = *(indexIter);
                        RealType prospectiveMedialness = this->m_MedialnessImage->GetPixel(prospectiveIndex);
                        if(prospectiveMedialness >= this->m_MedialnessThreshold )
                        {
                            averageGoodRadius += this->m_RadiusImage->GetPixel(prospectiveIndex);
                            ++numAboveAverage;
                        }
                    }
                    
                    if(numAboveAverage > 0)
                    {
                        averageGoodRadius = averageGoodRadius / RealType(numAboveAverage);
                    }
                    
                    indexIter = indexList.begin();

                    for(; indexIter != indexList.end(); ++indexIter)
                    {
                        IndexType prospectiveIndex = *(indexIter);
                        LabelType prospectiveCenterlineLabel = this->m_CenterlineImage->GetPixel(prospectiveIndex);
                        
                        if(prospectiveCenterlineLabel == 0)
                        {
                            this->m_CenterlineImage->SetPixel(prospectiveIndex, s_MedialCenterlinePredicted);
                            
                            RealType prospectiveMedialness = this->m_MedialnessImage->GetPixel(prospectiveIndex);
                            
                            if(prospectiveMedialness < this->m_MedialnessThreshold)
                            {
                                if(numAboveAverage > 0)
                                {
                                    this->m_RadiusImage->SetPixel(prospectiveIndex, averageGoodRadius);
                                    this->LabelSphereOnCenterlineAtIndex(prospectiveIndex, this->m_CenterlineSphereImage, s_MedialCenterlinePredicted, this->m_RadiusImage->GetPixel(prospectiveIndex), s_MedialSpherePredicted);
                                }
                            }
                            else
                            {
                                this->LabelSphereOnCenterlineAtIndex(prospectiveIndex, this->m_CenterlineSphereImage, s_MedialCenterlinePredicted, this->m_RadiusImage->GetPixel(prospectiveIndex), s_MedialSpherePredicted);
                            }
                        }
                        

                    }
                }
                
            }//end forward and backward loop
        
        }//end centerline condition
        
    }//end centerline loop
    
}

template <class TInputImage, class TOutputImage>                
typename MedialVesselFilter<TInputImage, TOutputImage>::IndexListType 
MedialVesselFilter<TInputImage, TOutputImage>::GetInDirectionIndices(
                const IndexType& centerIndex,
                const RealVectorType& vector, const RealType& angle) const
{
    IndexListType indexList;
    typename TInputImage::SpacingType imageSpacing = this->m_CenterlineImage->GetSpacing();
    
    typename TInputImage::PointType centerPoint;
    this->m_CenterlineImage->TransformIndexToPhysicalPoint(centerIndex, centerPoint);
    
    RealVectorType unitVector = vector;
    unitVector.Normalize();
 
    for(long int i = -1; i <= 1; i++)
    {
        for(long int j = -1; j <= 1; j++)
        {
            for(long int k = -1; k <= 1; k++)
            {               
                typename TInputImage::PointType movingPoint;
                movingPoint[0] = centerPoint[0] + i*imageSpacing[0];
                movingPoint[1] = centerPoint[1] + j*imageSpacing[1]; 
                movingPoint[2] = centerPoint[2] + k*imageSpacing[2];   
                
                RealVectorType movingVector((movingPoint-centerPoint).GetDataPointer());
                if(movingVector.GetNorm() != 0.0)
                {
                    movingVector.Normalize();
                    RealType movingDot = movingVector*unitVector;
                     
                    RealType movingAngle = acos(movingDot);    
                    
                    typename TInputImage::IndexType movingIndex;
                    bool insideImage = this->m_CenterlineImage->TransformPhysicalPointToIndex(movingPoint, movingIndex);
                    
                    typename IndexListType::const_iterator foundIter = std::find(indexList.begin(), indexList.end(), movingIndex);
                    
                    if(insideImage == true && movingIndex != centerIndex && movingAngle >= 0 && movingAngle <= angle && foundIter == indexList.end())
                    {
                        indexList.push_back(movingIndex);
                    }
                }
            }
        }
    }
    
    return(indexList);
}
                
                
template <class TInputImage, class TOutputImage>                
typename MedialVesselFilter<TInputImage, TOutputImage>::IndexListType 
MedialVesselFilter<TInputImage, TOutputImage>::GetInPlaneIndices(
                const IndexType& centerIndex,
                const RealVectorType& vectorPlane1, 
                const RealVectorType& vectorPlane2) const
{
    IndexListType indexList;
    
    typename TInputImage::SpacingType imageSpacing = this->m_CenterlineImage->GetSpacing();
    
    typename TInputImage::PointType centerPoint;
    this->m_CenterlineImage->TransformIndexToPhysicalPoint(centerIndex, centerPoint);

    int numAngles = 8;
    for(unsigned int i = 0; i < numAngles; i++)
    {
        RealType currentAngle = 2.0*M_PI*RealType(i)/RealType(numAngles);
        RealVectorType currentVector = vectorPlane1*cos(currentAngle) + vectorPlane2*sin(currentAngle);
        typename TInputImage::PointType offsetPoint;
        offsetPoint[0] = centerPoint[0] + imageSpacing[0]*currentVector[0];
        offsetPoint[1] = centerPoint[1] + imageSpacing[1]*currentVector[1];
        offsetPoint[2] = centerPoint[2] + imageSpacing[2]*currentVector[2];
        
        typename TInputImage::IndexType currentIndex;
        
        bool insideImage = this->m_CenterlineImage->TransformPhysicalPointToIndex(offsetPoint, currentIndex);
        typename IndexListType::const_iterator foundIter = std::find(indexList.begin(), indexList.end(), currentIndex);

        if(insideImage == true && currentIndex != centerIndex && foundIter == indexList.end())
        {
            indexList.push_back(currentIndex);
        }
    }
    
    return(indexList);
}

template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage, TOutputImage>::
ComputeCenterlineFromMedialness()
{
    typedef std::vector<IndexType>  IndexVectorType; 

    RealVectorType zeroVector;
    zeroVector.Fill(0.0);  
 
    typename RealScalarImageType::Pointer centerInputImage = this->m_LineCenterImage;
    typename RealVectorImageType::Pointer normalVectorImage = this->hessianVesselFilter->GetVesselVectorImage();

    //Get medialness neighboor 
    typename itk::ConstNeighborhoodIterator< RealScalarImageType >::RadiusType 
        neighborhoodRadius;
    unsigned int numNeighborsEachDirection = 1;
    unsigned int numNeighbors = (2*numNeighborsEachDirection+1)*
        (2*numNeighborsEachDirection+1)*(2*numNeighborsEachDirection+1);
    neighborhoodRadius.Fill(numNeighborsEachDirection);
    typename itk::ConstNeighborhoodIterator< RealScalarImageType > 
        neighIter(neighborhoodRadius, centerInputImage, 
        centerInputImage->GetRequestedRegion() );
    
    
    typename LinearInterpolateImageFunction<TInputImage>::Pointer interpolator =  LinearInterpolateImageFunction<TInputImage>::New();
    interpolator->SetInputImage(centerInputImage);
    
    typename TInputImage::SpacingType imageSpacing = this->GetInput()->GetSpacing();
    
    //Determine candidate voxels to compute medialness based on vesselness    
    for (neighIter.GoToBegin(); !neighIter.IsAtEnd(); ++neighIter)
    {   
        RealType currentValue = neighIter.GetCenterPixel();
        typename RealScalarImageType::IndexType currentIndex = 
            neighIter.GetIndex();  
            
        if(currentValue > this->m_MedialnessThreshold && neighIter.InBounds() == true )
        {
            //Get neighborhood around index        
            typename itk::Neighborhood<RealType, 
                itkGetStaticConstMacro(ImageDimension)> hood = 
                neighIter.GetNeighborhood();
            
            //store neighboring indices above medial threshold  
            IndexVectorType neighIndicesAbove;
            
            for (unsigned int i = 0; i < numNeighbors; i++)
            {
                RealType neighValue = hood.GetElement(i);
                IndexType neighIndex = neighIter.GetIndex(i);
                
                /*
                if(neighMedialness == 0.0)
                {
                    RealVectorType neighNormal = normalVectorImage->GetPixel(neighIndex);
                    this->medialFunction->SetNormal(neighNormal);
                    typename MedialFunctionType::TOutput neighMedial = this->medialFunction->EvaluateAtIndex(neighIndex);
                    this->m_MedialnessImage->SetPixel(neighIndex, neighMedial.first);
                    this->m_RadiusImage->SetPixel(neighIndex, neighMedial.second);
                    neighMedialness = this->m_MedialnessImage->GetPixel(neighIndex);               
                }
                */
                
                if(neighValue > this->m_MedialnessThreshold)
                {     
                    neighIndicesAbove.push_back(neighIndex);       
                }
            }
                       
            bool isInPlaneMax = true;
            RealVectorType planeNormal = normalVectorImage->GetPixel(currentIndex);
            
            IndexVectorType inPlaneIndices;
            
            if(planeNormal != zeroVector && neighIndicesAbove.size() >= 3 )
            {
                this->medialFunction->SetNormal(planeNormal);  
                RealVectorType vectorPlane1 = zeroVector; 
                RealVectorType vectorPlane2 = zeroVector;   
                this->medialFunction->ComputeInPlaneVector(vectorPlane1, vectorPlane2);
                
                if(vectorPlane1 == zeroVector || vectorPlane2 == zeroVector )
                {
                    isInPlaneMax = false;
                }
                else
                {
                    typename TInputImage::PointType currentPoint;
                    centerInputImage->TransformIndexToPhysicalPoint(currentIndex, currentPoint );

                    int numAngles = 8;
                    for(unsigned int i = 0; i < numAngles; i++)
                    {
                        RealType currentAngle = 2.0*M_PI*RealType(i)/RealType(numAngles);
                        RealVectorType currentVector = vectorPlane1*cos(currentAngle) + vectorPlane2*sin(currentAngle);
                        typename TInputImage::PointType inPlanePoint;
                        inPlanePoint[0] = currentPoint[0] + imageSpacing[0]*currentVector[0];
                        inPlanePoint[1] = currentPoint[1] + imageSpacing[1]*currentVector[1];
                        inPlanePoint[2] = currentPoint[2] + imageSpacing[2]*currentVector[2];
                        
                        RealType inPlaneValue= interpolator->Evaluate(inPlanePoint);
                        
                        IndexType inPlaneIndex;
                        centerInputImage->TransformPhysicalPointToIndex(inPlanePoint, inPlaneIndex );
                        inPlaneIndices.push_back(inPlaneIndex);
                        
                        if(currentValue < inPlaneValue)
                        {
                            isInPlaneMax = false;
                            break;
                        }
                    }
                }
            }
            
            unsigned int numAboveNotInPlane = 0;
            typename IndexVectorType::const_iterator neighIter = neighIndicesAbove.begin();
            for(; neighIter != neighIndicesAbove.end(); ++neighIter)
            {
                typename IndexVectorType::const_iterator findIter = std::find(inPlaneIndices.begin(), inPlaneIndices.end(), *neighIter);
                if(findIter == inPlaneIndices.end() )
                {
                    ++numAboveNotInPlane;
                }
            }
            
            if(isInPlaneMax == true && numAboveNotInPlane >= 1)
            {
                this->m_CenterlineImage->SetPixel(currentIndex,  s_MedialCenterlineUnlinked);
                this->LabelSphereOnCenterlineAtIndex(currentIndex, this->m_CenterlineSphereImage, s_MedialCenterlineUnlinked, this->m_RadiusImage->GetPixel(currentIndex), s_MedialSphereUnlinked);
            }

        }
    }//end medialness loop
}

template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage, TOutputImage>::ComputeMedialness() 
{
    RealVectorType zeroVector;
    zeroVector.Fill(0.0);  
    
    typename RealScalarImageType::Pointer medialnessInput = this->hessianVesselFilter->GetVesselnessImage();
    
    typename TInputImage::SpacingType pixelSpacing = this->GetInput()->GetSpacing();
    typedef std::vector<RealType> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(pixelSpacing[0]);
    spacingList.push_back(pixelSpacing[1]);
    spacingList.push_back(pixelSpacing[2]);
    
    RealType maxSpacing = *std::max_element(spacingList.begin(), spacingList.end(), std::less_equal<RealType>() );
    RealType minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<RealType>() );
    
    typename LinearInterpolateImageFunction<TInputImage>::Pointer imageInterpolator =  LinearInterpolateImageFunction<TInputImage>::New();
    imageInterpolator->SetInputImage(this->GetInput());

    //Get vesselness neighboor 
    typename itk::ConstNeighborhoodIterator< RealScalarImageType >::RadiusType 
        neighborhoodRadius;
    unsigned int numNeighborsEachDirection = 1;
    unsigned int numNeighbors = (2*numNeighborsEachDirection+1)*
        (2*numNeighborsEachDirection+1)*(2*numNeighborsEachDirection+1);
    neighborhoodRadius.Fill(numNeighborsEachDirection);
    typename itk::ConstNeighborhoodIterator< RealScalarImageType > 
        vesselnessNeighIter(neighborhoodRadius, medialnessInput, 
        medialnessInput->GetRequestedRegion() );

    const TInputImage* inputImage = this->GetInput();
    
    //Determine candidate voxels to compute medialness based on vesselness    
    for (vesselnessNeighIter.GoToBegin(); !vesselnessNeighIter.IsAtEnd(); 
        ++vesselnessNeighIter)
    {   
        RealType currentVesselness = vesselnessNeighIter.GetCenterPixel();
        typename RealScalarImageType::IndexType currentIndex = 
            vesselnessNeighIter.GetIndex();
            
        typename itk::Neighborhood<RealType, 
            itkGetStaticConstMacro(ImageDimension)> vesselHood = 
            vesselnessNeighIter.GetNeighborhood();
        
        //count number of neighbors with vesselness above threshold    
        unsigned int vesselnessCounter = 0;
        for (unsigned int i = 0; i < numNeighbors; i++)
        {
            RealType neighVesselness = vesselHood.GetElement(i);
            if(neighVesselness > this->m_VesselnessThreshold)
            {
                ++vesselnessCounter;              
            }
        }
        
        //vectors defining vessel plane at this voxel
        RealVectorType planeNormal = this->hessianVesselFilter->GetVesselVectorImage()->GetPixel(currentIndex);
             
        //Perform medialness calculation if current voxel has at least 2 
        //neighbors with vesselness above threshold
        if(vesselnessCounter >=3 
            && currentVesselness > this->m_VesselnessThreshold  
            && planeNormal != zeroVector 
            && inputImage->GetPixel(currentIndex) >= this->m_IntensityMinimum
            && inputImage->GetPixel(currentIndex) <= this->m_IntensityMaximum )
        {
                      
            this->medialFunction->SetNormal(planeNormal);
            typename MedialFunctionType::TOutput medialOutput = this->medialFunction->EvaluateAtIndex(currentIndex);
            RealType currentMedialness = medialOutput.first;
            RealType currentRadius = medialOutput.second;
           
            RealType currentIntensity = this->GetInput()->GetPixel(currentIndex);
            typename TInputImage::PointType currentPoint;
            this->GetInput()->TransformIndexToPhysicalPoint( currentIndex, currentPoint);
            typename MedialFunctionType::VectorType planeVector1;
            typename MedialFunctionType::VectorType planeVector2;
            this->medialFunction->ComputeInPlaneVector(planeVector1, planeVector2);
            planeVector1 = (currentRadius + 1.5*maxSpacing)*planeVector1;
            planeVector2 = (currentRadius + 1.5*maxSpacing)*planeVector2;
            
            int numAngles = int(ceil( 2.0*M_PI * (currentRadius/minSpacing) + 4));
            unsigned int numAboveDiff = 0;
            for(unsigned int i = 0; i < numAngles; i++)
            {
                RealType currentAngle = 2.0*M_PI*RealType(i)/RealType(numAngles);
                RealVectorType currentVector = planeVector1*cos(currentAngle) + planeVector2*sin(currentAngle);
                typename TInputImage::PointType inPlanePoint = currentPoint + currentVector;
               
                if(imageInterpolator->IsInsideBuffer(inPlanePoint) == true)
                {
                    RealType inPlaneIntensity = imageInterpolator->Evaluate(inPlanePoint);
                    
                    if(currentIntensity - inPlaneIntensity > this->m_IntensityDifference)
                    {
                        ++numAboveDiff;
                    }
                }
            }
            
            if(RealType(numAboveDiff)/RealType(numAngles) >= this->m_EdgeRatio)
            {
                this->m_MedialnessImage->SetPixel( currentIndex, currentMedialness);
                this->m_RadiusImage->SetPixel( currentIndex, currentRadius);    
                
//                std::cout << numAboveDiff << " " << numAngles << " " << currentIntensity << " " << currentMedialness << std::endl;
            } 
   
        }//end vesselness condition
        
    }//end vesselness loop   
    
    //foo
    std::cout << "Completed Medialness " << std::endl;

    typedef LineTransform3DFunction<TInputImage, RealType> LineFunctionType;
    typedef LineTransformMinVariance3DFunction<TInputImage, RealType> LineMinVarFunctionType;
    typedef VectorLineTransform3DFunction<RealVectorImageType, RealType> VectorLineFunctionType;
    typedef CorrelateVectorsOnLine3DFunction<RealVectorImageType, RealType> VectorLineCorrelateFunctionType;

    
    typename RealScalarImageType::Pointer startImage = this->m_MedialnessImage;
    typename RealScalarImageType::Pointer startRadiusImage = this->m_RadiusImage;

    typename LineFunctionType::Pointer lineFunction = LineFunctionType::New();
    lineFunction->SetInputImage(startImage);
    lineFunction->SetParameters(1*maxSpacing, 10);

    
    typename RealVectorImageType::Pointer lineVectorImage = RealVectorImageType::New();
    lineVectorImage->CopyInformation(this->GetInput());
    lineVectorImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
    lineVectorImage->SetRequestedRegionToLargestPossibleRegion();
    lineVectorImage->Allocate();
    lineVectorImage->FillBuffer(zeroVector); 
    
    
    typename RealScalarImageType::Pointer lineImage = RealScalarImageType::New();
    lineImage->CopyInformation(this->GetInput());
    lineImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
    lineImage->SetRequestedRegionToLargestPossibleRegion();
    lineImage->Allocate();
    lineImage->FillBuffer(0.0); 
    
    this->m_LineCenterImage = RealScalarImageType::New();
    this->m_LineCenterImage->CopyInformation(this->GetInput());
    this->m_LineCenterImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
    this->m_LineCenterImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_LineCenterImage->Allocate();
    this->m_LineCenterImage->FillBuffer(0.0); 
  
    typename itk::ConstNeighborhoodIterator< RealScalarImageType > 
        startIter(neighborhoodRadius, startImage, startImage->GetLargestPossibleRegion() );
        
     for (startIter.GoToBegin(); !startIter.IsAtEnd();++startIter)
    {   
        RealType currentValue = startIter.GetCenterPixel();
        typename RealScalarImageType::IndexType currentIndex = 
            startIter.GetIndex();
            
        typename itk::Neighborhood<RealType, 
            itkGetStaticConstMacro(ImageDimension)> hood = 
            startIter.GetNeighborhood();
        
        //count number of neighbors with vesselness above threshold    
        unsigned int counter = 0;
        for (unsigned int i = 0; i < numNeighbors; i++)
        {
            RealType neighValue = hood.GetElement(i);
            if(neighValue > this->m_MedialnessThreshold*this->m_MedialnessPercent)
            {
                ++counter;              
            }
        }
        
        if(counter >= 3)
        {
            RealType currentRadius = startRadiusImage->GetPixel(currentIndex);
            
            lineFunction->SetParameters(currentRadius/2.0, 10);
            typename LineFunctionType::TOutput lineOutput = lineFunction->EvaluateAtIndex( currentIndex );
            
            lineImage->SetPixel( currentIndex, lineOutput.GetNorm());
            //this->m_TestImage->SetPixel( currentIndex, lineOutput.GetNorm() );
            
            if(lineOutput.GetNorm() > 0.0)
            {
                lineOutput.Normalize();
            }
            
            lineVectorImage->SetPixel( currentIndex, lineOutput);
        }
    }

    typename VectorLineFunctionType::Pointer vectorLineFunction = VectorLineFunctionType::New();
    vectorLineFunction->SetInputImage(lineVectorImage);
    vectorLineFunction->SetParameters(1*maxSpacing, 10);
 
    itk::ImageRegionConstIterator<RealScalarImageType>  
        lineIter(lineImage, lineImage->GetRequestedRegion() );

    for (lineIter.GoToBegin(); !lineIter.IsAtEnd(); ++lineIter)
    { 
        if(lineIter.Get() >  this->m_MedialnessThreshold * this->m_MedialnessPercent)
        {
            typename RealScalarImageType::IndexType currentIndex = lineIter.GetIndex();
            RealType currentRadius = startRadiusImage->GetPixel(currentIndex);
            
            vectorLineFunction->SetParameters(currentRadius/2.0, 10);
            typename VectorLineFunctionType::TOutput vectorLineOutput = vectorLineFunction->EvaluateAtIndex( currentIndex );
            RealVectorType lineVector = vectorLineOutput.second;
            
            this->m_LineCenterImage->SetPixel( currentIndex, vectorLineOutput.first * lineImage->GetPixel(currentIndex));
            
        }
    }
 
}    


template <class TInputImage, class TOutputImage>
void
MedialVesselFilter<TInputImage, TOutputImage>::SetVesselFilter(HessianVesselFilterType* h )
{
    this->hessianVesselFilter = h;
}

template <class TInputImage, class TOutputImage>
void
MedialVesselFilter<TInputImage, TOutputImage>::SetMedialFunction(MedialFunctionType* m)
{
    this->medialFunction = m;
}


//Print
template <class TInputImage, class TOutputImage>
void MedialVesselFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

//  os << indent << "foo:  " << m_foo << std::endl;

}


 
}//end namespace                
