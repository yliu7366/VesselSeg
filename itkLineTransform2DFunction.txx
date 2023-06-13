#include <cmath>
#include <vector>
#include <functional>

#include "itkLineTransform2DFunction.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"


namespace itk
{


template < typename TInputImage, typename TCoordRep >
LineTransform2DFunction<TInputImage, TCoordRep>::LineTransform2DFunction()
{
    //typename LinearInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater = LinearInterpolateImageFunction<TInputImage>::New();
    typename NearestNeighborInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater = NearestNeighborInterpolateImageFunction<TInputImage>::New();
    defaultIterpolater->Register();
    this->interpolator = defaultIterpolater;
    
    this->radius = 0.0;
    this->maxNumAngles = 0;
    this->minSpacing = 0.0;
    
}

template < typename TInputImage, typename TCoordRep >
LineTransform2DFunction<TInputImage, TCoordRep>::~LineTransform2DFunction()
{
}

template < typename TInputImage, typename TCoordRep >
void LineTransform2DFunction<TInputImage, TCoordRep>::
SetInputImage(const TInputImage * ptr)
{
    Superclass::SetInputImage(ptr);
    this->interpolator->SetInputImage(this->m_Image);   
    
    typename TInputImage::SpacingType pixelSpacing = this->m_Image->GetSpacing();
    typedef std::vector<TCoordRep> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(pixelSpacing[0]);
    spacingList.push_back(pixelSpacing[1]);
    
    this->minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<TCoordRep>() );
}

template < typename TInputImage, typename TCoordRep >
void LineTransform2DFunction<TInputImage, TCoordRep>::
SetInterpolater(InterpolaterType* i)
{
    this->interpolator = i;
}

template < typename TInputImage, typename TCoordRep >
void LineTransform2DFunction<TInputImage, TCoordRep>::
SetParameters(TCoordRep radius, unsigned int maxNumAngles)
{
    this->radius = radius;
    this->maxNumAngles = maxNumAngles;
}


template < typename TInputImage, typename TCoordRep >
typename LineTransform2DFunction<TInputImage, TCoordRep>::TOutput LineTransform2DFunction<TInputImage, TCoordRep>
::Evaluate(const PointType &point) const
{
    VectorType lineResponse;
    lineResponse.Fill(0.0);

    int numRadii = int(ceil(this->radius/ TCoordRep(this->minSpacing)));
    if(numRadii < 1)
    {
        numRadii = 1;
    }
    int numAnglesA = int(ceil( 2.0*M_PI * (this->radius/this->minSpacing) + 2));
        
    if(numAnglesA > this->maxNumAngles)
    {
        numAnglesA = this->maxNumAngles;
    }
        
    for(int i = 0; i < numAnglesA; i++)
    {
        TCoordRep azimuth = 2.0*M_PI * TCoordRep(i) / TCoordRep(numAnglesA);
        
        TCoordRep currentLine = 0.0;
        
        for(int k = -numRadii; k <= numRadii; k++)
        {
            TCoordRep currentRadius = this->radius * TCoordRep(k) / TCoordRep(numRadii);
            
            VectorType offset;
            offset[0] =  currentRadius * cos(azimuth);
            offset[1] =  currentRadius * sin(azimuth);
            
            PointType currentPoint = point + offset;
        
            if(this->interpolator->IsInsideBuffer(currentPoint) == true )
            {
                TCoordRep currentValue = this->interpolator->Evaluate(currentPoint);
                currentLine += currentValue; 
            }
        }
        currentLine = currentLine / TCoordRep((2*numRadii+1));
        if(currentLine > lineResponse.GetNorm())
        {
            lineResponse[0] = currentLine * cos(azimuth);
            lineResponse[1] = currentLine * sin(azimuth);           
        }

    }
        
    return(lineResponse);
}

template < typename TInputImage, typename TCoordRep >
typename LineTransform2DFunction<TInputImage, TCoordRep>::TOutput LineTransform2DFunction<TInputImage, TCoordRep>::
EvaluateAtContinuousIndex (const ContinuousIndexType &index) const
{
    PointType point;
    this->m_Image->TransformContinuousIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}

template < typename TInputImage, typename TCoordRep >
typename LineTransform2DFunction<TInputImage, TCoordRep>::TOutput LineTransform2DFunction<TInputImage, TCoordRep>::
EvaluateAtIndex (const IndexType &index) const
{
    PointType point;
    this->m_Image->TransformIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}

//
template < typename TInputImage, typename TCoordRep >
typename LineTransformMinVariance2DFunction<TInputImage, TCoordRep>::TOutput LineTransformMinVariance2DFunction<TInputImage, TCoordRep>
::Evaluate(const PointType &point) const
{
    VectorType lineResponse;
    lineResponse.Fill(0.0);

    int numRadii = int(ceil(this->radius/ TCoordRep(this->minSpacing)));
    if(numRadii < 1)
    {
        numRadii = 1;
    }
    
    int numAnglesA = int(ceil( 2.0*M_PI * (this->radius/this->minSpacing) + 2));
        
    if(numAnglesA > this->maxNumAngles)
    {
        numAnglesA = this->maxNumAngles;
    }
    
    
    for(int i = 0; i < numAnglesA; i++)
    {
        TCoordRep azimuth = 2.0*M_PI * TCoordRep(i) / TCoordRep(numAnglesA);
      
        TCoordRep currentSum = 0.0;
        TCoordRep currentSumSq = 0.0;
        
        for(int k = -numRadii; k <= numRadii; k++)
        {
            TCoordRep currentRadius = this->radius * TCoordRep(k) / TCoordRep(numRadii);
            
            VectorType offset;
            offset[0] =  currentRadius * cos(azimuth);
            offset[1] =  currentRadius * sin(azimuth);
            
            PointType currentPoint = point + offset;
        
            if(this->interpolator->IsInsideBuffer(currentPoint) == true )
            {
                TCoordRep currentValue = this->interpolator->Evaluate(currentPoint);
                currentSum += currentValue; 
                currentSumSq += currentValue*currentValue;
            }
        }
        TCoordRep mean = currentSum / TCoordRep((2*numRadii+1));
        TCoordRep currentLine = (currentSumSq -currentSum*mean) / TCoordRep((2*numRadii+1));
        
        
        if(i==0)
        {
            lineResponse[0] = cos(azimuth);
            lineResponse[1] = sin(azimuth);           
        }
        else if(currentLine < lineResponse.GetNorm())
        {
            lineResponse[0] = cos(azimuth);
            lineResponse[1] = sin(azimuth);
        }
    
    }
        
    return(lineResponse);
}

//
template < typename TInputImage, typename TCoordRep >
VectorLineTransform2DFunction<TInputImage, TCoordRep>::VectorLineTransform2DFunction()
{
    //typename VectorLinearInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater = VectorLinearInterpolateImageFunction<TInputImage>::New();
    typename VectorNearestNeighborInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater = VectorNearestNeighborInterpolateImageFunction<TInputImage>::New();
    defaultIterpolater->Register();
    this->interpolator = defaultIterpolater;
    
    this->radius = 0.0;
    this->maxNumAngles = 0;
    this->minSpacing = 0.0;
    
}

template < typename TInputImage, typename TCoordRep >
VectorLineTransform2DFunction<TInputImage, TCoordRep>::~VectorLineTransform2DFunction()
{
}

template < typename TInputImage, typename TCoordRep >
void VectorLineTransform2DFunction<TInputImage, TCoordRep>::
SetInputImage(const TInputImage * ptr)
{
    Superclass::SetInputImage(ptr);
    this->interpolator->SetInputImage(this->m_Image);   
    
    typename TInputImage::SpacingType pixelSpacing = this->m_Image->GetSpacing();
    typedef std::vector<TCoordRep> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(pixelSpacing[0]);
    spacingList.push_back(pixelSpacing[1]);
    
    this->minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<TCoordRep>() );
}


template < typename TInputImage, typename TCoordRep >
void VectorLineTransform2DFunction<TInputImage, TCoordRep>::
SetInterpolater(InterpolaterType* i)
{
    this->interpolator = i;
}



template < typename TInputImage, typename TCoordRep >
void VectorLineTransform2DFunction<TInputImage, TCoordRep>::
SetParameters(TCoordRep radius, unsigned int maxNumAngles)
{
    this->radius = radius;
    this->maxNumAngles = maxNumAngles;
}


template < typename TInputImage, typename TCoordRep >
typename VectorLineTransform2DFunction<TInputImage, TCoordRep>::TOutput 
VectorLineTransform2DFunction<TInputImage, TCoordRep>
::Evaluate(const PointType &point) const
{
    TCoordRep maxCorrelation = 0.0;
    TCoordRep maxAzimuth = 0.0;
    VectorType maxOrientation;

    int numRadii = int(ceil(this->radius/ TCoordRep(this->minSpacing)));
    if(numRadii < 1)
    {
        numRadii = 1;
    }
    int numAnglesA = int(ceil( 2.0*M_PI * (this->radius/this->minSpacing) + 2));
        
    if(numAnglesA > this->maxNumAngles)
    {
        numAnglesA = this->maxNumAngles;
    }

    for(int i = 0; i < numAnglesA; i++)
{
        TCoordRep azimuth = 2.0*M_PI * TCoordRep(i) / TCoordRep(numAnglesA);
        
        VectorType orientation;
        orientation[0] =  cos(azimuth);
        orientation[1] =  sin(azimuth);
        
        TCoordRep currentLine = 0.0;
        
        for(int k = -numRadii; k <= numRadii; k++)
        {
            TCoordRep currentRadius = this->radius * TCoordRep(k) / TCoordRep(numRadii);
            
            VectorType offset;
            offset[0] =  currentRadius * cos(azimuth);
            offset[1] =  currentRadius * sin(azimuth);
            
            PointType currentPoint = point + offset;

            if(this->interpolator->IsInsideBuffer(currentPoint) == true )
            {
                VectorType currentValue = this->interpolator->Evaluate(currentPoint);
                currentLine += fabs(currentValue*orientation);
            }
        }
        currentLine = currentLine / TCoordRep((2*numRadii+1));
        if(currentLine > maxCorrelation)
        {
            maxCorrelation = currentLine;
            maxAzimuth = azimuth;
            maxOrientation[0] = this->radius;
            maxOrientation[1] = maxAzimuth;              
        }
    
    }
    
    TOutput output;
    output.first = maxCorrelation;
    output.second = maxOrientation;
    
    return(output);
}

template < typename TInputImage, typename TCoordRep >
typename VectorLineTransform2DFunction<TInputImage, TCoordRep>::TOutput 
VectorLineTransform2DFunction<TInputImage, TCoordRep>::
EvaluateAtContinuousIndex (const ContinuousIndexType &index) const
{
    PointType point;
    this->m_Image->TransformContinuousIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}

template < typename TInputImage, typename TCoordRep >
typename VectorLineTransform2DFunction<TInputImage, TCoordRep>::TOutput 
VectorLineTransform2DFunction<TInputImage, TCoordRep>::
EvaluateAtIndex (const IndexType &index) const
{
    PointType point;
    this->m_Image->TransformIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}


//foo
template < typename TInputImage, typename TCoordRep >
CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::CorrelateVectorsOnLine2DFunction()
{
    //typename VectorLinearInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater1 = VectorLinearInterpolateImageFunction<TInputImage>::New();
    typename VectorNearestNeighborInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater1 = VectorNearestNeighborInterpolateImageFunction<TInputImage>::New();
    defaultIterpolater1->Register();
    this->interpolator1 = defaultIterpolater1;
    
    //typename VectorLinearInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater2 = VectorLinearInterpolateImageFunction<TInputImage>::New();
    typename VectorNearestNeighborInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater2 = VectorNearestNeighborInterpolateImageFunction<TInputImage>::New();
    defaultIterpolater2->Register();
    this->interpolator2 = defaultIterpolater2;
    
    this->radius = 0.0;
    this->elevation = 0.0;
    this->minSpacing = 0.0;
    
}

template < typename TInputImage, typename TCoordRep >
CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::~CorrelateVectorsOnLine2DFunction()
{

}

template < typename TInputImage, typename TCoordRep >
void CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::
SetInputImages(const TInputImage * ptr1, const TInputImage * ptr2)
{
    Superclass::SetInputImage(ptr1);
    this->interpolator1->SetInputImage(ptr1); 
    
    this->inputImage2 = ptr2;  
    this->interpolator2->SetInputImage(ptr2);  
    
    typename TInputImage::SpacingType pixelSpacing = this->m_Image->GetSpacing();
    typedef std::vector<TCoordRep> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(pixelSpacing[0]);
    spacingList.push_back(pixelSpacing[1]);
    
    this->minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<TCoordRep>() );
}

template < typename TInputImage, typename TCoordRep >
void CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::
SetInterpolaters(InterpolaterType* i1, InterpolaterType* i2)
{
    this->interpolator1 = i1;
    this->interpolator2 = i2;
}

template < typename TInputImage, typename TCoordRep >
void CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::
SetParameters(TCoordRep radius, TCoordRep azimuth, TCoordRep elevation)
{
    this->radius = radius;
    this->azimuth = azimuth;
}


template < typename TInputImage, typename TCoordRep >
typename CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::TOutput CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>
::Evaluate(const PointType &point) const
{
    TOutput lineResponse = 0.0;

    int numRadii = int(ceil(this->radius/ TCoordRep(this->minSpacing)));
    
    if(numRadii < 1)
    {
        numRadii = 1;
    }
            
    for(int k = -numRadii; k <= numRadii; k++)
    {
        TCoordRep currentRadius = this->radius * TCoordRep(k) / TCoordRep(numRadii);
        
        VectorType offset;
        offset[0] =  currentRadius * cos(this->azimuth);
        offset[1] =  currentRadius * sin(this->azimuth);
        
        PointType currentPoint = point + offset;
    
        if(this->interpolator1->IsInsideBuffer(currentPoint) == true &&  this->interpolator2->IsInsideBuffer(currentPoint) == true)
        {
            typename TInputImage::PixelType currentValue1 = this->interpolator1->Evaluate(currentPoint);
            typename TInputImage::PixelType currentValue2 = this->interpolator2->Evaluate(currentPoint);
            if(currentValue1.GetNorm() > 0)
            {
                currentValue1.Normalize();
            }
            if(currentValue2.GetNorm() > 0)
            {
                currentValue2.Normalize();
            }
            lineResponse += fabs(currentValue1*currentValue2); 
        }
    }
    lineResponse = lineResponse / TCoordRep(2*numRadii+1);
    return(lineResponse);
}

template < typename TInputImage, typename TCoordRep >
typename CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::TOutput CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::
EvaluateAtContinuousIndex (const ContinuousIndexType &index) const
{
    PointType point;
    this->m_Image->TransformContinuousIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}

template < typename TInputImage, typename TCoordRep >
typename CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::TOutput CorrelateVectorsOnLine2DFunction<TInputImage, TCoordRep>::
EvaluateAtIndex (const IndexType &index) const
{
    PointType point;
    this->m_Image->TransformIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}


}//end namespace itk


