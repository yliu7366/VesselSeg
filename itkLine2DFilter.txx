#include "itkLine2DFilter.h"


#include "itkConstNeighborhoodIterator.h"

namespace itk
{
template <class TInputImage, class TOutputImage>
Line2DFilter<TInputImage,TOutputImage>::Line2DFilter()
{
    
}

template <class TInputImage, class TOutputImage>
Line2DFilter<TInputImage,TOutputImage>::~Line2DFilter()
{
    
}

template <class TInputImage, class TOutputImage>
void Line2DFilter<TInputImage, TOutputImage>::GenerateData()
{
}

template <class TInputImage, class TOutputImage>
void Line2DFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream& os, Indent indent) const
{
}

template <class TInputImage, class TOutputImage>
void Line2DFilter<TInputImage, TOutputImage>::SetLineFunction(LineFunctionType* l)
{
    this->lineFunction = l;
}
template <class TInputImage, class TOutputImage>
void Line2DFilter<TInputImage, TOutputImage>::SetVectorLineFunction(VectorLineFunctionType* l)
{
    this->vectorLineFunction = l;
}

template <class TInputImage, class TOutputImage>
void Line2DFilter<TInputImage, TOutputImage>::Allocate()
{
    typename TInputImage::ConstPointer inputImage = this->GetInput();
    
    this->m_LineImage = RealScalarImageType::New();
    this->m_LineImage->CopyInformation(inputImage);
    this->m_LineImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_LineImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_LineImage->Allocate(); 
    this->m_LineImage->FillBuffer(0.0);
    
    RealVectorType zeroVector;
    zeroVector.Fill(0.0);  
    this->m_LineVectorImage = RealVectorImageType::New();
    this->m_LineVectorImage->CopyInformation(inputImage);
    this->m_LineVectorImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_LineVectorImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_LineVectorImage->Allocate();
    this->m_LineVectorImage->FillBuffer(zeroVector); 
}

template <class TInputImage, class TOutputImage>
void Line2DFilter<TInputImage, TOutputImage>::ComputeLineFilter()
{
    itk::ImageRegionConstIterator<TInputImage>  
        imageIter(this->GetInput(), 
        this->GetInput()->GetRequestedRegion() );
    
    typename RealScalarImageType::Pointer intermediateLineImage = RealScalarImageType::New();
    intermediateLineImage->CopyInformation(this->GetInput());
    intermediateLineImage->SetRegions(this->GetInput()->GetLargestPossibleRegion() );
    intermediateLineImage->SetRequestedRegionToLargestPossibleRegion();
    intermediateLineImage->Allocate(); 
    intermediateLineImage->FillBuffer(0.0);
    
    for(imageIter.GoToBegin(); !imageIter.IsAtEnd(); ++imageIter)
    {
        typename TInputImage::ValueType currentIntensity = imageIter.Get();

        if(currentIntensity >= this->m_IntensityMinimum && currentIntensity <= this->m_IntensityMaximum)
        {
            typename TInputImage::IndexType currentIndex = imageIter.GetIndex();
            typename LineFunctionType::TOutput lineOutput = lineFunction->EvaluateAtIndex( currentIndex );
            
            //this->m_LineImage->SetPixel( currentIndex, lineOutput.GetNorm());
            intermediateLineImage->SetPixel(currentIndex, lineOutput.GetNorm());
            if(lineOutput.GetNorm() > 0)
            {
                lineOutput.Normalize();
            }
            this->m_LineVectorImage->SetPixel(currentIndex, lineOutput);
        }
    }
    
    vectorLineFunction->SetInputImage(this->m_LineVectorImage);
    
    itk::ImageRegionConstIterator<RealScalarImageType>  
        intermediateIter(intermediateLineImage, 
        intermediateLineImage->GetRequestedRegion() );
    
    for(intermediateIter.GoToBegin(); !intermediateIter.IsAtEnd(); ++intermediateIter)
    {
        typename TInputImage::ValueType currentValue = intermediateIter.Get();

        if(currentValue >= this->m_IntensityMinimum && currentValue <= this->m_IntensityMaximum)
        {
            typename TInputImage::IndexType currentIndex = intermediateIter.GetIndex();
            typename VectorLineFunctionType::TOutput lineVectorOutput = vectorLineFunction->EvaluateAtIndex( currentIndex );
            
            this->m_LineImage->SetPixel( currentIndex, lineVectorOutput.first);
        }
    }
}


}
