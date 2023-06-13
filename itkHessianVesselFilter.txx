#include <cmath>
#include <algorithm>

#include "itkHessianVesselFilter.h"

#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkDiscreteHessianGaussianImageFunction.h"
#include "itkConstNeighborhoodIterator.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"

#include "itkLinearInterpolateImageFunction.h"

namespace itk
{

//Constructor
template <class TInputImage, class TOutputImage>
HessianVesselFilter<TInputImage,TOutputImage>::HessianVesselFilter()
{
    
}

//Destructor
template <class TInputImage, class TOutputImage>
HessianVesselFilter<TInputImage,TOutputImage>::~HessianVesselFilter()
{
}

//GenerateData
template <class TInputImage, class TOutputImage>
void HessianVesselFilter<TInputImage, TOutputImage>::GenerateData()
{
    this->Allocate();
    this->ComputeVesselness();
//    this->GraftOutput(?);
//    this->GraftNthOutput(0, ?);
}



template <class TInputImage, class TOutputImage>
void HessianVesselFilter<TInputImage, TOutputImage>::Allocate()
{
    typename TInputImage::ConstPointer inputImage = this->GetInput();
    
    //Scalar images
    this->m_VesselnessImage = RealScalarImageType::New();
    this->m_VesselnessImage->CopyInformation(inputImage);
    this->m_VesselnessImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_VesselnessImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_VesselnessImage->Allocate(); 
    this->m_VesselnessImage->FillBuffer(0.0);
      
    this->m_ScaleImage = RealScalarImageType::New();
    this->m_ScaleImage->CopyInformation(inputImage);
    this->m_ScaleImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_ScaleImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_ScaleImage->Allocate(); 
    this->m_ScaleImage->FillBuffer(0.0);
    
    this->m_VesselMedialImage = RealScalarImageType::New();
    this->m_VesselMedialImage->CopyInformation(inputImage);
    this->m_VesselMedialImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_VesselMedialImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_VesselMedialImage->Allocate(); 
    this->m_VesselMedialImage->FillBuffer(0.0);
    
    //vector images 
    RealVectorType zeroVector;
    zeroVector.Fill(0.0);
    
    this->m_VesselVectorImage = RealVectorImageType::New();
    this->m_VesselVectorImage->CopyInformation(inputImage);
    this->m_VesselVectorImage->SetRegions(inputImage->GetLargestPossibleRegion() );
    this->m_VesselVectorImage->SetRequestedRegionToLargestPossibleRegion();
    this->m_VesselVectorImage->Allocate();
    this->m_VesselVectorImage->FillBuffer(zeroVector);
    
}

template <class TInputImage, class TOutputImage>
void HessianVesselFilter<TInputImage, TOutputImage>::ComputeVesselness() 
{
    //loop through all scales
    unsigned int numScales = this->m_Scales.size();
    for(unsigned int j = 0; j < numScales; j++)
    {
        //Filter to construct Hessian matrix
        typename itk::HessianRecursiveGaussianImageFilter<RealScalarImageType, 
            RealSymmetricMatrixImageType>::Pointer hessianFilter = 
            itk::HessianRecursiveGaussianImageFilter<TInputImage, 
            RealSymmetricMatrixImageType>::New();
    
        // Convolve with "normalized Gaussian" in each direction
        // 1/\sqrt(2s^2 \pi) \exp(-|x|^2 / (2 s^2) ) 
        hessianFilter->SetNormalizeAcrossScale(false);    
        hessianFilter->SetInput( this->GetInput() );
        RealType currentScale = this->m_Scales[j];
        hessianFilter->SetSigma( currentScale );
        hessianFilter->Update();
        typename RealSymmetricMatrixImageType::Pointer currentHessianImage = 
            hessianFilter->GetOutput();
        itk::ImageRegionConstIterator<RealSymmetricMatrixImageType>  
            currentHessianIter(currentHessianImage, 
            currentHessianImage->GetRequestedRegion() );
        
        for(currentHessianIter.GoToBegin(); !currentHessianIter.IsAtEnd();
            ++currentHessianIter)
        {
            //Get Hessian matrix for this voxel
            RealSymmetricMatrixType hessianMatrix = 
                currentHessianIter.Get();
            typename RealSymmetricMatrixImageType::IndexType index = 
                currentHessianIter.GetIndex(); 
            
            if( this->GetInput()->GetPixel(index) >= this->m_IntensityMinimum &&
                this->GetInput()->GetPixel(index) <= this->m_IntensityMaximum )
            {        
                RealVectorType eigenValues;
                RealMatrixType eigenVectors;
                
                //Compute EigenSystem
                this->ComputeEigenSystem( hessianMatrix, eigenValues, 
                    eigenVectors);
                
                //Scale eigenvalues by sigma squared   
                RealType smallestEValue = currentScale*currentScale*eigenValues[0];
                RealType middleEValue = currentScale*currentScale*eigenValues[1];
                RealType largestEValue = currentScale*currentScale*eigenValues[2];
                
                //Compute Vesselness
                RealType vesselness = this->VesselnessMeasure(smallestEValue, 
                    middleEValue, largestEValue, this->m_VesselnessAlpha, 
                    this->m_VesselnessBeta, this->m_VesselnessGamma);
                
                /*
                if(vesselness > 0)
                {
                    RealVectorType normalVector;
                    normalVector.SetVnlVector(eigenVectors.GetVnlMatrix().get_row(0));
                    this->medialFunction->SetRadii(0.5*currentScale, 1.5*currentScale, 5, 10);
                    this->medialFunction->SetNormal(normalVector);
                    typename MedialFunctionType::TOutput medialOutput = this->medialFunction->EvaluateAtIndex(index);
                    RealType medialness = medialOutput.first;
                    RealType radius = medialOutput.second;
          
                    vesselness = (vesselness+medialness)/(2.0);
                }
                */
                //std::cout << "Scale " << currentScale << " Index " << index << " Vector " << vesselVector << " Medialness " << vesselness << " Radius " << radius << std::endl;
                    
                
                //Save largest vesselness so far and associated scale and 
                //coordinate system            
                if(vesselness > this->m_VesselnessImage->GetPixel( index ) )
                {
                    this->m_VesselnessImage->SetPixel( index, vesselness );
                    this->m_ScaleImage->SetPixel( index, currentScale ); 
                    
                    //get vectors
                    RealVectorType vesselVector;
                    vesselVector.SetVnlVector(eigenVectors.GetVnlMatrix().get_row(0));
                
                    
                    //normalize vectors
                    if(vesselVector.GetNorm() > 0.0)
                    {
                        vesselVector.Normalize();
                        this->m_VesselVectorImage->SetPixel( index, vesselVector); 
                    }
                }
                
            }//end largest vesselness condition   
   
        }//end hessian loop
            
    }//end scales loop
}
    
template <class TInputImage, class TOutputImage>
void HessianVesselFilter<TInputImage, TOutputImage>::ComputeMedialVesselness() 
{    
    
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
    itk::ImageRegionConstIterator<RealScalarImageType>  
        vesselnessIter(this->m_VesselnessImage, 
        this->m_VesselnessImage->GetRequestedRegion() );

    for(vesselnessIter.GoToBegin(); !vesselnessIter.IsAtEnd();
        ++vesselnessIter)
    {
        RealType currentVesselness = vesselnessIter.Get();
        if(currentVesselness > 0)
        {
            typename TInputImage::IndexType currentIndex = vesselnessIter.GetIndex();
            RealType currentScale = this->m_ScaleImage->GetPixel(currentIndex);
            this->medialFunction->SetRadii(0.5*currentScale, 2*currentScale, 5, 20);
            
            RealVectorType normalVector = this->m_VesselVectorImage->GetPixel( currentIndex ); 
            this->medialFunction->SetNormal(normalVector);
            typename MedialFunctionType::TOutput medialOutput = this->medialFunction->EvaluateAtIndex(currentIndex);
            RealType currentMedialness = medialOutput.first;
            RealType currentRadius = medialOutput.second;
            
            //RealType currentRadius = currentScale;
            
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
            
            if(RealType(numAboveDiff)/RealType(numAngles) >= this->m_EdgeRatio )
            {
                this->m_VesselMedialImage->SetPixel( currentIndex, this->m_VesselnessImage->GetPixel( currentIndex )); 
            } 
            
        }
    }
    
}

//Compute vesselness measure
template <class TInputImage, class TOutputImage>
typename HessianVesselFilter<TInputImage, TOutputImage>::RealType 
    HessianVesselFilter<TInputImage, TOutputImage>::VesselnessMeasure(
    RealType lambda1, RealType lambda2, RealType lambda3, RealType alpha, 
    RealType beta, RealType gamma) const
{

    RealType vesselMeasure = 0.0;
    if (   lambda2 >= 0.0 ||  lambda3 >= 0.0 )
    {
        vesselMeasure = 0.0;
    }
    else
    {
        //Sato    
        //RealType arg1 = fabs(lambda2/lambda3);
        //RealType arg2 = fabs(lambda1/lambda2);
        //vesselMeasure = (1.0-arg2)*arg1;
        
        //Frangi like
        RealType arg1 = -(lambda2*lambda2) / (2.0*alpha*alpha*lambda3*lambda3);
        RealType arg2 = -(lambda1*lambda1) / (2.0*beta*beta*fabs(lambda2*lambda3));
        RealType arg3 = -( (lambda1*lambda1) + (lambda2*lambda2) + (lambda3*lambda3) ) / (2.0*gamma*gamma);
        vesselMeasure = (1.0 - exp(arg1))*exp(arg2)*(1.0 - exp(arg3));        
    }
    
    return(vesselMeasure);                  
}

//Compute Eigen system
template <class TInputImage, class TOutputImage>
void HessianVesselFilter<TInputImage, TOutputImage>::ComputeEigenSystem(
    const RealSymmetricMatrixType& symmetricMatrix, 
    RealVectorType& eigenValues, RealMatrixType& eigenVectors) const
{
    RealVectorType tempEigenValues;
    RealMatrixType tempEigenVectors;
    
    itk::SymmetricEigenAnalysis<RealSymmetricMatrixType, RealVectorType, 
        RealMatrixType> eigenSolver(itkGetStaticConstMacro(ImageDimension));
    //itk bug with sorting
    //eigenSolver.SetOrderEigenMagnitudes(true);
    unsigned int result = eigenSolver.ComputeEigenValuesAndVectors(
        symmetricMatrix, tempEigenValues, tempEigenVectors);

    //sort by hand 
    //we want |l1| <= |l2| <= |l3|        
    RealType l1 = tempEigenValues[0];
    RealType l2 = tempEigenValues[1];
    RealType l3 = tempEigenValues[2];
    
    vnl_vector<RealType> v1 = tempEigenVectors.GetVnlMatrix().get_row(0);
    vnl_vector<RealType> v2 = tempEigenVectors.GetVnlMatrix().get_row(1);
    vnl_vector<RealType> v3 = tempEigenVectors.GetVnlMatrix().get_row(2);            
    
    if ( fabs(l1) > fabs(l2)  ) 
    {
        std::swap(l1, l2);
        std::swap(v1, v2);
    }
    if ( fabs(l2) > fabs(l3)  )
    {        
        std::swap(l2, l3);
        std::swap(v2, v3);
    }
    if ( fabs(l1) > fabs(l2)  ) 
    {
        std::swap(l1, l2);
        std::swap(v1, v2);
    }
    
    vnl_matrix<RealType> sortedEigenVectors(itkGetStaticConstMacro(ImageDimension), itkGetStaticConstMacro(ImageDimension));
    sortedEigenVectors.set_row(0, v1);
    sortedEigenVectors.set_row(1, v2);
    sortedEigenVectors.set_row(2, v3);
    eigenVectors = sortedEigenVectors;
    
    eigenValues[0] = l1;
    eigenValues[1] = l2;
    eigenValues[2] = l3;    
}

template <class TInputImage, class TOutputImage>
void
HessianVesselFilter<TInputImage, TOutputImage>::ComputeEigenHessianAtIndex(typename TInputImage::IndexType index, RealType scale, RealType& vesselness, RealVectorType& vesselVector) 
{
    typename DiscreteHessianGaussianImageFunction<TInputImage, RealType>::Pointer hessianFunction = DiscreteHessianGaussianImageFunction<TInputImage, RealType>::New();
    hessianFunction->SetInputImage(this->GetInput() );
    hessianFunction->Initialize();
    hessianFunction->SetSigma(scale);
    
    RealSymmetricMatrixType hessianMatrix = hessianFunction->EvaluateAtIndex(index);
    
    RealVectorType eigenValues;
    RealMatrixType eigenVectors;
    
    this->ComputeEigenSystem(hessianMatrix, eigenValues, eigenVectors);
    
    vesselness = this->VesselnessMeasure(eigenValues[0], 
                    eigenValues[1], eigenValues[2], this->m_VesselnessAlpha, 
                    this->m_VesselnessBeta, this->m_VesselnessGamma);
                    
    vesselVector.SetVnlVector(eigenVectors.GetVnlMatrix().get_row(0));
    
}

template <class TInputImage, class TOutputImage>
void
HessianVesselFilter<TInputImage, TOutputImage>::SetScales(RealListType s)
{
    this->m_Scales = s;
} 
         
template <class TInputImage, class TOutputImage>
typename HessianVesselFilter<TInputImage, TOutputImage>::RealListType 
HessianVesselFilter<TInputImage, TOutputImage>::GetScales() const
{
    return(this->m_Scales);
}

//Print
template <class TInputImage, class TOutputImage>
void HessianVesselFilter<TInputImage,TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

//  os << indent << "foo:  " << m_foo << std::endl;

}

//foo
template <class TInputImage, class TOutputImage>
void
HessianVesselFilter<TInputImage, TOutputImage>::SetMedialFunction(MedialFunctionType* m)
{
    this->medialFunction = m;
}
 
}//end namespace                
