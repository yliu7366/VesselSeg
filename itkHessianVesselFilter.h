#ifndef __itkHessianVesselFilter_h
#define __itkHessianVesselFilter_h

#include <vector>

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkMatrix.h"
#include "itkSymmetricSecondRankTensor.h"
//foo
#include "itkMedialness3DFunction.h"


namespace itk
{
/**\class itkHessianVesselFilter
 *  \brief Compute the vesselnes of an image
 *
 *  Use the eigenvalues of the multi Hessian matrix to compute a vesselness
 *  measure.  For each voxel we compute the vesselness, the scale, and the
 *  direction of the vessel
 */
    template <typename TInputImage,typename TOutputImage >
    class ITK_EXPORT HessianVesselFilter 
    : public ImageToImageFilter< TInputImage, TOutputImage > 
    {
        public:
        
            //typedefs ITK
            typedef HessianVesselFilter Self;
            typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            
            typedef TInputImage InputImageType;
            typedef TOutputImage OuputImageType;
            
            //typedefs this class
            //define macro to get image dimension
            itkStaticConstMacro(ImageDimension, unsigned int, 
                ::itk::GetImageDimension<InputImageType>::ImageDimension);
            
            //typedefs for data
            typedef double RealType;  
            typedef std::vector<RealType> RealListType; 
            
            typedef itk::Image< RealType, 
                itkGetStaticConstMacro(ImageDimension) > RealScalarImageType;
            typedef itk::Vector< RealType, 
                itkGetStaticConstMacro(ImageDimension) > RealVectorType;
            typedef itk::Image< RealVectorType, 
                itkGetStaticConstMacro(ImageDimension) > RealVectorImageType;
            typedef itk::Matrix< RealType, 
                itkGetStaticConstMacro(ImageDimension),  
                itkGetStaticConstMacro(ImageDimension) > RealMatrixType;               
            typedef itk::SymmetricSecondRankTensor< RealType, 
                itkGetStaticConstMacro(ImageDimension) > 
                RealSymmetricMatrixType;
            typedef itk::Image< RealSymmetricMatrixType, 
                itkGetStaticConstMacro(ImageDimension) > 
                RealSymmetricMatrixImageType;
     
            typedef Medialness3DFunction < RealVectorImageType, RealType > MedialFunctionType;
                        
            //macros
            itkNewMacro(Self);

            //Getter and Setter macros
            //parameter in ComputeVesselness
            itkGetMacro(IntensityMaximum, RealType);
            itkSetMacro(IntensityMaximum, RealType);
            
             //parameter in ComputeVesselness
            itkGetMacro(IntensityMinimum, RealType);
            itkSetMacro(IntensityMinimum, RealType);
                       
            //parameter in ComputeVesselness
            itkGetMacro(VesselnessAlpha, RealType);
            itkSetMacro(VesselnessAlpha, RealType);
            
            //parameter in ComputeVesselness
            itkGetMacro(VesselnessBeta, RealType);
            itkSetMacro(VesselnessBeta, RealType);  
            
            //parameter in ComputeVesselness
            itkGetMacro(VesselnessGamma, RealType);
            itkSetMacro(VesselnessGamma, RealType); 
            
            //intensity difference between center and boundary of vessel
            itkGetMacro(IntensityDifference, RealType);
            itkSetMacro(IntensityDifference, RealType);
            
            //percentage of voxels in medial plane at vessel boundary
            // that are above intensity difference
            itkGetMacro(EdgeRatio, RealType);
            itkSetMacro(EdgeRatio, RealType);
            
            //set the scales used (in units of image spacing)
            void SetScales(RealListType);
            RealListType GetScales() const;
            
            //set the medialness function
            void SetMedialFunction(MedialFunctionType*);

            //Get intermediate images
            itkGetMacro(VesselnessImage, typename RealScalarImageType::Pointer);
            itkGetMacro(ScaleImage, typename RealScalarImageType::Pointer);
            itkGetMacro(VesselMedialImage, typename RealScalarImageType::Pointer);
            itkGetMacro(VesselVectorImage, typename RealVectorImageType::Pointer);
     
            /**
            *  Allocate memory used by filter. The memory usage is 
            * 6 RealScalarImageType
            */
            void Allocate();
            
            /**Compute the maximum vesselness response over given scales.
            *Each voxel of the input image gets assigned a vesselness response
            *which is stored in this->m_VesselnessImage.  Furthermore a local
            *coordinates system (three vectors) is assigned to each voxel.  
            *The vectors have Euclidean norm (2-norm) 1 or they are the zero
            *vector. The vector which points in the vessel direction is stored
            *in this->m_VesselVectorImage and the vectors which define the 
            *vessel plane are stored in this->m_VesselPlane1Image and
            *this->m_VesselPlane2Image.
            */
            virtual void ComputeVesselness();
            
            /**Must be called after this->ComputeVesselness.  For all non-zero
            * vesselness voxels we compute the medialness of that voxel by 
            * using the vesselness scale as a guide for the radius range used
            * in medialness calculation.  Furthermore, we determine if the 
            * contrast between the voxel and those at the radius which maximizes
            * the medialness is above this->_IntensityDifference.  For the voxels
            * whose percentage of voxels in the medial plane is above the
            * the this->m_EdgeRatio we record the vesselness in 
            * this->m_VesselMedialImage.  The idea here is that voxels near 
            * center of the vessel will have high contrast with those just outside
            * the vessel boundary.
            */
            virtual void ComputeMedialVesselness();
            
            /**Compute a vesselness measure from the eigenvalues of the Hessian
            *matrix. This measure should be a number between 0 and 1 and can
            *be interpretted as a probability.
            * /f[
            * |\lambda_{1}|\le|\lambda_{2}|\le|\lambda_{3}|
            * \left(1-e^{-\frac{\lambda_{2}^{2}}{2\alpha^{2}\lambda_{3}^{2}}}\right)e^{-\frac{\lambda_{1}^{2}}{2\beta^{2}|\lambda_{2}\lambda_{3}|}}
            * /f]  
            *
            *@param lambda1 smallest (in absolute value) eigenvalue of Hessian
            *@param lambda2 middle (in absolute value) eigenvalue of Hessian
            *@param lambda3 largest (in absolute value) eigenvalue of Hessian
            *@param alpha standard deviation for middle/largest ratio term
            *@param beta standard deviation for other ratio term
            *@param gamma standard deviation for magnitude term
            *@return the vesselness measure between 0 and 1
            */
            virtual RealType VesselnessMeasure(RealType lambda1, 
                RealType lambda2, RealType lambda3, RealType alpha, 
                RealType beta, RealType gamma) const;
            
            
            /**Compute eigenvalues and eigenvectors of a real symmetric matrix
            *with dimension given by input image dimension.
            *
            *@param symmetricMatrix (input) real symmetric matrix
            *@param eigenValues (output) Eigenvalues are ordered from smallest 
            *to largest in absolute value.
            *@param eigenVectors (output) matrix containing eigenvectors as row. 
            */
            virtual void ComputeEigenSystem(
                const RealSymmetricMatrixType& symmetricMatrix, 
                RealVectorType& eigenValues, RealMatrixType& eigenVectors) 
                const;
            
            /**Compute the vesselness given the desired
            * scale
            * @param index voxel coordinate of input image for which 
            *  we compute the hessian
            * @param scale scale in units of image spacing
            * @param vesselness vesselness at desired location
            * @param vesselVector vector that points in the direction
            *  of vessel (eigenvector associated with smallest eigenvalue)
            */           
            virtual void ComputeEigenHessianAtIndex(typename TInputImage::IndexType index, RealType scale, RealType& vesselness, RealVectorType& vesselVector);   
            

        protected:
            HessianVesselFilter();
            ~HessianVesselFilter();
            virtual void PrintSelf(std::ostream& os, Indent indent) const;
            virtual void GenerateData();

        private:

            //purposely not implemented
            HessianVesselFilter(const Self&);
            
            //Intermediate outputs
            typename RealScalarImageType::Pointer m_VesselnessImage;
            typename RealScalarImageType::Pointer m_VesselMedialImage;
            typename RealScalarImageType::Pointer m_ScaleImage;
            typename RealVectorImageType::Pointer m_VesselVectorImage;
          
            
            //parameters
            RealType m_IntensityMaximum;
            RealType m_IntensityMinimum;
            RealType m_IntensityDifference;
            RealType m_EdgeRatio;
            RealListType m_Scales; 
            
            RealType m_VesselnessAlpha;
            RealType m_VesselnessBeta;
            RealType m_VesselnessGamma;
            
           
            MedialFunctionType* medialFunction;

    };//end HessianVesselFilter class
}//end itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHessianVesselFilter.txx"
#endif

#endif
