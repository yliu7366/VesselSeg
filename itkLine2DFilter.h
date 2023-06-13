#ifndef __itkLine2DFilter_h
#define __itkLine2DFilter_h


#include "itkImageToImageFilter.h"
#include "itkImage.h"

#include "itkLineTransform2DFunction.h"

namespace itk
{
/**\class itkLine2DFilter
 * 
 */
    template <typename TInputImage,typename TOutputImage >
    class ITK_EXPORT Line2DFilter 
    : public ImageToImageFilter< TInputImage, TOutputImage > 
    {
        public:
        
            //typedefs ITK
            typedef Line2DFilter Self;
            typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            
            typedef TInputImage InputImageType;
            typedef TOutputImage OuputImageType;
            
            //typedefs this class
            //define macro to get image dimension
            itkStaticConstMacro(ImageDimension, unsigned int, 
                ::itk::GetImageDimension<InputImageType>::ImageDimension);
               
            
            typedef double RealType;
            typedef itk::Image< RealType, itkGetStaticConstMacro(ImageDimension) > RealScalarImageType;
            typedef itk::Vector< RealType, itkGetStaticConstMacro(ImageDimension) > RealVectorType;
            typedef itk::Image< RealVectorType, itkGetStaticConstMacro(ImageDimension) > RealVectorImageType;
            
            typedef itk::LineTransform2DFunction<InputImageType, RealType> LineFunctionType;
            typedef itk::VectorLineTransform2DFunction<RealVectorImageType, RealType> VectorLineFunctionType;
            
                        
            //macros
            itkNewMacro(Self);

            //Getters and Setters
            //parameter in ComputeVesselness
            itkGetMacro(IntensityMaximum, RealType);
            itkSetMacro(IntensityMaximum, RealType);
            
             //parameter in ComputeVesselness
            itkGetMacro(IntensityMinimum, RealType);
            itkSetMacro(IntensityMinimum, RealType);
                       

            //Get intermediate images
            itkGetMacro(LineImage, typename RealScalarImageType::Pointer)   
            itkGetMacro(LineVectorImage, typename RealVectorImageType::Pointer);

            /**
            *Call after SetInput.  Allocates and initializes
            *all the memory for intermediates
            */
            virtual void Allocate();
         
            virtual void ComputeLineFilter();
           
            void SetLineFunction(LineFunctionType*);
            void SetVectorLineFunction(VectorLineFunctionType*);
           
   
        protected:
            Line2DFilter();
            ~Line2DFilter();
            virtual void PrintSelf(std::ostream& os, Indent indent) const;
            virtual void GenerateData();

        private:

            //purposely not implemented
            Line2DFilter(const Self&);
            
            //Intermediate outputs
            typename RealScalarImageType::Pointer m_LineImage;
            typename RealVectorImageType::Pointer m_LineVectorImage;

            LineFunctionType* lineFunction;
            VectorLineFunctionType* vectorLineFunction;
    
            //parameters
            RealType m_IntensityMaximum;
            RealType m_IntensityMinimum;


    };//end Line2DFilter class
}//end itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLine2DFilter.txx"
#endif

#endif
