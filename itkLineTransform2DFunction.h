#ifndef __itkLineTransform2DFunction_h
#define __itkLineTransform2DFunction_h

#include <utility>

#include "itkImageFunction.h"
#include "itkImageBase.h"
#include "itkInterpolateImageFunction.h"
#include "itkVector.h"
#include "itkVectorInterpolateImageFunction.h"

namespace itk
{
/**\class LineTransformFunction

 */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT LineTransform2DFunction 
        : public ImageFunction< TInputImage, Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  >, TCoordRep > 
    {
        public:
            typedef Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > VectorType;
            typedef VectorType TOutput;
            typedef LineTransform2DFunction Self;
            typedef ImageFunction< TInputImage, TOutput, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
            typedef typename Superclass::IndexType IndexType;
            typedef InterpolateImageFunction< TInputImage, TCoordRep > InterpolaterType;
            
            itkNewMacro(Self);
            
            
            /**
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
            /**
            */
            virtual TOutput EvaluateAtContinuousIndex (const ContinuousIndexType &index) const;
            
            /**
            */
            virtual TOutput EvaluateAtIndex (const IndexType &index) const;
            
            /**
            */
            void SetInterpolater(InterpolaterType*);  
            
            /**
            */
            virtual void SetInputImage(const TInputImage * ptr);
            
            /**
            */
            virtual void SetParameters(TCoordRep radius, unsigned int maxNumAngles);
            

   
   
        protected:
            LineTransform2DFunction ();
            ~LineTransform2DFunction ();
            //virtual void PrintSelf(std::ostream& os, Indent indent) const;
            
            InterpolaterType* interpolator;
            
            unsigned int maxNumAngles;
            TCoordRep radius;
            TCoordRep minSpacing;
            
            
        private:

            //purposely not implemented
            LineTransform2DFunction (const Self&);
            void operator=(const Self);

    };//end class LineTransform2DFunction 
    
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT LineTransformMinVariance2DFunction 
        : public LineTransform2DFunction< TInputImage,TCoordRep > 
    {
    
        public:
            typedef LineTransformMinVariance2DFunction Self;
            typedef LineTransform2DFunction< TInputImage, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::VectorType VectorType;
            typedef typename Superclass::TOutput TOutput;

            itkNewMacro(Self);
            
            /**
            */
            virtual TOutput Evaluate (const PointType &point) const;
    };//end class LineTransformMinVariance2DFunction
    
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT VectorLineTransform2DFunction 
        : public ImageFunction< TInputImage, std::pair<TCoordRep, Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > >, TCoordRep > 
    {
        public:
            typedef Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > VectorType;
            typedef std::pair<TCoordRep, VectorType> TOutput;
            typedef VectorLineTransform2DFunction Self;
            typedef ImageFunction< TInputImage, TOutput, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
            typedef typename Superclass::IndexType IndexType;
            typedef VectorInterpolateImageFunction< TInputImage, TCoordRep > InterpolaterType;
            
            
            itkNewMacro(Self);
            
            
            /**
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
            /**
            */
            virtual TOutput EvaluateAtContinuousIndex (const ContinuousIndexType &index) const;
            
            /**
            */
            virtual TOutput EvaluateAtIndex (const IndexType &index) const;
            
            /**
            */
            void SetInterpolater(InterpolaterType*);  
            
            /**
            */
            virtual void SetInputImage(const TInputImage * ptr);                        
            
            /**
            */
            virtual void SetParameters(TCoordRep radius, unsigned int maxNumAngles);
            

   
   
        protected:
            VectorLineTransform2DFunction ();
            ~VectorLineTransform2DFunction ();
            //virtual void PrintSelf(std::ostream& os, Indent indent) const;
            
            InterpolaterType* interpolator;
            
            unsigned int maxNumAngles;
            TCoordRep radius;
            TCoordRep minSpacing;
            
            
        private:

            //purposely not implemented
            VectorLineTransform2DFunction (const Self&);
            void operator=(const Self);

    };//end class VectorLineTransform2DFunction 
    
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT CorrelateVectorsOnLine2DFunction 
        : public ImageFunction< TInputImage, TCoordRep, TCoordRep > 
    {
        public:
            typedef TCoordRep TOutput;
            typedef CorrelateVectorsOnLine2DFunction Self;
            typedef ImageFunction< TInputImage, TOutput, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
            typedef typename Superclass::IndexType IndexType;
            typedef VectorInterpolateImageFunction< TInputImage, TCoordRep > InterpolaterType;
            typedef Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > VectorType;
            
            itkNewMacro(Self);
            
            /**
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
            /**
            */
            virtual TOutput EvaluateAtContinuousIndex (const ContinuousIndexType &index) const;
            
            /**
            */
            virtual TOutput EvaluateAtIndex (const IndexType &index) const;
            
            /**
            */
            void SetInterpolaters(InterpolaterType*, InterpolaterType*);  
        
            virtual void SetInputImages(const TInputImage * , const TInputImage * );
            virtual void SetParameters(TCoordRep radius, TCoordRep azimuth, TCoordRep elevation);
        
         protected:
            CorrelateVectorsOnLine2DFunction ();
            ~CorrelateVectorsOnLine2DFunction ();
            //virtual void PrintSelf(std::ostream& os, Indent indent) const;
            
            const TInputImage * inputImage2;
            InterpolaterType* interpolator1;
            InterpolaterType* interpolator2;
            
            TCoordRep radius;
            TCoordRep azimuth;
            TCoordRep minSpacing;
            
        private:

            //purposely not implemented
            CorrelateVectorsOnLine2DFunction (const Self&);
            void operator=(const Self);
    };//end class CorrelateVectorsOnLine2DFunction
    
}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLineTransform2DFunction.txx"
#endif

#endif
