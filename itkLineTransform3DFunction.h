#ifndef __itkLineTransform3DFunction_h
#define __itkLineTransform3DFunction_h

#include <utility>

#include "itkImageFunction.h"
#include "itkImageBase.h"
#include "itkInterpolateImageFunction.h"
#include "itkVector.h"
#include "itkVectorInterpolateImageFunction.h"

namespace itk
{
    /**\class LineTransformFunction
     * \brief A function which computes the line transform given a scalar volume 
     *
     * The line transform integrates image values along a line of specified
     * length and orientation.  The output is obtained when the maximum 
     * over all orientations is obtained.   
     */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT LineTransform3DFunction 
        : public ImageFunction< TInputImage, Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  >, TCoordRep > 
    {
        public:
            typedef Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > VectorType;
            typedef VectorType TOutput;
            typedef LineTransform3DFunction Self;
            typedef ImageFunction< TInputImage, TOutput, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
            typedef typename Superclass::IndexType IndexType;
            typedef InterpolateImageFunction< TInputImage, TCoordRep > InterpolaterType;
            
            itkNewMacro(Self);
            
            
            /**
            *  Evaluate line transform
            *  @param point physical location to compute line transform
            *  @return a vector pointing in the direction of the largest
            *    response and whose norm is the line transform value in [0,1]
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
            /**
            *  Evaluate line transform
            *  @param index image location to compute line transform
            *  @return a vector pointing in the direction of the largest
            *    response and whose norm is the line transform value in [0,1]
            */
            virtual TOutput EvaluateAtContinuousIndex (const ContinuousIndexType &index) const;
            
            /**
            *  Evaluate line transform
            *  @param index image location to compute line transform
            *  @return a vector pointing in the direction of the largest
            *    response and whose norm is the line transform value in [0,1]
            */
            virtual TOutput EvaluateAtIndex (const IndexType &index) const;
            
            /**
            * Set the type of interpolation used on input image
            */
            void SetInterpolater(InterpolaterType*);  
            
            /**
            * Set input image
            */
            virtual void SetInputImage(const TInputImage * ptr);
            
            /**
            * Set parameters used in line transform
            * @param radius half of line length (in units of image spacing)
            * @param maxNumAngles maximum number of angles to use
            */
            virtual void SetParameters(TCoordRep radius, unsigned int maxNumAngles);
            

   
   
        protected:
            LineTransform3DFunction ();
            ~LineTransform3DFunction ();
            //virtual void PrintSelf(std::ostream& os, Indent indent) const;
            
            InterpolaterType* interpolator;
            
            unsigned int maxNumAngles;
            TCoordRep radius;
            TCoordRep minSpacing;
            
            
        private:

            //purposely not implemented
            LineTransform3DFunction (const Self&);
            void operator=(const Self);

    };//end class LineTransform3DFunction 
    
    /**\class LineTransformFunction
     * \brief A function which computes the line transform given a scalar volume 
     *
     * The line transform considers image values along a line of specified
     * length and orientation.  The output is obtained when the mininum variance 
     * over all orientations is obtained.   
     */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT LineTransformMinVariance3DFunction 
        : public LineTransform3DFunction< TInputImage,TCoordRep > 
    {
    
        public:
            typedef LineTransformMinVariance3DFunction Self;
            typedef LineTransform3DFunction< TInputImage, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::VectorType VectorType;
            typedef typename Superclass::TOutput TOutput;

            itkNewMacro(Self);
            
            /**
            *  Evaluate line transform
            *  @param point physical location to compute line transform
            *  @return a vector pointing in the direction of the largest
            *    response and whose norm is the line transform value in [0,1]
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
    };//end class LineTransformMinVariance3DFunction
    
    
    /**\class VectorLineTransform3DFunction
     * \brief A function which computes the line transform given a vector volume 
     *
     * The line transform integrates image values along a line of specified
     * length and orientation.  The output is obtained when the maximum 
     * over all orientations is obtained.   
     */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT VectorLineTransform3DFunction 
        : public ImageFunction< TInputImage, std::pair<TCoordRep, Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > >, TCoordRep > 
    {
        public:
            typedef Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > VectorType;
            typedef std::pair<TCoordRep, VectorType> TOutput;
            typedef VectorLineTransform3DFunction Self;
            typedef ImageFunction< TInputImage, TOutput, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
            typedef typename Superclass::IndexType IndexType;
            typedef VectorInterpolateImageFunction< TInputImage, TCoordRep > InterpolaterType;
            
            
            itkNewMacro(Self);
            
            
            /**
            *  Evaluate line transform
            *  @param point physical location to compute line transform
            *  @return the orientation that correspond to the direction when
            *    the vector field is maximally aligned and the magnitude
            *    in that direction
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
            /**
            *  Evaluate line transform
            *  @param index image location to compute line transform
            *  @return the orientation that correspond to the direction when
            *    the vector field is maximally aligned and the magnitude
            *    in that direction
            */
            virtual TOutput EvaluateAtContinuousIndex (const ContinuousIndexType &index) const;
            
            /**
            *  Evaluate line transform
            *  @param index image location to compute line transform
            *  @return the orientation that correspond to the direction when
            *    the vector field is maximally aligned and the magnitude
            *    in that direction
            */
            virtual TOutput EvaluateAtIndex (const IndexType &index) const;
            
            /**
            * Set the type of interpolation used on input image
            */
            void SetInterpolater(InterpolaterType*);  
            
            /**
            * Set input image
            */
            virtual void SetInputImage(const TInputImage * ptr);                        
            
            /**
            * Set parameters used in line transform
            * @param radius half of line length (in units of image spacing)
            * @param maxNumAngles maximum number of angles to use
            */
            virtual void SetParameters(TCoordRep radius, unsigned int maxNumAngles);

        protected:
            VectorLineTransform3DFunction ();
            ~VectorLineTransform3DFunction ();
            //virtual void PrintSelf(std::ostream& os, Indent indent) const;
            
            InterpolaterType* interpolator;
            
            unsigned int maxNumAngles;
            TCoordRep radius;
            TCoordRep minSpacing;
                    
        private:

            //purposely not implemented
            VectorLineTransform3DFunction (const Self&);
            void operator=(const Self);

    };//end class VectorLineTransform3DFunction 
    
    /**\class CorrelateVectorsOnLine3DFunction
     * \brief A function which computes the correlation between two vector volumes  
     */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT CorrelateVectorsOnLine3DFunction 
        : public ImageFunction< TInputImage, TCoordRep, TCoordRep > 
    {
        public:
            typedef TCoordRep TOutput;
            typedef CorrelateVectorsOnLine3DFunction Self;
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
            * Find the correlation between two vector fields using the
            * specified orientation and distance
            *  @param point physical location to compute correlation
            *  @return correlation value in [0,1]
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
            /**
            * Find the correlation between two vector fields using the
            * specified orientation and distance
            *  @param index image location to compute correlation
            *  @return correlation value in [0,1]
            */
            virtual TOutput EvaluateAtContinuousIndex (const ContinuousIndexType &index) const;
            
            /**
            * Find the correlation between two vector fields using the
            * specified orientation and distance
            *  @param index image location to compute correlation
            *  @return correlation value in [0,1]
            */
            virtual TOutput EvaluateAtIndex (const IndexType &index) const;
            
            /**
            * Set the type of interpolation used on input image
            */
            void SetInterpolaters(InterpolaterType*, InterpolaterType*); 
             
            /**
            * Set input image
            */
            virtual void SetInputImages(const TInputImage * , const TInputImage * );
            
            /**
            * Set parameters used in correlation
            * @param radius half of line length (in units of image spacing)
            * @param azimuth azimuthal angle
            * @param elevation elevation angle
            */
            virtual void SetParameters(TCoordRep radius, TCoordRep azimuth, TCoordRep elevation);
        
         protected:
            CorrelateVectorsOnLine3DFunction ();
            ~CorrelateVectorsOnLine3DFunction ();
            //virtual void PrintSelf(std::ostream& os, Indent indent) const;
            
            const TInputImage * inputImage2;
            InterpolaterType* interpolator1;
            InterpolaterType* interpolator2;
            
            TCoordRep radius;
            TCoordRep azimuth;
            TCoordRep elevation;
            TCoordRep minSpacing;
            
        private:

            //purposely not implemented
            CorrelateVectorsOnLine3DFunction (const Self&);
            void operator=(const Self);
    };//end class CorrelateVectorsOnLine3DFunction
    
}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLineTransform3DFunction.txx"
#endif

#endif
