#ifndef __itkMedialness3DFunction_h
#define __itkMedialness3DFunction_h

#include <vector>
#include <utility>

#include "itkImageFunction.h"
#include "itkImageBase.h"
#include "itkVector.h"
#include "itkVectorInterpolateImageFunction.h"

namespace itk
{
/**\class Medialness3DFunction
 * \brief A function which computes the medialness given a volume of three
 * dimension gradient vectors.
 * 
 * Medialness is a value from [0,1] which describes the about of "circleness"
 * in a given plane.  The output of this function is a pair of numbers, 
 * the first representing the largest medialness over the user specified 
 * radius range and the radius at which it occured.
 * The user may call either EvaluateWithPlane with the desired point and 
 * vectors which define the cross sectional plane or EvaluateWithNormal with the
 * planes normal vector.  The user also may call Evaluate, 
 * EvaluateAtContinuousIndex, or EvaluateAtIndex after setting the plane normal
 * vector at the specified point or index.  The computation of medialness
 * occurs in physical space and is thus appropriate even if voxel spacing 
 * is anisotropic. To override the functionality of this class simply override 
 * the method EvaluateWithPlane.  Since interpolation of the gradient vectors
 * is required the user may set the type of interpolation.
 *
 */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT Medialness3DFunction 
        : public ImageFunction< TInputImage, std::pair<TCoordRep, TCoordRep>, TCoordRep > 
    {
        public:
        
            typedef std::pair<TCoordRep, TCoordRep> TOutput; 
            typedef Medialness3DFunction Self;
            typedef ImageFunction< TInputImage, TOutput, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
            typedef typename Superclass::IndexType IndexType;
            typedef VectorInterpolateImageFunction< TInputImage, TCoordRep > VectorInterpolaterType;
            typedef Vector< TCoordRep, GetImageDimension<TInputImage>::ImageDimension  > VectorType;
            typedef std::vector< TCoordRep > ListType;
            
            itkNewMacro(Self);
            
            /**
            * Compute medialness at the specified point in the plane defined by
            * the specified unit vectors
            * @param point physical location to compute medialness
            * @param inPlane1 vector which defines medial plane
            * @param inPlane2 vector which defines medial plane
            * @return medialness and radius
            */
            virtual TOutput EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const;
            
            /**
            * Compute medialness at the specified point in the plane defined by
            * the specified unit normal vector
            * @param point physical location to compute medialness
            * @param normal vector which defines medial plane normal
            * @return medialness and radius
            */
            virtual TOutput EvaluateWithNormal(const PointType &point, const VectorType &normal);
            
            /**
            * Compute medialness at the specified point in the plane defined by
            * the normal vector set by calling SetNormal
            * @param point physical location to compute medialness
            * @return medialness and radius
            */
            virtual TOutput Evaluate (const PointType &point) const;
            
            /**
            * Compute medialness at the specified continuous index in the plane 
            * defined by the normal vector set by calling SetNormal
            * @index point image location to compute medialness
            * @return medialness and radius
            */
            virtual TOutput EvaluateAtContinuousIndex (const ContinuousIndexType &index) const;
            
            /**
            * Compute medialness at the specified index in the plane 
            * defined by the normal vector set by calling SetNormal
            * @index point image location to compute medialness
            * @return medialness and radius
            */
            virtual TOutput EvaluateAtIndex (const IndexType &index) const;
            
            /**
            * Set the unit normal vector.
            */
            void virtual SetNormal(const VectorType& normal);
            
            /**
            * Compute the unit vectors within the plane defined by the normal 
            * vector
            * @param inPlane1 vector which defines medial plane
            * @param inPlane2 vector which defines medial plane
            */
            void virtual ComputeInPlaneVector(VectorType& inPlane1, VectorType& inPlane2) const;
            
            /**
            * Set the type of interpolation used on the gradient vectors
            */
            void SetInterpolater(VectorInterpolaterType*);  
            
            /**
            * Set the gradient field
            */
            virtual void SetInputImage(const TInputImage * ptr);
            
            /**
            * Set radii for which to compute the medialness.  Also specify 
            * the maximum number of angles to use.
            */
            void SetRadii(TCoordRep start, TCoordRep stop, unsigned int numRadii, unsigned int maxNumAngles);
            
            /**
            * Set radii for which to compute the medialness.  Also specify 
            * the maximum number of angles to use.
            */
            void SetRadii(ListType radii,  unsigned int maxNumAngles);
   
   
        protected:
            Medialness3DFunction ();
            ~Medialness3DFunction ();
            virtual void PrintSelf(std::ostream& os, Indent indent) const;
            
            VectorInterpolaterType* gradientInterpolator;
            TCoordRep minSpacing;
            TCoordRep maxSpacing;
            unsigned int maxNumAngles;
            ListType radii;
            VectorType normal;
            

        private:

            //purposely not implemented
            Medialness3DFunction (const Self&);
            void operator=(const Self);

    };//end class Medialness3DFunction 
    
    /**\class MedialnessMin3DFunction
    * \brief A function which computes the medialness given a volume of three
    * dimension gradient vectors.
    * Here the medialness is computed by taking the minimum of opposing 
    * gradient vectors.
    */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT MedialnessMin3DFunction 
        : public Medialness3DFunction< TInputImage, TCoordRep > 
    {
        public:
         
            typedef MedialnessMin3DFunction Self;
            typedef Medialness3DFunction< TInputImage, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::TOutput TOutput; 
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::VectorType VectorType;
     
            itkNewMacro(Self);
            
            /**
            * Compute medialness at the specified point in the plane defined by
            * the specified unit vectors
            * @param point physical location to compute medialness
            * @param inPlane1 vector which defines medial plane
            * @param inPlane2 vector which defines medial plane
            * @return medialness and radius
            */
            virtual TOutput EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const;
        
        protected:
            MedialnessMin3DFunction ();
            ~MedialnessMin3DFunction ();
            
    };//end class
    
    /**\class MedialnessVarRadius3DFunction 
    * \brief A function which computes the medialness given a volume of three
    * dimension gradient vectors.
    * Here the medialness is computed by looping through angles first.  So for
    * each angle we get a medialness associated with each radius.  We determine
    * the medialness by getting the maximum medialness at each radius
    * and multipling by the variance in the radius
    */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT MedialnessVarRadius3DFunction 
        : public Medialness3DFunction< TInputImage, TCoordRep > 
    {
        public:
         
            typedef MedialnessVarRadius3DFunction Self;
            typedef Medialness3DFunction< TInputImage, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::TOutput TOutput; 
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::VectorType VectorType;
            
            itkNewMacro(Self);
            
            /**
            * Compute medialness at the specified point in the plane defined by
            * the specified unit vectors
            * @param point physical location to compute medialness
            * @param inPlane1 vector which defines medial plane
            * @param inPlane2 vector which defines medial plane
            * @return medialness and radius
            */
            virtual TOutput EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const;
        
        protected:
            MedialnessVarRadius3DFunction ();
            ~MedialnessVarRadius3DFunction ();
            
    };//end class
    
    /**\class MedialnessVarRadius3DFunction 
    * \brief A function which computes the medialness given a volume of three
    * dimension gradient vectors.
    * Here the medialness is computed by computing the dot product between
    * the gradient vectors and the vector normal to the medial plane (in the
    * direction of the vessel)
    */
    template < typename TInputImage, typename TCoordRep >
    class ITK_EXPORT MedialnessNormal3DFunction 
        : public Medialness3DFunction< TInputImage, TCoordRep > 
    {
        public:
         
            typedef MedialnessNormal3DFunction Self;
            typedef Medialness3DFunction< TInputImage, TCoordRep > Superclass;
            typedef SmartPointer< Self > Pointer;
            typedef SmartPointer< const Self > ConstPointer;
            typedef typename Superclass::TOutput TOutput; 
            typedef typename Superclass::PointType PointType;
            typedef typename Superclass::VectorType VectorType;
     
            itkNewMacro(Self);
            
            /**
            * Compute medialness at the specified point in the plane defined by
            * the specified unit vectors
            * @param point physical location to compute medialness
            * @param inPlane1 vector which defines medial plane
            * @param inPlane2 vector which defines medial plane
            * @return medialness and radius
            */
            virtual TOutput EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const;
        
        protected:
            MedialnessNormal3DFunction ();
            ~MedialnessNormal3DFunction ();
            
    };//end class
    
    
}//end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMedialness3DFunction.txx"
#endif

#endif
