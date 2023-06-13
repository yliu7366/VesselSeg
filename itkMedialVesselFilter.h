#ifndef __itkMedialVesselFilter_h
#define __itkMedialVesselFilter_h

#include <vector>
#include <map>
#include <list>
#include <utility>

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVector.h"
#include "itkIndex.h"
#include "itkOffset.h"


#include "itkMedialness3DFunction.h"
#include "itkHessianVesselFilter.h"


namespace itk
{
/**\class itkMedialVesselFilter
 * \brief Compute the medialness of an image and use it to generate label
 * map of vasculature
 *
 * This class uses the vesselness calculation and associated coordinate 
 * system to compute medialness of an image.  After post processing the 
 * medialness a label of the vasulcature is constructed
 */
    template <typename TInputImage,typename TOutputImage >
    class ITK_EXPORT MedialVesselFilter 
    : public ImageToImageFilter< TInputImage, TOutputImage > 
    {
        public:
        
            //typedefs ITK
            typedef MedialVesselFilter Self;
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
            typedef itk::Image< RealType, 
                itkGetStaticConstMacro(ImageDimension) > RealScalarImageType;
            typedef itk::Vector< RealType, 
                itkGetStaticConstMacro(ImageDimension) > RealVectorType;
            typedef itk::Image< RealVectorType, 
                itkGetStaticConstMacro(ImageDimension) > RealVectorImageType;
            typedef itk::Index< itkGetStaticConstMacro(ImageDimension) > 
                IndexType;
            typedef itk::Offset< itkGetStaticConstMacro(ImageDimension) > 
                OffsetType;
            typedef unsigned int LabelType;
            typedef itk::Image< LabelType, 
                itkGetStaticConstMacro(ImageDimension) > LabelImageType;     
            
            typedef HessianVesselFilter < TInputImage, TOutputImage > HessianVesselFilterType;
            typedef Medialness3DFunction < RealVectorImageType, RealType > MedialFunctionType;
            
            typedef std::list<IndexType> IndexListType; 
            
            //Label map values     
            static const LabelType s_MedialCenterlineUnlinked = 1;
            static const LabelType s_MedialSphereUnlinked = 2;
            static const LabelType s_MedialCenterlineLinked = 3;
            static const LabelType s_MedialSphereLinked = 4;
            static const LabelType s_MedialCenterlinePredicted = 5;
            static const LabelType s_MedialSpherePredicted = 6;
            
                        
            //macros
            itkNewMacro(Self);

            //Getters and Setters
            //parameter in ComputeVesselness
            itkGetMacro(IntensityMaximum, RealType);
            itkSetMacro(IntensityMaximum, RealType);
            
             //parameter in ComputeVesselness
            itkGetMacro(IntensityMinimum, RealType);
            itkSetMacro(IntensityMinimum, RealType);
                       
            //only voxels with vesselness above this threshold are consider
            //by medialness
            itkGetMacro(VesselnessThreshold, RealType);
            itkSetMacro(VesselnessThreshold, RealType);
            
            //medial threshold
            itkGetMacro(MedialnessThreshold, RealType);
            itkSetMacro(MedialnessThreshold, RealType);
            
            //medial Factor
            itkGetMacro(MedialnessPercent, RealType);
            itkSetMacro(MedialnessPercent, RealType);

            //alginment parameter when connecting centerlines
            itkGetMacro(DotProductThreshold, RealType);
            itkSetMacro(DotProductThreshold, RealType);
            
            //parameter in LinkMedialnessPlanes
            //Maximum number of prediction steps allowed
            itkGetMacro(MaxNumMedialPrediction, unsigned int);
            itkSetMacro(MaxNumMedialPrediction, unsigned int);
            
            //intensity difference between center and boundary of vessel
            itkGetMacro(IntensityDifference, RealType);
            itkSetMacro(IntensityDifference, RealType);
            
            //percentage of voxels in medial plane at vessel boundary
            // that are above intensity difference
            itkGetMacro(EdgeRatio, RealType);
            itkSetMacro(EdgeRatio, RealType);

            //Get intermediate images
            itkGetMacro(MedialnessImage, typename RealScalarImageType::Pointer);
            itkGetMacro(RadiusImage, typename RealScalarImageType::Pointer);
            itkGetMacro(CenterlineImage, typename LabelImageType::Pointer);
            itkGetMacro(CenterlineSphereImage, typename LabelImageType::Pointer);
            
            itkGetMacro(LineCenterImage, typename RealScalarImageType::Pointer);
            itkGetMacro(LinePredictImage, typename RealScalarImageType::Pointer);
            

            /**
            * Get the voxel indices in 3x3x3 neighborhood that are in the
            * direction of the specified vector
            * @param center image location
            * @param vector the desired direction
            * @param angle the allowable angle from which we can deviate 
            *   from the desired direction
            * @return a list of indices that are within the cone specified
            *   by the vector and angle
            */
            virtual IndexListType GetInDirectionIndices(
                const IndexType& center,
                const RealVectorType& vector,
                const RealType& angle) const;
                
            /**
            * Get the voxel indices in 3x3x3 neighborhood that are in the
            * plane specified by the vectors
            * @param center image location
            * @param vectorPlane1 vector describing plane
            * @param vectorPlane2 vector describing plane
            * @return a list of indices that are in the plane
            */
            virtual IndexListType GetInPlaneIndices(
                const IndexType& center,
                const RealVectorType& vectorPlane1, 
                const RealVectorType& vectorPlane2) const;
           
                
            /**
            * Create a label of the segmented vasculature using the centerline
            * and radius at centerline location to draw a sphere
            * @param centerlineImage binary image depicting centerline voxels
            * @param centerlineLabel value of centerline voxels in label map 
            * @param radiusImage image with the radius of the vessel at each voxel
            * @param radiusLabel value of non-centerline voxels in label map
            * @return image with vasculature labeled
            */
            virtual typename LabelImageType::Pointer LabelSphereOnCenterline(
                typename LabelImageType::Pointer centerlineImage, LabelType centerlineLabel, 
                typename RealScalarImageType::Pointer radiusImage, LabelType radiusLabel) const;
                
            /**
            * Label the specified location on the specified image with a sphere
            * of the specified radius
            * @param index image location of sphere origin
            * @param labelImage image to draw sphere
            * @param centerlineLabel value of centerline voxels in label map 
            * @param radius radius of sphere
            * @param radiusLabel value of non-centerline voxels in label map
            */
            virtual void LabelSphereOnCenterlineAtIndex(
                const IndexType& index,
                typename LabelImageType::Pointer labelImage, LabelType centerlineLabel, 
                RealType radius, LabelType radiusLabel) const;
           
            
            /**
            * Convert centerline to binary image
            * @return image with 1 denoting centerline and 0 otherwise
            */
            virtual typename LabelImageType::Pointer GetCenterlineAsBinary() const;
            
            
            /**
            *Call after SetInput.  Allocates and initializes
            *all the memory for intermediates
            */
            virtual void Allocate();
            
            
            
            /**
            *Link medial planes by looking in the neighborhood of each center
            *point.  If there are other centerpoints around we re-labed
            *the vessel planes in this->m_LabeledPlanesImage with centerpoints
            *labeled with s_MedialCenterlineLinked
            *and voxels in plane given a value of s_MedialPlaneLinked.
            */
            virtual void LinkMedialnessPlanes();

           
            
            /**
            *Predict centerlines by using the vessel vector to predict candidate
            *centerline voxels. If enough predicted adjacent voxels are in plane
            *medialness maximum then add them to this->m_CenterlineImage.  The
            *alogirhtm is parameters are this->MedialnessThreshold,
            *this->m_MedialnessPercent, this->m_MaxNumMedialPrediction,   
            *this->m_NumMustBeMedialConnected.
            */
            virtual void PredictCenterlineFromMedialness();
            
            /**
            *Compute estimated centerlines using medialness and local
            *coordinate system.  An estimated centerline is defined by
            *a local maximum of medialness in the cross-sectional plane.
            *The results is stored in this->m_CenterlineImage
            */
            virtual void ComputeCenterlineFromMedialness();
            
            /**Compute medialness of all voxels whose have a neighborhood
            *of voxels with vesselness above the threshold 
            *this->m_VesselnessThreshold. Medialness is calculated in the
            *plane defined by the vectors in this->m_VesselPlane1Image and
            *this->m_VesselPlane2Image and at the radii specified in the list
            *this->m_Radii.  We store the maximum medialness (over radii) in
            *this->m_MedialnessImage and the radii it occur at in
            *this->m_RadiusImage.
            * /f[
            * M_{r}(\mathbf{x})=\frac{1}{N}\sum_{i=0}^{N-1}\min\left\{ -\nabla I(\mathbf{x}+r\mathbf{v}_{\alpha_{i}})\cdot\mathbf{v}_{\alpha_{i}},\,\,-\nabla I(\mathbf{x}+r\mathbf{v}_{\alpha_{i}+\pi})\cdot\mathbf{v}_{\alpha_{i}+\pi}\right\}
            * \alpha_{i}=\pi i/N$,~~ $i=0,\dots,N-1$,~~$N=\left\lceil 2\pi r+1\right\rceil 
            * r^{*}=\arg\max_{r}\, M_{r}(\mathbf{x})
            * \mathcal{M}_{r^{*}}(\mathbf{x})=\frac{1}{N}\sum_{i=0}^{N-1}\min\left\{ \frac{-\nabla I(\mathbf{x}+r^{*}\mathbf{v}_{\alpha_{i}})}{\left\Vert \nabla I(\mathbf{x}+r^{*}\mathbf{v}_{\alpha_{i}})\right\Vert }\cdot\mathbf{v}_{\alpha_{i}},\,\,\frac{-\nabla I(\mathbf{x}+r^{*}\mathbf{v}_{\alpha_{i}+\pi})}{\left\Vert \nabla I(\mathbf{x}+r^{*}\mathbf{v}_{\alpha_{i}+\pi})\right\Vert }\cdot\mathbf{v}_{\alpha_{i}+\pi}\right\}
            * /f] 
            */
            virtual void ComputeMedialness();
             
            /**
            *  Set the medial function
            */
           void SetVesselFilter(HessianVesselFilterType*);
           
           /**
           *  Set the medial function
           */
           void SetMedialFunction(MedialFunctionType*);
   
           
   
        protected:
            MedialVesselFilter();
            ~MedialVesselFilter();
            virtual void PrintSelf(std::ostream& os, Indent indent) const;
            virtual void GenerateData();

        private:

            //purposely not implemented
            MedialVesselFilter(const Self&);
            
            //Intermediate outputs
            typename RealScalarImageType::Pointer m_MedialnessImage;
            typename RealScalarImageType::Pointer m_RadiusImage;
            typename LabelImageType::Pointer m_CenterlineImage;
            typename LabelImageType::Pointer m_CenterlineSphereImage;
            
            typename RealScalarImageType::Pointer m_LineCenterImage;
            typename RealScalarImageType::Pointer m_LinePredictImage;
            
            HessianVesselFilterType* hessianVesselFilter;
            MedialFunctionType* medialFunction;
          
            
            //parameters
            RealType m_IntensityMaximum;
            RealType m_IntensityMinimum;
            RealType m_VesselnessThreshold;
            RealType m_MedialnessThreshold;
            RealType m_MedialnessPercent;
            RealType m_DotProductThreshold;
            unsigned int m_MaxNumMedialPrediction;
            RealType m_IntensityDifference;
            RealType m_EdgeRatio;

    };//end MedialVesselFilter class
}//end itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMedialVesselFilter.txx"
#endif

#endif
