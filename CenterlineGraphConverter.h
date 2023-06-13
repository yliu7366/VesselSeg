#ifndef __CenterlineGraphConverter_h
#define __CenterlineGraphConverter_h

#include <iostream>
#include <vector>

#include "itkImage.h"
#include "itkVector.h"

#include "VesselGraph.h"
#include "VesselSegment.h"
#include "VesselNode.h"
#include "FixedArray.h"


/**\class CenterlineGraphConverter
 * This class simply contains methods which convert a binary centerline image
 * to a graph.
 */ 
template<typename TLabelType, typename TCoordType, unsigned int VDimension>
class CenterlineGraphConverter
{
    public:      
    
    typedef itk::Image<TLabelType, VDimension> LabelImageType;
    typedef itk::Image<TCoordType, VDimension> CoordImageType;
    typedef itk::Vector<TCoordType, VDimension> VectorType;
    typedef itk::Image< VectorType, VDimension> VectorImageType;
    typedef typename LabelImageType::IndexType::IndexValueType IDType;
    typedef VesselGraph< VDimension, TCoordType > VesselGraphType;
    
    /**A helper method to convert an image voxel index to a unique id using 
    *  the image dimensions.
    *@param image the image from which the index comes from
    *@param index the index to be transformed into an id
    *@return the id 
    */
    IDType GetIDFromIndex(typename LabelImageType::Pointer image, const typename LabelImageType::IndexType &index) const;
    
    /**Convert an image to a graph.  The graph attributes are populated with
    *  the specified images and the node positions are computed by converting
    *  the centerline voxels to physical (world) coordinates.
    *@param centerlineImage a binary image with non-zero voxels corresponding to
    *  vessel centerlines
    *@param radiusImage image whose elements which correspond with centerline
    *  voxels represent the radius of the vessel
    *@param probabilityImage image whose elements which correspond with centerline
    *  voxels represent the probability of this voxel being a vessel
    *@param tangentImage image whose elements which correspond with centerline
    *  voxels represent the tangent vector of the centerline
    *@return a vessel graph with position, radius, probability, and tangent 
    *  attributes populated
    */
    VesselGraphType ImageToGraph(typename LabelImageType::Pointer centerlineImage, typename CoordImageType::Pointer radiusImage, typename CoordImageType::Pointer probabilityImage, typename VectorImageType::Pointer tangentImage) const;
    
    /**Convert an image to a graph.  The graph attributes are computed by 
    *  transforming the specified vessel label map into a centerline using a 
    *  skeletonization algorithm, the radii are computed using a distance 
    *  using a distance transform, the probabilities are fixed and the tangents
    *  are computed using the node positions
    *@param labelMapImage a binary image with non-zero voxels corresponding to
    *  vessels
    *@return a vessel graph with position, radius, probability, and tangent 
    *  attributes populated
    */
    VesselGraphType ImageToGraph(typename LabelImageType::Pointer labelMapImage) const;
    
    /**Convert a graph to an image.   
    * @param centerlineImage a binary image with non-zero voxels corresponding to
    *  vessel centerlines.  This image is used to determine the properties of
    *  the output image.  It can be the image from with the graph was created
    *  or an image whose meta data is populated.
    * @param oversamplingFactor factor by which the image resolution is increased.
    *  For example, if the specified factor is 2 then number of voxels in 
    *  each direction increases by 2 and the image spacing is halved.
    * @param graph the graph from which the image is created
    * @return a binary label depicting the vessels.  Here the nodes radius 
    *  are used to specify the vessel width at each centerline location
    */
    typename CoordImageType::Pointer GraphToImage(typename LabelImageType::Pointer centerlineImage, int oversamplingFactor, const VesselGraphType& graph) const;

        
    private:

};



#include "CenterlineGraphConverter.txx"

#endif
