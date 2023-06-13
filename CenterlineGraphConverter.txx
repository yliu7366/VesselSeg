
#include "CenterlineGraphConverter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkPoint.h"
#include "itkSphereSpatialFunction.h"
#include "itkResampleImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
//not in ITK 
#include "itkBinaryThinningImageFilter3D.h"

template<typename TLabelType, typename TCoordType, unsigned int VDimension>
typename CenterlineGraphConverter<TLabelType, TCoordType, VDimension>::IDType
CenterlineGraphConverter<TLabelType, TCoordType, VDimension>::GetIDFromIndex( typename LabelImageType::Pointer image, const typename LabelImageType::IndexType &index) const
{
    typename LabelImageType::SizeType size = image->GetLargestPossibleRegion().GetSize();
    typename LabelImageType::IndexType::IndexValueType id = index[0];
    for(unsigned int i = 1; i < VDimension; i++)
    {
        typename LabelImageType::SizeValueType multiplyFactor = size[0];
        for(unsigned int j = 1; j < i; j++)
        {
            multiplyFactor *= size[j];
        }    
        
        id += multiplyFactor*index[i];
    }

    return(id);
}

template<typename TLabelType, typename TCoordType, unsigned int VDimension>
typename CenterlineGraphConverter<TLabelType, TCoordType, VDimension>::VesselGraphType 
CenterlineGraphConverter<TLabelType, TCoordType, VDimension>::ImageToGraph(typename LabelImageType::Pointer centerlineImage, typename CoordImageType::Pointer radiusImage, typename CoordImageType::Pointer probabilityImage, typename VectorImageType::Pointer tangentImage) const
{

    //neighborhood iterator with zero boundary condition
     typename itk::ConstNeighborhoodIterator< LabelImageType >::RadiusType 
        neighborhoodRadius;
    unsigned int numNeighborsEachDirection = 1;
    unsigned int numNeighbors = (2*numNeighborsEachDirection+1)*
        (2*numNeighborsEachDirection+1)*(2*numNeighborsEachDirection+1);
    neighborhoodRadius.Fill(numNeighborsEachDirection);
    typename itk::ConstNeighborhoodIterator< LabelImageType > 
        centerlineNeighIter(neighborhoodRadius, centerlineImage, 
        centerlineImage->GetRequestedRegion() );
    itk::ConstantBoundaryCondition<LabelImageType> zeroBC;
    zeroBC.SetConstant(0);
    centerlineNeighIter.OverrideBoundaryCondition(&zeroBC);    

    VesselGraphType vesselGraph;
 
    for(centerlineNeighIter.GoToBegin(); !centerlineNeighIter.IsAtEnd(); ++centerlineNeighIter)
    {
        typename LabelImageType::PixelType currentCenterline = centerlineNeighIter.GetCenterPixel();
        typename LabelImageType::IndexType currentIndex = centerlineNeighIter.GetIndex();
        
        if(currentCenterline != 0)
        {           
            //create node
            typename VesselGraphType::NodeType currentNode;
            typename VesselGraphType::NodeType::IDType currentID = this->GetIDFromIndex(centerlineImage, currentIndex);
            currentNode.SetID(currentID);
            
            //physical location of node
            itk::Point<TCoordType, VDimension> currentPhysicalPoint;
            centerlineImage->TransformIndexToPhysicalPoint(currentIndex, currentPhysicalPoint);
            typename VesselGraphType::NodeType::PositionType currentPosition(currentPhysicalPoint.GetDataPointer());
            currentNode.SetPosition(currentPosition);
            
            //node properties
            currentNode.SetRadius(radiusImage->GetPixel(currentIndex));
            currentNode.SetMedialness(probabilityImage->GetPixel(currentIndex));
            typename VesselGraphType::NodeType::VectorType currentTangent(tangentImage->GetPixel(currentIndex).GetDataPointer());
            currentNode.SetTangent(currentTangent);
                
            for (unsigned int i = 0; i < numNeighbors; i++)
            {
                typename LabelImageType::PixelType neighCenterline = centerlineNeighIter.GetPixel(i);

                if(neighCenterline != 0)
                {
                    typename LabelImageType::IndexType neighIndex = centerlineNeighIter.GetIndex(i);
                    if(currentIndex != neighIndex)
                    {
                        typename VesselGraphType::NodeType::IDType neighID = this->GetIDFromIndex(centerlineImage, neighIndex);
                        currentNode.AddNeighbor(neighID);                          
                    }
                }
            }
            
            vesselGraph.AddNode(currentNode);
        }//end centerline condition
    }//end neighborhood loop

    vesselGraph.FormGraph();
    
    return(vesselGraph);
}

template<typename TLabelType, typename TCoordType, unsigned int VDimension>
typename CenterlineGraphConverter<TLabelType, TCoordType, VDimension>::VesselGraphType 
CenterlineGraphConverter<TLabelType, TCoordType, VDimension>:: ImageToGraph(typename LabelImageType::Pointer labelMapImage) const
{
    typename itk::BinaryThinningImageFilter3D< LabelImageType, LabelImageType >::Pointer thinningFilter = itk::BinaryThinningImageFilter3D< LabelImageType, LabelImageType >::New();
    thinningFilter->SetInput(labelMapImage);
    thinningFilter->Update();
    typename LabelImageType::Pointer centerlineImage = thinningFilter->GetOutput();
   
    typename itk::SignedMaurerDistanceMapImageFilter< LabelImageType, CoordImageType >::Pointer distanceFilter = itk::SignedMaurerDistanceMapImageFilter< LabelImageType, CoordImageType >::New();
    distanceFilter->SetInput(labelMapImage);
    distanceFilter->SetSquaredDistance(false);
    distanceFilter->SetUseImageSpacing(true);
    distanceFilter->Update();

    typename itk::MultiplyByConstantImageFilter<CoordImageType, TCoordType, CoordImageType>::Pointer multFilter = itk::MultiplyByConstantImageFilter<CoordImageType, TCoordType, CoordImageType>::New();
    multFilter->SetConstant(-1.0);
    multFilter->SetInput(distanceFilter->GetOutput());
    multFilter->Update();
    
    typename itk::AddConstantToImageFilter<CoordImageType, TCoordType, CoordImageType>::Pointer addFilter = itk::AddConstantToImageFilter<CoordImageType, TCoordType, CoordImageType>::New();
    addFilter->SetConstant(1.0);
    addFilter->SetInput(multFilter->GetOutput());
    addFilter->Update();
    
    typename CoordImageType::Pointer radiusImage = addFilter->GetOutput();
 
/*  
    typename itk::DanielssonDistanceMapImageFilter< LabelImageType, CoordImageType >::Pointer distanceFilter = itk::DanielssonDistanceMapImageFilter< LabelImageType, CoordImageType >::New();
    distanceFilter->SetInput(labelMapImage);
    distanceFilter->SetSquaredDistance(false);
    distanceFilter->SetUseImageSpacing(true);
    distanceFilter->InputIsBinaryOn(); 
    distanceFilter->Update();
    typename CoordImageType::Pointer radiusImage = distanceFilter->GetOutput();
*/  
    typename LabelImageType::SpacingType spacing = labelMapImage->GetSpacing();
    typename LabelImageType::PointType origin = labelMapImage->GetOrigin();
    typename LabelImageType::RegionType region = labelMapImage->GetLargestPossibleRegion();
    

    typename CoordImageType::Pointer probabilityImage = CoordImageType::New();
    probabilityImage->SetSpacing( spacing );
    probabilityImage->SetOrigin( origin );
    probabilityImage->SetRegions( region );
    probabilityImage->Allocate(); 
    probabilityImage->FillBuffer(1.0);
    
    VectorType zeroVector;
    zeroVector.Fill(0.0);
    typename VectorImageType::Pointer tangentImage = VectorImageType::New();
    tangentImage->SetSpacing( spacing );
    tangentImage->SetOrigin( origin );
    tangentImage->SetRegions( region );
    tangentImage->Allocate(); 
    tangentImage->FillBuffer(zeroVector);
    
    VesselGraphType vesselGraph = this->ImageToGraph(centerlineImage, radiusImage,  probabilityImage, tangentImage);
    
    vesselGraph.SetTangentByPosition();
    vesselGraph.ComputeNodeParameters();
    vesselGraph.ComputeSegmentParameters();
    
    return(vesselGraph);
}


template<typename TLabelType, typename TCoordType, unsigned int VDimension>
typename CenterlineGraphConverter<TLabelType, TCoordType, VDimension>::CoordImageType::Pointer 
CenterlineGraphConverter<TLabelType, TCoordType, VDimension>::GraphToImage(typename LabelImageType::Pointer centerlineImage, int oversamplingFactor, const VesselGraphType& graph) const
{
    typename CoordImageType::Pointer upSampledLabelImage = CoordImageType::New();
    upSampledLabelImage->CopyInformation(centerlineImage);
    typename CoordImageType::RegionType upSampledRegion = centerlineImage->GetLargestPossibleRegion();
    typename CoordImageType::SpacingType upSampledSpacing = centerlineImage->GetSpacing();
    
    for(unsigned int i = 0; i < VDimension; i++)
    {
        upSampledRegion.SetSize(i, oversamplingFactor*upSampledRegion.GetSize(i));
        upSampledSpacing[i] = upSampledSpacing[i]/oversamplingFactor;
    }
    
    upSampledLabelImage->SetSpacing(upSampledSpacing);
    upSampledLabelImage->SetRegions( upSampledRegion );
    upSampledLabelImage->SetRequestedRegionToLargestPossibleRegion();
    upSampledLabelImage->Allocate(); 
    upSampledLabelImage->FillBuffer(0);
    typename LabelImageType::SpacingType imageSpacing = upSampledLabelImage->GetSpacing();
    
    typename VesselGraphType::ConstIterator segmentIter = graph.Begin();
    for(; segmentIter != graph.End(); ++segmentIter)
    {
        typename VesselGraphType::SegmentType currentSegment = segmentIter->second;
        typename VesselGraphType::SegmentType::ConstIterator nodeIter = currentSegment.Begin();
        for(; nodeIter != currentSegment.End(); ++nodeIter)
        {
            typename VesselGraphType::NodeIDType currentNodeID = *nodeIter;
            typename VesselGraphType::NodeType currentNode = graph.GetNode(currentNodeID);
            typename VesselGraphType::NodeType::PositionType currentPosition = currentNode.GetPosition();
            typename VesselGraphType::NodeType::RealType currentRadius = currentNode.GetRadius();
            
            int radiusFactor = 3;

            int xRadiusOffset = radiusFactor * (int(ceil(currentRadius/imageSpacing[0])) + 1);
            int yRadiusOffset = radiusFactor * (int(ceil(currentRadius/imageSpacing[1])) + 1);
            int zRadiusOffset = radiusFactor * (int(ceil(currentRadius/imageSpacing[2])) + 1);
            
            typename itk::SphereSpatialFunction<3, typename LabelImageType::PointType>::Pointer currentSphere = itk::SphereSpatialFunction<3, typename LabelImageType::PointType>::New();
            
            typename LabelImageType::PointType currentPoint;
            currentPoint[0] = currentPosition[0];
            currentPoint[1] = currentPosition[1];
            currentPoint[2] = currentPosition[2];
            currentSphere->SetCenter(currentPoint);
            currentSphere->SetRadius(currentRadius);

            for(int i = -xRadiusOffset; i <= xRadiusOffset; i++)
            {
                for(int j = -yRadiusOffset; j <= yRadiusOffset; j++)
                {
                    for(int k = -zRadiusOffset; k <= zRadiusOffset; k++)
                    {
                        typename LabelImageType::PointType movingPoint;
                        movingPoint[0] = currentPoint[0] + (TCoordType(i)/TCoordType(radiusFactor))*imageSpacing[0];
                        movingPoint[1] = currentPoint[1] + (TCoordType(j)/TCoordType(radiusFactor))*imageSpacing[1];
                        movingPoint[2] = currentPoint[2] + (TCoordType(k)/TCoordType(radiusFactor))*imageSpacing[2];
                    
                        bool insideSphere = currentSphere->Evaluate(movingPoint);
                        typename LabelImageType::IndexType movingIndex;
                        bool insideImage = upSampledLabelImage->TransformPhysicalPointToIndex(movingPoint, movingIndex);
                        
                        if(insideSphere == true && insideImage == true)
                        {
                            upSampledLabelImage->SetPixel(movingIndex, /*upSampledLabelImage->GetPixel(movingIndex) +*/ 1);
                        }
                    }
                }
            }
            
        }
            
    }
    

    typename itk::ResampleImageFilter< CoordImageType, CoordImageType, TCoordType>::Pointer resampler = itk::ResampleImageFilter< CoordImageType, CoordImageType, TCoordType>::New();
    resampler->SetInput(upSampledLabelImage);
    resampler->SetOutputParametersFromImage(centerlineImage);
    resampler->Update();
    typename CoordImageType::Pointer smoothedLabelMap = resampler->GetOutput();
    
    
    typename CoordImageType::Pointer downSampledLabelImage = CoordImageType::New();
    downSampledLabelImage->CopyInformation(centerlineImage);
    downSampledLabelImage->SetRegions(centerlineImage->GetLargestPossibleRegion() );
    downSampledLabelImage->SetRequestedRegionToLargestPossibleRegion();
    downSampledLabelImage->Allocate(); 
    downSampledLabelImage->FillBuffer(0);
    
    itk::ImageRegionIterator<CoordImageType> downIter( downSampledLabelImage, downSampledLabelImage->GetLargestPossibleRegion() );   
    typename itk::ConstNeighborhoodIterator< CoordImageType >::RadiusType neighborhoodRadius;
    unsigned int numNeighbors = (2*oversamplingFactor+1)*
        (2*oversamplingFactor+1)*(2*oversamplingFactor+1);
    neighborhoodRadius.Fill(oversamplingFactor);
    typename itk::ConstNeighborhoodIterator< CoordImageType > upNeighIter(neighborhoodRadius, upSampledLabelImage, 
        upSampledLabelImage->GetRequestedRegion() );
    
    for(downIter.GoToBegin(); !downIter.IsAtEnd(); ++downIter)
    {
        typename CoordImageType::IndexType currentIndex = downIter.GetIndex();
        typename CoordImageType::PointType currentPoint;
        downSampledLabelImage->TransformIndexToPhysicalPoint(currentIndex, currentPoint);
        typename CoordImageType::IndexType currentUpIndex;
        bool isInImage = upSampledLabelImage->TransformPhysicalPointToIndex(currentPoint, currentUpIndex);
        upNeighIter.SetLocation(currentUpIndex);
        typename itk::Neighborhood<TCoordType, VDimension> upHood = upNeighIter.GetNeighborhood();
        TCoordType currentValue = 0.0;    
        for (unsigned int i = 0; i < numNeighbors; i++)
        {
            TCoordType currentNeigh = upHood.GetElement(i);
            currentValue += currentNeigh;
        }
        downIter.Set(currentValue/numNeighbors);
    }
    
    //return(downSampledLabelImage);
    return(upSampledLabelImage);
}

