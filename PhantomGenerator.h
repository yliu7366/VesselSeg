#ifndef __PhantomGenerator_h
#define __PhantomGenerator_h

#include <iostream>
#include <vector>

#include "VesselGraph.h"
#include "itkImage.h"

/**\class PhantomGenerator
 * \brief Generate vessel tree phantom 
 *
 * Use L-system approach to construct simple vessel trees.
 */ 
template<typename TCoordType, unsigned int VDimension>
class PhantomGenerator
{
    public:      
    
    typedef VesselGraph< VDimension, TCoordType > GraphType;
    typedef typename GraphType::NodeType NodeType;
    typedef typename NodeType::PositionType PointType;
    typedef typename NodeType::VectorType VectorType;
    
    typedef itk::Image<TCoordType, VDimension> CoordImageType;
    
    /**
    * Generate graph using L-system approach and specified parameters.
    * This method starts from a root location grows a segment and then
    * bifurcates.  This process is repeated
    * @return graph
    */
    GraphType GenerateGraph() const;
    
    /**
    * Randomly perturb the node locations of a specified graph
    * @param g input graph
    * @param randomness the amount perturbation
    * @return perturbed graph
    */
    GraphType PerturbGraph(const GraphType& g, TCoordType randomness) const;
    
    /**
    * Interpolate the node locations of a specified graph
    * @param g input graph
    * @param oversamplingFactor factor which determines the number of additional
    *   nodes on each segment
    * @return interpolated graph
    */
    GraphType InterpolateGraph(const GraphType& g, TCoordType oversamplingFactor) const;
    
    /**
    * When constructing the segments in this method generates the next node
    * @param currentNode the node to grow the next node from
    * @param previousNode the node before currentNode
    * @param newID the id of the node to generated 
    * @param radialScaling we multiple this factor by the currentNodes radius
    *   to generate the radius used for the generated node
    * @param lengthScaling we multiple this factor by the distance between 
    *   previous nodes to determine the spacing
    * @param angle1 angle defining about of rotate 
    * @param angle2 angle defining about of rotate 
    * @return a node with all attributes populated except coordinate system
    *  and update the coordinate system on the currentNode
    */
    NodeType GenerateNextNode(NodeType& currentNode, NodeType& previousNode, typename NodeType::IDType newID, TCoordType radialScaling, TCoordType lengthScaling, TCoordType angle1, TCoordType angle2) const;
    
    /**
    * Generate the node after the root node
    * @param currentNode root node
    * @param newID the id of the node to generated 
    * @param radialScaling we multiple this factor by the currentNodes radius
    *   to generate the radius used for the generated node
    * @param lengthScaling we multiple this factor by the distance between  
    * @return first node with all attributes populated except coordinate system
    *   and update the coordinate system on the root node
    */
    NodeType GenerateFirstNode(NodeType& currentNode, typename NodeType::IDType newID, TCoordType radialScaling, TCoordType lengthScaling) const;
    
    /**
    * Generate a node which will begin a bifurcation
    * @param currentNode the node to grow the next node from
    * @param previousNode the node before currentNode
    * @param newID the id of the node to generated 
    * @param radialScaling we multiple this factor by the currentNodes radius
    *   to generate the radius used for the generated node
    * @param lengthScaling we multiple this factor by the distance between 
    *   previous nodes to determine the spacing
    * @param angle1 angle defining about of rotate 
    * @param angle2 angle defining about of rotate 
    * @return a node with all attributes populated except coordinate system
    *  and update the coordinate system on the currentNode
    */
    NodeType GenerateBifurcationNode(NodeType& currentNode, NodeType& previousNode, typename NodeType::IDType  newID, TCoordType lengthScaling, TCoordType radialScaling, TCoordType angle1, TCoordType angle2) const;
    
    /**
    * Generate a node just after a bifurcation
    * @param currentNode the node to grow the next node from
    * @param previousNode the node before currentNode
    * @param newID the id of the node to generated 
    * @param radialScaling we multiple this factor by the currentNodes radius
    *   to generate the radius used for the generated node
    * @param lengthScaling we multiple this factor by the distance between 
    *   previous nodes to determine the spacing
    * @param angle1 angle defining about of rotate 
    * @param angle2 angle defining about of rotate 
    * @return a node with all attributes populated except coordinate system
    *  and update the coordinate system on the currentNode
    */
    NodeType GenerateAfterBifurcationNode(NodeType& currentNode, NodeType& previousNode, typename NodeType::IDType  newID, TCoordType lengthScaling, TCoordType radialScaling, TCoordType angle1, TCoordType angle2) const;
    
    /**
    * Set the parameters 
    * @param leftRadialScaling radial scaling factor for left branches
    * @param rightRadialScaling radial scaling factor for right branches
    * @param straightRadialScaling radial scaling factor along a segment
    * @param leftLengthScaling spatial scaling for left branches
    * @param rightLengthScaling spatial scaling for right branches
    * @param straightLengthScaling spatial scaling factor along a segment
    * @param angleLeftIn angle 1 for left branches
    * @param angleLeftOut angle 2 for left branches
    * @param angleRightIn angle 1 for right branches
    * @param angleRightOut angle 2 for right branches
    * @param angleStraightIn angle 1 for segment
    * @param angleStraightOut angle 2 for segment
    * @param numScales number of generations in tree
    * @param maxSegmentLength largest number of nodes in a segment
    * @param minSegmentLength smallest number of nodes in a segment
    * @param root root node
    * @param randomnessFactor degree of randomness when generating tree
    */
    void SetParameters(TCoordType leftRadialScaling, TCoordType rightRadialScaling, TCoordType straightRadialScaling, TCoordType leftLengthScaling, TCoordType rightLengthScaling, TCoordType straightLengthScaling, TCoordType angleLeftIn, TCoordType angleLeftOut, TCoordType angleRightIn, TCoordType angleRightOut, TCoordType angleStraightIn, TCoordType angleStraightOut, unsigned int numScales, unsigned int maxSegmentLength, unsigned int minSegmentLength, NodeType root, TCoordType randomnessFactor);
    
    /**
    * Convert the tree representation to an image.  Map the node positions
    * to image coordinates using the specified spacing
    * @param g tree
    * @param border amount of spacing around border of image
    * @param spacing image spacing
    * @return an image representation of the tree
    */
    typename CoordImageType::Pointer GraphToImage(const GraphType& g, TCoordType border, const VectorType& spacing) const;

    /**
    * Rotate a vector in plane and out of plane by the specified amoumts.  
    * Generate local coordinate system using v1 as tangent.  We compute
    * normal and binormal vector.  
    * @param v1 vector which defines the tangent
    * @param v2 vector which is not colinear with v1.  It will be orthogonalized
    *   to a become the normal vector
    * @param a1 angle to rotate in tangent/normal plane
    * @param a2 angle to rotate in normal/binormal plane
    * @return rotated vector
    */
    VectorType Rotate(VectorType v1, VectorType v2, TCoordType a1, TCoordType a2) const;
       
    private:
    
    /**
    * Get random number
    */
    TCoordType GetRandom() const;
    
    /**
    * Get random number using specified see
    */
    TCoordType GetRandom(TCoordType) const;
    
    /**
    * Generate normal and binormal vectors
    */
    void GetNormalBiNormal(VectorType, VectorType, VectorType&, VectorType&) const;
    
    TCoordType leftRadialScaling;
    TCoordType rightRadialScaling;
    TCoordType straightRadialScaling;
    TCoordType leftLengthScaling;
    TCoordType rightLengthScaling;
    TCoordType straightLengthScaling;
    TCoordType angleLeftIn;
    TCoordType angleLeftOut;
    TCoordType angleRightIn;
    TCoordType angleRightOut;
    TCoordType angleStraightIn;
    TCoordType angleStraightOut;
    unsigned int numScales;
    unsigned int maxSegmentLength;
    unsigned int minSegmentLength;
    TCoordType randomnessFactor;
    NodeType root;
    

};



#include "PhantomGenerator.txx"

#endif
