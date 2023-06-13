#include <cmath>
#include <vector>
#include <utility>

#include "PhantomGenerator.h"
#include "CenterlineGraphConverter.h"

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator<TCoordType, VDimension>::NodeType 
PhantomGenerator< TCoordType, VDimension>::GenerateNextNode(NodeType& currentNode, NodeType& previousNode, typename NodeType::IDType newID,  TCoordType radialScaling, TCoordType lengthScaling, TCoordType angle1, TCoordType angle2) const
{
    NodeType newNode;
    newNode.SetID(newID);
    
    PointType currentPoint = currentNode.GetPosition();
    
    PointType previousPoint = previousNode.GetPosition();
    VectorType previousTangent = previousNode.GetTangent();
    VectorType previousNormal = previousNode.GetNormal();
    
    VectorType currentTangent = this->Rotate(previousTangent, previousNormal, angle1, angle2);
    currentTangent = currentTangent.Normalize();
    
    PointType newPosition = currentPoint + lengthScaling*currentTangent;
    
    newNode.SetPosition(newPosition);
    newNode.SetRadius(radialScaling * currentNode.GetRadius());
    newNode.SetMedialness(1);
    newNode.AddNeighbor(currentNode.GetID());

    currentNode.SetTangent(currentTangent);
    currentNode.AddNeighbor(newNode.GetID());
    
    VectorType normal;
    VectorType binormal;
    TCoordType dot = dot_product(currentTangent, previousTangent);
    if(fabs(dot) >= 0.999)
    {
        this->GetNormalBiNormal(currentTangent, previousNormal, normal, binormal);
    }
    else
    {
        this->GetNormalBiNormal(currentTangent, previousTangent, normal, binormal);
    }

    currentNode.SetNormal(normal);
    currentNode.SetBiNormal(binormal);
  
    return(newNode);
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator<TCoordType, VDimension>::NodeType 
PhantomGenerator< TCoordType, VDimension>::GenerateAfterBifurcationNode(NodeType& currentNode, NodeType& previousNode, typename NodeType::IDType  newID, TCoordType lengthScaling, TCoordType radialScaling, TCoordType angle1, TCoordType angle2) const
{
    NodeType newNode;
    newNode.SetID(newID);
    
    PointType currentPoint = currentNode.GetPosition();
    PointType previousPoint = previousNode.GetPosition(); 

    VectorType previousTangent = previousNode.GetTangent();
    VectorType currentPreviousVector = (currentPoint - previousPoint).Normalize();
    VectorType previousNormal;
    VectorType previousBinormal;
    this->GetNormalBiNormal(currentPreviousVector, previousTangent, previousNormal, previousBinormal);
    
    VectorType currentTangent = this->Rotate(currentPreviousVector, previousTangent, angle1, angle2);
    currentTangent = currentTangent.Normalize();
    
    PointType newPosition = currentPoint + lengthScaling*currentTangent;
    
    newNode.SetPosition(newPosition);
    newNode.SetRadius(radialScaling * currentNode.GetRadius());
    newNode.SetMedialness(1);
    newNode.AddNeighbor(currentNode.GetID());

    currentNode.SetTangent(currentTangent);
    currentNode.AddNeighbor(newNode.GetID());
    
    VectorType normal;
    VectorType binormal;
    VectorType v12 = currentNode.GetTangent();
    VectorType v22 = (currentPoint - previousPoint).Normalize();
    
    TCoordType dot = dot_product(currentTangent, currentPreviousVector);
    
    this->GetNormalBiNormal(currentTangent, previousTangent, normal, binormal); 
    
    currentNode.SetNormal(normal);
    currentNode.SetBiNormal(binormal);
    
    return(newNode);
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator<TCoordType, VDimension>::NodeType 
PhantomGenerator< TCoordType, VDimension>::GenerateFirstNode(NodeType& currentNode, typename NodeType::IDType newID, TCoordType radialScaling, TCoordType lengthScaling) const
{
    NodeType newNode;
    newNode.SetID(newID);
    
    PointType currentPoint = currentNode.GetPosition();
    VectorType currentTangent = currentNode.GetTangent();
    currentTangent = lengthScaling * currentTangent;
    
    PointType newPosition = currentPoint + currentTangent;

    newNode.SetPosition(newPosition);
    newNode.SetRadius(radialScaling * currentNode.GetRadius());
    newNode.SetMedialness(1);
    newNode.AddNeighbor(currentNode.GetID());
    
    currentNode.AddNeighbor(newNode.GetID());
    
    return(newNode);
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator<TCoordType, VDimension>::NodeType
PhantomGenerator< TCoordType, VDimension>::GenerateBifurcationNode(NodeType& currentNode, NodeType& previousNode, typename NodeType::IDType newID, TCoordType radialScaling,  TCoordType lengthScaling, TCoordType angle1, TCoordType angle2) const
{
    NodeType newNode;
    newNode.SetID(newID);

    PointType currentPoint = currentNode.GetPosition();
    PointType previousPoint = previousNode.GetPosition();

    VectorType previousTangent = previousNode.GetTangent();
    VectorType previousNormal = previousNode.GetNormal();
    
    VectorType currentTangent = this->Rotate(previousTangent, previousNormal, angle1, angle2);
    currentTangent = currentTangent.Normalize();

    PointType newPosition = currentPoint + lengthScaling * currentTangent;

    newNode.SetPosition(newPosition);
    newNode.SetRadius(radialScaling * currentNode.GetRadius());
    newNode.SetMedialness(1);
    newNode.AddNeighbor(currentNode.GetID());
   
    currentNode.AddNeighbor(newNode.GetID());
    currentTangent.Normalize();
    currentNode.SetTangent((currentPoint-previousPoint).Normalize());
    
    return(newNode);
}

template<typename TCoordType, unsigned int VDimension>
void PhantomGenerator< TCoordType, VDimension>::SetParameters(TCoordType leftRadialScaling, TCoordType rightRadialScaling, TCoordType straightRadialScaling, TCoordType leftLengthScaling, TCoordType rightLengthScaling, TCoordType straightLengthScaling, TCoordType angleLeftIn, TCoordType angleLeftOut, TCoordType angleRightIn, TCoordType angleRightOut, TCoordType angleStraightIn, TCoordType angleStraightOut, unsigned int numScales, unsigned int maxSegmentLength, unsigned int minSegmentLength, NodeType root , TCoordType randomnessFactor )
{
    this->leftRadialScaling = leftRadialScaling;
    this->rightRadialScaling = rightRadialScaling;
    this->straightRadialScaling = straightRadialScaling;
    this->leftLengthScaling = leftLengthScaling;
    this->rightLengthScaling = rightLengthScaling;
    this->straightLengthScaling = straightLengthScaling;
    this->angleLeftIn = angleLeftIn;
    this->angleLeftOut = angleLeftOut;
    this->angleRightIn = angleRightIn;
    this->angleRightOut = angleRightOut;
    this->angleStraightIn = angleStraightIn;
    this->angleStraightOut = angleStraightOut;
    this->numScales = numScales;
    this->maxSegmentLength = maxSegmentLength;
    this->minSegmentLength = minSegmentLength;
    this->root = root;
    this->randomnessFactor = randomnessFactor;
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator< TCoordType, VDimension>::GraphType PhantomGenerator< TCoordType, VDimension>::GenerateGraph() const
{

    typedef std::pair<NodeType, NodeType> NodePairType;
    typedef std::vector<NodePairType> NodePairListType;

    GraphType graph;
    typename NodeType::IDType nodeIDCounter = this->root.GetID();
    typename GraphType::SegmentType::IDType segmentIDCounter = 1;
    
    typename GraphType::NodeType::VectorType zeroVector;
    zeroVector.Fill(0.0);
    
    typename GraphType::SegmentType firstSegment;
    firstSegment.SetID(segmentIDCounter);
    ++segmentIDCounter;

    typename GraphType::NodeType firstNode = this->root;
    ++nodeIDCounter;
    firstNode.AddSegment(firstSegment.GetID());
    firstSegment.AddNode(firstNode.GetID());   
    
    //Get node after root
    typename GraphType::NodeType secondNode = GenerateFirstNode(firstNode, nodeIDCounter, this->straightRadialScaling, this->straightLengthScaling);
    ++nodeIDCounter;
    secondNode.AddSegment(firstSegment.GetID());
    firstSegment.AddNode(secondNode.GetID());
    
    graph.AddNode(firstNode);
    
    //Build first segment
    typename GraphType::NodeType previousNode = firstNode;
    typename GraphType::NodeType currentNode = secondNode;
    
    unsigned int numNodesFirstSegment = this->minSegmentLength + rand() % (this->maxSegmentLength-this->minSegmentLength) + 1;
    for(unsigned int i = 0; i < numNodesFirstSegment; i++)
    {
        TCoordType angle1 = (1.0/previousNode.GetRadius())*  this->GetRandom() * this->angleStraightIn;
        TCoordType angle2 = (1.0/previousNode.GetRadius())*  this->GetRandom() * this->angleStraightOut;
        NodeType newNode = this->GenerateNextNode(currentNode, previousNode, nodeIDCounter, this->straightRadialScaling, this->straightLengthScaling, angle1, angle2);
        ++nodeIDCounter;
        newNode.AddSegment(firstSegment.GetID());
        firstSegment.AddNode(newNode.GetID());
        
        graph.AddNode(currentNode);
        
        previousNode = currentNode;
        currentNode = newNode;      
    }
    
    
    currentNode.SetTangent(zeroVector);
    graph.AddNode(currentNode);
    graph.AddSegment(firstSegment);
    
    
    //Build branches of tree
    NodePairType firstPair(currentNode, previousNode); 
    NodePairListType currentJunctions;
    currentJunctions.push_back(firstPair);
    
    for(unsigned int j = 1; j <= this->numScales; j++)
    {
        int numJunctions = pow(2, j-1);
        NodePairListType newJunctions;
        for(unsigned int k = 1; k <= numJunctions; k++)
        {
            NodePairType currentPair = currentJunctions[k-1];
            
            NodeType currentJunctionNode = currentPair.first;
            NodeType previousJunctionNode = currentPair.second;
            
            //left branch
            TCoordType angle1L = this->GetRandom() + this->angleLeftIn;
            TCoordType angle2L = this->GetRandom() + this->angleLeftOut;

            NodeType newLeftNode = this->GenerateBifurcationNode(currentJunctionNode, previousJunctionNode, nodeIDCounter, this->leftRadialScaling, this->leftLengthScaling, angle1L, angle2L);
            ++nodeIDCounter;
                    
            typename GraphType::SegmentType leftSegment;
            leftSegment.SetID(segmentIDCounter);
            ++segmentIDCounter;
            
            currentJunctionNode.AddNeighbor(newLeftNode.GetID());
            currentJunctionNode.AddSegment(leftSegment.GetID());
            leftSegment.AddNode(currentJunctionNode.GetID()); 
            newLeftNode.AddSegment(leftSegment.GetID());
            leftSegment.AddNode(newLeftNode.GetID());
            
            TCoordType angle1 = (1.0/newLeftNode.GetRadius())*  this->GetRandom() * this->angleStraightIn;
            TCoordType angle2 = (1.0/newLeftNode.GetRadius())*  this->GetRandom() * this->angleStraightOut;
            
            NodeType newLeftNextNode = this->GenerateAfterBifurcationNode(newLeftNode, currentJunctionNode, nodeIDCounter, this->leftRadialScaling, this->leftLengthScaling, angle1, angle2);
            ++nodeIDCounter;
            
            newLeftNode.AddNeighbor(newLeftNextNode.GetID());
            newLeftNextNode.AddNeighbor(newLeftNode.GetID());
            newLeftNextNode.AddSegment(leftSegment.GetID());
            leftSegment.AddNode(newLeftNextNode.GetID()); 
            
            graph.AddNode(newLeftNode);

            typename GraphType::NodeType previousLeftNode = newLeftNode;
            typename GraphType::NodeType currentLeftNode = newLeftNextNode;
            
            unsigned int numNodesLeftSegment = this->minSegmentLength + rand() % (this->maxSegmentLength-this->minSegmentLength) + 1;
            for(unsigned int i = 0; i < numNodesLeftSegment; i++)
            {
                angle1 = (1.0/previousLeftNode.GetRadius())* this->GetRandom() * this->angleStraightIn;
                angle2 = (1.0/previousLeftNode.GetRadius())* this->GetRandom() * this->angleStraightOut;
            
                NodeType newLeftStraightNode = this->GenerateNextNode(currentLeftNode, previousLeftNode, nodeIDCounter, this->straightRadialScaling, this->straightLengthScaling, angle1, angle2);
                ++nodeIDCounter;
                newLeftStraightNode.AddSegment(leftSegment.GetID());
                leftSegment.AddNode(newLeftStraightNode.GetID());
                
                graph.AddNode(currentLeftNode);
                
                previousLeftNode = currentLeftNode;
                currentLeftNode = newLeftStraightNode;      
            }
            
            currentLeftNode.SetTangent(zeroVector);
            graph.AddNode(currentLeftNode);
            graph.AddSegment(leftSegment);
            
            NodePairType leftPair(currentLeftNode, previousLeftNode); 
            newJunctions.push_back(leftPair);
            
            //right
            TCoordType angle1R = this->GetRandom() + this->angleRightIn;
            TCoordType angle2R = this->GetRandom() + this->angleRightOut;
            
            NodeType newRightNode = this->GenerateBifurcationNode(currentJunctionNode, previousJunctionNode, nodeIDCounter, this->rightRadialScaling, this->rightLengthScaling, angle1R, angle2R);
            ++nodeIDCounter;
 
            typename GraphType::SegmentType rightSegment;
            rightSegment.SetID(segmentIDCounter);
            segmentIDCounter++;

            currentJunctionNode.AddNeighbor(newRightNode.GetID());
            currentJunctionNode.AddSegment(rightSegment.GetID());
            rightSegment.AddNode(currentJunctionNode.GetID()); 
            newRightNode.AddSegment(rightSegment.GetID());
            rightSegment.AddNode(newRightNode.GetID());
            
            angle1 = -(1.0/newRightNode.GetRadius())*  this->GetRandom() * this->angleStraightIn;
            angle2 = -(1.0/newRightNode.GetRadius())*  this->GetRandom() * this->angleStraightOut;
                
            NodeType newRightNextNode = this->GenerateAfterBifurcationNode(newRightNode, currentJunctionNode, nodeIDCounter, this->leftRadialScaling, this->leftLengthScaling, angle1, angle2);
            ++nodeIDCounter;
            
            newRightNode.AddNeighbor(newRightNextNode.GetID());
            newRightNextNode.AddNeighbor(newRightNode.GetID());
            newRightNextNode.AddSegment(rightSegment.GetID());
            rightSegment.AddNode(newRightNextNode.GetID()); 
            
            graph.AddNode(newRightNode);
            
            typename GraphType::NodeType previousRightNode = newRightNode;
            typename GraphType::NodeType currentRightNode = newRightNextNode;
            
            unsigned int numNodesRightSegment = this->minSegmentLength + rand() % (this->maxSegmentLength-this->minSegmentLength) + 1;
            for(unsigned int i = 0; i < numNodesRightSegment; i++)
            {
                angle1 = -(1.0/previousRightNode.GetRadius())*  this->GetRandom() * this->angleStraightIn;
                angle2 = -(1.0/previousRightNode.GetRadius())*  this->GetRandom() * this->angleStraightOut;
                
                NodeType newRightStraightNode = this->GenerateNextNode(currentRightNode, previousRightNode, nodeIDCounter, this->straightRadialScaling, this->straightLengthScaling, angle1, angle2);
                ++nodeIDCounter;
                newRightStraightNode.AddSegment(rightSegment.GetID());
                rightSegment.AddNode(newRightStraightNode.GetID());
                
                graph.AddNode(currentRightNode);
                
                previousRightNode = currentRightNode;
                currentRightNode = newRightStraightNode;      
            }
            
            graph.AddNode(currentRightNode);
            graph.AddSegment(rightSegment);
            
            graph.AddNode(currentJunctionNode);
            
            NodePairType rightPair(currentRightNode, previousRightNode); 
            newJunctions.push_back(rightPair);
            
        }
        currentJunctions = newJunctions;
    }
    

    graph.ComputeFrenetParameters();
    graph.ComputeSegmentParameters();

    return(graph);
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator< TCoordType, VDimension>::GraphType
PhantomGenerator< TCoordType, VDimension>::PerturbGraph(const GraphType& g, TCoordType randomness) const
{
    GraphType perturbGraph = g;
    
    typename GraphType::ConstIterator segmentIter = perturbGraph.Begin();
    for(; segmentIter != perturbGraph.End(); ++segmentIter)
    {
        typename GraphType::SegmentType currentSegment = segmentIter->second;
        typename GraphType::NodeIDType currentFrontNodeID = currentSegment.Front();
        typename GraphType::NodeIDType currentBackNodeID = currentSegment.Back();
        typename GraphType::SegmentType::ConstIterator nodeIter = currentSegment.Begin();
        typename GraphType::SegmentType::ConstIterator nodeAheadIter = currentSegment.Begin();
        nodeAheadIter++;
        for(; nodeAheadIter != currentSegment.End(); ++nodeIter, ++nodeAheadIter)
        {
            typename GraphType::NodeIDType currentNodeID = *nodeIter;
            typename GraphType::NodeType currentNode = perturbGraph.GetNode(currentNodeID);
            typename GraphType::NodeType oldCurrentNode = g.GetNode(currentNodeID);
            typename GraphType::NodeType::PositionType currentPosition = currentNode.GetPosition();
            
            typename GraphType::NodeIDType nextNodeID = *nodeAheadIter;
            typename GraphType::NodeType nextNode = perturbGraph.GetNode(nextNodeID);
            typename GraphType::NodeType::PositionType nextPosition = nextNode.GetPosition();
            

            //change this junction point position
            TCoordType xPert = this->GetRandom(randomness);
            TCoordType yPert = this->GetRandom(randomness);
            TCoordType zPert = this->GetRandom(randomness);
            nextPosition[0] = nextPosition[0] + xPert;
            nextPosition[1] = nextPosition[1] + yPert;
            nextPosition[2] = nextPosition[2] + zPert;
            
            nextNode.SetPosition(nextPosition);
            perturbGraph.AddNode(nextNode);
            
            if(currentNodeID != currentFrontNodeID)
            {
                typename GraphType::NodeType::VectorType theTangent = nextPosition-currentPosition;
                theTangent.Normalize();
                currentNode.SetTangent(theTangent);
                perturbGraph.AddNode(currentNode);       
            }    
        }
    }
    
    perturbGraph.ComputeNodeParameters();
    perturbGraph.ComputeSegmentParameters();

    return(perturbGraph);
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator< TCoordType, VDimension>::GraphType
PhantomGenerator< TCoordType, VDimension>::InterpolateGraph(const GraphType& graph, TCoordType oversamplingFactor) const
{
    GraphType interpGraph = graph.Interpolate(oversamplingFactor);
    
    typename GraphType::ConstIterator segmentIter = interpGraph.Begin();
    for(; segmentIter != interpGraph.End(); ++segmentIter)
    {
        typename GraphType::SegmentType currentSegment = segmentIter->second;
        typename GraphType::SegmentType::ConstIterator nodeIter = currentSegment.Begin();
        typename GraphType::SegmentType::ConstIterator nodeAheadIter = currentSegment.Begin();
        nodeAheadIter++;
        
        typename GraphType::SegmentType::ConstIterator lastNode = currentSegment.End();
        --lastNode;
        for(; nodeAheadIter != lastNode; ++nodeIter, ++nodeAheadIter)
        {
            typename GraphType::NodeIDType previousNodeID = *nodeIter;
            typename GraphType::NodeType previousNode = interpGraph.GetNode(previousNodeID);
            typename GraphType::NodeType::PositionType previousPosition = previousNode.GetPosition();
            typename GraphType::NodeType::VectorType previousTangent = previousNode.GetTangent();
            typename GraphType::NodeType::VectorType previousNormal = previousNode.GetNormal();
            
            typename GraphType::NodeIDType currentNodeID = *nodeAheadIter;
            typename GraphType::NodeType currentNode = interpGraph.GetNode(currentNodeID);
            typename GraphType::NodeType::PositionType currentPosition = currentNode.GetPosition();
            
            typename GraphType::NodeType::VectorType currentTangent = (currentPosition - previousPosition).Normalize();
            
            VectorType normal;
            VectorType binormal;
            TCoordType dot = dot_product(currentTangent, previousTangent);
            if(fabs(dot) >= 0.999)
            {
                this->GetNormalBiNormal(currentTangent, previousNormal, normal, binormal);
            }
            else
            {
                this->GetNormalBiNormal(currentTangent, previousTangent, normal, binormal);
            }

            currentNode.SetNormal(normal);
            currentNode.SetBiNormal(binormal);
            
            interpGraph.AddNode(currentNode);
        }
    }
    
    interpGraph.ComputeNodeParameters();
    interpGraph.ComputeSegmentParameters();
        
    return(interpGraph);
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator< TCoordType, VDimension>::CoordImageType::Pointer 
PhantomGenerator< TCoordType, VDimension>::GraphToImage(const GraphType& graph, TCoordType border, const VectorType& spacing) const
{
    typename CoordImageType::Pointer inputImage = CoordImageType::New();
    typedef CenterlineGraphConverter<TCoordType, TCoordType, VDimension> ConverterType;
    ConverterType converter;
    
    PointType corner1;
    PointType corner2; 
    graph.GetBoundingPoints(corner1, corner2);
    
    for(unsigned int i = 0; i < VDimension; i++)
    {
        corner1[i] = corner1[i] - border;
        corner2[i] = corner2[i] + border;
    }                
    
    typename CoordImageType::PointType imageOrigin;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        imageOrigin[i] = corner1[i];
    }
    inputImage->SetOrigin(imageOrigin);
    
    typename CoordImageType::SpacingType imageSpacing;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        imageSpacing[i] = spacing[i];
    }
    inputImage->SetSpacing(imageSpacing);
    
    typename CoordImageType::SizeType imageSize;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        imageSize[i] = typename CoordImageType::SizeType::SizeValueType( (ceil(corner2[i]- corner1[i])/imageSpacing[i]));
    }
    
    typename CoordImageType::IndexType imageStart;
    imageStart.Fill(0);
    typename CoordImageType::RegionType imageRegion;
    imageRegion.SetSize( imageSize );
    imageRegion.SetIndex( imageStart );
    inputImage->SetRegions( imageRegion );
    inputImage->Allocate();
    inputImage->FillBuffer( 0 );
    
    
    
    typename CoordImageType::Pointer outputImage = converter.GraphToImage(inputImage, 1, graph);
    return(outputImage);
    
}

template<typename TCoordType, unsigned int VDimension>
typename PhantomGenerator< TCoordType, VDimension>::VectorType 
PhantomGenerator< TCoordType, VDimension>::Rotate(VectorType v1, VectorType v2, TCoordType angle1, TCoordType angle2) const
{
    VectorType tangent = v1.Normalize();
    VectorType normal;
    VectorType binormal;
    
    this->GetNormalBiNormal(v1, v2, normal, binormal);

    VectorType inPlaneVector = cos(angle1)*tangent + sin(angle1)*normal;
    VectorType rotateVector = cos(angle2)*inPlaneVector + sin(angle2)*binormal;
    rotateVector.Normalize();
    
    return(rotateVector);
}

template<typename TCoordType, unsigned int VDimension>
TCoordType PhantomGenerator< TCoordType, VDimension>::GetRandom() const
{
    int num = int(ceil(1000.0*double(this->randomnessFactor)));
    if(num < 1000)
    {
        num = 1000;
    }
    TCoordType randNum = this->randomnessFactor*(double(rand() % (2*num)) - double(num)) / double(num);
    return(randNum);
}

template<typename TCoordType, unsigned int VDimension>
TCoordType PhantomGenerator< TCoordType, VDimension>::GetRandom(TCoordType randomness) const
{
    int num = int(ceil(1000.0*double(randomness)));
    if(num < 1000)
    {
        num = 1000;
    }
    TCoordType randNum = randomness*((double(rand() % (2*num)) - double(num)) / double(num));
    
    //std::cout << randNum << " " << ((double(rand() % (2*num)) - double(num)) / double(num))<< ", ";
    return(randNum);
}

template<typename TCoordType, unsigned int VDimension>
void
PhantomGenerator< TCoordType, VDimension>::GetNormalBiNormal(VectorType v1, VectorType v2, VectorType& normal, VectorType& binormal) const
{
    VectorType tangent = v1.Normalize();
    v2 = v2.Normalize();
    
    TCoordType dot = dot_product(tangent, v2);
    
    VectorType temp;
    temp.Fill(0.0);
    
    if(fabs(dot) > 0.999)
    { 
        if(fabs(tangent[2]) <= fabs(tangent[1]) && fabs(tangent[2]) <= fabs(tangent[0]))
        {
            temp[2] = 1.0;
        }
        else if(fabs(tangent[1]) <= fabs(tangent[2]) && fabs(tangent[1]) <= fabs(tangent[0]))
        {
            temp[1] = 1.0;
        }
        else 
        {
            temp[0] = 1.0;
        }

        TCoordType proj = dot_product(temp, tangent);
        normal = temp - proj*tangent;
        normal = normal.Normalize();

        binormal = cross_product(tangent, normal);    
        binormal = binormal.Normalize();
    }
    else
    {
        TCoordType proj = dot_product(v2, tangent);
        normal = v2 - proj*tangent;
        
        normal = normal.Normalize();
        binormal = cross_product(tangent, normal); 
        binormal = binormal.Normalize();
    }
    
    
    TCoordType proj = dot_product(v2, tangent);
    normal = v2 - proj*tangent;
    
    normal = normal.Normalize();
    binormal = cross_product(tangent, normal); 
    binormal = binormal.Normalize();
    
}

