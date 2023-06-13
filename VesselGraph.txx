#include <cmath>
#include <sstream>
#include <list>
#include <vector>
#include <functional>
#include <iterator>

#include "VesselGraph.h"
#include "Interpolater.h"


template<unsigned int VIndexDimension, typename TCoordinateType>
VesselGraph<VIndexDimension, TCoordinateType>::VesselGraph()
{

}

template<unsigned int VIndexDimension, typename TCoordinateType>
VesselGraph<VIndexDimension, TCoordinateType>::~VesselGraph()
{

}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::Clear()
{
    this->nodeTable.clear();
    this->segmentTable.clear();
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::AddNode(NodeType node)
{
    this->nodeTable[node.GetID()] = node;
    //std::cout << "Added Node ID " << node.GetID() << std::endl;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::NodeType 
VesselGraph<VIndexDimension, TCoordinateType>::GetNode(const NodeIDType& id) const
{
    typename NodeTableType::const_iterator foundIter = this->nodeTable.find(id);
    if(foundIter != this->nodeTable.end())
    {
        return(foundIter->second);
    }
    else
    {
        std::stringstream message;
        message << "Node " << id << " not in Graph" << std::endl;
        //std::cout << message.str();
        throw(message.str());
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::AddSegment(SegmentType segment)
{
    this->segmentTable[segment.GetID()] = segment;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::SegmentType 
VesselGraph<VIndexDimension, TCoordinateType>::GetSegment(const SegmentIDType& id) const
{
    typename SegmentTableType::const_iterator foundIter = this->segmentTable.find(id);
    if(foundIter != this->segmentTable.end())
    {
        return(foundIter->second);
    }
    else
    {
        std::stringstream message;
        message << "Segment " << id << " not in Graph" << std::endl;
        //std::cout << message.str();
        throw(message.str());
    }
}
template<unsigned int VIndexDimension, typename TCoordinateType>
bool VesselGraph<VIndexDimension, TCoordinateType>::IsSegmentConnected(const SegmentIDType& id) const
{
    bool isConnected = true;
    
    typename SegmentTableType::const_iterator foundIter = this->segmentTable.find(id);
    try
    {
        if(foundIter != this->segmentTable.end())
        {
            SegmentType segment = foundIter->second;
            NodeIDType frontID = segment.Front();
            NodeType frontNode = this->GetNode(frontID);
            NodeIDType backID = segment.Back();
            NodeType backtNode = this->GetNode(backID);
            
            if( (frontNode.GetType() == NodeType::TERMINAL && backtNode.GetType() == NodeType::TERMINAL) || (frontNode.GetType() == NodeType::SINGLE && backtNode.GetType() == NodeType::SINGLE) )
            {
                isConnected = false;
            }
        }
    }
    catch(...)
    {
    
    }
    
    return(isConnected);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::RemoveUnconnectedSegments(unsigned int numNodes)
{
    typedef std::vector<NodeIDType> NodeIDListType;
    typedef std::vector<SegmentIDType> SegmentIDListType;
    
    NodeIDListType removeNodeIDList;
    SegmentIDListType removeSegmentIDList;
    
    try
    {
        ConstIterator segmentIter = this->Begin();
        for(; segmentIter != this->End(); ++segmentIter)
        {
            SegmentIDType currentSegmentID = segmentIter->first;
            SegmentType currentSegment = segmentIter->second;
            if(currentSegment.NumberOfNodes() <= numNodes && this->IsSegmentConnected(currentSegmentID) == false)
            {  
                typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
                for(; nodeIter != currentSegment.End(); ++nodeIter)
                {
                    NodeIDType currentNodeID = *nodeIter;

                    if(std::find( removeNodeIDList.begin(), removeNodeIDList.end(), currentNodeID) == removeNodeIDList.end() )
                    {
                        removeNodeIDList.push_back(currentNodeID);
                    }
                }
                
                if(std::find( removeSegmentIDList.begin(), removeSegmentIDList.end(), currentSegmentID) == removeSegmentIDList.end() )
                {
                    removeSegmentIDList.push_back(currentSegmentID);
                }
            }
        }
        
        typename NodeIDListType::const_iterator removeNodeIter = removeNodeIDList.begin();
        for(; removeNodeIter != removeNodeIDList.end(); ++removeNodeIter)
        {
            this->nodeTable.erase(*removeNodeIter);
        }
        typename SegmentIDListType::const_iterator removeSegmentIter = removeSegmentIDList.begin();
        for(; removeSegmentIter != removeSegmentIDList.end(); ++removeSegmentIter)
        {
            this->segmentTable.erase(*removeSegmentIter);
        }
        
    }
    catch(...)
    {
    
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Self 
VesselGraph<VIndexDimension, TCoordinateType>::Interpolate(RealType oversamplingFactor) const
{
    typename NodeType::VectorType zeroVector;
    zeroVector.Fill(0.0);
    NodeIDType newNodeIDCoutner = 0;

    Self graph;
    typedef std::vector< typename NodeType::PositionType > PositionListType;
    typedef Interpolater<TCoordinateType, VIndexDimension> InterpolaterType;
    InterpolaterType interpolater;
    typedef Interpolater<TCoordinateType, 1> OneInterpolaterType;
    OneInterpolaterType oneInterpolater;
    int order = 3;
    
    typename NodeTableType::const_iterator nodeOnlyIter = this->nodeTable.begin();
    for(; nodeOnlyIter != this->nodeTable.end(); nodeOnlyIter++)
    {
        NodeIDType currentID = nodeOnlyIter->first;
        if(currentID > newNodeIDCoutner)
        {
            newNodeIDCoutner = currentID;
        }
    }
    newNodeIDCoutner++;
    try
    {
        ConstIterator segmentIter = this->Begin();
        for(; segmentIter != this->End(); ++segmentIter)
        {
            SegmentType currentSegment = segmentIter->second;
            SegmentType currentInterpolatedSegment;
            currentInterpolatedSegment.SetID(currentSegment.GetID());
            
            if(currentSegment.NumberOfNodes() > 2*order)
            {
                order = 3;
            }
            else if(currentSegment.NumberOfNodes() <= 2*order && currentSegment.NumberOfNodes() > 2)
            {
                order = int(floor( currentSegment.NumberOfNodes() / 2));
            }

            if(currentSegment.NumberOfNodes() > 2*order)
            {
                //get old nodes they will become control points for interp
                PositionListType positionList;
                typename OneInterpolaterType::PointListType radiusList;
                typename OneInterpolaterType::PointListType medialnessList;
                typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
                for(; nodeIter != currentSegment.End(); ++nodeIter)
                {
                    NodeIDType currentNodeID = *nodeIter;
                    NodeType currentNode = this->GetNode(currentNodeID);
                    typename NodeType::PositionType currentPosition = currentNode.GetPosition();
                    positionList.push_back(currentPosition);
                    
                    typename NodeType::RealType currentRadius = currentNode.GetRadius();
                    typename OneInterpolaterType::PointType currentRadiusPoint;
                    currentRadiusPoint[0] = currentRadius;
                    radiusList.push_back(currentRadiusPoint);
                    
                    typename NodeType::RealType currentMedialness = currentNode.GetMedialness();
                    typename OneInterpolaterType::PointType currentMedialnessPoint;
                    currentMedialnessPoint[0] = currentMedialness;
                    medialnessList.push_back(currentMedialnessPoint);
                }
                
                //build knot vector
                int numControlPoints = positionList.size();
                int numKnots = numControlPoints + order + 1;
                typename InterpolaterType::TValueListType knotVector = interpolater.GetEqualKnots(order, numKnots);
                //typename InterpolaterType::TValueListType knotVector = interpolater.GetEndWeightKnots(order, numKnots);

                //keep first node the same
                NodeIDType currentFrontNodeID = currentSegment.Front();
                NodeType currentFrontNode;
                try
                {
                    currentFrontNode = graph.GetNode(currentFrontNodeID);
                }
                catch(...)
                {
                    currentFrontNode = this->GetNode(currentFrontNodeID);
                }
                
                RealType currentAverageRadius = currentSegment.GetAverageRadius();
                
                //remove old neighbors and add new ones
                typename SegmentType::ConstIterator currentFrontNextIter = currentSegment.Begin();
                std::advance(currentFrontNextIter,1);
                NodeIDType currentFrontNextNodeID = *(currentFrontNextIter);
                currentFrontNode.RemoveNeighbor(currentFrontNextNodeID);
                currentFrontNode.AddNeighbor(newNodeIDCoutner);
                //currentFrontNode.SetRadius(currentAverageRadius);
                //currentFrontNode.SetTangent(zeroVector);
                graph.AddNode(currentFrontNode);
                
                currentInterpolatedSegment.AddNode(currentFrontNode.GetID());
                
                //Do Interpolation for interior points only
                int numInterpPoints = int(ceil(oversamplingFactor*numControlPoints));
                for(int i = 1; i < numInterpPoints; i++)
                {
                    RealType t = RealType(i)/RealType(numInterpPoints);
                    typename InterpolaterType::PointType currentPoint = interpolater.EvaluateCurvePoint( numControlPoints-1, order, t, knotVector, positionList);
                    typename InterpolaterType::PointType currentDerivative = interpolater.EvaluateDerivatePoint( numControlPoints-1, order, t, knotVector, positionList);
                    typename OneInterpolaterType::PointType currentRadius = oneInterpolater.EvaluateCurvePoint(numControlPoints-1, order, t, knotVector, radiusList);
                    typename OneInterpolaterType::PointType currentMedialness = oneInterpolater.EvaluateCurvePoint(numControlPoints-1, order, t, knotVector, medialnessList);
                    
                    NodeType newNode;
                    newNode.SetID(newNodeIDCoutner);
                    newNode.SetPosition(currentPoint);
                    newNode.SetRadius(currentRadius[0]);
                    //newNode.SetRadius(currentAverageRadius);
                    newNode.SetMedialness(currentMedialness[0]);
                    newNode.SetTangent(zeroVector);
                    if(fabs(currentDerivative.Norm()) > 1.0e-6)
                    {
                        newNode.SetTangent(currentDerivative.Normalize());
                    }
                    newNode.AddNeighbor(newNodeIDCoutner-1);
                    newNode.AddNeighbor(newNodeIDCoutner+1);
                    newNode.AddSegment(currentInterpolatedSegment.GetID());
                    graph.AddNode(newNode);
                    currentInterpolatedSegment.AddNode(newNode.GetID());
                    newNodeIDCoutner++;
                }
                NodeIDType currentBackNodeID = currentSegment.Back();

                NodeType currentBackNode;
                try
                {
                    currentBackNode = graph.GetNode(currentBackNodeID);
                }
                catch(...)
                {
                    currentBackNode = this->GetNode(currentBackNodeID);
                }
                //remove old neighbors and add new one
                typename SegmentType::ConstIterator currentBackPreviousIter = currentSegment.Begin();
                std::advance(currentBackPreviousIter,currentSegment.NumberOfNodes()-2);
                NodeIDType currentBackPreviousNodeID = *(currentBackPreviousIter);
                currentBackNode.RemoveNeighbor(currentBackPreviousNodeID);
                currentBackNode.AddNeighbor(newNodeIDCoutner-1);
                //currentBackNode.SetRadius(currentAverageRadius);
                //currentBackNode.SetTangent(zeroVector);
                graph.AddNode(currentBackNode);
                currentInterpolatedSegment.AddNode(currentBackNode.GetID());
                
                graph.AddSegment(currentInterpolatedSegment);
               
            }
            else
            {
                typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
                for(; nodeIter != currentSegment.End(); ++nodeIter)
                {
                    NodeIDType currentNodeID = *nodeIter;
                    NodeType currentNode;
                    try
                    {
                        currentNode = graph.GetNode(currentNodeID);
                    }
                    catch(...)
                    {
                        currentNode = this->GetNode(currentNodeID);
                    }
                    graph.AddNode(currentNode);
                }
                graph.AddSegment(currentSegment);
            }
        }
        
    }
    catch(...)
    {
        
    }
 
    graph.ComputeNodeParameters();
    graph.ComputeSegmentParameters();
         
    return(graph);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Self 
VesselGraph<VIndexDimension, TCoordinateType>::GetDownSampledByCurvature(RealType curvatureThreshold) const
{
    Self graph;
    
    typedef std::list<NodeIDType> NodeIDListType;
    try
    {
        ConstIterator segmentIter = this->Begin();
        for(; segmentIter != this->End(); ++segmentIter)
        {
            SegmentType currentSegment = segmentIter->second;
            SegmentType currentDownSampledSegment;
            currentDownSampledSegment.SetID(currentSegment.GetID());
            
            
            if(currentSegment.NumberOfNodes() > 2)
            {
                NodeIDListType removeNodeList;
                
                //don't remove ends
                typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
                NodeIDType currentStartNodeID = *nodeIter;
                NodeType currentStartNode = this->GetNode(currentStartNodeID);
                currentDownSampledSegment.AddNode(currentStartNodeID);
                graph.AddNode(currentStartNode);
                ++nodeIter;
                
                NodeIDType currentStopNodeID = currentSegment.Back();

                for(; nodeIter != currentSegment.End(); ++nodeIter)
                {
                    NodeIDType currentNodeID = *nodeIter;
                    NodeType currentNode = this->GetNode(currentNodeID);
                    typename NodeType::RealType currentCurvature = currentNode.GetCurvature();
                    
                    if(fabs(currentCurvature) <= curvatureThreshold && (currentNodeID != currentStopNodeID) )
                    {
                        removeNodeList.push_back(currentNodeID);
                    }
                    else
                    {
                        //get ready to add nodes to graph but make sure they
                        //are not there already
                        
                        try
                        {
                            currentNode = graph.GetNode(currentNodeID);
                        }
                        catch(...)
                        {
                        }
                        
                        try
                        {
                            currentStartNode = graph.GetNode(currentStartNodeID);
                        }
                        catch(...)
                        {
                        }
                        
                        //remove neighbors which may be in remove list
                        typename NodeIDListType::const_iterator removeIter = removeNodeList.begin();
                        for(; removeIter != removeNodeList.end(); ++removeIter)
                        {
                            NodeIDType currentRemoveNodeID = *removeIter;
                            currentNode.RemoveNeighborIfFound(currentRemoveNodeID);
                            currentStartNode.RemoveNeighborIfFound(currentRemoveNodeID);
                        }
                        
                        currentStartNode.AddNeighbor(currentNodeID);
                        currentNode.AddNeighbor(currentStartNodeID);
                        
                        currentDownSampledSegment.AddNode(currentNodeID);
                        
                        graph.AddNode(currentNode);
                        
                        currentStartNodeID = currentNodeID;
                        currentStartNode = currentNode;
                        removeNodeList.clear();
                    }
                }
                //add last node
                graph.AddSegment(currentDownSampledSegment);
            }
            else
            {
                typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
                for(; nodeIter != currentSegment.End(); ++nodeIter)
                {
                    NodeIDType currentNodeID = *nodeIter;
                    NodeType currentNode;
                    try
                    {
                        currentNode = graph.GetNode(currentNodeID);
                    }
                    catch(...)
                    {
                        currentNode = this->GetNode(currentNodeID);
                    }
                    graph.AddNode(currentNode);
                }
                graph.AddSegment(currentSegment);
            }
            
        }
    }
    catch(...)
    {
    
    }

    //we could also use the old segment parameters
    graph.ComputeSegmentParameters();
    
    return(graph);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Self 
VesselGraph<VIndexDimension, TCoordinateType>::GetGraphInSphere(const NodeIDType& centerNodeID, RealType radius) const
{
    Self graph;
    typedef std::list<NodeType> NodeListType;
    
    try
    {
        NodeType centerNode = this->GetNode(centerNodeID);
        
        //find nodes close to center node
        NodeListType closeNodeList;
        typename NodeTableType::const_iterator nodeIter = nodeTable.begin();
        for(; nodeIter != nodeTable.end(); ++nodeIter)
        {
            NodeType currentNode = nodeIter->second;
            RealType currentDistance = centerNode.Distance(currentNode);
            if(currentDistance <= radius)
            {
                closeNodeList.push_back(currentNode);
            }
        }

        //form graph with close nodes       
        typename NodeListType::const_iterator closeNodeIter = closeNodeList.begin();
        for(; closeNodeIter != closeNodeList.end(); ++closeNodeIter)
        {
            NodeType currentNode = *closeNodeIter;
            typename NodeType::SegmentIDListType currentSegmentIDList = currentNode.GetSegments();
            typename NodeType::SegmentIDListType::const_iterator currentSegmentIDIter = currentSegmentIDList.begin();
            for(; currentSegmentIDIter != currentSegmentIDList.end(); ++currentSegmentIDIter)
            {
                SegmentIDType currentSegmentID = *currentSegmentIDIter;
                SegmentType currentSegment = this->GetSegment(currentSegmentID);
        
                if(graph.IsSegmentUnique(currentSegment) == true)
                {
                    graph.AddSegment(currentSegment);
                    typename SegmentType::ConstIterator currentSegmentIter = currentSegment.Begin();
                    for(; currentSegmentIter != currentSegment.End(); ++currentSegmentIter)
                    {
                        NodeIDType currentNodeID = *currentSegmentIter;
                        NodeType currentNode = this->GetNode(currentNodeID);
                        graph.AddNode(currentNode);
                    }
                }
                
            }
        }
    }
    catch(...)
    {
    
    }

    return(graph);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Self 
VesselGraph<VIndexDimension, TCoordinateType>::GetSubGraph(const SegmentIDType& id) const
{
    Self subGraph;
    
    this->AddSegmentAndNeighbors(id, subGraph);

    return(subGraph);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Self 
VesselGraph<VIndexDimension, TCoordinateType>::RemoveSubGraph(const SegmentIDType& id) const
{
    Self graphNoSub;
    
    Self subGraph = this->GetSubGraph(id);
    try
    {
        ConstIterator segmentIter = this->Begin();
        for(; segmentIter != this->End(); ++segmentIter)
        {
            SegmentIDType currentSegmentID = segmentIter->first;
            SegmentType currentSegment = segmentIter->second;
            
            try
            {
                subGraph.GetSegment(currentSegmentID);
            }
            catch(...)
            {
                graphNoSub.AddSegment(currentSegment);
                typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
                for(; nodeIter != currentSegment.End(); ++nodeIter)
                {
                    NodeIDType currentNodeID = *nodeIter;
                    NodeType currentNode = this->GetNode(currentNodeID);
                    graphNoSub.AddNode(currentNode);
                }
                
            }
        }
    }
    catch(...)
    {
    
    }

    return(graphNoSub);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::RemoveSegment(const SegmentIDType& id)
{
    try
    {   
        SegmentType segment = this->GetSegment(id);
        NodeIDType frontNodeID = segment.Front();
        NodeType frontNode = this->GetNode(frontNodeID);
        NodeIDType backNodeID = segment.Back();
        NodeType backNode = this->GetNode(backNodeID);
        
        //remove segment and all nodes unless we have junction points
        this->segmentTable.erase(segment.GetID());
        typename SegmentType::ConstIterator nodeIter = segment.Begin();
        
        for(; nodeIter != segment.End(); ++nodeIter)
        {
            NodeIDType currentNode = nodeIter->second;
            if(currentNode.GetType() != NodeType::JUNCTION)
            {
                this->nodeTable.erase(currentNode.GetID());
            }
        }
        
        //is front node a junction
        //if so then remove deleted neighbor and join other segments
        nodeIter = segment.Begin();
        if(frontNode.GetType() == NodeType::JUNCTION)
        {
            std::advance(nodeIter, 1);
            frontNode.RemoveNeighbor(*nodeIter);
            frontNode.RemoveSegment(id);
            this->AddNode(frontNode);
            
            //Get neighbors and joint
            if(frontNode.GetNumberOfNeighbors() == 3)
            {
                typename NodeType::SegmentIDListType joinedSegments;
                typename NodeType::SegmentIDListType currentSegments = frontNode.GetSegments();
                typename NodeType::SegmentIDListType::const_iterator nodeSegIter = currentSegments.begin();
                for(; nodeSegIter != currentSegments.end(); ++nodeSegIter)
                {
                    if(*nodeSegIter != id)
                    {
                        joinedSegments.push_back(*nodeSegIter);
                    }
                }
                
                this->JoinSegments(*(joinedSegments.begin()), *(std::advance(joinedSegments.begin(),1)));
            }
        }
        
        nodeIter = segment.End();
        if(backNode.GetType() == NodeType::JUNCTION)
        {
            std::advance(nodeIter, -2);
            backNode.RemoveNeighbor(*nodeIter);
            backNode.RemoveSegment(id);
            this->AddNode(backNode);
            
            if(backNode.GetNumberOfNeighbors() == 3)
            {
                typename NodeType::SegmentIDListType joinedSegments;
                typename NodeType::SegmentIDListType currentSegments = backNode.GetSegments();
                typename NodeType::SegmentIDListType::const_iterator nodeSegIter = currentSegments.begin();
                for(; nodeSegIter != currentSegments.end(); ++nodeSegIter)
                {
                    if(*nodeSegIter != id)
                    {
                        joinedSegments.push_back(*nodeSegIter);
                    }
                }
                
                this->JoinSegments(*(joinedSegments.begin()), *(std::advance(joinedSegments.begin(),1)));
            }
        }
        
        
    }
    catch(...)
    {
    
    }
    this->ComputeNodeParameters();
    this->ComputeSegmentParameters();
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::JoinSegments(const SegmentIDType& id1, const SegmentIDType& id2)
{
    try
    {
        SegmentType segment1 = this->GetSegment(id1);
        NodeIDType front1NodeID = segment1.Front();
        NodeType front1Node = this->GetNode(front1NodeID);
        NodeIDType back1NodeID = segment1.Back();
        NodeType back1Node = this->GetNode(back1NodeID);
        
        SegmentType segment2 = this->GetSegment(id1);
        NodeIDType front2NodeID = segment2.Front();
        NodeType front2Node = this->GetNode(front2NodeID);
        NodeIDType back2NodeID = segment2.Back();
        NodeType back2Node = this->GetNode(back2NodeID);
        
        RealType ff = front1Node.Distance(front2Node);
        RealType fb = front1Node.Distance(back2Node);
        RealType bf = back1Node.Distance(front2Node);
        RealType bb = back1Node.Distance(back2Node);
        

        if(ff <= fb && ff <= bf && ff <= bb)
        {
            front1Node.AddNeighbor(front2NodeID);
            this->AddNode(front1Node);
            front2Node.AddNeighbor(front1NodeID);
            this->AddNode(front2Node);
        }
        else if(fb <= ff && fb <= bf && fb <= bb)
        {
            front1Node.AddNeighbor(back2NodeID);
            this->AddNode(front1Node);
            back2Node.AddNeighbor(front1NodeID);
            this->AddNode(back2Node);
        }
        else if(fb <= ff && bf <= fb && bf <= bb)
        {
            back1Node.AddNeighbor(front2NodeID);
            this->AddNode(back1Node);
            front2Node.AddNeighbor(back1NodeID);
            this->AddNode(front2Node);
        }
        else
        {
            back1Node.AddNeighbor(back2NodeID);
            this->AddNode(back1Node);
            back2Node.AddNeighbor(back1NodeID);
            this->AddNode(back2Node);
        }
      
        typename SegmentType::ConstIterator nodeIter = segment2.Begin();
        for(; nodeIter != segment2.End(); ++nodeIter)
        {
            NodeType currentNode = nodeIter->second;
            currentNode.RemoveSegment(segment2.GetID());
            currentNode.AddSegment(segment1.GetID());
            this->AddNode(currentNode);
        }
        this->segmentTable.erase(segment2);
        
        
    }
    catch(...)
    {
    
    }
    
    this->ComputeNodeParameters();
    this->ComputeSegmentParameters();
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::AddSegmentAndNeighbors(const typename SegmentType::IDType& segmentID, Self& graph) const
{
    //add segment and nodes to graph
    SegmentType currentSegment = this->GetSegment(segmentID);
    
    if(graph.IsSegmentUnique(currentSegment) == true)
    {
        graph.AddSegment(currentSegment);
        typename SegmentType::ConstIterator currentSegmentIter = currentSegment.Begin();
        for(; currentSegmentIter != currentSegment.End(); ++currentSegmentIter)
        {
            NodeIDType currentNodeID = *currentSegmentIter;
            NodeType currentNode = this->GetNode(currentNodeID);
            graph.AddNode(currentNode);
        }
        
        //go to front of segment and see if there are new segments 
        typename SegmentType::NodeIDType frontNodeID = currentSegment.Front();
        NodeType frontNode = this->GetNode(frontNodeID);
        
        if(frontNode.GetType() == NodeType::JUNCTION)
        {
            typename NodeType::SegmentIDListType frontNeighborList = frontNode.GetNeighbors();
            typename NodeType::SegmentIDListType::const_iterator frontNeighborIter = frontNeighborList.begin();
            for(; frontNeighborIter != frontNeighborList.end(); ++frontNeighborIter)
            {
                typename NodeType::IDType currentNeigborID = *frontNeighborIter;
                NodeType currentNeighbor = this->GetNode(currentNeigborID);

                typename NodeType::SegmentIDListType currentSegmentIDList = currentNeighbor.GetSegments();
                typename NodeType::SegmentIDListType::const_iterator currentSegmentIter = currentSegmentIDList.begin();
                for(; currentSegmentIter != currentSegmentIDList.end(); ++currentSegmentIter)
                {
                    this->AddSegmentAndNeighbors(*currentSegmentIter, graph);
                }
            
            }
        }
        
        //go to back of segment
        typename SegmentType::NodeIDType backNodeID = currentSegment.Back();
        NodeType backNode = this->GetNode(backNodeID);
        
        if(backNode.GetType() == NodeType::JUNCTION)
        {
            typename NodeType::SegmentIDListType backNeighborList = backNode.GetNeighbors();
            typename NodeType::SegmentIDListType::const_iterator backNeighborIter = backNeighborList.begin();
            for(; backNeighborIter != backNeighborList.end(); ++backNeighborIter)
            {
                typename NodeType::IDType currentNeigborID = *backNeighborIter;
                NodeType currentNeighbor = this->GetNode(currentNeigborID);
                
                typename NodeType::SegmentIDListType currentSegmentIDList = currentNeighbor.GetSegments();
                typename NodeType::SegmentIDListType::const_iterator currentSegmentIter = currentSegmentIDList.begin();
                for(; currentSegmentIter != currentSegmentIDList.end(); ++currentSegmentIter)
                {
                    this->AddSegmentAndNeighbors(*currentSegmentIter, graph);
                }
            
            }
        }
        
        
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Iterator 
VesselGraph<VIndexDimension, TCoordinateType>::Begin()
{
    return(this->segmentTable.begin());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Iterator 
VesselGraph<VIndexDimension, TCoordinateType>::End()
{
    return(this->segmentTable.end());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::ConstIterator 
VesselGraph<VIndexDimension, TCoordinateType>::Begin() const
{
    return(this->segmentTable.begin());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::ConstIterator 
VesselGraph<VIndexDimension, TCoordinateType>::End() const
{
    return(this->segmentTable.end());
}


template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::FormGraph()
{
   
    //Get segments
    //segment is defined as a list of indices where the beginning and end
    //points have either only one neighbor or strictly more than
    //two neighbors.
    //terminal to terminal
    //terminal to junction (vica versa)
    //junction to junction
    //This algorithm will miss closed loops
    
    //remove neighbor ids that are not in list
    typename NodeTableType::const_iterator nodeIter = nodeTable.begin();
    for(; nodeIter != nodeTable.end(); ++nodeIter)
    {
        NodeType currentNode = nodeIter->second;
        typename NodeType::NodeIDListType currentNeighborIDList = currentNode.GetNeighbors();
        typename NodeType::NodeIDListType::const_iterator neighborIter = currentNeighborIDList.begin();
        for(; neighborIter != currentNeighborIDList.end(); ++neighborIter)
        {
            NodeIDType currentNeighborID = *(neighborIter);
            try
            {
                this->GetNode(currentNeighborID);
            }
            catch(...)
            {
                currentNode.RemoveNeighbor(currentNeighborID);
                this->nodeTable[currentNode.GetID()] = currentNode;
            }
        }
    }
    
    unsigned long newSegmentCounter = 1;
    
    nodeIter = nodeTable.begin();
    for(; nodeIter != nodeTable.end(); ++nodeIter)
    {
        NodeType candidateNode = nodeIter->second;
        NodeIDType candidateID = candidateNode.GetID();
        typename NodeType::NodeType candidateType = candidateNode.GetType();
        
        //process if terminal or junction
        if(candidateType == NodeType::JUNCTION || candidateType == NodeType::TERMINAL)
        {
            typename NodeType::NodeIDListType candidateNeighborIDList = candidateNode.GetNeighbors();
            typename NodeType::NodeIDListType::const_iterator neighborIter = candidateNeighborIDList.begin();
            for(; neighborIter != candidateNeighborIDList.end(); ++neighborIter)
            {
                SegmentType candidateSegment;
                candidateSegment.AddNode(candidateNode.GetID());
                
                NodeIDType previousID = candidateNode.GetID();
                NodeIDType currentID = *(neighborIter);
                bool keepFindingNode = true; 
                
                for(unsigned int i = 0; i < 10000 && keepFindingNode == true; i++)
                {
                    try
                    {
                        NodeType currentNode = this->GetNode(currentID);
                        candidateSegment.AddNode(currentNode.GetID());
                        typename NodeType::NodeType currentType = currentNode.GetType();
                        
                        if(currentType == NodeType::MIDSEGMENT)
                        {
                            //figure out next index by avoiding where you came from
                            NodeIDType nextID;
                            if(currentNode.GetNeighbors().front() != previousID)
                            {
                                nextID = currentNode.GetNeighbors().front();
                            }
                            else
                            {
                                nextID = currentNode.GetNeighbors().back();
                            }
                            
                            previousID = currentID;
                            currentID = nextID;
                            keepFindingNode = true;
                        }
                        else
                        {
                            keepFindingNode = false;
                        }
                    }
                    catch(...)
                    {
                        //couldn't find node
                        std::cout << "HELP!" << std::endl;
                        break;
                    }
                }
                
                if(this->IsSegmentUnique(candidateSegment) == true)
                {
                    candidateSegment.SetID(newSegmentCounter);
                    this->AddSegment(candidateSegment);
                    this->UpdateSegmentIDOnNodes(newSegmentCounter);
                    newSegmentCounter++;
                }
                
            }
        }
        else if(candidateType == NodeType::SINGLE)
        {
            SegmentType candidateSegment;
            candidateSegment.AddNode(candidateNode.GetID());
            
            if(this->IsSegmentUnique(candidateSegment) == true)
            {
                candidateSegment.SetID(newSegmentCounter);
                this->AddSegment(candidateSegment);
                this->UpdateSegmentIDOnNodes(newSegmentCounter);
                newSegmentCounter++;
            }
            
        }
        //else midsegment do nothing
    }
    
    //clean up graph
    //this->FixPoorRadii();
    //this->RemoveTriangleSegments();
    this->RemoveExcessNodes();
    this->ComputeNodeParameters();
    this->ComputeSegmentParameters();
    
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::Self VesselGraph<VIndexDimension, TCoordinateType>::RemoveRedundentJunctions() const 
{
    Self newGraph;
            
    NodeTableType newNodeTable(this->nodeTable);
    SegmentTableType newSegmentTable(this->segmentTable);
    
    SegmentIDType newSegmentID = 0;
    typename SegmentTableType::const_iterator segmentTableIter = this->segmentTable.begin();
    for(; segmentTableIter != this->segmentTable.begin(); ++segmentTableIter)
    {
        SegmentIDType currentSegmentID = segmentTableIter->first;
        if(currentSegmentID > newSegmentID)
        {
            newSegmentID = currentSegmentID;
        }
    }
    newSegmentID++;
       
    try
    {
        NodeIDType newNodeID = 0;
        typename NodeTableType::const_iterator nodeIter = this->nodeTable.begin();
        for(; nodeIter != this->nodeTable.end(); ++nodeIter)
        {
            NodeIDType currentID = nodeIter->first;
            if(currentID >= newNodeID)
            {
                newNodeID = currentID;
            }
        }
        ++newNodeID;
        
        

        
        nodeIter = this->nodeTable.begin();
        for(; nodeIter != this->nodeTable.end(); ++nodeIter)
        {
            NodeType currentNode = nodeIter->second;
             
            NodeTableType currentJunctionNeighbors;
            
            //get all junction points in a cluster
            this->GetNeighborJunctions(currentNode, currentJunctionNeighbors);
           
            if(currentJunctionNeighbors.size() > 1)
            {   
                //find new center 
                typename NodeType::PositionType averagePosition;
                averagePosition.Fill(0.0);
                typename NodeTableType::const_iterator neighborIter = currentJunctionNeighbors.begin();
                for(; neighborIter != currentJunctionNeighbors.end(); ++neighborIter)
                {
                    averagePosition = averagePosition + (neighborIter->second).GetPosition();
                }
                averagePosition = 1.0/RealType(currentJunctionNeighbors.size()) * averagePosition;
                NodeType newNode;
                newNode.SetID(newNodeID);
                ++newNodeID;
                newNode.SetPosition(averagePosition);
                newNode.SetRadius(1.0);
                
                
                neighborIter = currentJunctionNeighbors.begin();
                for(; neighborIter != currentJunctionNeighbors.end(); ++neighborIter)
                {
                    NodeIDType currentNeighborNodeID = neighborIter->first;
                    NodeType currentNeighborNode = neighborIter->second;
                    typename NodeType::SegmentIDListType currentNeighborSegmentList = currentNeighborNode.GetSegments();
                    typename NodeType::SegmentIDListType::const_iterator neighSegIter = currentNeighborSegmentList.begin();
                    for(; neighSegIter != currentNeighborSegmentList.end(); ++neighSegIter)
                    {
                        SegmentIDType currentNeighSegID = *neighSegIter;
                        SegmentType  currentNeighSeg = this->GetSegment(currentNeighSegID);
                        NodeIDType frontID = currentNeighSeg.Front();
                        NodeType frontNode = this->GetNode(frontID);
                        NodeIDType backID = currentNeighSeg.Back();
                        NodeType backNode = this->GetNode(backID);
                        //remove segment
                        if(currentNeighSeg.NumberOfNodes() == 2 && currentNeighborNodeID == frontID && backNode.GetType() == NodeType::JUNCTION)
                        {
                            newSegmentTable.erase(currentNeighSegID);  
                            newNodeTable.erase(currentNeighborNodeID);    
                            newNodeTable.erase(backID);                     
                        }
                        else if(currentNeighSeg.NumberOfNodes() == 2 && currentNeighborNodeID == backID && frontNode.GetType() == NodeType::JUNCTION)
                        {
                            newSegmentTable.erase(currentNeighSegID);     
                            newNodeTable.erase(currentNeighborNodeID); 
                            newNodeTable.erase(frontID);                     
                        }
                    }
                    
                    typename NodeType::NodeIDListType currentNeighborNodeList = currentNeighborNode.GetNeighbors();
                    typename NodeType::NodeIDListType::const_iterator neighNodeIter = currentNeighborNodeList.begin();
                    for(; neighNodeIter != currentNeighborNodeList.end(); ++neighNodeIter)
                    {
                        NodeIDType currentMidNodeID = *(neighNodeIter);
                        NodeType currentMidNode = this->GetNode(currentMidNodeID);
                        if(currentMidNode.GetType() == NodeType::MIDSEGMENT || currentMidNode.GetType() == NodeType::TERMINAL)
                        {
                            SegmentIDType currentMidSegmentID = *(currentMidNode.GetSegments().begin());
                            SegmentType currentMidSegment = newSegmentTable[currentMidSegmentID];
                            currentMidSegment.ReplaceNode(currentNeighborNodeID, newNode.GetID());
                            newSegmentTable[currentMidSegment.GetID()] = currentMidSegment;
                            
                            currentMidNode.RemoveNeighbor(currentNeighborNodeID);
                            currentMidNode.AddNeighbor(newNode.GetID());
                            newNodeTable[currentMidNode.GetID()] = currentMidNode;
                            
                            
                            newNode.AddNeighbor(currentMidNode.GetID());
                            newNode.AddSegment(currentMidSegmentID); 
                        }
                        
                    }

                }
                
                //remove old cluster points nodes
                //add later after neighbors and segment ids have been set
                newNodeTable[newNode.GetID()] = newNode;
                
            }
            
            
        }
        
    }
    catch(...)
    {
    
    }
    
    newGraph.nodeTable = newNodeTable;
    newGraph.segmentTable = newSegmentTable;

    newGraph.RemoveExcessNodes();
    return(newGraph);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::RemoveExcessNodes()
{
    typedef std::list<NodeIDType> NodeIDListType;
    NodeIDListType nodeIDList;
    
    ConstIterator segmentIter = this->Begin();
    for(; segmentIter != this->End(); ++segmentIter)
    {
        typename SegmentType::ConstIterator nodeIter = (segmentIter->second).Begin();
        for(; nodeIter != (segmentIter->second).End(); ++nodeIter)
        {
            NodeIDType currentNodeID = *nodeIter;
            typename NodeIDListType::const_iterator findIter = std::find(nodeIDList.begin(), nodeIDList.end(), currentNodeID);
            if(findIter == nodeIDList.end())
            {
                nodeIDList.push_back(currentNodeID);
            }
        }
    }
    
    NodeTableType newNodeTable(this->nodeTable);
    typename NodeTableType::const_iterator nodeIter = this->nodeTable.begin();
    for(; nodeIter != this->nodeTable.end(); nodeIter++)
    {
        NodeIDType currentNodeID = nodeIter->first;
        typename NodeIDListType::const_iterator findIter = std::find(nodeIDList.begin(), nodeIDList.end(), currentNodeID);
        if(findIter == nodeIDList.end())
        {
            newNodeTable.erase(currentNodeID);
        }
    }
    this->nodeTable = newNodeTable;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::GetNeighborJunctions(const NodeType node, NodeTableType& nodeList) const
{

    try
    {
        if(node.GetType() == NodeType::JUNCTION)
        {
            nodeList[node.GetID()] = node;
            
            typename NodeType::NodeIDListType neighbors = node.GetNeighbors();
            typename NodeType::NodeIDListType::const_iterator neighIter = neighbors.begin();
            for(; neighIter != neighbors.end(); ++neighIter)
            {
                NodeIDType neighborNodeID = *neighIter;
                
                if(nodeList.find(neighborNodeID) == nodeList.end())
                {
                    NodeType neighborNode = this->GetNode(neighborNodeID);
                    this->GetNeighborJunctions(neighborNode, nodeList);
                }
            }
        }
    }
    catch(...)
    {
    
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::FixPoorRadii(RealType medialnessThreshold)
{
    Iterator segmentIter = this->Begin();
    for(; segmentIter != this->End(); ++segmentIter)
    {
        //get max/min medialness and associated radius in segment
        RealType maxMedialness = 0.0;
        RealType radiusAtMaxMedialness = 0.0;
        typename SegmentType::Iterator nodeIter = (segmentIter->second).Begin();
        for(; nodeIter != (segmentIter->second).End(); ++nodeIter)
        {
            NodeIDType currentNodeID = *nodeIter;
            NodeType currentNode = this->GetNode(currentNodeID);
            typename NodeType::RealType currentMedialness = currentNode.GetMedialness();
            typename NodeType::RealType currentRadius = currentNode.GetRadius();
            if(currentMedialness >= maxMedialness)
            {
                maxMedialness = currentMedialness;
                radiusAtMaxMedialness = currentRadius;
            }
        }
        nodeIter = (segmentIter->second).Begin();
        for(; nodeIter != (segmentIter->second).End(); ++nodeIter)
        {
            NodeIDType currentNodeID = *nodeIter;
            NodeType currentNode = this->GetNode(currentNodeID);
            typename NodeType::RealType currentMedialness = currentNode.GetMedialness();
            if(currentMedialness < medialnessThreshold)
            {
                typename SegmentType::ConstIterator nodeForwardIter = nodeIter;
                std::advance(nodeForwardIter, 1);
                typename SegmentType::ConstIterator nodeBackwardIter = nodeIter;
                std::advance(nodeBackwardIter, -1);
                
                typename NodeType::RealType forwardMedialness = 0.0;
                typename NodeType::RealType backwardMedialness = 0.0;
                typename NodeType::RealType forwardRadius = 0.0;
                typename NodeType::RealType backwardRadius = 0.0;
                bool foundIt = false;
                unsigned int i = 0;
                for(; i < (segmentIter->second).NumberOfNodes(); i++)
                {
                    if(nodeForwardIter != (segmentIter->second).End())
                    {
                        NodeIDType forwardNodeID = *nodeForwardIter;
                        NodeType forwardNode = this->GetNode(forwardNodeID);
                        forwardMedialness = forwardNode.GetMedialness();
                        forwardRadius = forwardNode.GetRadius(); 
                        std::advance(nodeForwardIter, 1);
                    }
                    typename SegmentType::ConstIterator rbeginIter= (segmentIter->second).Begin();
                    std::advance(rbeginIter, -1);                    
                    if(nodeBackwardIter != rbeginIter)
                    {
                        NodeIDType backwardNodeID = *nodeBackwardIter;
                        NodeType backwardNode = this->GetNode(backwardNodeID);
                        backwardMedialness = backwardNode.GetMedialness();
                        backwardRadius = backwardNode.GetRadius();
                        std::advance(nodeBackwardIter, -1);
                    }
                    if(forwardMedialness >= backwardMedialness && forwardMedialness > medialnessThreshold)
                    {
                        currentNode.SetRadius(forwardRadius);
                        this->nodeTable[currentNode.GetID()] = currentNode;
                        foundIt = true;
                        break;
                    }
                    else if(backwardMedialness > forwardMedialness && backwardMedialness > medialnessThreshold)
                    {
                        currentNode.SetRadius(backwardRadius); 
                        this->nodeTable[currentNode.GetID()] = currentNode;
                        
                        foundIt = true;
                        break;
                    }
                    
                }
                if(foundIt == true)
                {
                    currentNode.SetRadius(radiusAtMaxMedialness);
                    this->nodeTable[currentNode.GetID()] = currentNode;
                }
                else
                {
                    //consider removing segment
                }
                
            }
        }
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::UpdateSegmentIDOnNodes(SegmentIDType id)
{
    try
    {
        SegmentType segment = this->GetSegment(id);
        typename SegmentType::Iterator iter = segment.Begin();
        for(; iter != segment.End(); ++iter)
        {
            NodeType currentNode = this->GetNode(*iter);
            currentNode.AddSegment(id);
            this->nodeTable[currentNode.GetID()] = currentNode;
        }
    }
    catch(...)
    {
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
unsigned int VesselGraph<VIndexDimension, TCoordinateType>::NumberOfSegments() const
{
    return(this->segmentTable.size());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
unsigned int VesselGraph<VIndexDimension, TCoordinateType>::NumberOfJunctions() const
{
    unsigned int numJunctions = 0;
    typename NodeTableType::const_iterator nodeIter = nodeTable.begin();
    for(; nodeIter != nodeTable.end(); ++nodeIter)
    {
        NodeType currentNode = nodeIter->second;

        if(currentNode.GetType() == NodeType::JUNCTION)
        {
            ++numJunctions;
        }
    }
    
    return(numJunctions);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
unsigned int VesselGraph<VIndexDimension, TCoordinateType>::NumberOfNJunctions(unsigned int n) const
{
    unsigned int numJunctions = 0;
    
    typename NodeTableType::const_iterator nodeIter = nodeTable.begin();
    for(; nodeIter != nodeTable.end(); ++nodeIter)
    {
        if(nodeIter->second.GetNumberOfNeighbors() == n)
        {
            ++numJunctions;
        }
    }
    
    return(numJunctions);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
unsigned int VesselGraph<VIndexDimension, TCoordinateType>::NumberOfNodes(const SegmentIDType& id) const
{
    unsigned int numNodes = 0;
    try
    {
        SegmentType segment = this->GetSegment(id);
        numNodes = segment.NumberOfNodes();
    }
    catch(...)
    {
    }
    
    return(numNodes);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentLength(SegmentType& segment) const
{
    RealType distance = 0.0;
    try
    {
        typename SegmentType::Iterator iter = segment.Begin();
        for(; iter != segment.End(); ++iter)
        {
            NodeIDType currentID = *(iter);
            NodeType currentNode = this->GetNode(currentID);

            typename SegmentType::Iterator nextIter = iter;
            std::advance(nextIter, 1);
            if(nextIter != segment.End())
            {
                NodeIDType nextID = *(nextIter);
                NodeType nextNode = this->GetNode(nextID);
                distance += currentNode.Distance(nextNode);
            }
        }
    }
    catch(...)
    {
    }
    
    segment.SetLength(distance);
}

template<unsigned int VIndexDimension, typename TCoordinateType>   
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentChordLength( SegmentType& segment ) const
{
    RealType distance = 0.0;
    try
    {
        NodeIDType frontID = segment.Front();
        NodeType frontNode = this->GetNode(frontID);
        NodeIDType backID = segment.Back();
        NodeType backNode = this->GetNode(backID);
        distance = frontNode.Distance(backNode); 
    }
    catch(...)
    {
    }
    
    segment.SetChordLength(distance);
}

template<unsigned int VIndexDimension, typename TCoordinateType>   
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentSumOfTotalAngle( SegmentType& segment ) const
{
    RealType sumOfTotalAngle = 0.0;
    try
    {
        typename SegmentType::Iterator iter = segment.Begin();
        for(; iter != segment.End(); ++iter)
        {
            NodeIDType currentID = *(iter);
            NodeType currentNode = this->GetNode(currentID);
            sumOfTotalAngle += currentNode.GetTotalAngle();
        }
        double length = segment.GetLength();
        if(length > 0)
        {
            sumOfTotalAngle = sumOfTotalAngle/length;
        }
    }
    catch(...)
    {
    }
    
    segment.SetSumOfAngleMetric(sumOfTotalAngle);
}

template<unsigned int VIndexDimension, typename TCoordinateType>   
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentAverageRadius( SegmentType& segment ) const
{
    RealType averageRadius = 0.0;
    try
    {
        typename SegmentType::Iterator iter = segment.Begin();
        for(; iter != segment.End(); ++iter)
        {
            NodeIDType currentID = *(iter);
            NodeType currentNode = this->GetNode(currentID);
            averageRadius += currentNode.GetRadius();
        }
        unsigned int numNodes = segment.NumberOfNodes();
        if(numNodes > 0)
        {
            averageRadius = averageRadius/numNodes;
        }
    }
    catch(...)
    {
    }
    
    segment.SetAverageRadius(averageRadius);
}

template<unsigned int VIndexDimension, typename TCoordinateType>   
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentSumMagnitudeOfCurvature( SegmentType& segment) const
{
    RealType sumCurvature = 0.0;
    RealType sumCurvatureSquared = 0.0;
    try
    {
        typename SegmentType::Iterator iter = segment.Begin();
        for(; iter != segment.End(); ++iter)
        {
            NodeIDType currentID = *(iter);
            NodeType currentNode = this->GetNode(currentID);
            RealType currentCurvature = currentNode.GetCurvature();
            sumCurvature += fabs(currentCurvature);
            sumCurvatureSquared +=  currentCurvature*currentCurvature;
        }
    }
    catch(...)
    {
    }
    
    segment.SetSumMagnitudeOfCurvature(sumCurvature);
    segment.SetSumMagnitudeOfCurvatureSquared(sumCurvatureSquared);
}

template<unsigned int VIndexDimension, typename TCoordinateType>   
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentSumMagnitudeOfDerivativeCurvature(SegmentType& segment) const
{   
    RealType sumDCurvature = 0.0;
    
    try
    {
        typename SegmentType::Iterator iter = segment.Begin();
        for(; iter != segment.End(); ++iter)
        {
            NodeIDType currentID = *(iter);
            NodeType currentNode = this->GetNode(currentID);
            RealType currentCurvature = currentNode.GetCurvature();
            
            typename SegmentType::Iterator nextIter = iter;
            std::advance(nextIter, 1);
            
            if(nextIter != segment.End())
            {
                NodeIDType nextID = *(nextIter);
                NodeType nextNode = this->GetNode(nextID);
                RealType nextCurvature = nextNode.GetCurvature();
                
                sumDCurvature += fabs(currentCurvature-nextCurvature);
            }
            
        }
    }
    catch(...)
    {
    }
    
    segment.SetSumMagnitudeOfDervativeCurvature(sumDCurvature);
}


template<unsigned int VIndexDimension, typename TCoordinateType>   
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentTortuosityDensity( SegmentType& segment) const
{
    RealType tortDensity = 0.0;
   
    try
    {
        RealType subSegmentLength = 0.0;
        int numSubSegments = 0;
        if(segment.NumberOfNodes() >= 2)
        {
            RealType segmentLength = segment.GetLength();
            
            typename SegmentType::Iterator iter = segment.Begin();
            NodeType startNode = this->GetNode(*(iter));
            NodeType stopNode;
            
            for(; iter != segment.End(); ++iter)
            {
                NodeIDType currentID = *(iter);
                NodeType currentNode = this->GetNode(currentID);
                RealType currentCurvature = currentNode.GetCurvature();
                bool currentCurvaturePositive;
                if(currentCurvature >= 0)
                {
                    currentCurvaturePositive = true;
                }
                else
                {
                    currentCurvaturePositive = false;
                }
                
                typename SegmentType::Iterator nextIter = iter;
                std::advance(nextIter, 1);
                
                if(nextIter != segment.End())
                {
                    NodeIDType nextID = *(nextIter);
                    NodeType nextNode = this->GetNode(nextID);
                    RealType nextCurvature = nextNode.GetCurvature();
                    
                    bool nextCurvaturePositive;
                    if(nextCurvature >= 0)
                    {
                        nextCurvaturePositive = true;
                    }
                    else
                    {
                        nextCurvaturePositive = false;
                    }
                    
                    subSegmentLength += currentNode.Distance(nextNode);
                    
                    if(currentCurvaturePositive != nextCurvaturePositive)
                    {
                        stopNode = nextNode;
                        RealType chordLength = startNode.Distance(stopNode);
                        tortDensity += (subSegmentLength/chordLength) - 1.0;
                        
                        startNode = stopNode;
                        subSegmentLength = 0.0;
                        ++numSubSegments;
                    }
                }
                
            }
            
            tortDensity *= (numSubSegments-1)/segmentLength;
        }
    }
    catch(...)
    {
    }
    
    segment.SetTortuosityDensity(tortDensity);
}

template<unsigned int VIndexDimension, typename TCoordinateType>   
void VesselGraph<VIndexDimension, TCoordinateType>::SetSegmentInflectionCount( SegmentType& segment) const
{
    RealType inflectionCount = 0.0;
    
    try
    {
        typename SegmentType::Iterator iter = segment.Begin();
        for(; iter != segment.End(); ++iter)
        {
            NodeIDType currentID = *(iter);
            NodeType currentNode = this->GetNode(currentID);
            typename NodeType::VectorType currentNormal = currentNode.GetNormal();
            
            typename SegmentType::Iterator nextIter = iter;
            std::advance(nextIter, 1);

            if(nextIter != segment.End())
            {
                NodeIDType nextID = *(nextIter);
                NodeType nextNode = this->GetNode(nextID);
                typename NodeType::VectorType nextNormal = nextNode.GetNormal();
                
                typename NodeType::VectorType normalDiff = nextNormal - currentNormal;
                
                inflectionCount += normalDiff.Norm();
            }
            
        }
    }
    catch(...)
    {
    }
    inflectionCount = inflectionCount/RealType(segment.NumberOfNodes());
    segment.SetInflectionCountMetric(inflectionCount);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
bool VesselGraph<VIndexDimension, TCoordinateType>::IsSegmentUnique(SegmentType aSegment) const
{
    ConstIterator segmentIter = this->Begin();
    for(; segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;
        
        typename SegmentType::ConstIterator searchIterForward = std::search(currentSegment.Begin(), currentSegment.End(), aSegment.Begin(), aSegment.End());
        aSegment.Reverse();
        typename SegmentType::ConstIterator searchIterBackward = std::search(currentSegment.Begin(), currentSegment.End(), aSegment.Begin(), aSegment.End());
        
        if(searchIterForward != currentSegment.End() || searchIterBackward != currentSegment.End())
        {
            return(false);
        }
    }
                
    return(true);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::AlignTangentVector()
{
    typename NodeType::VectorType zeroVector;
    zeroVector.Fill(0.0);
    
    for(Iterator segmentIter = this->Begin(); segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;
        NodeIDType currentBackID = currentSegment.Back();
        typename NodeType::RealType currentBackRadius = GetNode(currentBackID).GetRadius();
        NodeIDType currentFrontID = currentSegment.Front();
        typename NodeType::RealType currentFrontRadius = GetNode(currentFrontID).GetRadius();
        
        if( currentBackRadius > currentFrontRadius )
        {
            currentSegment.Reverse();
        }
        this->segmentTable[currentSegment.GetID()] = currentSegment;
        
        //align tangent vectors
        typename SegmentType::Iterator nodeIter = currentSegment.Begin();
        ++nodeIter;
        for(; nodeIter != currentSegment.End(); ++nodeIter)
        {
            NodeIDType currentID = *(nodeIter);
            NodeType currentNode = this->GetNode(currentID);
            typename NodeType::VectorType currentTangent = currentNode.GetTangent();
            
            typename SegmentType::Iterator iterBehind = nodeIter;
            --iterBehind;
            NodeIDType previousID = *(iterBehind);
            typename NodeType::VectorType previousTangent = this->GetNode(previousID).GetTangent();
            
            for(; previousTangent == zeroVector && iterBehind != currentSegment.Begin(); --iterBehind)
            {
                previousID = *(iterBehind);
                previousTangent = this->GetNode(previousID).GetTangent();
            }
            
            RealType dot = dot_product(currentTangent, previousTangent);
            if(dot < 0)
            {
                currentTangent = -currentTangent;
                currentNode.SetTangent(currentTangent);
                this->nodeTable[currentNode.GetID()] = currentNode;
            }
        }
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::ComputeFrenetFrame()
{
   typename NodeType::VectorType zeroVector;
    zeroVector.Fill(0.0);
    
    for(Iterator segmentIter = this->Begin(); segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;
                
        //frenet frame
        typename SegmentType::Iterator nodeIter = currentSegment.Begin();
        ++nodeIter;
        for(; nodeIter != currentSegment.End(); ++nodeIter)
        {
            NodeIDType currentID = *(nodeIter);
            NodeType currentNode = this->GetNode(currentID);
            typename NodeType::VectorType currentTangent = currentNode.GetTangent();
            
            typename SegmentType::Iterator nextIter = nodeIter;
            std::advance(nextIter, 1);
            typename SegmentType::Iterator previousIter = nodeIter;
            std::advance(previousIter, -1);
            
            if(nextIter != currentSegment.End())
            {
                NodeIDType nextID = *(nextIter);
                NodeType nextNode = this->GetNode(nextID);
                typename NodeType::VectorType nextTangent = nextNode.GetTangent();
                
                NodeIDType previousID = *(previousIter);
                NodeType previousNode = this->GetNode(previousID);
                typename NodeType::VectorType previousTangent = previousNode.GetTangent();
                
                typename NodeType::VectorType planeVector = nextTangent - previousTangent;
                planeVector = planeVector.Normalize();
                
                //RealType dotProduct = dot_product(currentTangent, nextTangent);
                RealType dotProduct = dot_product(currentTangent, planeVector);
                
            
                if(fabs(dotProduct) <= 0.999999)
                {
                    //typename NodeType::VectorType currentNormal = (nextTangent - dotProduct*currentTangent).Normalize();
                    typename NodeType::VectorType currentNormal = (planeVector - dotProduct*currentTangent).Normalize();
                    currentNode.SetNormal(currentNormal);
                         
                    typename NodeType::VectorType currentBiNormal= cross_product(currentTangent, currentNormal);
                    if(currentBiNormal.Norm() >= 1.0e-6)
                    {
                        currentBiNormal.Normalize();
                        currentNode.SetBiNormal(currentBiNormal);
                        
                        if(nextIter == currentSegment.Begin())
                        {
                            NodeIDType firstID = *(nextIter);
                            NodeType firstNode = this->GetNode(firstID);
                            firstNode.SetNormal(currentNormal);
                            firstNode.SetBiNormal(currentBiNormal);
                            this->nodeTable[firstNode.GetID()] = firstNode;
                        }
                    }
                    
                    this->nodeTable[currentNode.GetID()] = currentNode;
                    
                }
            }
        }    
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::ComputeFrenetParameters()
{
   typename NodeType::VectorType zeroVector;
    zeroVector.Fill(0.0);
    
    for(Iterator segmentIter = this->Begin(); segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;
     
        //curvature, torsion, angles
        typename SegmentType::Iterator nodeIter = currentSegment.Begin();
        nodeIter++;
        for(; nodeIter != currentSegment.End(); ++nodeIter)
        {
            NodeIDType currentID = *(nodeIter);
            NodeType currentNode = this->GetNode(currentID);
            typename NodeType::VectorType currentTangent = currentNode.GetTangent();
            typename NodeType::VectorType currentNormal = currentNode.GetNormal();
            typename NodeType::VectorType currentBiNormal = currentNode.GetBiNormal();
            
            typename SegmentType::Iterator nextIter = nodeIter;
            std::advance(nextIter, 1);
            
            if(nextIter != currentSegment.End())
            {
                NodeIDType nextID = *(nextIter);
                NodeType nextNode = this->GetNode(nextID);
                typename NodeType::VectorType nextTangent = nextNode.GetTangent();
                typename NodeType::VectorType nextNormal = nextNode.GetNormal();
                typename NodeType::VectorType nextBiNormal = nextNode.GetBiNormal();
                
                typename NodeType::VectorType currentTangentDeriv = nextTangent - currentTangent;
                typename NodeType::RealType currentCurvature = dot_product(currentTangentDeriv, currentNormal);
                
                typename NodeType::VectorType currentNormalDeriv = nextNormal - currentNormal;
                typename NodeType::RealType currentTorsion = dot_product(currentNormalDeriv, currentBiNormal);
                
                typename NodeType::RealType currentInPlaneAngle = 0.0;
                typename NodeType::RealType currentInPlaneDot = dot_product(currentTangent, nextTangent);
                if(fabs(currentInPlaneDot) <= 1.0 && currentTangent.Norm() > 0.0 && nextTangent.Norm() > 0.0)
                {
                    currentInPlaneAngle = acos(currentInPlaneDot);
                }
                
                typename NodeType::RealType currentTorsionalAngle = 0.0;
                typename NodeType::RealType currentTorsionalDot = dot_product(currentNormal, nextNormal);
                if(fabs(currentTorsionalDot) <= 1.0 && currentNormal.Norm() > 0.0 && nextNormal.Norm() > 0.0)
                {
                    currentTorsionalAngle = acos(currentTorsionalDot);
                }
                
                currentNode.SetCurvature(currentCurvature);
                currentNode.SetTorsion(currentTorsion);
                currentNode.SetInPlaneAngle(currentInPlaneAngle);
                currentNode.SetTorsionalAngle(currentTorsionalAngle);
                
                this->nodeTable[currentNode.GetID()] = currentNode;
            }
        }    
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::ComputeNodeParameters()
{
    this->AlignTangentVector();
    this->ComputeFrenetFrame();
    this->ComputeFrenetParameters();
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::HistogramType 
VesselGraph<VIndexDimension, TCoordinateType>::GetHistogram() const
{
    int maxNumberConnections = 0;

    typename NodeTableType::const_iterator nodeIter = this->nodeTable.begin();
    for(; nodeIter != nodeTable.end(); ++nodeIter)
    {
        NodeType currentNode = nodeIter->second;
        unsigned int numNeigh = currentNode.GetNumberOfNeighbors();
        if(numNeigh > maxNumberConnections)
        {
            maxNumberConnections = numNeigh;
        }
    }
    
    maxNumberConnections = maxNumberConnections+1;
    
    HistogramAxisType xAxis(maxNumberConnections);
    HistogramAxisType yAxis(maxNumberConnections);
    for(unsigned int i = 0; i < maxNumberConnections; i++)
    {
        xAxis[i] = i;
        yAxis[i] = 0;
    }
    
    nodeIter = this->nodeTable.begin();
    for(; nodeIter != nodeTable.end(); ++nodeIter)
    {
        NodeType currentNode = nodeIter->second;
        unsigned int numNeigh = currentNode.GetNumberOfNeighbors();
        yAxis[numNeigh] = yAxis[numNeigh] + 1;

    }
    
    HistogramType histogram;
    histogram.first = xAxis;
    histogram.second = yAxis;
    
    return(histogram);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::ComputeSegmentParameters()
{
    Iterator segmentIter = this->Begin();
    for(; segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;

        this->SetSegmentLength(currentSegment);
        this->SetSegmentChordLength(currentSegment);        
        this->SetSegmentSumOfTotalAngle(currentSegment);        
        this->SetSegmentAverageRadius(currentSegment);        
        this->SetSegmentSumMagnitudeOfCurvature(currentSegment);        
        this->SetSegmentSumMagnitudeOfDerivativeCurvature(currentSegment);        
        this->SetSegmentTortuosityDensity(currentSegment);
        this->SetSegmentInflectionCount(currentSegment);

        this->segmentTable[currentSegment.GetID()] = currentSegment;
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::SetTangentByPosition()
{
    typename NodeType::VectorType zeroVector;
    zeroVector.Fill(0.0);
    for(Iterator segmentIter = this->Begin(); segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;
        NodeIDType currentBackID = currentSegment.Back();
        typename NodeType::RealType currentBackRadius = GetNode(currentBackID).GetRadius();
        NodeIDType currentFrontID = currentSegment.Front();
        typename NodeType::RealType currentFrontRadius = GetNode(currentFrontID).GetRadius();
        
        if( currentBackRadius > currentFrontRadius )
        {
            currentSegment.Reverse();
        }
        this->segmentTable[currentSegment.GetID()] = currentSegment;
        
        //align tangent vectors
        typename SegmentType::Iterator currentNodeIter = currentSegment.Begin();
        typename SegmentType::Iterator nextNodeIter = currentSegment.Begin();
        nextNodeIter++;
        for(; nextNodeIter != currentSegment.End() && currentNodeIter != currentSegment.End(); ++nextNodeIter, ++currentNodeIter)
        {
            NodeIDType currentNodeID = *currentNodeIter;
            NodeIDType nextNodeID = *nextNodeIter;
            
            NodeType currentNode = this->GetNode(currentNodeID);
            typename NodeType::PositionType currentPosition = currentNode.GetPosition();
            
            NodeType nextNode = this->GetNode(nextNodeID);
            typename NodeType::PositionType nextPosition = nextNode.GetPosition();
            
            typename NodeType::VectorType currentTangent = nextPosition-currentPosition;
            currentTangent = currentTangent.Normalize();
            if(nextNode.GetType() != NodeType::JUNCTION || currentNode.GetType() != NodeType::JUNCTION)
            {
                currentNode.SetTangent(currentTangent);
            }
            else
            {
                currentNode.SetTangent(zeroVector);
            }
            this->nodeTable[currentNode.GetID()] = currentNode;            
        }
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselGraph<VIndexDimension, TCoordinateType>::GetBoundingPoints(PointType& corner1, PointType& corner2) const
{
    try
    {
        corner1 = (this->nodeTable.begin())->second.GetPosition();
        corner2 = corner1;
    
        ConstIterator segmentIter = this->Begin();
        for(; segmentIter != this->End(); ++segmentIter)
        {
            SegmentIDType currentSegmentID = segmentIter->first;
            SegmentType currentSegment = segmentIter->second;
            typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
            for(; nodeIter != currentSegment.End(); ++nodeIter)
            {
                NodeIDType currentNodeID = *nodeIter;
                NodeType currentNode = this->GetNode(currentNodeID);
                PointType currentPosition = currentNode.GetPosition();
                
                for(unsigned int i = 0; i < VIndexDimension; i++)
                {
                    if(currentPosition[i] < corner1[i])
                    {
                        corner1[i] = currentPosition[i];
                    }
                    if(currentPosition[i] > corner2[i])
                    {
                        corner2[i] = currentPosition[i];
                    }
                }
                    
            }
                
                
        }
          
    }
    catch(...)
    {
    
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::ParameterTableType 
VesselGraph<VIndexDimension, TCoordinateType>::GetNodeParameters() const
{   
    ParameterTableType parameterTable;

    //First pass through
    typename NodeTableType::const_iterator nodeIter = nodeTable.begin();

    NodeType currentNode = nodeIter->second;
    typename NodeType::PropertyTable propertyTable = currentNode.GetProperties();
    
    typename NodeType::PropertyTable::const_iterator propertyIter = propertyTable.begin();
    for(; propertyIter != propertyTable.end(); ++propertyIter)
    {
        typename NodeType::PropertyTable::key_type currentKey = propertyIter->first;
        RealType currentValue = atof(propertyIter->second.c_str());
        ParameterRangeType firstRange(currentValue, currentValue);
        parameterTable[currentKey] = firstRange;
    }

    //Go through the rest of the nodes
    for(; nodeIter != nodeTable.end(); ++nodeIter)
    {
        currentNode = nodeIter->second;
        typename NodeType::PropertyTable propertyTable = currentNode.GetProperties();
        
        typename NodeType::PropertyTable::const_iterator propertyIter = propertyTable.begin();
        for(; propertyIter != propertyTable.end(); ++propertyIter)
        {
            typename NodeType::PropertyTable::key_type currentKey = propertyIter->first;
            RealType currentValue = atof(propertyIter->second.c_str());
            ParameterRangeType currentRange = parameterTable[currentKey];
            
            if(currentValue < currentRange.first)
            {
                currentRange.first = currentValue;
            }
            else if(currentValue > currentRange.second)
            {
                currentRange.second = currentValue;
            }
            parameterTable[currentKey] = currentRange;
        }
        
    }
  
    return(parameterTable);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselGraph<VIndexDimension, TCoordinateType>::ParameterTableType 
VesselGraph<VIndexDimension, TCoordinateType>::GetSegmentParameters() const
{
    ParameterTableType parameterTable;

    //First pass through
    typename SegmentTableType::const_iterator segmentIter = segmentTable.begin();

    SegmentType currentSegment = segmentIter->second;
    typename SegmentType::PropertyTable propertyTable = currentSegment.GetProperties();
    
    typename SegmentType::PropertyTable::const_iterator propertyIter = propertyTable.begin();
    for(; propertyIter != propertyTable.end(); ++propertyIter)
    {
        typename SegmentType::PropertyTable::key_type currentKey = propertyIter->first;
        RealType currentValue = atof(propertyIter->second.c_str());
        ParameterRangeType firstRange(currentValue, currentValue);
        parameterTable[currentKey] = firstRange;
    }

    //Go through the rest of the nodes
    for(; segmentIter != segmentTable.end(); ++segmentIter)
    {
        currentSegment = segmentIter->second;
        typename SegmentType::PropertyTable propertyTable = currentSegment.GetProperties();
        
        typename SegmentType::PropertyTable::const_iterator propertyIter = propertyTable.begin();
        for(; propertyIter != propertyTable.end(); ++propertyIter)
        {
            typename SegmentType::PropertyTable::key_type currentKey = propertyIter->first;
            RealType currentValue = atof(propertyIter->second.c_str());
            ParameterRangeType currentRange = parameterTable[currentKey];
            
            if(currentValue < currentRange.first)
            {
                currentRange.first = currentValue;
            }
            else if(currentValue > currentRange.second)
            {
                currentRange.second = currentValue;
            }
            parameterTable[currentKey] = currentRange;
        }
        
    }
    
    return(parameterTable);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& VesselGraph<VIndexDimension, TCoordinateType>::ReportMetricsAsPython(std::ostream &os) const
{
    unsigned int numPoints = this->nodeTable.size();
    unsigned int numSegments = this->NumberOfSegments();
    unsigned int numJunctions = this->NumberOfJunctions();

    //Segment Properties
    std::list<SegmentIDType> segmentIds;
    std::list<unsigned int> numNodes;
    std::list<RealType> lengths;
    std::list<RealType> chordLengths;
    std::list<RealType> distanceMetrics;
    std::list<RealType> soams;
    std::list<RealType> averageRadii;
    std::list<RealType> curvPerLengths;
    std::list<RealType> curvSqPerLengths;
    std::list<RealType> derivCurvPerLengths;
    std::list<RealType> tortDenisty;
    std::list<RealType> inflections;
    std::list<RealType> curvMag;

    ConstIterator segmentIter = this->Begin();
    for(; segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;
        SegmentIDType currentSegmentID = currentSegment.GetID();
        segmentIds.push_back(currentSegmentID);
    
        numNodes.push_back(currentSegment.NumberOfNodes());
        lengths.push_back(currentSegment.GetLength());
        chordLengths.push_back(currentSegment.GetChordLength());
        distanceMetrics.push_back(currentSegment.GetDistanceMetric());
        soams.push_back(currentSegment.GetSumOfAngleMetric());
        averageRadii.push_back(currentSegment.GetAverageRadius());
        curvPerLengths.push_back(currentSegment.GetCurvatureOverLengthMetric());
        curvSqPerLengths.push_back(currentSegment.GetCurvatureSquaredTimesLength());
        derivCurvPerLengths.push_back(currentSegment.GetDerivateOfCurvatureOverLengthMetric());
        tortDenisty.push_back(currentSegment.GetTortuosityDensity());
        inflections.push_back(currentSegment.GetInflectionCountMetric());
        curvMag.push_back(currentSegment.GetSumMagnitudeOfCurvature());
    }
    
    //node properties
    typename NodeTableType::const_iterator nodeIter = this->nodeTable.begin();
    std::list<NodeIDType> nodeIds;
    std::list<RealType> radii;
    std::list<RealType> medialness;
    std::list<RealType> curvatures;
    std::list<RealType> torsions;
    std::list<RealType> inPlaneAngles;
    std::list<RealType> torsionalAngles;
    std::list<RealType> totalAngles;
    
    for(; nodeIter != this->nodeTable.end(); ++nodeIter)
    {
        NodeType currentNode = nodeIter->second;
        nodeIds.push_back(currentNode.GetID());
        radii.push_back(currentNode.GetRadius());
        medialness.push_back(currentNode.GetMedialness());
        curvatures.push_back(currentNode.GetCurvature());
        torsions.push_back(currentNode.GetTorsion());
        inPlaneAngles.push_back(currentNode.GetInPlaneAngle());
        torsionalAngles.push_back(currentNode.GetTorsionalAngle());
        totalAngles.push_back(currentNode.GetTotalAngle());
    }
    
    
    os << "import numpy" << std::endl << std::endl;
    
    ValueToPython(os, "numNodes", numPoints);
    ValueToPython(os, "numSegments", numSegments);
    ValueToPython(os, "numJunctions", numJunctions);
    
    //histogram
    HistogramType histogram = this->GetHistogram();
    HistogramAxisType histX = histogram.first;
    HistogramAxisType histY = histogram.second;
    ArrarToPython(os, "xDegreeHistogram", histX.begin(), histX.end());
    ArrarToPython(os, "yDegreeHistogram", histY.begin(), histY.end());
    
    //segment 
    ArrarToPython(os, "segmentIds", segmentIds.begin(), segmentIds.end());
    ArrarToPython(os, "numNodes", numNodes.begin(), numNodes.end());
    ArrarToPython(os, "lengths", lengths.begin(), lengths.end());
    ArrarToPython(os, "chordLengths", chordLengths.begin(), chordLengths.end());
    ArrarToPython(os, "distanceMetrics", distanceMetrics.begin(), distanceMetrics.end());
    ArrarToPython(os, "soams", soams.begin(), soams.end());
    ArrarToPython(os, "averageRadii", averageRadii.begin(), averageRadii.end());
    ArrarToPython(os, "curvatureMetrics", curvPerLengths.begin(), curvPerLengths.end());
    ArrarToPython(os, "curvatureSquaredMetrics", curvSqPerLengths.begin(), curvSqPerLengths.end());
    ArrarToPython(os, "dervivativeCurvatureMetrics", derivCurvPerLengths.begin(), derivCurvPerLengths.end());
    ArrarToPython(os, "tortuosityDensity", tortDenisty.begin(), tortDenisty.end());
    ArrarToPython(os, "inflectionCount", inflections.begin(), inflections.end());
    ArrarToPython(os, "curvatureMagnitude", curvMag.begin(), curvMag.end());
    
    //node
    ArrarToPython(os, "nodeIds", nodeIds.begin(), nodeIds.end());
    ArrarToPython(os, "radii", radii.begin(), radii.end());
    ArrarToPython(os, "medialness", medialness.begin(), medialness.end());
    ArrarToPython(os, "curvature", curvatures.begin(), curvatures.end());
    ArrarToPython(os, "torsions", torsions.begin(), torsions.end());
    ArrarToPython(os, "inPlaneAngles", inPlaneAngles.begin(), inPlaneAngles.end());
    ArrarToPython(os, "torsionalAngles", torsionalAngles.begin(), torsionalAngles.end());
    ArrarToPython(os, "totalAngles", totalAngles.begin(), totalAngles.end());

}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& VesselGraph<VIndexDimension, TCoordinateType>::Write(std::ostream& os) const
{
    os << std::endl;
    os << segmentTable.size() << std::endl;
    
    ConstIterator segmentIter = this->Begin();
    for(; segmentIter != this->End(); ++segmentIter)
    {
        SegmentType currentSegment = segmentIter->second;
        os << currentSegment << std::endl;
        
        typename SegmentType::ConstIterator nodeIter = currentSegment.Begin();
        for(; nodeIter != currentSegment.End(); ++nodeIter)
        {
            try
            {
                NodeType currentNode = this->GetNode(*nodeIter);
                os << currentNode;
            }
            catch(...)
            {
                std::cout << "Writing graph.  Couldn't find Node " << *nodeIter << " on Segment " << currentSegment.GetID() << std::endl;
            }
        }
    }

    return(os);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& VesselGraph<VIndexDimension, TCoordinateType>::WriteGraphML(std::ostream& os) const
{
    typedef unsigned long GMLIDType;
    typedef std::map<NodeIDType, GMLIDType> NodeIDToGMLIDType;
    NodeIDToGMLIDType nodeIDToGMLIDMap;

    GMLIDType currentNodeID = 1;

    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << std::endl;
    os << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"" << std::endl;  
    os << "    xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" << std::endl;
    os << "    xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns " << std::endl; 
    os << "    http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">" << std::endl;
    os << std::endl;
    
    //node attributes
    os << "    <key id=\"na1\" for=\"node\" attr.name=\"label\" attr.type=\"long\"/>" << std::endl;
    os << "    <key id=\"na2\" for=\"node\" attr.name=\"positionX\" attr.type=\"double\"/>" << std::endl;
    os << "    <key id=\"na3\" for=\"node\" attr.name=\"positionY\" attr.type=\"double\"/>" << std::endl;
    os << "    <key id=\"na4\" for=\"node\" attr.name=\"positionZ\" attr.type=\"double\"/>" << std::endl;
    os << "    <key id=\"na5\" for=\"node\" attr.name=\"radius\" attr.type=\"double\"/>" << std::endl;
    os << "    <key id=\"na6\" for=\"node\" attr.name=\"curvature\" attr.type=\"double\"/>" << std::endl;
    os << "    <key id=\"na7\" for=\"node\" attr.name=\"torsion\" attr.type=\"double\"/>" << std::endl;

    os << std::endl;
    os << "    <graph id=\"G\" edgedefault=\"undirected\">" << std::endl;


    //write nodes
    typename NodeTableType::const_iterator nodeIter = this->nodeTable.begin();
    for(; nodeIter != this->nodeTable.end(); nodeIter++)
    {
        NodeType currentNode = nodeIter->second;
        
        os << "        <node id=\"" << currentNodeID << "\">" << std::endl;
        os << "            <data key=\"na1\">" <<  currentNode.GetID() << "</data>" << std::endl;
        os << "            <data key=\"na2\">" <<  currentNode.GetPosition()[0] << "</data>" << std::endl;
        os << "            <data key=\"na3\">" <<  currentNode.GetPosition()[1] << "</data>" << std::endl;
        os << "            <data key=\"na4\">" <<  currentNode.GetPosition()[2] << "</data>" << std::endl;
        os << "            <data key=\"na5\">" <<  currentNode.GetRadius() << "</data>" << std::endl;
        os << "            <data key=\"na6\">" <<  currentNode.GetCurvature() << "</data>" << std::endl;
        os << "            <data key=\"na7\">" <<  currentNode.GetTorsion() << "</data>" << std::endl;
        os << "        </node>" << std::endl;
        
        nodeIDToGMLIDMap[currentNode.GetID()] = currentNodeID;
        currentNodeID++;
    }
    
    typedef std::pair<GMLIDType, GMLIDType> GMLIDPairType;
    typedef std::list<GMLIDPairType> GMLIDPairListType;
    GMLIDPairListType gmlIDPairList;
    //write edges
    GMLIDType currentEdgeID = 1;
    nodeIter = this->nodeTable.begin();
    for(; nodeIter != this->nodeTable.end(); nodeIter++)
    {
        NodeType currentNode = nodeIter->second;
        GMLIDType currentGMLID = nodeIDToGMLIDMap[currentNode.GetID()];
        typename NodeType::NodeIDListType currentNeighbors = currentNode.GetNeighbors();
        typename NodeType::NodeIDListType::const_iterator neighIter = currentNeighbors.begin();
        for(; neighIter != currentNeighbors.end(); ++neighIter)
        {
            NodeIDType neigborNodeID = *neighIter;
            if(nodeIDToGMLIDMap.find(neigborNodeID) != nodeIDToGMLIDMap.end())
            {
                GMLIDType neigborGMLID = nodeIDToGMLIDMap[neigborNodeID];
                
                GMLIDPairType pairF(currentGMLID, neigborGMLID);
                GMLIDPairType pairB(neigborGMLID, currentGMLID);
                if(std::find(gmlIDPairList.begin(), gmlIDPairList.end(), pairF) == gmlIDPairList.end() &&
                   std::find(gmlIDPairList.begin(), gmlIDPairList.end(), pairB) == gmlIDPairList.end() )
                {   
                    os << "        <edge id=\"" << currentEdgeID << "\" source=\"" << currentGMLID << "\" target=\"" << neigborGMLID <<"\"/>" << std::endl;
                    
                    gmlIDPairList.push_back(pairF);
                    currentEdgeID++;
                }
            }
        }
    }
    
    os << "     </graph>"<< std::endl;
    os << "</graphml>"<< std::endl;
    os << std::endl;
    
    return(os);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& VesselGraph<VIndexDimension, TCoordinateType>::WriteGMLJunctionOnly(std::ostream& os) const
{
    typedef unsigned long GMLIDType;
    typedef std::map<NodeIDType, GMLIDType> NodeIDToGMLIDType;
    NodeIDToGMLIDType nodeIDToGMLIDMap;
    
    typedef std::pair<NodeIDType, NodeIDType> NodeIDPairType;
    typedef std::list<NodeIDPairType> NodeIDPairListType;
    NodeIDPairListType nodeIDPairList;

    GMLIDType currentNodeID = 1;

    os << "graph [" << std::endl;
    
    for(ConstIterator segmentIter = this->Begin(); segmentIter != this->End(); ++segmentIter)
    {
        SegmentIDType currentSegmentID = segmentIter->first;
        SegmentType currentSegment = segmentIter->second;
        
        NodeIDType frontID = currentSegment.Front();
        NodeType frontNode = this->GetNode(frontID);
        
        NodeIDType backID = currentSegment.Back();
        NodeType backNode = this->GetNode(backID);
        
        if(nodeIDToGMLIDMap.find(frontID) == nodeIDToGMLIDMap.end())
        {
            os << "  node [" << std::endl;
            os << "    id " << currentNodeID << std::endl;
            os << "    label \"" << frontNode.GetID() <<"\"" << std::endl;
            os << "    position [" << std::endl;
            typename NodeType::PositionType frontPosition = frontNode.GetPosition();
            for(unsigned int i = 0; i <  VIndexDimension; i++)
            {
                os << "      x" << i << " " << frontPosition[i] << std::endl;
            }
            os << "    ]" << std::endl;
            os << "  ]" << std::endl;
            
            nodeIDToGMLIDMap[frontID] =currentNodeID;
            currentNodeID++;
        }
        if(nodeIDToGMLIDMap.find(backID) == nodeIDToGMLIDMap.end())
        {                                                              
            os << "  node [" << std::endl;
            os << "    id " << currentNodeID << std::endl;
            os << "    label \"" << backNode.GetID() <<"\"" << std::endl;
            os << "    position [" << std::endl;
            typename NodeType::PositionType backPosition = backNode.GetPosition();
            for(unsigned int i = 0; i <  VIndexDimension; i++)
            {
                os << "      x" << i << " " << backPosition[i] << std::endl;
            }
            os << "    ]" << std::endl;
            os << "  ]" << std::endl;
            
            nodeIDToGMLIDMap[backID] = currentNodeID;
            currentNodeID++;
        }
        
        NodeIDPairType currentEdge(frontID, backID);
        if( std::find(nodeIDPairList.begin(), nodeIDPairList.end(), currentEdge) == nodeIDPairList.end() && frontID != backID)
        {
            os << "  edge [" << std::endl;
            os << "    source " << nodeIDToGMLIDMap[frontID] << std::endl;
            os << "    target " << nodeIDToGMLIDMap[backID] << std::endl;
            os << "    weight " << frontNode.Distance(backNode) << std::endl;
            os << "  ]" << std::endl;
            nodeIDPairList.push_back(currentEdge);    
        }

    }
    
    os << "]"<< std::endl;
    os << std::endl;
    
    return(os);
}
template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& VesselGraph<VIndexDimension, TCoordinateType>::WriteGML(std::ostream& os) const
{
    typedef unsigned long GMLIDType;
    typedef std::map<NodeIDType, GMLIDType> NodeIDToGMLIDType;
    NodeIDToGMLIDType nodeIDToGMLIDMap;

    GMLIDType currentNodeID = 1;

    os << "graph [" << std::endl;

    //write nodes
    typename NodeTableType::const_iterator nodeIter = this->nodeTable.begin();
    for(; nodeIter != this->nodeTable.end(); nodeIter++)
    {
        NodeType currentNode = nodeIter->second;
        os << "  node [" << std::endl;
        os << "    id " << currentNodeID << std::endl;
        os << "    label \"" << currentNode.GetID() <<"\"" << std::endl;
        os << "    position [" << std::endl;
        typename NodeType::PositionType currentPosition = currentNode.GetPosition();
        for(unsigned int i = 0; i <  VIndexDimension; i++)
        {
            os << "      x" << i << " " << currentPosition[i] << std::endl;
        }
        os << "    ]" << std::endl;
        os << "    radius " << currentNode.GetRadius() << std::endl;
        os << "    curvature " << currentNode.GetCurvature() << std::endl;
        os << "    torsion " << currentNode.GetTorsion() << std::endl;
        os << "  ]" << std::endl;
        
        nodeIDToGMLIDMap[currentNode.GetID()] = currentNodeID;
        currentNodeID++;
    }
    
    typedef std::pair<GMLIDType, GMLIDType> GMLIDPairType;
    typedef std::list<GMLIDPairType> GMLIDPairListType;
    GMLIDPairListType gmlIDPairList;
    //write edges
    nodeIter = this->nodeTable.begin();
    for(; nodeIter != this->nodeTable.end(); nodeIter++)
    {
        NodeType currentNode = nodeIter->second;
        GMLIDType currentGMLID = nodeIDToGMLIDMap[currentNode.GetID()];
        typename NodeType::NodeIDListType currentNeighbors = currentNode.GetNeighbors();
        typename NodeType::NodeIDListType::const_iterator neighIter = currentNeighbors.begin();
        for(; neighIter != currentNeighbors.end(); ++neighIter)
        {
            NodeIDType neigborNodeID = *neighIter;
            if(nodeIDToGMLIDMap.find(neigborNodeID) != nodeIDToGMLIDMap.end())
            {
                GMLIDType neigborGMLID = nodeIDToGMLIDMap[neigborNodeID];
                
                GMLIDPairType pairF(currentGMLID, neigborGMLID);
                GMLIDPairType pairB(neigborGMLID, currentGMLID);
                if(std::find(gmlIDPairList.begin(), gmlIDPairList.end(), pairF) == gmlIDPairList.end() &&
                   std::find(gmlIDPairList.begin(), gmlIDPairList.end(), pairB) == gmlIDPairList.end() )
                {
                    os << "  edge [" << std::endl;
                    os << "    source " << currentGMLID << std::endl;
                    os << "    target " << neigborGMLID << std::endl;
                    os << "  ]" << std::endl;
                    gmlIDPairList.push_back(pairF);
                }
            }
            
        }
        
    }
    
    os << "]"<< std::endl;
    os << std::endl;
    
    return(os);
}


template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream& VesselGraph<VIndexDimension, TCoordinateType>::Read(std::istream& is) 
{
    unsigned int numSegments;
   
    is >> numSegments;
    for(unsigned int i = 0; i < numSegments; i++)
    {
        SegmentType s;
        is >> s;
        this->AddSegment(s);
        
        unsigned int numNodes = s.NumberOfNodes();
        
        for(unsigned int j = 0; j < numNodes; j++)
        {
            NodeType n;
            is >> n;
            this->AddNode(n);
        }   
    }
    
    this->ComputeSegmentParameters();
    return(is);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream & operator<<(std::ostream &os, const VesselGraph<VIndexDimension, TCoordinateType> &vg)
{
    return(vg.Write(os));
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream & operator>>(std::istream &is, VesselGraph<VIndexDimension, TCoordinateType> &vg)
{
    return(vg.Read(is));
}

template<typename T>
std::ostream& ValueToPython(std::ostream &os, std::string name, T param)
{
    os << name << " = " << param << std::endl;
    return(os);
}

template<typename RandomAccessIterator>
std::ostream& ArrarToPython(std::ostream &os, std::string name, RandomAccessIterator begin, RandomAccessIterator end)
{
    os << name << " = numpy.array([";
    for(; begin != end; ++begin)
    {
        os << *begin << ", ";
    }
    os << "])" << std::endl;
    return(os);
}

