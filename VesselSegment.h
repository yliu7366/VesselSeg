#ifndef __VesselSegment_h
#define __VesselSegment_h

#include <list>
#include <iostream>
#include <map>
#include <string>

#include "VesselNode.h"

/**\class VesselSegment
 * \brief A simple container holding the properties of a segment in a vessel network
 *
 * This class keeps track of the nodes belonging to the segment and 
 * morphological properties of the segment
 */
template<unsigned int VIndexDimension, typename TCoordinateType>  
class VesselSegment
{
    public:
        typedef double RealType;
    
        typedef VesselNode< VIndexDimension, TCoordinateType> NodeType;
        typedef typename NodeType::IDType NodeIDType;
        typedef std::list<NodeIDType> NodeIDListType; 
        typedef typename NodeIDListType::iterator Iterator; 
        typedef typename NodeIDListType::const_iterator ConstIterator; 
        typedef NodeIDType IDType;
        typedef std::map<std::string, std::string> PropertyTable;
        
        VesselSegment();
        virtual ~VesselSegment();
        
        /**
        * Set segment id
        * @param id segment id
        */
        void SetID(IDType id);
        
        /**
        * Get segment id
        */
        IDType GetID() const;
        
        /**
        * Add node
        * @param id node id to add
        */
        void AddNode(NodeIDType id);
        
        /**
        * Replace node
        * @param id1 node id to remove
        * @param id2 node id to add
        */
        void ReplaceNode(NodeIDType id1, NodeIDType id2);
        
        /**
        * Remove node
        */
        void RemoveNode(NodeIDType);
        
        /**
        * Remove all nodes
        */
        void RemoveAllNodes();
        
        /**
        * Get number of nodes
        */
        unsigned int NumberOfNodes() const;
    
        /**
        * Get first node in segment 
        */
        NodeIDType Front() const;
        
        /**
        * Get last node in segment 
        */
        NodeIDType Back() const;
        
        /**
        * Get beginning of node collection
        * @return node collection iterator
        */
        Iterator Begin();
        
        /**
        * Get end of node collection
        * @return node collection iterator
        */
        Iterator End();
        
        /**
        * Get beginning of node collection
        * @return node collection iterator
        */
        ConstIterator Begin() const;
        
        /**
        * Get end of node collection
        * @return node collection iterator
        */
        ConstIterator End() const;
        
        /**
        * Reverse order of nodes
        */
        void Reverse();
       
        /**
        * Write segment in .vessel format
        */
        std::ostream& Write(std::ostream& os) const;
        
        /**
        * Read segment in .vessel format
        */
        std::istream& Read(std::istream& is);
        
        RealType GetLength() const;
        void SetLength(RealType);
        
        RealType GetChordLength() const;
        void SetChordLength(RealType);
        
        RealType GetAverageRadius() const;
        void SetAverageRadius(RealType);
        
        RealType GetSumOfAngleMetric() const;
        void SetSumOfAngleMetric(RealType);
        
        RealType GetSumMagnitudeOfCurvature() const;
        void SetSumMagnitudeOfCurvature(RealType);
        
        RealType GetSumMagnitudeOfCurvatureSquared() const;
        void SetSumMagnitudeOfCurvatureSquared(RealType);
        
        RealType GetSumMagnitudeOfDervativeCurvature() const;
        void SetSumMagnitudeOfDervativeCurvature(RealType);
        
        RealType GetTortuosityDensity() const;
        void SetTortuosityDensity(RealType);
        
        RealType GetInflectionCountMetric() const;
        void SetInflectionCountMetric(RealType);
        

        RealType GetDistanceMetric() const;
        RealType GetCurvatureOverLengthMetric() const;
        RealType GetCurvatureSquaredTimesLength() const;
        RealType GetDerivateOfCurvatureOverLengthMetric() const;
        
        PropertyTable GetProperties() const;
        
    protected:
        NodeIDListType nodeList;
        IDType id;
        
        RealType length;
        RealType chordLength;
        RealType averageRadius;
        RealType soam;
        RealType sumMagnitudeOfCurvature;
        RealType sumMagnitudeOfCurvatureSquared;
        RealType sumMagnitudeOfDervativeCurvature;
        RealType tortuosityDensity;
        RealType inflectionPoint;
    
};

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream & operator<<(std::ostream &os, const VesselSegment<VIndexDimension, TCoordinateType> &vs);

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream & operator>>(std::istream &os, VesselSegment<VIndexDimension, TCoordinateType> &vs);

#include "VesselSegment.txx"

#endif
