#ifndef __VesselNode_h
#define __VesselNode_h

#include <list>
#include <iostream>
#include <map>
#include <string>

#include "FixedArray.h"
//#include "VesselSegment.h"

/**\class VesselNode
 * \brief A simple container holding the properties of a node in a vessel network
 * 
 * This class keeps track of node attributes like position, radius, etc ...
 * , it neighboring nodes and which vessel segments it belongs to 
 */ 
template<unsigned int VIndexDimension, typename TCoordinateType>
class VesselNode
{
    public:
        typedef unsigned long IDType;
        typedef FixedArray< TCoordinateType, VIndexDimension > PositionType;
        typedef std::list< IDType > NodeIDListType;
        //typedef typename VesselSegment<VIndexDimension, TCoordinateType>::SegmentID SegmentIDType;
        typedef IDType SegmentIDType;
        typedef std::list< SegmentIDType > SegmentIDListType;
        typedef double RealType;
        typedef FixedArray< RealType, VIndexDimension > VectorType;
        typedef std::map<std::string, std::string> PropertyTable;


        VesselNode();
        virtual ~VesselNode();
        
        /**
        * Set node id
        * @param id node id
        */
        void SetID(IDType id);
        
        /**
        * Get node id
        */
        IDType GetID() const;
        
        /**
        * Node Type
        * TERMINAL has 1 neighbor
        * JUNCTION has 3 or more neighbors
        * MIDSEGMENT has 2 neighbors
        * SINGLE has no neighbors
        */
        enum NodeType{TERMINAL, JUNCTION, MIDSEGMENT, SINGLE};
        
        /**
        * Get node type
        */
        NodeType GetType() const;
        
        /**
        * Get nodes position
        */
        PositionType GetPosition() const;
        
        /**
        * Set nodes position in phyical coordinates
        */
        void SetPosition(PositionType );
        
        /**
        * Add node neighbor
        */
        void AddNeighbor(IDType );
        
        /**
        * Get number of neighbors
        */
        unsigned int GetNumberOfNeighbors() const;
        
        /**
        * Get node neighbors
        */
        NodeIDListType GetNeighbors() const;
        
        /**
        * Remove node neighbor
        */
        void RemoveNeighbor(IDType);
        
        /**
        * Remove node neighbor if its really a neighbor
        */
        void RemoveNeighborIfFound(IDType);
        
        /**
        * Add an assoication to a vessel segment
        */
        void AddSegment(SegmentIDType);
        
        /**
        * Get all segments associated with this node
        */
        SegmentIDListType GetSegments() const;
        
        /**
        * Remove assocation with this node and specified segment
        */
        void RemoveSegment(SegmentIDType);
        
        
        RealType GetRadius() const;
        void SetRadius(RealType );
        
        RealType GetMedialness() const;
        void SetMedialness(RealType );
        
        VectorType GetTangent() const;
        void SetTangent(VectorType );
        
        VectorType GetNormal() const;
        void SetNormal(VectorType );
        
        VectorType GetBiNormal() const;
        void SetBiNormal(VectorType );
        
        RealType GetCurvature() const;
        void SetCurvature(RealType );
        
        RealType GetTorsion() const;
        void SetTorsion(RealType );
        
        RealType GetInPlaneAngle() const;
        void SetInPlaneAngle(RealType );
        
        RealType GetTorsionalAngle() const;
        void SetTorsionalAngle(RealType );
        
        RealType GetTotalAngle() const;      
        
        RealType Distance(const VesselNode& ) const;
        
        PropertyTable GetProperties() const;
        
        std::ostream& Write(std::ostream& ) const;
        std::istream& Read(std::istream& );
        
    protected:
        IDType id;
        PositionType position;
        SegmentIDListType segments;
        NodeIDListType neighbors;
        RealType radius;
        RealType medialness;
        VectorType tangent;
        VectorType normal;
        VectorType biNormal;
        RealType curvature;
        RealType torsion;
        RealType inPlaneAngle;
        RealType torsionalAngle;
};

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream & operator<<(std::ostream &os, const VesselNode<VIndexDimension, TCoordinateType> &vn);

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream & operator>>(std::istream &is, VesselNode<VIndexDimension, TCoordinateType> &vn);

#include "VesselNode.txx"

#endif
