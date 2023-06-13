#include <iterator>
#include <algorithm>
#include <limits>
#include <sstream> 

#include "VesselNode.h"

template<unsigned int VIndexDimension, typename TCoordinateType>
VesselNode<VIndexDimension, TCoordinateType>::VesselNode()
{
    this->id = 0;
    this->position.Fill(0.0);
    this->tangent.Fill(0.0);
    this->normal.Fill(0.0);
    this->biNormal.Fill(0.0);
    this->radius = 0.0;
    this->curvature = 0.0;
    this->torsion = 0.0;
    this->inPlaneAngle = 0.0;
    this->torsionalAngle = 0.0;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
VesselNode<VIndexDimension, TCoordinateType>::~VesselNode()
{

}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetID(IDType anID)
{
    this->id = anID;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::IDType 
VesselNode<VIndexDimension, TCoordinateType>::GetID() const
{
    return(this->id);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::NodeType 
VesselNode<VIndexDimension, TCoordinateType>::GetType() const
{
    NodeType nodeType;
    if(this->neighbors.size() == 0 )
    {
        nodeType = SINGLE;
        return(nodeType);
    }
    else if(this->neighbors.size() == 1 )
    {
        nodeType = TERMINAL;
        return(nodeType);
    }
    else if(this->neighbors.size() == 2 )
    {
        nodeType = MIDSEGMENT;
        return(nodeType);
    }
    else
    {
        nodeType = JUNCTION;
        return(nodeType);
    }
}


template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::PositionType 
VesselNode<VIndexDimension, TCoordinateType>::GetPosition() const
{
    return(this->position);
}


template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetPosition(PositionType p) 
{
    this->position = p;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::AddNeighbor(IDType id)
{
    if(std::find(this->neighbors.begin(), this->neighbors.end(), id) == this->neighbors.end())
    {
        this->neighbors.push_back(id);
    }
}


template<unsigned int VIndexDimension, typename TCoordinateType>
unsigned int VesselNode<VIndexDimension, TCoordinateType>::GetNumberOfNeighbors() const
{
    return(this->neighbors.size());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::NodeIDListType 
VesselNode<VIndexDimension, TCoordinateType>::GetNeighbors() const
{
    return(this->neighbors);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::RemoveNeighbor(IDType id)
{
    this->neighbors.remove(id);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::RemoveNeighborIfFound(IDType id)
{
    typename NodeIDListType::iterator foundIter = std::find(this->neighbors.begin(), this->neighbors.end(), id);
    if(foundIter != this->neighbors.end())
    {
        this->neighbors.erase(foundIter);
    }
}


template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::AddSegment(SegmentIDType segmentID)
{
    if(std::find(this->segments.begin(), this->segments.end(), segmentID) == this->segments.end())
    {
        this->segments.push_back(segmentID);
    }
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::SegmentIDListType 
VesselNode<VIndexDimension, TCoordinateType>::GetSegments() const
{
    return(this->segments);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::RemoveSegment(SegmentIDType id)
{
    this->segments.remove(id);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::GetRadius() const
{
    return(this->radius);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetRadius(RealType r) 
{
    this->radius = r;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::GetMedialness() const
{
    return(this->medialness);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetMedialness(RealType m) 
{
    this->medialness = m;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::VectorType 
VesselNode<VIndexDimension, TCoordinateType>::GetTangent() const
{
    return(this->tangent);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetTangent(VectorType v)
{
    this->tangent = v;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::VectorType 
VesselNode<VIndexDimension, TCoordinateType>::GetNormal() const
{
    return(this->normal);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetNormal(VectorType v)
{
    this->normal = v;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::VectorType 
VesselNode<VIndexDimension, TCoordinateType>::GetBiNormal() const
{
    return(this->biNormal);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetBiNormal(VectorType v)
{
    this->biNormal = v;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::GetCurvature() const
{
    return(this->curvature);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetCurvature(RealType c) 
{
    this->curvature = c;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::GetTorsion() const
{
    return(this->torsion);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetTorsion(RealType t) 
{
    this->torsion = t;
}  

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::GetInPlaneAngle() const
{
    return(this->inPlaneAngle);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetInPlaneAngle(RealType a) 
{
    this->inPlaneAngle = a;
}  

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::GetTorsionalAngle() const
{
    return(this->torsionalAngle);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselNode<VIndexDimension, TCoordinateType>::SetTorsionalAngle(RealType a) 
{
    this->torsionalAngle = a;
} 

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::GetTotalAngle() const
{
    return(this->inPlaneAngle + this->torsionalAngle);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::RealType 
VesselNode<VIndexDimension, TCoordinateType>::Distance(const VesselNode& vn) const
{
    PositionType i1 = this->GetPosition();
    PositionType i2 = vn.GetPosition();
    
    RealType distance = i1.Distance(i2);
    
    return(distance);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselNode<VIndexDimension, TCoordinateType>::PropertyTable 
VesselNode<VIndexDimension, TCoordinateType>::GetProperties() const
{
    PropertyTable propertyTable;
    std::ostringstream ss;
    ss << this->GetRadius();
    propertyTable["Node Radius"] = ss.str();
    ss.str("");

    ss << this->GetMedialness();
    propertyTable["Node Medialness"] = ss.str();
    ss.str("");
    
    ss << this->GetCurvature();
    propertyTable["Node Curvature"] = ss.str();
    ss.str("");
    
    ss << this->GetTorsion();
    propertyTable["Node Torsion"] = ss.str();
    ss.str("");
    
    ss << this->GetInPlaneAngle();
    propertyTable["Node In Plane Angle"] = ss.str();
    ss.str("");
    
    ss << this->GetTorsionalAngle();
    propertyTable["Node Torsional Angle"] = ss.str();
    ss.str("");
    
    ss << this->GetTotalAngle();
    propertyTable["Node Total Angle"] = ss.str();
    ss.str("");
/*    
    typename PropertyTable::const_iterator iter = propertyTable.begin();
    for(; iter != propertyTable.end(); ++iter)
    {
        std::cout << "Key " << iter->first << " Value " << iter->second << std::endl;
    }
*/    
    return(propertyTable);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& VesselNode<VIndexDimension, TCoordinateType>::Write(std::ostream& os) const
{
    //os.precision(std::numeric_limits< RealType >::digits10+5);
    os.precision(7);
    os << std::endl;
    os << this->id << std::endl;
    os << this->position << std::endl;
    
    os << this->neighbors.size() << std::endl;
    std::ostream_iterator<IDType> osIterNeigh(os, " ");
    std::copy(this->neighbors.begin(), this->neighbors.end(), osIterNeigh);
    os << std::endl;
    
    os << this->segments.size() << std::endl;
    std::ostream_iterator<SegmentIDType> osIterSeg(os, " ");
    std::copy(this->segments.begin(), this->segments.end(), osIterSeg);
    os << std::endl;

    os << this->radius << std::endl;
    os << this->medialness << std::endl;
    os << this->tangent << std::endl;
    os << this->normal << std::endl;
    os << this->biNormal << std::endl;
    os << this->curvature << std::endl;
    os << this->torsion << std::endl;
    os << this->inPlaneAngle << std::endl;
    os << this->torsionalAngle << std::endl;
    
    return(os);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream& VesselNode<VIndexDimension, TCoordinateType>::Read(std::istream& is) 
{
    //is.precision(std::numeric_limits< RealType >::digits10+5);
    is.precision(7);
    is >> this->id;
    is >> this->position;

    unsigned int numNeighbors;  
    is >> numNeighbors;
    for(unsigned int i = 0; i < numNeighbors; i++)
    {        
        IDType currentID;
        is >> currentID;
        
        this->AddNeighbor(currentID);
    }
    
    unsigned int numSegments;  
    is >> numSegments;
    for(unsigned int i = 0; i < numSegments; i++)
    {        
        SegmentIDType currentID;
        is >> currentID;
        
        this->AddSegment(currentID);
    }

    is >> this->radius;
    is >> this->medialness;
    is >> this->tangent;
    is >> this->normal;
    is >> this->biNormal;
    is >> this->curvature;
    is >> this->torsion;
    is >> this->inPlaneAngle;
    is >> this->torsionalAngle;
    
    return(is);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream & operator<<(std::ostream &os, const VesselNode<VIndexDimension, TCoordinateType> &vn)
{
    return(vn.Write(os));
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream & operator>>(std::istream &is, VesselNode<VIndexDimension, TCoordinateType> &vn)
{
    return(vn.Read(is));
}


