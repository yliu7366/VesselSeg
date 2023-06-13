#include <iterator>
#include <algorithm>
#include <sstream> 

#include "VesselSegment.h"

template<unsigned int VIndexDimension, typename TCoordinateType>
VesselSegment<VIndexDimension, TCoordinateType>::VesselSegment()
{
    this->id=0;
    this->length = 0.0;
    this->chordLength = 0.0;
    this->averageRadius = 0.0;
    this->soam = 0.0;
    this->sumMagnitudeOfCurvature = 0.0;
    this->sumMagnitudeOfCurvatureSquared = 0.0;
    this->sumMagnitudeOfDervativeCurvature = 0.0;
    this->tortuosityDensity = 0.0;
    this->inflectionPoint = 0.0;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
VesselSegment<VIndexDimension, TCoordinateType>::~VesselSegment()
{
     
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetID(IDType segmentID)
{
    this->id = segmentID;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::IDType 
VesselSegment<VIndexDimension, TCoordinateType>::GetID() const
{
    return(this->id);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::AddNode(NodeIDType nid)
{
    this->nodeList.push_back(nid);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
unsigned int VesselSegment<VIndexDimension, TCoordinateType>::NumberOfNodes() const
{
    return(this->nodeList.size());
}
template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::RemoveNode(NodeIDType nid)
{
    this->nodeList.remove(nid);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::RemoveAllNodes()
{
    this->nodeList.clear();
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::ReplaceNode(NodeIDType old, NodeIDType neww)
{
    typename NodeIDListType::iterator findIter = std::find(this->nodeList.begin(), this->nodeList.end(), old);
    if(findIter != this->nodeList.end())
    {
        this->nodeList.insert(findIter, neww);
    }
    this->nodeList.remove(old);   
}


template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::NodeIDType 
VesselSegment<VIndexDimension, TCoordinateType>::Front() const
{
    return(this->nodeList.front());
}
        
template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::NodeIDType 
VesselSegment<VIndexDimension, TCoordinateType>::Back() const
{
    return(this->nodeList.back());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::Iterator 
VesselSegment<VIndexDimension, TCoordinateType>::Begin()  
{
    return(this->nodeList.begin());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::Iterator 
VesselSegment<VIndexDimension, TCoordinateType>::End()  
{
    return(this->nodeList.end());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::ConstIterator 
VesselSegment<VIndexDimension, TCoordinateType>::Begin() const
{
    return(this->nodeList.begin());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::ConstIterator 
VesselSegment<VIndexDimension, TCoordinateType>::End() const
{
    return(this->nodeList.end());
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::Reverse()
{
    this->nodeList.reverse();
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& VesselSegment<VIndexDimension, TCoordinateType>::Write(std::ostream& os) const
{
    os << std::endl;
    os << this->id << std::endl;
    os << this->nodeList.size() << std::endl;
    std::ostream_iterator<NodeIDType> osIter(os, " ");
    std::copy(this->nodeList.begin(), this->nodeList.end(), osIter);
    os << std::endl;
    
    return(os);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream& VesselSegment<VIndexDimension, TCoordinateType>::Read(std::istream& is) 
{
    is >> this->id;
    unsigned int numNodes;
    is >> numNodes;
    for(unsigned int i = 0; i < numNodes; i++)
    {
        NodeIDType n;
        is >> n;
        this->AddNode(n);
    }
    return(is);
}


template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetLength() const
{
    return(this->length);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetLength(RealType l)
{
    this->length = l;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetChordLength() const
{
    return(this->chordLength);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetChordLength(RealType cl)
{
    this->chordLength = cl;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetAverageRadius() const
{
    return(this->averageRadius);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetAverageRadius(RealType r)
{
    this->averageRadius = r;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetSumOfAngleMetric() const
{
    return(this->soam);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetSumOfAngleMetric(RealType s)
{
    this->soam = s;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetSumMagnitudeOfCurvature() const
{
    return(this->sumMagnitudeOfCurvature);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetSumMagnitudeOfCurvature(RealType c)
{
    this->sumMagnitudeOfCurvature = c;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetSumMagnitudeOfCurvatureSquared() const
{
    return(this->sumMagnitudeOfCurvatureSquared);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetSumMagnitudeOfCurvatureSquared(RealType c)
{
    this->sumMagnitudeOfCurvatureSquared = c;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetSumMagnitudeOfDervativeCurvature() const
{
    return(this->sumMagnitudeOfDervativeCurvature);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetSumMagnitudeOfDervativeCurvature(RealType c)
{
    this->sumMagnitudeOfDervativeCurvature = c;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetTortuosityDensity() const
{
    return(this->tortuosityDensity);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetTortuosityDensity(RealType t)
{
    this->tortuosityDensity = t;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetInflectionCountMetric() const
{
    return(this->inflectionPoint);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
void VesselSegment<VIndexDimension, TCoordinateType>::SetInflectionCountMetric(RealType i)
{
    this->inflectionPoint = i;
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType VesselSegment<VIndexDimension, TCoordinateType>::GetDistanceMetric() const
{
    RealType dm = 0.0;
    if(this->chordLength != 0.0 && this->length != 0.0)
    {
        dm = this->chordLength / this->length;
    }
    
    return(dm);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType 
VesselSegment<VIndexDimension, TCoordinateType>::GetCurvatureOverLengthMetric() const
{
    RealType colm = 0;
    if(this->length != 0.0)
    {
        colm = (this->sumMagnitudeOfCurvature * this->sumMagnitudeOfCurvature) / this->length;
    }
    return(colm);
}


template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::PropertyTable 
VesselSegment<VIndexDimension, TCoordinateType>::GetProperties() const
{
    PropertyTable propertyTable;
    std::stringstream ss;
    
    ss << this->NumberOfNodes();
    propertyTable["Segment Number of Nodes"] = ss.str();
    ss.str("");
    
    ss << this->GetLength();
    propertyTable["Segment Length"] = ss.str();
    ss.str("");
    
    ss << this->GetChordLength();
    propertyTable["Segment Chord Length"] = ss.str();
    ss.str("");
    
    ss << this->GetAverageRadius();
    propertyTable["Segment Average Radius"] = ss.str();
    ss.str("");
    
    ss << this->GetSumOfAngleMetric();
    propertyTable["Segment SOAM"] = ss.str();
    ss.str("");
    
    ss << this->GetSumMagnitudeOfCurvature();
    propertyTable["Segment Curvature"] = ss.str();
    ss.str("");
    
    ss << this->GetSumMagnitudeOfCurvatureSquared();
    propertyTable["Segment Curvature Squared"] = ss.str();
    ss.str("");
    
    ss << this->GetSumMagnitudeOfDervativeCurvature();
    propertyTable["Segment Derivature Curvature"] = ss.str();
    ss.str("");
    
    ss << this->GetTortuosityDensity();
    propertyTable["Segment Tortuosity Density"] = ss.str();
    ss.str("");
    
    ss << this->GetInflectionCountMetric();
    propertyTable["Segment Inflection Count"] = ss.str();
    ss.str("");
    
    ss << this->GetDistanceMetric();
    propertyTable["Segment Distance Metric"] = ss.str();
    ss.str("");
    
    ss << this->GetCurvatureOverLengthMetric();
    propertyTable["Segment Curvature Over Length"] = ss.str();
    ss.str("");
    
    ss << this->GetCurvatureSquaredTimesLength(); 
    propertyTable["Segment Curvature Squared Times Length"] = ss.str();
    ss.str("");
    
    ss << this->GetDerivateOfCurvatureOverLengthMetric();
    propertyTable["Segment Curvature Derivative Over Length"] = ss.str();
    ss.str("");
    
    
    return(propertyTable);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType 
VesselSegment<VIndexDimension, TCoordinateType>::GetCurvatureSquaredTimesLength() const
{
    RealType csolm = 0;
    if(this->length != 0.0)
    {
        csolm = (this->sumMagnitudeOfCurvature * this->sumMagnitudeOfCurvatureSquared) * RealType(this->GetChordLength());
    }
    csolm = sqrt(csolm);
    return(csolm);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
typename VesselSegment<VIndexDimension, TCoordinateType>::RealType 
VesselSegment<VIndexDimension, TCoordinateType>::GetDerivateOfCurvatureOverLengthMetric() const
{
    RealType dcolm = 0.0;
    if(this->length != 0.0)
    {    
        dcolm = (this->sumMagnitudeOfDervativeCurvature* this->sumMagnitudeOfDervativeCurvature) / this->length;
    }
    return(dcolm);
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream & operator<<(std::ostream &os, const VesselSegment<VIndexDimension, TCoordinateType> &vs)
{
    return(vs.Write(os));
}

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream & operator>>(std::istream &is, VesselSegment<VIndexDimension,TCoordinateType> &vs)
{
    return(vs.Read(is));
}

