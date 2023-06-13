#ifndef __VesselGraph_h
#define __VesselGraph_h

#include <list>
#include <map>
#include <string>
#include <utility>
#include <iterator>
#include <iostream>
#include <vector>

#include "VesselNode.h"
#include "VesselSegment.h"

/**\class VesselGraph
 * \brief A class which represents a vessel graph
 *
 * To define the graph we have a node table and a segment table.
 * The node table simply maps node ids to node objects.  We define
 * a segment as a collection of nodes where the beginning and end of the segment
 * contain either a end or junction node.  An end node has only one
 * node neighboor while a junction node has three or more.  The segment table
 * maps a segment id to a collection of node ids.
 * This class also provides some functionality for graph topology and 
 * morphology
 */
template<unsigned int VIndexDimension, typename TCoordinateType>  
class VesselGraph
{
    public:
        typedef VesselGraph<VIndexDimension, TCoordinateType> Self;
        typedef double RealType;
    
        typedef VesselNode<VIndexDimension, TCoordinateType> NodeType;
        typedef typename NodeType::IDType NodeIDType;
        typedef std::map<NodeIDType, NodeType> NodeTableType;
        typedef typename NodeType::PositionType PointType;
        
        typedef VesselSegment<VIndexDimension, TCoordinateType> SegmentType;
        typedef typename SegmentType::IDType SegmentIDType;
        typedef std::map<SegmentIDType, SegmentType> SegmentTableType;
        typedef typename SegmentTableType::iterator Iterator;
        typedef typename SegmentTableType::const_iterator ConstIterator;
        
        typedef std::pair<RealType, RealType> ParameterRangeType;
        typedef std::map<std::string, ParameterRangeType> ParameterTableType;
        
        typedef std::vector<int> HistogramAxisType;
        typedef std::pair<HistogramAxisType, HistogramAxisType> HistogramType;
        
       
        VesselGraph();
        virtual ~VesselGraph();
        
        /**
        * Remove all nodes and segments from the graph
        */
        void Clear();

        /**
        * Add a node to the graph
        * @param n the node
        */
        void AddNode(NodeType n);
        
        /**
        * Get a node from the graph
        * @param nid  node id
        * @return the node
        */
        NodeType GetNode(const NodeIDType& nid) const;
        
        /**
        * Add a segment to the graph
        * @param s the node
        */
        void AddSegment(SegmentType s);
        
        /**
        * Get a segment from the graph
        * @param sid  segment id
        * @return the segment
        */
        SegmentType GetSegment(const SegmentIDType& sid) const;
        
        /**
        * Remove all segments (and associated nodes) whose number of nodes
        * is less than the specified and that are not connected to any other
        * segments
        * @param n number of nodes 
        */
        void RemoveUnconnectedSegments(unsigned int n);
        
        /**
        * Query a segment to determine if it is connected to any other segments
        * @param sid  segment id
        * @return true is segment is connected to another segment
        */
        bool IsSegmentConnected(const SegmentIDType& sid) const;
        
        /**
        * Remove segment from graph and make appropriate changes to nodes
        * @param sid  segment id
        */
        void RemoveSegment(const SegmentIDType& sid);
        
        /**
        * Join two segments together and make appropriate changes to nodes
        * @param sid1  segment id
        * @param sid2  segment id
        */
        void JoinSegments(const SegmentIDType& sid1, const SegmentIDType& sid2);
        
        /**
        * Remove all segments connected to specified segment
        * @param sid  segment id
        * @return the new graph
        */
        Self RemoveSubGraph(const SegmentIDType& sid) const;
        
        /**
        * Get the graph containing all segments connected to specified segment
        * @param sid  segment id
        * @return the new graph
        */
        Self GetSubGraph(const SegmentIDType& sid) const;
        
        /**
        * Get the graph containing all segments withing the specified radius
        * @param sid  segment id
        * @param r radius
        * @return the new graph
        */
        Self GetGraphInSphere(const NodeIDType& sid, RealType r) const;
        
        /**
        * Remove nodes in graph if there curvature is below the specified
        * @param c curvature
        * @return the new graph
        */
        Self GetDownSampledByCurvature(RealType c) const;
        
        /**
        * Interpolate all graph segments
        * @param os oversampling factor.  Determines the amount of new nodes
        * @return the new graph
        */
        Self Interpolate(RealType os) const;

        /**
        * Remove all redundant junction points ie junction points that 
        * are neighbors to junction points because of discretization errors
        * @return the new graph
        */
        Self RemoveRedundentJunctions() const;

        /**
        * Get beginning of segment table
        * @return segment table iterator
        */
        Iterator Begin();
        
        /**
        * Get end of segment table
        * @return segment table iterator
        */
        Iterator End();
        
        /**
        * Get beginning of segment table
        * @return segment table iterator
        */
        ConstIterator Begin() const;
        
        /**
        * Get end of segment table
        * @return segment table iterator
        */
        ConstIterator End() const;
        
        /**
        * Form graph after node table has been populated
        */
        void FormGraph();
        
        /**
        * Compute the graph attributes associated with the nodes
        */
        void ComputeNodeParameters();
        
        /**
        * Align the tangent vectors along a segment.  Have vectors point
        * from thick end to thin end
        */
        void AlignTangentVector();
        
        /**
        * Compute local coordinate system for each node 
        * i.e. tangent, normal,binormal vectors
        */
        void ComputeFrenetFrame();
        
        /**
        * Compute parameters associated with Frenet frame 
        */
        void ComputeFrenetParameters();
        
        /**
        * Compute parameters associated with segments
        */
        void ComputeSegmentParameters();
        
        /**
        * Compute the tangent by considering the positions neighboring nodes
        */
        void SetTangentByPosition();

        /**
        * Build a table of all node attributes
        */
        ParameterTableType GetNodeParameters() const;
        
        /**
        * Build a table of all segment attributes
        */
        ParameterTableType GetSegmentParameters() const;
        
        /**
        * Return the number of segments in the graph
        */
        unsigned int NumberOfSegments() const;
        
        /**
        * Return the number of junction points in the graph
        */
        unsigned int NumberOfJunctions() const;
        
        /**
        * Return the number of junction points in the graph with specified
        * number of neighbors
        */
        unsigned int NumberOfNJunctions(unsigned int n) const;
        
        /**
        * Return the number of nodes in the graph
        */
        unsigned int NumberOfNodes(const SegmentIDType&) const;  
        
        /**
        * Return the degree histogram of the graph
        */
        HistogramType GetHistogram() const;

        /**
        * Dump all graph metrics to a python file for plotting
        */
        std::ostream& ReportMetricsAsPython(std::ostream &os) const;
        
        /**
        * Remove the radii whose vessel probability is below the specified
        * threshold and replace them with a more sensible value
        * @param p probability threshold
        */
        void FixPoorRadii(RealType p);
        
        /**
        * Get points which define the corners of a box that the graph would
        * fit in
        */
        void GetBoundingPoints(PointType&, PointType&) const;
        
        /**
        * Write graph in .vessel format
        */
        std::ostream& Write(std::ostream& os) const;
        
        /**
        * Read graph in .vessel format
        */
        std::istream& Read(std::istream& is);
        
        /**
        * Write graph in .gml format
        */
        std::ostream& WriteGML(std::ostream& os) const;
        
        /**
        * Write graph defined by junction points only and edges as 
        * distance between junction nodes in .gml format
        */
        std::ostream& WriteGMLJunctionOnly(std::ostream& os) const;
        
        /**
        * Write graph in .graphml format
        */
        std::ostream& WriteGraphML(std::ostream& os) const;
    
    protected:
        /**
        * Fix possible problems in the node table
        */
        void RemoveExcessNodes();
        
        /**
        * Determine if graph contains segment 
        */
        bool IsSegmentUnique(SegmentType) const;
        
        /**
        * Fix associations between nodes and segment ids
        */
        void UpdateSegmentIDOnNodes(SegmentIDType);
        
        /**
        * Method used in sub graph recursion
        */
        void AddSegmentAndNeighbors(const typename SegmentType::IDType&, Self& ) const;
        
        /**
        * Get all neighboring nodes of the specified node
        * @param n node to find neighbors of
        * @param nn collection of neighbors
        */
        void GetNeighborJunctions(const NodeType n, NodeTableType& nn) const;

        /**
        * Set the segment length on the specified segment
        */ 
        void SetSegmentLength( SegmentType&) const;   
        
        /**
        * Set the segment chord length on the specified segment
        */  
        void SetSegmentChordLength( SegmentType&) const;
        
        /**
        * Set the sum of total angle on the specified segment
        */ 
        void SetSegmentSumOfTotalAngle( SegmentType&) const;
        
        /**
        * Set the average radius on the specified segment
        */
        void SetSegmentAverageRadius( SegmentType&) const;
        
        /**
        * Set the accumlated curvature on the specified segment
        */
        void SetSegmentSumMagnitudeOfCurvature( SegmentType&) const;
        
        /**
        * Set the accumlated curvature derivative on the specified segment
        */
        void SetSegmentSumMagnitudeOfDerivativeCurvature( SegmentType&) const;
        
        /**
        * Set the tortuoisty density on the specified segment
        */
        void SetSegmentTortuosityDensity( SegmentType&) const;
        
        /**
        * Set the number of inflection points on the specified segment
        */
        void SetSegmentInflectionCount( SegmentType&) const;
        
        //data structures
        SegmentTableType segmentTable;
        NodeTableType nodeTable;
};

template<unsigned int VIndexDimension, typename TCoordinateType>
std::ostream& operator<<(std::ostream &os, const VesselGraph<VIndexDimension, TCoordinateType> &vg);

template<unsigned int VIndexDimension, typename TCoordinateType>
std::istream& operator>>(std::istream &os, VesselSegment<VIndexDimension, TCoordinateType> &vg);

template<typename T>
std::ostream& ValueToPython(std::ostream &os, std::string, T);

template<typename RandomAccessIterator>
std::ostream& ArrarToPython(std::ostream &os, std::string, RandomAccessIterator, RandomAccessIterator);

#include "VesselGraph.txx"

#endif
