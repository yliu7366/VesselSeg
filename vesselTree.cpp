
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <limits>

#include "VesselNode.h"
#include "VesselSegment.h"
#include "VesselGraph.h"

#include "itkIndex.h"

#include <boost/config.hpp>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/wavefront.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/connected_components.hpp>

typedef double CoordinateType;
const unsigned int Dimension = 3; 
typedef VesselGraph<Dimension, CoordinateType> VesselGraphType; 

namespace boost {
    enum edge_id_t { edge_id };
    enum vertex_id_t { vertex_id };
    enum vertex_centralbetween_t { vertex_centralbetween };

    BOOST_INSTALL_PROPERTY(edge, id);
    BOOST_INSTALL_PROPERTY(vertex, id);
    BOOST_INSTALL_PROPERTY(vertex, centralbetween);
}

int main(int argc, char ** argv)
{

    std::string inputName(argv[1]);
    std::ifstream inputfile;
    inputfile.open(inputName.c_str(), std::ios::in);
    VesselGraphType vesselGraph;
    inputfile >> vesselGraph;
    inputfile.close();

   
    //Form Boost Graph    
    unsigned numSegments = vesselGraph.NumberOfSegments();
    //std::cout << "Number of Segments " << numSegments << std::endl;
    typedef boost::property < boost::vertex_in_degree_t, int > VertexDegreeType;
    typedef boost::property < boost::vertex_id_t, long int> VertexPropertyType;
    
    typedef boost::property < boost::edge_weight_t, double >Weight;
    typedef boost::property<boost::edge_id_t, long int, Weight> EdgePropertyType;

    typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS, VertexPropertyType, EdgePropertyType > UndirectedGraph;
    UndirectedGraph undigraph;//(numSegments+1);
    
    long int nodeCounter = 0;
    long int numEndPoints = 0;
    
    for(VesselGraphType::ConstIterator segmentIter = vesselGraph.Begin(); segmentIter != vesselGraph.End(); ++segmentIter)
    {
        VesselGraphType::SegmentIDType currentSegmentID = segmentIter->first;
        VesselGraphType::SegmentType currentSegment = segmentIter->second;
        VesselGraphType::NodeIDType frontID = currentSegment.Front();
        VesselGraphType::NodeType frontNode = vesselGraph.GetNode(frontID);
        
        VesselGraphType::NodeIDType backID = currentSegment.Back();
        VesselGraphType::NodeType backNode = vesselGraph.GetNode(backID);
        
        if(frontNode.GetType() == VesselGraphType::NodeType::SINGLE)
        {
            //boost::graph_traits < UndirectedGraph >::vertex_descriptor u = boost::vertex(frontID, undigraph);
            //boost::add_vertex(u, undigraph);
        }
        else
        {
        
            boost::graph_traits<UndirectedGraph>::vertex_iterator vertexIter;
            boost::graph_traits<UndirectedGraph>::vertex_iterator vertexEnd;
            boost::tie(vertexIter, vertexEnd) = boost::vertices(undigraph);
            boost::graph_traits < UndirectedGraph >::vertex_descriptor u; 
            boost::graph_traits < UndirectedGraph >::vertex_descriptor v;
            bool foundFront = false;
            bool foundBack = false;
            
            
            for(; vertexIter != vertexEnd; ++vertexIter)
            {
                boost::graph_traits < UndirectedGraph >::vertex_descriptor currentVertex = *vertexIter;
            
                long vertexIDProprety = boost::get(boost::vertex_id, undigraph, currentVertex);
                if(vertexIDProprety == frontID)
                {
                    u = currentVertex;
                    foundFront = true;
                }
                if(vertexIDProprety == backID)
                {
                    v = currentVertex;
                    foundBack = true;
                }
                
            }
            if(foundFront == false)
            {
                boost::add_vertex(undigraph);
                u = boost::vertex(nodeCounter, undigraph);
                boost::put(boost::vertex_id, undigraph, u, frontID);
                nodeCounter++;
            }
            if(foundBack == false)
            {
                boost::add_vertex(undigraph);
                v = boost::vertex(nodeCounter, undigraph);
                boost::put(boost::vertex_id, undigraph, v, backID);
                nodeCounter++;
            }
            
            //boost::graph_traits < UndirectedGraph >::vertex_descriptor u = boost::vertex(frontID, undigraph);
            //boost::graph_traits < UndirectedGraph >::vertex_descriptor v = boost::vertex(backID, undigraph);
            
            boost::graph_traits < UndirectedGraph >::edge_descriptor addedDescrip;
            bool didAdd;
            boost::tie(addedDescrip, didAdd) = boost::add_edge( u, v, undigraph);
            boost::graph_traits<UndirectedGraph>::vertex_descriptor source = boost::source(addedDescrip, undigraph);
            boost::graph_traits<UndirectedGraph>::vertex_descriptor target = boost::target(addedDescrip, undigraph);
            
            boost::graph_traits < UndirectedGraph >::edge_descriptor descrip;
            bool isInGraph;
            boost::tie( descrip, isInGraph ) = boost::edge(u,v,undigraph);
            boost::put(boost::edge_weight, undigraph, descrip, currentSegment.GetLength());
            //boost::put(boost::edge_weight, undigraph, descrip, 1.0);
            boost::put(boost::edge_id, undigraph, descrip, currentSegmentID);
        }
        if(frontNode.GetType() == VesselGraphType::NodeType::TERMINAL)
        {
            ++numEndPoints;    
        }
        if(backNode.GetType() == VesselGraphType::NodeType::TERMINAL)
        {
            ++numEndPoints;        
        }
        for(VesselGraphType::SegmentType::ConstIterator nodeIter = currentSegment.Begin(); nodeIter != currentSegment.End(); ++nodeIter)
        {
            VesselGraphType::SegmentType::NodeType::IDType currentNodeID = *(nodeIter);
            VesselGraphType::SegmentType::NodeType currentNode = vesselGraph.GetNode(currentNodeID);
        }
    }
    
    int numEdges = boost::num_edges(undigraph);
    long int V = boost::num_vertices(undigraph);
    //std::cout << "Number of Vertices " << V << std::endl;
    std::vector < double > aRow(V, 0.0);
    std::vector < std::vector<double> > matrix(V, aRow);
    boost::johnson_all_pairs_shortest_paths(undigraph, matrix);
    
    double averagePath = 0.0;
    int numPaths = 0;
    //std::cout << "Shortest Path " << std::endl;
    for(unsigned int i = 0; i < V; i++)
    {
        for(unsigned int j = 0; j < i; j++)
        {
            double currentPathLength = matrix[i][j];
            //std::cout << "(" << i << "," << j << ")=" << currentPathLength << " ";
            if(currentPathLength < std::numeric_limits<double>::max())
            {
                averagePath += currentPathLength;
                numPaths++;
            }
        }
        
    }
    //std::cout << std::endl;
    averagePath = averagePath/numPaths;

    std::cout << "Average Path Length " << averagePath << std::endl;
    
    std::vector<double> centrality(V);
    boost::brandes_betweenness_centrality(undigraph, &centrality[0]);
    std::vector<double> absoluteCentrality(centrality);
    boost::relative_betweenness_centrality(undigraph, &centrality[0]);
    std::vector<double> relativeCentrality(centrality);
    double centralDominance = boost::central_point_dominance(undigraph, &relativeCentrality[0]);
    
    std::cout << "Central Point Dominance = " << centralDominance << std::endl;
    
    for(unsigned int i = 0; i < V; i++)
    {
        //std::cout << "Vertex index "<< i << " Abs Cent " << absoluteCentrality[i] << " Rel Cent " << relativeCentrality[i] << std::endl;
    }

    
    
    std::vector< boost::graph_traits < UndirectedGraph >::vertex_descriptor > art_points;
    boost::articulation_points(undigraph, std::back_inserter(art_points));
    if(art_points.size() == 0)
    {
        //std::cout << "No Art Points " << std::endl;
    }
    for(unsigned int i = 0; i < art_points.size(); i++)
    {
        //std::cout << "Art Point " << art_points[i] << std::endl;
    }   
    
    std::cout << "Number of Articulation Points " << art_points.size() << std::endl;
    double articulationPercent = double(art_points.size())/double(V) ;
    std::cout << "Percent of Articulation Points " << articulationPercent << std::endl;
    std::cout << "Number of Terminal Points " << numEndPoints << std::endl;
    
    std::vector<int> component(V);
    int numComponents = boost::connected_components(undigraph, &component[0]);
    std::cout << "Number of Components " << numComponents << std::endl;
    double relativeComponents = double(numComponents) / double(V);
    
    int bandwidth = boost::bandwidth(undigraph);
    std::cout << "Bandwidth " << bandwidth << std::endl;
    int maxWavefront = boost::max_wavefront(undigraph);
    std::cout << "Max Wavefront " << maxWavefront << std::endl;
    int avgWavefront = boost::aver_wavefront(undigraph);
    std::cout << "Avg Wavefront " << avgWavefront << std::endl;
    int rmsWavefront = boost::rms_wavefront(undigraph);
    std::cout << "RMS Wavefront " << rmsWavefront << std::endl;
    
    double averageHubWeight =0.0;

    std::vector<int> degree;
    boost::graph_traits<UndirectedGraph>::vertex_iterator vertexIter;
    boost::graph_traits<UndirectedGraph>::vertex_iterator vertexEnd;
    boost::tie(vertexIter, vertexEnd) = boost::vertices(undigraph);
    for(; vertexIter != vertexEnd; ++vertexIter)
    { 
        boost::graph_traits<UndirectedGraph>::vertex_descriptor currentVertex = *vertexIter;
        int currentDegreeOut = boost::out_degree(currentVertex, undigraph);
        int currentDegreeIn = boost::in_degree(currentVertex, undigraph);
        degree.push_back(currentDegreeOut);
        
        if(currentDegreeOut >= 3)
        {
            boost::graph_traits<UndirectedGraph>::out_edge_iterator currentEdgeIter;
            boost::graph_traits<UndirectedGraph>::out_edge_iterator currentEdgeEnd;
            boost::tie(currentEdgeIter, currentEdgeEnd) = boost::out_edges(currentVertex, undigraph);
            std::vector<double> currentVertexWeights;
            for(; currentEdgeIter != currentEdgeEnd; ++currentEdgeIter)
            {
                boost::graph_traits<UndirectedGraph>::edge_descriptor currentEdge = *currentEdgeIter;
                double currentWeight = boost::get(boost::edge_weight, undigraph, currentEdge);
                if(currentWeight > 0)
                {                
                    currentVertexWeights.push_back(currentWeight);
                }
            }
            double minWeight = *std::min_element(currentVertexWeights.begin(), currentVertexWeights.end(), std::less_equal<double>() );
            averageHubWeight += double(currentDegreeOut) / (minWeight);
        }
        //std::cout << "Vertex ID " << currentVertex << " Degree Out " << currentDegreeOut << " Degree In " << currentDegreeIn << std::endl;
    }
    
    averageHubWeight = averageHubWeight/double(V);
    std::cout << "Hub Weight " << averageHubWeight << std::endl;
    int hubDegree = *std::max_element(degree.begin(), degree.end(), std::less_equal<int>() );
    std::cout << "Hub Degree " << hubDegree << std::endl;
    
    //Cluster Coeff
    double clusterCoeff = 0.0;
    std::vector<double> cluster(V);
    boost::tie(vertexIter, vertexEnd) = boost::vertices(undigraph);
    for(; vertexIter != vertexEnd; ++vertexIter)
    { 
        boost::graph_traits<UndirectedGraph>::vertex_descriptor currentVertex = *vertexIter;
        int currentDegree = boost::out_degree(currentVertex, undigraph);
        
        std::vector<int> currentNeighbors;
        boost::graph_traits<UndirectedGraph>::out_edge_iterator currentEdgeIter;
        boost::graph_traits<UndirectedGraph>::out_edge_iterator currentEdgeEnd;
        boost::tie(currentEdgeIter, currentEdgeEnd) = boost::out_edges(currentVertex, undigraph);
        for(; currentEdgeIter != currentEdgeEnd; ++currentEdgeIter)
        {
            boost::graph_traits<UndirectedGraph>::edge_descriptor currentEdge = *currentEdgeIter;
            boost::graph_traits<UndirectedGraph>::vertex_descriptor target = boost::target(currentEdge, undigraph);
            currentNeighbors.push_back(target);
        }

        int currentClusterCoeff = 0;
        boost::tie(currentEdgeIter, currentEdgeEnd) = boost::out_edges(currentVertex, undigraph);
        for(; currentEdgeIter != currentEdgeEnd; ++currentEdgeIter)
        {
            boost::graph_traits<UndirectedGraph>::edge_descriptor currentEdge = *currentEdgeIter;
            boost::graph_traits<UndirectedGraph>::vertex_descriptor target2 = boost::target(currentEdge, undigraph);
            boost::graph_traits<UndirectedGraph>::out_edge_iterator currentEdgeIter2;
            boost::graph_traits<UndirectedGraph>::out_edge_iterator currentEdgeEnd2;
            boost::tie(currentEdgeIter2, currentEdgeEnd2) = boost::out_edges(target2, undigraph);
            for(; currentEdgeIter2 != currentEdgeEnd2; ++currentEdgeIter2)
            {
                boost::graph_traits<UndirectedGraph>::edge_descriptor currentEdge2 = *currentEdgeIter2;
                boost::graph_traits<UndirectedGraph>::vertex_descriptor target3= boost::target(currentEdge2, undigraph);
                std::vector<int>::iterator findIter = std::find(currentNeighbors.begin(), currentNeighbors.end(), target3);
                if(findIter != currentNeighbors.end())
                {
                    //std::cout << " found it " << " (" << currentVertex << "," << target2 << ","<< target3 << ") " << std::endl;
                    ++currentClusterCoeff;
                }
            }
            
        }
        //std::cout << "Vertex ID " << currentVertex << " Num triangles " << currentClusterCoeff << std::endl;
        if(currentDegree > 1)
        {
            double currentValue = double(currentClusterCoeff)/(double(currentDegree)*double(currentDegree-1));
            clusterCoeff += currentValue;
            cluster[currentVertex] = currentValue;
        }
        else
        {
            cluster[currentVertex] = 0.0;
        }
    }
    clusterCoeff = clusterCoeff / double(V);
    std::cout << "Cluster Coeff " << clusterCoeff << std::endl;

    
    std::vector<int> degreeHistogram(hubDegree+1, 0);
    for(unsigned int i = 0; i < degree.size(); i++)
    {
        int currentDegree = degree[i];
        degreeHistogram[currentDegree] =  degreeHistogram[currentDegree] + 1;
    }

    std::cout << "Degree Histogram ";
    std::copy(degreeHistogram.begin(), degreeHistogram.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    
    double entropy = 0.0;
    for(unsigned int i = 0; i < degreeHistogram.size(); i++)
    {
        int currentNumber = degreeHistogram[i];
        if(currentNumber != 0)
        {
            double currentFaction = double(currentNumber)/double(V);
            entropy += currentFaction*log(currentFaction);
        }
    }
    entropy = -entropy;
    std::cout << "Entropy " << entropy << std::endl;
    
    std::cout << "[" << averagePath << ", " << centralDominance << ", " << articulationPercent << ", " << bandwidth << ", " << maxWavefront << ", " << avgWavefront << ", " << rmsWavefront << ", " << hubDegree << ", " << averageHubWeight << ", " << clusterCoeff << ", " << entropy  << ", " << relativeComponents << "]" << std::endl;
    
    std::string outputName(argv[2]);
    std::ofstream outputfile;
    outputfile.open(outputName.c_str(), std::ios::out);
    outputfile << "import numpy" << std::endl;
    outputfile << "metrics = numpy.array([" << averagePath << ", " << centralDominance << ", " << articulationPercent << ", " << bandwidth << ", " << maxWavefront << ", " << avgWavefront << ", " << rmsWavefront << ", " << hubDegree << ", " << averageHubWeight << ", " << clusterCoeff << ", " << entropy << ", " << relativeComponents << "])" << std::endl;
    outputfile.close();


/*
    //Write graph
    std::ofstream fout("TestGraphWrite.dot");
    fout << "graph A {\n"
    << "  rankdir=LR\n"
    << "size=\"8,8\"\n"
    << "ratio=\"fill\"\n"
    << "edge[style=\"bold\"]\n" << "node[shape=\"circle\"]\n";
    
    
    boost::graph_traits<UndirectedGraph>::edge_iterator edgeIter;
    boost::graph_traits<UndirectedGraph>::edge_iterator edgeEnd;
    boost::tie(edgeIter, edgeEnd) = boost::edges(undigraph);
    for(; edgeIter != edgeEnd; ++edgeIter)
    {   
        boost::graph_traits<UndirectedGraph>::edge_descriptor currentEdge = *edgeIter;
        boost::graph_traits<UndirectedGraph>::vertex_descriptor source = boost::source(currentEdge, undigraph);
        boost::graph_traits<UndirectedGraph>::vertex_descriptor target = boost::target(currentEdge, undigraph);
        long sourceID = boost::get(boost::vertex_id, undigraph, source);
        long targetID = boost::get(boost::vertex_id, undigraph, target);
        double currentWeight = boost::get(boost::edge_weight, undigraph, currentEdge);
        
        //std::cout << std::endl;
        //std::cout << "Source " << source << " Target " << target << std::endl;
        //std::cout << "Source ID " << sourceID << " Target ID " << targetID << std::endl;
        //std::cout << "Source Degree " << boost::in_degree(source, undigraph) << std::endl;
        //std::cout << "Target Degree " << boost::in_degree(target, undigraph) << std::endl;
        //std::cout << "Edge Weight " << currentWeight << std::endl;

        
        //fout << sourceID << " -- " << targetID << "[label=" << currentWeight << "]\n";
        fout << source << " -- " << target << "[label=" << currentWeight << "]\n";
    }
    fout << "}\n";
    fout.close();
    std::string outputFilename = "TestGraphWrite.dot";
    std::ofstream outputFile;
    outputFile.open( outputFilename.c_str(), std::ios::out);

    boost::write_graphviz(outputFile, undigraph);

    outputFile.close();
    */
    
    
    return 0;
}
