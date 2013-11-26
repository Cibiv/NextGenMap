#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() 
{
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
	TSize numEdges = 10;
	TVertexDescriptor edges[] = {0,1, 0,4, 1,5, 2,5, 2,6, 2,3, 3,6, 3,7, 5,6, 6,7};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
	String<char> nameMap;
	char names[] = {'r', 's', 't', 'u', 'v', 'w', 'x', 'y'};
	assignVertexMap(g,nameMap, names);
	String<unsigned int> predMap;
	String<unsigned int> distMap;
	breadthFirstSearch(g, 1, predMap, distMap);
	std::cout << "Breadth-First search: " << std::endl;
	typedef Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	while(!atEnd(it)) {
		std::cout << "Vertex " << getProperty(nameMap, getValue(it)) << ": ";
		if (getProperty(distMap, getValue(it))== _getInfinityDistance(distMap)) {
			std::cout << "Not reachable!";
		} else {
			std::cout << "Level = " << getProperty(distMap, getValue(it));
		}
		typedef Value<String<unsigned int> >::Type TPredVal;
		TPredVal pre = getProperty(predMap, getValue(it));
		if (pre != getNil<TVertexDescriptor>()) {
			std::cout << ", Predecessor = " << getProperty(nameMap, pre) << std::endl;
		} else {
			std::cout << ", Predecessor = nil" << std::endl;
		}
		goNext(it);
	}
	return 0;
}
