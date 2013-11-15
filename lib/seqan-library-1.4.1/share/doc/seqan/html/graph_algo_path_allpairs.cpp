#include <iostream>
#include <seqan/graph_algorithms.h>

using namespace seqan;

int main() {
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;
	TSize numEdges = 9;
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	TGraph g;
	addEdges(g, edges, numEdges);
	std::cout << g << std::endl;
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};
	String<int> weightMap;
	assignEdgeMap(g,weightMap, weights);
	String<int> distMat;
	String<TVertexDescriptor> predMat;
	allPairsShortestPath(g,weightMap, distMat, predMat);
	unsigned int len = (unsigned int) std::sqrt((double) length(distMat));
	for (TSize row=0;row < len;++row) {
		for (TSize col=0;col < len;++col) {
			std::cout << row << "," << col << " (Distance=" << getValue(distMat, row*len + col) << "): "; 
			_printAllPairsShortestPath(g,predMat,row,col);
			std::cout << std::endl;
		}
	}
	return 0;
}
