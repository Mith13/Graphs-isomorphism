#include "Topology.hpp"
#include "Graph.h"
#include <map>
#include <memory>
using std::cout;
using std::endl;


int main()
{

  Graph g1(6), g2(6);

  g1.insertEdge(0,1);
  g1.insertEdge(2,3);
  g1.insertEdge(4,5);
  g1.insertEdge(1,3);
  g1.insertEdge(3,5);
  g1.insertEdge(5,0);
  g1.insertEdge(0,2);
  g1.insertEdge(2,4);
  g1.insertEdge(4,1);

  g2.insertEdge(0,4);
  g2.insertEdge(4,5);
  g2.insertEdge(5,2);
  g2.insertEdge(2,1);
  g2.insertEdge(1,3);
  g2.insertEdge(3,0);
  g2.insertEdge(1,4);
  g2.insertEdge(0,2);
  g2.insertEdge(3,5);
  
   if(g1.compare(g2))std::cout << "graph is isomorphic" << std::endl;
	 else std::cout << "graph is not isomorphic" << std::endl;

}
