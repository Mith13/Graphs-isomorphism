/* This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include "Graph.h"

int main()
{
  Graph g1(7), g2(7);
//////////GRAPH 1////////////////
  g1.insertEdge(1,2);
  g1.insertEdge(1,4);
  g1.insertEdge(2,5);

  g1.insertEdge(3,1);
  g1.insertEdge(4,5);
  g1.insertEdge(5,7);

  g1.insertEdge(6,3);
  g1.insertEdge(6,7);
  g1.insertEdge(4,7);
  g1.printAdjacencyMatrix();
//////////GRAPH 2////////////////
  g2.insertEdge(1,2);
  g2.insertEdge(1,6);
  g2.insertEdge(1,5);

  g2.insertEdge(2,3); 
  g2.insertEdge(3,7);
  g2.insertEdge(7,4);

  g2.insertEdge(7,6);
  g2.insertEdge(4,5);
  g2.insertEdge(5,6);
  
   if(g1.isIsomoprhic(g2))std::cout << "graph is isomorphic" << std::endl;
	 else std::cout << "graph is not isomorphic" << std::endl;
  
}
