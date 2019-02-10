# Graphs-isomorphism

Testing implementation of searching for all graph isomorphisms of Hugenholtz diagrams.
Uses ideas of McKay and his Practical Graph Isomorphism algorithm without any automorphism pruning
see (<a href="https://arxiv.org/pdf/1301.1493.pdf" target="_blank">`arXiv:1301.1493`</a>).
One difference is the estimation of number of degrees because Hugenholtz diagrams can have multiple edges for the same
pair of vertices. The PseudoEdge class was created because of this see (Graph.h).

## Example 

```c++

Graph g(7); //five vertices

g1.insertEdge(1,2);
g1.insertEdge(1,4);
g1.insertEdge(2,5);
g1.insertEdge(3,1);
g1.insertEdge(4,5);
g1.insertEdge(5,7);
g1.insertEdge(6,3);
g1.insertEdge(6,7);
g1.insertEdge(4,7);

if(g1.isIsomoprhic(g1))std::cout << "graph is isomorphic" << std::endl;
else std::cout << "graph is not isomorphic" << std::endl;
```
## Author
Jakub Lang

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**
