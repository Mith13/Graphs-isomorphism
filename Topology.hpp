/******************************************************************************
 *                                                                            *
 * Copyright (C) 2002-2005 Michal Czakon                                      *
 *                                                                            *
 ******************************************************************************/

#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP 1

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>
#include <string>
#include <stack>
#include <list>
#include <set>
using std::ostream;
using std::vector;
using std::string;
using std::stack;
using std::list;
using std::multiset;
using std::set;
using std::pair;

#ifdef DEBUG
#include "Debug.hpp"
#endif

class AdjacentEdges;
class AdjacentNodes;
class TopologyComponent;

/******************************************************************************
 *                                                                            *
 * Topology                                                                   *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The QFT topology is considered to be a directed pseudo-graph (multiple edges
 * and loops, where both the source and target nodes of an edge are the same,
 * are allowed). Nodes are elements of the set {0, ..., n_nodes-1}, which must
 * be specified at construction. Edges may carry a momentum, but changing its
 * value is not considered to modify the topology. Node labelling and topology
 * comparison can be done either with equivalent external nodes or with their
 * ordering defined by integer comparison.
 *
 */

class Topology
{
/******************************************************************************
 ******************************************************************************
 **                                                                          **
 ** Public interface part                                                    **
 **                                                                          **
 ******************************************************************************
 ******************************************************************************/

public:

  /****************************************************************************
   *                                                                          *
   * construction and destruction                                             *
   *                                                                          *
   ****************************************************************************/
   
  explicit Topology(int n_nodes);
  
  virtual ~Topology();

  /****************************************************************************
   *                                                                          *
   * modification                                                             *
   *                                                                          *
   ****************************************************************************/
   
  virtual void
  insert_edge(int source_node, int target_node);
  
  virtual void
  insert_edge(int source_node,
	      int target_node,
	      const string& momentum,
	      int direction = +1);

  virtual void
  erase_edge(int);

  virtual void
  erase_edges(int first_node, int second_node);

  /****************************************************************************
   *                                                                          *
   * node and edge properties                                                 *
   *                                                                          *
   ****************************************************************************/
   
  int
  n_nodes() const;

  int
  n_external_nodes() const;

  int
  n_edges() const;

  int
  n_external_edges() const;

  vector<int>
  external_nodes() const;

  bool
  is_external_node(int) const;

  AdjacentEdges
  external_edges() const;
  
  bool
  is_external_edge(int) const;

  AdjacentEdges
  adjacent_edges(int node) const;

  AdjacentNodes
  adjacent_nodes(int edge) const;
  
  int
  adjacency(int first_node, int second_node) const;

  /****************************************************************************
   *                                                                          *
   * momenta                                                                  *
   *                                                                          *
   ****************************************************************************/
   
  void
  assign_momenta(const string& external_momentum_prefix = "p",
		 const string& loop_momentum_prefix = "k") const;
  
  void
  assign_momentum(int edge, const string& momentum, int direction = +1) const;
  
  void
  assign_momentum(int edge, const vector<int>& momentum) const;

  void
  assign_external_momentum(int edge,
			   const string& momentum,
			   int direction) const;
  
  void
  copy_external_momentum(int edge, const Topology& t, int t_edge) const;

  void
  clear_momenta() const;

  vector<int>
  momentum(int edge) const;
  
  void
  set_momentum_basis(const vector<string>&) const;

  vector<string>
  get_momentum_basis() const;

  /****************************************************************************
   *                                                                          *
   * topological properties                                                   *
   *                                                                          *
   ****************************************************************************/
   
  int
  n_cycles() const;

  bool
  is_connected() const;

  bool
  is_one_particle_irreducible() const;

  bool
  has_tadpoles() const;

  bool
  has_self_energies() const;

  bool
  is_on_shell() const;

  vector<TopologyComponent>
  biconnected_components() const;

  vector<vector<int> >
  node_symmetry_group() const;

  vector<vector<int> >
  internal_node_symmetry_group() const;

  vector<vector<int> >
  independent_external_node_permutations() const;
  
  /****************************************************************************
   *                                                                          *
   * comparison                                                               *
   *                                                                          *
   ****************************************************************************/
   
  void
  fix_external_nodes() const;

  void
  free_external_nodes() const;

  bool
  fixed_external_nodes() const;

  vector<int>
  node_labelling() const; 

  virtual int
  compare(const Topology& t) const;

  /****************************************************************************
   *                                                                          *
   * output                                                                   *
   *                                                                          *
   ****************************************************************************/
   
  virtual void
  print(ostream& output = std::cout) const;

  virtual void
  print_edge_list(ostream& output = std::cout) const;

  virtual void
  postscript_print(const string& file_name) const;
  
  void
  print_adjacency_matrix(ostream& output = std::cout) const;

  void
  print_adjacency_list(ostream& output = std::cout) const;

  void
  print_momentum(ostream& output, int edge) const;

  void
  print_momentum(ostream& output, int edge, const string& index) const;

  void
  print_momentum(ostream& output, const vector<int>& momentum) const;

  void
  print_momentum(ostream& output,
		 const vector<int>& momentum,
		 const string& index) const;

  list<set<int> >
	  free_equitable_partition() const;

/******************************************************************************
 ******************************************************************************
 **                                                                          **
 ** Private implementation part                                              **
 **                                                                          **
 ******************************************************************************
 ******************************************************************************/

protected:

  vector<AdjacentEdges> _node;

  vector<AdjacentNodes> _edge;

private:

  /// Topology must be constructed with a known number of nodes
  Topology();

  vector<int>
  subgraph(int starting_node, const vector<bool>& forbidden_edge) const;

  vector<vector<int> >
  symmetry_group(const list<set<int> >& partition) const;

  struct BCCStatus
  {
    BCCStatus(int n_nodes, int n_edges);

    vector<TopologyComponent> _complete;

    /// the first element is a node or edge, while the second is true for nodes.
    stack<pair<int, bool> >   _incomplete;
    
    vector<bool>               _articulation_point;

    int                        _n_external_nodes;

    int                       _count;

    vector<int>               _parent;

    vector<int>               _index;

    vector<int>               _min_index;

    vector<bool>              _used_edge;

    vector<bool>              _visited_node;
  };
  
  void
  bcc(BCCStatus& status, int node) const;

  /// the generalized degree of a node with respect to a set of nodes
  class NodeDegree
  {
  public:
    
    NodeDegree();
    
    /// n should be the adjacency of two nodes, or minus adjacency for loops
    void
    insert_adjacency(int n);
    
    bool
    operator<(const NodeDegree&) const;

	int           _total_degree;
  private:
    
    /// graph theoretical degree of the node
    
    int           _key_size;
    
    /// adjacencies
    multiset<int> _key;
  };

  list<set<int> >
  fixed_equitable_partition() const;

  bool
  refine(list<set<int> >& partition, list<set<int> >& active) const;

  struct PartitionStatus
  {
    list<set<int> >           _partition;
    
    list<set<int> >::iterator _cell;
    
    set<int>::iterator        _node;
  };
  
  void
  init_partition_stack(stack<PartitionStatus>& partition_stack,
		       const list<set<int> >& partition) const;

  vector<int>
  next_leaf(stack<PartitionStatus>& partition_stack) const;

  /// compares the adjacency matrix transformed by two node permutations
  int
  compare_adjacency_matrix(const vector<int>& first_permutation,
			   const vector<int>& second_permutation) const;
  
private:

  vector<vector<int> >         _adjacency_matrix;

  mutable vector<vector<int> > _momentum;

  mutable vector<string>       _momentum_basis;

  mutable bool                 _fixed_external_nodes;

  // needed for consistency of different versions of momentum assignement
  mutable bool                 _user_assigned_momenta;

  // implementation speedup, since labelling is needed for topology comparison
  mutable vector<int>          _node_labelling;
};

/******************************************************************************
 *                                                                            *
 * AdjacentEdges                                                              *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Wrapper for a pair of vectors with the edges adjacent to a node. Ingoing
 * edges come first. Direction is either 0 (ingoing) or 1 (outgoing).
 *
 */

struct AdjacentEdges
{
  vector<int> _edges[2];
};

/******************************************************************************
 *                                                                            *
 * AdjacentNodes                                                              *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Wrapper for the pair of nodes adjacent to an edge. The source node comes
 * first. Direction is either 0 (source) or 1 (target).
 *
 */

struct AdjacentNodes
{
  AdjacentNodes(int source_node, int target_node);

  int _node[2];
};

/******************************************************************************
 *                                                                            *
 * TopologyComponent                                                          *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * A topology component is a subgraph of the topology. The attribute _vaccum is
 * true iff there is no external momentum flow through the component.
 *
 */

struct TopologyComponent
{
  TopologyComponent();

  vector<int> _nodes;

  vector<int> _edges;

  vector<int> _articulation_points;

  bool        _vacuum;
};

/******************************************************************************
 *                                                                            *
 * Inlines                                                                    *
 *                                                                            *
 ******************************************************************************/

inline
Topology::Topology(int n_nodes) :
  _node(n_nodes),
  _adjacency_matrix(n_nodes, vector<int>(n_nodes)),
  _fixed_external_nodes(false),
  _user_assigned_momenta(true)
{
  // If the implementation of the vector container provides a small default
  // space then during edge insertions the container is reallocated. There is a
  // small gain of performance if we reserve enough space at the beginning. We
  // use an upper bound of twice the number of nodes, which is based on the fact
  // that this is the number of edges that corresponds to a topology on
  // quadruple nodes only.

  _edge.reserve(2*n_nodes);
}

inline
Topology::~Topology()
{}

inline
int
Topology::n_nodes() const
{
  return _node.size();
}

inline
int
Topology::n_external_nodes() const
{
  return external_nodes().size();
}

inline
int
Topology::n_edges() const
{
  return _edge.size();
}

inline
int
Topology::n_external_edges() const
{
  return external_nodes().size();
}

inline
bool
Topology::is_external_node(int n) const
{
  return _node[n]._edges[0].size()+_node[n]._edges[1].size() == 1;
}

inline
bool
Topology::is_external_edge(int e) const
{
  return is_external_node(_edge[e]._node[0]) ||
    is_external_node(_edge[e]._node[1]);
}

inline
AdjacentEdges
Topology::adjacent_edges(int n) const
{
  return _node[n];
}

inline
AdjacentNodes
Topology::adjacent_nodes(int e) const
{
  return _edge[e];
}

inline
int
Topology::adjacency(int n1, int n2) const
{
  return _adjacency_matrix[n1][n2];
}

inline
vector<int>
Topology::momentum(int e) const
{
  return _momentum[e];
}

inline
void
Topology::set_momentum_basis(const vector<string>& basis) const
{
  clear_momenta();
  _momentum_basis = basis;
}

inline
vector<string>
Topology::get_momentum_basis() const
{
  return _momentum_basis;
}

inline
bool
Topology::is_connected() const
{
  return subgraph(0, vector<bool>(n_edges(), false)).size() == n_nodes();
}

inline
vector<vector<int> >
Topology::node_symmetry_group() const
{
  return symmetry_group(free_equitable_partition());
}

inline
vector<vector<int> >
Topology::internal_node_symmetry_group() const
{
  return symmetry_group(fixed_equitable_partition());
}

inline
bool
Topology::fixed_external_nodes() const
{
  return _fixed_external_nodes;
}

inline
Topology::BCCStatus::BCCStatus(int n_nodes, int n_edges) :
  _articulation_point(n_nodes, false),
  _n_external_nodes(0),
  _count(0),
  _parent(n_nodes),
  _index(n_nodes),
  _min_index(n_nodes),
  _used_edge(n_edges, false),
  _visited_node(n_nodes, false)
{}

inline
AdjacentNodes::AdjacentNodes(int source_node, int target_node)
{
  _node[0] = source_node;
  _node[1] = target_node;
}

inline
TopologyComponent::TopologyComponent() : _vacuum(false)
{}

inline
Topology::NodeDegree::NodeDegree() : _total_degree(0), _key_size(0)
{}

inline
void
Topology::NodeDegree::insert_adjacency(int n)
{
  if (!n) return;

  if (n > 0) _total_degree += n;
  else _total_degree -= n;
  
  ++_key_size;
  _key.insert(n);
}

inline
bool
Topology::NodeDegree::operator<(const NodeDegree& d) const
{
  if (_total_degree != d._total_degree) return _total_degree < d._total_degree;
  if (_key_size != d._key_size) return _key_size < d._key_size;
  return _key < d._key;
}

inline
bool
isomorphic(const Topology& t1, const Topology& t2)
{
  return t1.compare(t2) == 0;
}

inline
bool
operator<(const Topology& t1, const Topology& t2)
{
  return t1.compare(t2) < 0;
}

inline
ostream&
operator<<(ostream& output, const Topology& t)
{
  t.print(output);
  return output;
}

inline
ostream&
operator<<(ostream& output, const vector<int>& v)
{
  copy(v.begin(), v.end(), std::ostream_iterator<int>(output, " "));
  output << "\n";
  return output;
}

inline
ostream&
operator<<(ostream& output, const vector<vector<int> >& m)
{
  for (vector<vector<int> >::const_iterator r = m.begin(); r != m.end(); ++r)
    output << *r;
  return output;
}

#endif
