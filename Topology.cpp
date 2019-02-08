/******************************************************************************
 *                                                                            *
 * Copyright (C) 2002-2005 Michal Czakon                                      *
 *                                                                            *
 ******************************************************************************/

#include <iterator>
#include <sstream>
#include <fstream>
#include <map>

#ifdef _WIN32
#include <io.h>
#endif
#ifdef linux
#include <unistd.h>
#endif

#include "Topology.hpp"
using namespace std;
void printSet(list<set<int>> set, std::string s = "")
{
	std::cout << s;
	for (auto s : set) {
		std::cout << std::endl;
		for (auto e : s) {
			std::cout << " " << e + 1;
		}
	}
	std::cout << std::endl;
}
/******************************************************************************
 *                                                                            *
 * insert_edge                                                                *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Inserting an edge without explicit momentum assignement invalidates the
 * node_labelling and the momentum distribution if it has been automatically
 * determined with assign_momenta.
 *
 */

void Topology::insert_edge(int n1, int n2)
{
  if (n1 >= n_nodes() || n1 < 0 || n2 >= n_nodes() || n2 < 0)
  {
    cerr << "Topology::insert_edge: at least one of the nodes " << n1 << ", "
         << n2 << " does not belong to a topology on " << n_nodes()
         << " nodes\n";

    exit(1);
  }

  //if (!_user_assigned_momenta) clear_momenta();
  _node_labelling.clear();

  _node[n1]._edges[1].push_back(n_edges());
  _node[n2]._edges[0].push_back(n_edges());

  _edge.push_back(AdjacentNodes(n1, n2));

  ++_adjacency_matrix[n1][n2];
  ++_adjacency_matrix[n2][n1];
}

/**
 *
 * Inserting an edge with momentum assignement invalidates the node_labelling
 * and the momentum distribution if it has been automatically determined with
 * assign_momenta. The direction of the momentum is defined with respect to the
 * edge.
 *
 */

void Topology::insert_edge(int n1, int n2, const string &p, int dir)
{
  insert_edge(n1, n2);

  if (p != "0")
  {
    _momentum.resize(n_edges());
    _momentum.back().resize(n_edges());
    _momentum.back().back() = dir;
  }
  _momentum_basis.resize(n_edges());
  _momentum_basis.back() = p;
}

/******************************************************************************
 *                                                                            *
 * erase_edge                                                                 *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Erasing an edge invalidates the momentum distribution and the node labelling
 *
 */

void Topology::erase_edge(int e)
{
  if (e >= n_edges() || e < 0)
  {
    cerr << "Topology::erase_edge: the edge " << e
         << " does not belong to a topology on " << n_edges() << " edges\n";

    exit(1);
  }

  //clear_momenta();
  _node_labelling.clear();

  --_adjacency_matrix[_edge[e]._node[0]][_edge[e]._node[1]];
  --_adjacency_matrix[_edge[e]._node[1]][_edge[e]._node[0]];

  _edge.erase(_edge.begin() + e);

  for (vector<AdjacentEdges>::iterator n = _node.begin(); n != _node.end(); ++n)
    for (int d = 0; d < 2; ++d)
      for (vector<int>::iterator ep = n->_edges[d].begin();
           ep != n->_edges[d].end(); ++ep)
        if (*ep == e)
          n->_edges[d].erase(ep--);
        else if (*ep > e)
          --*ep;
}

/******************************************************************************
 *                                                                            *
 * erase_edges                                                                *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Erases all the edges between the two nodes irrespective of their direction.
 *
 */

void Topology::erase_edges(int n1, int n2)
{
  if (n1 >= n_nodes() || n1 < 0 || n2 >= n_nodes() || n2 < 0)
  {
    cerr << "Topology::erase_edges: at least one of the nodes " << n1 << ", "
         << n2 << " does not belong to a topology on " << n_nodes()
         << " nodes\n";

    exit(1);
  }

  for (;;)
  {
    bool erased_edge = false;

    for (int d = 0; d < 2; ++d)
      for (vector<int>::iterator e = _node[n1]._edges[d].begin();
           e != _node[n1]._edges[d].end(); ++e)
        if (_edge[*e]._node[d] == n2)
        {
          erase_edge(*e);
          erased_edge = true;
          break;
        }

    if (!erased_edge)
      return;
  }
}

/******************************************************************************
 *                                                                            *
 * external_nodes                                                             *
 *                                                                            *
 ******************************************************************************/

vector<int>
Topology::external_nodes() const
{
  vector<int> nodes;

  for (int n = 0; n < n_nodes(); ++n)
    if (is_external_node(n))
      nodes.push_back(n);

  return nodes;
}

/******************************************************************************
 *                                                                            *
 * external_edges                                                             *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Similarly to edges adjacent to a node, edges ingoing with respect to the
 * topology come first.
 *
 */

AdjacentEdges
Topology::external_edges() const
{
  AdjacentEdges external;

  for (vector<AdjacentEdges>::const_iterator n = _node.begin();
       n != _node.end(); ++n)
    if (n->_edges[0].size() == 1 && n->_edges[1].size() == 0)
      external._edges[1].push_back(n->_edges[0].back());
    else if (n->_edges[0].size() == 0 && n->_edges[1].size() == 1)
      external._edges[0].push_back(n->_edges[1].back());

  return external;
}

/******************************************************************************
 *                                                                            *
 * n_cycles                                                                   *
 *                                                                            *
 ******************************************************************************/

int Topology::n_cycles() const
{
  int count = 0;

  // depth first search on each of the connected components

  vector<bool> used_edge(n_edges(), false);
  vector<bool> visited_node(n_nodes(), false);
  stack<int> nodes_to_visit;

  for (;;)
  {
    int starting_node;
    for (starting_node = 0; starting_node < n_nodes(); ++starting_node)
    {
      if (!visited_node[starting_node])
      {
        nodes_to_visit.push(starting_node);
        visited_node[starting_node] = true;
        break;
      }
    }
    if (starting_node == n_nodes())
      break;

    while (!nodes_to_visit.empty())
    {
      int n = nodes_to_visit.top();
      nodes_to_visit.pop();

      for (int d = 0; d < 2; ++d)
        for (vector<int>::const_iterator e = _node[n]._edges[d].begin();
             e != _node[n]._edges[d].end(); ++e)
        {
          if (used_edge[*e])
            continue;
          else
            used_edge[*e] = true;

          int next = _edge[*e]._node[d];
          if (visited_node[next])
            ++count;
          else
          {
            nodes_to_visit.push(next);
            visited_node[next] = true;
          }
        }
    }
  }

  return count;
}

/******************************************************************************
 *                                                                            *
 * is_one_particle_irreducible                                                *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The algorithm is a straightforward application of the definition, i.e.
 * the topology is one particle irreducible iff cutting any of its internal
 * edges doesn't lead to a disconnected graph.
 *
 */

bool Topology::is_one_particle_irreducible() const
{
  // check for connectedness after cutting an edge. edges that belong to a
  // multiple edge set or a loop are not cut

  bool cut_at_least_once = false;
  vector<bool> forbidden_edge(n_edges(), false);

  for (int e = 0; e < n_edges(); ++e)
  {
    if (!is_external_edge(e) &&
        adjacency(_edge[e]._node[0], _edge[e]._node[1]) == 1)
    {
      cut_at_least_once = true;
      forbidden_edge[e] = true;

      if (subgraph(0, forbidden_edge).size() != n_nodes())
        return false;

      forbidden_edge[e] = false;
    }
  }

  // it may happen that no edge will qualify to be cut, then run connected()

  if (!cut_at_least_once)
    return is_connected();
  else
    return true;
}

vector<int>
Topology::subgraph(int starting_node, const vector<bool> &forbidden_edge) const
{
  vector<int> g(1, starting_node);

  vector<bool> visited_node(n_nodes(), false);
  stack<int> nodes_to_visit;

  nodes_to_visit.push(starting_node);
  visited_node[starting_node] = true;

  while (!nodes_to_visit.empty())
  {
    int n = nodes_to_visit.top();
    nodes_to_visit.pop();

    for (int d = 0; d < 2; ++d)
      for (vector<int>::const_iterator e = _node[n]._edges[d].begin();
           e != _node[n]._edges[d].end(); ++e)
      {
        if (forbidden_edge[*e])
          continue;

        int next = _edge[*e]._node[d];
        if (!visited_node[next])
        {
          nodes_to_visit.push(next);
          visited_node[next] = true;
          g.push_back(next);
        }
      }
  }

  sort(g.begin(), g.end());
  return g;
}

/******************************************************************************
 *                                                                            *
 * biconnected_components                                                     *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Implementation of the recursive biconnected component algorithm from
 *
 * K. Melhorn, Algorithms and Data Structures,
 *
 */

vector<TopologyComponent>
Topology::biconnected_components() const
{
  vector<TopologyComponent> components;
  vector<int> external = external_nodes();
  BCCStatus status(n_nodes(), n_edges());

  for (int n = 0; n < n_nodes(); ++n)
    if (!status._visited_node[n])
    {
      bcc(status, n);

      for (vector<TopologyComponent>::iterator c = status._complete.begin();
           c != status._complete.end(); ++c)
      {
        for (vector<int>::iterator n = c->_nodes.begin();
             n != c->_nodes.end(); ++n)
          if (status._articulation_point[*n])
            c->_articulation_points.push_back(*n);

        if (!c->_vacuum)
        {
          c->_vacuum = true;

          if (status._n_external_nodes > 1)
          {
            vector<bool> forbidden_edge(n_edges(), false);
            for (vector<int>::iterator e = c->_edges.begin();
                 e != c->_edges.end(); ++e)
              forbidden_edge[*e] = true;

            for (vector<int>::iterator n = c->_articulation_points.begin();
                 n != c->_articulation_points.end(); ++n)
            {
              vector<int> seen = subgraph(*n, forbidden_edge);

              vector<int> intersection;

              set_intersection(seen.begin(), seen.end(),
                               external.begin(), external.end(),
                               back_inserter(intersection));

              if (intersection.size() != status._n_external_nodes &&
                  intersection.size() != 0)
              {
                c->_vacuum = false;
                break;
              }
            }
          }
        }

        components.push_back(*c);
      }

      // different connectedness components are independent

      status._complete.clear();
      status._n_external_nodes = 0;
      status._count = 0;
    }

  return components;
}

/******************************************************************************
 *                                                                            *
 * bcc                                                                        *
 *                                                                            *
 ******************************************************************************/

void Topology::bcc(BCCStatus &status, int n) const
{
  status._incomplete.push(make_pair(n, true));
  if (is_external_node(n))
    ++status._n_external_nodes;
  status._index[n] = status._count++;
  status._min_index[n] = status._index[n];
  status._visited_node[n] = true;

  for (int d = 0; d < 2; ++d)
    for (vector<int>::const_iterator e = _node[n]._edges[d].begin();
         e != _node[n]._edges[d].end(); ++e)
    {
      const int next = _edge[*e]._node[d];

      if (!status._used_edge[*e])
      {
        status._used_edge[*e] = true;

        if (next != n)
          status._incomplete.push(make_pair(*e, false));
        else
        {
          TopologyComponent component;
          component._nodes.push_back(n);
          component._edges.push_back(*e);
          component._vacuum = true;
          status._articulation_point[n] = true;
          ;
          status._complete.push_back(component);
          continue;
        }
      }

      if (status._visited_node[next])
        status._min_index[n] = min(status._min_index[n], status._index[next]);
      else
      {
        status._parent[next] = n;

        bcc(status, next);

        status._min_index[n] =
            min(status._min_index[n],
                status._min_index[next]);
      }
    }

  if ((status._index[n] >= 1) &&
      (status._min_index[n] == status._index[status._parent[n]]))
  {
    TopologyComponent component;

    pair<int, bool> object;
    do
    {
      object = status._incomplete.top();
      status._incomplete.pop();

      if (object.second)
        component._nodes.push_back(object.first);
      else
        component._edges.push_back(object.first);
    } while (!object.second || object.first != n);

    object = status._incomplete.top();
    status._incomplete.pop();

    const int parent = status._parent[n];
    component._nodes.push_back(parent);
    component._edges.push_back(object.first);

    // it cannot be decided whether the parent is an articulation point if
    // we are in the first component.
    if (status._index[n] > 1)
      status._articulation_point[parent] = true;
    ;
    status._complete.push_back(component);
  }
}

/******************************************************************************
 *                                                                            *
 * independent_external_node_permutations                                     *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Two external node permutations are independent when they cannot be
 * transformed into each other through a topology's symmetry transformation.
 * The current algorithm is ineffective for large n and should be improved in
 * future. It consists of filtering all permutations of the external nodes
 * through the symmetry transformations, and rejecting those that are not
 * lexicographically minimal.
 *
 */

vector<vector<int>>
Topology::independent_external_node_permutations() const
{
  vector<vector<int>> permutations;

  vector<int> external = external_nodes();
  const vector<vector<int>> group = node_symmetry_group();
  do
  {
    for (vector<vector<int>>::const_iterator permutation = group.begin();
         permutation != group.end(); ++permutation)
    {
      for (vector<int>::const_iterator n = external.begin();
           n != external.end(); ++n)
      {
        if (*n < (*permutation)[*n])
          break;
        if (*n > (*permutation)[*n])
          goto NotMinimal;
      }
    }

    permutations.push_back(external);

  NotMinimal:;
  } while (next_permutation(external.begin(), external.end()));

  return permutations;
}

/******************************************************************************
 *                                                                            *
 * symmetry_group                                                             *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The algorithm used is a linear search among the allowed permutations of the
 * search tree (see next_leaf()). Just as node_labelling(), it could be
 * made faster.
 *
 */

vector<vector<int>>
Topology::symmetry_group(const list<set<int>> &partition) const
{
  vector<int> identity;
  for (int n = 0; n < n_nodes(); ++n)
    identity.push_back(n);

  vector<vector<int>> group(1, identity);

  // if the partition is discrete, then the symmetry group is trivial

  if (partition.size() == n_nodes())
    return group;

  // if the initial partition is not discrete, then perform a search
  // among the allowed permutations, for all those that leave the adjacency
  // matrix invariant

  stack<PartitionStatus> partition_stack;
  init_partition_stack(partition_stack, partition);

  vector<int> node_index(n_nodes());
  const vector<int> reference_permutation = next_leaf(partition_stack);
  for (int n = 0; n < n_nodes(); ++n)
    node_index[reference_permutation[n]] = n;

  vector<int> permutation;
  while (!(permutation = next_leaf(partition_stack)).empty())
    if (compare_adjacency_matrix(permutation, reference_permutation) == 0)
    {
      vector<int> symmetry_permutation;
      for (int n = 0; n < n_nodes(); ++n)
        symmetry_permutation.push_back(permutation[node_index[n]]);

      group.push_back(symmetry_permutation);
    }

  return group;
}

/******************************************************************************
 *                                                                            *
 * fix_external_nodes                                                         *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * If external nodes are fixed, then node_labelling() and compare() assume that
 * their permutations are not allowed. Obviously, changing the bahaviour of
 * external nodes invalidates node_labelling().
 *
 */

void Topology::fix_external_nodes() const
{
  if (!_fixed_external_nodes)
  {
    _node_labelling.clear();
    _fixed_external_nodes = true;
  }
}

/******************************************************************************
 *                                                                            *
 * free_external_nodes                                                        *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Opposite of fix_external_nodes().
 *
 */

void Topology::free_external_nodes() const
{
  if (_fixed_external_nodes)
  {
    _node_labelling.clear();
    _fixed_external_nodes = false;
  }
}

/******************************************************************************
 *                                                                            *
 * node_labelling                                                             *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The node labelling gives a unique representation of the adjacency matrix with
 * respect to permutations of the nodes. External nodes may or may not be
 * permuted according to fixed_external_nodes().
 * 
 * The algorithm consists in searching for the lexicographically minimal
 * adjacency matrix obtained by permutations given by the leaves of the search
 * tree (see next_leaf()). Less comparisons would result in the full algorithm
 * of
 *
 * B. D. McKay, Practical Graph Isomorphism,
 * Congressus Numerantium 30 (1981) 45.
 * 
 */

vector<int>
Topology::node_labelling() const
{
  // return if the node labelling has already been determined

  if (!_node_labelling.empty())
    return _node_labelling;

  list<set<int>> partition;

  if (_fixed_external_nodes)
    partition = fixed_equitable_partition();
  else
    partition = free_equitable_partition();

  // if the equitable partition is discrete, then the node permutation
  // is as in the partition

  if (partition.size() == n_nodes())
  {
    for (list<set<int>>::iterator cell = partition.begin();
         cell != partition.end(); ++cell)
      _node_labelling.push_back(*cell->begin());

    return _node_labelling;
  }

  // if the equitable partition is not discrete, then perform a search among
  // the allowed permutations to find the permutation minimizing the
  // adjacency matrix

  stack<PartitionStatus> partition_stack;
  init_partition_stack(partition_stack, partition);

  _node_labelling = next_leaf(partition_stack);
  for (auto i : _node_labelling)cout << i << " ";
  cout << endl;
  //return _node_labelling;
  vector<int> permutation;
  while (!(permutation = next_leaf(partition_stack)).empty()){
    for(auto i:permutation)cout<<i<<" ";
    cout<<endl;	
  }
    //if (compare_adjacency_matrix(permutation, _node_labelling) < 0)
    //  _node_labelling = permutation;

  return _node_labelling;
}

/******************************************************************************
 *                                                                            *
 * compare                                                                    *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * Lexicographical comparison of adjacency matrices in the node_labelling().
 * The comparison concerns only the topological structure of the object and has
 * nothing to do with other attributes, like momenta. If the external nodes
 * are fixed then their permutations are not allowed in the determination of
 * the labelling, which thus influences the comparison.
 *
 */

int Topology::compare(const Topology &t) const
{
  if (n_nodes() != t.n_nodes())
    return n_nodes() - t.n_nodes();

  const vector<int> label = node_labelling();
  const vector<int> t_label = t.node_labelling();

  for (int n1 = 0; n1 < n_nodes(); ++n1)
    for (int n2 = 0; n2 <= n1; ++n2)
    {
      const int first_adjacency = adjacency(label[n1], label[n2]);
      const int second_adjacency = t.adjacency(t_label[n1], t_label[n2]);

      if (first_adjacency != second_adjacency)
        return first_adjacency - second_adjacency;
    }

  return 0;
}

/******************************************************************************
 *                                                                            *
 * fixed_equitable_partition                                                  *
 *                                                                            *
 ******************************************************************************/

list<set<int>>
Topology::fixed_equitable_partition() const
{
  list<set<int>> partition;

  set<int> last_cell;
  for (int n = 0; n < n_nodes(); ++n)
    if (!is_external_node(n))
      last_cell.insert(n);
    else
    {
      set<int> singleton_cell;
      singleton_cell.insert(n);
      partition.push_back(singleton_cell);
    }
  partition.push_back(last_cell);

  list<set<int>> active(partition);
  refine(partition, active);

  return partition;
}

/******************************************************************************
 *                                                                            *
 * free_equitable_partition                                                   *
 *                                                                            *
 ******************************************************************************/

list<set<int>>
Topology::free_equitable_partition() const
{
  set<int> all_nodes;
  for (int i = 0; i < n_nodes(); ++i)
    all_nodes.insert(i);

  list<set<int>> partition, active;
  partition.push_back(all_nodes);
  active.push_back(all_nodes);

  refine(partition, active);
  printSet(partition, "Equitable");
  return partition;
}

/******************************************************************************
 *                                                                            *
 * refine                                                                     *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The procedure refine is an implementation of the algorithm 2.5 of
 *
 * B. D. McKay, Practical Graph Isomorphism,
 * Congressus Numerantium 30 (1981) 45.
 * 
 * If partition and active are equal, then partition is modified to contain 
 * the coarsest equitable partition finer than the original. A second use is 
 * when partition was obtained from an equitable partition by splitting a cell 
 * into an element and the rest, and active contains this element. The result 
 * is then again an equitable partition. The procedure differs from the
 * original of B. D. McKay in the way the degree of a node is defined. Since
 * we treat here pseudographs, where loops and multiple edges are possible
 * the degree is a structure \see NodeDegree.
 *
 * \param partition the original partition
 * \param active list of active cells
 * \return true if the equitable partition is discrete, false otherwise
 *
 */

bool Topology::refine(list<set<int>> &partition, list<set<int>> &active) const
{
  // a discrete partition is equitable

  if (partition.size() == n_nodes())
    return true;

  // the main loops
  //std::cout << "New refinement\n";
  //printSet(partition, "Partition START");
  //printSet(active, "Active START");
  for (list<set<int>>::iterator active_cell = active.begin();
       active_cell != active.end(); ++active_cell)
  {
    for (list<set<int>>::iterator partition_cell = partition.begin();
         partition_cell != partition.end(); ++partition_cell)
    {

	//printSet(partition, "Partition ");
	//	printSet(active, "Active");
      // cells containing a single element can't be further divided.

      if (partition_cell->size() == 1)
        continue;

      // the cell is divided according to the scheme defined by the
      // NodeDegree structure.

      map<NodeDegree, set<int>> divided_partition_cell;
      for (set<int>::const_iterator partition_node = partition_cell->begin();
           partition_node != partition_cell->end(); ++partition_node)
      {
        NodeDegree d;
        for (set<int>::const_iterator active_node = active_cell->begin();
             active_node != active_cell->end(); ++active_node)
        {
          int n = adjacency(*partition_node, *active_node);
          if (*partition_node == *active_node)
            n = -n;
          d.insert_adjacency(n);
        }
        divided_partition_cell[d].insert(*partition_node);
      }

	  for (auto i : divided_partition_cell) {
		 // std::cout << "Key: " << i.first._total_degree << " ";
	//	  for (auto e : i.second) {
//			  std::cout << e << " ";
//		  }
//		  std::cout << std::endl;
	  }
      // in case the cell was not divided, there is nothing to update.

      if (divided_partition_cell.size() == 1)
        continue;

      // to perform the updates, proceed as in the original algorithm,
      // that is first find the largest subcell of the smallest degree.
      // the reason is most probably performance, since that will lead to
      // the highest number of subdivisions in the next pass.

      int max_subcell_size = 0;
      map<NodeDegree, set<int>>::const_iterator subcell, max_subcell;
      for (subcell = divided_partition_cell.end();
           subcell-- != divided_partition_cell.begin();)
      {
        int subcell_size = subcell->second.size();

        if (subcell_size > max_subcell_size)
        {
          max_subcell_size = subcell_size;
          max_subcell = subcell;
        }
      }
	  //for (auto i : max_subcell->second)cout << i<<" ";
	  //std::cout << std::endl;
      // if there is a copy of the partition cell in active, then it gets
      // replaced by the largest subcell. remark that such a copy will
      // always be there in the case partition and active are equal at the
      // beginnig. however, in the second use of the procedure, there will
      // be no match at least in the first pass.

      replace(active_cell, active.end(),
              *partition_cell, max_subcell->second);

      // the rest of the subcells are added at the end of active. on the
      // contrary, all of the cells replace the partition cell that was
      // divided. remark that we have to do some iterator magic so that
      // at the end partition cell is the last subcell in the partition.

      list<set<int>>::iterator last_subcell;
      for (subcell = divided_partition_cell.begin();
           subcell != divided_partition_cell.end(); ++subcell)
      {
        if (subcell != max_subcell)
          active.push_back(subcell->second);

        last_subcell = partition.insert(partition_cell, subcell->second);
      }
      partition.erase(partition_cell);
      partition_cell = last_subcell;
    }

    // check if the partition is discrete.

    if (partition.size() == n_nodes())
      return true;
  }

  return false;
}

/******************************************************************************
 *                                                                            *
 * init_partition_stack                                                       *
 *                                                                            *
 ******************************************************************************/

void Topology::init_partition_stack(stack<PartitionStatus> &partition_stack,
                                    const list<set<int>> &partition) const
{
  partition_stack.push(PartitionStatus());
  PartitionStatus &root = partition_stack.top();

  root._partition = partition;
  for (root._cell = root._partition.begin(); root._cell->size() == 1; ++root._cell)
    ;
  root._node = root._cell->begin();
}

/******************************************************************************
 *                                                                            *
 * next_leaf                                                                  *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The procedure next_leaf is an implementation of the algorithm 2.12 of
 *
 * B. D. McKay, Practical Graph Isomorphism,
 * Congressus Numerantium 30 (1981) 45.
 * 
 * It generates the leaves of the search tree. Contrary to 2.12, at any
 * moment only a branch of the tree is in memory in form of the partition_stack.
 * Remark also that the target cell at any level of the tree is chosen as in
 * the new version of nauty, that is as the first non singleton cell and
 * not as in the original algorithm, as the first non singleton cell of the
 * smallest size. There is also a misprint in the article, since as is the
 * algorithm 2.12 would not generate the whole tree, since it only takes the
 * first node in the target cell.
 *
 * \param partition_stack the stack of partitions. it must be allocated by the
 *        user and contains on the first entry the equitable partition with a
 *        pointer to the first node in the first non singleton cell.
 * \return the next allowed permutation or an empty vector if no further
 *         permutations could be found.
 *
 */

vector<int>
Topology::next_leaf(stack<PartitionStatus> &partition_stack) const
{
  vector<int> permutation;

  while (!partition_stack.empty())
  {
    PartitionStatus &parent = partition_stack.top();

    if (parent._node == parent._cell->end())
    {
      partition_stack.pop();
      continue;
    }

    // allocate the child partition on the stack.

    partition_stack.push(PartitionStatus());
    PartitionStatus &child = partition_stack.top();

    // split the target cell by taking out one of its elements.

    set<int> singleton_cell;
    singleton_cell.insert(*parent._node);
    for (list<set<int>>::iterator parent_cell = parent._partition.begin();
         parent_cell != parent._partition.end(); ++parent_cell)
    {
      child._partition.push_back(*parent_cell);
      if (parent_cell == parent._cell)
      {
        child._partition.insert(--child._partition.end(), singleton_cell);
        child._partition.back().erase(*parent._node);
      }
    }

    // as the active partition of the refine procedure use the
    // partition containing only a singleton cell.

    list<set<int>> singleton_partition(1, singleton_cell);

    // this is necessary, since in the next pass we want to split on
    // the next element.

    ++parent._node;

    // refine to get the equitable partition. the data on the stack is
    // completed only if the partition is not discrete.

    if (refine(child._partition, singleton_partition))
    {
      for (list<set<int>>::iterator child_cell = child._partition.begin();
           child_cell != child._partition.end(); ++child_cell)
        permutation.push_back(*child_cell->begin());
      partition_stack.pop();
      break;
    }
    else
    {
      for (child._cell = child._partition.begin();
           child._cell->size() == 1; ++child._cell)
        ;

      child._node = child._cell->begin();
    }
  }

  return permutation;
}

/******************************************************************************
 *                                                                            *
 * compare_adjacency_matrix                                                   *
 *                                                                            *
 ******************************************************************************/

int Topology::compare_adjacency_matrix(const vector<int> &first_permutation,
                                       const vector<int> &second_permutation) const
{
  for (int n1 = 0; n1 < n_nodes(); ++n1)
    for (int n2 = 0; n2 <= n1; ++n2)
    {
      const int a1 = adjacency(first_permutation[n1], first_permutation[n2]);
      const int a2 = adjacency(second_permutation[n1], second_permutation[n2]);

      if (a1 != a2)
        return a1 - a2;
    }

  return 0;
}

/******************************************************************************
 *                                                                            *
 * print                                                                      *
 *                                                                            *
 ******************************************************************************/

void Topology::print(ostream &output) const
{
  output << "Adjacency matrix:\n\n";
  print_adjacency_matrix(output);

  output << "\nAdjacency list:\n\n";
  print_adjacency_list(output);

  output << "\nEdge list:\n\n";
  print_edge_list(output);
}

/******************************************************************************
 *                                                                            *
 * print_edge_list                                                            *
 *                                                                            *
 ******************************************************************************/

void Topology::print_edge_list(ostream &output) const
{
  for (int e = 0; e < n_edges(); ++e)
  {
    output << e << ": " << _edge[e]._node[0] << " -> " << _edge[e]._node[1];

    if (!_momentum_basis.empty())
    {
      output << ", ";
      print_momentum(output, e);
    }
    output << "\n";
  }
}

/******************************************************************************
 *                                                                            *
 * postscript_print                                                           *
 *                                                                            *
 ******************************************************************************/

/**
 *
 * The topology is drawn using neato from the graphviz package, which can be
 * found at
 * 
 * http://www.research.att.com/sw/tools/graphviz/
 *
 * Currently the package doesn't support multiple edges, therefore every 
 * multiple edge is drawn with a bold line and has a label with the number of
 * edges.
 *
 */

void Topology::postscript_print(const string &file_name) const
{
  ostringstream name;
  name << "output.tmp";

  string temp_file_name(name.str());
  ofstream temp_file(temp_file_name.c_str());

  // the header defining the basic graph properties. currently nodes are
  // numbered

  //temp_file << "Graph G {\n"
  //    << "node [shape=point, style=filled, color=black, label=\"\"];\n";

  temp_file << "Graph G {\n"
            << "node [height=.1, width=.1, style=filled];\n";

  // loop over the nodes

  for (int n1 = 0; n1 < n_nodes(); ++n1)
    for (int n2 = 0; n2 <= n1; ++n2)
    {
      int value = adjacency(n1, n2);
      if (value == 0)
        continue;

      if (n1 == n2)
        value /= 2;

      temp_file << n1 << " -- " << n2;

      if (value > 1)
        temp_file << " [style=bold, label=" << value << "]";

      temp_file << ";\n";
    }

  // trailer

  temp_file << "}\n";

  temp_file.close();

  // prepare the command to run neato

  ostringstream command;
  command << "neato -Tps " << temp_file_name << " -o " << file_name;

  system(command.str().c_str());

  //unlink(temp_file_name.c_str());
}

/******************************************************************************
 *                                                                            *
 * print_adjacency_matrix                                                     *
 *                                                                            *
 ******************************************************************************/

void Topology::print_adjacency_matrix(ostream &output) const
{
  for (vector<vector<int>>::const_iterator i = _adjacency_matrix.begin();
       i != _adjacency_matrix.end(); ++i)
  {
    copy(i->begin(), i->end(), ostream_iterator<int>(output, " "));
    output << "\n";
  }
}

/******************************************************************************
 *                                                                            *
 * print_adjacency_list                                                       *
 *                                                                            *
 ******************************************************************************/

void Topology::print_adjacency_list(ostream &output) const
{
  for (int n = 0; n < n_nodes(); ++n)
  {
    vector<int> adjacent_nodes;

    for (int d = 0; d < 2; ++d)
      for (vector<int>::const_iterator e = _node[n]._edges[d].begin();
           e != _node[n]._edges[d].end(); ++e)
        adjacent_nodes.push_back(_edge[*e]._node[d]);

    sort(adjacent_nodes.begin(), adjacent_nodes.end());

    output << n << ": ";
    copy(adjacent_nodes.begin(),
         adjacent_nodes.end(),
         ostream_iterator<int>(output, " "));
    output << "\n";
  }
}

/******************************************************************************
 *                                                                            *
 * print_momentum                                                             *
 *                                                                            *
 ******************************************************************************/

void Topology::print_momentum(ostream &output, int e) const
{
  if (_momentum_basis.empty())
    return;

  bool first_output = true;
  vector<string>::const_iterator s = _momentum_basis.begin();
  for (vector<int>::const_iterator p = _momentum[e].begin();
       p != _momentum[e].end(); ++p, ++s)
  {
    switch (*p)
    {
    case 1:
      if (!first_output)
        output << "+";
      else
        first_output = false;
      output << *s;
      break;
    case -1:
      if (first_output)
        first_output = false;
      output << "-" << *s;
      break;
    case 0:
      break;
    default:
      cerr << "Topology::print_momentum: incorrect momentum "
           << _momentum[e] << endl;

      exit(1);
    }
  }
  if (first_output)
    output << 0;
}

void Topology::print_momentum(ostream &output, int e, const string &index) const
{
  if (_momentum_basis.empty())
    return;

  bool first_output = true;
  vector<string>::const_iterator s = _momentum_basis.begin();
  for (vector<int>::const_iterator p = _momentum[e].begin();
       p != _momentum[e].end(); ++p, ++s)
  {
    switch (*p)
    {
    case 1:
      if (!first_output)
        output << "+";
      else
        first_output = false;
      output << *s << "(" << index << ")";
      break;
    case -1:
      if (first_output)
        first_output = false;
      output << "-" << *s << "(" << index << ")";
      break;
    case 0:
      break;
    default:
      cerr << "Topology::print_momentum: incorrect momentum "
           << _momentum[e] << endl;

      exit(1);
    }
  }
  if (first_output)
    output << 0;
}

void Topology::print_momentum(ostream &output, const vector<int> &momentum) const
{
  if (_momentum_basis.empty())
    return;

  bool first_output = true;
  vector<string>::const_iterator s = _momentum_basis.begin();
  for (vector<int>::const_iterator p = momentum.begin();
       p != momentum.end() || s != _momentum_basis.end(); ++p, ++s)
  {
    switch (*p)
    {
    case 1:
      if (!first_output)
        output << "+";
      else
        first_output = false;
      output << *s;
      break;
    case -1:
      if (first_output)
        first_output = false;
      output << "-" << *s;
      break;
    case 0:
      break;
    default:
      cerr << "Topology::print_momentum: incorrect momentum "
           << momentum << endl;

      exit(1);
    }
  }
  if (first_output)
    output << 0;
}

void Topology::print_momentum(ostream &output,
                              const vector<int> &momentum,
                              const string &index) const
{
  if (_momentum_basis.empty())
    return;

  bool first_output = true;
  vector<string>::const_iterator s = _momentum_basis.begin();
  for (vector<int>::const_iterator p = momentum.begin();
       p != momentum.end() || s != _momentum_basis.end(); ++p, ++s)
  {
    switch (*p)
    {
    case 1:
      if (!first_output)
        output << "+";
      else
        first_output = false;
      output << *s << "(" << index << ")";
      break;
    case -1:
      if (first_output)
        first_output = false;
      output << "-" << *s << "(" << index << ")";
      break;
    case 0:
      break;
    default:
      cerr << "Topology::print_momentum: incorrect momentum "
           << momentum << endl;

      exit(1);
    }
  }
  if (first_output)
    output << 0;
}
