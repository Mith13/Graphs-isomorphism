#pragma once
#include <vector>
#include <set>
#include <stack>
#include <iostream>
#include <memory>
#include <iostream>

using std::vector;
using std::set;
using partition=vector<set<int>>;

// Class for pseudographs
// As we can have multiple edges from same node, we need to distinquish 
// two edges from same node and two edges from different edges
class PseudoEdge
{
public:
	PseudoEdge() :m_degree(0), m_nodes(0),m_node_types() {};
	int m_degree;
	int m_nodes;
	set<int> m_node_types;
	PseudoEdge& operator +=(const int& e) {
		m_degree += e;
		m_nodes++;
		if (e != 0) m_node_types.insert(e);
		return *this;
	}	
};

// Node class for search tree
class TreeNode
{
public:
	TreeNode(partition& data) :m_data(data), m_child(), m_it_target_vertex(), m_it_target_set() {};
	TreeNode(partition&& data) :m_data(std::forward<partition>(data)), m_child(), m_it_target_vertex(), m_it_target_set() {};
	TreeNode() = delete;
	partition m_data;
	set<int>::iterator m_it_target_vertex;
	partition::iterator m_it_target_set;
	vector<std::shared_ptr<TreeNode>>m_child;
};

class Graph
{
public:
	Graph() = delete;
	Graph(int n_nodes) :m_adjacent_list(n_nodes), m_adjacency_matrix(n_nodes, vector<int>(n_nodes)),m_canonical_label(n_nodes){};
	Graph(const Graph& g);
	Graph(const Graph&& g);

	void insertEdge(int a1, int a2);
	void createCanonicalLabel();
	void printAdjacencyMatrix() const;

	bool isIsomoprhic(Graph& rhs);

	int getVertices() const { return m_adjacent_list.size(); }
	int getAdjacency(int i, int j) const { return m_adjacency_matrix[i][j];}

	vector<int> m_canonical_label;
	~Graph();
private:
	bool refineColouring(partition& colouring_out, partition& colouring_wrk) const;
	
	int lower_permutation(const vector<int>& perm1,const vector<int>& perm2) const;

	vector<vector<int>> findIsomorphisms() const;
	vector<int> next_permutation(std::stack<std::shared_ptr<TreeNode>>& tree) const;
	vector<set<int>> equipartition() const;
	
	vector<vector<int>> m_adjacency_matrix;
	vector<set<int>> m_adjacent_list;
};
