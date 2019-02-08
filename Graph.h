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

class TreeNode
{
public:
	TreeNode(partition data) :m_data(data), m_child(),m_it_target_vertex(),m_it_target_set() {};
	TreeNode() = delete;
	partition m_data;
	partition::iterator m_it_target_set;
	set<int>::iterator m_it_target_vertex;
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

	bool compare(Graph& rhs);

	int getVertices() const { return m_adjacent_list.size(); }
	
	vector<int> m_canonical_label;
	~Graph();
private:
	bool refineColouring(partition& colouring_out, partition& colouring_wrk) const;
	
	int lesser_permutation(const vector<int>& perm1,const vector<int>& perm2) const;

	vector<vector<int>> findIsomorphisms() const;
	vector<int> next_permutation(std::stack<std::shared_ptr<TreeNode>>& tree) const;
	vector<set<int>> equipartition() const;
	
	vector<vector<int>> m_adjacency_matrix;
	vector<set<int>> m_adjacent_list;
};
