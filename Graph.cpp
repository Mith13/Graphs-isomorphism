#include "Graph.h"
#include <map>
#include <algorithm>
#include <iterator>
#include <string>
using std::map;

// Prints an entire parttion
inline
void printSet(partition set, std::string s = "")
{
	std::cout << s<<std::endl;
	for (auto s : set) {
		std::cout << "{";
		for (auto e : s) {
			std::cout << " " << e + 1;
		}
		std::cout << " }";
	}
	std::cout << std::endl;
}
// Prints only a set inside some partition 
inline
void printSet(set<int> set, std::string s = "")
{
	std::cout << s<<std::endl;
	for (auto s : set) {
			std::cout << " " << s + 1;
	}
	std::cout << std::endl;
}
inline
void printSet2(vector<int> set, std::string s = "")
{
	std::cout << s;
	for (auto s : set) {
			std::cout << " " << s + 1;
	}
	std::cout << std::endl;
}
// Copy ctor 
Graph::Graph(const Graph & g)
{
	m_adjacency_matrix=g.m_adjacency_matrix;
	m_adjacent_list=g.m_adjacent_list;
	m_canonical_label=g.m_canonical_label;
}
// Move ctor
Graph::Graph(const Graph && g)
{
	m_adjacency_matrix=std::move(g.m_adjacency_matrix);
	m_adjacent_list=std::move(g.m_adjacent_list);
	m_canonical_label=std::move(g.m_canonical_label);
}

// Inserts an edge into the graph
// Canonical label is always cleared after a new edge is added
void Graph::insertEdge(int a1, int a2)
{
	//If any vertex is outside for boundaries, abort
	if (a1 > getVertices() || a2 > getVertices() || a1 < 0 || a2 < 0) {
		std::cerr << "Problem";
		abort();
	}
	a1--;
	a2--;
	m_adjacent_list[a1].insert(a2);
	m_adjacent_list[a2].insert(a1);
	m_adjacency_matrix[a1][a2]++;
	m_adjacency_matrix[a2][a1]++;
	m_canonical_label.clear();

}

// Compares a graph with another
// Returns true if two graphs are isomorphic,
// otherwise false
bool Graph::isIsomoprhic(Graph& rhs){
	createCanonicalLabel();
	rhs.createCanonicalLabel();
	for (int i = 0; i < getVertices(); i++) {
		for (int j = 0; j <= i; j++) {
			//Do permutation of adjacency_matrix and then compare them
			int a1 = getAdjacency(m_canonical_label[i],m_canonical_label[j]);
			int a2 = rhs.getAdjacency(rhs.m_canonical_label[i],rhs.m_canonical_label[j]);
			if (a1 != a2)return false;
		}
	}
	return true;
}

// Creates canonical label according to McKay alghorithm
// see. Practical Graph Isomorphishm II. https://arxiv.org/abs/1301.1493
void Graph::createCanonicalLabel()
{
	//If canonical label is created skip the creation of new one
	if(m_canonical_label.empty()){
		vector<vector<int>> permutations=findIsomorphisms();
		vector<int> perm=permutations[0];
		for(auto it_perm=permutations.begin()+1;it_perm!=permutations.end();it_perm++){
			// canonical label is choosen according to the lowest adjacency matrix elements
			if(lower_permutation(*it_perm,perm)<0){
				perm=*it_perm;
			}
		}
		m_canonical_label=perm;
	}
	printSet2(m_canonical_label,"Canonical");
}

// Compares permuted elements of the adjacency matrix and returns the difference between the two
int Graph::lower_permutation(const vector<int>& perm1,const vector<int>& perm2)const 
{
	for(int i=0;i<getVertices();i++){
		for(int j=0;j<=i;j++){
			int a1=getAdjacency(perm1[i],perm1[j]);
			int a2=getAdjacency(perm2[i],perm2[j]);
			if(a1!=a2)return a1-a2;
		}
	}
	return 0;
}

// Returns an initial equitable splitting of graph vertices
partition Graph::equipartition() const
{
	// In this case final and working partition are the same
	set<int> all_nodes;
	for (int i = 0; i < getVertices(); i++) all_nodes.insert(i);
	partition colouring_out, colouring_wrk;
	colouring_out.push_back(all_nodes);
	colouring_wrk.push_back(all_nodes);
	refineColouring(colouring_out, colouring_wrk);
	printSet(colouring_out,"Equitable");
	//abort();
	return colouring_out;
}

// Returns all possible isomorphic permutations or vertices
vector<vector<int>> Graph::findIsomorphisms() const
{

	vector<vector<int>> permutations; // container for all permutations
	vector<int> permutation;
	// Do first splitting with equipartition()
	partition equitable_partition = equipartition(); 
	if (equitable_partition.size() == getVertices()) {
		for(const auto& set:equitable_partition)
			for(const auto& vertex:set)
			permutation.push_back(vertex);
		permutations.push_back(permutation);
		return permutations;
	}//if there are too many nodes, do iterative approach
	//else we can do it recursively.
	if (getVertices() > 2) {
		std::stack<std::shared_ptr<TreeNode>> tree; //stack for iterative solution
		auto root = std::make_shared<TreeNode>(equitable_partition);
		//Set up the root node
		// set m_it_target_set to first set with more than one vertex
		// and set m_it_target_vertex to first vertex of said set
		for (root->m_it_target_set = root->m_data.begin(); root->m_it_target_set->size() == 1;
			++root->m_it_target_set++) {}
		root->m_it_target_vertex = root->m_it_target_set->begin();
		tree.push(root);
		while(!(permutation=next_permutation(tree)).empty()){
			permutations.push_back(permutation);
			for(auto i:permutation)std::cout<<i<<" ";
			std::cout<<std::endl;
		}
		return permutations;
	}
	else {

	}
	return permutations;
}

// Finds one permutation at the time
// This is based on the Mckay tree search without automorphism pruning.
// Thus, the search for canonical label is much slower than in nauty alghorithm.
// In this case we keep only one branch of tree in memory. 
// Therefore, the need for two iterators in the TreeNode class.
vector<int> Graph::next_permutation(std::stack<std::shared_ptr<TreeNode>>& tree)const
{
	vector<int> permutation;
	while (!tree.empty()) {
		//std::cout << "Next while" << std::endl;
		//Get parent from stack
		auto parent = tree.top();
		// If we are at the end of target set, there will be any child pushed onto stack 
		// and we can pop out the parent
		if(parent->m_it_target_vertex==parent->m_it_target_set->end()){
				tree.pop();
				continue;
		}
			partition colouring_out; // colouring_out equals \pi partition in McKay algh.
			//for (auto it_colour = parent->m_data.begin(); it_colour!=parent->m_it_target_set; it_colour++) {
			//	colouring_out.push_back(*it_colour);
			//}

			// Copy all expect the target set to current \pi partition.
			// Later we will insert splitted target cell on the position of target set
			std::copy_if(parent->m_data.begin(), parent->m_data.end(),
				std::inserter(colouring_out, colouring_out.end()), [&](set<int> s)
                    {return s != *parent->m_it_target_set; });
                    
			// Split target set to singleton partition containing target vertex and rest of set
			set<int>tmp;
			if(parent->m_it_target_set->size()>1){            
				std::copy_if(parent->m_it_target_set->begin(),parent->m_it_target_set->end(),
					std::inserter(tmp, tmp.end()),[&](const int i)
                    	{return i !=*(parent->m_it_target_vertex);});
			}
            set<int> single_colour_set{*parent->m_it_target_vertex};
			partition single_colour_partition;
			single_colour_partition.push_back(single_colour_set);
			// Insert splitted target set into current partition
			colouring_out.insert(colouring_out.begin()+(parent->m_it_target_set-parent->m_data.begin()),{ {single_colour_set}, {tmp} });
			//printSet(colouring_out);
			
			//Move current partition to child member partition
			std::shared_ptr<TreeNode> child=std::make_shared<TreeNode>(std::move(colouring_out));
			
			// Incremenet the target vertex of parent, so in next run we extract a correct vertex 
			parent->m_it_target_vertex++;

			// Now do a refinement of partition
			// single_colour_partition contains target vertex with which help we introduce asymmetricity
			// during the refinement process
			if(refineColouring(child->m_data,single_colour_partition)){				
				for(auto colour_set=child->m_data.begin();colour_set!=child->m_data.end();colour_set++){
					for(auto vertex=colour_set->begin();vertex!=colour_set->end();++vertex){
						 permutation.push_back(*vertex);
					}
				}
				return permutation;
			}
			else{
				for (child->m_it_target_set = child->m_data.begin(); child->m_it_target_set->size() == 1;
					++child->m_it_target_set) {
				////		std::cout << "searchin for nonsingleton \n";
				//		printSet(*child->m_it_target_set, "");
				}
				child->m_it_target_vertex = child->m_it_target_set->begin();
				//std::cout << "Target set ";
				//printSet(*child->m_it_target_set, "");
				//std::cout << "Target vertex " << *child->m_it_target_vertex+1 <<"/////////////////////\n"<< std::endl;
				tree.push(child);
			}
			
		}
		return permutation;
}

// Prints an adjacency matrix
void Graph::printAdjacencyMatrix() const
{
	for (int i = 0; i < getVertices(); i++) {
		for (int j = 0; j < getVertices(); j++) {
			std::cout << m_adjacency_matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

// Refinement function
// Splits the partition according the colouring_wrk
// colouring_out equals to \alpha in McKay algh
// colouring_wkr equals to \pi
bool Graph::refineColouring(partition& colouring_out, partition& colouring_wrk) const
{
	//std::cout << "\nNew refinement\n";
	//printSet(colouring_out, "Partition START");
	//printSet(colouring_wrk, "Active START");

	// If there is n_vertices sets inside partition we have the final partitioning
	if (colouring_out.size() == getVertices()) return true;
	
	while(!colouring_wrk.empty()){
		
		//printSet(colouring_out, "Partition ");
		//printSet(colouring_wrk, "Active ");
		//printSet(spliting_set, "splitting set ");
		// According to McKay paper its best to chose single vertex set 
		/*set<int> spliting_set;
		if(colouring_wrk.size()==1){
			spliting_set=colouring_wrk[0];
			colouring_wrk.pop_back();
		}
		else {			
			partition::iterator it_set;
			for(it_set=colouring_wrk.begin();it_set->size()==1;it_set++){
				spliting_set.insert(it_set->begin(),it_set->end());
			}			
			colouring_wrk.erase(it_set);
		}*/
		// Old code		
		auto spliting_set=colouring_wrk.back();
		colouring_wrk.pop_back();


		for(int i=0;i<colouring_out.size();i++){
			auto& split_set=colouring_out[i];

			//if split set size is one, we do not have anything to 
			if (split_set.size() > 1) {
				//printSet(split_set, "Splitting this set");
				// map sets according their edge degree
				map<PseudoEdge, set<int>> quasicolouring;
				for (const auto& vertex_in_split : split_set) {		
					PseudoEdge degree;
					for (const auto& splitting_vertex : spliting_set) {
						degree +=(m_adjacency_matrix[vertex_in_split][splitting_vertex]);						
					}
					quasicolouring[degree].insert(vertex_in_split);
				}
				/*for (auto i : quasicolouring) {
					std::cout << "Key: [" << i.first << "] ";
					for (auto e : i.second) {
						std::cout << e+1<<" ";
					}
					std::cout << std::endl;
				}*/
				// Continue, if set was not splitted
				if (quasicolouring.size() == 1) continue;
				
				// Find the biggest set from all splitted sets
				auto biggest_set = std::max_element(
					quasicolouring.rbegin(), quasicolouring.rend(),
					[](auto& lhs, auto& rhs) {return lhs.second.size() < rhs.second.size();}
				);

				//now replace splitting_set with the biggest set from partitioned split_set
				std::replace(colouring_wrk.begin(), colouring_wrk.end(),
					split_set, biggest_set->second);
				//std::cout << "Bigges set is ";
				//printSet(biggest_set->second, "biggest set");
					
				partition tmp;
				for(auto it=quasicolouring.rbegin();it!=quasicolouring.rend();it++){
					tmp.push_back(it->second);
					if (biggest_set->second != it->second) {
						//add rest of parttioned split_set to splitting_sets
						colouring_wrk.push_back(it->second);
					}
				}
				
				// replace split_set with its partitioned sets
				// this invalidates a pointer to split_set, but we are at the end
				// of the loop anyway and the pointer will be rejuvinated.
				colouring_out.insert(colouring_out.begin()+i,tmp.begin(),tmp.end());
				i+=(tmp.size()-1);//update index
				colouring_out.erase(colouring_out.begin()+i+1);
			}
		}
		if (colouring_out.size() == getVertices()) {
			//std::cout << "Equitable" << std::endl;
			return true;
		}
	}
	return false;
}


Graph::~Graph()
{
}

//Custom comparator for PseudoEdge class
bool operator <(const PseudoEdge& lhs, const PseudoEdge& rhs) noexcept{
		if (lhs.m_nodes != rhs.m_nodes)return lhs.m_nodes<rhs.m_nodes;
		if (lhs.m_degree != rhs.m_degree) return lhs.m_degree<rhs.m_degree;
		return lhs.m_node_types<rhs.m_node_types;
}