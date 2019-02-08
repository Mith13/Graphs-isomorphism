#include "Graph.h"
#include <map>
#include <algorithm>
#include <iterator>
#include <string>
using std::map;

void printSet(vector<set<int>> set, std::string s = "")
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
void printSet(set<int> set, std::string s = "")
{
	std::cout << s;
	for (auto s : set) {
			std::cout << " " << s + 1;
	}
	std::cout << std::endl;
}
Graph::Graph(const Graph & g)
{
}

Graph::Graph(const Graph && g)
{
}

void Graph::insertEdge(int a1, int a2)
{
	if (a1 > getVertices() || a2 > getVertices() || a1 < 0 || a2 < 0) {
		std::cerr << "Problem";
		abort();
	}
	m_adjacent_list[a1].insert(a2);
	m_adjacent_list[a2].insert(a1);
	m_adjacency_matrix[a1][a2]++;
	m_adjacency_matrix[a2][a1]++;
	m_canonical_label.clear();

}
bool Graph::compare(Graph& rhs) const{
	return(lesser_permutation(m_canonical_label,rhs.m_canonical_label)==0);
}
void Graph::createCanonicalLabel()
{
	vector<vector<int>> permutations(std::move(findIsomorphisms()));
	vector<int> perm=permutations[0];
	for(auto it_perm=permutations.begin()+1;it_perm!=permutations.end();it_perm++){
		if(lesser_permutation(*it_perm,perm)<0){
			perm=*it_perm;
		}
	}
	m_canonical_label=perm;
}
int Graph::lesser_permutation(const vector<int>& perm1,const vector<int>& perm2)const 
{
	for(int i=0;i<perm1.size();i++){
		for(int j=0;i<=j;j++){
			int a1=m_adjacency_matrix[perm1[i]][perm1[j]];
			int a2=m_adjacency_matrix[perm2[i]][perm2[j]];
			if(a1!=a2)return a1-a2;
		}
	}
	return 0;
}
vector<set<int>> Graph::equipartition() const
{
	set<int> all_nodes;
	for (int i = 0; i < getVertices(); i++) all_nodes.insert(i);
	partition colouring_tmp, colouring_wrk;
	colouring_tmp.push_back(all_nodes);
	colouring_wrk.push_back(all_nodes);
	refineColouring(colouring_tmp, colouring_wrk);
	return colouring_tmp;
}

vector<vector<int>> Graph::findIsomorphisms() const
{
	partition equitable_partition = equipartition();
	vector<vector<int>> permutations;
	//if there are too many nodes, do iterative approach
	//else we can do it recursively.
	if (getVertices() > 6) {
		std::stack<std::shared_ptr<TreeNode>> tree;
		auto root = std::make_shared<TreeNode>(equitable_partition);
		root->m_it_target_set = root->m_data.begin();
		for (auto it_first_nosingle_colour = root->m_data.begin(); it_first_nosingle_colour->size() == 1;
			it_first_nosingle_colour++) {
         	    root->m_it_target_set = it_first_nosingle_colour;
		}
		root->m_it_target_vertex = root->m_it_target_set->begin();
		tree.push(root);
		vector<int> permutation;
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
vector<int> Graph::next_permutation(std::stack<std::shared_ptr<TreeNode>>& tree)const
{
	vector<int> permutation;
	while (!tree.empty()) {
		//std::cout << "Next while" << std::endl;
		auto parent = tree.top();
		if(parent->m_it_target_vertex==parent->m_it_target_set->end()){
				tree.pop();
				continue;
		}
			partition colouring_out;
			//for (auto it_colour = parent->m_data.begin(); it_colour!=parent->m_it_target_set; it_colour++) {
			//	colouring_out.push_back(*it_colour);
			//}

			std::copy_if(parent->m_data.begin(), parent->m_data.end(),
				std::inserter(colouring_out, colouring_out.end()), [&](set<int> s)
                    {return s != *parent->m_it_target_set; });
                    
            set<int>tmp;
			std::copy_if(parent->m_it_target_set->begin(),parent->m_it_target_set->end(),
				std::inserter(tmp, tmp.end()),[&](const int i)
                    {return i !=*(parent->m_it_target_vertex);});
            set<int> single_colour_set{*parent->m_it_target_vertex};
			vector<set<int>> single_colour_partition;
			single_colour_partition.push_back(single_colour_set);

						colouring_out.insert(colouring_out.begin()+(parent->m_it_target_set-parent->m_data.begin()),{ {single_colour_set}, {tmp} });
			std::shared_ptr<TreeNode> child=std::make_shared<TreeNode>(colouring_out);
			
			parent->m_it_target_vertex++;

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

void Graph::printAdjacencyMatrix() const
{
	for (int i = 0; i < getVertices(); i++) {
		for (int j = 0; j < getVertices(); j++) {
			std::cout << m_adjacency_matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}
bool Graph::refineColouring(partition& colouring_out, partition& colouring_wrk) const
{
	//std::cout << "\nNew refinement\n";
	//printSet(colouring_out, "Partition START");
	//printSet(colouring_wrk, "Active START");

	if (colouring_out.size() == getVertices()) return true;
	while(!colouring_wrk.empty()){
	//for (auto& spliting_set : colouring_wrk) {
		auto spliting_set=colouring_wrk.back();

		//printSet(colouring_out, "Partition ");
		//printSet(colouring_wrk, "Active ");
		//printSet(spliting_set, "splitting set ");
		colouring_wrk.pop_back();
		//for (vector<set<int>>::iterator split_set=colouring_out.begin();split_set!=colouring_out.end();split_set++) {
		for(int i=0;i<colouring_out.size();i++){
			auto& split_set=colouring_out[i];
			if (split_set.size() > 1) {
				//printSet(split_set, "Splitting this set");
				map<int, set<int>> quasicolouring;

				for (const auto& node_in_split : split_set) {		
					int degree=0;
					for (const auto& splitting_node : spliting_set) {
						degree += m_adjacency_matrix[node_in_split][splitting_node];						
						if(node_in_split==splitting_node)degree-= m_adjacency_matrix[node_in_split][splitting_node];
					}
					quasicolouring[degree].insert(node_in_split);
				}
				/*for (auto i : quasicolouring) {
					std::cout << "Key: [" << i.first << "] ";
					for (auto e : i.second) {
						std::cout << e+1<<" ";
					}
					std::cout << std::endl;
				}*/
				if (quasicolouring.size() == 1) continue;
				
				auto biggest_set = std::max_element(
					quasicolouring.begin(), quasicolouring.end(),
					[](auto& lhs, auto& rhs) {return lhs.second.size() < rhs.second.size();}
				);

				//now replace splitting_set with the biggest set from partitioned split_set
				std::replace(colouring_wrk.begin(), colouring_wrk.end(),
					split_set, biggest_set->second);
				//std::cout << "Bigges set is ";
				//printSet(biggest_set->second, "biggest set");
					
				partition tmp;
				for(auto it=quasicolouring.begin();it!=quasicolouring.end();it++){
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
			std::cout << "Equitable" << std::endl;
			return true;
		}
	}
	return false;
}


Graph::~Graph()
{
}
