/*
 	Main application
	Author : Bulat, 2015

*/
#include "graph.h"


int main()
{
	Graph graph;
	graph.init_test_graph(); // Reading graph from the file in COO format
	graph.print_coo_graph();
	graph.convert_to_CSR();
	graph.print_csr_graph();	

	// BFS 
	//graph.single_bfs(2);
	graph.form_full_level_graph();
	//graph.print_csr_graph();
	graph.calc_L_opacity();
	graph.print_opacity_matrix();

	
	return 0;
}
