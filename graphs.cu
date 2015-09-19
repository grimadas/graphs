/*
 	Main application
	Author : Bulat, 2015

*/
#include "graph.h"
#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

int main()
{
	Graph graph;
	graph.init_test_graph(); // Reading graph from the file in COO format
	graph.print_coo_graph();
	graph.convert_to_CSR();
	graph.print_csr_graph();
	
	return 0;
}
