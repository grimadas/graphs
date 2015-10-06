/*
 	Main application
	Author : Bulat, 2015

*/
//#include "graph.cuh"
#include "apsp.cuh"
#include "graph.cuh"

/*

*/


/************************************************
* 												*
*  CUDA graph representation					*
*  author: Bulat 								*
*  graph.h 										*
*************************************************
*		la - look ahed parametr
*		L - opacity threshold (1, 2) - small value
*		num_vertex - number of vetexes in a graph
*		n(n-1)/2 combinations pairs
*		2 variants : edge removal and remove (adding)
*		T - types calc distance to
*		0) Reading initial graph
*			a) In Edge format
*			b) Convert to Adj list (CSR)
*		1) Distance Matrix calculation
*		 - How to store ?
*		  * n arrays for different level opacity
*		  * special CSR format with different levels - 2 1D flatten matrix, since we know L and num_vertex
*			Example:
*				L1			L2 (same structure)
*				| 1| 2| 4|
*				| | | | | | |  Memory: O(L*|V| + L* 2*|E|) (what about dublication? ) two choices
*		  * COO format (same situation ?)
*							   Memory: O(2*|E|)

*/



int main()
{

	//Graph graph;
	init_test_graph(); // Reading graph from the file in COO format
	print_coo_graph();
	convert_to_CSR();
//	graph.print_csr_graph();

	UINT wTimerRes = 0;
	bool init = InitMMTimer(wTimerRes);
	DWORD startTime = timeGetTime();

	//test_funct<<<1, 1>>>(thrust::raw_pointer_cast(graph.full_vertex_array.data()), 4);

	// BFS
	//graph.single_bfs(2);
//	graph.form_full_level_graph();
//	graph.print_csr_graph();
//	graph.calc_L_opacity();
//	graph.print_opacity_matrix();

	unsigned int endTime = timeGetTime();
	unsigned int gpu_time = unsigned int(endTime - startTime);
	printf("GPU Timing(including all device-host, host-device copies, device allocations and freeing of device memory): %dms\n\n", gpu_time);
	DestroyMMTimer(wTimerRes, init);



	return 0;
}
