/*
 	Main application
	Author : Bulat, 2015

*/
//#include "graph.cuh"
#include "apsp.cuh"
#include "graph.cuh"
#include "headers.h"




__global__  void expander(device_ptr<vertex> previous, device_ptr<vertex> current,
	device_ptr<vertex> current_vertex, device_ptr<vertex> temp_from, device_ptr<vertex> temp_to,
	device_ptr<vertex> full_vertex_array, device_ptr<vertex> full_edge_array,	int number_of_vertex
	)
{

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	int start = 0;
	if (idx != 0)
	{
		start = full_vertex_array[current_vertex[idx - 1]];
	}
	int end = full_vertex_array[current_vertex[idx]];

	printf("Expander ok 0 \n");

 	current =
	thrust::copy_if(thrust::device,
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_from[idx])),
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_to[idx])),
		current,
		unique_edge());

		/*
	current_vertex_ptr =
	thrust::copy(thrust::device,
		thrust::make_constant_iterator(current_vertex),
		thrust::make_constant_iterator(current_vertex) + thrust::distance(previous, current),
		current_vertex_ptr);
 	previos = current;
	printf("Expander ok 1 \n");
	*/

/*

	printf("Expander ok 2 \n");
	current[current_vertex[idx]] = thrust::remove(thrust::device, previous[current_vertex[idx]], current[current_vertex[idx]], current_vertex[idx] + 1);
	printf("Expander ok 3 \n");

	thrust::sort(thrust::device, previous[current_vertex[idx]], current[current_vertex[idx]]);
	printf("Expander ok 4 \n");

	current[current_vertex[idx]] = thrust::unique(thrust::device, previous[current_vertex[idx]], current[current_vertex[idx]]);
	printf("Expander ok 5 \n");

	// Depends on degree !
	for (auto j = full_edge_array + start; j != full_edge_array + end; j++)
	{
		current[current_vertex[idx]] = thrust::remove(thrust::device, previous[current_vertex[idx]], current[current_vertex[idx]], *j);
	}
	printf("Expander ok 6 \n");

	full_vertex_array[current_vertex[idx] + number_of_vertex] = thrust::distance(previous[current_vertex[idx]], current[current_vertex[idx]]);
	previous[current_vertex[idx]] = current[current_vertex[idx]];
	printf("Expander ok 7 \n");


*/
}



/*
*	By finding shortest paths, form to L_VALUE level
*/
void form_full_level_graph(Graph graph)
{

	for (int current_level = 1; current_level < graph.L_VALUE; current_level++)
	{
	int starting_point = 0;
	int ending_point = graph.full_vertex_array[current_level * graph.number_of_vertex - 1];
	if (current_level != 1)
	{
		starting_point = graph.full_vertex_array[(current_level - 1) * graph.number_of_vertex - 1]; // previous last element
	}
	vertex number_edges_to_process =
					ending_point- starting_point;


	device_ptr<vertex> temp_to =  device_malloc<vertex>(graph.number_of_edges * 2);
	device_ptr<vertex> temp_from =  device_malloc<vertex>(graph.number_of_edges * 2);

	/* Form temp to an temp from vector from edge arrray */
	thrust::copy(thrust::device, graph.full_edge_array + starting_point,
	graph.full_edge_array + ending_point,
	temp_to)	;

	thrust::copy(thrust::device,
		graph.full_edge_array + starting_point, graph.full_edge_array + ending_point,
		 temp_from);
		//
	thrust::transform(thrust::device,
		temp_from, temp_from + number_edges_to_process,
		temp_from, previous_el(current_level * graph.number_of_vertex + 1));

	/* Store begining and ending */
	thrust::copy(
		thrust::device,
		thrust::make_permutation_iterator(graph.full_vertex_array, temp_to),
		thrust::make_permutation_iterator(graph.full_vertex_array, temp_to+ number_edges_to_process),
		temp_to);

	thrust::copy(
		thrust::device,
		thrust::make_permutation_iterator(graph.full_vertex_array,
			temp_from),
		thrust::make_permutation_iterator(graph.full_vertex_array,
			temp_from + number_edges_to_process), temp_from);


	device_ptr<vertex>  process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
	/* Find all breaking points */
	thrust::transform(
		thrust::device,
		thrust::make_counting_iterator<vertex>(0),
		thrust::make_counting_iterator<vertex>(number_edges_to_process),
		process_vetxes,
		replacer(graph.full_vertex_array, graph.number_of_vertex));
	/* Sum all previous results */
		thrust::inclusive_scan(thrust::device, process_vetxes,
									process_vetxes + number_edges_to_process, process_vetxes);
	/* Number of edges to process*/
		device_ptr<vertex>  process_number = device_malloc<vertex>(number_edges_to_process+1);

		thrust::transform(
			thrust::device,
			make_zip_iterator(thrust::make_tuple(temp_from, temp_to)),
			make_zip_iterator(thrust::make_tuple(temp_from + number_edges_to_process, temp_to + number_edges_to_process)),
			process_number,
			counter()
			);

		thrust::inclusive_scan(thrust::device, process_number,
													 process_number + number_edges_to_process, process_number);


	// number_edges_to_process in a nutshell ?
	// int NUM = thrust::distance(temp_to.begin(), temp_to.end());

	/* 	Graph init. */


	device_ptr<vertex> tempo_vector =  device_malloc<vertex>(graph.number_of_vertex * graph.number_of_vertex);
	device_ptr<vertex> previous = tempo_vector;
	device_ptr<vertex> current = previous;


//	domain _temp_from = thrust::raw_pointer_cast(temp_from.data());
//	domain _temp_to = thrust::raw_pointer_cast(temp_to.data());
//	domain _full_vertex = graph.full_vertex_array;
//	domain _full_edge = graph.full_edge_array;
//	domain _vertex_data = thrust::raw_pointer_cast(process_vetxes.data());

	expander <<<1, 2 >>> (previous, current, process_vetxes , temp_from, temp_to,
		graph.full_vertex_array, graph.full_edge_array, graph.number_of_vertex);

			cout << "The okey 6" << endl;
	/*
	// Update vertex
	thrust::inclusive_scan(full_vertex_array.begin() + (number_of_vertex - 1), full_vertex_array.begin() + 2 * (number_of_vertex), full_vertex_array.begin() + (number_of_vertex - 1));
	// Update edges
	thrust::copy(tempo.begin(), current, full_edge_array.begin() + full_vertex_array[number_of_vertex-1]);


	print_csr_graph();
	*/
}
}


__global__ void sorter(thrust::device_ptr<vertex> full_edge_array, thrust::device_ptr<vertex> full_vertex_array)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int starting_point = 0;
	if (idx != 0)
	{
		starting_point = full_vertex_array[idx - 1] ;
	}
	int ending_point = full_vertex_array[idx];

	thrust::sort(thrust::device, full_edge_array+ starting_point, full_edge_array + ending_point);

}

void ordering_function(Graph graph)
{
	sorter<<<1, graph.number_of_vertex>>>(graph.full_edge_array, graph.full_vertex_array);
}


/*********************************
*	L opacity matrix calculation
*********************************/
/*
void calc_L_opacity(Graph graph)
{
	for (int i = 1; i < graph.L_VALUE; i++)
	{
		// full_edge_array - here we store all adjasent

		vertex N = graph.full_vertex_array[graph.number_of_vertex - 1];
		cout << "N+ " << N << endl;
		thrust::device_vector<vertex> from(N);
		// Forming indexes (from vertex)
		int starting_point = 0;
		int ending_point = graph.full_vertex_array[(i)*graph.number_of_vertex - 1];

		if (i != 1)

		{
			starting_point = graph.full_vertex_array[(i - 1)* graph.number_of_vertex - 1];
		}
		/*
		*	Expanding indexes. Finding break points
		*	Example: 0 1 2 3 4 .. 20 => 0 0 1 0 0 0 1 ...
		*/
		/*
		thrust::transform(
			thrust::device,
			thrust::make_counting_iterator<vertex>(starting_point),
			thrust::make_counting_iterator<vertex>(ending_point),
			from.begin(), replacer(thrust::raw_pointer_cast(graph.full_vertex_array.data()), graph.number_of_vertex)
			);

		// debug print
		cout << endl << "From degrees ";
		for (auto iter : from)
		{
			cout << "  " << iter;
		}


		//	from[0] = full_vertex_array[(number_of_vertex-1)*(i-1)];
		/*
		*	Transorming into indexes:
		*	Example:	0 0 1 0 0 0 1 => 0 0 1 1 1 1 2 2 2 ..
		*/
		/*
		thrust::inclusive_scan(thrust::device, from.begin(), from.end(), from.begin());

		cout << endl << "From degrees ";
		for (auto iter : from)
		{
			cout << "  " << iter;
		}


		/*
		*	Transforming from indexes into degrees:
		*	Example:  0 0 1 1 1 1 2 2 2.. => 2 2 4 4 4 4 ...
		*/
		/*
		thrust::transform(
			thrust::device,
			thrust::make_permutation_iterator(graph.vertex_degrees.begin(), from.begin()),
			thrust::make_permutation_iterator(graph.vertex_degrees.begin(), from.end()),
			from.begin(), thrust::identity<vertex>());

		cout << endl << "From degrees ";
		for (auto iter : from)
		{
			cout << "  " << iter;
		}


		/*
		*	To vector. Transform edge list into degree list =>  similar techno
		*
		*/

		/*
		thrust::device_vector<vertex> to(N);
		//	auto iter_begin = thrust::make_transform_iterator(full_edge_array.begin(), minus_one());
		//	auto iter_end =   thrust::make_transform_iterator(full_edge_array.begin() + N, minus_one());

		thrust::copy(thrust::device, graph.full_edge_array.begin(), graph.full_edge_array.begin() + N, to.begin());

		cout << endl << "TO degrees ";
		for (auto iter : to)
		{
			cout << "  " << iter;
		}


		thrust::transform(
			thrust::device,
			thrust::make_permutation_iterator(graph.vertex_degrees.begin(), to.begin()),
			thrust::make_permutation_iterator(graph.vertex_degrees.begin(), to.end()),
			to.begin(), thrust::identity<vertex>());

		cout << endl << "TO degrees ";
		for (auto iter : to)
		{
			cout << "  " << iter;
		}

		/*
		*  Find max and min in zip iterator of to - from pairs
		*/
		/*
		thrust::transform(
			thrust::device,
			thrust::make_zip_iterator(thrust::make_tuple(from.begin(), to.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(from.end(), to.end())),
			thrust::make_zip_iterator(thrust::make_tuple(from.begin(), to.begin())),
			min_max_transform()
			);
		cout << endl << "From degrees ";
		for (auto iter : from)
		{
			cout << "  " << iter;
		}

		cout << endl << "TO degrees ";
		for (auto iter : to)
		{
			cout << "  " << iter;
		}

		/*
		* 	Opacity  matrix forming. Now it is n^ 2 memory TODO: IN PARALLEL using cuda kernel
		* 	Assumptions !!: Not optimum for undericted (div 2).
		* 	Problem with same degree. Example: {4 = > 4} - must count only degree of one.
		*/
		/*
		for (int i = 0; i < N; i++)
		{
			double min = graph.degree_count[from[i] - 1] * graph.degree_count[to[i] - 1];
			if (graph.degree_count[from[i] - 1] == graph.degree_count[to[i] - 1])
				min = graph.degree_count[from[i] - 1];
			cout << endl << "FOR " << from[i] << " " << to[i] << " = " << min;
			graph.opacity_matrix[graph.max_degree*(from[i] - 1) + (to[i] - 1)] += 1.0 / (2.0*min);
		}


		/*
		* Sort by key. Indexes (values) and degrees (keys)
		*/

		/*
		* 	Reduce by key. Count how many pairs we have. 0 0 3 3
		*/
/*
	}
}
*/




int main()
{

	Graph graph;
	graph.init_test_graph(); // Reading graph from the file in COO format
	graph.print_coo_graph();
	graph.convert_to_CSR();
	graph.print_csr_graph();
  //form_full_level_graph(graph);
	ordering_function(graph);
	graph.print_csr_graph();
	/*
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

	*/

	return 0;
}
