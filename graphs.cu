/*
 	Main application
	Author : Bulat, 2015

*/
//#include "graph.cuh"
#include "apsp.cuh"
#include "graph.cuh"
#include "headers.h"



__global__  void expander(
	device_ptr<vertex> current_vertex, device_ptr<vertex> temp_from, device_ptr<vertex> temp_to,
	device_ptr<vertex> full_vertex_array, device_ptr<vertex> full_edge_array,	
	device_ptr<vertex> position_in_array, 
	device_ptr<vertex> expanded_array, device_ptr<vertex> from_vertex_array, 
	int number_edges_to_process
	)
{

	int idx = blockIdx.x*blockDim.x + threadIdx.x;

	int start_point_in_edge_list = 0;
	int offset_to_put_exp_array = 0;
	/*
	*	For searching in full_edge_list if already discovered
	*/
	if (idx != 0)
	{
		if (current_vertex[idx]!= 0)
			start_point_in_edge_list = full_vertex_array[current_vertex[idx]-1];
		offset_to_put_exp_array = position_in_array[idx - 1]; // reserved position in expanded array
	}
	int end_point_in_edge_list = full_vertex_array[current_vertex[idx]];
	// DEbug print TODO: remove
	printf("Expander ok 0 \n");
	/*
	*	Copy to expanded array if the edge is unique (was not discovered previouslly, not equal to vertex itself)
	*	Result:			1 2 1 .... (expanded edges)
	*/
	thrust::device_ptr<vertex> current_position =
 	thrust::copy_if(thrust::device,
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_from[idx])),
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_to[idx])),
		expanded_array + offset_to_put_exp_array,
		unique_edge(full_edge_array + start_point_in_edge_list,full_edge_array + end_point_in_edge_list,
					current_vertex[idx]));

	int planed_size = temp_to[idx] - temp_from[idx];
	int real_size = thrust::distance(expanded_array + offset_to_put_exp_array, current_position);
	int starting = current_vertex[idx];
	/*
	*	Expand the current processing vertex to the size *real size*;
	*			Result : 0 0 0 1 1 ... (the vertex from expanded)
	*/
	thrust::copy(thrust::device,
		thrust::make_constant_iterator(starting),
		thrust::make_constant_iterator(starting) + real_size,
		from_vertex_array + offset_to_put_exp_array);

	if (planed_size != real_size)
		if (idx != 0)
	{
		thrust::transform(thrust::device, position_in_array + idx - 1,
				position_in_array + number_edges_to_process, position_in_array + idx -1, minus_value(planed_size - real_size));
	}
		else
		{
			thrust::transform(thrust::device, position_in_array,
				position_in_array + number_edges_to_process, position_in_array, minus_value(planed_size - real_size));
		}
		
	printf("Expander ok 6 \n");
}


/*
*	
*
*/
__global__ void sorter(thrust::device_ptr<vertex> full_edge_array,
						thrust::device_ptr<vertex> full_vertex_array)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int starting_point = 0;
	if (idx != 0)
	{
		starting_point = full_vertex_array[idx - 1];
	}
	int ending_point = full_vertex_array[idx];
	// sort 
	thrust::sort(thrust::device, full_edge_array + starting_point, full_edge_array + ending_point);
	
}
/*
*
*
*/


__global__ void unifier(
	device_ptr<vertex> expanded_array,
	device_ptr<vertex> positions_vertex_current_level,
	device_ptr<vertex> current_ending_offset
	)
	{
		int idx = blockIdx.x*blockDim.x + threadIdx.x;
		int start_point = 0;
	
		if (idx != 0)
		{
			start_point = positions_vertex_current_level[idx - 1];
		
		}
		int end_point = positions_vertex_current_level[idx];

		printf("So far so good \n");
		thrust::sort(thrust::device, expanded_array + start_point, expanded_array + end_point);
		printf("Boom fuck you sort \n");
		// remove dublicates
		thrust::device_ptr<vertex> current_position =
			thrust::unique(thrust::device, expanded_array + start_point, expanded_array + end_point);
		printf("Unique fuck you sort \n");
		vertex real_size = thrust::distance(expanded_array + start_point, current_position);
		current_ending_offset[idx] = real_size;
		
	}


__global__ void edge_copier(
	device_ptr<vertex> expanded_array,
	device_ptr<vertex> positions_vertex_current_level,
	device_ptr<vertex> current_ending_offset,
	device_ptr<vertex> full_vertex_array,
	device_ptr<vertex> full_edge_array,
	int L_VALUE,
	int number_of_vertex
	)

{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int start_point = 0;

	if (idx != 0)
	{
		start_point = positions_vertex_current_level[idx - 1];

	}
	int end_point = start_point + current_ending_offset[idx];

	int edge_put_list_start = full_vertex_array[L_VALUE *number_of_vertex + idx - 1];
	
	thrust::copy(thrust::device, expanded_array + start_point, expanded_array + end_point, full_edge_array + edge_put_list_start);
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
	int added_offset = 0;
	if (current_level != 1)
		added_offset = graph.full_vertex_array[(current_level - 1) * graph.number_of_vertex - 1];

	device_ptr<vertex> temp_to =  device_malloc<vertex>(graph.number_of_edges * 2);
	device_ptr<vertex> temp_from =  device_malloc<vertex>(graph.number_of_edges * 2);

	/* Form temp to an temp from vector from edge arrray */
	thrust::copy(thrust::device, graph.full_edge_array + starting_point,
	graph.full_edge_array + ending_point,
	temp_to);

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

	/*
	*	Array of vertex, from which we will expand. Proces vertex
	*/
	device_ptr<vertex>  process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
	/* Find all breaking points. 0 0 1 0 0 0 1 ... */
	thrust::transform(
		thrust::device,
		thrust::make_counting_iterator<vertex>(0),
		thrust::make_counting_iterator<vertex>(number_edges_to_process),
		process_vetxes,
		replacer(graph.full_vertex_array, graph.number_of_vertex));
	/*
		Sum all previous results (sum of breaking points).
		Result: 0 0 1 1 1 1 2 2 2 (vertex array)
	*/
		thrust::inclusive_scan(thrust::device, process_vetxes,
									process_vetxes + number_edges_to_process, process_vetxes);
	/*
		Offset array, step 1:
		Result <- Temp_TO - Temp_FROM
	*/
		device_ptr<vertex>  position_in_array = device_malloc<vertex>(number_edges_to_process+1);

		thrust::transform(
			thrust::device,
			make_zip_iterator(thrust::make_tuple(temp_from, temp_to)),
			make_zip_iterator(thrust::make_tuple(temp_from + number_edges_to_process,
												temp_to + number_edges_to_process)),
			position_in_array,
			counter());
	/*
		Forming offset array from process number, step 2:
		2 4 4 => 2 6 10
	*/
	thrust::inclusive_scan(thrust::device, position_in_array,
					 		    position_in_array + number_edges_to_process,
							    position_in_array);



	
	device_ptr<vertex> expanded_array = device_malloc<vertex>(position_in_array[number_edges_to_process - 1]);
	thrust::fill(thrust::device, expanded_array, expanded_array + position_in_array[number_edges_to_process - 1], -1);
	// process number contains the maximum needed memory to store if all vertexes are unique
	device_ptr<vertex> from_vertex_array = device_malloc<vertex>(position_in_array[number_edges_to_process - 1]);
	thrust::fill(thrust::device, from_vertex_array, from_vertex_array + position_in_array[number_edges_to_process - 1], -1);

	int prev_max_position = position_in_array[number_edges_to_process - 1];
	//  Expand array on one level
	//	Can contain non unique values
	expander<<< 1, number_edges_to_process >>>(
		process_vetxes, temp_from, temp_to,
		graph.full_vertex_array, graph.full_edge_array,
		position_in_array, 
		expanded_array, from_vertex_array,
		number_edges_to_process
		);
	//cudaThreadSynchronize();
	cudaDeviceSynchronize();
	cout << endl << "Expander finished work" << endl;
	/*
	*	Remove empty, non used data
	*/
	thrust::remove(thrust::device, expanded_array, expanded_array + prev_max_position, -1);
	thrust::remove(thrust::device, from_vertex_array, from_vertex_array + prev_max_position, -1);
	cout << endl << "Removing ok" << endl;
	/*
	*	Form vertex offset list
	*/
	
	thrust::device_ptr<vertex> positions_vertex_current_level = thrust::device_malloc<vertex>(graph.number_of_vertex);
	thrust::reduce_by_key(thrust::device,
		from_vertex_array, from_vertex_array + position_in_array[number_edges_to_process-1],
		thrust::make_constant_iterator(1),
		thrust::make_discard_iterator(),
		positions_vertex_current_level);
	// STEP 2: Forming offset
	thrust::inclusive_scan(thrust::device, positions_vertex_current_level,
							positions_vertex_current_level + graph.number_of_vertex, positions_vertex_current_level);


	thrust::device_ptr<vertex> vertex_ending_offsets = thrust::device_malloc<vertex>(graph.number_of_vertex);
	
	cout << "OMG it works untill now !!" << endl;


	
	unifier <<<1, graph.number_of_vertex >>>( expanded_array, positions_vertex_current_level, vertex_ending_offsets);
	cudaDeviceSynchronize();
	cout << " I think it works" << endl;
	
	thrust::device_ptr<vertex> position_in_edge_list = thrust::device_malloc<vertex>(graph.number_of_vertex);
	cout << "Yep here we go" << endl;
	
	thrust::copy(thrust::device, vertex_ending_offsets, vertex_ending_offsets + graph.number_of_vertex, 
		graph.full_vertex_array + current_level * graph.number_of_vertex);
		
	
	cout << "BAM THAT IS IT" << endl;
	thrust::inclusive_scan(thrust::device, graph.full_vertex_array + current_level * graph.number_of_vertex - 1,
		graph.full_vertex_array + (current_level + 1) * graph.number_of_vertex, graph.full_vertex_array + current_level * graph.number_of_vertex - 1);
	

	int* a = new int[100];
	thrust::copy(graph.full_vertex_array, graph.full_vertex_array + 14, a);

	for (int i = 0; i < 14; i++)
	{
		cout << a[i] << " ";
	}
	cout << endl;


	
	cout << "OMG o_______O !!" << endl;

	edge_copier<<<1, graph.number_of_vertex>>>(
		expanded_array,
		positions_vertex_current_level,
		vertex_ending_offsets,
		graph.full_vertex_array,
		graph.full_edge_array,
		current_level,
		graph.number_of_vertex);
	cudaDeviceSynchronize(); 


	cout << "The okey 6" << endl;
	}
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
	ordering_function(graph);
	form_full_level_graph(graph);
	
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
