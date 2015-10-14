/*
 	Main application
	Author : Bulat, 2015

*/
//#include "graph.cuh"
#include "apsp.cuh"
#include "graph.cuh"
#include "headers.h"


/*
*	Expand array according to it's from and to values
*	Input : device_ptr<vertex> expanded_array
*			device_ptr<vertex> position_current_level
*			device_ptr<vertex> current_ending_offset
*	Out: 	sorted and normalized expanded_array
*/
__global__  void expander(
	device_ptr<vertex> current_vertex, device_ptr<vertex> temp_from, device_ptr<vertex> temp_to,
	device_ptr<vertex> full_vertex_array, device_ptr<vertex> full_edge_array,
	device_ptr<vertex> position_in_array,
	device_ptr<vertex> expanded_array, device_ptr<vertex> from_vertex_array,
	int number_edges_to_process, int number_of_vertex,
	int current_level,
	device_ptr<vertex> vertex_offset
	)
{

	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	int offset_to_put_exp_array = 0;
	/*
	*	For searching in full_edge_list if already discovered
	*/
	if (idx != 0)
	{
		offset_to_put_exp_array = position_in_array[idx - 1]; // reserved position in expanded array
	}
	/*
	*	Copy to expanded array if the edge is unique (was not discovered previously, not equal to vertex itself)
	*	Result:			1 2 1 .... (expanded edges)
	*/
	thrust::device_ptr<vertex> current_position =
 	thrust::copy_if(thrust::device,
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_from[idx])),
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_to[idx])),
		expanded_array + offset_to_put_exp_array,
		unique_edge(full_vertex_array, full_edge_array, current_vertex[idx],
		number_of_vertex, current_level));

	int planed_size = temp_to[idx] - temp_from[idx];
	int real_size = thrust::distance(expanded_array + offset_to_put_exp_array, current_position);
	int starting = current_vertex[idx];
	// TODO : check real size value

	/*
	*	Expand the current processing vertex to the size *real size*;
	*			Result : 0 0 0 1 1 ... (the vertex from expanded)
	*/
	thrust::copy(thrust::device,
		thrust::make_constant_iterator(starting),
		thrust::make_constant_iterator(starting) + real_size,
		from_vertex_array + offset_to_put_exp_array);
	vertex* k = thrust::raw_pointer_cast(vertex_offset + starting);
	atomicAdd(k, real_size);

	if (planed_size != real_size)
		if (idx != 0)
	{
			thrust::transform(thrust::device, position_in_array + idx - 1,
			position_in_array + number_edges_to_process,
			position_in_array + idx -1, minus_value(planed_size - real_size));
	}
		else
		{
			thrust::transform(thrust::device, position_in_array,
				position_in_array + number_edges_to_process, position_in_array, minus_value(planed_size - real_size));
		}
}


/*
*	Sorting function for each vertex offset position. Initial sort (not for L > 1)
*	Input: device_ptr full_edge_array
*		   device_ptr full_vertex_array
*	Out : Sorted full_edge_array
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
*	Sort and remove duplicates in edge_array
*	Input : device_ptr<vertex> expanded_array
*			device_ptr<vertex> position_current_level
*			device_ptr<vertex> current_ending_offset
*	Out: 	sorted and normalized expanded_array
*/
__global__ void unifier(
	device_ptr<vertex> expanded_array,
	device_ptr<vertex> positions_vertex_current_level,
	device_ptr<vertex> current_ending_offset)
	{
		int idx = blockIdx.x*blockDim.x + threadIdx.x;
		int start_point = 0;

		if (idx != 0)
		{
			start_point = positions_vertex_current_level[idx - 1];

		}
		int end_point = positions_vertex_current_level[idx];
		if (end_point > start_point)
		{
			thrust::sort(thrust::device, expanded_array + start_point, expanded_array + end_point);
			// remove dublicates
			thrust::device_ptr<vertex> current_position =
				thrust::unique(thrust::device, expanded_array + start_point, expanded_array + end_point);
			vertex real_size = thrust::distance(expanded_array + start_point, current_position);
			current_ending_offset[idx] = real_size;
		}
		else
		{
			current_ending_offset[idx] = 0;
		}


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

__global__ void opacity_former(
	device_ptr<vertex> from,
	device_ptr<vertex> to,
	device_ptr<vertex> degree_count,
	device_ptr<opacity> opacity_matrix,
	int max_degree
	)
	{
		// TODO: Not sure about this
		int i = blockIdx.x*blockDim.x + threadIdx.x;
		opacity min = degree_count[from[i] - 1] * degree_count[to[i] - 1];
		if (degree_count[from[i] - 1] == degree_count[to[i] - 1])
			min = degree_count[from[i] - 1];
		opacity* k = thrust::raw_pointer_cast(opacity_matrix + max_degree*(from[i] - 1) + (to[i] - 1));
		opacity added_value = 1.0/ (2.0 * min);
		atomicAdd(k, added_value);

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

	vertex number_edges_to_process = ending_point- starting_point;
	// Value to add as an offset (previous end - current begin)
	int added_offset = starting_point;

	device_ptr<vertex> temp_to =  device_malloc<vertex>(number_edges_to_process);
	device_ptr<vertex> temp_from =  device_malloc<vertex>(number_edges_to_process);

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
		thrust::make_permutation_iterator(graph.full_vertex_array, temp_to + number_edges_to_process),
		temp_to);

	thrust::copy(
		thrust::device,
		thrust::make_permutation_iterator(graph.full_vertex_array, temp_from),
		thrust::make_permutation_iterator(graph.full_vertex_array,
		temp_from + number_edges_to_process), temp_from);


	/*
	*	Array of vertex, from which we will expand. Proces vertex
	*/
	device_ptr<vertex>  process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
	/* Find all breaking points. 0 0 1 0 0 0 1 ... */
	thrust::transform(
		thrust::device,
		thrust::make_counting_iterator<vertex>(added_offset),
		thrust::make_counting_iterator<vertex>(number_edges_to_process + added_offset),
		process_vetxes,
		replacer(graph.full_vertex_array + (current_level-1)* graph.number_of_vertex, graph.number_of_vertex));
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
		device_ptr<vertex>  position_in_array = device_malloc<vertex>(number_edges_to_process);

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

	int grid_size = number_edges_to_process;
	thrust::device_ptr<vertex> positions_vertex_current_level = thrust::device_malloc<vertex>(graph.number_of_vertex);
	thrust::fill(thrust::device, positions_vertex_current_level, positions_vertex_current_level + graph.number_of_vertex, 0);

	expander<<< 1, grid_size >>>(
		process_vetxes, temp_from, temp_to,
		graph.full_vertex_array, graph.full_edge_array,
		position_in_array,
		expanded_array, from_vertex_array,
		number_edges_to_process,
		graph.number_of_vertex,
		current_level,
		positions_vertex_current_level);
	//cudaThreadSynchronize();
	cudaDeviceSynchronize();

	device_free(temp_from);
	device_free(temp_to);
	device_free(process_vetxes);

	/*
	*	Remove empty, non used data
		*/
	thrust::remove(thrust::device, expanded_array, expanded_array + prev_max_position, -1);
	thrust::remove(thrust::device, from_vertex_array, from_vertex_array + prev_max_position, -1);

	/*
	*	Form vertex offset list
	*/
	//	int gridsize = position_in_array[number_edges_to_process - 1];

	// STEP 2: Forming offset
	thrust::inclusive_scan(thrust::device, positions_vertex_current_level,
							positions_vertex_current_level + graph.number_of_vertex, positions_vertex_current_level);

	thrust::device_ptr<vertex> vertex_ending_offsets = thrust::device_malloc<vertex>(graph.number_of_vertex);
	unifier <<<1, graph.number_of_vertex >>>( expanded_array, positions_vertex_current_level, vertex_ending_offsets);


	thrust::device_ptr<vertex> position_in_edge_list = thrust::device_malloc<vertex>(graph.number_of_vertex);

	cout << endl;
	cudaDeviceSynchronize();

	thrust::copy(thrust::device, vertex_ending_offsets, vertex_ending_offsets + graph.number_of_vertex,
		graph.full_vertex_array + current_level * graph.number_of_vertex);



	thrust::inclusive_scan(thrust::device, graph.full_vertex_array + current_level * graph.number_of_vertex - 1,
		graph.full_vertex_array + (current_level + 1) * graph.number_of_vertex, graph.full_vertex_array + current_level * graph.number_of_vertex - 1);


	grid_size = graph.number_of_vertex;

	edge_copier<<<1, grid_size>>>(
		expanded_array,
		positions_vertex_current_level,
		vertex_ending_offsets,
		graph.full_vertex_array,
		graph.full_edge_array,
		current_level,
		graph.number_of_vertex);
	cudaDeviceSynchronize();
	device_free(expanded_array);
	device_free(positions_vertex_current_level);
	device_free(vertex_ending_offsets);


	}
}



void ordering_function(Graph graph)
{
	sorter<<<1, graph.number_of_vertex>>>(graph.full_edge_array, graph.full_vertex_array);
}


/*********************************
*	L opacity matrix calculation
*********************************/

void calc_L_opacity(Graph graph)
{
	for (int i = 1; i <= graph.L_VALUE; i++)
	{
		// full_edge_array - here we store all adjasent

		// Forming indexes (from vertex)
		int starting_point = 0;
		int ending_point = graph.full_vertex_array[(i)*graph.number_of_vertex - 1];

		if (i != 1)

		{
			starting_point = graph.full_vertex_array[(i - 1)* graph.number_of_vertex - 1];
		}

		vertex N = ending_point - starting_point;
		device_ptr<vertex> from = device_malloc<vertex>(N);

		/*
		*	Expanding indexes. Finding break points
		*	Example: 0 1 2 3 4 .. 20 => 0 0 1 0 0 0 1 ...
		*/

		thrust::transform(
			thrust::device,
			thrust::make_counting_iterator<vertex>(starting_point),
			thrust::make_counting_iterator<vertex>(ending_point),
			from, replacer(graph.full_vertex_array + (i-1)* graph.number_of_vertex, graph.number_of_vertex)
			);



		//	from[0] = full_vertex_array[(number_of_vertex-1)*(i-1)];
		/*
		*	Transorming into indexes:
		*	Example:	0 0 1 0 0 0 1 => 0 0 1 1 1 1 2 2 2 ..
		*/

		thrust::inclusive_scan(thrust::device, from, from + N , from);


		/*
		*	Transforming from indexes into degrees:
		*	Example:  0 0 1 1 1 1 2 2 2.. => 2 2 4 4 4 4 ...
		*/

		thrust::transform(
			thrust::device,
			thrust::make_permutation_iterator(graph.vertex_degrees, from),
			thrust::make_permutation_iterator(graph.vertex_degrees, from + N),
			from, thrust::identity<vertex>());

		/*
		*	To vector. Transform edge list into degree list =>  similar techno
		*
		*/


		thrust::device_ptr<vertex> to = device_malloc<vertex>(N);
		//	auto iter_begin = thrust::make_transform_iterator(full_edge_array.begin(), minus_one());
		//	auto iter_end =   thrust::make_transform_iterator(full_edge_array.begin() + N, minus_one());

		thrust::copy(thrust::device, graph.full_edge_array + starting_point, graph.full_edge_array + ending_point, to);

		thrust::transform(
			thrust::device,
			thrust::make_permutation_iterator(graph.vertex_degrees, to),
			thrust::make_permutation_iterator(graph.vertex_degrees, to + N),
			to, thrust::identity<vertex>());

		/*
		*  Find max and min in zip iterator of to - from pairs
		*/


		thrust::transform(
			thrust::device,
			thrust::make_zip_iterator(thrust::make_tuple(from, to)),
			thrust::make_zip_iterator(thrust::make_tuple(from + N, to + N)),
			thrust::make_zip_iterator(thrust::make_tuple(from, to)),
			min_max_transform());

		/*
		* 	Opacity  matrix forming. Now it is n^ 2 memory TODO: IN PARALLEL using cuda kernel
		* 	Assumptions !!: Not optimum for undericted (div 2).
		* 	Problem with same degree. Example: {4 = > 4} - must count only degree of one.
		*/
		int gridsize = N;
		opacity_former<<<1, gridsize>>>(from, to, graph.degree_count, graph.opacity_matrix, graph.max_degree);

		/*
		* Sort by key. Indexes (values) and degrees (keys)
		*/

		/*
		* 	Reduce by key. Count how many pairs we have. 0 0 3 3
		*/

	}
}


int main(int argc, char* argv[])
{

	Graph graph;
	cout << "Converting to " << endl;
	int l_value = std::atoi(argv[1]);
	cout << l_value << " L value" << endl;
	graph.L_VALUE = l_value;
	graph.init_test_graph(); // Reading graph from the file in COO format
	graph.print_coo_graph();
	graph.convert_to_CSR();
	graph.print_csr_graph();
	ordering_function(graph);



	UINT wTimerRes = 0;
	bool init = InitMMTimer(wTimerRes);
	DWORD startTime = timeGetTime();

		form_full_level_graph(graph);

	// BFS
	//graph.single_bfs(2);
//	graph.form_full_level_graph();
//	graph.print_csr_graph();
		calc_L_opacity(graph);
  	graph.print_opacity_matrix();

	unsigned int endTime = timeGetTime();
	unsigned int gpu_time = unsigned int(endTime - startTime);
	printf("GPU Timing(including all device-host, host-device copies, device allocations and freeing of device memory): %dms\n\n", gpu_time);
	DestroyMMTimer(wTimerRes, init);

	graph.print_csr_graph();


	return 0;
}
