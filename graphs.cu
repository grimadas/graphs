/*
 	Main application
	Author : Bulat, 2015

*/
//#include "graph.cuh"
#include "apsp.cuh"
#include "graph.cuh"
#include "headers.h"



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
	device_ptr<vertex> temp_vertexes =  device_malloc<vertex>(graph.number_of_vertex);
	device_ptr<vertex> process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
	thrust::copy(thrust::device, graph.full_vertex_array + (current_level-1)*graph.number_of_vertex, graph.full_vertex_array + (current_level)*graph.number_of_vertex , temp_vertexes);
	thrust::adjacent_difference(thrust::device, temp_vertexes, temp_vertexes + graph.number_of_vertex, temp_vertexes);
	temp_vertexes[0] = temp_vertexes[0] - starting_point;

	expand(temp_vertexes, temp_vertexes + graph.number_of_vertex, thrust::make_counting_iterator<vertex>(0), process_vetxes);
	device_free(temp_vertexes);

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
		device_ptr<vertex> from_vertex = device_malloc<vertex>(graph.number_of_vertex);
		device_ptr<vertex> from = device_malloc<vertex>(N);

		/*
		*	Expanding indexes. Finding break points
		*	Example: 0 1 2 3 4 .. 20 => 0 0 1 0 0 0 1 ...
		*/
		thrust::copy(thrust::device, graph.full_vertex_array + (i-1)*graph.number_of_vertex, graph.full_vertex_array + (i)*graph.number_of_vertex , from_vertex);
		thrust::adjacent_difference(thrust::device, from_vertex, from_vertex + graph.number_of_vertex, from_vertex);

		from_vertex[0] = from_vertex[0] - starting_point;

   	expand(from_vertex, from_vertex + graph.number_of_vertex, thrust::make_counting_iterator<vertex>(0), from);
/*
		thrust::transform(
			thrust::device,
			thrust::make_counting_iterator<vertex>(starting_point),
			thrust::make_counting_iterator<vertex>(ending_point),
			from, replacer(graph.full_vertex_array + (i-1)* graph.number_of_vertex, graph.number_of_vertex)
			);
			*/



		//	from[0] = full_vertex_array[(number_of_vertex-1)*(i-1)];
		/*
		*	Transorming into indexes:
		*	Example:	0 0 1 0 0 0 1 => 0 0 1 1 1 1 2 2 2 ..
		*/

	//	thrust::inclusive_scan(thrust::device, from, from + N , from);


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
//		thrust::fill(thrust::device, graph.opacity_matrix, graph.opacity_matrix + graph.max_degree * graph.max_degree, 0);
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
	graph.convert_to_CSR();
	ordering_function(graph);

	UINT wTimerRes = 0;
	bool init = InitMMTimer(wTimerRes);
	DWORD startTime = timeGetTime();

	form_full_level_graph(graph);
	calc_L_opacity(graph);

	unsigned int endTime = timeGetTime();
	unsigned int gpu_time = unsigned int(endTime - startTime);
	printf("GPU Timing(including all device-host, host-device copies, device allocations and freeing of device memory): %dms\n\n", gpu_time);
	DestroyMMTimer(wTimerRes, init);

	graph.print_opacity_matrix();

	graph.print_csr_graph();


	return 0;
}
