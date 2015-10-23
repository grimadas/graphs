/*
 	Main application
	Author : Bulat, 2015

*/
//#include "apsp.cuh"
#include "graph.cuh"




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
			copy(device, graph.full_edge_array + starting_point,
			graph.full_edge_array + ending_point,
			temp_to);

			copy(device,
				graph.full_edge_array + starting_point, graph.full_edge_array + ending_point,
				 temp_from);
			std::cout << "Ok in 1" << std::endl;
				//
				transform(device,
					temp_from, temp_from + number_edges_to_process,
					temp_from, previous_el(current_level * graph.number_of_vertex + 1));

			/* Store begining and ending */
			copy(
				device,
				make_permutation_iterator(graph.full_vertex_array, temp_to),
				make_permutation_iterator(graph.full_vertex_array, temp_to + number_edges_to_process),
				temp_to);

			copy(
				device,
				make_permutation_iterator(graph.full_vertex_array, temp_from),
				make_permutation_iterator(graph.full_vertex_array,
				temp_from + number_edges_to_process), temp_from);
				std::cout << "Ok in 2" << std::endl;

			/*
			*	Array of vertex, from which we will expand. Proces vertex
			*/
			device_ptr<vertex> temp_vertexes =  device_malloc<vertex>(graph.number_of_vertex);
			device_ptr<vertex> process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
			copy(device, graph.full_vertex_array + (current_level-1)*graph.number_of_vertex, graph.full_vertex_array + (current_level)*graph.number_of_vertex , temp_vertexes);
			adjacent_difference(device, temp_vertexes, temp_vertexes + graph.number_of_vertex, temp_vertexes);
			temp_vertexes[0] = temp_vertexes[0] - starting_point;

			expand(temp_vertexes, temp_vertexes + graph.number_of_vertex, make_counting_iterator<vertex>(0), process_vetxes);
			device_free(temp_vertexes);
			std::cout << "Ok in 3" << std::endl;
			/*
				Offset array, step 1:
				Result <- Temp_TO - Temp_FROM
			*/
			device_ptr<vertex>  position_in_array = device_malloc<vertex>(number_edges_to_process);

			transform(
					device,
					make_zip_iterator(make_tuple(temp_from, temp_to)),
					make_zip_iterator(make_tuple(temp_from + number_edges_to_process,
														temp_to + number_edges_to_process)),
					position_in_array,
					counter());
				std::cout << "Ok in 4" << std::endl;
				/*
					Forming offset array from process number, step 2:
					2 4 4 => 2 6 10
					*/
				inclusive_scan(device, position_in_array,
							 		    position_in_array + number_edges_to_process,
									    position_in_array);

			//TODO: check position in array
			domain prints = new vertex[number_edges_to_process];
			copy(position_in_array, position_in_array + number_edges_to_process, prints);
			std::cout << "Position In array:  " << std::endl;
			for(int i=0; i < number_edges_to_process; i++)
			{
				std::cout << prints[i]<<" ";
			}

			std::cout << "Ok in 5" << std::endl;

			device_ptr<vertex> expanded_array = device_malloc<vertex>(position_in_array[number_edges_to_process - 1]);
			fill(device, expanded_array, expanded_array + position_in_array[number_edges_to_process - 1], -1);
			// process number contains the maximum needed memory to store if all vertexes are unique
			device_ptr<vertex> from_vertex_array = device_malloc<vertex>(position_in_array[number_edges_to_process - 1]);
			fill(device, from_vertex_array, from_vertex_array + position_in_array[number_edges_to_process - 1], -1);

			int prev_max_position = position_in_array[number_edges_to_process - 1];
			//  Expand array on one level
			//	Can contain non unique values
			std::cout << "Ok in 6" << std::endl;

			device_ptr<vertex> positions_vertex_current_level = device_malloc<vertex>(graph.number_of_vertex);
			fill(device, positions_vertex_current_level, positions_vertex_current_level + graph.number_of_vertex, 0);

			int grid_size = (number_edges_to_process + BLOCK_SIZE - 1) / BLOCK_SIZE;
			// Expand here and put to expanded_array
			expander<<< grid_size, BLOCK_SIZE >>>(
				process_vetxes, temp_from, temp_to,
				graph.full_vertex_array, graph.full_edge_array,
				position_in_array,
				expanded_array, from_vertex_array,
				number_edges_to_process,
				graph.number_of_vertex,
				current_level,
				positions_vertex_current_level);

			std::cout << "Ok in 7" << std::endl;
			//cudaThreadSynchronize();
			cudaDeviceSynchronize();
			std::cout << "Ok in 9" << std::endl;
			device_free(temp_from);
			device_free(temp_to);
			device_free(process_vetxes);

			/*
			*	Remove empty, non used data
				*/
			remove(device, expanded_array, expanded_array + prev_max_position, -1);
			remove(device, from_vertex_array, from_vertex_array + prev_max_position, -1);
			std::cout << "Ok in 10" << std::endl;

			prints = new vertex[position_in_array[number_edges_to_process-1]];
			copy(expanded_array, expanded_array + position_in_array[number_edges_to_process-1], prints);
			std::cout << "After expanded : " << std::endl;
			for(int i=0; i <  position_in_array[number_edges_to_process-1]; i++)
			{
				std::cout << prints[i]<<" ";
			}
			/*
			*	Form vertex offset list
			*/
			//	int gridsize = position_in_array[number_edges_to_process - 1];

			// STEP 2: Forming offset
			inclusive_scan(device, positions_vertex_current_level,
									positions_vertex_current_level + graph.number_of_vertex, positions_vertex_current_level);


			copy(positions_vertex_current_level, positions_vertex_current_level + graph.number_of_vertex, prints);
			std::cout << "Position after expanded : " << std::endl;
			for(int i=0; i < graph.number_of_vertex; i++)
			{
					std::cout << prints[i]<<" ";
			}


			device_ptr<vertex> vertex_ending_offsets = device_malloc<vertex>(graph.number_of_vertex);
			grid_size =  (graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE;
			unifier <<<grid_size, BLOCK_SIZE >>>( expanded_array, positions_vertex_current_level, vertex_ending_offsets, graph.number_of_vertex);
			cudaDeviceSynchronize();

			copy(vertex_ending_offsets, vertex_ending_offsets + graph.number_of_vertex, prints);
			std::cout << "vertex_ending_offsets : " << std::endl;
			for(int i=0; i < graph.number_of_vertex; i++)
			{
					std::cout << prints[i]<<" ";
			}

			copy(expanded_array, expanded_array + position_in_array[number_edges_to_process - 1], prints);
			std::cout <<std::endl << "Exapnded array now " << std::endl;
			for(int i=0; i < position_in_array[number_edges_to_process - 1]; i++)
			{
					std::cout << prints[i]<<" ";
			}

			std::cout << "Ok in 11" << std::endl;

			copy(device, vertex_ending_offsets, vertex_ending_offsets + graph.number_of_vertex,
				graph.full_vertex_array + current_level * graph.number_of_vertex);

			inclusive_scan(device, graph.full_vertex_array + current_level * graph.number_of_vertex - 1,
				graph.full_vertex_array + (current_level + 1) * graph.number_of_vertex, graph.full_vertex_array + current_level * graph.number_of_vertex - 1);

			std::cout << "Ok in 12" << std::endl;
			grid_size =  (graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE;
			edge_copier<<<grid_size, BLOCK_SIZE>>>(
				expanded_array,
				positions_vertex_current_level,
				vertex_ending_offsets,
				graph.full_vertex_array,
				graph.full_edge_array,
				current_level,
				graph.number_of_vertex);

			cudaDeviceSynchronize();
			std::cout << "Ok in 13" << std::endl;
			device_free(expanded_array);
			device_free(positions_vertex_current_level);
			device_free(vertex_ending_offsets);
			device_free(position_in_array);
			delete prints;
			std::cout << "Ok in 14" << std::endl;


	}
}



void ordering_function(Graph graph)
{
	int	gridsize = (graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE;
	sorter<<<gridsize, BLOCK_SIZE>>>(graph.full_edge_array, graph.full_vertex_array, graph.number_of_vertex);
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
		copy(device, graph.full_vertex_array + (i-1)*graph.number_of_vertex, graph.full_vertex_array + (i)*graph.number_of_vertex , from_vertex);
		adjacent_difference(device, from_vertex, from_vertex + graph.number_of_vertex, from_vertex);

		from_vertex[0] = from_vertex[0] - starting_point;
		std::cout << "Before Expanding opacity";
   	expand(from_vertex, from_vertex + graph.number_of_vertex, make_counting_iterator<vertex>(0), from);
/*
		transform(
			device,
			make_counting_iterator<vertex>(starting_point),
			make_counting_iterator<vertex>(ending_point),
			from, replacer(graph.full_vertex_array + (i-1)* graph.number_of_vertex, graph.number_of_vertex)
			);
			*/



		//	from[0] = full_vertex_array[(number_of_vertex-1)*(i-1)];
		/*
		*	Transorming into indexes:
		*	Example:	0 0 1 0 0 0 1 => 0 0 1 1 1 1 2 2 2 ..
		*/

	//	inclusive_scan(device, from, from + N , from);


		/*
		*	Transforming from indexes into degrees:
		*	Example:  0 0 1 1 1 1 2 2 2.. => 2 2 4 4 4 4 ...
		*/

		transform(
			device,
			make_permutation_iterator(graph.vertex_degrees, from),
			make_permutation_iterator(graph.vertex_degrees, from + N),
			from, identity<vertex>());

		/*
		*	To vector. Transform edge list into degree list =>  similar techno
		*
		*/


		device_ptr<vertex> to = device_malloc<vertex>(N + 1);
		//	auto iter_begin = make_transform_iterator(full_edge_array.begin(), minus_one());
		//	auto iter_end =   make_transform_iterator(full_edge_array.begin() + N, minus_one());
		std::cout << "Before Expanding opacity 2 statt = " << starting_point << " Ending "  << ending_point;
		copy(device, graph.full_edge_array + starting_point, graph.full_edge_array + ending_point, to);
		std::cout << "Before Expanding opacity 3 ";
		transform(
			device,
			make_permutation_iterator(graph.vertex_degrees, to),
			make_permutation_iterator(graph.vertex_degrees, to + N),
			to, identity<vertex>());

			std::cout << "Before Expanding opacity 4 ";

		/*
		*  Find max and min in zip iterator of to - from pairs
		*/


		transform(
			device,
			make_zip_iterator(make_tuple(from, to)),
			make_zip_iterator(make_tuple(from + N, to + N)),
			make_zip_iterator(make_tuple(from, to)),
			min_max_transform());
			std::cout << "Before Expanding opacity 5";

		/*
		* 	Opacity  matrix forming. Now it is n^ 2 memory TODO: IN PARALLEL using cuda kernel
		* 	Assumptions !!: Not optimum for undericted (div 2).
		* 	Problem with same degree. Example: {4 = > 4} - must count only degree of one.
		*/
//		fill(device, graph.opacity_matrix, graph.opacity_matrix + graph.max_degree * graph.max_degree, 0);

		int gridsize =(N + BLOCK_SIZE - 1) / BLOCK_SIZE;
		opacity_former<<<gridsize, BLOCK_SIZE>>>(from, to, graph.degree_count, graph.opacity_matrix, graph.max_degree, N);
		std::cout << "Before Expanding opacity 6";

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
	unsigned long long t1_time, t2_time, t3_time;
	t1_time = dtime_usec(0);
	std::cout << "Converting to " << std::endl;
	int l_value = std::atoi(argv[1]);
	std::cout << l_value << " L value" << std::endl;
	graph.L_VALUE = l_value;
	graph.init_test_graph(); // Reading graph from the file in COO format
	std::cout << "Init ready " << std::endl;
	graph.convert_to_CSR();
	ordering_function(graph);
	t1_time = dtime_usec(t1_time);

//	UINT wTimerRes = 0;
//	bool init = InitMMTimer(wTimerRes);
//	DWORD startTime = timeGetTime();
	t2_time = dtime_usec(0);
	form_full_level_graph(graph);
	t2_time = dtime_usec(t2_time);
	t3_time = dtime_usec(0);
	calc_L_opacity(graph);
	t3_time = dtime_usec(t3_time);

//	unsigned int endTime = timeGetTime();
//	unsigned int gpu_time = unsigned int(endTime - startTime);
//	printf("GPU Timing(including all device-host, host-device copies, device allocations and freeing of device memory): %dms\n\n", gpu_time);
//	DestroyMMTimer(wTimerRes, init);


	graph.print_csr_graph();
	graph.print_opacity_matrix();

	std::cout << "Graph preprocessing " << t1_time/(float)USECPSEC << " Graph L - apsp " << t2_time/(float)USECPSEC
	<< " Opacity matrix calculation " << t3_time/(float)USECPSEC <<std::endl;

	return 0;
}
