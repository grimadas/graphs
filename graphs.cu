/*
 	Main application
	Author : Bulat, 2015

*/
//#include "apsp.cuh"
#include "graph.cuh"


void order_graph(Graph& graph)
{
	int	gridsize = (graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE;
	sorter<<<gridsize, BLOCK_SIZE>>>(graph.full_edge_array, graph.full_vertex_array, 0, graph.number_of_vertex);
	cudaDeviceSynchronize();
}



/********************************************************************
*	L opacity matrix calculation
*
*	 Input: Graph with non empty full_edge_array and full_vertex_array
*	 Output: opacity_matrix
*
*********************************************************************/
void calc_L_opacity(Graph &graph,
	device_ptr<int> from_vertex_array,
	device_ptr<vertex> expanded_array,
	device_ptr<int>  new_degrees,int total_size)
{
		// Forming indexes (from vertex)
		int N = total_size;
		//
		//	Expanding indexes. Finding break points
		//	Example: 0 1 2 3 4 .. 20 => 0 0 1 0 0 0 1 ...
		//

		device_ptr<int> from_degree = device_malloc<int>(N);

		//
		//	Transforming from indexes into degrees:
		//	Example:  0 0 1 1 1 1 2 2 2.. => 2 2 4 4 4 4 ...
		//

		transform(
			device,
			make_permutation_iterator(graph.initial_vertex_degrees, from_vertex_array),
			make_permutation_iterator(graph.initial_vertex_degrees, from_vertex_array + N),
			from_degree, identity<vertex>());

		device_ptr<int> to_degree = device_malloc<int>(N);

		transform(
			device,
			make_permutation_iterator(graph.initial_vertex_degrees, expanded_array),
			make_permutation_iterator(graph.initial_vertex_degrees, expanded_array + N),
			to_degree, identity<vertex>());


		/*
		*  Find max and min in zip iterator of to - from pairs
		*/


		transform(
			device,
			make_zip_iterator(make_tuple(from_degree, to_degree)),
			make_zip_iterator(make_tuple(from_degree + N, to_degree + N)),
			make_zip_iterator(make_tuple(from_degree, to_degree)),
			min_max_transform());

		//
		// 	Opacity  matrix forming. Now it is n^ 2 memory TODO: IN PARALLEL using cuda kernel
		// 	Assumptions !!: Not optimum for undericted (div 2).
		// 	Problem with same degree. Example: {4 = > 4} - must count only degree of one.
		//

		device_ptr<int> opacity_index = device_malloc<int>(N);
		if (debug)
			std::cout << "Opacity  started ok";

		int gridsize =(N + BLOCK_SIZE - 1) / BLOCK_SIZE;
		opacity_former<<<gridsize, BLOCK_SIZE>>>(from_degree,
																						 to_degree,
																						 opacity_index,
																						 graph.degree_count,
																						 graph.opacity_matrix,
																						 graph.lessL_matrix,
																						 from_vertex_array
																						 expanded_array,
																						 new_degrees,
																						 graph.max_degree,
																						 graph.threshold,
																						 N);


		cudaDeviceSynchronize();

		// Remove all edges more than threshold
		if (debug)
		{
			int num = N;
			domain temp_array = new vertex[num];
			std::cout << "After opacity removal " << std::endl;
			std::cout << "Opacity index: ";
			copy(opacity_index, opacity_index + num, temp_array);
			for(int i=0; i < num; i++)
				std::cout << temp_array[i] << " ";
			std::cout << std::endl;
			std::cout << " Expanded : ";
			copy(expanded_array, expanded_array + num, temp_array);
			for(int i=0; i < num; i++)
				std::cout << temp_array[i] << " ";
			std::cout << std::endl;
			delete temp_array;

		}



	}







/********************************************************************
*	L opacity matrix calculation
*
*	 Input: Graph with non empty full_edge_array and full_vertex_array
*	 Output: opacity_matrix
*
*********************************************************************/
void calc_L_opacity(Graph &graph,int start_level ,int end_level)
{
	for (int i = start_level; i <= end_level; i++)
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

		//
		//	Expanding indexes. Finding break points
		//	Example: 0 1 2 3 4 .. 20 => 0 0 1 0 0 0 1 ...
		//
		copy(device, graph.full_vertex_array + (i-1)*graph.number_of_vertex, graph.full_vertex_array + (i)*graph.number_of_vertex , from_vertex);
		adjacent_difference(device, from_vertex, from_vertex + graph.number_of_vertex, from_vertex);

		//
		//	Transorming into indexes:
		//	Example:	0 0 1 0 0 0 1 => 0 0 1 1 1 1 2 2 2 ..
		//
		from_vertex[0] = from_vertex[0] - starting_point;
	 	expand(from_vertex, from_vertex + graph.number_of_vertex, make_counting_iterator<vertex>(0), from);


		device_ptr<int> from_degree = device_malloc<int>(N);

		//
		//	Transforming from indexes into degrees:
		//	Example:  0 0 1 1 1 1 2 2 2.. => 2 2 4 4 4 4 ...
		//

		transform(
			device,
			make_permutation_iterator(graph.initial_vertex_degrees, from),
			make_permutation_iterator(graph.initial_vertex_degrees, from + N),
			from_degree, identity<vertex>());

		/*
		*	To vector. Transform edge list into degree list =>  similar techno
		*
		*/


		device_ptr<vertex> to = device_malloc<vertex>(N+1);
		device_ptr<int> to_degree = device_malloc<int>(N+1);
		//	auto iter_begin = make_transform_iterator(full_edge_array.begin(), minus_one());
		//	auto iter_end =   make_transform_iterator(full_edge_array.begin() + N, minus_one());
		copy(device, graph.full_edge_array + starting_point, graph.full_edge_array + ending_point, to);

		transform(
			device,
			make_permutation_iterator(graph.initial_vertex_degrees, to),
			make_permutation_iterator(graph.initial_vertex_degrees, to + N),
			to_degree, identity<vertex>());


		/*
		*  Find max and min in zip iterator of to - from pairs
		*/


		transform(
			device,
			make_zip_iterator(make_tuple(from_degree, to_degree)),
			make_zip_iterator(make_tuple(from_degree + N, to_degree + N)),
			make_zip_iterator(make_tuple(from_degree, to_degree)),
			min_max_transform());





		//
		// 	Opacity  matrix forming. Now it is n^ 2 memory TODO: IN PARALLEL using cuda kernel
		// 	Assumptions !!: Not optimum for undericted (div 2).
		// 	Problem with same degree. Example: {4 = > 4} - must count only degree of one.
		//

		device_ptr<int> opacity_index = device_malloc<int>(N);
		if (debug)
			std::cout << "Opacity  started ok";

		int gridsize =(N + BLOCK_SIZE - 1) / BLOCK_SIZE;
		opacity_former<<<gridsize, BLOCK_SIZE>>>(from_degree,
																						 to_degree,
																						 opacity_index,
																						 graph.degree_count,
																						 graph.opacity_matrix,
																						 graph.lessL_matrix,
																						 graph.max_degree, N);


		cudaDeviceSynchronize();
		// Initial remove

		if(i==1)
		{
			transform(
				device,
				make_zip_iterator(make_tuple(from, to)),
				make_zip_iterator(make_tuple(from+N, to+N)),
				make_zip_iterator(make_tuple(from,to)),
				min_max_transform_int());

			stable_sort_by_key(device,
				from, from + N,
			make_zip_iterator(make_tuple(to, opacity_index)));
			int max_value = from[N-1] + 1;
			if(debug)
				std::cout << "Size is " << max_value;


			gridsize = (max_value + BLOCK_SIZE - 1) / BLOCK_SIZE;
			fill(device,
				make_zip_iterator( make_tuple(graph.opacity_index,
					 														graph.from_array, graph.to_array)),
				make_zip_iterator( make_tuple(graph.opacity_index + N,
																			graph.from_array + N, graph.to_array + N)),
				 make_tuple(-1, -1, -1));

			coo_unifier<<<1, max_value>>>(
																				from,
																		  	to,
																		  	opacity_index,
																		    graph.from_array,
																		    graph.to_array,
																		    graph.opacity_index,
																				N,
																		    max_value);
			cudaDeviceSynchronize();

			remove(device,
				make_zip_iterator( make_tuple(graph.opacity_index,
																			graph.from_array, graph.to_array)),
			  make_zip_iterator( make_tuple(graph.opacity_index + N,
																			graph.from_array + N, graph.to_array + N)),
				make_tuple(-1, -1, -1));
		}

		if(debug)
			std::cout << "Opacity finished ok";


	}
}


void edge_removal(Graph& graph)
{
	// Degree value in opacity matrix
		counting_iterator<vertex> degree_value(0);
	int N = graph.max_degree;
	if (debug)
	{
		std::cout << "Edge removal start" << std::endl;
	}

	// Form possible candidates for removal or insertion

	int possible_size = N*N;
	device_ptr<int> remove_edges_count = device_malloc<int>(N*N);
	device_ptr<int> opacity_index = device_malloc<int>(N*N);
	if (debug)
	{
		std::cout << "Sorting in edge removal " << std::endl;
	}
	domain a = new vertex[graph.number_of_edges];
	copy(graph.opacity_index, graph.opacity_index + graph.number_of_edges, a);
	std::cout << "Opacity index: ";
	for (int i =0; i < graph.number_of_edges; i++)
	{
			std::cout << a[i];
	}
	std::cout << std::endl;

	stable_sort_by_key(device,
		graph.opacity_index, graph.opacity_index + graph.number_of_edges,
	make_zip_iterator(make_tuple(graph.from_array, graph.to_array)));


	fill(device,
		make_zip_iterator(
			make_tuple(remove_edges_count,
								 opacity_index)),
		make_zip_iterator(
			make_tuple(remove_edges_count + possible_size,
									 opacity_index + possible_size)),
		make_tuple(-1,-1));

	transform_if(
		device,
		make_zip_iterator(
			make_tuple(graph.opacity_matrix,
			    			 graph.lessL_matrix,
								 degree_value)),
		make_zip_iterator(
			make_tuple(graph.opacity_matrix + N*N,
								 graph.lessL_matrix + N*N,
								 degree_value + N*N)),
		make_zip_iterator(
			make_tuple(opacity_index,
								 remove_edges_count)),
		form_remove_candidates(graph.threshold),
		more_than_threshold(graph.threshold));

	// Calculate vertex size :
	device_ptr<int> last_el=
	remove(opacity_index,
		opacity_index + possible_size,
		-1);

	remove(remove_edges_count,
		remove_edges_count + possible_size,
		-1);

	int real_size = distance(opacity_index, last_el);

	domain pr = new vertex[real_size];
	domain pr2 = new vertex[real_size];
	copy(opacity_index, opacity_index+ real_size, pr);
	copy(remove_edges_count, remove_edges_count+ real_size, pr2);

	if (debug)
	{
		for(int i=0; i < real_size; i++)
		{
			std::cout << "Opacity " << pr[i] << " " << " remove" << pr2[i] << std::endl;
		}
	}

	device_ptr<int> index_array = device_malloc<int>(graph.number_of_edges);
	fill(device, index_array, index_array + graph.number_of_edges, -1);
	if (debug)
	{
		std::cout << "Edge removal kernel" << std::endl;
	}

	edge_remover<<<1, real_size>>>(
		opacity_index,
		remove_edges_count,
		graph.opacity_index,
		index_array,
		graph.number_of_edges,
		real_size);

	device_ptr<int> index_size_ptr =
	remove(device, index_array, index_array+graph.number_of_edges, -1);
	int grid_size = distance(index_array, index_size_ptr);
	if (debug)
	{
		int size = grid_size;
		domain temp_array = new  vertex[size];
		copy(index_array, index_array + size, temp_array);
		std::cout << " \n Should be removed: ";
		for(int i = 0; i < size; i++)
		{
			std::cout << " " << temp_array[i];
		}
		std::cout << "\n";
	}

	remove_by_index<<< 1, grid_size>>> (index_array,
											graph.degree_count,graph.opacity_matrix,graph.lessL_matrix,
											graph.max_degree,
											graph.from_array, graph.to_array, graph.opacity_index);

	device_ptr<vertex> last_elto  =
	remove(device,
				graph.from_array, graph.from_array + graph.number_of_edges, -1);

	remove(device,
				make_zip_iterator(make_tuple(graph.to_array, graph.opacity_index)),
				make_zip_iterator(make_tuple(graph.to_array + graph.number_of_edges,
																		graph.opacity_index + graph.number_of_edges)),
			  make_tuple(-1,-1));
	graph.number_of_edges = distance(graph.from_array,last_elto);
	if (debug)
	{
		std::cout << "Edge removal done" << std::endl;
	}

}


/*
*	By finding shortest paths, form to L_VALUE level
*/
void form_full_level_graph(Graph &graph)
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

			device_ptr<vertex> inter_vertex = device_malloc<vertex>(number_edges_to_process);
			device_ptr<vertex> temp_to =  device_malloc<vertex>(number_edges_to_process);
			device_ptr<vertex> temp_from =  device_malloc<vertex>(number_edges_to_process);

			/* Form temp to an temp from vector from edge arrray */
			copy(device, graph.full_edge_array + starting_point,
									 graph.full_edge_array + ending_point,
									 inter_vertex);

			copy(device,
				graph.full_edge_array + starting_point,
				 graph.full_edge_array + ending_point,
				 temp_from);

			 std::cout << "Ok in 1" << std::endl;
				//
				transform(
					device,
					temp_from,
					temp_from + number_edges_to_process,
					temp_from,
					previous_el(current_level * graph.number_of_vertex + 1));

			/* Store begining and ending */
			copy(
				device,
				make_permutation_iterator(graph.full_vertex_array, inter_vertex),
				make_permutation_iterator(graph.full_vertex_array, inter_vertex + number_edges_to_process),
				temp_to);

			copy(
				device,
				make_permutation_iterator(graph.full_vertex_array, temp_from),
				make_permutation_iterator(graph.full_vertex_array,
				temp_from + number_edges_to_process), temp_from);
				// std::cout << "Ok in 2" << std::endl;

			/*
			*	Array of vertex, from which we will expand. Proces vertex
			*/
			device_ptr<vertex> temp_vertexes =  device_malloc<vertex>(graph.number_of_vertex);
			device_ptr<vertex> process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
			copy(device,
				   graph.full_vertex_array + (current_level-1)*graph.number_of_vertex,
					 graph.full_vertex_array + (current_level)*graph.number_of_vertex ,
					temp_vertexes);

			adjacent_difference(device, temp_vertexes, temp_vertexes + graph.number_of_vertex, temp_vertexes);
			temp_vertexes[0] = temp_vertexes[0] - starting_point;

			expand(temp_vertexes, temp_vertexes + graph.number_of_vertex,
						 make_counting_iterator<vertex>(0), process_vetxes);
			device_free(temp_vertexes);
			// std::cout << "Ok in 3" << std::endl;
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
				// std::cout << "Ok in 4" << std::endl;
				/*
					Forming offset array from process number, step 2:
					2 4 4 => 2 6 10
					*/
				inclusive_scan(device, position_in_array,
							 		    position_in_array + number_edges_to_process,
									    position_in_array);

		//	domain prints = new vertex[number_edges_to_process];
		//	copy(position_in_array, position_in_array + number_edges_to_process, prints);
		//	std::cout << "Position In array:  " << std::endl;
		//	for(int i=0; i < number_edges_to_process; i++)
		//	{
		//		std::cout << prints[i]<<" ";
		//	}

			//// std::cout << "Ok in 5" << std::endl;
			int full_size = position_in_array[number_edges_to_process - 1];
			device_ptr<vertex> expanded_array = device_malloc<vertex>(full_size);
			device_ptr<vertex> from_vertex_array = device_malloc<vertex>(full_size);
			device_ptr<vertex> inter_vertex_array = device_malloc<vertex>(full_size);

			//fill(device, expanded_array, expanded_array + full_size, -1);
			// process number contains the maximum needed memory to store if all vertexes are unique

			//fill(device, from_vertex_array, from_vertex_array + full_size, -1);



			fill(device,
				make_zip_iterator(
					make_tuple(expanded_array, from_vertex_array, inter_vertex_array)),
				make_zip_iterator(
					make_tuple(expanded_array + full_size,
										 from_vertex_array + full_size,
										 inter_vertex_array + full_size)),
					make_tuple(-1,-1,-1));


			int prev_max_position = full_size;

			//  Expand array on one level
			//	Can contain non unique values

			device_ptr<vertex> positions_vertex_current_level = device_malloc<vertex>(graph.number_of_vertex);
			fill(device, positions_vertex_current_level, positions_vertex_current_level + graph.number_of_vertex, 0);

			dim3 grid_size((number_edges_to_process + BLOCK_SIZE - 1) / BLOCK_SIZE);

		  //if(debug)
				std::cout << "Ok in 6 " << "Number of edges to process is " << grid_size.x <<std::endl;
			// Expand here and put to expanded_array


			expander<<< grid_size, BLOCK_SIZE >>>(
				process_vetxes, inter_vertex,
				temp_from, temp_to,
				graph.full_vertex_array, graph.full_edge_array,
				position_in_array,
				expanded_array,
				from_vertex_array,
				inter_vertex_array,
				number_edges_to_process,
				graph.number_of_vertex,
				current_level,
				positions_vertex_current_level);

			cudaDeviceSynchronize();
		//	if (debug)
			  std::cout << "Ok in 9" << std::endl;
			device_free(temp_from);
			device_free(temp_to);
			device_free(process_vetxes);
			cudaFree(0);


			//
			//	Remove empty, non used data
			//

			remove(device,
				make_zip_iterator(
					make_tuple(expanded_array,
										 from_vertex_array,
										 inter_vertex_array)),
				make_zip_iterator(
					make_tuple(expanded_array + full_size,
				   					 from_vertex_array + full_size,
										 inter_vertex_array + full_size)),
				make_tuple(-1,-1,-1));

			if (debug)
			{
				int num = full_size;
				domain temp_array = new vertex[num];
				std::cout << "Expanded candidates " << std::endl;
				copy(from_vertex_array, from_vertex_array + num, temp_array);
				for(int i=0; i < num; i++)
					std::cout << temp_array[i] << " ";
				std::cout << std::endl;

				copy(inter_vertex_array, inter_vertex_array + num, temp_array);
				for(int i=0; i < num; i++)
					std::cout << temp_array[i] << " ";
				std::cout << std::endl;

				copy(expanded_array, expanded_array + num, temp_array);
				for(int i=0; i < num; i++)
					std::cout << temp_array[i] << " ";
				std::cout << std::endl;


				delete temp_array;

			}

		//	if(debug)
				std::cout << "Ok in 10" << std::endl;


			// STEP 2: Forming offset for unifier function
			inclusive_scan(device, positions_vertex_current_level,
									           positions_vertex_current_level + graph.number_of_vertex,
									           positions_vertex_current_level);

//			if(debug)
				std::cout << "Position vertex " << positions_vertex_current_level[graph.number_of_vertex - 1] << std::endl;
			int ze_size = positions_vertex_current_level[graph.number_of_vertex - 1];
			device_ptr<vertex> vertex_ending_offsets = device_malloc<vertex>(graph.number_of_vertex);
		  //fill(device, vertex_ending_offsets, vertex_ending_offsets + graph.number_of_vertex, 0);
			//	positions_vertex_current_level[graph.number_of_vertex - 1]+1)
			device_ptr<vertex> expanded_array_copy = device_malloc<vertex>(ze_size);
			device_ptr<vertex> temp_expanded = device_malloc<vertex>(ze_size);


			copy(device, expanded_array, expanded_array + ze_size, temp_expanded);

			device_ptr<vertex> from_vertex_array_copy = device_malloc<vertex>(ze_size);
			device_ptr<vertex> inter_vertex_array_copy = device_malloc<vertex>(ze_size);

			device_ptr<int> indices = device_malloc<int>(ze_size);

			// Remove all dublicates. Calculate the opacity, remove if needed


		//	if (debug)
				std::cout << "Starting unification " << std::endl;
			int grid_size_uni = ((graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE);
			unifier <<<grid_size_uni, BLOCK_SIZE >>>(	expanded_array,
																								expanded_array_copy,
																								temp_expanded,
																							 	from_vertex_array,
																								from_vertex_array_copy,
																								inter_vertex_array,
																								inter_vertex_array_copy,
																								indices,
																								positions_vertex_current_level,
																								vertex_ending_offsets,
																								graph.number_of_vertex);

			cudaDeviceSynchronize();
			// Now we have unique expanded_copy_array
			// vertex_ending_offsets[graph.number_of_vertex - 1];
			int proceed_size =  // vertex_ending_offsets[graph.number_of_vertex - 1];
			reduce(device, vertex_ending_offsets,
										 vertex_ending_offsets + graph.number_of_vertex);

			if (debug)
				std::cout << "proceed_size " << proceed_size << std::endl;
			// Remove shlack

			if(debug)
				std::cout << "unification done " << std::endl;
			copy(device, vertex_ending_offsets,
				vertex_ending_offsets + graph.number_of_vertex,
				graph.full_vertex_array + current_level * graph.number_of_vertex);


			inclusive_scan(device,
				graph.full_vertex_array + current_level * graph.number_of_vertex - 1,
				graph.full_vertex_array + (current_level + 1) * graph.number_of_vertex,
				graph.full_vertex_array + current_level * graph.number_of_vertex - 1);

			if(debug)
					std::cout << "Starting copier " << std::endl;

			device_free(temp_expanded);
			temp_expanded = device_malloc<vertex>(proceed_size);
			device_ptr<vertex> temp_from_array = device_malloc<vertex>(proceed_size);



			 temp_edge_copier<<<grid_size_uni, BLOCK_SIZE>>>
		 			(expanded_array, from_vertex_array,
			  	positions_vertex_current_level, vertex_ending_offsets,
					graph.full_vertex_array,
			  	temp_from_array, temp_expanded,
			  	graph.number_of_vertex, current_level);

			cudaDeviceSynchronize();
			calc_L_opacity(graph, temp_from_array, temp_expanded,
														vertex_ending_offsets, proceed_size);
			// Renew vertex_offset array :
			copy(device, vertex_ending_offsets, vertex_ending_offsets + graph.number_of_vertex,
									 graph.full_vertex_array + (current_level)*graph.number_of_vertex);
			inclusive_scan(device,
							graph.full_vertex_array + (current_level)*graph.number_of_vertex - 1,
							graph.full_vertex_array + (current_level+1)*graph.number_of_vertex,
							graph.full_vertex_array + (current_level)*graph.number_of_vertex - 1);
			copy(device, vertex_ending_offsets,
									 vertex_ending_offsets + graph.number_of_vertex,
									 graph.real_vertex_degrees + (current_level) * graph.number_of_vertex);
			// Copy all non negative values to full_edge_array
			remove_copy_if(device,
										temp_expanded,
										temp_expanded + proceed_size,
										graph.full_edge_array +
										graph.full_vertex_array[current_level*graph.number_of_vertex - 1],
										negative());
			// Remove all positive values
			zip_iterator< tuple < device_ptr<int>,  device_ptr<int> > > new_ending =
			remove_if(device,
								make_zip_iterator(make_tuple(temp_expanded, temp_from_array)),
								make_zip_iterator(make_tuple(temp_expanded + proceed_size, temp_from_array + proceed_size)),
								non_negative());
			device_ptr<int> temp_ptr = get<0>(new_ending.get_iterator_tuple());
			int new_size = distance(temp_expanded, temp_ptr);
			// Translate removing array back to positive
			transform(device, temp_expanded, temp_expanded + new_size, to_positive_plus_one());
			// Search position in copied arrays
			device_ptr<unsigned int> start_points = device_malloc<int>(new_size);
			lower_bound(   from_vertex_array_copy, from_vertex_array_copy + ze_size,
      							 temp_from_array,  temp_from_array + new_size,
      						   start_points);

			device_ptr<unsigned int> end_points = device_malloc<int>(new_size);
			upper_bound(  from_vertex_array_copy, from_vertex_array_copy + ze_size,
      							temp_from_array,  temp_from_array + new_size,
      						  end_points);
			// STEP 6. Preparing the removal candidates
			dim3 grid_size_remove1(new_size );
			form_remove_array<<<grid_size_remove1, blocksize>>>(
											start_points,
											end_points,
											temp_expanded,
											from_vertex_array_copy,
											inter_vertex_array_copy,
											expanded_array_copy,
											int N);
			new_ending =
			remove_if(device,
								make_zip_iterator(make_tuple(expanded_array_copy, inter_vertex_array_copy)),
								make_zip_iterator(make_tuple(expanded_array_copy + ze_size, inter_vertex_array_copy + ze_size)),
								non_negative());
			temp_ptr = get<0>(new_ending.get_iterator_tuple());
			int remove_array_size = distance(expanded_array_copy, temp_ptr);
				transform(device, expanded_array_copy,
													expanded_array_copy + remove_array_size,
													to_positive_plus_one());
			// STEP 7. Remove from inter_vertex all expanded arrays

			dim3 grid_size_remove1();
			first_level_remove_via_sort<<<grid_size_remove1, blocksize>>>(
											graph.full_vertex_array,
											graph.full_edge_array,
											from_vertex_array_copy,
											inter_vertex_array_copy,
											expanded_array_copy,
											graph.real_vertex_degrees,
											remove_array_size);

			cudaDeviceSynchronize();
			// Calculate degrees via adjacent difference
			// Copy first level
			// STEP 8. Removal function
			for(int level = 1; level < current_level; level++)
			{
						int starter = 0;
						if (level!=1)
						{
							starter = graph.full_vertex_array[(level-1)*graph.number_of_vertex - 1];
						}
						int ender = graph.full_vertex_array[level*graph.number_of_vertex - 1];
						device_ptr<int> indexes = device_malloc<int>(ender - starter);
						sequence(device, indexes, indexes + ender - starter);
						// Now we store the from vertex
						transform(device, indexes, indexes + ender - starter, indexes, to_from_vertex(graph.full_vertex_array));

						device_ptr<int> values = device_malloc<int>(ender - starter);
						copy(device, graph.full_edge_array + starter, graph.full_edge_array + ender, values);
						sort_by_key(device, values, values + ender - starter,indexes);

						device_ptr<unsigned int> temp_start_points = device_malloc<int>(remove_array_size);
						lower_bound(   values, values + ender - starter,
			      							 inter_vertex_array_copy,  inter_vertex_array_copy + remove_array_size,
			      						   temp_start_points);

						device_ptr<unsigned int> temp_end_points = device_malloc<int>(remove_array_size);
						upper_bound(  values, values + ender - starter,
			      							inter_vertex_array_copy,  inter_vertex_array_copy + remove_array_size,
			      						  temp_end_points);

						int blocksize = std::min(new_size,BLOCK_SIZE);
						dim3 grid_size_remove((new_size + blocksize - 1) / blocksize);
						edge_remove_by_vertex<<<grid_size_remove, blocksize>>>
																							(	temp_start_points,
																								temp_end_points,
																								expanded_array_copy,
																								temp_expanded, temp_from_array,
																							 from_vertex_array_copy,
																							 inter_vertex_array_copy,
																							 expanded_array_copy,
																							 graph.full_vertex_array,
																							 graph.full_edge_array,
																							 cur_vert,
																							 new_size);
							cudaDeviceSynchronize();
							sorter<<<grid_size_uni, BLOCK_SIZE>>>(graph.full_edge_array,
													           								graph.full_vertex_array,
																										level,
							                      								graph.number_of_vertex);
							}


					}

			}


			/***************************************************
			*			Break point
			*			Recalculate at each level, or apsp
			****************************************************/

				if (debug)
				{
					int num = graph.number_of_vertex;
					int new_num = proceed_size;
					domain temp_array = new vertex[graph.number_of_vertex];
					std::cout << "After unification " << std::endl;
					std::cout << "Vertex pos: ";
					copy(positions_vertex_current_level, positions_vertex_current_level + num, temp_array);
					for(int i=0; i < num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;
					std::cout << "Final offsets : ";
					copy(vertex_ending_offsets, vertex_ending_offsets + num, temp_array);
					for(int i=0; i < num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;
					delete temp_array;

					temp_array = new vertex[new_num];
					std::cout << "From vertex array copy: ";
					copy(from_vertex_array, from_vertex_array + new_num, temp_array);
					for(int i=0; i < new_num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;

					std::cout << "Inter array copy: ";
					copy(inter_vertex_array_copy, inter_vertex_array_copy + new_num, temp_array);
					for(int i=0; i < new_num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;

					std::cout << "Expanded array copy: ";
					copy(expanded_array, expanded_array + new_num, temp_array);
					for(int i=0; i < new_num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;

					std::cout << "Temporal expanded: ";
					copy(temp_expanded, temp_expanded + new_num, temp_array);
					for(int i=0; i < new_num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;


					std::cout << "Temporal from: ";
					copy(temp_from_array, temp_from_array + new_num, temp_array);
					for(int i=0; i < new_num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;

					std::cout << "Indices : ";
					copy(indices, indices + new_num, temp_array);
					for(int i=0; i < new_num; i++)
						std::cout << temp_array[i] << " ";
					std::cout << std::endl;
					delete temp_array;
				}

			 std::cout << "Ok in 12" << std::endl;

			 std::cout << "Ok in 13" << std::endl;
		//	calc_L_opacity(graph, current_level + 1, current_level + 1);

			device_free(positions_vertex_current_level);
			device_free(vertex_ending_offsets);
			device_free(position_in_array);
			device_free(expanded_array_copy);

			//cudaFree(0);
			 std::cout << "Ok in 14" << std::endl;


	}
}


int main(int argc, char* argv[])
{

	Graph graph;
	unsigned long long read_time,t1_time, t2_time, t3_time,t4_time;

	int l_value = std::atoi(argv[1]);
	float threshold = std::atof(argv[2]);
	graph.L_VALUE = l_value;
	graph.threshold = threshold;

	read_time = dtime_usec(0);
	graph.init_test_graph(); // Reading graph from the file in COO format
	read_time = dtime_usec(read_time);
	std::cout << "Graph is readed "<< read_time/(float)USECPSEC  << " seconds" << std::endl;

	size_t free_mem;
	size_t total_mem;
	cudaMemGetInfo(&free_mem,&total_mem);
	std::cout << "Memory 1: " << (long)free_mem/1024/1024 << " Total " << (long)total_mem/1024/1024 << std::endl;


	// Initial convertion
	t1_time = dtime_usec(0);
	graph.convert_to_CSR(true, true);

	order_graph(graph);

	cudaMemGetInfo(&free_mem,&total_mem);
	std::cout << "Memory 2: " << (long)free_mem/1024/1024 << " Total " << (long)total_mem/1024/1024 << std::endl;

	calc_L_opacity(graph, 1, 1);
	graph.print_opacity_matrix();

	// Something must be removed
	if (graph.opacity_index[0] != -1)
	{
		edge_removal(graph);
		cudaMemGetInfo(&free_mem,&total_mem);
		std::cout << "Memory 3: " << (long)free_mem/1024/1024 << " Total " << (long)total_mem/1024/1024 << std::endl;
		graph.convert_to_CSR(false, false);
		order_graph(graph);
	}
	t1_time = dtime_usec(t1_time);
	std::cout << "Graph is ready for apsp. Initial removal done in  "<< t1_time/(float)USECPSEC  << " seconds" << std::endl;

	cudaMemGetInfo(&free_mem,&total_mem);
	std::cout << "Memory 4: " << (long)free_mem/1024/1024 << " Total " << (long)total_mem/1024/1024 << std::endl;

	t2_time = dtime_usec(0);
	form_full_level_graph(graph);
	t2_time = dtime_usec(t2_time);
	std::cout << " \n Graph L - apsp time : " << t2_time/(float)USECPSEC << std::endl;
	t3_time = dtime_usec(0);

	t3_time = dtime_usec(t3_time);

	cudaMemGetInfo(&free_mem,&total_mem);
	std::cout << "Memory 5: " << (long)free_mem/1024/1024 << " Total " << (long)total_mem/1024/1024 << std::endl;

	std::cout << " Opacity matrix calculation time " << t3_time/(float)USECPSEC <<std::endl;
	t4_time = dtime_usec(0);
//	edge_removal(graph);
	t4_time = dtime_usec(t4_time);
	std::cout << " Edge removal function time" << t4_time/(float)USECPSEC <<std::endl;

	graph.print_csr_graph();
	graph.print_opacity_matrix();

	return 0;
}
