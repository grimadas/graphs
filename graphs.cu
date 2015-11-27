/*
 	Main application
	Author : Bulat, 2015

*/
//#include "apsp.cuh"
#include "graph.cuh"


/**
*		Filling the arrays that will be further expanded.
*		PostCondition: filled from, to, inter
*/
void form_copy_arrays(Graph& graph,
											device_ptr<vertex> temp_from,
											device_ptr<vertex> inter_vertex,
											device_ptr<vertex> temp_to,
											int current_level,
											int starting_point,
											int ending_point)
{
	/* Form temp to an temp from vector from edge arrray */
	int number_edges_to_process = ending_point - starting_point;
	copy(device, graph.full_edge_array + starting_point,
							 graph.full_edge_array + ending_point,
							 inter_vertex);

	copy(device,
		graph.full_edge_array + starting_point,
		 graph.full_edge_array + ending_point,
		 temp_from);

	 std::cout << "Forming copy arrays " << std::endl;
		//Form from vertex array
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
	 std::cout << "Copy candidates are formed " << std::endl;
}


/**
*	  Caluclate position in new array for each vertex
*		Postconditon: position in array filled
*/
void calc_position_in_new_array(device_ptr<int> temp_from,
																device_ptr<int> temp_to,
																device_ptr<int> position_in_array,
																int number_edges_to_process)
{
		transform(
			device,
			make_zip_iterator(make_tuple(temp_from, temp_to)),
			make_zip_iterator(make_tuple(temp_from + number_edges_to_process,
												temp_to + number_edges_to_process)),
			position_in_array,
			counter());

		 /*
			Forming offset array from process number, step 2:
			2 4 4 => 2 6 10
			*/
		inclusive_scan(device, position_in_array,
									position_in_array + number_edges_to_process,
									position_in_array);
}


/**
*	Renew vertex offset in graph in current level
*
*/
void renew_graph_offset(Graph& graph,
												device_ptr<int> new_offset_array,
												int current_level)
{
			// Renew vertex_offset array :
  		copy(device, new_offset_array, new_offset_array + graph.number_of_vertex,
						 graph.full_vertex_array + (current_level)*graph.number_of_vertex);
			inclusive_scan(device,
					graph.full_vertex_array + (current_level)*graph.number_of_vertex - 1,
					graph.full_vertex_array + (current_level+1)*graph.number_of_vertex,
	  			graph.full_vertex_array + (current_level)*graph.number_of_vertex - 1);

			copy(device, new_offset_array,
					 new_offset_array + graph.number_of_vertex,
					 graph.real_vertex_degrees + (current_level) * graph.number_of_vertex);
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
																						 from_vertex_array,
																						 expanded_array,
																						 new_degrees,
																						 graph.number_of_vertex + 1,
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


/**
*		Remove edges that are non_violating rule.
*		All violated values are negative
*		Return new size of an array
*/
int remove_non_violated_edges(
			device_ptr<int> temp_expanded,
			device_ptr<int> temp_from_array,
			int number_of_vertex,
			int true_unique_size)
{
	zip_iterator< tuple < device_ptr<int>,  device_ptr<int> > > new_ending =
	remove_if(device,
						make_zip_iterator(make_tuple(temp_expanded, temp_from_array)),
						make_zip_iterator(make_tuple(temp_expanded + true_unique_size,
																				 temp_from_array + true_unique_size)),
						non_negative());
	device_ptr<int> temp_ptr = get<0>(new_ending.get_iterator_tuple());
	int new_size = distance(temp_expanded, temp_ptr);
	// Translate removing array back to positive
	transform(device, temp_expanded,
										temp_expanded + new_size,
										temp_expanded,
										to_positive(number_of_vertex+1));
	return new_size;
}


/**
*		Generate removing candiadets using inter_vertex_array_copy
*/
int transform_removing_edges(
													 device_ptr<vertex> inter_vertex_array_copy,
													 device_ptr<vertex> expanded_array_copy,
													 device_ptr<int> positions_vertex_current_level,
													 int true_size,
													 device_ptr<vertex> temp_from_array,
													 device_ptr<vertex> temp_expanded,
													 int number_of_vertex,
													 int true_unique_size)
{


	int new_size = remove_non_violated_edges(temp_expanded, temp_from_array,
																						number_of_vertex, true_unique_size);
	if(debug)
	{
	 	std::cout << " Temp expanded ";
		for(int i=0; i< new_size; i++)
			std::cout << temp_expanded[i] << " ";
		std::cout << std::endl;

		std::cout << " Temp from ";
		for(int i=0; i< new_size; i++)
			std::cout << temp_from_array[i] << " ";
		std::cout << std::endl;
	}

	int grid_size_remove1 = (new_size + BLOCK_SIZE - 1)/BLOCK_SIZE;
	form_remove_array<<<grid_size_remove1, BLOCK_SIZE>>>(
																	expanded_array_copy,
		 															positions_vertex_current_level,
																	temp_from_array,
		 															temp_expanded,
																	number_of_vertex,
																	new_size);
		new_size =
				remove_non_violated_edges(expanded_array_copy, inter_vertex_array_copy,
			 	  												number_of_vertex,
																	true_size);
		return new_size;
}



/**
*		At current level :
*			1. copy all vertexes
*				calc from_vertex (indexes)
*			2. sort them -> values
*			3. find begin and end in values all inter_vertex_array_copy(that should be removed)
*			4.
*/
void edge_remove_step(Graph& graph,
											device_ptr<vertex> inter_vertex_array_copy,
											device_ptr<vertex> expanded_array_copy,
											int remove_array_size, int level)
{
	int starter = 0;
	int adding_offset = 0;
	if (level!=0)
	{

		adding_offset = (level)*graph.number_of_vertex - 1;
		starter = graph.full_vertex_array[adding_offset];
	}
	int ender = graph.full_vertex_array[(level+1)*graph.number_of_vertex - 1];
	device_ptr<int> indexes = device_malloc<int>(ender - starter);
	copy(device, make_counting_iterator(starter),
							 make_counting_iterator(ender),
							 indexes);
	// sequence(device, indexes, indexes + ender - starter);
	// Now we store the from vertex


	transform(device, indexes, indexes + ender - starter,
					 indexes,
					 to_from_vertex(graph.full_vertex_array + adding_offset,
					  							graph.number_of_vertex));

	device_ptr<int> values = device_malloc<int>(ender - starter);
	copy(device, graph.full_edge_array + starter, graph.full_edge_array + ender, values);
	sort_by_key(device, values, values + ender - starter,indexes);

	device_ptr<unsigned int> temp_start_points = device_malloc<unsigned int>(remove_array_size);
	lower_bound(  device, values, values + ender - starter,
								 inter_vertex_array_copy,  inter_vertex_array_copy + remove_array_size,
								 temp_start_points);

	device_ptr<unsigned int> temp_end_points = device_malloc<unsigned int>(remove_array_size);
	upper_bound( device, values, values + ender - starter,
								inter_vertex_array_copy,  inter_vertex_array_copy + remove_array_size,
								temp_end_points);
	device_ptr<int> temp_offsets = device_malloc<int>(remove_array_size);
	transform(device, make_zip_iterator(
											make_tuple(temp_start_points, temp_end_points)),
										make_zip_iterator(
											make_tuple(temp_start_points + remove_array_size,
																			temp_end_points + remove_array_size)),
										temp_offsets,
										counter());
	inclusive_scan(device, temp_offsets, temp_offsets + remove_array_size, temp_offsets);
	int total_size = temp_offsets[remove_array_size - 1];
	device_ptr<int> temp_removing_array = device_malloc<int>(total_size);
	fill(device, temp_removing_array, temp_removing_array + total_size, -1);
	dim3 grid_size_remove((remove_array_size + BLOCK_SIZE - 1) / BLOCK_SIZE);
	edge_remove_by_vertex<<<grid_size_remove, BLOCK_SIZE>>> ( expanded_array_copy,
																													 temp_start_points,
																													 temp_end_points,
																													 indexes,
																													 temp_offsets,
																													 temp_removing_array,
																													 graph.full_vertex_array,
																													 graph.full_edge_array,
																													 graph.real_vertex_degrees + graph.number_of_vertex*(level+1),
																													 graph.number_of_vertex,
																													 level, remove_array_size);
	cudaDeviceSynchronize();
	dim3 grid_size_uni((graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE);
	//sorter<<<grid_size_uni, BLOCK_SIZE>>>(graph.full_edge_array,
	//																			graph.full_vertex_array,
	//																			level + 1,
	//																			graph.number_of_vertex);
	device_free(indexes);
	device_free(values);
	device_free(temp_start_points);
	device_free(temp_end_points);
	cudaDeviceSynchronize();
}


void expand_step(Graph &graph, int current_level)
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

	form_copy_arrays(graph, temp_from, inter_vertex, temp_to,
											 current_level,
											 starting_point, ending_point);
	/*
	*	Array of vertex, from which we will expand. Proces vertex
	*/
	device_ptr<vertex> temp_vertexes =  device_malloc<vertex>(graph.number_of_vertex);
	device_ptr<vertex> process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
	copy(device,
			 graph.real_vertex_degrees + (current_level-1)*graph.number_of_vertex,
			 graph.real_vertex_degrees + (current_level)*graph.number_of_vertex,
			 temp_vertexes);

	expand(temp_vertexes, temp_vertexes + graph.number_of_vertex,
				 make_counting_iterator<vertex>(0), process_vetxes);
	device_free(temp_vertexes);

 device_ptr<vertex>  position_in_array = device_malloc<vertex>(number_edges_to_process);

 calc_position_in_new_array(temp_from, temp_to, position_in_array,
																				 number_edges_to_process);

 int full_size = position_in_array[number_edges_to_process - 1];
 device_ptr<vertex> expanded_array = device_malloc<vertex>(full_size);


 // Fill with -1 initially
 fill(device,
		 expanded_array,
		 expanded_array + full_size,
		 -1);


 device_ptr<vertex> positions_vertex_current_level = device_malloc<vertex>(graph.number_of_vertex);
 fill(device, positions_vertex_current_level, positions_vertex_current_level + graph.number_of_vertex, 0);
 // TODO Change here
 int grid_size = (number_edges_to_process + BLOCK_SIZE - 1) / BLOCK_SIZE;


 //************************ Starting expander function **********************************************/

 if(debug)
	 std::cout << "Starting main expand once kernel " << " Number of edges to process is " << grid_size <<std::endl;

 simple_expander<<< grid_size, BLOCK_SIZE >>>(
	 process_vetxes,
	 temp_from, temp_to,
	 graph.full_vertex_array, graph.full_edge_array,
	 position_in_array,
	 expanded_array,
	 number_edges_to_process,
	 graph.number_of_vertex,
	 current_level,
	 positions_vertex_current_level);

 cudaDeviceSynchronize();
 device_free(temp_from);
 device_free(temp_to);
 device_free(process_vetxes);
 device_free(inter_vertex);


 if (debug)
	 std::cout << "Expand kernel finished working " <<std::endl;


 /*
 *	Remove empty, non used data
 */
 thrust::remove(device,
	 				expanded_array,
	 				expanded_array + full_size,
					-1);


 /*********************************    Expander kernel finished ************************************/

 /************************************ Unifier kernel prepare **************************************/


 inclusive_scan(device, positions_vertex_current_level,
											positions_vertex_current_level + graph.number_of_vertex,
											positions_vertex_current_level);

 if(debug)
	 std::cout << "Position vertex " << positions_vertex_current_level[graph.number_of_vertex - 1] << std::endl;

 int true_size = positions_vertex_current_level[graph.number_of_vertex - 1];
 device_ptr<vertex> vertex_ending_offsets = device_malloc<vertex>(graph.number_of_vertex);


 // Remove all dublicates. Calculate the opacity, remove if needed
 if (debug)
	 std::cout << "Starting unification " << std::endl;

 /******************************************* Unifier kernel start *************************************/

 grid_size = ((graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE);

 simple_unifier <<<grid_size, BLOCK_SIZE >>>(	expanded_array,
																			 positions_vertex_current_level,
																			 vertex_ending_offsets,
																			 graph.number_of_vertex);

 cudaDeviceSynchronize();
 // Now we have unique expanded_array

 int proceed_size = vertex_ending_offsets[graph.number_of_vertex - 1];

 int true_unique_size =
 reduce(device, vertex_ending_offsets,
								vertex_ending_offsets + graph.number_of_vertex);

 /******************************************************* Unifier kernel ready ***********************************************/

 renew_graph_offset(graph, vertex_ending_offsets, current_level);

 if(debug)
		 std::cout << "Starting copier " << std::endl;

	 /******************************************************* Edge copy kernel start  ********************************************/

	 edge_copier<<<grid_size, BLOCK_SIZE>>>(
 		expanded_array,
 		positions_vertex_current_level,
 		vertex_ending_offsets,
 		graph.full_vertex_array,
 		graph.full_edge_array,
 		current_level,
 		graph.number_of_vertex );
	 cudaDeviceSynchronize();
	 device_free(expanded_array);
	 device_free(positions_vertex_current_level);
	 device_free(vertex_ending_offsets);


	 /****************************************************** Edge copied finished*********************************************************/

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

				 form_copy_arrays(graph, temp_from, inter_vertex, temp_to,
															current_level,
															starting_point, ending_point);
				 /*
				 *	Array of vertex, from which we will expand. Proces vertex
				 */
				 device_ptr<vertex> temp_vertexes =  device_malloc<vertex>(graph.number_of_vertex);
				 device_ptr<vertex> process_vetxes = device_malloc<vertex>(number_edges_to_process+1);
				 copy(device,
							graph.real_vertex_degrees + (current_level-1)*graph.number_of_vertex,
							graph.real_vertex_degrees + (current_level)*graph.number_of_vertex,
							temp_vertexes);

				 expand(temp_vertexes, temp_vertexes + graph.number_of_vertex,
 				 				make_counting_iterator<vertex>(0), process_vetxes);
				 device_free(temp_vertexes);

			  device_ptr<vertex>  position_in_array = device_malloc<vertex>(number_edges_to_process);

			  calc_position_in_new_array(temp_from, temp_to, position_in_array,
																								number_edges_to_process);

				int full_size = position_in_array[number_edges_to_process - 1];
				device_ptr<vertex> expanded_array = device_malloc<vertex>(full_size);
				device_ptr<vertex> inter_vertex_array = device_malloc<vertex>(full_size);

				// Fill with -1 initially
				fill(device,
						expanded_array,
					  expanded_array + full_size,
						-1);


				device_ptr<vertex> positions_vertex_current_level = device_malloc<vertex>(graph.number_of_vertex);
				fill(device, positions_vertex_current_level, positions_vertex_current_level + graph.number_of_vertex, 0);
				// TODO Change here
				int grid_size = (number_edges_to_process + BLOCK_SIZE - 1) / BLOCK_SIZE;


				//************************ Starting expander function **********************************************/

			  if(debug)
					std::cout << "Starting main expand once kernel " << " Number of edges to process is " << grid_size <<std::endl;

				expander<<< grid_size, BLOCK_SIZE >>>(
					process_vetxes, inter_vertex,
					temp_from, temp_to,
					graph.full_vertex_array, graph.full_edge_array,
					position_in_array,
					expanded_array,
					inter_vertex_array,
					number_edges_to_process,
					graph.number_of_vertex,
					current_level,
					positions_vertex_current_level);

				cudaDeviceSynchronize();
				device_free(temp_from);
				device_free(temp_to);
				device_free(process_vetxes);
				device_free(inter_vertex);
				device_free(position_in_array);
				cudaFree(0);

				if (debug)
					std::cout << "Expand kernel finished working " <<std::endl;


				/*
				*	Remove empty, non used data
				*/
				remove_if(device,
					make_zip_iterator(
						make_tuple(expanded_array,
											 inter_vertex_array)),
					make_zip_iterator(
						make_tuple(expanded_array + full_size,
					   					 inter_vertex_array + full_size)),
					is_minus_one());
				/*********************************    Expander kernel finished ************************************/

				if (debug)
				{
					int num = full_size;
					domain temp_array = new vertex[num];
					std::cout << "Expanded candidates " << std::endl;

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

				/************************************ Unifier kernel prepare **************************************/


		  	inclusive_scan(device, positions_vertex_current_level,
									           positions_vertex_current_level + graph.number_of_vertex,
									           positions_vertex_current_level);

 	  		if(debug)
  				std::cout << "Position vertex " << positions_vertex_current_level[graph.number_of_vertex - 1] << std::endl;

				int true_size = positions_vertex_current_level[graph.number_of_vertex - 1];
				device_ptr<vertex> vertex_ending_offsets = device_malloc<vertex>(graph.number_of_vertex);

				device_ptr<vertex> from_vertex_array = device_malloc<vertex>(true_size);
				device_ptr<vertex> expanded_array_copy = device_malloc<vertex>(true_size);
				device_ptr<vertex> temp_expanded = device_malloc<vertex>(true_size);
				copy(device, expanded_array, expanded_array + true_size, temp_expanded);

				device_ptr<vertex> inter_vertex_array_copy = device_malloc<vertex>(true_size);

				device_ptr<int> indices = device_malloc<int>(true_size);

				// Remove all dublicates. Calculate the opacity, remove if needed
				if (debug)
					std::cout << "Starting unification " << std::endl;

				/******************************************* Unifier kernel start *************************************/

				grid_size = ((graph.number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE);

				unifier <<<grid_size, BLOCK_SIZE >>>(	expanded_array,
																							expanded_array_copy,
																							temp_expanded,
																							from_vertex_array,
																							inter_vertex_array,
																							inter_vertex_array_copy,
																							indices,
																							positions_vertex_current_level,
																							vertex_ending_offsets,
																							graph.number_of_vertex);

		    cudaDeviceSynchronize();
				// Now we have unique expanded_array

				int proceed_size = vertex_ending_offsets[graph.number_of_vertex - 1];
				if (debug)
				{
					 int num = graph.number_of_vertex;
					 int new_num = true_size;
					 domain temp_array = new vertex[graph.number_of_vertex];
					 std::cout << "After unification " << std::endl;
					 copy(positions_vertex_current_level, positions_vertex_current_level + num, temp_array);
					 for(int i=0; i < num; i++)
						 std::cout << temp_array[i] << " ";
					 std::cout << std::endl;

					 copy(vertex_ending_offsets, vertex_ending_offsets + num, temp_array);
					 for(int i=0; i < num; i++)
						 std::cout << temp_array[i] << " ";
					 std::cout << std::endl;
					 delete temp_array;
					 temp_array = new vertex[new_num];

					 std::cout << "\n Expanded array copy: " << std::endl;
					 copy(expanded_array_copy, expanded_array_copy + new_num, temp_array);
					 for(int i=0; i < new_num; i++)
						 std::cout << temp_array[i] << " ";
					 std::cout << std::endl;

					 copy(indices, indices + new_num, temp_array);
					 for(int i=0; i < new_num; i++)
						 std::cout << temp_array[i] << " ";
					 std::cout << std::endl;
					 delete temp_array;
				}

				int true_unique_size =
				reduce(device, vertex_ending_offsets,
											 vertex_ending_offsets + graph.number_of_vertex);

				device_free(temp_expanded);
				device_free(inter_vertex_array);

				if (debug)
					std::cout << "True unique size  " << true_unique_size << std::endl;

				if(debug)
					std::cout << "unification done " << std::endl;

				/******************************************************* Unifier kernel ready ***********************************************/

				renew_graph_offset(graph, vertex_ending_offsets, current_level);

				if(debug)
						std::cout << "Starting copier " << std::endl;

				  /******************************************************* Edge copy kernel start  ********************************************/

					temp_expanded = device_malloc<vertex>(true_unique_size);
					device_ptr<vertex> temp_from_array = device_malloc<vertex>(true_unique_size);

					 temp_edge_copier<<<grid_size, BLOCK_SIZE>>>
				 			(expanded_array,
					  	positions_vertex_current_level, vertex_ending_offsets,
							graph.full_vertex_array,
					  	temp_from_array, temp_expanded,
					  	graph.number_of_vertex, current_level);

					cudaDeviceSynchronize();
					device_free(expanded_array);
					//device_free(positions_vertex_current_level);

					if (debug)
					 {
						 int new_num = true_unique_size;

						 std::cout << "After copier " << std::endl;

						 std::cout << " Temporal expand array: ";
						 for(int i=0; i < new_num; i++)
							 std::cout << temp_expanded[i] << " ";
						 std::cout << std::endl;

						 std::cout << " Temporal from array: ";
						 for(int i=0; i < new_num; i++)
							 std::cout << temp_from_array[i] << " ";
						 std::cout << std::endl;

					 }


					/****************************************************** Edge copied finished*********************************************************/

					// Calculate how they will effect on the system
					calc_L_opacity(graph, temp_from_array, temp_expanded, vertex_ending_offsets, true_unique_size);

					renew_graph_offset(graph, vertex_ending_offsets, current_level);

					// COPY TO EDGE LIST LAST LEVEL

					int temp_offset = graph.full_vertex_array[current_level*graph.number_of_vertex - 1];

					remove_copy_if(device,
										temp_expanded,
										temp_expanded + true_unique_size,
										graph.full_edge_array + temp_offset,
										negative());

					/******************* Removing candidates former *************/

					// Create removing candidates
					if (debug)
					{
						std::cout << "Forming remove Candidates" << std::endl;
					}
					int remove_array_size = transform_removing_edges(
															 inter_vertex_array_copy,
															 expanded_array_copy,
															 positions_vertex_current_level,
															 true_size,
															 temp_from_array,
															 temp_expanded,
															 graph.number_of_vertex,
															 true_unique_size);


					if (debug)
					 {
						 int new_num = true_size;
						 std::cout << "Candidates to remove are formed " << std::endl;

						 std::cout << "\n From array copy: ";

 						 for(int i=0; i < new_num; i++)
	 							 std::cout << from_vertex_array[i] << " ";
						 std::cout << std::endl;

						 std::cout << "Inter array copy: ";

 						 for(int i=0; i < new_num; i++)
	 							 std::cout << inter_vertex_array_copy[i] << " ";
						 std::cout << std::endl;


 						 std::cout << " Expanded array copy: ";

 						 for(int i=0; i < new_num; i++)
	 							 std::cout << expanded_array_copy[i] << " ";
						 std::cout << std::endl;

					 }

					device_free(temp_expanded);
					device_free(temp_from_array);
					device_free(vertex_ending_offsets);
					device_free(positions_vertex_current_level);

					/********************************* First level remove kernel  ********************************************/
					if(debug)
						std::cout << "First level remove kernel started " << std::endl;

					grid_size = (remove_array_size+BLOCK_SIZE-1)/BLOCK_SIZE;
					first_level_remove<<<grid_size, BLOCK_SIZE>>>(
																				graph.full_vertex_array,
																				graph.full_edge_array,
																				inter_vertex_array_copy,
																				expanded_array_copy,
																				remove_array_size);
					 cudaDeviceSynchronize();
					 grid_size = (graph.number_of_vertex+BLOCK_SIZE-1)/BLOCK_SIZE;
					 calc_new_offsets<<<grid_size, BLOCK_SIZE>>>(
						 												graph.full_vertex_array,
					 													graph.full_edge_array,
																		graph.real_vertex_degrees,
																		graph.number_of_vertex);


					 std::cout << "Ok in 13" << std::endl;

					 /********************************* First level  kernel remove finished  ********************************************/
				/*	 remove(device, graph.full_edge_array,
													graph.full_edge_array + graph.full_vertex_array[graph.number_of_vertex-1],
													-1);
					 inclusive_scan(device, graph.real_vertex_degrees,
																	graph.real_vertex_degrees + graph.number_of_vertex,
																	graph.full_vertex_array);

					 fill(make_zip_iterator(make_tuple(graph.opacity_matrix, graph.lessL_matrix)),
					 			make_zip_iterator(make_tuple(graph.opacity_matrix + graph.max_degree*graph.max_degree,
																						 graph.lessL_matrix + graph.max_degree*graph.max_degree)),
								make_tuple(0.0,0.0));

						for(int i=1;i<= current_level; i++)
						{
							expand_step(graph, i);
						} */




					 for(int level = 0; level < current_level; level++)
						{
							edge_remove_step(graph,
															 inter_vertex_array_copy,
															 expanded_array_copy,
															 remove_array_size, level);
						}

						/********************************* All removing finished  ********************************************/
						/******************************* True removal *******************************************/
						remove(device, graph.full_edge_array,
													 graph.full_edge_array + graph.full_vertex_array[graph.number_of_vertex*(current_level+1)-1],
													 -1);
						inclusive_scan(device, graph.real_vertex_degrees,
																	 graph.real_vertex_degrees + graph.number_of_vertex*(current_level + 1),
																	 graph.full_vertex_array);



				  device_free(expanded_array_copy);
					device_free(inter_vertex_array_copy);

				//cudaFree(0);
				 std::cout << "Ok in 14" << std::endl;


	}
}



void order_graph(Graph graph)
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
void calc_L_opacity(Graph graph,int start_level ,int end_level)
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

		/*
		*	Expanding indexes. Finding break points
		*	Example: 0 1 2 3 4 .. 20 => 0 0 1 0 0 0 1 ...
		*/
		copy(device, graph.full_vertex_array + (i-1)*graph.number_of_vertex, graph.full_vertex_array + (i)*graph.number_of_vertex , from_vertex);
		adjacent_difference(device, from_vertex, from_vertex + graph.number_of_vertex, from_vertex);

		/*
		*	Transorming into indexes:
		*	Example:	0 0 1 0 0 0 1 => 0 0 1 1 1 1 2 2 2 ..
		*/
		from_vertex[0] = from_vertex[0] - starting_point;
	 	expand(from_vertex, from_vertex + graph.number_of_vertex, make_counting_iterator<vertex>(0), from);


		device_ptr<int> from_degree = device_malloc<int>(N);

		/*
		*	Transforming from indexes into degrees:
		*	Example:  0 0 1 1 1 1 2 2 2.. => 2 2 4 4 4 4 ...
		*/

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


		transform(
			device,
			make_zip_iterator(make_tuple(from, to)),
			make_zip_iterator(make_tuple(from+N, to+N)),
			make_zip_iterator(make_tuple(from,to)),
			min_max_transform_int());


		//
		// 	Opacity  matrix forming. Now it is n^ 2 memory TODO: IN PARALLEL using cuda kernel
		// 	Assumptions !!: Not optimum for undericted (div 2).
		// 	Problem with same degree. Example: {4 = > 4} - must count only degree of one.
		//

		device_ptr<int> opacity_index = device_malloc<int>(N);
		if (debug)
			std::cout << "Opacity ok";

		int gridsize =(N + BLOCK_SIZE - 1) / BLOCK_SIZE;
		opacity_former<<<gridsize, BLOCK_SIZE>>>(from_degree,
																						 to_degree,
																						 opacity_index,
																						 graph.degree_count,
																						 graph.opacity_matrix,
																						 graph.lessL_matrix,
																						 graph.max_degree, N);


		cudaDeviceSynchronize();

		if(i==1)
		{
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
			std::cout << "Copy ok";


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
//	graph.print_opacity_matrix();
//	graph.print_csr_graph();

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
//	calc_L_opacity(graph, 2 ,graph.L_VALUE);
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
