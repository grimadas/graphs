/*
 	Main application
	Author : Bulat, 2015

*/
//#include "graph.cuh"
#include "apsp.cuh"
#include "graph.cuh"
#include "headers.h"




__global__  void expander(thrust::device_vector<vertex>::iterator* previous,
	thrust::device_vector<vertex>::iterator* current,
	domain current_vertex, domain temp_from, domain temp_to, domain full_vertex_array, domain full_edge_array,
	int number_of_vertex
	)
{

	int idx = blockIdx.x*blockDim.x + threadIdx.x;
	printf("Expander = current vertex  %u \n", current_vertex[idx]);
	printf("Expander = idx  %u \n", idx);
	temp_from[idx];
	temp_from[idx];
	current[current_vertex[idx]];
	printf("First in edge %u \n", *full_edge_array);
	printf("Ex[ everything fine ] \n");

	current[current_vertex[idx]] =
	thrust::copy(thrust::device,
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_from[idx])),
		thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_to[idx])),
		current[current_vertex[idx]]);

	printf("Expander ok 1 \n");
	int start = 0;
	if (idx != 0)
	{
		start = full_vertex_array[current_vertex[idx - 1]];
	}
	int end = full_vertex_array[current_vertex[idx]];
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


	// Change it to paralel version
	// expander<<<1, NUM>>(c, full_vertex_array, full_edge_array)
	// {



		// Put a value into vertex array



	/*
	for (int i = 0; i < NUM; i++)
	{
		if (c[i] != current_index)

	{
			// We finish expanding for current index
			int start = 0;
			if (current_index != 0)
			{
				start = full_vertex_array[current_index - 1];
			}
			int end = full_vertex_array[current_index];
			//thrust::device_vector<vertex> temp_vect(number_of_edges);

			//thrust::merge(previous, current, full_edge_array.begin() + start, full_edge_array.begin() + end, temp_vect.begin());
			// Remove a self vertex
			//	current = thrust::remove(previous, current, 0);

			current = thrust::remove(previous, current, current_index + 1);
			// Remove copies

			domain from = thrust::raw_pointer_cast(&full_edge_array[start]);
			domain to = thrust::raw_pointer_cast(&full_edge_array[end]);

			cout << "Tempo " << current_index << ": ";

			//		pred_if(from, to));

		//	thrust::device_vector<vertex> temp_vector(number_of_edges);
		//	temp_vec
			thrust::sort(previous, current);
			current = thrust::unique(previous, current);
			for (auto j = full_edge_array.begin() + start; j != full_edge_array.begin() + end; j++)
			{
				current = thrust::remove(previous, current, *j);
			}

		//	current = thrust::remove_if(thrust::device, previous, current,
		//		pred_if(thrust::raw_pointer_cast(full_edge_array.data()), thrust::raw_pointer_cast(full_vertex_array.data()), current_index));

			thrust::for_each(previous, current, printer());
			cout << endl;
			full_vertex_array[current_index + number_of_vertex] = thrust::distance(previous, current);
			previous = current;
			current_index = c[i];

			// Put a value into vertex array

		}
		//thrust::make_transform_iterator(thrust::make_counting_iterator<vertex>(0), fun())
		current = thrust::copy(thrust::seq,
			thrust::make_permutation_iterator(full_edge_array.begin(), thrust::make_counting_iterator<vertex>(temp_from[i])),
			thrust::make_permutation_iterator(full_edge_array.begin(), thrust::make_counting_iterator<vertex>(temp[i])),
			current);



	}

	int start = full_vertex_array[current_index - 1];
	int end = full_vertex_array[current_index];
	// Remove current vertex (avoid loops)
	current = thrust::remove(previous, current, current_index + 1);

	cout <<endl << "Tempo " << current_index << ": ";
	// sort connected edges
	thrust::sort(previous, current);
	// remove dublicates
	current = thrust::unique(previous, current);
	// Remove all vertex that are already discovered
	for (auto j = full_edge_array.begin() + start; j != full_edge_array.begin() + end; j++)
	{
		current = thrust::remove(previous, current, *j);
	}

	// Print all edges from current vertex
	thrust::for_each(previous, current, printer());
	full_vertex_array[current_index + number_of_vertex] = thrust::distance(previous, current);
	*/
}


/*
*	By finding shortest paths, form to L_VALUE level
*/
void form_full_level_graph(Graph graph)
{

	// Number of edges in a nutshell TODO: change to cycle
	vertex number_edges_to_process =
					graph.full_vertex_array[graph.number_of_vertex - 1]
					- 0;


	thrust::device_vector<vertex> temp_to(graph.number_of_edges * 2);
	thrust::device_vector<vertex> temp_from(graph.number_of_edges * 2);

	/* Form temp to an temp from vector from edge arrray */
	thrust::copy(thrust::device, graph.full_edge_array,
	graph.full_edge_array + number_edges_to_process,
	temp_to.begin());

	thrust::copy(thrust::device,
		graph.full_edge_array, graph.full_edge_array + number_edges_to_process,
		 temp_from.begin());

	// TODO Change index get previos index
	thrust::transform(thrust::device,
		temp_from.begin(), temp_from.end(),
		temp_from.begin(), previous_el(graph.number_of_vertex + 1));
	cout << "The okey 1" << endl;
	/* Store begining and ending */
	thrust::copy(
		thrust::device,
		thrust::make_permutation_iterator(graph.full_vertex_array, temp_to.begin()),
		thrust::make_permutation_iterator(graph.full_vertex_array, temp_to.end()),
		temp_to.begin());

	thrust::copy(
		thrust::device,
		thrust::make_permutation_iterator(graph.full_vertex_array,
			temp_from.begin()),
		thrust::make_permutation_iterator(graph.full_vertex_array,
			temp_from.end()), temp_from.begin());

		cout << "The okey 2" << endl;
	thrust::device_vector<vertex>  process_vetxes(number_edges_to_process+1);
	/* Find all breaking points */
	thrust::transform(
		thrust::device,
		thrust::make_counting_iterator<vertex>(0),
		thrust::make_counting_iterator<vertex>(number_edges_to_process),
		process_vetxes.begin(),
		replacer(graph.full_vertex_array, graph.number_of_vertex));
	/* Sum all previous results */
	thrust::inclusive_scan(thrust::device, process_vetxes.begin(),
									process_vetxes.end(), process_vetxes.begin());
	/**/
	cout << "The okey 3" << endl;

	int NUM = thrust::distance(temp_to.begin(), temp_to.end());

	/* 	Graph init. */


	thrust::device_vector<vertex>::iterator* previous = new thrust::device_vector<vertex>::iterator[graph.number_of_vertex];
	thrust::device_vector<vertex>::iterator* current = new thrust::device_vector<vertex>::iterator[graph.number_of_vertex];
	thrust::device_vector<vertex>* tempo_vector = new thrust::device_vector<vertex>[graph.number_of_vertex];
	cout << "The okey 4" << endl;
	for (int i=0; i< graph.number_of_vertex; i++)
	{
		tempo_vector[i].reserve(20 * graph.number_of_vertex);
		current[i] = tempo_vector[i].begin();
		previous[i] = current[i];
	}

	domain _temp_from = thrust::raw_pointer_cast(temp_from.data());
	domain _temp_to = thrust::raw_pointer_cast(temp_to.data());
	domain _full_vertex = thrust::raw_pointer_cast(graph.full_vertex_array.data());
	domain _full_edge = thrust::raw_pointer_cast(&graph.full_edge_array[0]);
	domain _vertex_data = thrust::raw_pointer_cast(process_vetxes.data());
	cout << "The okey 5" << endl;

	expander << <1, 1 >> > (previous, current, _vertex_data , _temp_from, _temp_to,
		_full_vertex, _full_edge, graph.number_of_vertex);

			cout << "The okey 6" << endl;
	/*
	// Update vertex
	thrust::inclusive_scan(full_vertex_array.begin() + (number_of_vertex - 1), full_vertex_array.begin() + 2 * (number_of_vertex), full_vertex_array.begin() + (number_of_vertex - 1));
	// Update edges
	thrust::copy(tempo.begin(), current, full_edge_array.begin() + full_vertex_array[number_of_vertex-1]);


	print_csr_graph();
	*/

}




/*********************************
*	L opacity matrix calculation
*********************************/
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

	}
}





int main()
{

	Graph graph;
	graph.init_test_graph(); // Reading graph from the file in COO format
	graph.print_coo_graph();
	graph.convert_to_CSR();
	graph.print_csr_graph();
	form_full_level_graph(graph);
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
