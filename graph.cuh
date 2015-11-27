/************************************************
* 												*
*  CUDA graph representation					*
*  author: Bulat 								*
*  graph.h 										*
*************************************************
*		la - look ahead parametr
*		L (L_VALUE) - opacity threshold (1, 2) - small value
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
*				| | | | | | |  Memory: O(L*|V| + L* 2*|E|) (what about duplication? ) two choices
*		  * COO format (same situation ?)
*							   Memory: O(2*|E|)

*/

#ifndef GRAPH_USED
#define GRAPH_USED



#include "headers.h"
#include "functors.cuh"
#include "kernels.cuh"



class Graph {

public:
// CSR with levels graph format
// Distance matrix in a nuttshell
device_ptr<vertex> full_vertex_array;
device_ptr<vertex> full_edge_array;

int L_VALUE;
float threshold;

// Current
// domain vertex_current_end; We don't need this ?
// int current_end;

	// COO graph format (coordinate list)
domain from_array_host;
domain to_array_host;
domain from_to_host_matrix;
// Copy of arrays in the device
device_ptr<vertex> from_array;
device_ptr<vertex> to_array;
device_ptr<int> opacity_index;
int size_from_to;

device_ptr<int> remove_count;
device_ptr<int> removing_opacity_index;
int size_check;

// Additional arrays
device_ptr<int> initial_vertex_degrees;
device_ptr<int> real_vertex_degrees;
device_ptr<int> degree_count;
field max_degree;

device_ptr<opacity> opacity_matrix;
device_ptr<opacity> lessL_matrix;

int number_of_vertex;
int number_of_edges;

bool directed;

	int* num_edges()
	{
		return &number_of_edges;
	}

	/**********************************************
	*		Read graph in Edge list format (COO)
	*		input:
	*					string file_name
	***********************************************/
	void read_COO_format(const char* file_name)
	{

			std::ifstream myfile;
			myfile.open(file_name);
			myfile >> number_of_vertex >> number_of_edges;
		//	number_of_edges = 0;
			// reserve maximum value needed
			from_to_host_matrix = new vertex[number_of_vertex * number_of_vertex];
			fill(from_to_host_matrix, from_to_host_matrix + number_of_vertex*number_of_vertex, 0);
			// Read a pair of vertex - vertex forming an edge
			int a, b;
			number_of_edges = 0;
			while (myfile >> a >> b)
			{
				from_to_host_matrix[number_of_vertex*min(a,b) + max(a,b)] = 1;
				number_of_edges++;
			}
			from_array_host = new vertex[number_of_edges];
			to_array_host = new vertex[number_of_edges];
			number_of_edges = 0;
			for(int i =0; i< number_of_vertex; i++)
			{
				for(int j = i+1; j< number_of_vertex; j++)
				{
					if (from_to_host_matrix[i*number_of_vertex + j] == 1)
					{
						from_array_host[number_of_edges] = i;
						to_array_host[number_of_edges] = j;
						number_of_edges++;
					}
				}
			}
			// Reading from file
			delete from_to_host_matrix;
			if (debug)
				printf("Graph parametrs %d %d \n", number_of_vertex, number_of_edges);

			myfile.close();
	}

	/**********************************
	*	Print full graph CSR format
	* Require:
	*					initial_vertex_degrees != NULL
	*					full_vertex_array != NULL
	*					after_array_init
	************************************/
	void print_csr_graph()
	{
		std::cout << "\n Vertex degrees :";
		domain a = new vertex[number_of_vertex];
		copy(initial_vertex_degrees, initial_vertex_degrees + number_of_vertex, a);
		for(int i=0; i < number_of_vertex; i++)
		{
			 std::cout << a[i] << " ";
		}

		std::cout << "\n Real vertex degree :";
		delete a;
		a = new int[L_VALUE*number_of_vertex];
		copy(real_vertex_degrees, real_vertex_degrees + L_VALUE*number_of_vertex, a);
		for(int i=0; i < L_VALUE*number_of_vertex; i++)
		{
			 std::cout << a[i] << " ";
		}

		std::cout << "\n Degree count :";
		domain b= new vertex[max_degree];
		copy(degree_count, degree_count + max_degree, b);
		for(int i=0; i < max_degree; i++)
		{
			 std::cout << b[i] << " ";
		}

		std::cout << "\n Vertex offset: ";
		domain c = new vertex[L_VALUE* number_of_vertex];
		copy(full_vertex_array, full_vertex_array +  L_VALUE*number_of_vertex, c);
		for(int i=0; i < L_VALUE*number_of_vertex; i++)
		{
			 std::cout << c[i] << " ";
		}

		std::cout << "\n Connected Edges ";
		int size_to_print = full_vertex_array[L_VALUE*number_of_vertex - 1];
		domain d = new vertex[size_to_print];
		copy(full_edge_array, full_edge_array + size_to_print, d);
		for(int i=0; i < size_to_print; i++)
		{
			 std::cout << d[i] << " ";
		}
		delete a,b,c,d;
	}

	/**********************************************
	*	Print opacity matrix ON_HOST
	*	Require:
	*						opacity_matrix != NULL
	*
	***********************************************/
	void print_opacity_matrix()
	{
		printf("\n Opacity : \n");
		opacity* a = new opacity[max_degree * max_degree];
		copy(opacity_matrix, opacity_matrix + max_degree*max_degree, a);
		for (int i = 0; i< max_degree; i++)
		{
			for (int j = 0; j< max_degree; j++)
			{
				std::cout << a[i*max_degree + j] << " ";
			}
			std::cout << std::endl;
		}
		printf("\n Number of edges less L : \n");
		copy(lessL_matrix, lessL_matrix + max_degree*max_degree, a);
		for (int i = 0; i< max_degree; i++)
		{
			for (int j = 0; j< max_degree; j++)
			{
				std::cout << a[i*max_degree + j] << " ";
			}
			std::cout << std::endl;
		}

		delete a;
	}


	/********************************************
	*	Print graph in (one layer, initial state)
	* COO format (edge list) ON_HOST
	*
	*********************************************/
	void print_coo_graph()
	{
		int total_size = number_of_edges;
		std::cout << "EDGES " << number_of_edges;
		std::cout << " From: \n";
		domain a = new vertex[total_size];
		copy(from_array, from_array + total_size, a);
		for(int i=0; i < total_size; i++)
		{
			 std::cout << a[i] << " ";
		}

		std::cout << "\n To: \n";
		domain b= new vertex[total_size];
		copy(to_array, to_array + total_size, b);
		for(int i=0; i < total_size; i++)
		{
			 std::cout << b[i] << " ";
		}
		std::cout << std::endl << " Opacity index: \n";
		copy(opacity_index, opacity_index + total_size, b);
		for(int i=0; i < total_size; i++)
		{
			 std::cout << b[i] << " ";
		}
		std::cout << std::endl;
		delete a,b;
	}

	/****************************************************
	* 	Reading test graph presented in the paper "L-opacity"
	*******************************************************/
	void init_test_graph()
	{
		// COO format

		read_COO_format("graph.txt");
		std::cout << "Reading finidhed" << std::endl;
		std::cout << "Vertex " << number_of_vertex << std::endl;
		std::cout << "Edges " << number_of_edges << std::endl;

	}

	/******************************************
	*	Init arrays  via after to arrays
	*****************************************/
	void init_arrays()
	{

		from_array = device_malloc<vertex>(2*number_of_edges);
		to_array = device_malloc<vertex>(2*number_of_edges);
		opacity_index = device_malloc<int>(2*number_of_edges);


			/* Copy arrays to device */
		copy(from_array_host, from_array_host + number_of_edges, from_array);
		copy(to_array_host, to_array_host + number_of_edges, to_array);
		delete from_array_host, to_array_host;

		initial_vertex_degrees = device_malloc<vertex>(number_of_vertex);
<<<<<<< HEAD
		real_vertex_degrees = device_malloc<vertex>(L_VALUE*number_of_vertex);
=======

>>>>>>> origin/Save_Thrust

		int num_vertex=L_VALUE*number_of_vertex;
		//	if (!directed)
		//			num_edges *= 2; // double edges
		full_vertex_array = device_malloc<vertex>(num_vertex);
		real_vertex_degrees = device_malloc<vertex>(num_vertex);
		fill(device, full_vertex_array, full_vertex_array + num_vertex, 0);
	}

	/********************************************************************
	*  Converting from COO (edge list) format to CSR (adjaceny list) format
	*  Run it after something is in COO list (from and to).
	*		Require:
	*						directed = False
	*						from_array contains values
	*						to_array contains values
	*		Input:  Create arrays ? : bool
	*						change degree properties : ? bool
	********************************************************************/
	void convert_to_CSR(bool create_arrays, bool change_properties)
	{
		/*
		* First combine and sort data from and to array - this will be our new edge_list acording to their indexes
		*/
		if (debug)
			std::cout << "Starting convertion ";

		if (create_arrays)
			init_arrays();
		else
			fill(device, full_vertex_array, full_vertex_array+ L_VALUE*number_of_vertex, 0);

		device_ptr<vertex> temp_indx = device_malloc<vertex>(2*number_of_edges);
		device_ptr<vertex> temp_indx2 = device_malloc<vertex>(2*number_of_edges);

		//wrap raw pointer with a device_ptr to use with Thrust functions
		fill(device, temp_indx, temp_indx + 2*number_of_edges, 0);
		counting_iterator<vertex> index_from(0);
		counting_iterator<vertex> index_to(number_of_edges);

		//	Merging from and to arrays are keys,
		//	indexes are (0..number_edges) and (num_edges to 2*number_edges)
			// Copy from to values
			copy(device, from_array, from_array + number_of_edges, temp_indx);
			copy(device, to_array, to_array + number_of_edges, temp_indx + number_of_edges);
			// Copy indexes
			copy(device, index_from, index_from + number_of_edges, temp_indx2);
			copy(device, index_to, index_to + number_of_edges, temp_indx2 + number_of_edges);
			if (debug)
				std::cout << "Merge ok : ";

			sort_by_key(device,
			temp_indx, temp_indx + 2*number_of_edges,
			temp_indx2);

			if(debug)
				std::cout << "Sort ok: ";

		/*
		*	Form vertex offset list
		*/

		int gridsize =(2*number_of_edges + BLOCK_SIZE - 1) / BLOCK_SIZE;
		degree_count_former<<<gridsize, BLOCK_SIZE>>>(temp_indx, full_vertex_array,2*number_of_edges,0);
		cudaDeviceSynchronize();

		copy(device,
			full_vertex_array,
			full_vertex_array + number_of_vertex,
			real_vertex_degrees);

		if(change_properties)
		{
			//
			//	Form degree vector.
			//	Each vertex has degree
			//	Total size: number_of_vertex
			//
			copy(device,
				full_vertex_array,
				full_vertex_array + number_of_vertex,
				initial_vertex_degrees);
			if(debug)
				std::cout << "Copy ok";

		// Find maximum degree value
			max_degree = reduce(device, initial_vertex_degrees,
																initial_vertex_degrees + number_of_vertex,
																0, maximum<vertex>());
		if(debug)
			std::cout << "Degree ok";

		// Opacity matrix should be create again
		opacity_matrix = device_malloc<opacity>(max_degree*max_degree);
		fill(device, opacity_matrix, opacity_matrix + max_degree*max_degree, 0);

		// Malloc lessL_matrix in memory: n^2
	 	lessL_matrix= device_malloc<opacity>(max_degree*max_degree);
		fill(device,  lessL_matrix, lessL_matrix + max_degree*max_degree, 0);

		degree_count = device_malloc<vertex>(max_degree);
			fill(device, degree_count, degree_count + max_degree, 0);
			gridsize = (number_of_vertex + BLOCK_SIZE - 1) / BLOCK_SIZE;
			// Offset is 1
			degree_count_former<<<gridsize, BLOCK_SIZE>>>
																(initial_vertex_degrees, degree_count,number_of_vertex, 1);
			cudaDeviceSynchronize();
		 }

		//
		//	Form vertex offset array
		//	Result: vertex offser array => 2 4 10 ...
		//
		inclusive_scan(device,
			 full_vertex_array,
			 full_vertex_array+number_of_vertex,
			 full_vertex_array);

		if(debug)
			std::cout << "Inclusive scan ok";
		// Clean temporal array
		device_free(temp_indx);


		int num_edges=number_of_vertex*max_degree*L_VALUE;
		if (num_edges > number_of_vertex*number_of_vertex)
		{
			num_edges = number_of_vertex*number_of_vertex;
		}
		device_free(full_edge_array);
		full_edge_array = device_malloc<vertex>(num_edges);
		fill(device, full_edge_array, full_edge_array + num_edges, -1);


		/*
		*	Transform the edge list array according to they paired edge.
		*	Form edge list combined by vertexes
		*/

		transform(device,
			temp_indx2, temp_indx2 + 2*number_of_edges,
			full_edge_array,
			coo_to_csr_converter(from_array, to_array,
				number_of_edges));


		device_free(temp_indx2);

		if(debug)
			std::cout<< std::endl << "Converting to CSR done " << std::endl;
	}
};
#endif
