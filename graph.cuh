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

#include "headers.h"
#include "functors.cuh"

class Graph {

public:
// CSR with levels graph format
// Distance matrix in a nuttshell
domain full_vertex_array;
domain full_edge_array;

int L_VALUE = 2;

// Current
// domain vertex_current_end; We don't need this ?
// int current_end;

	// COO graph format (coordinate list)
domain from_array_host;
domain to_array_host;
// Copy of arrays in the device
domain from_array;
domain to_array;

// Additional arrays
domain vertex_degrees;
domain degree_count;
field max_degree;

opacity* opacity_matrix;

int number_of_vertex;
int number_of_edges;

bool directed = False;

	/************************************
	*		Read graph in Edge list format (COO)
	*		input:
	*					string file_name
	**************************************/
	void read_COO_format(string file_name)
	{
			ifstream myfile;
			myfile.open(file_name);
			myfile >> number_of_vertex >> number_of_edges;
		//	number_of_edges = 0;
			printf("Graph parametrs %d %d \n", number_of_vertex, number_of_edges);
			// reserve maximum value needed
			from_array_host = new vertex[number_of_vertex];
			to_array_host = new vertex[number_of_vertex];
			// Read a pair of vertex - vertex forming an edge
			int a, b;
			number_of_edges = 0;
			while (myfile >> a >> b)
			{
				from_array[i] = a;
				to_array[i] = b;
				number_edges++;
			}
			// Reading from file
			myfile.close();
	}

	/**********************************
	*	Print full graph CSR format
	* Require:
	*					vertex_degrees != NULL
	*					full_vertex_array != NULL
	*					after_array_init
	************************************/
	void print_csr_graph()
	{
		cout << "Vertex degrees ";
		for (auto iter : vertex_degrees)
		{
			cout << "  " << iter;
		}
		cout << endl;

		cout << "Degree count ";
		for (auto iter : degree_count)
		{
			cout << "  " << iter;
		}
		cout << endl;

		cout << "Vertex offset: ";
		for (auto iter : full_vertex_array)
		{
			cout << "  " << iter;
		}
		cout << endl;

		cout << "Connected Edges ";
		for (auto iter : full_edge_array)
		{
			cout << "  " << iter;
		}
		cout << endl;
	}

	/**********************************
	*	Print opacity matrix ON_HOST
	*	Require:
	*						opacity_matrix != NULL
	*
	*********************************/
	void print_opacity_matrix()
	{
		cout<< endl  << "Opacity : " << endl;
		for (auto i = 0; i < max_degree; i++)
		{
			for (auto j = 0; j < max_degree; j++)
				cout <<" " << opacity_matrix[max_degree*i + j];
				//	printf(" %f",  opacity_matrix[max_degree*i + j]);
			cout << endl;
		}
		cout << endl;
	}


	/**************************************
	*	Print graph in (one layer, initial state)
	* COO format (edge list) ON_HOST
	*
	*************************************/
	void print_coo_graph()
	{
		cout << "From ";
		for (auto iter : from_array_host)
		{
			cout << "  " << iter;
		}
		cout << endl;
		cout << "To   ";
		for (auto iter : to_array_host)
		{
			cout << "  " << iter;
		}
		cout << endl;


	}

	/*******************************************
	* 	Reading test graph presented in the paper "L-opacity"
	********************************************/
	void init_test_graph()
	{
		// COO format
		read_COO_format("graph.txt");
	}

	/******************************************
	*	Init arrays  via after to arrays
	*****************************************/
	void init_arrays()
	{
		int num_edges=L_VALUE*number_of_edges*sizeof(vertex);
		int num_vertex=L_VALUE*number_of_vertex*sizeof(vertex);
		if (!directed)
				num_edges *= 2; // double edges
		cudaMalloc((domain*)&full_edge_array,num_edges);
		cudaMalloc((domain*)&full_vertex_array,num_vertex);
		/* Copy arrays to device */
		cudaMemcpy(from_array,from_array_host,number_of_edges*sizeof(vertex),_HTD);
		cudaMemcpy(to_array,to_array_host,number_of_edges*sizeof(vertex),_HTD);
		// Init degree vector
		cudaMalloc((domain*)&vertex_degrees, number_of_vertex * sizeof(vertex));
		cudaMalloc((domain*)&degree_count, number_of_vertex * sizeof(vertex));
		// Init opacity matrix TODO: memory if n^2
		cudaMalloc((opacity**)&opacity_matrix, number_of_vertex* number_of_vertex * sizeof(opacity));
	}

	/*********************************************
	*  Converting from COO (edge list) format to CSR (adjaceny list) format
	*  Run it after something is in COO list (from and to).
	*		Require:
	*						directed = False
	********************************************/
	void convert_to_CSR()
	{
		/*
		* First combine and sort data from and to array - this will be our new edge_list acording to their indexes
		*/
		init_arrays();
		thrust::device_vector<vertex> temp_indx(2 * number_of_edges);

		thrust::counting_iterator<vertex> index_from(0);
		thrust::counting_iterator<vertex> index_to(number_of_edges);

		//	Merging from and to arrays are keys,
		//	indexes are (0..number_edges) and (num_edges to 2*number_edges)
		thrust::merge_by_key(thrust::device,
			from_array, from_array + number_of_edges,
			to_array, to_array + number_of_edges,
			index_from, index_to,
			temp_indx.begin(),
			full_edge_array
			);

		thrust::sort_by_key(thrust::device,
			temp_indx.begin(), temp_indx.end(),
			full_edge_array;

			/*
			*	Form vertex offset list
			*/


		thrust::reduce_by_key(thrust::device,
			temp_indx.begin(), temp_indx.end(),
			thrust::make_constant_iterator(1),
			temp_indx.begin(),
			full_vertex_array);

		/*
		*	Form degree vector
		*/

		thrust::copy(thrust::device,
			full_vertex_array,
			full_vertex_array + number_of_vertex,
			vertex_degrees);

		thrust::copy(thrust::device,
			vertex_degrees, vertex_degrees + number_of_vertex,
			degree_count);
		thrust::sort(thrust::device,
			degree_count, degree_count + number_of_vertex;
		max_degree = degree_count[number_of_vertex - 1];

		thrust::reduce_by_key(thrust::device,
			degree_count, degree_count + number_of_vertex,
			thrust::make_constant_iterator(1),
			thrust::make_discard_iterator(),
			degree_count);


		thrust::inclusive_scan(thrust::device,
			 full_vertex_array,
			 full_vertex_array+number_of_vertex,
			 full_vertex_array);

		// Clean temporal arrays

		temp_indx.clear();
		temp_indx.shrink_to_fit();

		/*
		*	Transform the edge list array according to they paired edge.
		*	Form edge list combined by vertexes
		*/
		
		thrust::transform(thrust::device,
			full_edge_array, full_edge_array + 2*number_of_edges,
			full_edge_array,
			coo_to_csr_converter(from_array, to_array,
				number_of_edges));
	}


	/***
	*  Converting from COO (edge list) format to CSR (adjacency list) format
	*  Run it after someting is in COO list (from and to).
	*/
	void convert_to_CSR_no_doubles()
	{
		/*
		* First combine and sort data from and to array - this will be our new edge_list according to their indexes
		*/
		init_arrays();
		thrust::device_vector<vertex> temp_indx(2 * number_of_edges);
		thrust::fill(temp_indx.begin(), temp_indx.end(), 0);


		thrust::counting_iterator<vertex> index_from(0);
		thrust::counting_iterator<vertex> index_to(number_of_edges);

		//	Merging from and to arrays are keys,
		//	indexes are (0..number_edges) and (num_edges to 2*number_edges)
		thrust::merge_by_key(from_array.begin(), from_array.end(),
			to_array.begin(), to_array.end(),
			index_from, index_to,
			temp_indx.begin(),
			full_edge_array.begin()
			);

		thrust::sort_by_key(temp_indx.begin(), temp_indx.end(), full_edge_array.begin());


		/*
		*	Form vertex offset list
		*/


		thrust::reduce_by_key(temp_indx.begin(), temp_indx.end(),
			thrust::make_constant_iterator(1), temp_indx.begin(), full_vertex_array.begin());

		thrust::inclusive_scan(full_vertex_array.begin(), full_vertex_array.begin() + number_of_vertex, full_vertex_array.begin());

		// Clean temporal arrays

		temp_indx.clear();
		temp_indx.shrink_to_fit();

		/*
		*	Transform the edge list array according to they paired edge.
		*	Form edge list combined by vertexes
		*/
		domain a = thrust::raw_pointer_cast(from_array.data());
		domain b = thrust::raw_pointer_cast(to_array.data());

		thrust::transform(full_edge_array.begin(), full_edge_array.begin() + 2 * number_of_edges, full_edge_array.begin(),
			coo_to_csr_converter(a, b, number_of_edges));
	}
};
