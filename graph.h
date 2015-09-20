/*

*/


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


#ifndef graph_h
#define graph_h

// STL includes
#include <stdio.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
// Thrust includes
#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/find.h>
#include <thrust/count.h>
#include <thrust/reduce.h>
#include <thrust/merge.h>
#include <thrust/sequence.h>
#include <thrust\sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
#include <thrust\iterator\counting_iterator.h>
#include <thrust\iterator\permutation_iterator.h>


using namespace std;

int a = 0;

class Graph
{

#define vertex  unsigned int
#define edge  unsigned int

#define domain vertex*
#define field  vertex


	int L_VALUE = 2;
public:

	// CSR with levels graph format
	// Distance matrix in a nuttshell
	thrust::device_vector<vertex> full_vertex_array;
	thrust::device_vector<edge> full_edge_array;

	// Current
	thrust::device_vector<vertex>::iterator vertex_current_end;
	int current_end;

	// COO graph format (coordinate list)
	thrust::device_vector<vertex> from_array;
	thrust::device_vector<vertex> to_array;
	// Distance oracle
	// ?
	// Additional arrays
	thrust::device_vector<int> vertex_degrees;


	unsigned int number_of_vertex;
	unsigned int number_of_edges;

	/** ** ** **
	*		Read graph in Edge list format (COO)
	*		input:
	*					string file_name
	*/
	void read_COO_format(string file_name)
	{
		ifstream myfile;
		myfile.open(file_name);
		myfile >> number_of_vertex >> number_of_edges;
		printf("%d %d \n", number_of_vertex, number_of_edges);
		// reserve maximum value needed
		from_array.reserve(number_of_edges);
		to_array.reserve(number_of_vertex);
		// Read a pair of vertex - vertex forming an edge
		int a, b;
		while (myfile >> a >> b)
		{
			from_array.push_back(a);
			to_array.push_back(b);
		}
		// Reading from file
		myfile.close();
	}

	/**
	*	Print full graph CSR format
	*
	*/
	void print_csr_graph()
	{
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

	/**
	*	Print graph in (one layer, initial state) COO format (edge list)
	*
	*/
	void print_coo_graph()
	{
		cout << "From ";
		for (auto iter : from_array)
		{
			cout << "  " << iter;
		}
		cout << endl;
		cout << "To   ";
		for (auto iter : to_array)
		{
			cout << "  " << iter;
		}
		cout << endl;


	}

	/**
	* 	Reading test graph presented in the paper "L-opacity"
	*/
	void init_test_graph()
	{
		// COO format
		read_COO_format("graph.txt");

	}
	// ----------------------------------------------------------------
	/**
	*	 Converter functor,
	*	 INPUT:
	*				_a - from_array
	*				_b - to_array
	*				_size - number_of_edges
	*/
	struct coo_to_csr_converter
	{
		__host__ __device__
		coo_to_csr_converter(domain _a, domain _b, int _size) : a(_a), b(_b), size(_size){}

		__host__ __device__
			field operator()(field x)
		{

				if (x < size)
				{
					return b[x];
				}
				else
				{
					return a[x - size];
				}


			}

		domain a;
		domain b;
		int size;
	};


	void init_arrays()
	{

	//	full_edge_array.reserve(2 * L_VALUE * number_of_edges);
	//	full_vertex_array.reserve(L_VALUE*number_of_vertex);

		thrust::device_vector<vertex> temp_indx(2 * L_VALUE* number_of_edges);
		thrust::fill(temp_indx.begin(), temp_indx.end(), 0);
		full_edge_array = temp_indx;
		temp_indx.resize(L_VALUE*number_of_vertex);
		temp_indx.shrink_to_fit();

		full_vertex_array = temp_indx;
		temp_indx.clear();
		temp_indx.shrink_to_fit();

		current_end = 2 * number_of_edges;

	}

	/***
	*  Converting from COO (edge list) format to CSR (adjaceny list) format
	*  Run it after someting is in COO list (from and to).
	*/
	void convert_to_CSR()
	{
		/*
		* First combine and sort data from and to array - this will be our new edge_list acording to their indexes
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

		thrust::inclusive_scan(full_vertex_array.begin(), full_vertex_array.begin()+number_of_vertex, full_vertex_array.begin());

		// Clean temporal arrays

		temp_indx.clear();
		temp_indx.shrink_to_fit();

		/*
		*	Transform the edge list array according to they paired edge.
		*	Form edge list combined by vertexes
		*/
		domain a = thrust::raw_pointer_cast(from_array.data());
		domain b = thrust::raw_pointer_cast(to_array.data());

		thrust::transform(full_edge_array.begin(), full_edge_array.begin() + 2*number_of_edges, full_edge_array.begin(),
			coo_to_csr_converter(a, b, number_of_edges));



	}

	struct test_func
	{
		__host__ __device__
		test_func(Graph* graph) : 	gr(graph) {}

		__host__ __device__
			void operator()(field x)
		{
	//			gr->
	//			single_bfs(x);
		}

		Graph* gr;
	};


	struct  copier
	{
	  __host__ __device__
	 copier(domain _a, int _n) : a(_a), n(_n) {}

		__host__ __device__
		void test()
		{
			thrust::unique(thrust::device,a, a+n);
		}

		__host__ __device__
			void operator()(field x)
		{
				 test();
		}
		domain a;
		int n;
	};


	struct  equal
	{

		__host__ __device__
		field operator()(field x)
		{
			
			return x - 1;
		}

	};

	struct  prev
	{
		__host__ __device__
		prev(int _nums) : nums(_nums) {}

		__host__ __device__
		field operator()(field x)
		{
			if (x < 2)
			{
				return nums;
			}
			return x - 2;
		}
		int nums;

	};

	struct pred_if
	{
		
		__host__ __device__
			bool operator()(vertex t)
		{
				return true;
		}
	};


	struct  calcus
	{
		__host__ __device__
		calcus(domain& c) : edge_array(c) 
		{

		}

		template <typename Tuple>
		__host__ __device__
			void operator()(Tuple t)
		{
				// Device vector temporal array (candidate)
				thrust::device_vector<vertex> temp_vector(thrsut::get<1>(t) - thrust::get<0>(t));
				// copy indexes to temp 
				thrust::copy(thrust::device, thrust::counting_iterator(thrust::get<0>(t)),
					thrust::counting_iterator(thrust::counting_iterator(thrust::get<1>(t))), temp_vector.begin());
				
				current = thrust::copy_if(
					thrust::make_permutation_iterator(full_vertex_array.begin(), temp_vector.begin()),
					thrust::make_permutation_iterator(full_vertex_array.end(), temp_vector.end()), 
					current, 
					pred_if());
					
		}
		domain edge_array;
		thrust::device_ptr<vertex> current;
		

	};


	/*
	*	Perform a single best first search
	*	Input: vertex num_vertex
	*	Require: num_vertex > 0 && num_vertex <= number_vertex
	*/
	void single_bfs(vertex num_vertex)
	{
	//	full_edge_array.resize(2 * L_VALUE * number_of_edges);
		int begin = 0;
		int end;
		if (num_vertex > 1)
		{
			begin = full_vertex_array[num_vertex - 2];
		}
		end = full_vertex_array[num_vertex - 1];


		auto past_pos = thrust::find(full_edge_array.begin() + current_end, full_edge_array.end(), 0);

		thrust::copy(full_edge_array.begin() + begin, full_edge_array.begin() + end, past_pos);

		thrust::copy
			(thrust::make_permutation_iterator(full_edge_array.begin(), full_vertex_array.begin()),
			thrust::make_permutation_iterator(full_edge_array.end(), full_vertex_array.end()),
			full_edge_array.begin());

		 // get all indexes of edge list
	//	thrust::make_permutation_iterator(full_vertex_array.begin(), full_edge_array.begin())
		// transform edge list with minus one
		//thrust::make_transform_iterator(full_edge_array.begin(), minus_one());
		//thrust::make_permutation_iterator(full_vertex_array.begin(), thrust::make_transform_iterator(full_edge_array.begin(), minus_one()));


	}

	void form_local(vertex cur_vertex, domain a)
	{


//		thrust::device_vector<vertex> temp_arr(number_of_vertex);
	//	thrust::fill(temp_arr.begin(), temp_arr.end(), 0);

		// For each vertex adjacent with cur_vertex - form candiate
		//thrust::for_each(thrust::device, a + begin, a + end, form_candidate);

	//	auto ending = thrust::copy(a, a + 5, temp_arr.begin());


	}


	/*
	*	By finding shortest paths, form to L_VALUE level
	*/

	void form_full_level_graph()
	{
	//	domain a = thrust::raw_pointer_cast(full_edge_array.data());

	//	thrust::for_each(thrust::device, full_vertex_array.begin(), full_vertex_array.begin() + number_of_vertex, copier(a, 5));
		//thrust::make_transform_iterator(full_edge_array.begin(), minus_one())
		thrust::device_vector<vertex> temp(number_of_edges*2);
	

		thrust::device_vector<vertex> temp_fin(number_of_edges * number_of_edges);
		thrust::device_vector<vertex> temp_fin2(number_of_edges * number_of_edges);
		thrust::device_vector<vertex> temp_from(number_of_edges * 2 );

		thrust::copy(full_edge_array.begin(), full_edge_array.begin() + 2 * number_of_edges, temp.begin());
		thrust::transform(temp.begin(), temp.end(), temp.begin(), equal());
	
		thrust::copy(full_edge_array.begin(), full_edge_array.begin() + 2 * number_of_edges, temp_from.begin());
		thrust::transform(temp_from.begin(), temp_from.end(), temp_from.begin(), prev(number_of_vertex + 1));
	
		thrust::copy(
		
			thrust::make_permutation_iterator(full_edge_array.begin(), temp.begin()),
			thrust::make_permutation_iterator(full_edge_array.end(), temp.end()),  temp.begin());

		thrust::copy(

			thrust::make_permutation_iterator(full_vertex_array.begin(), temp_from.begin()),
			thrust::make_permutation_iterator(full_vertex_array.end(), temp_from.end()), temp_from.begin());
	      
		
		cout << "From ";
		for (auto iter : temp_from)
		{
			cout << "  " << iter;
		}
		cout << endl;

		cout << "To   ";
		for (auto iter : temp)
		{
			cout << "  " << iter;
		}
		cout << endl;

		//thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(temp_from.begin(), temp.begin())),
		//				 thrust::make_zip_iterator(thrust::make_tuple(temp_from.end(), temp.end())), 
//)

	}


};
#endif
