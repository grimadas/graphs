/*

*/


/************************************************
* 													 										*
*  CUDA graph representation										*
*  author: Bulat 																*
*  graph.h 																			*
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


using namespace std;

#define vertex  unsigned int
#define edge  unsigned int

#define domain vertex*
#define field  vertex


	int L_VALUE = 1;

  // CSR graph format
  thrust::device_vector<vertex> vertex_array;
  thrust::device_vector<edge> edge_array;
  // CSR with levels graph format
  // Distance matrix in a nuttshell
  thrust::device_vector<vertex> full_vertex_array;
  thrust::device_vector<edge> full_edge_array;

  // COO graph format (coordinate list)
  thrust::device_vector<vertex> from_array;
  thrust::device_vector<vertex> to_array;
  // Distance oracle

  // Additional arrays
  thrust::device_vector<int> vertex_degrees;

  // Storing shortest path in COO format matrix
  thrust::device_vector<vertex> from_array_SP;
  thrust::device_vector<vertex> to_array_SP;

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
	*	Print graph CSR format
	*	INLINE Grpah function
	*/
  void print_csr()
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
 		 cout <<"  " << iter;
 	 }
 	 cout << endl;
 	 cout << "To   ";
 	 for (auto iter : to_array)
 	 {
 		 cout << "  "<< iter;
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
		converter(domain _a, domain _b, int _size) : a(_a), b(_b), size(_size){}

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


   void convert()
   {

	   thrust::counting_iterator<vertex> first(0);
	   thrust::counting_iterator<vertex> last = first + number_of_edges;
	   /*
	   * First combine and sort data from and to array - this will be our new edge_list acording to their indexes
	   */

	   thrust::device_vector<vertex> temp_indx(2 * L_VALUE* number_of_edges);
	   thrust::fill(temp_indx.begin(), temp_indx.end(), 0);

	   full_edge_array = temp_indx;

	   thrust::device_vector<vertex> index_from(number_of_edges);
	   thrust::device_vector<vertex> index_to(number_of_edges);

	   thrust::sequence(index_from.begin(), index_from.end(), 0);
	   thrust::sequence(index_to.begin(), index_to.end(), number_of_edges);

	   cout << "Full edge array  ";
	   for (auto iter : full_edge_array)
	   {
		   cout << "  " << iter;
	   }
	   cout << endl;

	   thrust::merge_by_key(from_array.begin(), from_array.end(),
		   to_array.begin(), to_array.end(),
		   index_from.begin(), index_to.begin(),
		   temp_indx.begin(),
		   full_edge_array.begin()
		   );

	   thrust::sort_by_key(temp_indx.begin(), temp_indx.end(), full_edge_array.begin());

	   thrust::device_vector<vertex> temp1(2 * number_of_edges);
	   thrust::fill(temp1.begin(), temp1.end(), 1);

	   thrust::device_vector<vertex> temp_arr(number_of_vertex);

	   thrust::reduce_by_key(temp_indx.begin(), temp_indx.end(), temp1.begin(), temp_indx.begin(), temp_arr.begin());
	   full_vertex_array.reserve(L_VALUE * number_of_vertex);
	   full_vertex_array = temp_arr;

	   cout << "Temp indx  ";
	   for (auto iter : temp_indx)
	   {
		   cout << "  " << iter;
	   }
	   cout << endl;

	   cout << "Temp 1   ";
	   for (auto iter : temp1)
	   {
		   cout << "  " << iter;
	   }
	   cout << endl;


	   /*
	   *	Transform the edge list array according to they
	   */
	   domain a = thrust::raw_pointer_cast(from_array.data());
	   domain b = thrust::raw_pointer_cast(to_array.data());

	   thrust::transform(full_edge_array.begin(), full_edge_array.end(), full_edge_array.begin(), converter(a, b, number_of_edges));

   }


	void test_func()
	{


	 convert();

	 cout << "Full vertex  ";
	 for (auto iter : full_vertex_array)
	 {
		 cout << "  " << iter;
	 }
	 cout << endl;


	 cout << "Full edge  ";
	 for (auto iter : full_edge_array)
	 {
		 cout << "  " << iter;
	 }
	 cout << endl;

	}

	/***
  *  Converting from COO (edge list) format to CSR (adjaceny list) format
  *  ! Run it after someting is in COO list (from and to)
  */
  void convert_to_CSR()
  {

	  full_vertex_array.reserve(2*L_VALUE * number_of_edges);
	  full_edge_array.reserve(2 * L_VALUE * number_of_vertex);


  }
#endif
