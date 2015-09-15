/*
*  CUDA graph representation
*  author: Bulat
*  graph.h
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

using namespace std;

#define vertex  unsigned int
#define edge  unsigned int

#define domain vertex*
#define field  vertex

/*
* Class Graph
* representation
*/
class Graph
{

private:

	int L_VALUE = 1;

	// CSR graph format
  thrust::device_vector<vertex> vertex_array;
  thrust::device_vector<edge> edge_array;
	// CSR with levels graph format
	thrust::device_vector<vertex> vertex_array_exp;
  thrust::device_vector<edge> edge_array_exp;

  // COO graph format (coordinate list)
  thrust::device_vector<vertex> from_array;
  thrust::device_vector<vertex> to_array;
	// Distance oracle

	// Vertex adjacency

  // Additional arrays
  thrust::device_vector<int> vertex_degrees;

  // Storing shortest path in COO format matrix
  thrust::device_vector<vertex> from_array_SP;
  thrust::device_vector<vertex> to_array_SP;

  unsigned int number_of_vertex;
  unsigned int number_of_edges;
  // __device__ int d1;
  bool directed = false;


public:
  /*
  * Std constructor
  */
  Graph()
  {

  }



  void random_graph()
  {

  }

  double get_density()
  {
      // For directed
      return directed? (double)number_of_edges/ (number_of_vertex * (number_of_vertex - 1)) :  (double)2 * number_of_edges/ (number_of_vertex * (number_of_vertex - 1));
  }

  void test_func()
  {
	
		  // thrust::merge(from_array.begin(), from_array.end(), to_array.begin(), to_array.end(), )
		  from_array.insert(from_array.end(), to_array.begin(), to_array.end());
		  
		  int i = thrust::count(from_array.begin(), from_array.end(), 2);
		  printf("2 degree %d ", i);
  }	


  void calc_degree()
  {
//	  thrust::make_counting_iterator()
	  // #TODO try to make directly
	//  thrust::device_vector<vertex> merged = thrust::merge(from_array.begin(), from_array.end(), to_array.begin(), to_array.end(), merged);
	  vertex_degrees.reserve(number_of_vertex);

	 int d =
	  thrust::count(from_array.begin(), from_array.end(), 0);
//	 thrust::for_each(thrust::constant_iterator<vertex>(0), thrust::constant_iterator<vertex>(number_of_vertex),
//	thrust::count(from_array.begin(), from_array.end(), count_functor(from_array, vertex_degrees));
//	 thrust::transform(from_array.begin(), from_array.end(), vertex_degrees.begin(), thrust::constant_iterator<vertex>(2));
	 cout << "The degree " << d;
	  //thrust::transform();
  }

  /**
  * Generate a random graph with given parameters
  * Input:
  *      - int vertex_num
  *      - int edge_num
  */

  void random_sparse_derict_graph(int vertex_num)
  {
     int max_num = 5;
     int offset = 0;

     directed = true;

    //number_of_edges = edge_num;
    number_of_vertex = vertex_num;
    srand(time(NULL));
    // debug
//    cout << vertex_num << " " << edge_num << endl;
    vertex_array.reserve(number_of_vertex);
    for (int i =0; i < number_of_vertex; i++ )
    {
          int rand_num = rand()% max_num;
          vertex_array.push_back(offset + rand_num);

          //cout << vertex_array.at(i);

          offset += rand_num;

    }
    number_of_edges = vertex_array[number_of_vertex - 1];



    edge_array.reserve(number_of_edges);




    for (int i=0; i < number_of_edges; i++)
    {
      edge_array.push_back(rand()% vertex_num);


    }
    //cout << "Vector check " << edge_array[0] << " \n";



  }

  int get_first_free(int* array, int size)
  {
      for (int i = 0; i < size; i++)
     {    if (array[i] < 0)
         return i;
     }
     return size;

  }
  /**
  * Generating random graph in COO format
  *
  */
  void random_coo_graph(int vertex_num, int max_num_per_vertex)
  {

	  int offset = 1;

	  directed = false;

	  //number_of_edges = edge_num;
	  number_of_vertex = vertex_num;
	  from_array.reserve(max_num_per_vertex*number_of_vertex);
	  to_array.reserve(max_num_per_vertex*number_of_vertex);


	  // Random seed
	  srand(time(NULL));

	  for (int current_vertex = 0; current_vertex < vertex_num-2; current_vertex++)
	  {
		  // Generating data for all the vertex

		  int number_of_edges = rand() % max_num_per_vertex;
		  for (int j = 0; j < number_of_edges; j++)
		  {
			  int rand_vertex = current_vertex+1 + rand() % (vertex_num-current_vertex);

			  if (rand_vertex < vertex_num)
			  //  if (thrust::find(to_array.begin(), to_array.end(), rand_vertex) != to_array.end())
			  {
				  from_array.push_back(current_vertex);
				  to_array.push_back(rand_vertex);
			  }

		  }
	  }


  }

  /**
  * Generating random graph in COO format
  * EXPEREMETNAL
  */
  void random_coo_graph_expr(int vertex_num, int max_num_per_vertex)
  {

	  int offset = 1;

	  directed = false;

	  //number_of_edges = edge_num;
	  number_of_vertex = vertex_num;
//	  from_to_array.reserve(2*max_num_per_vertex*number_of_vertex);


	  // Random seed
	  srand(time(NULL));

	  for (int current_vertex = 0; current_vertex < vertex_num - 2; current_vertex++)
	  {
		  // Generating data for all the vertex

		  int number_of_edges = rand() % max_num_per_vertex;
		  for (int j = 0; j < number_of_edges; j++)
		  {
			  int rand_vertex = current_vertex + 1 + rand() % (vertex_num - current_vertex);

			  if (rand_vertex < vertex_num)
				  //  if (thrust::find(to_array.begin(), to_array.end(), rand_vertex) != to_array.end())
			  {
			//	  cout << current_vertex << " ok" << endl;
	//			  from_to_array[current_vertex + j] = current_vertex;
		//		  from_to_array[2 * current_vertex + j] = rand_vertex;
			//	  cout << current_vertex << " still ok" << endl;
			  }

		  }



	  }


  }
  /*
  *  Print graph coo format
  *  Print is on  HOST
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

 
  void init_test_graph()
  {
      // COO format
	  read_COO_format("graph.txt");

      // CSR format
//	  read_CSR_graph("graph.txt");

  }


  void read_COO_format(string file_name)
  {
	  ifstream myfile;
	  int vertex_num, edge_num;
	  myfile.open(file_name);
	  myfile >> vertex_num >> edge_num;
	  printf("%d %d \n", vertex_num, edge_num);
	  // reserve maximum value needed
	  from_array.reserve(edge_num);
	  to_array.reserve(edge_num);
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


  /*
  *	Reading function
  * Input: FileName
  */
  void read_CSR_graph(string file_name)
  {
	  ifstream myfile;
	  int vertex_num, edge_num;
	  myfile.open(file_name);
	  myfile >> vertex_num >> edge_num;
	  // reserve maximum value needed
	  vertex_array_exp.reserve(vertex_num*(L_VALUE + 1));
	  edge_array_exp.reserve(vertex_num*vertex_num); // TODO more efficient space wasting assumption
	  // Read a pait of vertex - vertex forming an edge
	  int a, b;
	  while (myfile >> a >> b)
	  {
		  // Increment offset
		  vertex_array_exp.push_back(a);
		  // add new connection to edge array
		  edge_array_exp.push_back(b);
	  }
	  // Reading from file
	  myfile.close();

  }

  /**
  *
  */
  void random_sparse_underict_graph(int vertex_num)
  {
     int max_num = 5;
     int offset = 0;

     directed = false;

    //number_of_edges = edge_num;
    number_of_vertex = vertex_num;
    // Generae temporal array to store edges
     int** temp_array = new int*[number_of_vertex];
      for (int i = 0; i < number_of_vertex; ++i) {
          temp_array[i] = new int[max_num];
          for (int j = 0; j< max_num; j++)
          {
              temp_array[i][j] = -1;
          }
      }
   // Random seed
    srand(time(NULL));


    vertex_array.reserve(number_of_vertex);
    // Main loop for generating
    for (int vertex_indx =0; vertex_indx < number_of_vertex; vertex_indx++ )
    {
          // Generate random number of edges
          int rand_num = rand()% max_num;
          // add the offset and number of edges to vector
          vertex_array.push_back(offset + rand_num);


          offset += rand_num;
          // strating index is not always zero
          int edges_generated = get_first_free(temp_array[vertex_indx], max_num);
          if (edges_generated > rand_num && edges_generated < max_num)
          {
              vertex_array[vertex_indx] += edges_generated-rand_num;
              offset += edges_generated-rand_num;
          }
          while (edges_generated < rand_num)

          {
              int rand_vertex = rand() % vertex_num;
              if (rand_vertex != vertex_indx) // no cycles allowed
               {
                   int index_to_put = get_first_free(temp_array[rand_vertex], max_num);
                   cout << "Index to put " << index_to_put << endl;
                   if (index_to_put < max_num) // there should be place to put
                   {
                       temp_array[vertex_indx][edges_generated] = rand_vertex;
                       temp_array[rand_vertex][index_to_put] = vertex_indx;
                       edges_generated++;
                   }

               }

          }


    }
    number_of_edges = vertex_array[number_of_vertex - 1];


      for (int i = 0; i < number_of_vertex; ++i) {

          for (int j = 0; j < max_num; j++) {
              cout << " " <<temp_array[i][j];
          }
          cout << endl;
      }

    edge_array.reserve(number_of_edges);


    for (int vertexs = 0; vertexs < number_of_vertex; vertexs++)
    {
        int i = 0;
        while (temp_array[vertexs][i] >= 0)
        {
            edge_array.push_back(temp_array[vertexs][i]);
            i++;
        }
    }

    //cout << "Vector check " << edge_array[0] << " \n";



  }

  /*
  *  Print graph Vertex and Edges
  */
  void print_graph()
  {
    cout << "Vertex ";
    for(auto iter: vertex_array)
    {
      cout << iter << " ";
    }
    cout << endl;
    cout << "Edges ";
    for(auto iter: edge_array)
    {
      cout << iter << " ";
    }
    cout << endl;

  }

  /**
  *     Print graph in list maneer
  */
  void print_graph_in_list()
  {
        cout << "Graph:  " << endl;
        int starting_index = 0;



        for (int i = 0; i < number_of_vertex; i++)
        {

            cout << "  " << i << " : ";
            ///for (int edge = starting_index; edge < (int)vertex_array[i]; edge++)
            int count = (int)vertex_array[i];
            for (int j = starting_index; j <= count; j++) {

                cout << "-> " << edge_array[j];
            }
            starting_index = vertex_array[i]+1;
            cout << endl;

        }
  }

};
#endif
