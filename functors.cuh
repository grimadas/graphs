#include "headers.h"

/****************************
*	 Converter functor,
*	 INPUT:
*				_a - from_array
*				_b - to_array
*				_size - number_of_edges
*       field X - element: index of element in to (b[x]) and element in from array a[x - size]
*/
struct coo_to_csr_converter
{
  __host__ __device__
  coo_to_csr_converter(thrust::device_ptr<vertex> _a, thrust::device_ptr<vertex> _b, int _size) : a(_a), b(_b), size(_size){}

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

  thrust::device_ptr<vertex> a;
  thrust::device_ptr<vertex> b;
  int size;
};


/***************************
*   Equaluty functor
*   Returns element - 1 :
*   x - 1
****************************/
struct  previous_el
{
  __host__ __device__
  previous_el(field _nums) : nums(_nums) {}

  __host__ __device__
  field operator()(field x)
  {
    if (x < 1)
    {
      return nums;
    }
    return x - 1;
  }

  field nums;
};

/***********************************************
*   If element exists in start - end of current
*   Input: domain a - starting vector
*          domain b - ending vector
*          int vert - cureent vertex
************************************************/
struct if_exists
{
  __host__ __device__
  if_exists(thrust::device_ptr<vertex> a, thrust::device_ptr<vertex> b, int vert) : start(a), end(b), current(vert)
  {

  }

  __host__ __device__
    bool operator()(vertex x)
  {
    //	printf("SEARCHING distance %i \n", thrust::distance(start, end));
      int from = 0;
      if (current != 0)
      {
        from = end[current - 1];
      }
      int to = end[current];

      return thrust::binary_search(thrust::device, start+from, start+to, x);
  }

  thrust::device_ptr<vertex> start;
  thrust::device_ptr<vertex> end;
  int current;
};

/***************************************************
*   Replace with 1 if found in cidt - cidt + size
*   Input:  cidt - vector Input
*           size - size until we search
****************************************************/
struct  replacer
{
		__host__ __device__
		replacer(thrust::device_ptr<vertex> c, int _size) : cidt(c), size(_size)
		{

		}


		__host__ __device__
			vertex operator()(vertex t)
		{
				// Device vector temporal array (candidate)

				if (thrust::binary_search(thrust::device, cidt, cidt + size, t))
					return 1;
				return 0;
		}


		thrust::device_ptr<vertex> cidt;
		int size;
};

/****************************************************
*   Printer functor
*   Input:  each vertex t (unsigned int)
****************************************************/
struct  printer
{

  __host__ __device__
    void operator()(vertex t)
  {
      printf("%u ", t);
  }

};

/****************************************************
*   Printer functor
*   Input:  each vertex t (opacity)
****************************************************/
struct  printer_opacity
{

  __host__ __device__
    void operator()(opacity t)
  {
      printf("%f ", t);
  }
};


/***********************************************************
*   Transform min and max in pair of arrays
*   Input:   For each tuple<double double> retun (min, max )
**********************************************************/
struct min_max_transform
{

  __host__ __device__
  thrust::pair<double, double> operator()(thrust::tuple<double, double> t)
  {

      double min = thrust::get<0>(t) - thrust::get<1>(t) < 0? thrust::get<0>(t) : thrust::get<1>(t);
      double max = thrust::get<0>(t) -thrust::get<1>(t) > 0 ? thrust::get<0>(t) : thrust::get<1>(t);

      return thrust::make_pair(min, max);
  }
};


/***********************************************************
*   Transform min and max in pair of arrays
*   Input:   For each tuple<double double> retun (min, max )
**********************************************************/
struct counter
{

  __host__ __device__
  vertex operator()(thrust::tuple<vertex, vertex> t)
  {

    return thrust::get<1>(t) - thrust::get<0>(t);
  }
};


/***********************************************
*
*
***********************************************/
struct unique_edge
{
  __host__ __device__
  unique_edge(thrust::device_ptr<vertex> _start_point, thrust::device_ptr<vertex> _end_point,
              thrust::device_ptr<vertex> _vertex_edge_list, thrust::device_ptr<vertex> _full_edge_list,
              vertex _current_vertex,vertex _additional_vertex_search_index) : cidt(c), size(_size)
  {

  }


  __host__ __device__
    bool operator()(vertex t)
  {
      // Device vector temporal array (candidate)
      if (t == current_vertex)
        return false;
      bool vertex_previously_found = thrust::binary_search(thrust::device, start_point, end_point, t);
      if   
    //  if (thrust::binary_search(thrust::device, start_point, end_point, t)))
    //    return false;
      return true;
    //  if (additional_vertex_search_index !=0 && thrust::binary_search(thrust::device, ) )
  }

  thrust::device_ptr<vertex> start_point;
  thrust::device_ptr<vertex> end_point;
  int additional_vertex_search_index;
  vertex current_vertex;
  thrust::device_ptr<vertex> full_edge_list;
};
