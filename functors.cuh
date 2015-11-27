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
  coo_to_csr_converter(device_ptr<vertex> _a, device_ptr<vertex> _b, int _size) : a(_a), b(_b), size(_size){}

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

  device_ptr<vertex> a;
  device_ptr<vertex> b;
  int size;
};


/***************************
*   Equality functor
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
*          int vert - current vertex
************************************************/
struct if_exists
{
  __host__ __device__
  if_exists(device_ptr<vertex> a, device_ptr<vertex> b, int vert) : start(a), end(b), current(vert)
  {

  }

  __host__ __device__
    bool operator()(vertex x)
  {
    //	printf("SEARCHING distance %i \n", distance(start, end));
      int from = 0;
      if (current != 0)
      {
        from = end[current - 1];
      }
      int to = end[current];

      return binary_search(device, start+from, start+to, x);
  }

  device_ptr<vertex> start;
  device_ptr<vertex> end;
  int current;
};

/***************************************************
*   Replace with 1 if found in cidt - cidt + size;
 *   Find breaking points
*   Input:  cidt - vector Input
*           size - size until we search
****************************************************/
struct  replacer
{
		__host__ __device__
		replacer(device_ptr<vertex> c, int _size) : cidt(c), size(_size)
		{

		}


		__host__ __device__
			vertex operator()(vertex t)
		{
				// Device vector temporal array (candidate)

				if (binary_search(device, cidt, cidt + size, t))
					return 1;
				return 0;
		}


		device_ptr<vertex> cidt;
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
      printf("%d ", t);
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
*   Input:   For each tuple<double double> return (min, max )
**********************************************************/

struct min_max_transform
{

  __host__ __device__
    tuple<double, double> operator()(tuple<double, double> t)
  {

    double min = get<0>(t) - get<1>(t) < 0? get<0>(t) : get<1>(t);
    double max = get<0>(t) -get<1>(t) > 0 ? get<0>(t) : get<1>(t);
    return make_tuple(min,max);
  }
};

/***********************************
* Min max transoform for integers
*
***********************************/
struct min_max_transform_int
{
__host__ __device__
  tuple<int, int> operator()(tuple<int, int> t)
{
    int min = get<0>(t) - get<1>(t) < 0? get<0>(t) : get<1>(t);
    int max = get<0>(t) -get<1>(t) > 0 ? get<0>(t) : get<1>(t);
    return make_tuple(min,max);
}

};
/*********************************
* Removing function
**********************************/
struct removing_funct
{
  __device__
  removing_funct(
    device_ptr<int> a_full_vertex_array,
    device_ptr<vertex> a_full_edge_array,
    int a_number_of_vertex) :
  full_vertex_array(a_full_vertex_array),
  full_edge_array(a_full_edge_array),
  number_of_vertex(a_number_of_vertex)
  {

  }


  __device__
  void operator()(tuple<bool, vertex, vertex> t)
  // should_remove, where_remove, what_remove
  {
    if(get<0>(t))
    {
      int level = 1;
      int vertex_number = get<1>(t); // vertex_number
      int starting_offset =
        full_vertex_array[number_of_vertex*level + vertex_number - 1];
      int ending_offset =
        full_vertex_array[number_of_vertex*level + vertex_number];
      device_ptr<vertex> remove_value =
          lower_bound(device,
            full_edge_array + starting_offset,
            full_edge_array + ending_offset, get<2>(t));

      int dist = distance(full_edge_array + starting_offset, remove_value);
      __syncthreads();
      full_edge_array[starting_offset + dist] = -1;
    }
  }
  device_ptr<int> full_vertex_array;
  device_ptr<vertex> full_edge_array;
  int number_of_vertex;
};


/*********************************
*
*
*********************************/
struct to_minus_one
{
  __host__ __device__
  vertex operator()(vertex t)
  {
    return -1;
  }
};

struct is_minus_one
{
  template <typename Tuple>
  __host__ __device__
  bool operator()(Tuple t)
  {
    return get<0>(t)==-1? true:false;
  }
};


struct found_in
{
  __host__ __device__
  found_in(
    device_ptr<int> a_expanded_array,
    int a_new_size) :
  expanded_array(a_expanded_array),
  new_size(a_new_size)
  {

  }

  __host__ __device__
  bool operator()(vertex t)
  {

    return binary_search(device, expanded_array,
                        expanded_array + new_size, t);
  }
  device_ptr<vertex> expanded_array;
  int new_size;
};

/*****************************************************
*
*
*****************************************************/
struct to_from_vertex
{
  __host__ __device__
  to_from_vertex(device_ptr<int> a_search_array, int a_size) :
  search_array(a_search_array), size(a_size)
  {

  }


  __host__ __device__
  vertex operator()(int index)
  {

      device_ptr<int> temp_ptr = lower_bound(device, search_array,
                                                     search_array + size,
                                                     index);
      int index_vertex = distance(search_array, temp_ptr);
      return index_vertex;
  }
  device_ptr<int> search_array;
  int size;
  };

/***********************************************************
*   Form candiadets where to and what to remove
*   Input:
*       opacity_matrix, lessL_matrix, index in opacity_matrix
*   Output:  For each tuple<opacity, opacity, int>
*                return (where_start, where_end,
*                        whar_start, what_end)
*
**********************************************************/
struct form_remove_candidates
{
  __host__ __device__
  form_remove_candidates(opacity a_threshold) :
  threshold(a_threshold)
  {

  }


  __host__ __device__
  tuple<int, int> operator()(
    tuple<opacity, opacity, int> t)
  {

      opacity opacity_value = get<0>(t);
      opacity count_less_L = get<1>(t);
      // From degree value
      int index = get<2>(t);
      // To degree value

      int to_delete = ceil(count_less_L - threshold * count_less_L /opacity_value);


     return make_tuple(index, to_delete);
  }
  opacity threshold;
};
/*******
* If opacity value is more than threshold
*/
struct more_than_threshold
{
  __host__ __device__
  more_than_threshold(float a_threshold) : threshold(a_threshold)
  {

  }

  __host__ __device__
  bool operator()(tuple<opacity, opacity, int> t)
  {
    //	printf("SEARCHING distance %i \n", distance(start, end));


      return get<0>(t) > threshold;
  }

  float threshold;
};
/*******************************
* Initial remove in coo format
*
**********************************/
struct remove_index_former
{
  __host__ __device__
  remove_index_former(device_ptr<int> a_degree_count,
                      device_ptr<opacity> a_opacity_matrix,
                      device_ptr<opacity> a_lessL_matrix,
                      int a_max_degree,
                      device_ptr<vertex> a_from,
                      device_ptr<vertex> a_to,
                      device_ptr<vertex> a_index) :
  from(a_from),
  to(a_to),
  opacity_index(a_index),
  degree_count(a_degree_count),
  opacity_matrix(a_opacity_matrix),
  lessL_matrix(a_lessL_matrix),
  max_degree(a_max_degree)
  {

  }


  __host__ __device__
  void operator()(int removal)
  {
    int from_degree = (opacity_index[removal] / max_degree) + 1;
    int to_degree =  (opacity_index[removal] % max_degree) + 1;

    opacity min = degree_count[from_degree - 1] * degree_count[to_degree - 1];

     if (from_degree  == to_degree)
     {
       opacity k = degree_count[from_degree-1];
       min = k * (k-1);

      }

    //min = min * 2.0;
    opacity added_value = 1.0/ (min);
    opacity_matrix[opacity_index[removal]] = opacity_matrix[opacity_index[removal]] - added_value;
    lessL_matrix[opacity_index[removal]] = lessL_matrix[opacity_index[removal]] - 1;

    from[removal] = -1;
    to[removal] = -1;
    opacity_index[removal] = -1;
  }
  device_ptr<vertex> from;
  device_ptr<vertex> to;
  device_ptr<int> opacity_index;
  device_ptr<int> degree_count;
  device_ptr<opacity> opacity_matrix;
  device_ptr<opacity> lessL_matrix;
  int max_degree;
};

/******************************
* Calucate difference in tuple
********************************/
struct size_counter
{

  __host__ __device__
  tuple<vertex, vertex> operator()(
    tuple<vertex, vertex, vertex, vertex> t)
  {
      int from_start = get<0>(t);
      int from_end = get<1>(t);

      int to_start = get<2>(t);
      int to_end = get<3>(t);

      return make_tuple(from_end - from_start, to_end - to_start );

  }
};

struct to_negative
{
  __host__ __device__
  to_negative(int an_except_state):
        except_state(an_except_state)
  {

  }
  __host__ __device__
  vertex operator()(vertex t)
  {
      if (t==0)
        return -except_state;
      return -(t);
  }
  int except_state;
};


struct to_positive
{
  __host__ __device__
  to_positive(int an_except_state):
        except_state(an_except_state)
  {

  }
  __host__ __device__
  vertex operator()(vertex t)
  {
      if (t==-except_state)
        return 0;
      return -(t);
  }
  int except_state;
};

/****************************************
* Functor for 4tuples. positive all values
*****************************************/
struct positive
{
  __host__ __device__
  bool operator()(
    tuple<vertex, vertex, vertex, vertex> t)
  {
      return get<0>(t) >= 0 && get<1>(t) >= 0 && get<2>(t) >= 0 && get<3>(t) >= 0;

  }

};
/************************************************
* Functor for general Tuples.
Checking only if the first is non_negative
*****************************************/
struct non_negative
{

  template <typename Tuple>
  __host__ __device__
  bool operator()(Tuple t)
  {
      return get<0>(t) >= 0;

  }
};

struct negative
{

  __host__ __device__
  bool operator()(vertex t)
  {
      return t < 0;

  }
};

/**************************************************************
*    Search in full_edge_array and delete in next level
*    returns current vertex - remove_value
****************************************************************/
struct seek_and_destroy
{
   __device__
  seek_and_destroy(
      int a_what_to_search,
      int a_what_to_remove,
      int a_current_level,
      device_ptr<int> a_full_vertex_array,
      device_ptr<vertex> a_full_edge_array,
      int a_number_of_vertex,
      device_ptr<int> a_new_vertex_array ):
       what_to_search(a_what_to_search), what_to_remove(a_what_to_remove), current_level(a_current_level),
       full_edge_array(a_full_edge_array), full_vertex_array(a_full_vertex_array), new_vertex_array(a_new_vertex_array),
        number_of_vertex(a_number_of_vertex) { }

  __device__
  void operator()(int vertex_number)
  {

      // Device vector temporal array (candidate)

      int start_offset = 0;
      if(current_level!=0 || vertex_number != 0 )
      {
        start_offset = full_vertex_array[current_level*number_of_vertex + vertex_number-1];
      }
      int ending_offset = full_vertex_array[current_level*number_of_vertex + vertex_number];

      // Search the presented from_vertex (what should be removed)
      if (binary_search(device,
                        full_edge_array + start_offset,
                        full_edge_array + ending_offset, what_to_search)){

                            // Search from where it could be expanded
                            start_offset = full_vertex_array[(current_level+1)*number_of_vertex + vertex_number - 1];
                            ending_offset = full_vertex_array[(current_level+1)*number_of_vertex + vertex_number];
                            if (binary_search(device, full_edge_array + start_offset,
                                              full_edge_array + ending_offset, what_to_remove))
                                              {
                                                device_ptr<int> temp_ptr = lower_bound(device,
                                                                      full_edge_array + start_offset,
                                                                      full_edge_array + ending_offset, what_to_remove);

                                                // index where it should be removed
                                                int index_to_remove = distance(full_edge_array + start_offset, temp_ptr);
                                                __syncthreads();

                                                full_edge_array[start_offset + index_to_remove] = -1;

                                                int* k = raw_pointer_cast(new_vertex_array + vertex_number);
                                                atomicSub(k, 1);

                                              }
                          }

        }
  int what_to_search, what_to_remove, current_level;
  device_ptr<int> new_vertex_array;
  device_ptr<int> full_vertex_array;
  device_ptr<vertex> full_edge_array;
  int number_of_vertex;
};


/***********************************************************
*   Transform min and max in pair of arrays
*   Input:   For each tuple<double double> return (min, max )
**********************************************************/
struct counter
{

  __host__ __device__
  vertex operator()(tuple<vertex, vertex> t)
  {

    return get<1>(t) - get<0>(t);
  }
};



/***********************************************
* Check if the vertex was previously discovered by current_vertex in each level until current_level
* Input :  device_ptr<vertex> _full_vertex_array,
           device_ptr<vertex> _full_edge_array,
            vertex _current_vertex,
            int _number_of_vertex,
            int _current_level,
* Result : bool true - if unique value
***********************************************/
struct unique_edge
{
  __host__ __device__
  unique_edge(
    device_ptr<vertex> _full_vertex_array,
    device_ptr<vertex> _full_edge_array,
    vertex _current_vertex,
    int _number_of_vertex,
    int _current_level):   full_vertex_array(_full_vertex_array),
         full_edge_array(_full_edge_array),
         current_vertex(_current_vertex),
         number_of_vertex(_number_of_vertex),
         current_level(_current_level)
  {

  }


  __host__ __device__
    bool operator()(vertex t)
  {
      // Device vector temporal array (candidate)
      if (t == current_vertex)
        return false;

      for (int i= 0; i < current_level; i++)
      {
        int starting = 0;
        if (current_vertex != 0 || i!=0)
        {
          starting = full_vertex_array[i*number_of_vertex + current_vertex - 1];
        }


        int ending = full_vertex_array[i*number_of_vertex + current_vertex];

        bool vertex_previously_found = binary_search(device,
          full_edge_array + starting, full_edge_array + ending, t);
        if  (vertex_previously_found)
          return false;
      }


      return true;
    }

  device_ptr<vertex> full_vertex_array;
  device_ptr<vertex> full_edge_array;
  vertex current_vertex;
  int number_of_vertex;
  int current_level;

};

/**************************************************************
*     Minus data
*    returns current vertex - remove_value
****************************************************************/
struct minus_value
{
  __host__ __device__
  minus_value(int i): remove_value(i)
  {

  }
  __host__ __device__
    vertex operator()(vertex t)
  {
      // Device vector temporal array (candidate)
      if (remove_value > 0)
        return t - remove_value;
      else
        return t;
  }
  int remove_value;
};

/**************************************************************
*     Zero elimenter
*    returns replace_value against zero
****************************************************************/
struct add_offset
{
  __host__ __device__
  add_offset(vertex to_add): value_to_add(to_add)  {  }

  __host__ __device__
    vertex operator()(vertex t)
  {
      // Device vector temporal array (candidate)
      return t + value_to_add;
  }
  vertex value_to_add;
};
