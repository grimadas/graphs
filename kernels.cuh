/*
*   File with device kernels
*   @author = Bulat Nasrulin
*/
#include "headers.h"



/*
*   Count occurance of degree in vertex_degrees and puts data in degree_count.
*   Puts max_degree in value.
*   Input:    device_ptr<vertex> vertex_degrees (from_where),
*             device_ptr<vertex> degree_count (where_to_add),
*              int added_value - offset value.
*
*   Output:   degree_count - with number of occurance
*/
__global__ void degree_count_former
(
  device_ptr<vertex> vertex_degrees,
  device_ptr<vertex> degree_count,
  int edges_to_process,
  int added_value)
{
  	int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < edges_to_process)
    {
      vertex current_pos = vertex_degrees[i];
      vertex* k = raw_pointer_cast(degree_count + current_pos - added_value);
      atomicAdd(k, 1);
    }

}


/*
*	Expand array according to it's from and to values
*	Input : device_ptr<vertex> expanded_array
*			device_ptr<vertex> position_current_level
*			device_ptr<vertex> current_ending_offset
*	Out: 	sorted and normalized expanded_array
*/
__global__  void expander(
	device_ptr<vertex> current_vertex,
  device_ptr<vertex> inter_vertex,
  device_ptr<vertex> temp_from, device_ptr<vertex> temp_to,
	device_ptr<vertex> full_vertex_array, device_ptr<vertex> full_edge_array,
	device_ptr<vertex> position_in_array,
	device_ptr<vertex> expanded_array,
  device_ptr<vertex> inter_vertex_array,
	int number_edges_to_process, int number_of_vertex,
	int current_level,
  device_ptr<vertex> vertex_offset)
{

  for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
             idx < number_edges_to_process;
             idx += blockDim.x * gridDim.x)

  {
    int offset_to_put_exp_array = 0;
  	/*
  	*	For searching in full_edge_list if already discovered
  	*/

  	if (idx != 0)
  	{
  		offset_to_put_exp_array = position_in_array[idx - 1]; // reserved position in expanded array
  	}

  	/*
  	*	Copy to expanded array if the edge is unique (was not discovered previously, not equal to vertex itself)
  	*	Result:			1 2 1 .... (expanded edges)
  	*/
  	  device_ptr<vertex> current_position =
   	  copy_if(device,
    		make_permutation_iterator(full_edge_array, make_counting_iterator<vertex>(temp_from[idx])),
    		make_permutation_iterator(full_edge_array, make_counting_iterator<vertex>(temp_to[idx])),
    		expanded_array + offset_to_put_exp_array,
    		unique_edge(full_vertex_array, full_edge_array,
                    current_vertex[idx],
    		            number_of_vertex, current_level));

    	int planed_size = temp_to[idx] - temp_from[idx];
    	int real_size = distance(expanded_array + offset_to_put_exp_array, current_position);
    	int starting = current_vertex[idx];

      int in_vert = inter_vertex[idx];
      copy(device,
          make_constant_iterator(in_vert),
          make_constant_iterator(in_vert) + real_size,
          inter_vertex_array + offset_to_put_exp_array);

    	vertex* k = raw_pointer_cast(vertex_offset + starting);
    	atomicAdd(k, real_size);
    }
}

/*
*	Expand array according to it's from and to values
*	Input : device_ptr<vertex> expanded_array
*			device_ptr<vertex> position_current_level
*			device_ptr<vertex> current_ending_offset
*	Out: 	sorted and normalized expanded_array
*/
__global__  void simple_expander(
	device_ptr<vertex> current_vertex,
  device_ptr<vertex> temp_from, device_ptr<vertex> temp_to,
	device_ptr<vertex> full_vertex_array, device_ptr<vertex> full_edge_array,
	device_ptr<vertex> position_in_array,
	device_ptr<vertex> expanded_array,
  int number_edges_to_process, int number_of_vertex,
	int current_level,
  device_ptr<vertex> vertex_offset)
{

  for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
             idx < number_edges_to_process;
             idx += blockDim.x * gridDim.x)

  {
    int offset_to_put_exp_array = 0;
  	/*
  	*	For searching in full_edge_list if already discovered
  	*/

  	if (idx != 0)
  	{
  		offset_to_put_exp_array = position_in_array[idx - 1]; // reserved position in expanded array
  	}

  	/*
  	*	Copy to expanded array if the edge is unique (was not discovered previously, not equal to vertex itself)
  	*	Result:			1 2 1 .... (expanded edges)
  	*/
  	  device_ptr<vertex> current_position =
   	  copy_if(device,
    		make_permutation_iterator(full_edge_array, make_counting_iterator<vertex>(temp_from[idx])),
    		make_permutation_iterator(full_edge_array, make_counting_iterator<vertex>(temp_to[idx])),
    		expanded_array + offset_to_put_exp_array,
    		unique_edge(full_vertex_array, full_edge_array,
                    current_vertex[idx],
    		            number_of_vertex, current_level));

    	int planed_size = temp_to[idx] - temp_from[idx];
    	int real_size = distance(expanded_array + offset_to_put_exp_array, current_position);
    	int starting = current_vertex[idx];

    	vertex* k = raw_pointer_cast(vertex_offset + starting);
    	atomicAdd(k, real_size);
    }
}


/*
*	Sorting function for each vertex offset position. Initial sort (not for L > 1)
*	Input: device_ptr full_edge_array
*		   device_ptr full_vertex_array
*	Out : Sorted full_edge_array
*/
__global__ void sorter(device_ptr<vertex> full_edge_array,
						          device_ptr<vertex> full_vertex_array,
                      int level,
<<<<<<< HEAD
                      int number_of_vertex)
=======
                      int process_number)
>>>>>>> origin/Save_Thrust
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < number_of_vertex)
  {
	int starting_point = 0;
	if (idx != 0 || level!=0)
	{
<<<<<<< HEAD
		starting_point = full_vertex_array[level*number_of_vertex + idx - 1];
	}
	int ending_point = full_vertex_array[level*number_of_vertex + idx];
=======
		starting_point = full_vertex_array[level*graph.number_of_vertex + idx - 1];
	}
	int ending_point = full_vertex_array[level*graph.number_of_vertex + idx];
>>>>>>> origin/Save_Thrust
	// sort
	sort(device, full_edge_array + starting_point, full_edge_array + ending_point);
  }
}


/*

*	Sort and remove duplicates in edge_array
*	Input : device_ptr<vertex> expanded_array
*			device_ptr<vertex> position_current_level
*			device_ptr<vertex> current_ending_offset
*	Out: 	sorted and normalized expanded_array
*/
__global__ void unifier(
	device_ptr<vertex> expanded_array,
  device_ptr<vertex> expanded_array_copy,
  device_ptr<vertex> temp_expanded,
  device_ptr<vertex> from_vertex_array,
  device_ptr<vertex> inter_vertex_array,
  device_ptr<vertex> inter_vertex_array_copy,
  device_ptr<int> indices,
	device_ptr<int> positions_vertex_current_level,
	device_ptr<int> current_ending_offset,
  int N)
{
		int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx < N)
    {
      int start_point = 0;

  		if (idx != 0)
  		{
  			start_point = positions_vertex_current_level[idx - 1];

<<<<<<< HEAD
  		}
  		int end_point = positions_vertex_current_level[idx];
      //  printf("ZZZ = s", start_point, end_point);
      copy(device, make_constant_iterator(idx) + start_point,
                   make_constant_iterator(idx) + end_point, from_vertex_array + start_point);
      // ----------------------------- Sorting --------------------------
      sequence(device, indices + start_point, indices + end_point);
      sort_by_key(device, expanded_array + start_point, expanded_array + end_point, indices + start_point);
      // remove dublicates
      // Now reorder the ID arrays using the sorted indices

      thrust::gather(device, indices + start_point, indices + end_point, inter_vertex_array + start_point, inter_vertex_array_copy + start_point);
      thrust::gather(device, indices + start_point, indices + end_point, temp_expanded + start_point, expanded_array_copy + start_point);
      // ---------------------------- Unification ------------------------------

      device_ptr<vertex> current_position =
  		unique(device, expanded_array + start_point, expanded_array + end_point);

      int real_size = distance(expanded_array + start_point, current_position);
  	  current_ending_offset[idx] = real_size;
=======
		}
		int end_point = positions_vertex_current_level[idx];

    sequence(device, indices + start_point, indices + end_point);
    sort_by_key(device, expanded_array + start_point, expanded_array + end_point, indices + start_point);

    // remove dublicates
    // Now reorder the ID arrays using the sorted indices
    thrust::gather(device, indices + start_point, indices + end_point, from_vertex_array + start_point, from_vertex_array_copy + start_point);
    thrust::gather(device, indices + start_point, indices + end_point, inter_vertex_array + start_point, inter_vertex_array_copy + start_point);
    thrust::gather(device, indices + start_point, indices + end_point, temp_expanded + start_point, expanded_array_copy + start_point);
    //thrust::copy(device, expanded_array + start_point, expanded_array + end_point, expanded_array_copy + start_point);
    copy(device, from_vertex_array_copy + start_point, from_vertex_array_copy + end_point, from_vertex_array + start_point);

    // device_ptr<vertex> current_position =
    // 	unique(device, expanded_array + start_point, expanded_array + end_point);
      //thrust::tuple<thrust::device_ptr<int>, thrust::device_ptr<int> > current_position =
    int* l1 = raw_pointer_cast(expanded_array + start_point);
    int* l2 = raw_pointer_cast(from_vertex_array + start_point);

   thrust::pair<int*,int*> new_end =
      unique_by_key(device,
          l1, l1 + end_point - start_point, l2);
    vertex* k = raw_pointer_cast(expanded_array + start_point);
    //vertex real_size = distance(expanded_array + start, current_position);
    vertex real_size = distance(l1, new_end.first);

    //int real_size = end_point - start_point;
    current_ending_offset[idx] = real_size;
>>>>>>> origin/Save_Thrust


    }
}


__global__ void simple_unifier(
  device_ptr<vertex> expanded_array,
  device_ptr<int> positions_vertex_current_level,
	device_ptr<int> current_ending_offset,
  int N)
{
		int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx < N)
    {
      int start_point = 0;

  		if (idx != 0)
  		{
  			start_point = positions_vertex_current_level[idx - 1];

  		}
  		int end_point = positions_vertex_current_level[idx];
      // ----------------------------- Sorting --------------------------

      sort(device, expanded_array + start_point, expanded_array + end_point);
      // ---------------------------- Unification ------------------------------

      device_ptr<vertex> current_position =
  		  unique(device, expanded_array + start_point, expanded_array + end_point);

      int real_size = distance(expanded_array + start_point, current_position);
  	  current_ending_offset[idx] = real_size;


    }
}



  __global__ void temp_edge_copier(
    	device_ptr<vertex> expanded_array,
    	device_ptr<vertex> positions_vertex_current_level,
    	device_ptr<vertex> current_ending_offset,
      device_ptr<int> full_vertex_array,
    	device_ptr<vertex> from_array,
    	device_ptr<vertex> to_array,
    	int number_of_vertex, int L_VALUE )
      {
      	int idx = blockIdx.x*blockDim.x + threadIdx.x;
       	 if (idx < number_of_vertex)
        	{
          int start_point = 0;
          if (idx != 0)
        	{
        		start_point = positions_vertex_current_level[idx - 1];

        	}

        	int end_point = start_point + current_ending_offset[idx];
        	int edge_put_list_start = full_vertex_array[L_VALUE *number_of_vertex + idx - 1] - full_vertex_array[L_VALUE *number_of_vertex - 1];

        	copy(device,
            make_zip_iterator(make_tuple(make_constant_iterator(idx) + start_point, expanded_array + start_point)),
            make_zip_iterator(make_tuple(make_constant_iterator(idx) + end_point, expanded_array + end_point)),
            make_zip_iterator(make_tuple(from_array + edge_put_list_start, to_array + edge_put_list_start)));
        }
      }



  __global__ void temp_edge_copier(
  	device_ptr<vertex> expanded_array,
    device_ptr<vertex> from_vertex_array,
  	device_ptr<vertex> positions_vertex_current_level,
  	device_ptr<vertex> current_ending_offset,
    device_ptr<int> full_vertex_array,
  	device_ptr<vertex> from_array,
  	device_ptr<vertex> to_array,
  	int number_of_vertex, int L_VALUE )

  {
  	int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx < number_of_vertex)
    {
      int start_point = 0;
      if (idx != 0)
    	{
    		start_point = positions_vertex_current_level[idx - 1];

    	}

    	int end_point = start_point + current_ending_offset[idx];
    	int edge_put_list_start = full_vertex_array[L_VALUE *number_of_vertex + idx - 1] - full_vertex_array[L_VALUE *number_of_vertex - 1];

    	copy(device,
        make_zip_iterator(make_tuple(from_vertex_array + start_point, expanded_array + start_point)),
        make_zip_iterator(make_tuple(from_vertex_array + end_point, expanded_array + end_point)),
        make_zip_iterator(make_tuple(from_array + edge_put_list_start, to_array + edge_put_list_start)));
    }
  }

__global__ void edge_copier(
	device_ptr<vertex> expanded_array,
	device_ptr<vertex> positions_vertex_current_level,
	device_ptr<vertex> current_ending_offset,
	device_ptr<vertex> full_vertex_array,
	device_ptr<vertex> full_edge_array,
	int L_VALUE,
	int number_of_vertex )
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < number_of_vertex)
  {
    int start_point = 0;
    if (idx != 0)
  	{
  		start_point = positions_vertex_current_level[idx - 1];

  	}

  	int end_point = start_point + current_ending_offset[idx];

  	int edge_put_list_start = full_vertex_array[L_VALUE *number_of_vertex + idx - 1];

  	copy(device, expanded_array + start_point, expanded_array + end_point, full_edge_array + edge_put_list_start);
  }
}

<<<<<<< HEAD





__global__ void opacity_former(
      device_ptr<vertex> from,
      device_ptr<vertex> to,
      device_ptr<int> opacity_index,
      device_ptr<vertex> degree_count,
      device_ptr<opacity> opacity_matrix,
      device_ptr<opacity> lessL_matrix,
      device_ptr<vertex> from_vertex_array,
      device_ptr<vertex> violater,
      device_ptr<int> new_degrees,
      int except_state,
      int max_degree,
      float threshold,
      int N)
{
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if (i < N)
        {
        int index = max_degree*(from[i] - 1) + (to[i] - 1);
        opacity min = degree_count[from[i] - 1] * degree_count[to[i] - 1];
        opacity* k2 = raw_pointer_cast(lessL_matrix + index);
     		if (from[i]  == to[i])
    		{
          opacity k = degree_count[from[i]-1];
          min = k * (k-1)/2;
          atomicAdd(k2, 1.0/2.0);
        }
        min = min * 2.0;
    		opacity* k = raw_pointer_cast(opacity_matrix + index);
        opacity_index[i] = index;
    		opacity added_value = 1.0/ (min);
        //  if (1.0 - *k > 0.001)
          atomicAdd(k, added_value);
          atomicAdd(k2, 1.0/2.0);

          if (*k > threshold)
          {
            // Can be zero - so subsrit 1
            if (violater[i] == 0)
                violater[i] = except_state;
            violater[i] = -violater[i];
            atomicAdd(k, -added_value);
            atomicAdd(k2, -1.0/2.0);
            int* kz = raw_pointer_cast(new_degrees + from_vertex_array[i]);
            atomicAdd(kz, -1);
          }

        }
}












=======
__global__ void form_remove_array(
  device_ptr<unsigned int> start_points,  device_ptr<unsigned int> end_points,
  device_ptr<vertex> temp_expanded,
  device_ptr<vertex> from_vertex_array_copy,
  device_ptr<vertex> inter_vertex_array_copy,
  device_ptr<vertex> expanded_array_copy,
  device_ptr<int>  full_vertex_array,
  device_ptr<vertex> full_edge_array,
  int N)
  )
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < N)
    {
      unsigned int start_point = start_points[idx];
      unsigned int end_point = end_points[idx];
      int exp_remove = temp_expanded[idx];
      thrust::pair< device_ptr<int>, device_ptr<int> > search_seq =
        equal_range(device, expanded_array_copy + start_point,
                            expanded_array_copy + end_point, exp_remove);
      int start_offset = distance(expanded_array_copy + start, search_seq.first);
      int  ending_offset = distance(expanded_array_copy +start, search_seq.second);
      __syncthreads();
      transform(device, expanded_array_copy + start_offset,
                        expanded_array_copy + ending_offset,
                        expanded_array_copy + start_offset,
                        to_negative_minus_one());

    }
  }
>>>>>>> origin/Save_Thrust

__global__ void opacity_former(
	device_ptr<vertex> from,
	device_ptr<vertex> to,
  device_ptr<int> opacity_index,
  device_ptr<vertex> degree_count,
	device_ptr<opacity> opacity_matrix,
  device_ptr<opacity> lessL_matrix,
  int max_degree,
  int N)
	{

		int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < N)
    {
    int index = max_degree*(from[i] - 1) + (to[i] - 1);
    opacity min = degree_count[from[i] - 1] * degree_count[to[i] - 1];
    opacity* k2 = raw_pointer_cast(lessL_matrix + index);
 		if (from[i]  == to[i])
		{
      opacity k = degree_count[from[i]-1];
      min = k * (k-1)/2;
      atomicAdd(k2, 1.0/2.0);
    }
    min = min * 2.0;
		opacity* k = raw_pointer_cast(opacity_matrix + index);
    opacity_index[i] = index;
		opacity added_value = 1.0/ (min);
    //  if (1.0 - *k > 0.001)
    atomicAdd(k, added_value);
    atomicAdd(k2, 1.0/2.0);
    }

	}



    __global__ void opacity_former(
      device_ptr<vertex> from,
      device_ptr<vertex> to,
      device_ptr<int> opacity_index,
      device_ptr<vertex> degree_count,
      device_ptr<opacity> opacity_matrix,
      device_ptr<opacity> lessL_matrix,
      device_ptr<vertex> from_vertex_array,
      device_ptr<vertex> violater,
      device_ptr<int> new_degrees;
      int max_degree,
      float threshold,
      int N)
    	{
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if (i < N)
        {
        int index = max_degree*(from[i] - 1) + (to[i] - 1);
        opacity min = degree_count[from[i] - 1] * degree_count[to[i] - 1];
        opacity* k2 = raw_pointer_cast(lessL_matrix + index);
     		if (from[i]  == to[i])
    		{
          opacity k = degree_count[from[i]-1];
          min = k * (k-1)/2;
          atomicAdd(k2, 1.0/2.0);
        }
        min = min * 2.0;
    		opacity* k = raw_pointer_cast(opacity_matrix + index);
        opacity_index[i] = index;
    		opacity added_value = 1.0/ (min);
        //  if (1.0 - *k > 0.001)
          atomicAdd(k, added_value);
          atomicAdd(k2, 1.0/2.0);

          if (*k > threshold)
          {
            // Can be zero - so subsrit 1
            violater[i] = -violater[i] - 1;
            atomicAdd(k, -added_value);
            atomicAdd(k2, -1.0/2.0);
            opacity* kz = raw_pointer_cast(new_degrees + from_vertex_array[i]);
            atomicAdd(kz, -1);
          }

        }


    	}





  __global__ void remove_by_index(
                      device_ptr<int> index_array,
                      device_ptr<int> degree_count,
                      device_ptr<opacity> opacity_matrix,
                      device_ptr<opacity> lessL_matrix,
                      int max_degree,
                      device_ptr<vertex> from,
                      device_ptr<vertex> to,
                      device_ptr<vertex> opacity_index)
      {
          int idx = blockIdx.x*blockDim.x + threadIdx.x;
          int removal = index_array[idx];

          int from_degree = (opacity_index[removal] / max_degree) + 1;
          int to_degree =  (opacity_index[removal] % max_degree) + 1;

          from_degree = from_degree < to_degree?from_degree:to_degree;
          to_degree = to_degree > (opacity_index[removal] / max_degree) + 1?
                      to_degree:(opacity_index[removal] / max_degree) + 1;

          __syncthreads();
          opacity min = degree_count[from_degree - 1] * degree_count[to_degree - 1];

           if (from_degree  == to_degree)
           {
             opacity k = degree_count[from_degree-1];
             min = k * (k-1);

            }

          //min = min * 2.0;
          opacity added_value = 1.0/ (min);
          opacity* opacity_matrix_pointer = raw_pointer_cast(opacity_matrix + opacity_index[removal]);
          opacity* lessL_matrix_pointer = raw_pointer_cast(lessL_matrix + opacity_index[removal]);

          atomicAdd(opacity_matrix_pointer, -added_value);
          atomicAdd(lessL_matrix_pointer, -1);

          from[removal] = -1;
          to[removal] = -1;
          opacity_index[removal] = -1;

      }



        __global__ void coo_unifier(
        	device_ptr<vertex> from,
        	device_ptr<vertex> to,
        	device_ptr<vertex> opacity_index,
          device_ptr<vertex> real_from,
          device_ptr<vertex> real_to,
          device_ptr<vertex> real_opacity,
          int size,
          int N)
        	{
        		int idx = blockIdx.x*blockDim.x + threadIdx.x;
            if (idx < N)
            {
              device_ptr<int> start = lower_bound(device, from, from + size, idx);
              device_ptr<int> end = upper_bound(device, from, from + size, idx);

              int starting = distance(from, start);
              int ending = distance(from, end);
              //printf("REal_size  %d %d\n", starting,idx);

              stable_sort_by_key(device,
                      to + starting, to + ending, opacity_index + starting);

            	// remove dublicates

            thrust::pair<device_ptr<vertex>, device_ptr<vertex> > new_end =
            	unique_by_key(device, to + starting, to + ending, opacity_index + starting);
            int real_size = distance(to+starting, new_end.first);
          //  printf(" ");
            copy(device,
              make_zip_iterator(make_tuple(from + starting, to + starting, opacity_index + starting)),
              make_zip_iterator(make_tuple(from + starting + real_size,
                                           to + starting + real_size,
                                           opacity_index + starting + real_size)),
              make_zip_iterator(make_tuple(real_from + starting, real_to + starting, real_opacity + starting)));



            }
        	}


    __device__ int remove_candidate()
    {
      // 0 -remove from inter_vertex expanded value
      return 0;
    }

    __device__ void remove_from_edge_array(
      device_ptr<int>   from_vertex_remove,
      device_ptr<vertex> expanded_vertex_remove,
      device_ptr<int> full_vertex_array,
      device_ptr<vertex> full_edge_array,
      device_ptr<vertex> from_vertex_array_copy,
      device_ptr<vertex> inter_vertex_array_copy,
      device_ptr<vertex> expanded_array_copy,
      device_ptr<int> start_points,
      device_ptr<int> end_points,
      int number_of_vertex,
      int current_level, int cur_vert,int N)
    {
      int idx = blockIdx.x*blockDim.x + threadIdx.x;
      if (idx < N)
      {
        int start_offset = 0;
        if (cur_vert != 0 || current_level!=1)
        {
            start_offset = full_vertex_array[(current_level-1)*number_of_vertex + cur_vert - 1];

        }
        int ending_offset = full_vertex_array[(current_level-1)*number_of_vertex + cur_vert];
        int start_point = start_points[idx];
        int end_point = end_points[idx];
        // Find all values to remove (could be more than 1)
        thrust::pair< device_ptr<int>, device_ptr<int> > search_seq =
          equal_range(device, expanded_array_copy + start_point, expanded_array_copy + end_point, exp_remove);
        int start_offset = distance(expanded_array_copy +start, search_seq.first);
        int  ending_offset = distance(expanded_array_copy +start, search_seq.second);
        for (int removing_offset = start_offset; removing_offset < ending_offset; removing_offset++)
        {
          // Here could be  possibilities
          if (remove_candidate == 0)
          {

            int inter_vertex = inter_vertex_array_copy[removing_offset];
            int expanded_vertex = expanded_array_copy[removing_offset];
            // First level
            int start_edge_offset = 0;
            if(inter_vertex != 0)
            {
              start_edge_offset = full_vertex_array[inter_vertex-1];
            }
            int ending_offset = full_vertex_array[inter_vertex];


    __global__ void first_level_remove_via_sort(
      device_ptr<int> full_vertex_array,
      device_ptr<vertex> full_edge_array,
      device_ptr<vertex> from_vertex_array_copy,
      device_ptr<vertex> inter_vertex_array_copy,
      device_ptr<vertex> expanded_array_copy,
      device_ptr<int> new_vertex_array,
      int N)
      {
        int idx = blockIdx.x*blockDim.x + threadIdx.x;
        if (idx < N)
        {
          int inter_vertex = inter_vertex_array_copy[idx];
          int expanded_vertex = expanded_array_copy[idx];

          int start_offset = 0;
          if (inter_vertex != 0)
          {
              start_offset = full_vertex_array[inter_vertex - 1];

          }
          int ending_offset = full_vertex_array[inter_vertex];
          device_ptr<int> start = lower_bound(device,
                                    full_edge_array + start_offset,
                                    full_edge_array + ending_offset, expanded_vertex);

          int index_to_remove = distance(full_edge_array + start_offset,
                                                              start);
          int* k = raw_pointer_cast(new_vertex_array + inter_vertex);
          atomicSub(k, 1);
          __syncthreads();
          full_edge_array[start_offset + index_to_remove] = -1;
          sort(device, full_edge_array + start_offset,
                        full_edge_array + ending_offset);
        }

      }


    __global__ void first_level_remove(
      device_ptr<int> full_vertex_array,
      device_ptr<vertex> full_edge_array,
      device_ptr<vertex> from_vertex_array_copy,
      device_ptr<vertex> inter_vertex_array_copy,
      device_ptr<vertex> expanded_array_copy,
      device_ptr<int> new_vertex_array,
      int N)
    {
      int idx = blockIdx.x*blockDim.x + threadIdx.x;
      if (idx < N)
      {
        int inter_vertex = inter_vertex_array_copy[idx];
        int expanded_vertex = expanded_array_copy[idx];

        int start_offset = 0;
        if (inter_vertex != 0)
        {
            start_offset = full_vertex_array[inter_vertex - 1];

        }
        int ending_offset = full_vertex_array[inter_vertex];
        device_ptr<int> start = lower_bound(device,
                                  full_edge_array + start_offset,
                                  full_edge_array + ending_offset, expanded_vertex);

        int index_to_remove = distance(full_edge_array + start_offset,
                                                            start);
        int* k = raw_pointer_cast(new_vertex_array + inter_vertex);
        atomicSub(k, 1);
        __syncthreads();
        full_edge_array[start_offset + index_to_remove] = -1;


      }
    }


    __global__ void edge_remove_by_vertex(
      device_ptr<vertex> expanded_array_copy,
      device_ptr<vertex> start_points,
      device_ptr<vertex> end_points,
      device_ptr<int> vertex_index_array,
      device_ptr<int> full_vertex_array,
      device_ptr<vertex> full_edge_array,
      device_ptr<int> new_vertex_array,
      int number_of_vertex,
      int current_level, int N)
    {

      int idx = blockIdx.x*blockDim.x + threadIdx.x;
      if (idx < N)
      {
        int expanded_vertex = expanded_array_copy[idx];
        int removing_vertexes[end_points[idx] - start_points[idx]];
        for (int i=start_points[idx]; i < end_points[idx]; i++)
        {
          removing_vertexes[i - start_points[idx]] = -1;
          int found_vertex = vertex_index_array[i];
          int start_offset  = full_vertex_array[number_of_vertex * (current_level + 1) - 1 + found_vertex];

          int ending_offset = full_vertex_array[number_of_vertex * (current_level + 1) + found_vertex];
          device_ptr<int> start = lower_bound(device,
                                    full_edge_array + start_offset,
                                    full_edge_array + ending_offset, expanded_vertex);
          int index_to_remove = distance(full_edge_array + start_offset,
                                                              start);
          if (full_edge_array[start_offset + index_to_remove] == expanded_vertex)
          {
            removing_vertexes[i - start_points[idx]] = start_offset + index_to_remove;
          }

        }
        __syncthreads();
        for (int i= 0 ; i < end_points[idx] - start_points[idx]; i++)
        {
          if (removing_vertexes != -1)
          {
              // remove procedure
              full_edge_array[removing_vertexes[i]] = -1;
              int* k = raw_pointer_cast(new_vertex_array + vertex_index_array[i + start_points[idx]]);
              atomicSub(k, 1);

          }

        }


          // Here could be  possibilities
          if (remove_candidate == 0)
          {

            int inter_vertex = inter_vertex_array_copy[removing_offset];
            int expanded_vertex = expanded_array_copy[removing_offset];
            // First level
            int start_edge_offset = 0;
            if(inter_vertex != 0)
            {
              start_edge_offset = full_vertex_array[inter_vertex-1];
            }
            int ending_offset = full_vertex_array[inter_vertex];
            device_ptr<int> start = lower_bound(device,
                                      full_edge_array + start_offset,
                                      full_edge_array + ending_offset, expanded_vertex);


<<<<<<< HEAD
__global__ void form_remove_array(
          device_ptr<vertex> expanded_array_copy,
          device_ptr<int> positions_vertex_current_level,
          device_ptr<int> temp_from_array,
          device_ptr<vertex> temp_expanded,
          int number_of_vertex,
          int N)
{
              int idx = blockIdx.x*blockDim.x + threadIdx.x;
              if (idx < N)
              {

                int search_vertex = temp_from_array[idx];
                int start_point = 0;
                if (search_vertex!=0)
                  start_point = positions_vertex_current_level[search_vertex - 1];
                int end_point = positions_vertex_current_level[search_vertex];
                int exp_remove = temp_expanded[idx];
                thrust::pair< device_ptr<int>, device_ptr<int> > search_seq =
                  equal_range(device, expanded_array_copy + start_point,
                                      expanded_array_copy + end_point, exp_remove);
                int start_offset = distance(expanded_array_copy + start_point, search_seq.first);
                int  ending_offset = distance(expanded_array_copy +start_point, search_seq.second);
                __syncthreads();
                transform(device, expanded_array_copy + start_point + start_offset,
                                  expanded_array_copy + start_point + ending_offset,
                                  expanded_array_copy + start_point + start_offset,
                                  to_negative(number_of_vertex+1));

              }
}

=======
            int index_to_remove = distance(full_edge_array + start_offset,
                                                                start);

             __syncthreads();
             full_edge_array[start_offset + index_to_remove] = -1;
             int* k = raw_pointer_cast(new_vertex_array + inter_vertex);
             atomicSub(k, 1);



          }

        }

      }
>>>>>>> origin/Save_Thrust

    }

    __global__ void remover(
    device_ptr<int>   from_vertex_remove,
    device_ptr<vertex> expanded_vertex_remove,
    device_ptr<vertex> from_vertex_array,
    device_ptr<vertex> inter_vertex_array,
    device_ptr<vertex> expanded_array,
    device_ptr<int> new_vertex_array,
    int total_size, int N)
  {
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx < N)
    {
        int from_remove = from_vertex_remove[idx];
        int exp_remove = expanded_vertex_remove[idx];
        from start_point
        from ending_point
        thrust::pair< device_ptr<int>, device_ptr<int> > search_seq =
        equal_range(expanded_array_copy + start, expanded_array_copy + end_point, exp_remove, )

        start_offset = distance(expanded_array_copy +start, search_seq.first);
        ending_offset = distance(expanded_array_copy +start, search_seq.second);

        // Two possibilities:
        full_edge_array
        inter_vertex_array_copy + start_offset

        lower_bound(expanded_array_copy, expanded_array_copy )


        device_ptr<vertex> start = lower_bound(device,
                                  from_vertex_array,
                                  from_vertex_array + total_size, expanded_value);


        int from_vertex = from_vertex_array[idx];
        int inter_vertex = inter_vertex_array[idx];
        int expanded_value = expanded_array[idx];
        if (random_remove() == 0)
        {
          int start_offset = 0;
          if(inter_vertex != 0)
          {
            start_offset = full_vertex_array[inter_vertex-1];
          }
          int ending_offset = full_vertex_array[inter_vertex];
          device_ptr<int> start = lower_bound(device,
                                    full_edge_array + start_offset,
                                    full_edge_array + ending_offset, expanded_value);


          int index_to_remove = distance(full_edge_array + start_offset,
                                                              start);
          __syncthreads();
          full_edge_array[start_offset + index_to_remove] = -1;

          int* k = raw_pointer_cast(new_vertex_array + inter_vertex);
          atomicSub(k, 1);
      }
     }
    else{
      int inter_vertex = inter_vertex_array[idx];
      int expanded_value = expanded_array[idx];

      if(random_remove() == 0){
        // TODO: current_level
        for_each( device,
                 full_vertex_array + current_level*number_of_vertex,
                 full_vertex_array + (current_level+1)*number_of_vertex,
                 seek_and_destroy(
                   inter_vertex, expanded_value, current_level,
                   full_vertex_array,
                   full_edge_array,
                   number_of_vertex,
                   new_vertex_array));

        for(int i=0; i<number_of_vertex; i++)
        {
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
                device_ptr<int> temp_ptr = lower_bound(device, full_edge_array + start_offset,
                                                               full_edge_array + ending_offset, what_to_remove);
                // The connection indeed exists
                if(*temp_ptr == what_to_remove) {
                      // index where it should be removed
                      int index_to_remove = distance(full_edge_array + start_offset, temp_ptr);
                      int start_offset_copy = start_offset - full_vertex_array[(current_level+1)*number_of_vertex - 1];
                      full_edge_array_copy[start_offset_copy + index_to_remove] = -1;
                      int* k = raw_pointer_cast(new_vertex_array + vertex_number);
                      atomicAdd(k, -1);
                  }

                }
            }
            __syncthreads();
            sort(device, full_edge_array + start_offset, full_edge_array + ending_offset);
          }
        }
      }




  __global__ via_vertex_search()
  {
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
                      }

  }

        __global__ void edge_remover(
          device_ptr<int> remove_opacity_index,
          device_ptr<int> remove_count,
          device_ptr<int> opacity_index,
          device_ptr<int> index_array,
          int number_of_edges, int N)
        {
          int i = blockIdx.x*blockDim.x + threadIdx.x;
          if (i < N)
          {
            int count = remove_count[i];
            int remove_value = remove_opacity_index[i];

            device_ptr<int> start = lower_bound(device, opacity_index, opacity_index+number_of_edges, remove_value);

            int offset = distance(opacity_index, start);
            for(int i=0; i< count; i++)
              index_array[offset+i] = offset + i;

          }
        }


/**
*   inter_vertex - where find
*   expanded_vertex - what find
**/
__device__ int find_to_remove(
                    device_ptr<vertex> full_edge_array,
                    device_ptr<int> full_vertex_array,
                    vertex inter_vertex, vertex expanded_vertex)
{
  int start_offset = 0;
  if (inter_vertex != 0)
  {
      start_offset = full_vertex_array[inter_vertex - 1];

  }
  int ending_offset = full_vertex_array[inter_vertex];
  device_ptr<int> start = lower_bound(device,
                            full_edge_array + start_offset,
                            full_edge_array + ending_offset, expanded_vertex);

   return start_offset + distance(full_edge_array + start_offset, start);

}

__global__ void first_level_remove(
      device_ptr<int> full_vertex_array,
      device_ptr<vertex> full_edge_array,
      device_ptr<vertex> inter_vertex_array_copy,
      device_ptr<vertex> expanded_array_copy,
      int N)
  {
        int idx = blockIdx.x*blockDim.x + threadIdx.x;
        if (idx < N)
        {
          int inter_vertex = inter_vertex_array_copy[idx];
          int expanded_vertex = expanded_array_copy[idx];

          int direct_remove = find_to_remove(full_edge_array,
                                             full_vertex_array,
                                             inter_vertex, expanded_vertex);
          int inverse_remove = find_to_remove(full_edge_array, full_vertex_array,
                                              expanded_vertex, inter_vertex);
          __syncthreads();
          full_edge_array[direct_remove] = -1;
          full_edge_array[inverse_remove] = -1;

        }

  }

__global__ void calc_new_offsets(
              device_ptr<int> full_vertex_array,
              device_ptr<vertex> full_edge_array,
              device_ptr<int> new_vertex_array,
              int number_of_vertex)
  {
    for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
               idx < number_of_vertex;
               idx += blockDim.x * gridDim.x)
      {
        int start_offset = 0;
        if (idx != 0)
        {
            start_offset = full_vertex_array[idx - 1];

        }
        int ending_offset = full_vertex_array[idx];
        sort(device, full_edge_array + start_offset,
                     full_edge_array + ending_offset);
        int added_value = count(device, full_edge_array + start_offset,
                      full_edge_array + ending_offset, -1);
        int* k = raw_pointer_cast(new_vertex_array + idx);
        atomicAdd(k, -added_value);
      }

  }


        /**
        *
        *   Input: - Expanded vertex value - what should be found
        *          - start_points - starting offset where to search
        *          - end_points - ending offset where to search
        *          - vertex_index_array - vertex index where inter_vertex was found
        *          - full_vertex_array, full_edge_array - graph properties
        *          - new_vertex_array - new array degree offset
        *          - current_level - current search level, N - number of to process
        *     Total time:  N/p * d * (log(d) + 1)
        *     Total memory: O(N*d)
        *     Functionality:
        *              For each offset :
        *                 create array T[](O(max_degree)?);
        *                 find in graph(level + 1) expanded value
        *                 if found: put in array T
        *              __sync
        *              For each t in  T:
        *                    if t not null - mark edge to delete
        *
        **/
      __global__ void edge_remove_by_vertex(
        device_ptr<vertex> expanded_array_copy,
        device_ptr<unsigned int> start_points,
        device_ptr<unsigned int> end_points,
        device_ptr<int> vertex_index_array,
        device_ptr<int> removing_offsets,
        device_ptr<int> removing_vertexes,
        device_ptr<int> full_vertex_array,
        device_ptr<vertex> full_edge_array,
        device_ptr<int> new_vertex_array,
        int number_of_vertex,
        int current_level,

        int N)
      {

        int idx = blockIdx.x*blockDim.x + threadIdx.x;
        if (idx < N)
        {
          int expanded_vertex = expanded_array_copy[idx];
          const int total_size = end_points[idx] - start_points[idx];
          //int* const l = &total_size;
          int added_offset = 0;
          if (idx!=0)
              added_offset = removing_offsets[idx-1];

          for (int i=start_points[idx]; i < end_points[idx]; i++)
          {

            int found_vertex = vertex_index_array[i];
            int start_offset  = full_vertex_array[number_of_vertex * (current_level + 1) - 1 + found_vertex];

            int ending_offset = full_vertex_array[number_of_vertex * (current_level + 1) + found_vertex];
            device_ptr<int> start = lower_bound(device,
                                      full_edge_array + start_offset,
                                      full_edge_array + ending_offset, expanded_vertex);
            int index_to_remove = distance(full_edge_array + start_offset,
                                                                start);
            if (full_edge_array[start_offset + index_to_remove] == expanded_vertex)
            {
              removing_vertexes[added_offset + i - start_points[idx]] = start_offset + index_to_remove;
            }

          }
          __syncthreads();
          for (int i= 0 ; i < end_points[idx] - start_points[idx]; i++)
          {
            if (removing_vertexes[added_offset + i] != -1)
            {
                // remove procedure
                full_edge_array[removing_vertexes[added_offset + i]] = -1;
                int* k = raw_pointer_cast(new_vertex_array + vertex_index_array[i + start_points[idx]]);
                atomicAdd(k, -1);

            }
          }


          }
      }




  // This example demonstrates how to expand an input sequence by
  // replicating each element a variable number of times. For example,
  //
  //   expand([2,2,2],[A,B,C]) -> [A,A,B,B,C,C]
  //   expand([3,0,1],[A,B,C]) -> [A,A,A,C]
  //   expand([1,3,2],[A,B,C]) -> [A,B,B,B,C,C]
  //
  // The element counts are assumed to be non-negative integers

  template <typename InputIterator1,
            typename InputIterator2,
            typename OutputIterator>  OutputIterator expand(InputIterator1 first1,
                                      InputIterator1 last1,
                                      InputIterator2 first2,
                                      OutputIterator output)
  {
    typedef typename iterator_difference<InputIterator1>::type difference_type;

    difference_type input_size  = distance(first1, last1);
    difference_type output_size = reduce(device, first1, last1);

    // scan the counts to obtain output offsets for each input element
    device_vector<difference_type> output_offsets(input_size, 0);
    exclusive_scan(device, first1, last1, output_offsets.begin());

    // scatter the nonzero counts into their corresponding output positions
    device_vector<difference_type> output_indices(output_size, 0);
    scatter_if
      (device, counting_iterator<difference_type>(0),
       counting_iterator<difference_type>(input_size),
       output_offsets.begin(),
       first1,
       output_indices.begin());

    // compute max-scan over the output indices, filling in the holes
    inclusive_scan
      (device, output_indices.begin(),
       output_indices.end(),
       output_indices.begin(),
       maximum<difference_type>());

    // gather input values according to index array (output = first2[output_indices])
    OutputIterator output_end = output; advance(output_end, output_size);
    gather(device, output_indices.begin(),
                   output_indices.end(),
                   first2,
                   output);

    // return output + output_size
    advance(output, output_size);
    return output;
  }
