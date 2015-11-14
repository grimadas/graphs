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
	device_ptr<vertex> expanded_array, device_ptr<vertex> from_vertex_array,
  device_ptr<vertex> inter_vertex_array,
	int number_edges_to_process, int number_of_vertex,
	int current_level,
  device_ptr<vertex> vertex_offset)
{

	int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < number_edges_to_process)
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
    	// TODO : check real size value
    	/*
    	*	Expand the current processing vertex to the size *real size*;
    	*			Result : 0 0 0 1 1 ... (the vertex from expanded)
    	*/
    	copy(device,
    		make_constant_iterator(starting),
    		make_constant_iterator(starting) + real_size,
    		from_vertex_array + offset_to_put_exp_array);

      int in_vert = inter_vertex[idx];
      copy(device,
          make_constant_iterator(in_vert),
          make_constant_iterator(in_vert) + real_size,
          inter_vertex_array + offset_to_put_exp_array);


    	vertex* k = raw_pointer_cast(vertex_offset + starting);
    	atomicAdd(k, real_size);
      /*
    	if (planed_size != real_size)
    		if (idx != 0)
    	{
    			transform(device, position_in_array + idx - 1,
    			position_in_array + number_edges_to_process,
    			position_in_array + idx -1, minus_value(planed_size - real_size));
    	}
  		else
  		{
  			transform(device, position_in_array,
  				position_in_array + number_edges_to_process, position_in_array, minus_value(planed_size - real_size));
  		}
      */
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
                      int process_number)
{
	int idx = blockIdx.x*blockDim.x + threadIdx.x;
  if (idx < process_number)
  {
	int starting_point = 0;
	if (idx != 0)
	{
		starting_point = full_vertex_array[idx - 1];
	}
	int ending_point = full_vertex_array[idx];
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
  device_ptr<vertex> from_vertex_array,
  device_ptr<vertex> from_vertex_array_copy,
  device_ptr<vertex> inter_vertex_array,
  device_ptr<vertex> inter_vertex_array_copy,
  device_ptr<int> indices,
	device_ptr<vertex> positions_vertex_current_level,
	device_ptr<vertex> current_ending_offset,
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
    //  printf("ZZZ = s", start_point, end_point);
    //sort(device, expanded_array + start_point, expanded_array + end_point);
    sequence(device, indices + start_point, indices + end_point);
    sort_by_key(device, expanded_array + start_point, expanded_array + end_point, indices + start_point);
    // remove dublicates
    // Now reorder the ID arrays using the sorted indices
    thrust::gather(device, indices + start_point, indices + end_point, from_vertex_array + start_point, from_vertex_array_copy + start_point);
    thrust::gather(device, indices + start_point, indices + end_point, inter_vertex_array + start_point, inter_vertex_array_copy + start_point);

    device_ptr<vertex> current_position =
		unique(device, expanded_array + start_point, expanded_array + end_point);

      vertex real_size = distance(expanded_array + start_point, current_position);
		  current_ending_offset[idx] = real_size;


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
		// TODO: Not sure about this
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
            //  printf("REal_size  %d %d\n", starting,idx);

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
