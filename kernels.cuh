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
	device_ptr<vertex> current_vertex, device_ptr<vertex> temp_from, device_ptr<vertex> temp_to,
	device_ptr<vertex> full_vertex_array, device_ptr<vertex> full_edge_array,
	device_ptr<vertex> position_in_array,
	device_ptr<vertex> expanded_array, device_ptr<vertex> from_vertex_array,
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
    		unique_edge(full_vertex_array, full_edge_array, current_vertex[idx],
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
		if (end_point > start_point)
		{
      sort(device, expanded_array + start_point, expanded_array + end_point);
			// remove dublicates
			device_ptr<vertex> current_position =
				unique(device, expanded_array + start_point, expanded_array + end_point);
			vertex real_size = distance(expanded_array + start_point, current_position);
			current_ending_offset[idx] = real_size;
		}
		else
		{
			current_ending_offset[idx] = 0;
		}
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
	device_ptr<vertex> degree_count,
	device_ptr<opacity> opacity_matrix,
	int max_degree,
  int N)
	{
		// TODO: Not sure about this
		int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < N)
    {

    opacity min = degree_count[from[i] - 1] * degree_count[to[i] - 1];
 		if (from[i]  == to[i])
		{
      opacity k = degree_count[from[i]-1];
      min = k * (k-1)/2;
    }
      min = min * 2.0;
		opacity* k = raw_pointer_cast(opacity_matrix + max_degree*(from[i] - 1) + (to[i] - 1));
		opacity added_value = 1.0/ (min);
  //  if (1.0 - *k > 0.001)
      atomicAdd(k, added_value);
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
