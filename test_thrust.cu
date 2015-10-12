#include <thrust/device_malloc.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <iostream>
#include "headers.h"

using namespace thrust;
using namespace std;



struct unique_edge
{

  __host__ __device__
    bool operator()(vertex t)
  {
      // Device vector temporal array (candidate)
      if (t == 1)
        return false;
      return true;
  }

};

struct cleanup
{
  __host__ __device__
  cleanup(int i): remove_value(i)
  {

  }
  __host__ __device__
    vertex operator()(vertex t)
  {
      // Device vector temporal array (candidate)
      return t - remove_value;
  }
  int remove_value;
};





__global__ void someting(device_ptr<int> previos,  device_ptr<int> current,
                        device_ptr<int> full_edge_array, device_ptr<int> temp_from, device_ptr<int> temp_to,  device_ptr<int> current2,
                        device_ptr<int> vertex_array, device_ptr<int> position
      )
{
  	int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int offset = 0;
    if (idx != 0)
    {
      offset = position[idx - 1];
    }
    thrust::device_ptr<int> cur_pos =
      thrust::copy_if(thrust::device,
      thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_from[idx])),
      thrust::make_permutation_iterator(full_edge_array, thrust::make_counting_iterator<vertex>(temp_to[idx])),
      current + offset,
      unique_edge());

      //cur_pos =
      //thrust::remove(thrust::device, current + offset, cur_pos, -1);

      int planed_size = temp_to[idx] - temp_from[idx];
      int real_size = thrust::distance(current+offset, cur_pos);
      int starting = vertex_array[idx];
      thrust::copy(thrust::device,
      thrust::make_constant_iterator(starting),
      thrust::make_constant_iterator(starting) + real_size,
      current2 + offset);

            // reformat offset array
      if (planed_size != real_size)
        if (idx!=0)
        {thrust::transform(thrust::device, position + idx-1, position + 5, position + idx -1, cleanup(planed_size - real_size));}
        else
        {thrust::transform(thrust::device, position, position + 5, position, cleanup(planed_size - real_size));}
          __syncthreads();
        offset = 0;
        if (idx != 0)
        {
          offset = position[idx - 1];
        }





}

int main()
{
    // Initial host data
    int full2[40];
    int full[5] = {1, 2, 3, 4, 5};
    int vert[5] = {0, 0, 1, 1, 2};
    int _temp_from[5] = {0, 0, 1, 2, 3};
    int _temp_to[5] = {1, 2, 3, 4, 4};

    int _position[5] = {1,2,2,2,1};

    device_ptr<int> full_edge_array = device_malloc<int>(10);
    thrust::copy(full, full + 5, full_edge_array);

    device_ptr<int> vertex_array = device_malloc<int>(10);
    thrust::copy(vert, vert + 5, vertex_array);

    device_ptr<int> temp = device_malloc<int>(10);
    thrust::copy(vert, vert + 5, vertex_array);

    device_ptr<int> temp_from = device_malloc<int>(10);
    thrust::copy(_temp_from, _temp_from + 5, temp_from);

    device_ptr<int> temp_to = device_malloc<int>(10);
    thrust::copy(_temp_to, _temp_to + 5, temp_to);
    device_ptr<int> current = thrust::device_malloc<int>(100);
    thrust::fill(thrust::device, current, current + 100, -1);
    device_ptr<int> current2 = thrust::device_malloc<int>(100);

    thrust::device_ptr<int> position = device_malloc<int>(5);
    thrust::copy(_position, _position + 5, position);
    thrust::inclusive_scan(thrust::device, position, position + 5, position);

    thrust::fill(thrust::device, current2, current2 + 100, -1);
    device_ptr<int> current_begin = current;
    device_ptr<int> current_begin2 = current2;
    device_ptr<int> previos = current;
    someting<<<1,5>>>(previos, current,
            full_edge_array, temp_from,
            temp_to, current2,
                      vertex_array, position);

    thrust::copy(position, position + 5, _position);
    thrust::remove(thrust::device, current, current + 20, -1);
    thrust::remove(thrust::device, current2, current2 + 20, -1);

    for (int i=0; i< 5; i++)
    {
      full[i] = -1;
    }
    //thrust::remove(current_begin, current_begin + 40, -1);
    thrust::copy(current_begin2, current_begin2 + _position[4], full2);
    for (int i=0; i< _position[4]; i++)
    {
      cout << full2[i] << " ";
    }
    cout << endl;

    thrust::copy(current_begin, current_begin + _position[4], full2);
    for (int i=0; i< _position[4]; i++)
    {
      cout << full2[i] << " ";
    }
    cout << endl;

    // The slowest part TODO: in kernel
    /*
    for (int i=0; i< graph.number_of_vertex; i++)
    {
      tempo_vector[i] = device_malloc<vertex>(graph.number_of_vertex * graph.number_of_vertex);
      current[i] = tempo_vector[i];
      previous[i] = current[i];
    }
    */
    return 0;
}
