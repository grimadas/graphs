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
  thrust::device_vector<float> Similarities(8192*4320);
  thrust::device_vector<float> Similarities2(8192*4320);
  thrust::device_vector<float> Similarities3(2*8192*4320);
  Similarities[0] = 222222;
  Similarities[9830400]= 222223;
  std::cout << "Similarities[0] = " << Similarities[0] << std::endl;
  std::cout << "Similarities[9380400] = " << Similarities[9830400] << std::endl;

  thrust::merge(Similarities.begin(), Similarities.end(), Similarities2.begin(), Similarities2.end(), Similarities3.begin());
  std::cout << "Similarities3[9380400] = " << Similarities3[9830400] << std::endl;

  return 0;
}
