#include <thrust/device_malloc.h>
#include <thrust/device_ptr.h>
#include <thrust/copy.h>
#include <iostream>

using namespace thrust;

__global__ void someting(thrust::device_ptr<int>* a, device_ptr<int>* b)
{
  	int i = blockIdx.x*blockDim.x + threadIdx.x;
  a[i] = device_malloc<int>(20 * 20);
  b[i] = a[i];

}

int main()
{
    // Initial host data
    device_ptr<int>* tempo_vector = new device_ptr<int>[20];
    device_ptr<int>* previous = new device_ptr<int>[20];
    device_ptr<int>* current = new device_ptr<int>[20];
  
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
