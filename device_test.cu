#include <stdio.h>

__global__ void MyKernel(unsigned int *array, unsigned int arrayCount)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < arrayCount)
  {
    array[idx] *= array[idx];
  }
}

void launchMyKernel(unsigned int *array, unsigned int arrayCount)
{
  int blockSize;   // The launch configurator returned block size
  int minGridSize; // The minimum grid size needed to achieve the
                   // maximum occupancy for a full device launch
  int gridSize;    // The actual grid size needed, based on input size

  cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
                                      MyKernel, 0, 0);
  // Round up according to array size
  gridSize = (arrayCount + blockSize - 1) / blockSize;
  printf("GridSize %d. Block size: %d\n",
         gridSize,blockSize);
  MyKernel<<< gridSize, blockSize >>>(array, arrayCount);

  cudaDeviceSynchronize();

  // calculate theoretical occupancy
  int maxActiveBlocks;
  cudaOccupancyMaxActiveBlocksPerMultiprocessor( &maxActiveBlocks,
                                                 MyKernel, blockSize,
                                                 0);

  int device;
  cudaDeviceProp props;
  cudaGetDevice(&device);
  cudaGetDeviceProperties(&props, device);

  float occupancy = (maxActiveBlocks * blockSize / props.warpSize) /
                    (float)(props.maxThreadsPerMultiProcessor /
                            props.warpSize);

  printf("Launched blocks of size %d. Theoretical occupancy: %f\n",
         blockSize, occupancy);
}


int main()
{
  unsigned int arr_size = 1024 * 220000;
  unsigned int *dev_a, *dev_b, *dev_c;    // device копии of a, b, c
  unsigned int size = arr_size* sizeof( int );
  unsigned int* a = new unsigned int[arr_size];
  for(int i=0; i< arr_size; i++)
  {
    a[i] = i;
  }

  //выделяем память для device копий для a, b, c
  cudaMalloc( (void**)&dev_a, size );
  cudaMemcpy(dev_a, a, size, cudaMemcpyHostToDevice);
  launchMyKernel(dev_a, arr_size);
  cudaFree(dev_a);
  return 0;
}
