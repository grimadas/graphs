/*
 	Main application
	Author : Bulat, 2015

*/
#include "graph.h"
#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>



__global__ void CUDA_APSP()
{

}


__global__ void vector_sum(domain a, domain b, domain result, int size)
{

}

void simple_sum(domain a, domain b, domain result, int size)
{

}

void simple_zero_fill(domain a, int size)
{
	for (int i = 0; i < size; ++i)
	{
		a[i] = 0;
	}
}

void simple_fill(domain& a, int size)
{
	for (int i = 0; i < size; i++)
	{
		a[i] = 3;
	//	printf("Index i %f ", (field) i);
	}
}

__global__ void fill(domain a, int size)
{
//	unsigned int index;
//	index = blockIdx.x *blockDim.x + threadIdx.x;
//	c[index] = 1;
	//printf("The index is %d \n", index);
//	printf("%d \n", a[index]);
	//a[index] = 1;

	//printf("%d \n", *a);
	//printf("%d \n", a[index]);
}
__device__ unsigned int index2;



__global__ void my_sum(domain a, domain b, domain c)
	{
		unsigned int index;

	 index = blockIdx.x *blockDim.x + threadIdx.x;

	//	cout >> blockIdx.x >>" " >> blockDim.x >>" " >> threadIdx.x >> endl;
		//cout << "Index is " <<index;
	//	printf("%d\n",index );

		{
			c[index] = a[index] + b[index];
			printf(" %d  + %d = %d \n", a[index], b[index], c[index]);
		}
	}

 __global__ void cuda_change(domain a, int size)
	{

		int	index = blockIdx.x *blockDim.x + threadIdx.x;
		a[index] = (field)5;
	}




__global__ void cuda_print(domain a, int size)
{

	int	index = blockIdx.x *blockDim.x + threadIdx.x;
	printf("%d  = %u ",index, a[index]);
}

void print_vector(domain a,int size)
{
	for (int i = 0; i < size; ++i)
	{

		printf("%u ", a[i]);
	}
	printf("\n");
}

__shared__ domain  a;



void thrust_test()
{
	using namespace thrust;
	int major = THRUST_MAJOR_VERSION;
	int minor = THRUST_MINOR_VERSION;

	std::cout << "Thrust v" << major << "." << minor << std::endl;

	// H has storage for 4 integers
	thrust::host_vector<int> H(4);

	// initialize individual elements
	H[0] = 14;
	H[1] = 20;
	H[2] = 38;
	H[3] = 46;

	// H.size() returns the size of vector H
	std::cout << "H has size " << H.size() << std::endl;

	// print contents of H
	for(int i = 0; i < H.size(); i++)
	{
		std::cout << "H[" << i << "] = " << H[i] << std::endl;
	}

	// resize H
	H.resize(2);

	std::cout << "H now has size " << H.size() << std::endl;

	// Copy host_vector H to device_vector D
	device_vector<int> D = H;

	// elements of D can be modified
	// D[0] = 99;
    // 	D[1] = 88;

	// print contents of D
	for(int i = 0; i < D.size(); i++)
	{
		std::cout << "D[" << i << "] = " << D[i] << std::endl;
	}

	// H and D are automatically destroyed when the function returns



}


int main()
{
	Graph test_graph;
	
	test_graph.init_test_graph();
	test_graph.print_coo_graph();
	test_graph.test_func();
	test_graph.print_coo_graph();



	return 0;
}
