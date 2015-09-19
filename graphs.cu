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

__shared__ __device__ int iter_koef;

__global__ void a_iter(domain a, domain b, domain c, int k)
{
	int index =
		blockIdx.x *blockDim.x +
		threadIdx.x;
	if (a[index] == k)
	{
		c[iter_koef] = b[iter_koef];
		iter_koef++;
		__syncthreads();

	}

}


__global__ void b_iter(domain a, domain b, domain c, int k)
{
	int index =
		blockIdx.x *blockDim.x +
		threadIdx.x;
	if (b[index] == k)
	{
		c[iter_koef] = a[iter_koef];
		iter_koef++;
		__syncthreads();

	}

}

__global__ void k_iter(domain a, domain b, domain c, int size)
{
	
//	a_iter<<<1, size>>>(a, b, c, index);
//	b_iter<<<1, size>>>(a, b, c, index);
}
/*
void combine_edges()
{

	domain a = thrust::raw_pointer_cast(from_array.data());
	domain b = thrust::raw_pointer_cast(to_array.data());
	full_edge_array.reserve(L_VALUE * 2 * number_of_edges);
	domain c = thrust::raw_pointer_cast(&(full_edge_array[0]));


	k_iter << <1, number_of_vertex >> >(a,b,c,number_of_edges);
}*/

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


__global__ void calc_a()
{

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




int main()
{
	init_test_graph();
	print_test();
	test_func();
	return 0;
}
