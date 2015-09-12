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
			printf("index %d = %u  \n", index, c[index]);
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
//  const int size = 500;
/*	domain a;
	domain b;
	domain c;
 	 a = (domain)malloc(size*sizeof(field));
 	 b = (domain)malloc(size*sizeof(field));
 	 c = (domain)malloc(size*sizeof(field));
*/
	Graph test_graph;
	//
	int vertex_number = 50;
	int edges_per_vertex = 4;

	thrust::device_vector<vertex> sum(test_graph.from_array.size() + test_graph.to_array.size());
	domain f1 = 
	test_graph.random_coo_graph(vertex_number, edges_per_vertex);
	test_graph.print_coo_graph();
	test_graph.calc_degree();
	my_sum<<<100, 100>>>(thrust::raw_pointer_cast(test_graph.from_array),
		thrust::raw_pointer_cast(test_graph.to_array), 
		thrust::raw_pointer_cast(sum));


	//test_graph.random_coo_graph_expr(vertex_number, edges_per_vertex);
	//cout << "Must print ";
	//test_graph.print_coo_graph_expr();
	//cout << endl <<"Density of graph " <<test_graph.get_density();


  // --------------------------------

/*	domain cu_a;
	domain cu_b;
	domain cu_c;

	cudaMalloc((void**)&cu_a, size*sizeof(field));
	cudaMalloc((void**)&cu_b, size*sizeof(field));
	cudaMalloc((void**)&cu_c, size*sizeof(field)); */


	// -------------------------------
	/* Filling with zeros */
/* simple_fill(a, size);
	printf("A: \n");
	print_vector(a, size);
//	printf("Vector a: \n");
//	print_vector(a, size);
	simple_fill(b, size);
//	simple_fill(c, size);
	// ------------------------------
	cudaMemcpy(cu_a, a, size * sizeof(field), cudaMemcpyHostToDevice);
	cudaMemcpy(cu_b, b, size * sizeof(field), cudaMemcpyHostToDevice);

	dim3 grid_size(1);
  dim3 block_size(size);
	cuda_change<<<grid_size, block_size>>>(cu_a, size);
	cuda_change<<<grid_size, block_size>>>(cu_b, size);

	//cuda_print<<<grid_size, block_size>>>(cu_a, size);
	// main function
	clock_t begin = clock();
  my_sum<<<grid_size, block_size>>>(cu_a,
		 cu_b, cu_c, size);
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

	// Place memory back to CPU
	printf("C array \n");
  //cuda_print<<<grid_size, block_size>>>(cu_c, size);
 	cudaMemcpy(c, cu_c, size * sizeof(field), cudaMemcpyDeviceToHost);


	printf("%f \n", elapsed_secs);
	//print_vector(a, size);
//	print_vector(b, size);
  print_vector(c, size); */

//	thrust_test();

  return 0;
}
