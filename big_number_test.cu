#include <iostream>
#include <vector_types.h>
#include <stdlib.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/equal.h>
#include <thrust/fill.h>
#include <thrust/merge.h>
#include <algorithm>
#include <fstream>
#include <stdio.h>


#include <time.h>
#include <sys/time.h>
#define USECPSEC 1000000ULL

int number_of_edges;
int number_of_vertex;

using namespace thrust;
unsigned long long dtime_usec(unsigned long long start){

  timeval tv;
  gettimeofday(&tv, 0);
  return ((tv.tv_sec*USECPSEC)+tv.tv_usec)-start;
}

#define DSIZE (32*1048576)

struct sort_f4_w
{
  __host__ __device__
  bool operator()(const float4 &a, const float4 &b) const {
    return (a.w < b.w);}
};
// functor to extract the .w element from a float4
struct f4_to_fw : public thrust::unary_function<float4, float>
{
  __host__ __device__
  float operator()(const float4 &a) const {
    return a.w;}
};
// functor to extract the .x element from a float4
struct f4_to_fx : public thrust::unary_function<float4, float>
{
  __host__ __device__
  float operator()(const float4 &a) const {
    return a.x;}
};


bool validate(thrust::device_vector<float4> &d1, thrust::device_vector<float4> &d2){
  return thrust::equal(thrust::make_transform_iterator(d1.begin(), f4_to_fx()), thrust::make_transform_iterator(d1.end(), f4_to_fx()), thrust::make_transform_iterator(d2.begin(), f4_to_fx()));
}


void inline read_COO_format(char* file_name, int*& from_array_host, int*& to_array_host)
{


}




int main(){

  unsigned long long t1_time, t2_time;
    t1_time = dtime_usec(0);
    std::ifstream myfile;
    printf("Reading from file");
    myfile.open("graph.txt");
    myfile >> number_of_vertex >> number_of_edges;
  //	number_of_edges = 0;
    // reserve maximum value needed
    int* ar = new int[number_of_edges];
    int* br = new int[number_of_edges];
    int a,b;
    int i =0;
    while(myfile >> a >> b)
    {
      ar[i] = a;
      br[i] = b;
      i++;
    }
    number_of_edges = i-1;
    printf("Graph parametrs %d %d \n", number_of_vertex, number_of_edges);
    myfile.close();
    t1_time = dtime_usec(t1_time);
  // do once as a warm-up run, then report timings on second run

  device_ptr<int> d_from = device_malloc<int>(number_of_edges);

  device_ptr<int> d_to = device_malloc<int>(number_of_edges);

  copy(ar, ar+number_of_edges, d_from);
  copy(br, br+number_of_edges, d_to);
  device_ptr<int> temp_indx = device_malloc<int>(2*number_of_edges);
  //device_vector<int> temp_indx(20*number_of_edges);
  //wrap raw pointer with a device_ptr to use with Thrust functions
  //device_vector<int> rost(10*2*number_of_edges, -1);
  device_ptr<int> rost = device_malloc<int>(2*number_of_edges);
  counting_iterator<int> index_from(0);
  counting_iterator<int> index_to(number_of_edges);

  //	Merging from and to arrays are keys,
  //	indexes are (0..number_edges) and (num_edges to 2*number_edges)



  merge_by_key(device,
    d_from, d_from + number_of_edges,
    d_to, d_to + number_of_edges,
    index_from, index_to,
    temp_indx,
    rost);
    cudaDeviceSynchronize();
    thrust::copy(seq, rost, rost + number_of_edges, std::ostream_iterator<int>(std::cout, ","));
   std::cout << std::endl;

std::cout << "Merge ok";

/*    cudaDeviceSynchronize();
    sort_by_key(device,
      temp_indx, temp_indx + number_of_edges,
      rost);
      */
      std::cout << "Sort ok";

  // first time sort using typical thrust approach

  //  thrust::sort(d_data1.begin(), d_data1.end(), sort_f4_w());
//    cudaDeviceSynchronize();

  // now extract keys and create index values, sort, then rearrange
    t2_time = dtime_usec(0);

    t2_time = dtime_usec(t2_time);

  std::cout << "thrust t1 time: " << t1_time/(float)USECPSEC << "s, t2 time: " << t2_time/(float)USECPSEC << std::endl;
}
