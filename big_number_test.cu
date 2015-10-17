#include <iostream>
#include <vector_types.h>
#include <stdlib.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/equal.h>

#include <time.h>
#include <sys/time.h>
#define USECPSEC 1000000ULL

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


int main(){
  unsigned long long t1_time, t2_time;
  float4 *mydata = new float4[DSIZE];
  for (int i = 0; i < DSIZE; i++){
    mydata[i].x = i;
    mydata[i].y = i;
    mydata[i].z = i;
    mydata[i].w = rand()/(float)RAND_MAX;}

  thrust::host_vector<float4>   h_data(mydata, mydata+DSIZE);
  // do once as a warm-up run, then report timings on second run
  for (int i = 0; i < 2; i++){
    thrust::device_vector<float4> d_data1 = h_data;
    thrust::device_vector<float4> d_data2 = h_data;

  // first time sort using typical thrust approach
    t1_time = dtime_usec(0);
    thrust::sort(d_data1.begin(), d_data1.end(), sort_f4_w());
    cudaDeviceSynchronize();
    t1_time = dtime_usec(t1_time);
  // now extract keys and create index values, sort, then rearrange
    t2_time = dtime_usec(0);
    thrust::device_vector<float> keys(DSIZE);
    thrust::device_vector<int> vals(DSIZE);
    thrust::copy(thrust::make_transform_iterator(d_data2.begin(), f4_to_fw()), thrust::make_transform_iterator(d_data2.end(), f4_to_fw()), keys.begin());
    thrust::sequence(vals.begin(), vals.end());
    thrust::sort_by_key(keys.begin(), keys.end(), vals.begin());
    thrust::device_vector<float4> result(DSIZE);
    thrust::copy(thrust::make_permutation_iterator(d_data2.begin(), vals.begin()), thrust::make_permutation_iterator(d_data2.begin(), vals.end()), result.begin());
    cudaDeviceSynchronize();
    t2_time = dtime_usec(t2_time);
    if (!validate(d_data1, result)){
      std::cout << "Validation failure " << std::endl;
      }
    }
  std::cout << "thrust t1 time: " << t1_time/(float)USECPSEC << "s, t2 time: " << t2_time/(float)USECPSEC << std::endl;
}
