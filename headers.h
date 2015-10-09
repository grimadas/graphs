/*
*   Library headers
*   @author = Bulat Nasrulin
*/


// Thrust includes
#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/find.h>
#include <thrust/count.h>
#include <thrust/reduce.h>
#include <thrust/merge.h>
#include <thrust/sequence.h>
#include <thrust\sort.h>
#include <thrust/unique.h>
#include <thrust/execution_policy.h>
#include <thrust\iterator\counting_iterator.h>
#include <thrust\iterator\permutation_iterator.h>
#include <thrust\binary_search.h>
// Headers for floyd warshall algorithm
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <cstring>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <cuda.h>
#include <ctime>
#include <cassert>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define pb push_back
#define all(c) (c).begin(),(c).end()
#include <Windows.h>
#include <MMSystem.h>

#pragma comment(lib, "winmm.lib")
#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h> //to detect host memory leaks

#define _DTH cudaMemcpyDeviceToHost
#define _HTD cudaMemcpyHostToDevice

//typedef for vector used in path reconstruction
typedef pair<pair<int,int>,int> Piii;

using namespace std;


// Predefined constants
//these can be altered on user depending
// on data set and type of
// operation(random test, read from file etc)
#define BLOCK_SIZE 256
#define RANGE 1
#define RANDOM_GSIZE 2 * 500
#define FILE_GSIZE 8298 //the number of edges in Wiki-Vote.txt if the file test is run
#define INF (1<<22)
#define DO_TEST_RANDOM 0
#define DO_TEST_FROM_FILE 1
#define DENSITY 25

/* Defined fields for graph */
#define vertex  unsigned int
#define edge  unsigned int

#define domain vertex*
#define field  vertex

#define opacity double
