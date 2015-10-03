//The following code implemented by Oleg Konings in association with Morgan Hough and Gazzaley lab
//A simple implementation of the Floyd-Warshall all-pairs-shortest path algorithm with path reconstruction.
//This is indended to be used on directed graphs with no negative cycles
//The Adjacency Matrix is in Row-major format, and is implemented both in CUDA on a Nvidia GTX 680 2GB GPU,
// and in serial CPU code using an Intel i7-3770 3.9 ghz.
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <cstdio>
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
#include <crtdbg.h>//to detect host memory leaks
using namespace std;

#define _DTH cudaMemcpyDeviceToHost
#define _HTD cudaMemcpyHostToDevice

//these can be altered on user depending on data set and type of operation(random test, read from file etc)
#define BLOCK_SIZE 256
#define RANGE 1
#define RANDOM_GSIZE 2 * 500
#define FILE_GSIZE 8298//the number of edges in Wiki-Vote.txt if the file test is run
#define INF (1<<22)
#define DO_TEST_RANDOM 0
#define DO_TEST_FROM_FILE 1
#define DENSITY 25

//typedef for vector used in path reconstruction
typedef pair<pair<int,int>,int> Piii;

//forward function declarations
bool InitMMTimer(UINT wTimerRes);
void DestroyMMTimer(UINT wTimerRes, bool init);
void _CPU_Floyd(int *G,int *Gpath,int N);
void _showPath(int start,int end,const vector<Piii> &path,const int *D,const int N);
bool _getPath(int curEdge, int nxtEdge,vector<Piii> &path,const int *D, const int *Dpath,const int N);
void _get_full_paths(const int *D, const int *Dpath, const int N);

//CUDA GPU kernel/functions forward declaration
__global__ void _Wake_GPU(int reps);
__global__ void _GPU_Floyd_kernel(int k, int *G,int *P, int N);
void _GPU_Floyd(int *H_G, int *H_Gpath, const int N);

//other optional utility functions
int _read_from_file(int *G,const int N);
void _generateRandomGraph(int *G, int N, int range, int density);
void _generate_result_file(bool success, unsigned int cpu_time, unsigned int gpu_time, int N);
