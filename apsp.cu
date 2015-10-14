#include "apsp.cuh"


int main(int argc, char* argv[]){
	char ch;
	srand(time(NULL));

	int vertex_number = std::atoi(argv[1]);
	int l_value = std::atoi(argv[2]);
	const int N =  vertex_number;
	int L = l_value;
	cout << "Vertex num " << N << " l value " << 3;
	const int NumBytes = N*N*sizeof(int);
	//host allocations to create Adjancency matrix and result matrices with path matrices
	int *OrigGraph = (int *)malloc(NumBytes);//will be original Adjancency matrix, will NOT be changed
	int *H_G = (int *)malloc(NumBytes);
	int *H_Gpath = (int *)malloc(NumBytes);
	int *D_G = (int *)malloc(NumBytes);
	int *D_Gpath = (int *)malloc(NumBytes);

  thrust::fill(H_Gpath, H_Gpath + N*N, -1);
	thrust::fill(D_Gpath, D_Gpath + N*N, -1);
	thrust::fill(H_G, H_G + N*N, INF);
/*
	for (int i = 0; i<N*N; i++){//copy for use in computation
		//	H_G[i] = D_G[i] = OrigGraph[i];//copy for use in computation
		H_Gpath[i] = D_Gpath[i] = -1;//set to all negative ones for use in path construction
		H_G[i] = INF;
	}
*/
	//Read file from file
	cout << "Reading file with graph: \n ";
	_read_from_file(H_G, N);
	cout << "Assigning initial paths ";


	cout << "Graph " << endl;
	for (int i = 0; i<N; i++){//copy for use in computation
		{
				H_G[i*N + i] = 0;
			//			cout << H_G[i*N + j] << "  ";

		}
		//	cout << endl;
	}

	UINT wTimerRes = 0;
	bool init = InitMMTimer(wTimerRes);
	DWORD startTime = timeGetTime();
	cout << "\nFloyd-Warshall on GPU underway:\n";
	_Wake_GPU << <1, BLOCK_SIZE >> >(32);

	//call host function which will copy all info to device and run CUDA kernels


	_GPU_Floyd(H_G, H_Gpath, N, L);

	unsigned int endTime = timeGetTime();
	unsigned int gpu_time = unsigned int(endTime - startTime);
	printf("GPU Timing(including all device-host, host-device copies, device allocations and freeing of device memory): %dms\n\n", gpu_time);
	DestroyMMTimer(wTimerRes, init);
	// The result is stored in H_G (the level)
	cout << "H_PATH" << endl;
	for (int i = 0; i<N; i++){//copy for use in computation
		for (int j = 0; j < N; j++)
		{
				cout << H_Gpath[i*N + j] << "  ";
		}
		cout << endl;
	}

	cout << "H_G" << endl;
	for (int i = 0; i<N; i++){//copy for use in computation
		for (int j = 0; j < N; j++)
		{
					cout << H_G[i*N + j] << "  ";
		}
			cout << endl;
	}


	free(OrigGraph);
	free(H_G);
	free(H_Gpath);
	free(D_G);
	free(D_Gpath);
	//}
	//	*/

	cin >> ch;
	return 0;
}
