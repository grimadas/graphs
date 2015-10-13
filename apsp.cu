#include "apsp.cuh"

int main(){
	char ch;
	srand(time(NULL));

	const int N =  9;
	const int NumBytes = N*N*sizeof(int);
	//host allocations to create Adjancency matrix and result matrices with path matrices
	int *OrigGraph = (int *)malloc(NumBytes);//will be original Adjancency matrix, will NOT be changed
	int *H_G = (int *)malloc(NumBytes);
	int *H_Gpath = (int *)malloc(NumBytes);
	int *D_G = (int *)malloc(NumBytes);
	int *D_Gpath = (int *)malloc(NumBytes);


	for (int i = 0; i<N*N; i++){//copy for use in computation
		//	H_G[i] = D_G[i] = OrigGraph[i];//copy for use in computation
		H_Gpath[i] = D_Gpath[i] = -1;//set to all negative ones for use in path construction
		H_G[i] = INF;
	}

	//Read file from file
	cout << "Reading file with graph: \n ";
	_read_from_file(H_G, N);
	cout << "Assigning initial paths ";


	cout << "Graph " << endl;
	for (int i = 0; i<N; i++){//copy for use in computation
		for (int j = 0; j < N; j++)
		{
			if (i == j)
				H_G[i*N + j] = 0;
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

	int L = 10;
	_GPU_Floyd(H_G, H_Gpath, N, L);

	unsigned int endTime = timeGetTime();
	unsigned int gpu_time = unsigned int(endTime - startTime);
	printf("GPU Timing(including all device-host, host-device copies, device allocations and freeing of device memory): %dms\n\n", gpu_time);
	DestroyMMTimer(wTimerRes, init);

	cout << "H_PATH" << endl;
	for (int i = 0; i<N; i++){//copy for use in computation
		for (int j = 0; j < N; j++)
		{
			//	cout << H_Gpath[i*N + j] << "  ";
		}
		//cout << endl;
	}

	cout << "H_G" << endl;
	for (int i = 0; i<N; i++){//copy for use in computation
		for (int j = 0; j < N; j++)
		{
			//		cout << H_G[i*N + j] << "  ";
		}
		//	cout << endl;
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
