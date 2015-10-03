#include "apsp.cuh"

bool InitMMTimer(UINT wTimerRes){
	TIMECAPS tc;
	if (timeGetDevCaps(&tc, sizeof(TIMECAPS)) != TIMERR_NOERROR) {return false;}
	wTimerRes = min(max(tc.wPeriodMin, 1), tc.wPeriodMax);
	timeBeginPeriod(wTimerRes);
	return true;
}

void DestroyMMTimer(UINT wTimerRes, bool init){
	if(init)
		timeEndPeriod(wTimerRes);
}

void _CPU_Floyd(int *G,int *Gpath,int N){//standard N^3 algo
	for(int k=0;k<N;++k)for(int i=0;i<N;++i)for(int j=0;j<N;++j){
		int curloc=i*N+j,loca=i*N+k,locb=k*N+j;
		if(G[curloc]>(G[loca]+G[locb])){
			G[curloc]=(G[loca]+G[locb]);
			Gpath[curloc]=k;
		}
	}
}

void _showPath(int start,int end,const vector<Piii> &path,const int *D,const int N){
	cout<<"\nHere is the shortest cost path from "<<start<< " to "<<end<<", at a total cost of "<<D[start*N+end]<<".\n";
	for(int i=path.size()-1;i>=0;--i){
		cout<<"From "<<path[i].first.first<<" to "<<path[i].first.second<<" at a cost of "<<path[i].second<<'\n';
	}
	cout<<'\n';
}

bool _getPath(int curEdge, int nxtEdge,vector<Piii> &path,const int *D, const int *Dpath,const int N){
	int curIdx=curEdge*N+nxtEdge;
	if(D[curIdx]>=INF)return false;
	if(Dpath[curIdx]==-1){//end of backwards retracement
		path.push_back(make_pair(make_pair(curEdge,nxtEdge),D[curIdx]));
		return true;
	}else{//record las t edge cost and move backwards
		path.push_back(make_pair(make_pair(Dpath[curIdx],nxtEdge),D[Dpath[curIdx]*N+nxtEdge]));
		return _getPath(curEdge,Dpath[curIdx],path,D,Dpath,N);
	}
}

void _get_full_paths(const int *D, const int *Dpath, const int N){
	int start_vertex=-1,end_vertex=-1;
	vector<Piii> path;
	do{
		path.clear();
		cout<<"Enter start vertex #:";
		cin>>start_vertex;
		cout<<"Enter dest vertex(enter negative number to exit) #:";
		cin>>end_vertex;
		if(start_vertex<0 || start_vertex>=N || end_vertex<0 || end_vertex>=N)return;

		if(_getPath(start_vertex, end_vertex,path,D,Dpath,N)){
			_showPath(start_vertex,end_vertex,path,D,N);

		}else{
			cout<<"\nThere does not exist valid a path between "<<start_vertex<<" , and "<<end_vertex<<'\n';

		}
	}while(1);
}

__global__ void _Wake_GPU(int reps){
	int idx=blockIdx.x*blockDim.x + threadIdx.x;
	if(idx>=reps)return;
}

__global__ void _GPU_Floyd_kernel(int k, int *G,int *P, int N){//G will be the adjacency matrix, P will be path matrix
	int col=blockIdx.x*blockDim.x + threadIdx.x;
	if(col>=N)return;
	int idx=N*blockIdx.y+col;

	__shared__ int best;
	if(threadIdx.x==0)
		best=G[N*blockIdx.y+k];
	__syncthreads();
	if(best==INF || best > 10)return;
	int tmp_b=G[k*N+col];
	if(tmp_b==INF || tmp_b > 10)return;
//	if (cur > 1)
//		return;
	int cur = best + tmp_b;
	if(cur<G[idx]){
		G[idx]=cur;
		P[idx]=k;
	}
}
void _GPU_Floyd(int *H_G, int *H_Gpath, const int N){
	//allocate device memory and copy graph data from host
	int *dG,*dP;
	int numBytes=N*N*sizeof(int);
	cudaError_t err=cudaMalloc((int **)&dG,numBytes);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}
	err=cudaMalloc((int **)&dP,numBytes);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}
	//copy from host to device graph info
	err=cudaMemcpy(dG,H_G,numBytes,_HTD);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}
	err=cudaMemcpy(dP,H_Gpath,numBytes,_HTD);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}

	dim3 dimGrid((N+BLOCK_SIZE-1)/BLOCK_SIZE,N);

	for(int k=0;k<N;k++){//main loop

		_GPU_Floyd_kernel<<<dimGrid,BLOCK_SIZE>>>(k,dG,dP,N);
		err = cudaThreadSynchronize();
		if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}
	}
	//copy back memory
	err=cudaMemcpy(H_G,dG,numBytes,_DTH);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}
	err=cudaMemcpy(H_Gpath,dP,numBytes,_DTH);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}

	//free device memory
	err=cudaFree(dG);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}
	err=cudaFree(dP);
	if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);}
}

void _generateRandomGraph(int *G,int N,int range, int density){//density will be between 0 and 100, indication the % of number of directed edges in graph
	//range will be the range of edge weighting of directed edges
	int Prange=(100/density);
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			if(i==j){//set G[i][i]=0
				G[i*N+j]=0;
				continue;
			}
			int pr=rand()%Prange;
			G[i*N+j]= pr==0 ? ((rand()%range)+1):INF;//set edge random edge weight to random value, or to INF
		}
	}
}

int _read_from_file(int *G,const int N){//reads in edge list from file
	int num_edges=0;

	ifstream readfile;//enable stream for reading file
	readfile.open("Wiki-Vote.txt");
	cout << "Opening sucess \n";
	assert(readfile.good());//make sure it finds the file & file is
	string line;
	int v0,v1;
	while(getline(readfile,line)){
		istringstream linestream(line);
		linestream>>v0>>v1;
		cout << v0 <<"  "<< v1 <<  " readed ";
		G[v0*N+v1]=1;
		G[v1*N + v0] = 1;
		num_edges++;
	}
	readfile.close();
	cout << "File closing \n";
	return num_edges;
}

void _generate_result_file(bool success, unsigned int cpu_time,unsigned int gpu_time, int N){

	if(!success){
		cout<<"Error in calculation!\n";
		return;
	}else{
		ofstream myfile;
		myfile.open("Floyd-Warshall_result.txt");
		myfile<<"Success! The GPU Floyd-Warshall result and the CPU Floyd-Warshall results are identical(both final adjacency matrix and path matrix).\n\n";
		myfile<<"N= "<<N<<" , and the total number of elements(for Adjacency Matrix and Path Matrix) was "<<N*N<<" .\n";
		myfile<<"Matrices are int full dense format(row major) with a minimum of "<<(N*N)/4<<" valid directed edges.\n\n";
		myfile<<"The CPU timing for all was "<<float(cpu_time)/1000.0f<<" seconds, and the GPU timing(including all device memory operations(allocations,copies etc) ) for all was "<<float(gpu_time)/1000.0f<<" seconds.\n";
		myfile<<"The GPU result was "<<float(cpu_time)/float(gpu_time)<<" faster than the CPU version.\n";
		myfile.close();
	}
}


int main(){
	char ch;
	srand(time(NULL));

	const int N = 2 * 500;
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


	cout << "\nFloyd-Warshall on GPU underway:\n";
	_Wake_GPU << <1, BLOCK_SIZE >> >(32);

	//call host function which will copy all info to device and run CUDA kernels
	UINT wTimerRes = 0;
	bool init = InitMMTimer(wTimerRes);
	DWORD startTime = timeGetTime();

	_GPU_Floyd(H_G, H_Gpath, N);

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
	_CrtDumpMemoryLeaks();
	cin >> ch;
	return 0;
}
