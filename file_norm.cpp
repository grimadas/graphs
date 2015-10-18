#include<fstream>
#include<iostram>


int main(int argc, char* argv[])
{

	string a = std::atoi(argv[1]);
	std::cout <<" Reading file "<< a << std::endl;
	std::ifstream myfile;
	myfile.open(a);
	std::vector<std::pair<int, int>> items;
	int edges =0;
	int vertex = 0;
	int a, b;
	while (myfile >> a >> b)
	{
		items.
		number_of_edges++;
	}
	from_array_host = new vertex[number_of_edges];
	to_array_host = new vertex[number_of_edges];
	number_of_edges = 0;
	for(int i =0; i< number_of_vertex; i++)
	{
		for(int j = i+1; j< number_of_vertex; j++)
		{
			if (from_to_host_matrix[i*number_of_vertex + j] == 1)
			{
				from_array_host[number_of_edges] = i;
				to_array_host[number_of_edges] = j;
				number_of_edges++;
			}
		}
	}
	// Reading from file
	delete from_to_host_matrix;
	printf("Graph parametrs %d %d \n", number_of_vertex, number_of_edges);
	myfile.close();

	return 0;
}
