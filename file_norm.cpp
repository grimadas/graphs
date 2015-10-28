#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>
#include <algorithm>

struct value_comparer {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.second < right.second;
    }
};

struct key_comparer {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.first < right.first;
    }
};

int main(int argc, char* argv[])
{

    char* file_name = argv[1];
    std::cout <<" Reading file "<< file_name << std::endl;
    std::ifstream myfile;
    myfile.open(file_name);
    std::vector<std::pair<int, int> > items;
    int edges =0;
    int vertex = 0;
    int a,b;
    while (myfile >> a >> b)
    {
        if (a > vertex)
        {
            vertex = a;
        }
        if (b > vertex)
        {
            vertex = b;
        }
        items.push_back(std::make_pair(std::min(a,b),std::max(a,b)));
        edges++;
    }
    myfile.close();


    std::sort(items.begin(), items.end(), value_comparer());

    //sort by key using std::stable_sort
    std::stable_sort(items.begin(), items.end(), key_comparer());



    std::ofstream out_put_file;
    char* out_name = argv[2];
    out_put_file.open(out_name);
    out_put_file << vertex+1 << " " << edges << "\n";
    for (int i=0; i< edges; i++)
    {
        out_put_file << items[i].first << " " << items[i].second << "\n";

    }
    out_put_file.close();
    std::cout << "Done.";
    return 0;
}
