/*
 	Main application
	Author : Bulat, 2015

*/
#include "graph.h"
#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

int main()
{
	init_test_graph();
	print_test();
	test_func();
	return 0;
}
