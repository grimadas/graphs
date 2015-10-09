
#include "headers.h"
__global__ void printValue( int *value) {
printf("value %d",value[0]);
}

void hostFunction(int *value){

value[0]=1;
value[1]=2;
printValue<<< 1, 1 >>>(value);
cudaDeviceSynchronize();
cudaFree(value);
}

int main() {
int *value;
cudaMallocManaged(&value, 2 * sizeof(int));
hostFunction(value);
return 0;
}
