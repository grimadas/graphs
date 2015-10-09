#include "headers.h"

/****************************
*	 Converter functor,
*	 INPUT:
*				_a - from_array
*				_b - to_array
*				_size - number_of_edges
*/
struct coo_to_csr_converter
{
  __host__ __device__
  coo_to_csr_converter(domain _a, domain _b, int _size) : a(_a), b(_b), size(_size){}

  __host__ __device__
    field operator()(field x)
  {

      if (x < size)
      {
        return b[x];
      }
      else
      {
        return a[x - size];
      }


    }

  domain a;
  domain b;
  int size;
};
