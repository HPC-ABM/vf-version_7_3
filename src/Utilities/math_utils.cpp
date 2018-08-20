#include "math_utils.h"

int util::threadSafeRand(unsigned *seed_arr)
{
  int tid = 0;
#ifdef _OMP
  // Get thread id in order to access the seed that belongs to this thread
  tid = omp_get_thread_num();
#endif

  return rand_r(&(seed_arr[tid]));

}

float util::randFloatRange(float a, float b, unsigned *seed_arr) {
  float random = ((float) threadSafeRand(seed_arr)) / (float) RAND_MAX;
  float diff = b - a;
  float r = random * diff;
  return a + r;
}



bool util::divisible(long double a, long double b)
{
    bool result;
#if 1
    if(fabsl(((roundl(a/b)*b)- a)) <= (1E-9*b) ) {
        result = true;
    } else {
        result = false;
    }
#else
    if( fabsl(remainderl(a,b)) <= (1E-9*b ) ){
        result = true;
    } else {
        result = false;
    }
#endif
    // printf("divisible(%Lg, %Lg): %Lg, %Lg,%d\n", a, b, roundl(a/b), fabsl(((roundl(a/b)*b)-a)), result);
    return(result);
}

bool util::divisible(float a, float b)
{
    bool result;
#if 1
    if(fabs(((roundf(a/b)*b)- a)) <= (1E-6*b) ) {
        result = true;
    } else {
        result = false;
    }
#else
    if( fabs(remainderl(a,b)) <= (1E-6*b ) ){
        result = true;
    } else {
        result = false;
    }
#endif
    // printf("divisible(%Lg, %Lg): %Lg, %Lg,%d\n", a, b, roundl(a/b), fabsl(((roundl(a/b)*b)-a)), result);
    return(result);
}
