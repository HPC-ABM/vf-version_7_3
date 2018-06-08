#include "../World/Usr_World/woundHealingWorld.h"

#include <string>

class World;
class WHWorld;

using namespace std;
namespace util {


bool divisible(long double a, long double b)
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

bool divisible(float a, float b)
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



}
