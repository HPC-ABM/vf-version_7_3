
#ifndef MEM_UTILS_H_
#define MEM_UTILS_H_

#pragma once

#include <string>
#include <iostream>

//#include "error_utils.h"

using namespace std;
namespace util {

template<typename T>
void allocate(T** ptr, int n)
{
    cout << quote(ptr) << ": requesting for " <<
            (sizeof(T)*n)/(1024.0*1024.0*1024.0) << " GB\n" << endl;
    try {
        *ptr = new T [n];

//      if (util::ABMerror(
//              !((*ptr) = new T [n]),
//              "Mem alloc error!",
//              __FILE__,
//              __LINE__))
//          exit(1);

    } catch (std::bad_alloc& ba) {
        printf("Mem allocation error: %s\n", ba.what());
        exit(-1);
    }
}


}  // namespace util

#endif /* MEM_UTILS_H_ */

