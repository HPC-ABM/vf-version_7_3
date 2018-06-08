#ifndef PRINT_UTILS_H_
#define PRINT_UTILS_H_

#pragma once

#include <string>
#include "../common.h"

using namespace std;
namespace util {


void printWindow(float* a, int h, int w, int r){
#ifdef PRINT_KERNEL
	int mx = w/2;
	int my = h/2;
	int i_start	= my - r;
	int i_end	= my + r + 1;
	int j_start	= mx - r;
	int j_end	= mx + r + 1;
	if (i_start < 0 || i_end >= h || j_start < 0 || j_end >= w) {
		cout << "printWindowError: " << h << "x" << w << endl;
		cout << "	i_start:	" << i_start << endl;
		cout << "	i_end:		" << i_end << endl;
		cout << "	j_start:	" << j_start << endl;
		cout << "	j_end:		" << j_end << endl;
		exit(-1);
	}
	cout << "Matrix window: ---" << endl;
	for (int i = i_start; i < i_end; i++) {
		for (int j = j_start; j < j_end; j++) {
			cout << a[i*w + j] << " ";
		}
		cout << endl;
	}
#endif
}

void displayWindow(float *arr, int w, int h, int winsize)
{
    int halfw = w/2;
    int halfh = h/2;

    int halfs = winsize /2;

    int strt_i = halfw - halfs;
    int end_i  = halfw + halfs;
    int strt_j = halfh - halfs;
    int end_j  = halfh + halfs;

    if (winsize > w || winsize > h) {
        printf("Print window dimension errors\n");
        return;
    }

    printf("Window\n");

    for (int i = strt_i; i < end_i; i++)
    {
        for (int j = strt_j; j < end_j; j++)
        {
            printf("%f\t", arr[i*w + j]);
        }
        printf("\n");
    }

}


void displayWindowPlane(float *arr, int w, int h, int z, int winsize)
{
    int halfw = w/2;
    int halfh = h/2;

    int halfs = winsize /2;

    int strt_i = halfw - halfs;
    int end_i  = halfw + halfs;
    int strt_j = halfh - halfs;
    int end_j  = halfh + halfs;

    if (winsize > w || winsize > h) {
        printf("Print window dimension errors\n");
        return;
    }

    printf("Window Plane %d\n", z);

    for (int i = strt_i; i < end_i; i++)
    {
        for (int j = strt_j; j < end_j; j++)
        {
            printf("%f\t", arr[z*h*w + i*w + j]);
        }
        printf("\n");
    }

}

}  // namespace util

#endif /* PRINT_UTILS_H_ */
