/*
 * timer.h
 *
 * File Contents: Contains declarations for functions for timing code sections
 *
 * Created on: May 26, 2015
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef TIMER_H_
#define TIMER_H_

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>

#include <unistd.h>
#include <sys/time.h>

// Helper functions for CUDA
#include <helper_functions.h>
#include <helper_cuda.h>

using namespace std;

/*
 * Description:	Function to start timer
 *
 * Return: Time at the start
 *
 * Parameters: void
 */
struct timeval tic();

/*
 * Description:	Function to stop timer
 *
 * Return: void
 *
 * Parameters: timeval  -- Time at the start
 */
long toc(struct timeval);

/*
 * Description: Function to start GPU timer
 *
 * Return: Initialized and started timer
 *
 * Parameters: void
 */
StopWatchInterface* gputic();

/*
 * Description: Function to stop GPU timer
 *
 * Return: Total time
 *
 * Parameters: hTimer  -- Initialized and started timer
 */
double gputoc(StopWatchInterface* hTimer);

/*
 * Description:	Output time info
 *
 * Return: void
 *
 * Parameters: time        -- Elapsed time in millisecond
 *             stageName   -- The name of the stage to be displayed
 *             stageLevel  -- The stage level to be displayed
 */
void print_time(
		long time,
		const char* stageName,
		const char* stageLevel);

/*
 * Description:	A macro for timing a command/function of a major stage in go() and print the timing info
 *
 * Return: void
 *
 * Parameters: task    -- Command/function to perform
 *             sName   -- The name of the stage to be displayed
 *             sLevel  -- The stage level to be displayed
 *
 * Example Usage:
 * 			TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");
 */
#define TIME_STAGE(task, sName, sLevel) do { \
	struct timeval START = tic ();\
	task;\
	long ET = toc(START);\
	print_time(ET, sName, sLevel);\
	} while (0)

/*
 * Description: A macro for timing a command/function of a major stage in go() and print the timing info
 *
 * Return: void
 *
 * Parameters: task    -- Command/function to perform
 *             tvar    -- Name of the variable to store the time (double)
 *
 * Example Usage:
 *          TIME_GPU_STAGE(H2D(ic, this->chem_cctx, np), transerTime);
 */
#define TIME_GPU_STAGE(task, tvar) do { \
    StopWatchInterface* TIMER = gputic ();\
    task;\
    double ET = gputoc(TIMER);\
    memcpy(&tvar, &ET);\
    } while (0)


#endif /* TIMER_H_ */
