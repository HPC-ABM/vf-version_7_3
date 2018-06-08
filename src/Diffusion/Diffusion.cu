/*
 * Diffusion.cpp
 *
 *  Created on: May 18, 2016
 *      Author: NungnunG
 */

#include "../common.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
// Include CUDA runtime and CUFFT
#include <cuda_runtime.h>
#include <cufft.h>

// Helper functions for CUDA
#include <helper_functions.h>
#include <helper_cuda.h>

#include "../Utilities/output_utils.h"
#include "convolutionFFT_common.h"
//#include "convolutionFFT2D.cu"

using namespace std;

#ifdef GPU_DIFFUSE      // (*)

bool firstTransferCompleted = false;


////////////////////////////////////////////////////////////////////////////////
//// Helper functions
//////////////////////////////////////////////////////////////////////////////////
int snapTransformSize(int dataSize)
{
	int hiBit;
	unsigned int lowPOT, hiPOT;

	dataSize = iAlignUp(dataSize, 16);

	for (hiBit = 31; hiBit >= 0; hiBit--)
		if (dataSize & (1U << hiBit))
		{
			break;
		}

	lowPOT = 1U << hiBit;

	if (lowPOT == (unsigned int)dataSize)
	{
		return dataSize;
	}

	hiPOT = 1U << (hiBit + 1);

	if (hiPOT <= 1024)
	{
		return hiPOT;
	}
	else
	{
		//return iAlignUp(dataSize, 512);
		return pow(2, ceil(log(dataSize)/log(2)));
	}
}

void H2D(int ic, c_ctx* chem_cctx, int np)
{

#ifdef MODEL_3D
	int ig = chem_cctx->gpu_id[ic];
	cudaSetDevice(chem_cctx->dev_id[ic]);//ig);

#ifndef M40
	const int       fftD = chem_cctx->FFTD;
	const int       fftH = chem_cctx->FFTH;
	const int       fftW = chem_cctx->FFTW;

	size_t fssize_b	= fftD * fftH * (fftW / 2 + 1) * sizeof(fComplex);
	//    size_t fsize_b	= fftD		* fftH		* fftW		* sizeof(float);

	// Copy Kernel Spectrum from Host to Device
	checkCudaErrors(cudaMemcpy(chem_cctx->d_kernelspectrum_h[ic], chem_cctx->h_kernelspectrum[ic],
			fssize_b, cudaMemcpyHostToDevice));
#endif	// M40
#endif	// MODEL_3D

	float **h_ibuffs = chem_cctx->h_ibuffs;
#ifdef ASYNCCPY
	printf("async P2D copy: gpu[%d] chem [%d]\n", ig, ic);
	checkCudaErrors(cudaMemcpyAsync(chem_cctx->d_data[ic], h_ibuffs[ic], np*sizeof(float),
			cudaMemcpyHostToDevice, stream[ig]));
#else	// ASYNCCPY
	checkCudaErrors(cudaMemcpy(chem_cctx->d_data[ic], h_ibuffs[ic],
			np*sizeof(float), cudaMemcpyHostToDevice));
#endif	// ASYNCCPY
}

void D2H(int ic, c_ctx* chem_cctx, int np)
{

#ifdef MODEL_3D
	int ig = chem_cctx->gpu_id[ic];
	cudaSetDevice(chem_cctx->dev_id[ic]);//ig);
#endif	// MODEL_3D

	float **h_obuffs = chem_cctx->h_obuffs;

#ifdef ASYNCCPY
	printf("async D2P copy: gpu[%d] chem [%d]\n", i, ic);
	checkCudaErrors(cudaMemcpyAsync(h_obuffs[ic], chem_cctx->d_data[ic], np*sizeof(float),
			cudaMemcpyDeviceToHost, stream[ig]));
#else	// ASYNCCPY
	checkCudaErrors(cudaMemcpy(h_obuffs[ic], chem_cctx->d_data[ic],
			np*sizeof(float), cudaMemcpyDeviceToHost));
#endif	// ASYNCCPY
}

bool compareResults(float *h_ResultCPU, float *h_ResultGPU,
		int dataW, int dataH, int dataD, float eThreshold)
{
	double sum_delta2 = 0;
	double sum_ref2   = 0;
	double max_delta_ref = 0;

	double sum = 0;
	for (int z = 0; z < dataD; z++)
		for (int y = 0; y < dataH; y++)
			for (int x = 0; x < dataW; x++)
			{
				double  rCPU = (double)h_ResultCPU[z * dataH * dataW + y * dataW + x];
				double  rGPU = (double)h_ResultGPU[z * dataH * dataW + y * dataW  + x];
				double delta = (rCPU - rGPU) * (rCPU - rGPU);
				double   ref = rCPU * rCPU + rCPU * rCPU;

				if ((delta / ref) > max_delta_ref)
				{
					max_delta_ref = delta / ref;
				}

				//                if ((ref-0.0) > 0.000001){
					//                    if ((delta / ref) > max_delta_ref)
						//                    {
				//                        max_delta_ref = delta / ref;
				//                    }
				//                }

				sum_delta2 += delta;
				sum_ref2   += ref;

				sum += rGPU;

			}

	double L2norm = sqrt(sum_delta2 / sum_ref2);
	printf("rel L2 = %E (max delta = %E)\n", L2norm, sqrt(max_delta_ref));
	bool bRetVal = (L2norm < eThreshold) ? true : false;
	printf(bRetVal ? "L2norm Error OK\n" : "L2norm Error too high!\n");

	return bRetVal;


}

#ifndef MODEL_3D

/*
bool computeChemDiffusionCPU(
                float           *h_ChemOut,
                float           *h_ChemIn,
                float           *h_Kernel,
                c_ctx            cctx,
                int                      iter)
{

        bool bRetVal = 1;
        StopWatchInterface *hTimer = NULL;
        sdkCreateTimer(&hTimer);

        printf("Testing Chemical Diffusion CPU\n");
        const int    kernelH = cctx.KH;//7;
        const int    kernelW = cctx.KW;//6;
        const int    kernelY = cctx.KY;//3;
        const int    kernelX = cctx.KX;//4;
        const int      dataH = cctx.DH;//100;//1160;//2000;
        const int      dataW = cctx.DW;//100;//1660;//2000;
        const int outKernelH = cctx.DH;
        const int outKernelW = cctx.DW;
        const int       fftH = cctx.FFTH;
        const int       fftW = cctx.FFTW;


        sdkResetTimer(&hTimer);
        sdkStartTimer(&hTimer);



        fprintf(stderr,"...running reference CPU convolution\n");
        convolutionClampToBorderCPU(
                        h_ChemOut,
                        h_ChemIn,
                        h_Kernel,
                        dataH,
                        dataW,
                        kernelH,
                        kernelW,
                        kernelY,
                        kernelX
        );

        sdkStopTimer(&hTimer);
        double cpuTime = sdkGetTimerValue(&hTimer);
        printf("\t\tCPU chemical diffusion computation:\t%f MPix/s (%f ms)\n",
                        (double)dataH * (double)dataW * 1e-6 / (cpuTime * 0.001), cpuTime);

        sdkDeleteTimer(&hTimer);
        return bRetVal;
}
 */


void padData2D(
		pad_t pt,
		float *d_PaddedData,
		float *d_Data,
		int fftH,
		int fftW,
		int dataH,
		int dataW,
		int kernelH,
		int kernelW,
		int kernelY,
		int kernelX,
		int epiBoundary,
		float  baseChem)
{
	switch(pt)
	{
	case pClampToBorder:
		padDataClampToBorder(
				d_PaddedData,
				d_Data,
				fftH,
				fftW,
				dataH,
				dataW,
				kernelH,
				kernelW,
				kernelY,
				kernelX);
		break;
	case pRightWall:
		padDataRightWall(
				d_PaddedData,
				d_Data,
				fftH,
				fftW,
				dataH,
				dataW,
				kernelH,
				kernelW,
				kernelY,
				kernelX);
		break;
	case pMirror:
		padDataMirror(
				d_PaddedData,
				d_Data,
				fftH,
				fftW,
				dataH,
				dataW,
				kernelH,
				kernelW,
				kernelY,
				kernelX);
		break;
	case pConstantVF:
		padDataConstantVF(
				d_PaddedData,
				d_Data,
				fftH,
				fftW,
				dataH,
				dataW,
				kernelH,
				kernelW,
				kernelY,
				kernelX,
				epiBoundary,
				baseChem);
		break;

	}
}


bool computeKernel(
		float           *d_Window,
		int              kernelRadius,
		float            lambda,
		float            gamma,                                 // decay constant
		float            dt,
		c_ctx            cctx)
{

	/********************************************
	 * Declarations and allocations             *
	 ********************************************/
	float t = 0.0;

	int
	cpu_input  = 0,
	cpu_output = 1;

	float
	*h_Data,
	*h_Kernel,
	*h_Window,
	*h_ResultGPU,
	*h_ResultCPU[2];

	float
	*d_Data,
	*d_PaddedData,
	*d_Kernel,
	*d_PaddedKernel;

	fComplex
	*d_DataSpectrum,
	*d_KernelSpectrum;

	cufftHandle
	fftPlanFwd,
	fftPlanInv;

	bool bRetVal = true;
	StopWatchInterface *hTimer = NULL;
	sdkCreateTimer(&hTimer);

#ifdef PRINT_KERNEL
fprintf(stderr, "Testing kernel computation\n");
printf("Testing kernel computation\n");
fprintf(stderr, "\tBuilding filter kernel\n");
#endif
const int    kernelH = cctx.KH;//7;
const int    kernelW = cctx.KW;//6;
const int    kernelY = cctx.KY;//3;
const int    kernelX = cctx.KX;//4;
const int      dataH = cctx.DH;//100;//1160;//2000;
const int      dataW = cctx.DW;//100;//1660;//2000;
const int outKernelH = cctx.DH;
const int outKernelW = cctx.DW;
const int       fftH = cctx.FFTH;
const int       fftW = cctx.FFTW;
// Changed 2
const int    windowH = cctx.windowH;
const int    windowW = cctx.windowW;

#ifdef PRINT_KERNEL
fprintf(stderr, "...allocating memory\n");
#endif
h_Data                  = (float *)malloc(dataH   * dataW * sizeof(float));
h_Kernel                = (float *)malloc(kernelH * kernelW * sizeof(float));
h_Window                = (float *)malloc(windowH * windowW * sizeof(float));
h_ResultGPU             = (float *)malloc(dataH   * dataW * sizeof(float));
h_ResultCPU[cpu_input]  = (float *)malloc(dataH   * dataW * sizeof(float));
h_ResultCPU[cpu_output] = (float *)malloc(dataH   * dataW * sizeof(float));

checkCudaErrors(cudaMalloc((void **)&d_Data,   dataH   * dataW   * sizeof(float)));
checkCudaErrors(cudaMalloc((void **)&d_Kernel, kernelH * kernelW * sizeof(float)));

checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   fftH * fftW * sizeof(float)));
checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fftH * fftW * sizeof(float)));

checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum,   fftH * (fftW / 2 + 1) * sizeof(fComplex)));
checkCudaErrors(cudaMalloc((void **)&d_KernelSpectrum, fftH * (fftW / 2 + 1) * sizeof(fComplex)));

#ifdef PRINT_KERNEL
fprintf(stderr, "...generating 2D %d x %d kernel coefficients\n", kernelH, kernelW);
#endif

/********************************************
 * Initial kernel initialization            *
 ********************************************/
for (int i = 0; i < kernelH * kernelW; i++)
{
	h_Kernel[i] = 0;
}
h_Kernel[0 * kernelW + 1] = lambda;
h_Kernel[1 * kernelW + 0] = lambda;
h_Kernel[1 * kernelW + 2] = lambda;
h_Kernel[2 * kernelW + 1] = lambda;
h_Kernel[1 * kernelW + 1] = 1 - 4*lambda - gamma*dt;

for (int i = 0; i < dataH * dataW; i++)
{
	h_Data[i] = 0;
	h_ResultCPU[cpu_input][i] = 0;
}

// Copy kernel data to middle block of the input
int start_i = outKernelH/2 - kernelH/2;
int end_i   = outKernelH/2 + kernelH/2 + 1;
int start_j = outKernelW/2 - kernelW/2;
int end_j   = outKernelW/2 + kernelW/2 + 1;
int ki = 0, kj = 0;
for (int i = start_i; i < end_i; i++) {
	for (int j = start_j; j < end_j; j++) {
		h_Data                  [i * dataW + j] = h_Kernel[ki * kernelW + kj];
		h_ResultCPU [cpu_input] [i * dataW + j] = h_Kernel[ki * kernelW + kj];
#ifdef PRINT_KERNEL
printf("%d,%d -> %d,%d\n", ki, kj, i, j);
#endif
kj++;
	}
	ki++;
	kj = 0;
}

#ifdef PRINT_KERNEL
fprintf(stderr, "...creating R2C & C2R FFT plans for %i x %i\n", fftH, fftW);
#endif
checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));

#ifdef PRINT_KERNEL
fprintf(stderr, "...uploading to GPU and padding convolution kernel and input data\n");
#endif
sdkResetTimer(&hTimer);
sdkStartTimer(&hTimer);

checkCudaErrors(cudaMemcpy(d_Kernel, h_Kernel, kernelH * kernelW * sizeof(float), cudaMemcpyHostToDevice));
checkCudaErrors(cudaMemcpy(d_Data,   h_Data,   dataH   * dataW *   sizeof(float), cudaMemcpyHostToDevice));

sdkStopTimer(&hTimer);
double dataTransferTime = sdkGetTimerValue(&hTimer);
#ifdef PRINT_KERNEL
//      printf("\tData transfer: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
#endif
sdkResetTimer(&hTimer);
sdkStartTimer(&hTimer);

checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fftH * fftW * sizeof(float)));
checkCudaErrors(cudaMemset(d_PaddedData,   0, fftH * fftW * sizeof(float)));

padKernel(
		d_PaddedKernel,
		d_Kernel,
		fftH,
		fftW,
		kernelH,
		kernelW,
		kernelY,
		kernelX
);

padData2D(
		pClampToBorder,
		d_PaddedData,
		d_Data,
		fftH,
		fftW,
		dataH,
		dataW,
		kernelH,
		kernelW,
		kernelY,
		kernelX,
		-1,
		-1
);

sdkStopTimer(&hTimer);
double memsetPaddingTime = sdkGetTimerValue(&hTimer);
#ifdef PRINT_KERNEL
//      printf("\tMemset and padding: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (memsetPaddingTime * 0.001), memsetPaddingTime);
#endif

//Not including kernel transformation into time measurement,
//since convolution kernel is not changed very frequently
#ifdef PRINT_KERNEL
fprintf(stderr, "...transforming convolution kernel\n");
#endif

double buildKernelTimeTotalGPU = 0;
double buildKernelTimeTotalCPU = 0;
checkCudaErrors(cudaDeviceSynchronize());

/********************************************
 * Kernel Computation                       *
 ********************************************/
// d_KernelSpectrum = FFT(d_PaddedKernel)
checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum));
for (int iter = 0; iter < kernelRadius; iter++)
{

#ifdef PRINT_KERNEL
	fprintf(stderr, "...running GPU Kernel building iteration %d:\n", iter);
	printf("GPU Kernel building iteration %d:\n", iter);
#endif

	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);

	checkCudaErrors(cudaDeviceSynchronize());

	/********************************************
	 * Convolution                              *
	 ********************************************/
	checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedData, (cufftComplex *)d_DataSpectrum));
	modulateAndNormalize(d_DataSpectrum, d_KernelSpectrum, fftH, fftW, 1);
	checkCudaErrors(cufftExecC2R(fftPlanInv, (cufftComplex *)d_DataSpectrum, (cufftReal *)d_PaddedData));

	checkCudaErrors(cudaDeviceSynchronize());
	sdkStopTimer(&hTimer);
	double gpuTime = sdkGetTimerValue(&hTimer);
#ifdef PRINT_KERNEL
printf("\t\tGPU computation: %f MPix/s (%f ms)\n",
		(double)dataH * (double)dataW * 1e-6 / (gpuTime * 0.001), gpuTime);
#endif

sdkResetTimer(&hTimer);
sdkStartTimer(&hTimer);
#ifdef PRINT_KERNEL
fprintf(stderr, "...removing results padding\n");
#endif
unpadResult(
		d_Data,
		d_PaddedData,
		dataH,
		dataW,
		fftH,
		fftW
);

sdkStopTimer(&hTimer);
double unpadTime = sdkGetTimerValue(&hTimer);
#ifdef PRINT_KERNEL
printf("\t\textract results: %f MPix/s (%f ms)\n",
		(double)dataH * (double)dataW * 1e-6 / (unpadTime * 0.001), unpadTime);

fprintf(stderr, "...reading back GPU convolution results\n");
#endif
checkCudaErrors(cudaMemcpy(h_ResultGPU, d_Data, dataH * dataW * sizeof(float), cudaMemcpyDeviceToHost));

sdkResetTimer(&hTimer);
sdkStartTimer(&hTimer);

#ifdef PRINT_KERNEL
fprintf(stderr, "...running reference CPU convolution\n");
#endif

//              convolutionClampToBorderCPU(
//                              h_ResultCPU[cpu_output],
//                              h_ResultCPU[cpu_input],
//                              h_Kernel,
//                              dataH,
//                              dataW,
//                              kernelH,
//                              kernelW,
//                              kernelY,
//                              kernelX
//              );

sdkStopTimer(&hTimer);
double cpuTime = sdkGetTimerValue(&hTimer);
#ifdef PRINT_KERNEL
printf("\t\tCPU computation: %f MPix/s (%f ms)\n",
		(double)dataH * (double)dataW * 1e-6 / (cpuTime * 0.001), cpuTime);
#endif

buildKernelTimeTotalGPU += gpuTime;
buildKernelTimeTotalGPU += unpadTime;

buildKernelTimeTotalCPU += cpuTime;

t += dt;

// Update indices for CPU input/output
cpu_input  = (cpu_input  + 1) % 2;
cpu_output = (cpu_output + 1) % 2;

}

#ifdef PRINT_KERNEL
printf("...extract kernel window from center\n");
#endif
/********************************************************
 * Kernel center extraction - FINAL time domain results *
 ********************************************************/
extractCenter(
		d_Window,
		d_Data,
		dataH,
		dataW,
		windowH,
		windowW
);


#ifdef PRINT_KERNEL
printf("...reading back kernel center from GPU\n");
#endif
checkCudaErrors(cudaMemcpy(h_Window, d_Window, windowH * windowW * sizeof(float), cudaMemcpyDeviceToHost));


for (int i = 0; i < kernelH; i++) {
	for (int j = 0; j < kernelW; j++) {
#ifdef PRINT_KERNEL
		printf(", %f", h_Kernel[i*kernelW + j]);
#endif
	}
#ifdef PRINT_KERNEL
	printf("\n");
#endif
}

#ifdef PRINT_KERNEL
printf("\tData transfer:                %f MPix/s (%f ms)\n",
		(double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
printf("\tMemset and padding:   %f MPix/s (%f ms)\n",
		(double)dataH * (double)dataW * 1e-6 / (memsetPaddingTime * 0.001), memsetPaddingTime);
printf("\tTotal GPU time:               %f MPix/s (%f ms)\n",
		(double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalGPU * 0.001), buildKernelTimeTotalGPU);
printf("\tTotal CPU time:               %f MPix/s (%f ms)\n",
		(double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalCPU * 0.001), buildKernelTimeTotalCPU);


fprintf(stderr, "...comparing the results: ");
#endif


// Update indices for CPU input/output
cpu_input  = (cpu_input  + 1) % 2;
cpu_output = (cpu_output + 1) % 2;

//        bRetVal = compareResults(h_ResultCPU[cpu_output], h_ResultGPU,
//                                dataW, dataH, 1,
//                                1e-6);

#ifdef PRINT_KERNEL
fprintf(stderr, "...shutting down\n");
#endif
sdkDeleteTimer(&hTimer);

/********************************************
 * Pointer deallocations                    *
 ********************************************/

checkCudaErrors(cufftDestroy(fftPlanInv));
checkCudaErrors(cufftDestroy(fftPlanFwd));

checkCudaErrors(cudaFree(d_DataSpectrum));
checkCudaErrors(cudaFree(d_KernelSpectrum));
checkCudaErrors(cudaFree(d_PaddedKernel));
checkCudaErrors(cudaFree(d_Data));
checkCudaErrors(cudaFree(d_Kernel));
checkCudaErrors(cudaFree(d_PaddedData));

free(h_ResultCPU[0]);
free(h_ResultCPU[1]);
free(h_ResultGPU);
free(h_Window);
free(h_Data);
free(h_Kernel);

return bRetVal;
}


bool computeKernelSpectrum(
		fComplex        *d_KernelSpectrum,
		float           *d_Kernel,
		c_ctx           kernel_cctx,
		c_ctx           chem_cctx
)
{

	float
	*d_UnpaddedKernel,
	*d_PaddedKernel;

	// Changed
	cufftHandle
	//      fftPlanFwd,
	//      fftPlanInv;
	fftPlan;

	bool bRetVal = true;

#ifdef PRINT_KERNEL
	printf("Testing kernel spectrum computation\n");
	fprintf(stderr, "Testing kernel spectrum computation\n");
#endif
	const int    kernelH = chem_cctx.KH;//kernel_cctx.DH;//7;
	const int    kernelW = chem_cctx.KW;//kernel_cctx.DW;//6;
	const int    kernelY = chem_cctx.KY;//kernel_cctx.DH / 2;//3;
	const int    kernelX = chem_cctx.KX;//kernel_cctx.DW / 2;//4;
	const int       fftH = chem_cctx.FFTH;
	const int       fftW = chem_cctx.FFTW;

#ifdef PRINT_KERNEL
	printf("\tkernelH: %d\tkernelW: %d\n", kernelH, kernelW);
	printf("\tkernelX: %d\tkernelY: %d\n", kernelX, kernelY);
	printf("\tfftH: %d\tfftW: %d\n", fftH, fftW);

	fprintf(stderr,"...allocating memory\n");
#endif
	checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fftH * fftW * sizeof(float)));
	checkCudaErrors(cudaMalloc((void **)&d_UnpaddedKernel, kernelH * kernelW * sizeof(float)));

	// Changed
#ifdef PRINT_KERNEL
	//      printf("...creating R2C FFT plans for %i x %i\n", fftH, fftW);
#endif
	//      checkCudaErrors(cufftPlan2d(&fftPlanFwd, fftH, fftW, CUFFT_R2C));
	//      checkCudaErrors(cufftPlan2d(&fftPlanInv, fftH, fftW, CUFFT_C2R));
#ifdef PRINT_KERNEL
	printf("...creating C2C FFT plan for %i x %i\n", fftH, fftW / 2);
#endif
	checkCudaErrors(cufftPlan2d(&fftPlan, fftH, fftW / 2, CUFFT_C2C));

#ifdef PRINT_KERNEL
	fprintf(stderr,"...uploading to GPU and padding convolution kernel and input data\n");
#endif
	checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fftH * fftW * sizeof(float)));


	padKernel(
			d_PaddedKernel,
			d_Kernel,
			fftH,
			fftW,
			kernelH,
			kernelW,
			kernelY,
			kernelX
	);

	checkCudaErrors(cudaDeviceSynchronize());

	// Changed
	// d_KernelSpectrum = FFT(d_PaddedKernel)
	//      checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel,       (cufftComplex *)d_KernelSpectrum));

	//CUFFT_INVERSE works just as well...
	const int FFT_DIR = CUFFT_FORWARD;
#ifdef PRINT_KERNEL
	printf("...transforming convolution kernel\n");
#endif
	checkCudaErrors(cufftExecC2C(fftPlan, (cufftComplex *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum, FFT_DIR));

	checkCudaErrors(cudaDeviceSynchronize());

	// Changed
	//      checkCudaErrors(cufftDestroy(fftPlanFwd));
	checkCudaErrors(cufftDestroy(fftPlan));

	checkCudaErrors(cudaFree(d_PaddedKernel));
	checkCudaErrors(cudaFree(d_UnpaddedKernel));


	return bRetVal;
}


bool fftDiffuse2D(
		float           *d_Data,
		fComplex        *d_KernelSpectrum0,
		c_ctx            cctx,
		int              epiBoundary,
		float            baseChem)
{


	float
	*d_PaddedData;

	fComplex
	*d_DataSpectrum0;

	cufftHandle
	fftPlan;

	bool bRetVal = 1;


#ifdef PRINT_KERNEL
printf("Testing GPU chemical diffusion computation\n");
fprintf(stderr,"Testing GPU chemical diffusion computation\n");
#endif
const int    kernelH = cctx.KH;
const int    kernelW = cctx.KW;
const int    kernelY = cctx.KY;
const int    kernelX = cctx.KX;
const int      dataH = cctx.DH;
const int      dataW = cctx.DW;
const int       fftH = cctx.FFTH;
const int       fftW = cctx.FFTW;

#ifdef PRINT_KERNEL
printf("\tkernelH: %d\tkernelW: %d\n", kernelH, kernelW);
printf("\tkernelX: %d\tkernelY: %d\n", kernelX, kernelY);
printf("\tdataH: %d\tdataW: %d\n", dataH, dataH);
printf("\tfftH: %d\tfftW: %d\n", fftH, fftW);

fprintf(stderr,"...allocating memory\n");
#endif

checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   fftH * fftW * sizeof(float)));

checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum0,   fftH * (fftW / 2) * sizeof(fComplex)));

#ifdef PRINT_KERNEL
printf("...creating C2C FFT plan for %i x %i\n", fftH, fftW / 2);
#endif
checkCudaErrors(cufftPlan2d(&fftPlan, fftH, fftW / 2, CUFFT_C2C));

#ifdef PRINT_KERNEL
fprintf(stderr,"...uploading to GPU and padding input data\n");
#endif


checkCudaErrors(cudaMemset(d_PaddedData,   0, fftH * fftW * sizeof(float)));



// DEBUG
//      printf("------ before padding data\n");
//      for (int i = 0; i < 2; i++) {
//          for (int j = 0; j < dataW; j++) {
//              printf(" %f ",h_ChemIn[i * dataW + j]);
//          }
//          printf("\n");
//      }
//      printf("\n");
//        float *h_temp = (float *) malloc(fftH * fftH * sizeof(float));
//        checkCudaErrors(cudaMemcpy(h_temp, d_PaddedData, dataH * dataW * sizeof(float), cudaMemcpyDeviceToHost));
//        printf("------ before padding\n");
//        for (int i = 0; i < 2; i++) {
//            for (int j = 0; j < fftW; j++) {
//                printf(" %f ",h_temp[i * fftW + j]);
//            }
//            printf("\n");
//        }
//        printf("\n");

padData2D(
		pConstantVF,        // pRightWall,  //pMirror,  //pClampToBorder,
		d_PaddedData,
		d_Data,
		fftH,
		fftW,
		dataH,
		dataW,
		kernelH,
		kernelW,
		kernelY,
		kernelX,
		epiBoundary,
		baseChem
		//                        0.00
);

// DEBUG
//        checkCudaErrors(cudaMemcpy(h_temp, d_PaddedData, dataH * dataW * sizeof(float), cudaMemcpyDeviceToHost));
//        printf("------ after padding\n");
//        for (int i = 0; i < 2; i++) {
//            for (int j = 0; j < fftW; j++) {
//                printf(" %f ",h_temp[i * fftW + j]);
//            }
//            printf("\n");
//        }
//        free(h_temp);



#ifdef PRINT_KERNEL
fprintf(stderr,"...performing convolution\n");
#endif


// Changed : Added
//CUFFT_INVERSE works just as well...
const int FFT_DIR = CUFFT_FORWARD;

checkCudaErrors(cudaDeviceSynchronize());
// --------- Computing convolution ------------ begin

// d_DataSpectrum = FFT(d_PaddedData)
checkCudaErrors(cufftExecC2C(fftPlan, (cufftComplex *)d_PaddedData, (cufftComplex *)d_DataSpectrum0, FFT_DIR));

// d_DataSpectrum = d_DataSpectrum * d_KernelSpectrum
#ifdef PRINT_KERNEL
printf( "fftH: %d\tfftW: %d\n", fftH, fftW);
#endif
spProcess2D(d_DataSpectrum0, d_DataSpectrum0, d_KernelSpectrum0, fftH, fftW / 2, FFT_DIR);

// d_PaddedData = IFFT(d_DataSpectrum)                  <------- Output
checkCudaErrors(cufftExecC2C(fftPlan, (cufftComplex *)d_DataSpectrum0, (cufftComplex *)d_PaddedData, -FFT_DIR));

// --------- Computing convolution ------------ end
checkCudaErrors(cudaDeviceSynchronize());


#ifdef PRINT_KERNEL
fprintf(stderr,"...removing results padding\n");
#endif
unpadResult(
		d_Data,
		d_PaddedData,
		dataH,
		dataW,
		fftH,
		fftW
);


#ifdef PRINT_KERNEL
fprintf(stderr,"...reading back GPU convolution results\n");
#endif
checkCudaErrors(cufftDestroy(fftPlan));
#ifdef PRINT_KERNEL
printf("...freeing device pointers\n");
#endif
checkCudaErrors(cudaFree(d_DataSpectrum0));
checkCudaErrors(cudaFree(d_PaddedData));

#ifdef PRINT_KERNEL
printf("...returning to main()\n");
#endif
return bRetVal;
}


#else           // MODEL_3D


void my_sleep(unsigned usec) {
	struct timespec req, rem;
	int err;
	req.tv_sec = usec / 1000000;
	req.tv_nsec = (usec % 1000000) * 1000000000;
	while ((req.tv_sec != 0) || (req.tv_nsec != 0)) {
		if (nanosleep(&req, &rem) == 0)
			break;
		err = errno;
		// Interrupted; continue
		if (err == EINTR) {
			req.tv_sec = rem.tv_sec;
			req.tv_nsec = rem.tv_nsec;
		}
		// Unhandleable error (EFAULT (bad pointer), EINVAL (bad timeval in tv_nsec), or ENOSYS (function not supported))
		break;
	}
}

void reportMemUsageGPU()
{
	// show memory usage of GPU

	size_t free_byte ;

	size_t total_byte ;

	checkCudaErrors(cudaMemGetInfo( &free_byte, &total_byte )) ;


	double free_db = (double)free_byte ;

	double total_db = (double)total_byte ;

	double used_db = total_db - free_db ;

	printf("GPU memory usage: used = %f GB, free = %f GB, total = %f GB\n\n",
			used_db/1024.0/1024.0/1024.0, free_db/1024.0/1024.0/1024.0,
			total_db/1024.0/1024.0/1024.0);
}



bool computeKernel3DBatch(
		int      kernelRadius,
		float   *lambda,
		float   *gamma,
		float    dt,
		c_ctx    kern_cctx,
		c_ctx    chem_cctx
)
{

	int *gpu_id = chem_cctx.gpu_id;

	cufftHandle
	fftPlanFwd,
	fftPlanInv;

	const int       fftD = kern_cctx.FFTD;
	const int       fftH = kern_cctx.FFTH;
	const int       fftW = kern_cctx.FFTW;

	bool bRetVal;


	StopWatchInterface *hTimer = NULL;
	sdkCreateTimer(&hTimer);
	double kernelSpectrumComputationTime = 0.0;

	// Compute all kernels
	for (int ic = 0; ic < N_CHEM; ic++){
		int ig = gpu_id[ic];
		checkCudaErrors(cudaSetDevice(chem_cctx.dev_id[ic]));//ig));
		// Create FFT plans
		// TODO: Make plan plans reusable
		printf("...creating R2C & C2R 3D FFT plans for %i x %i x %i\n", fftD, fftH, fftW);
		printf("\tchem %d on device %d\n", ic, ig);

		checkCudaErrors(cufftPlan3d(&fftPlanFwd, fftD, fftH, fftW, CUFFT_R2C));
		checkCudaErrors(cufftPlan3d(&fftPlanInv, fftD, fftH, fftW, CUFFT_C2R));
		//  reportMemUsageGPU();

		/********************************************
		 * Kernel Computation                       *
		 ********************************************/
		bRetVal = computeKernel3D(
				kernelRadius,
				lambda[ic],
				gamma[ic],
				dt,
				kern_cctx,
				fftPlanFwd,
				fftPlanInv,
				ic
		);

		sdkResetTimer(&hTimer);
		sdkStartTimer(&hTimer);

		/********************************************
		 * Kernel Spectrum Computation              *
		 ********************************************/
		computeKernelSpectrum3D(
				kern_cctx,
				chem_cctx,
				ic);
		sdkStopTimer(&hTimer);
		kernelSpectrumComputationTime += sdkGetTimerValue(&hTimer);

		// Destroy reusable FFT plans
		checkCudaErrors(cufftDestroy(fftPlanInv));
		checkCudaErrors(cufftDestroy(fftPlanFwd));
	}

	printf("\tTotal ker spect computation: %f MPix/s (%f ms)\n",
			(double)chem_cctx.DD * (double)chem_cctx.DH * (double) chem_cctx.DW * 1e-6 /
			(kernelSpectrumComputationTime * 0.001),
			kernelSpectrumComputationTime);

	//  /********************************************
	//   * Pointer deallocations                    *
	//   ********************************************/
	//  for (int ig = 0; ig < N_GPU; ig++)
	//  {
	//    checkCudaErrors(cudaFree(kern_cctx.d_data[ig])); // [1] to [N_CHEM-1] is the same as [0]
	//    checkCudaErrors(cudaFree(kern_cctx.d_kernelspectrum_h[ig]));
	//  }

	printf("returning from compute kernel batch to main()\n");
	return bRetVal;
}

bool computeKernel3D(
		int		 kernelRadius,
		float		 lambda,
		float		 gamma,					// decay constant
		float		 dt,
		c_ctx 		 cctx,
		cufftHandle	 fftPlanFwd,
		cufftHandle	 fftPlanInv,
		short int	 ic)
{

#ifdef COMPUTE_COVERAGE
	int
	cpu_input  = 0,
	cpu_output = 1;
#endif  // COMPUTE_COVERAGE

	float	 *h_Window 			= cctx.h_data[ic];
	float	 *d_Window 			= cctx.d_data[ic];
	fComplex *d_DataSpectrum	= cctx.d_kernelspectrum_h[ic];	// Kernel result spectrum

	float
	*h_Data,
	*h_Kernel,
	*h_ResultGPU;
	//	*h_ResultCPU[2];

	float
	*d_Data,
	*d_PaddedData,
	*d_Kernel,
	*d_PaddedKernel;

	fComplex
	*d_KernelSpectrum;


	bool bRetVal = true;
	StopWatchInterface *hTimer = NULL;
	sdkCreateTimer(&hTimer);
#ifdef PRINT_KERNEL
	printf("Testing kernel computation\n");
	printf("\tBuilding filter kernel\n");
#endif	// PRINT_KERNEL
	const int    niter   = cctx.niter;
	const int    kernelD = cctx.KD;
	const int    kernelH = cctx.KH;//7;
	const int    kernelW = cctx.KW;//6;

	const int    kernelZ = cctx.KZ;
	const int    kernelY = cctx.KY;//3;
	const int    kernelX = cctx.KX;//4;

	const int      dataD = cctx.DD;
	const int      dataH = cctx.DH;//100;//1160;//2000;
	const int      dataW = cctx.DW;//100;//1660;//2000;

	const int outKernelD = cctx.DD;
	const int outKernelH = cctx.DH;
	const int outKernelW = cctx.DW;

	const int       fftD = cctx.FFTD;
	const int       fftH = cctx.FFTH;
	const int       fftW = cctx.FFTW;
	// Changed 2
	const int    windowD = cctx.windowD;
	const int    windowH = cctx.windowH;
	const int    windowW = cctx.windowW;

	int ksize	= kernelD	* kernelH	* kernelW;
	int dsize	= dataD		* dataH		* dataW;
	int fsize	= fftD		* fftH		* fftW;
	int wsize	= windowD	* windowH	* windowW;

	int ksize_b	= ksize * sizeof(float);
	int dsize_b = dsize * sizeof(float);
	int fsize_b = fsize * sizeof(float);
	int wsize_b = wsize * sizeof(float);
#ifdef PRINT_KERNEL
	printf("...allocating memory\n");
#endif	// PRINT_KERNEL
	h_Data      			= (float *)malloc(dsize_b);
	h_Kernel    			= (float *)malloc(ksize_b);

	// DEBUG
	reportMemUsageGPU();
	printf("allocating %d bytes, %d MB\n", dsize_b, dsize_b/(1024*1024));

	checkCudaErrors(cudaMalloc((void **)&d_Data,	dsize_b));
	checkCudaErrors(cudaMalloc((void **)&d_Kernel,	ksize_b));

	checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   fsize_b));
	checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fsize_b));
	printf("k: %p\n", d_PaddedKernel);

	checkCudaErrors(cudaMalloc((void **)&d_KernelSpectrum, fftD * fftH * (fftW / 2 + 1) * sizeof(fComplex)));


	/********************************************
	 * Initial kernel initialization            *
	 ********************************************/

	printf("...generating 3D %d x %d x %d kernel coefficients\n", kernelD, kernelH, kernelW);

	for (int i = 0; i < kernelD * kernelH * kernelW; i++)
	{
		h_Kernel[i] = 0;
	}

	int hStride = kernelW;
	int dStride = kernelH * kernelW;
	h_Kernel[0 * dStride + 1 * hStride + 1] = lambda;
	h_Kernel[2 * dStride + 1 * hStride + 1] = lambda;
	h_Kernel[1 * dStride + 0 * hStride + 1] = lambda;
	h_Kernel[1 * dStride + 2 * hStride + 1] = lambda;
	h_Kernel[1 * dStride + 1 * hStride + 0] = lambda;
	h_Kernel[1 * dStride + 1 * hStride + 2] = lambda;
	h_Kernel[1 * dStride + 1 * hStride + 1] = 1 - 6*lambda - gamma*dt;

	for (int i = 0; i < dataD * dataH * dataW; i++)
	{
		h_Data[i] = 0;
	}

	// Copy kernel data to middle block of the input
	int start_k = outKernelD/2 - kernelD/2;
	int end_k	= outKernelD/2 + kernelD/2 + 1;
	int start_i = outKernelH/2 - kernelH/2;
	int end_i	= outKernelH/2 + kernelH/2 + 1;
	int start_j = outKernelW/2 - kernelW/2;
	int end_j	= outKernelW/2 + kernelW/2 + 1;
	int kk = 0, ki = 0, kj = 0;
	int strideD		= dataH * dataW;
	int strideH		= dataW;
	int kstrideD	= kernelH * kernelW;
	int kstrideH	= kernelW;
	for (int k = start_k; k < end_k; k++) {
		for (int i = start_i; i < end_i; i++) {
			for (int j = start_j; j < end_j; j++) {
				h_Data					[k * strideD + i * strideH + j]	= h_Kernel[kk * kstrideD + ki * kstrideH + kj];
				//				printf("%d,%d,%d -> %d,%d,%d\n", kk, ki, kj, k, i, j);
				kj++;
			}
			ki++;
			kj = 0;
		}
		kk++;
		ki = 0;
	}

#ifdef CALC_MEM_RQ

	const size_t numGPUs = 1;
	size_t workSizeFwd[numGPUs];
	size_t workSizeInv[numGPUs];

	cufftGetSize3d(fftPlanFwd, fftD, fftH, fftW, CUFFT_R2C, workSizeFwd);
	cufftGetSize3d(fftPlanInv, fftD, fftH, fftW, CUFFT_C2R, workSizeInv);

	printf("Compute kernel forward size %d x %d x %d requires %d bytes\n", fftW, fftH, fftD, workSizeFwd[0]);
	printf("Compute kernel bckward size %d x %d x %d requires %d bytes\n", fftW, fftH, fftD, workSizeInv[0]);

	//return true;

#endif

	printf("...uploading to GPU and padding convolution kernel and input data\n");
	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);

	checkCudaErrors(cudaMemcpy(d_Kernel, h_Kernel, ksize_b, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_Data,   h_Data,   dsize_b, cudaMemcpyHostToDevice));

	sdkStopTimer(&hTimer);
	double dataTransferTime = sdkGetTimerValue(&hTimer);
	//	printf("\tData transfer: %f MPix/s (%f ms)\n", (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);

	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);

	checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fsize_b));
	checkCudaErrors(cudaMemset(d_PaddedData,   0, fsize_b));


	double memsetPaddingTime = 0, buildKernelTimeTotalGPU = 0, buildKernelTimeTotalCPU = 0;

	/********************************************
	 * Kernel Computation                       *
	 ********************************************/
	for (int filter_i = 0; filter_i < niter; filter_i++)
	{
		padKernel3D(
				d_PaddedKernel,
				d_Kernel,
				fftD,
				fftH,
				fftW,
				kernelD,
				kernelH,
				kernelW,
				kernelZ,
				kernelY,
				kernelX
		);

		checkCudaErrors(cudaDeviceSynchronize());

		padDataClampToBorder3D(
				d_PaddedData,
				d_Data,
				fftD,
				fftH,
				fftW,
				dataD,
				dataH,
				dataW,
				kernelD,
				kernelH,
				kernelW,
				kernelZ,
				kernelY,
				kernelX
		);

		sdkStopTimer(&hTimer);
		memsetPaddingTime += sdkGetTimerValue(&hTimer);

		//Not including kernel transformation into time measurement,
		//since convolution kernel is not changed very frequently
#ifdef PRINT_KERNEL
printf("...transforming convolution kernel k: %p\tkspec: %p\n",
		d_PaddedKernel, d_KernelSpectrum);
#endif	// PRINT_KERNEL

buildKernelTimeTotalGPU = 0;
buildKernelTimeTotalCPU = 0;
checkCudaErrors(cudaDeviceSynchronize());

// d_KernelSpectrum = FFT(d_PaddedKernel)
checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum));



sdkResetTimer(&hTimer);
sdkStartTimer(&hTimer);

checkCudaErrors(cudaDeviceSynchronize());
#ifdef PRINT_KERNEL
printf("HERE------- %p\t%p\n", d_PaddedData, d_DataSpectrum);
#endif	// PRINT_KERNEL

/********************************************
 * Convolution                              *
 ********************************************/
checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum));
checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedData, (cufftComplex *)d_DataSpectrum));

// Modulate WITHOUT scaling
int isFirstIter = (filter_i == 0)? 1 : 0;

//	    for (int iter = 1; iter < kernelRadius - isFirstIter; iter++)
//	    {
//	        modulate3D(d_DataSpectrum, d_KernelSpectrum, fftD, fftH, fftW, 1);
//	    }

	    complexPower(
	    		d_KernelSpectrum,
	    		fftD,
	    		fftH,
	    		fftW,
	    		1,
	    		kernelRadius - isFirstIter - 1
	    );


// Last iteration: Modulate AND scale
modulateAndNormalize3D(d_DataSpectrum, d_KernelSpectrum, fftD, fftH, fftW, 1);
checkCudaErrors(cufftExecC2R(fftPlanInv, (cufftComplex *)d_DataSpectrum, (cufftReal *)d_PaddedData));

checkCudaErrors(cudaDeviceSynchronize());

#ifdef PRINT_KERNEL
sdkStopTimer(&hTimer);
double gpuTime = sdkGetTimerValue(&hTimer);
printf("\t\tGPU computation: %f MPix/s (%f ms)\n",
		(double)dataD * (double)dataH * (double)dataW * 1e-6 / (gpuTime * 0.001), gpuTime);
#endif	// PRINT_KERNEL

sdkResetTimer(&hTimer);
sdkStartTimer(&hTimer);

unpadResult3D(
		d_Data,
		d_PaddedData,
		dataD,
		dataH,
		dataW,
		fftD,
		fftH,
		fftW
);

#ifdef PRINT_KERNEL
sdkStopTimer(&hTimer);
double unpadTime = sdkGetTimerValue(&hTimer);
printf("\t\tunpad results: %f MPix/s (%f ms)\n",
		(double)dataD * (double)dataH * (double)dataW * 1e-6 / (unpadTime * 0.001), unpadTime);

printf("...reading back GPU convolution results\n");
#endif	// PRINT_KERNEL

	}

	/********************************************************
	 * Kernel center extraction - FINAL time domain results *
	 ********************************************************/
	printf("...extract kernel window from center\n");
	extractCenter3D(
			d_Window,
			d_Data,
			dataD,
			dataH,
			dataW,
			windowD,
			windowH,
			windowW
	);

	printf("...reading back kernel center from GPU\n");
	checkCudaErrors(cudaMemcpy(h_Window, d_Window, wsize_b, cudaMemcpyDeviceToHost));


#ifdef TEST_KERNEL
	for (int k = 0; k < kernelD; k++) {
		for (int i = 0; i < kernelH; i++) {
			for (int j = 0; j < kernelW; j++) {
				printf(", %f", h_Kernel[k*kernelH*kernelD + i*kernelW + j]);
			}
			printf("\n");
		}
	}

	//	checkCudaErrors(cudaMemcpy(h_Data, d_Data, dsize_b, cudaMemcpyDeviceToHost));
	int testW  = 7;
	int xbegin = windowW/2 - testW/2;
	int xend   = windowW - xbegin;
	int ybegin = windowH/2 - testW/2;
	int yend   = windowH - ybegin;
	int zbegin = windowD/2 - testW/2;
	int zend   = windowD - zbegin;
	for (int z = zbegin; z < zend; z++)
		for (int y = ybegin; y < yend; y++)
			for (int x = xbegin; x < xend; x++)
			{
				double rGPU = (double)h_Window[z * windowH * windowW +
				                               y * windowW + x];
				printf("\t[%d,%d,%d] %.*f\n", x, y, z,
						OP_DBL_Digs + 6, rGPU);
			}

	// DEBUG rat
	// Check if filter add up to one
	double ksum = 0.0f;
	for (int z = 0; z < windowD; z++)
		for (int y = 0; y < windowH; y++)
			for (int x = 0; x < windowW; x++)
			{
				double rGPU = (double)h_Window[z * windowH * windowW + y * windowW + x];
				ksum += rGPU;
			}
	cout << "Kernel sum: " << ksum << endl;
	double diff = fabs(ksum - 1.0);
	double precision = 0.001;
	if (diff > precision){
		cout << "Error in kernel computation: Incorrect mass " << diff << " > " << precision << endl;
		exit(-2);
	}

	// print out kernel.vtk file
	util::outputDiffusionKernel(h_Window, windowW, windowH, windowD, "output/kernel.vtk");

#endif	// TEST_KERNEL

	printf("\tData transfer:		%f MPix/s (%f ms)\n",
			(double)dataD * (double)dataH * (double)dataW * 1e-6 / (dataTransferTime * 0.001), dataTransferTime);
	printf("\tMemset and padding:	%f MPix/s (%f ms)\n",
			(double)dataD * (double)dataH * (double)dataW * 1e-6 / (memsetPaddingTime * 0.001), memsetPaddingTime);
	printf("\tTotal GPU time:		%f MPix/s (%f ms)\n",
			(double)dataD * (double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalGPU * 0.001), buildKernelTimeTotalGPU);
	printf("\tTotal CPU time:		%f MPix/s (%f ms)\n",
			(double)dataD * (double)dataH * (double)dataW * 1e-6 / (buildKernelTimeTotalCPU * 0.001), buildKernelTimeTotalCPU);

#ifdef COMPUTE_COVERAGE
	printf("...comparing the results: ");


	// Update indices for CPU input/output
	cpu_input  = (cpu_input  + 1) % 2;
	cpu_output = (cpu_output + 1) % 2;

	printf("Results from GPU:\n");
	displayWindowPlane(h_ResultGPU, dataW, dataH, 15, 10);
	printf("Results from CPU:\n");
	displayWindowPlane(h_ResultCPU[cpu_output], dataW, dataH, 15, 10);



	//#ifdef COMPUTE_COVERAGE
	printf("...computing coverage\n");

	double sum_window;
	for (int z = 0; z < windowD; z++)
		for (int y = 0; y < windowH; y++)
			for (int x = 0; x < windowW; x++)
			{
				double rGPU = (double)h_Window[z * windowH * windowW +
				                               y * windowW + x];
				sum_window += rGPU;
			}

	int hlfW = windowW/2;
	int hlfH = windowH/2;
	int hlfD = windowD/2;
	printf("\tcoverage:\t\t%lf/%lf\t%lf\%\n", sum_window, sum, (sum_window/sum)*100.0);
	printf("\tzero threshold:\tx:[%.*f, %.*f]\n\t\t\ty:[%.*f, %.*f]\n\t\t\tz:[%.*f, %.*f]\n",
			OP_DBL_Digs, h_Window[hlfD * windowH * windowW + hlfH * windowW + 0],
			OP_DBL_Digs, h_Window[hlfD * windowH * windowW + hlfH * windowW + (windowW - 1)],
			OP_DBL_Digs, h_Window[hlfD * windowH * windowW + 0 * windowW + hlfW],
			OP_DBL_Digs, h_Window[hlfD * windowH * windowW + (windowH - 1) * windowW + hlfW],
			OP_DBL_Digs, h_Window[0 * windowH * windowW + hlfH * windowW + hlfW],
			OP_DBL_Digs, h_Window[(windowD - 1) * windowH * windowW + hlfH * windowW + hlfW]);

#endif	// COMPUTE_COVERAGE


	printf("...shutting down\n");
	sdkDeleteTimer(&hTimer);

	/********************************************
	 * Pointer deallocations                    *
	 ********************************************/

	checkCudaErrors(cudaFree(d_Data));
	checkCudaErrors(cudaFree(d_Kernel));
	checkCudaErrors(cudaFree(d_PaddedData));
	checkCudaErrors(cudaFree(d_PaddedKernel));

	// Free after all chems have used this
	checkCudaErrors(cudaFree(d_KernelSpectrum));

	free(h_Data);

	return bRetVal;
}

bool computeKernelSpectrum3D(
		c_ctx 		kernel_cctx,
		c_ctx		chem_cctx,
		short int	ic)
{

	float		*d_Kernel		= kernel_cctx.d_data[ic];
	fComplex	*d_KernelSpectrum	= chem_cctx.d_kernelspectrum_h[ic];
	fComplex	*h_KernelSpectrum	= chem_cctx.h_kernelspectrum[ic];		// permanent storage for padded spectrum on host

	float
	*d_PaddedKernel;

	cufftHandle
	fftPlanFwd;

	bool bRetVal = true;

	printf("Testing kernel spectrum computation\n");
	const int	 kernelD = chem_cctx.KD;
	const int    kernelH = chem_cctx.KH;//kernel_cctx.DH;//7;
	const int    kernelW = chem_cctx.KW;//kernel_cctx.DW;//6;

	const int	 kernelZ = chem_cctx.KZ;
	const int    kernelY = chem_cctx.KY;//kernel_cctx.DH / 2;//3;
	const int    kernelX = chem_cctx.KX;//kernel_cctx.DW / 2;//4;

	const int		fftD = chem_cctx.FFTD;
	const int       fftH = chem_cctx.FFTH;
	const int       fftW = chem_cctx.FFTW;

	size_t fsize_b	= fftD		* fftH		* fftW				* sizeof(float);
	size_t fssize_b = fftD		* fftH		* (fftW / 2 + 1)	* sizeof(fComplex);
	//	size_t ksize_b	= kernelD	* kernelH	* kernelW	* sizeof(float);

	printf("\tkernelD: %d\tkernelH: %d\tkernelW: %d\n", kernelD, kernelH, kernelW);
	printf("\tkernelX: %d\tkernelY: %d\tkernelZ: %d\n", kernelX, kernelY, kernelZ);
	printf("\tfftD: %d\tfftH: %d\tfftW: %d\n", fftD, fftH, fftW);

	printf("...allocating memory\n");
	checkCudaErrors(cudaMalloc((void **)&d_PaddedKernel, fsize_b));


	// TODO: Continue here
	// Changed
	printf("...creating R2C FFT plans for %i x %i x %i\n", fftD, fftH, fftW);
	checkCudaErrors(cufftPlan3d(&fftPlanFwd, fftD, fftH, fftW, CUFFT_R2C));

	printf("...uploading to GPU and padding convolution kernel and input data\n");
	checkCudaErrors(cudaMemset(d_PaddedKernel, 0, fsize_b));


	padKernel3D(
			d_PaddedKernel,
			d_Kernel,
			fftD,
			fftH,
			fftW,
			kernelD,
			kernelH,
			kernelW,
			kernelZ,
			kernelY,
			kernelX
	);

	checkCudaErrors(cudaDeviceSynchronize());


	// Changed
	printf("...transforming convolution kernel\n");
	// d_KernelSpectrum = FFT(d_PaddedKernel)
	checkCudaErrors(cufftExecR2C(fftPlanFwd, (cufftReal *)d_PaddedKernel, (cufftComplex *)d_KernelSpectrum));

	checkCudaErrors(cudaDeviceSynchronize());

	// Transfer data from device to host (permanent storage)
	checkCudaErrors(cudaMemcpy(h_KernelSpectrum,   d_KernelSpectrum, fssize_b, cudaMemcpyDeviceToHost));


	checkCudaErrors(cufftDestroy(fftPlanFwd));
	//checkCudaErrors(cudaFree(d_Kernel));		// Not used after spectrum is calculated
	checkCudaErrors(cudaFree(d_PaddedKernel));


	return bRetVal;
}


bool fftDiffuse3D(
		float    		*d_Data,
		fComplex 		*d_KernelSpectrum,
		cufftHandle  fftPlanFwd,
		cufftHandle  fftPlanInv,
		c_ctx	  		 cctx,
		int          epiBoundary,
		float        baseChem)
{

	int devID;

	StopWatchInterface *hTimer = NULL;
	sdkCreateTimer(&hTimer);
	cudaGetDevice(&devID);

	float
	*d_PaddedData;

	fComplex
	*d_DataSpectrum;

	/*	cufftHandle
	fftPlanFwd,
	fftPlanInv;
	 */
	bool bRetVal = true;
#ifdef PRINT_KERNEL
	printf("Testing GPU chemical diffusion computation\n");
#endif	// PRINT_KERNEL
	const int	 kernelD = cctx.KD;
	const int    kernelH = cctx.KH;
	const int    kernelW = cctx.KW;

	const int	 kernelZ = cctx.KZ;
	const int    kernelY = cctx.KY;
	const int    kernelX = cctx.KX;

	const int	   dataD = cctx.DD;
	const int      dataH = cctx.DH;
	const int      dataW = cctx.DW;

	const int		fftD = cctx.FFTD;
	const int       fftH = cctx.FFTH;
	const int       fftW = cctx.FFTW;

	int ksize	= kernelD	* kernelH	* kernelW;
	int dsize	= dataD		* dataH		* dataW;
	int fsize	= fftD		* fftH		* fftW;

	int ksize_b	= ksize * sizeof(float);
	int dsize_b = dsize * sizeof(float);
	int fsize_b = fsize * sizeof(float);
#ifdef PRINT_KERNEL
	printf("\tkernelD: %d\tkernelH: %d\tkernelW: %d\n", kernelD, kernelH, kernelW);
	printf("\tkernelX: %d\tkernelY: %d\tkernelZ: %d\n", kernelX, kernelY, kernelZ);
	printf("\tdataD: %d\tdataH: %d\tdataW: %d\n", dataD, dataH, dataW);
	printf("\tfftD: %d\tfftH: %d\tfftW: %d\n", fftD, fftH, fftW);

	printf("...allocating memory ------\n");
#endif	// PRINT_KERNEL

	checkCudaErrors(cudaMalloc((void **)&d_PaddedData,   	fsize_b));
	checkCudaErrors(cudaMalloc((void **)&d_DataSpectrum,	fftD * fftH * (fftW / 2 + 1) * sizeof(fComplex)));


#ifdef CALC_MEM_RQ

	const size_t numGPUs = 1;
	size_t workSizeFwd[numGPUs];
	size_t workSizeInv[numGPUs];

	checkCudaErrors(cufftEstimate3d(fftD, fftH, fftW, CUFFT_R2C, workSizeFwd));
	checkCudaErrors(cufftEstimate3d(fftD, fftH, fftW, CUFFT_C2R, workSizeInv));


	printf("Compute chem forward size %d x %d x %d requires %d bytes\n", fftW, fftH, fftD, workSizeFwd[0]);
	printf("Compute chem bckward size %d x %d x %d requires %d bytes\n", fftW, fftH, fftD, workSizeInv[0]);


	//return true;
#endif

#ifdef PRINT_KERNEL
	printf("...creating R2C & C2R FFT plans for %i x %i x %i\n", fftD, fftH, fftW);
#endif //PRINT_KERNEL
	//	checkCudaErrors(cufftPlan3d(&fftPlanFwd, fftD, fftH, fftW, CUFFT_R2C));
	//	checkCudaErrors(cufftPlan3d(&fftPlanInv, fftD, fftH, fftW, CUFFT_C2R));


#ifdef CALC_MEM_RQ


	cufftGetSize3d(fftPlanFwd, fftD, fftH, fftW, CUFFT_R2C, workSizeFwd);
	cufftGetSize3d(fftPlanInv, fftD, fftH, fftW, CUFFT_C2R, workSizeInv);

	printf("Compute kernel forward size %d x %d x %d requires %f GB\n", fftW, fftH, fftD,
			(float) workSizeFwd[0]/(1024.0*1024.0*1024.0));
	printf("Compute kernel bckward size %d x %d x %d requires %f GB\n", fftW, fftH, fftD,
			(float) workSizeInv[0]/(1024.0*1024.0*1024.0));

	//return true;

#endif


	firstTransferCompleted = true;
	checkCudaErrors(cudaMemset(d_PaddedData,    0,	fsize_b));

	/********************************************
	 * Pad data                                 *
	 ********************************************/

	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);

	//	padDataClampToBorder3D(
	//			d_PaddedData,
	//			d_Data,
	//			fftD,
	//			fftH,
	//			fftW,
	//			dataD,
	//			dataH,
	//			dataW,
	//			kernelD,
	//			kernelH,
	//			kernelW,
	//			kernelZ,
	//			kernelY,
	//			kernelX
	//	);


	// DEBUG rat
	padDataConstantVF3D(
			d_PaddedData,
			d_Data,
			fftD,
			fftH,
			fftW,
			dataD,
			dataH,
			dataW,
			kernelD,
			kernelH,
			kernelW,
			kernelZ,
			kernelY,
			kernelX,
			epiBoundary,
			baseChem
	);


#ifdef PRINT_KERNEL
	printf("...performing convolution\n");
#endif	// PRINT_KERNEL

	checkCudaErrors(cudaDeviceSynchronize());

	sdkStopTimer(&hTimer);
	double padTime = sdkGetTimerValue(&hTimer);
	// --------- Computing convolution ------------ begin
	/********************************************
	 * Compute FFT{data}                        *
	 ********************************************/
	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);

	// d_DataSpectrum = FFT{d_PaddedData}
	checkCudaErrors(cufftExecR2C(fftPlanFwd,
			(cufftReal *)d_PaddedData, (cufftComplex *)d_DataSpectrum));

	checkCudaErrors(cudaDeviceSynchronize());

	sdkStopTimer(&hTimer);
	double fftTime = sdkGetTimerValue(&hTimer);

	/********************************************
	 * Spectrum Point-wise Multiplication       *
	 ********************************************/
	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);
	// d_DataSpectrum = d_DataSpectrum * d_KernelSpectrum
	modulateAndNormalize3D(d_DataSpectrum, d_KernelSpectrum, fftD, fftH, fftW, 1);

	checkCudaErrors(cudaDeviceSynchronize());

	sdkStopTimer(&hTimer);
	double multTime = sdkGetTimerValue(&hTimer);

	/********************************************
	 * Compute IFFT{data_spectrum}              *
	 ********************************************/
	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);
	// d_PaddedData = IFFT{d_DataSpectrum}

	checkCudaErrors(cufftExecC2R(fftPlanInv, (cufftComplex *)d_DataSpectrum, (cufftReal *)d_PaddedData));

	checkCudaErrors(cudaDeviceSynchronize());

	sdkStopTimer(&hTimer);
	double ifftTime = sdkGetTimerValue(&hTimer);

	// --------- Computing convolution ------------ end

	/********************************************
	 * Unpad results                            *
	 ********************************************/
	sdkResetTimer(&hTimer);
	sdkStartTimer(&hTimer);
#ifdef PRINT_KERNEL
	printf("...removing result padding\n");
#endif	// PRINT_KERNEL
	unpadResult3D(
			d_Data,
			d_PaddedData,
			dataD,
			dataH,
			dataW,
			fftD,
			fftH,
			fftW
	);

	checkCudaErrors(cudaDeviceSynchronize());

	sdkStopTimer(&hTimer);
	double unpadTime = sdkGetTimerValue(&hTimer);

	/********************************************
	 * Execution Time Display                   *
	 ********************************************/
	printf("\t\t\t\t[%d] pad took:         %f ms\n", devID, padTime);
	printf("\t\t\t\t[%d] FFT took:         %f ms\n", devID, fftTime);
	printf("\t\t\t\t[%d] MULT took:        %f ms\n", devID, multTime);
	printf("\t\t\t\t[%d] IFFT took:        %f ms\n", devID, ifftTime);
	printf("\t\t\t\t[%d] unpad took:       %f ms\n", devID, unpadTime);


	/********************************************
	 * Deallocation of plans and memory         *
	 ********************************************/
	sdkDeleteTimer(&hTimer);

	//	checkCudaErrors(cufftDestroy(fftPlanInv));
	//	checkCudaErrors(cufftDestroy(fftPlanFwd));

#ifdef PRINT_KERNEL
	printf("...freeing device pointers\n");
#endif	// PRINT_KERNEL

	checkCudaErrors(cudaFree(d_PaddedData));
	checkCudaErrors(cudaFree(d_DataSpectrum));

#ifdef PRINT_KERNEL
	printf("...returning to main()\n");
#endif	// PRINT_KERNEL
	return bRetVal;
}


#endif  // MODEL_3D

#endif  // GPU_DIFFUSE (*)
