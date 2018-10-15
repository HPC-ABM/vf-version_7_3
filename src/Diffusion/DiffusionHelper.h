/*
 * DiffusionHelper.h
 *
 *  Created on: Oct 16, 2016
 *      Author: nseekhao
 */

#ifndef DIFFUSIONHELPER_H_
#define DIFFUSIONHELPER_H_

// CUDA runtime
#include <cuda_runtime.h>

// Helper functions for CUDA
#include <helper_functions.h>
#include <helper_cuda.h>


#include "../common.h"
#include "convolutionFFT_common.h"

#ifdef ECV_SAMPLE_CHEM
#include "../Visualization/3D/ChemBufferManager.h"
#endif	// ECV_SAMPLE_CHEM

#include "../FieldVariable/Usr_FieldVariables/WHChemical.h"

class WHChemical;

namespace diffusion_helper
{

void findGPUs();
void initKernelDimensions(c_ctx *kernel_cctx);
void allocKernelBuffers(c_ctx *kernel_cctx);
void deallocConvCtxBuffers(c_ctx *ctx);
void initChemDimensions(c_ctx *chem_cctx, c_ctx kernel_cctx, int nx, int ny, int nz);
void allocChemBuffers(c_ctx *chem_cctx, WHChemical *WHWorldChem, int nx, int ny, int nz);
void prepAndComputeKernelSpectrum(c_ctx* chem_cctx, c_ctx kernel_cctx, float ** h_dWindow);
void prepKernel(c_ctx* chem_cctx, c_ctx kernel_context, float *lambda, float *gamma, int numChem);
void printContext(c_ctx c);

}   // namespace diffusion_helper


#endif /* DIFFUSIONHELPER_H_ */
