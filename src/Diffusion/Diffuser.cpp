/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>

#include "Diffuser.h"


#ifdef GPU_DIFFUSE

Diffuser::Diffuser()
{
    this->nBaseChem     = 8;    // default number of chemical types
    this->chem_cctx     = NULL;
    this->kernel_cctx   = NULL;
    this->WHWorldChem   = WHWorldChem;
}

Diffuser::Diffuser(int nBaseChem, c_ctx *cc, c_ctx *kc, WHChemical *WHWorldChem)
{
    this->nBaseChem     = nBaseChem;
    this->chem_cctx     = cc;
    this->kernel_cctx   = kc;
    this->WHWorldChem   = WHWorldChem;

}

Diffuser::~Diffuser(){}

#ifdef OPT_PINNED_MEM
void Diffuser::diffuseChemGPU(int rightWallIndex){
    cerr << "Diffuse Chem (GPU)" << endl;


    int np          = this->WHWorldChem->np;
    int nBaseChem   = this->nBaseChem;

#ifdef PROFILE_GPU
    StopWatchInterface *hTimer = NULL;
    StopWatchInterface *hTimerTotal = NULL;
    sdkCreateTimer(&hTimer);
    sdkCreateTimer(&hTimerTotal);

    double
    unpackTime       = 0, packTime     = 0,
    dataTransferTime = 0, readbackTime = 0,
    convTime         = 0, totalTime    = 0;

    int dataW = this->chem_cctx->DW;
    int dataH = this->chem_cctx->DH;
    int dataD = this->chem_cctx->DD;

    int fftD = this->chem_cctx->FFTD;
    int fftH = this->chem_cctx->FFTH;
    int fftW = this->chem_cctx->FFTW;

    double gb      = 1024.0 * 1024.0 * 1024.0;
    double data_GB = (dataW * dataH * dataD * sizeof(float)) / gb;
    double fft_GB  = (fftD * fftH * (fftW / 2 + 1) * sizeof(fComplex)) / gb;
#endif


#ifdef PROFILE_GPU
        sdkResetTimer(&hTimerTotal);
        sdkStartTimer(&hTimerTotal);
#endif  // PROFILE_GPU


    // Loop over all types of chemical and perform convolution-based diffusion on GPU

    for (int ic = 0; ic < nBaseChem; ic++){
        cerr << "   Diffusing type " << ic << " of " << nBaseChem - 1 << endl;

        /********************************************
         * Unpack cytokine                          *
         ********************************************/
#ifdef OPT_CHEM
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer);
        sdkStartTimer(&hTimer);
#endif  // PROFILE_GPU

        // Preparing input chem buffer
        this->WHWorldChem->unpackChem(this->chem_cctx, ic);

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer);
        unpackTime = sdkGetTimerValue(&hTimer);
#endif  // PROFILE_GPU
#endif  // OPT_CHEM
#ifdef MODEL_3D
        int ig = this->chem_cctx->gpu_id[ic];
        checkCudaErrors(cudaSetDevice(this->chem_cctx->dev_id[ic]));//ig));
#endif      // MODEL_3D

        /********************************************
         * Pinned/System Host to Device Transfer    *
         ********************************************/
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer);
        sdkStartTimer(&hTimer);
#endif  // PROFILE_GPU
        // Copy from h_ibuffs[ic] to d_data[ic]
        H2D(ic, this->chem_cctx, np);

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer);
        dataTransferTime = sdkGetTimerValue(&hTimer);
#endif  // PROFILE_GPU

        /********************************************
         * Diffusion                                *
         ********************************************/
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer);
        sdkStartTimer(&hTimer);
#endif  // PROFILE_GPU

#ifdef MODEL_3D
// DEBUG rat
        fftDiffuse3D(
                this->chem_cctx->d_data[ic],
                this->chem_cctx->d_kernelspectrum_h[ic],
                *(this->chem_cctx),
                rightWallIndex,
                this->WHWorldChem->pbaseline[ic]);

#else   // MODEL_3D

        if (!fftDiffuse2D(
                this->chem_cctx->d_data[ic],
                this->chem_cctx->d_kernelspectrum_h[ic],
                *(this->chem_cctx),
                rightWallIndex,
                this->WHWorldChem->pbaseline[ic]))
        {
            // TODO: Error handling
            printf("Error: fftDiffuse2D()\n");
        }

#endif  // MODEL_3D

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer);
        convTime = sdkGetTimerValue(&hTimer);
#endif  // PROFILE_GPU

        /********************************************
         * Device to Pinned/System Host Transfer    *
         ********************************************/
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer);
        sdkStartTimer(&hTimer);
#endif  // PROFILE_GPU

        // Copy from d_data[ic] to h_obuffs[ic]
        D2H(ic, this->chem_cctx, np);

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer);
        readbackTime = sdkGetTimerValue(&hTimer);
#endif  // PROFILE_GPU


        /********************************************
         * Pack cytokine                            *
         ********************************************/
#ifdef OPT_CHEM
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer);
        sdkStartTimer(&hTimer);
#endif  // PROFILE_GPU

        this->WHWorldChem->packChem(this->chem_cctx, ic);

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer);
        packTime = sdkGetTimerValue(&hTimer);
#endif  // PROFILE_GPU

#ifdef PROFILE_GPU

        printf("\t\t\tunpack time by %d threads:        %f GB/s (%f ms)\n", PACKTH,
                data_GB / (readbackTime * 0.001), unpackTime);
        printf("\t\tdata transfer:                      %f GB/s (%f ms)\n",
                (data_GB + fft_GB) / (dataTransferTime * 0.001),
                dataTransferTime);
        printf("\t\tGPU chemical diffusion computation: %f MPix/s (%f ms)\n",
                (double)dataD * (double)dataH * (double)dataW * 1e-6 / (convTime * 0.001),
                convTime);
        printf("\t\tread back:                          %f GB/s (%f ms)\n",
                data_GB / (readbackTime * 0.001),
                readbackTime);
        printf("\t\t\tpack time by %d threads:          %f GB/s (%f ms)\n", PACKTH,
                data_GB / (packTime * 0.001), packTime);

#endif  // PROFILE_GPU
#endif  // OPT_CHEM
    }
#ifdef PROFILE_GPU
    sdkStopTimer(&hTimerTotal);
    totalTime = sdkGetTimerValue(&hTimerTotal);

    printf("==>\tTotal time for %d chemical types:\t%f ms\n", nBaseChem, totalTime);
    printf("\t\t\t\t\t\t\t%f ms per chem\n", totalTime/((double) nBaseChem));
#endif
        sdkDeleteTimer(&hTimer);
}

#else	// OPT_PINNED_MEM

void Diffuser::diffuseChemGPU(int rightWallIndex){
    cerr << "Diffuse Chem (GPU)" << endl;


    int np          = this->WHWorldChem->np;
    int nBaseChem   = this->nBaseChem;

#ifdef PROFILE_GPU
    StopWatchInterface *hTimer[N_GPU] = {0};
    for (int ig = 0; ig < N_GPU; ig++)
        sdkCreateTimer(&hTimer[ig]);

    StopWatchInterface *hTimerTotal = NULL;
    sdkCreateTimer(&hTimerTotal);

    double
    unpackTime              =  0 , packTime            =  0,
    dataTransferTime[N_GPU] = {0}, readbackTime[N_GPU] = {0},
    convTime        [N_GPU] = {0}, totalTime           =  0;

    int dataW = this->chem_cctx->DW;
    int dataH = this->chem_cctx->DH;
    int dataD = this->chem_cctx->DD;

    int fftD = this->chem_cctx->FFTD;
    int fftH = this->chem_cctx->FFTH;
    int fftW = this->chem_cctx->FFTW;

    double gb      = 1024.0 * 1024.0 * 1024.0;
    double data_GB = (dataW * dataH * dataD * sizeof(float)) / gb;
    double fft_GB  = (fftD * fftH * (fftW / 2 + 1) * sizeof(fComplex)) / gb;

    sdkResetTimer(&hTimerTotal);
    sdkStartTimer(&hTimerTotal);
#endif	// PROFILE_GPU

#ifdef MODEL_3D

    cufftHandle
    fftPlanFwd[N_GPU],
    fftPlanInv[N_GPU];

    /********************************************
     * Create reuseable FFT plans               *
     ********************************************/
    for (int ig = 0; ig < N_GPU; ig++)
    {
#ifdef PRINT_KERNEL
        printf("...creating R2C & C2R reusable FFT plans for %i x %i x %i on GPU %d\n",
			fftD, fftH, fftW, ig);
#endif //PRINT_KERNEL
        checkCudaErrors(cudaSetDevice(this->chem_cctx->dev_id[ig]));//ig));
        checkCudaErrors(cufftPlan3d(&(fftPlanFwd[ig]), fftD, fftH, fftW, CUFFT_R2C));
        checkCudaErrors(cufftPlan3d(&(fftPlanInv[ig]), fftD, fftH, fftW, CUFFT_C2R));
    }

#endif	// MODEL_3D


    /********************************************
     * Unpack all cytokines                     *
     ********************************************/
#ifdef OPT_CHEM
#ifdef PROFILE_GPU
    sdkResetTimer(&hTimer[0]);
    sdkStartTimer(&hTimer[0]);
#endif  // PROFILE_GPU

    // Preparing input chem buffer
    this->WHWorldChem->unpackAllChem(this->chem_cctx);

#ifdef PROFILE_GPU
    sdkStopTimer(&hTimer[0]);
    unpackTime = sdkGetTimerValue(&hTimer[0]);
#endif  // PROFILE_GPU
#endif  // OPT_CHEM


    // Loop over all types of chemical and perform convolution-based diffusion on GPU
#pragma omp parallel num_threads(N_GPU)
{
#ifdef _OMP
    int tid = omp_get_thread_num();
#else // _OMP
    int tid = 0;
#endif	// _OMP
    for (int ic = 0; ic < nBaseChem; ic++){
        
#ifdef MODEL_3D
        int ig = this->chem_cctx->gpu_id[ic];
	int devID = this->chem_cctx->dev_id[ic];
        if (ig != tid) continue;

        checkCudaErrors(cudaSetDevice(devID));//ig));

#ifdef GPU0_ONLY
	if (ig == 1) continue;
#elif defined(GPU1_ONLY)
	if (ig == 0) continue;
#endif	// GPU0_ONLY

#endif      // MODEL_3D
        cout << "   Diffusing type " << ic << " of " << nBaseChem - 1 << endl;
        cerr << "   Diffusing type " << ic << " of " << nBaseChem - 1 << endl;

        /********************************************
         * Pinned/System Host to Device Transfer    *
         ********************************************/
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer[tid]);
        sdkStartTimer(&hTimer[tid]);
#endif  // PROFILE_GPU
#ifndef M40
        if ((tid != 0) && (ic < N_GPU)) {
            usleep(600 * 1000);
        }
#endif

        // Copy from h_ibuffs[ic] to d_data[ic]
        H2D(ic, this->chem_cctx, np);

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer[tid]);
        dataTransferTime[tid] = sdkGetTimerValue(&hTimer[tid]);
#endif  // PROFILE_GPU

        /********************************************
         * Diffusion                                *
         ********************************************/
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer[tid]);
        sdkStartTimer(&hTimer[tid]);
#endif  // PROFILE_GPU

#ifdef MODEL_3D
        // DEBUG rat
        fftDiffuse3D(
                this->chem_cctx->d_data[ic],
                this->chem_cctx->d_kernelspectrum_h[ic],
                fftPlanFwd[ig],
                fftPlanInv[ig],
                *(this->chem_cctx),
                rightWallIndex,
                this->WHWorldChem->pbaseline[ic]);

#else   // MODEL_3D

        if (!fftDiffuse2D(
                this->chem_cctx->d_data[ic],
                this->chem_cctx->d_kernelspectrum_h[ic],
                *(this->chem_cctx),
                rightWallIndex,
                this->WHWorldChem->pbaseline[ic]))
        {
            // TODO: Error handling
            printf("Error: fftDiffuse2D()\n");
        }

#endif  // MODEL_3D

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer[tid]);
        convTime[tid] = sdkGetTimerValue(&hTimer[tid]);
#endif  // PROFILE_GPU

        /********************************************
         * Device to Pinned/System Host Transfer    *
         ********************************************/
#ifdef PROFILE_GPU
        sdkResetTimer(&hTimer[tid]);
        sdkStartTimer(&hTimer[tid]);
#endif  // PROFILE_GPU

        // Copy from d_data[ic] to h_obuffs[ic]
        D2H(ic, this->chem_cctx, np);

#ifdef PROFILE_GPU
        sdkStopTimer(&hTimer[tid]);
        readbackTime[tid] = sdkGetTimerValue(&hTimer[tid]);
#endif  // PROFILE_GPU


#ifdef PROFILE_GPU

        printf("\t\t[%d]data transfer:                      %f GB/s (%f ms)\n",
                tid,
                (data_GB + fft_GB) / (dataTransferTime[tid] * 0.001),
                dataTransferTime[tid]);
        printf("\t\t[%d]GPU chemical diffusion computation: %f MPix/s (%f ms)\n",
                tid,
                (double)dataD * (double)dataH * (double)dataW * 1e-6 / (convTime[tid] * 0.001),
                convTime[tid]);
        printf("\t\t[%d]read back:                          %f GB/s (%f ms)\n",
                tid,
                data_GB / (readbackTime[tid] * 0.001),
                readbackTime[tid]);

#endif  // PROFILE_GPU

    }
}

#ifdef MODEL_3D
    for (int ig = 0; ig < N_GPU; ig++)
    {
        checkCudaErrors(cufftDestroy(fftPlanInv[ig]));
        checkCudaErrors(cufftDestroy(fftPlanFwd[ig]));
    }
#endif
    /********************************************
     * Pack all cytokines                       *
     ********************************************/
#ifdef OPT_CHEM
#ifdef PROFILE_GPU
    sdkResetTimer(&hTimer[0]);
    sdkStartTimer(&hTimer[0]);
#endif  // PROFILE_GPU

    this->WHWorldChem->packAllChem(this->chem_cctx);

#ifdef PROFILE_GPU
    sdkStopTimer(&hTimer[0]);
    packTime = sdkGetTimerValue(&hTimerPack[0]);

    printf("\t\t\tunpack time by %d threads:        %f GB/s (%f ms)\n", PACKTH,
            data_GB / (unpackTime * 0.001), unpackTime);
    printf("\t\t\tpack time by %d threads:          %f GB/s (%f ms)\n", PACKTH,
            data_GB / (packTime * 0.001), packTime);
#endif  // PROFILE_GPU
#endif  // OPT_CHEM

#ifdef PROFILE_GPU


    sdkStopTimer(&hTimerTotal);
    totalTime = sdkGetTimerValue(&hTimerTotal);

    printf("==>\tTotal time for %d chemical types:\t%f ms\n", nBaseChem, totalTime);
    printf("\t\t\t\t\t\t\t%f ms per chem\n", totalTime/((double) nBaseChem));

    sdkDeleteTimer(&hTimerTotal);
    for (int ig = 0; ig < N_GPU; ig++)
        sdkDeleteTimer(&hTimer[ig]);
#endif
    cerr << "finished diffusion for all chems" << endl;
}
#endif	// OPT_PINNED_MEM

#else   // GPU_DIFFUSE

// Discretization of PDE Diffusion Equation using central difference approximation
// Note: diffuseChem() will very likely get replaced by a new correct version, thus this doesn't need comments just yet
void Diffuser::diffuseChemCPU(float dt){
#ifdef OPT_CHEM
    return;
#else       // OPT_CHEM
    //      cerr << " Diffuse Chem " << endl;

    int epiBoundary = this->bmcx-1;
    float coeff;
    // NS
    int wi_offset = nBaseChem;              // chem write index offset

    // For all chemical types
    for(int ichem=0; ichem<nBaseChem; ichem++){


        switch(ichem){  // cytokine specific diffusion coefficient
        case TNF:
            //                                coeff = 0.0018;                 // Diffusion coefficient (mm^2/min)
            //Examples of diffusion coefficients can be found at
            //http://www.math.ubc.ca/~ais/website/status/diffuse.html
            coeff = 0.0009;
            break;
        case TGF:
            //                                coeff = 0.00156;
            coeff = 0.00078;
            break;
        case FGF:
            //                                coeff = 0.00156;
            coeff = 0.00078;
            break;
        case MMP8:
            //                                coeff = 0.00156;
            coeff = 0.00078;
            break;
        case IL1beta:
            //                                coeff = 0.0018;
            coeff = 0.0009;
            break;
        case IL6:
            //                                coeff = 0.00162;
            coeff = 0.00081;
            break;
        case IL8:
            //                                coeff = 0.0018;
            coeff = 0.0009;
            break;
        case IL10:
            //                                coeff = 0.0018;
            coeff = 0.0009;
        }



        // Calculate change in concentration over dt at each patch
        float* tempPtr = new float[nx*ny*nz];
#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            if (tid == 0)
                *ntp = omp_get_num_threads();
            start_times[tid] = omp_get_wtime();
#pragma omp for nowait
#else
#ifdef V_a
#pragma omp parallel for schedule(dynamic)
#else
#pragma omp parallel
            {
#pragma omp for
#endif  // V_a or V_b
#endif



                // Calculate central difference approximation at each patch along each direction
                for(int yi=1; yi<ny-1; yi++){
#pragma omp simd
                    for(int xi=1; xi<epiBoundary; xi++){
#ifdef MODEL_3D
                        for(int zi=1; zi < nz-1; zi++){
#else
                            int zi = 0;
#endif
                            int index = xi+ yi*nx + zi*nx*ny;
                            REAL d;
                            int XPlusOne = (xi+1) + yi*nx + zi*nx*ny;
                            int XMinusOne = (xi-1) + yi*nx + zi*nx*ny;
                            int YPlusOne = xi + (yi+1)*nx + zi*nx*ny;
                            int YMinusOne = xi + (yi-1)*nx + zi*nx*ny;
#ifdef MODEL_3D
int ZPlusOne = xi + yi*nx + (zi+1)*nx*ny;
int ZMinusOne = xi + yi*nx + (zi-1)*nx*ny;
#endif

// Central Difference Approximation along x direction
float d2phi_dx2 = (this->chemAllocation[ichem][XPlusOne] - 2.*this->chemAllocation[ichem][index] + this->chemAllocation[ichem][XMinusOne])/(dx*dx);
// Central Difference Approximation along y direction
float d2phi_dy2 = (this->chemAllocation[ichem][YPlusOne] - 2.*this->chemAllocation[ichem][index] + this->chemAllocation[ichem][YMinusOne])/(dy*dy);
// 2D Central Differnce Approximation of Diffusion
// NS
tempPtr[index] = dt*coeff*(d2phi_dx2 + d2phi_dy2);
#ifdef MODEL_3D
// Central Difference Approximation along z direction
float d2phi_dz2 = (this->chemAllocation[ichem][ZPlusOne] - 2.*this->chemAllocation[ichem][index] + this->chemAllocation[ichem][ZMinusOne])/(dz*dz);
// 3D Central Differnce Approximation of Diffusion
tempPtr[index] = dt*coeff*(d2phi_dx2 + d2phi_dy2 + d2phi_dz2);
                        }
#endif
                    }
                }

#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
                end_times[tid] = omp_get_wtime();
                elapsed1[tid] = end_times[tid] - start_times[tid];
#pragma omp barrier
#endif


#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
                start_times[tid] = omp_get_wtime();
#pragma omp for nowait
#else
#ifdef V_a
#pragma omp parallel for schedule(dynamic)
#else
#pragma omp for
#endif  // V_a or V_b
#endif
                // NS
                // Update concentration from central difference scheme
                for(int yi=0; yi<ny; yi++){
#pragma omp simd
                    for(int xi=0; xi<epiBoundary+1; xi++){
#ifdef MODEL_3D
                        for(int zi=0; zi<nz; zi++){
#else
                            int zi = 0;
#endif
                            int index = xi+ yi*nx + zi*nx*ny;
                            if (yi == 0 || yi == ny-1 || xi == 0 || xi == epiBoundary) { // constant padding boundary condition
#ifdef MODEL_3D
                                if (zi == 0 || zi == nz-1) {
#endif
                                    int countTissue = WHWorld::initialTissue;
                                    this->chemAllocation[ichem][index] = this->baselineChem[ichem]/countTissue;
#ifdef MODEL_3D
                                }
#endif
                            } else {
                                this->chemAllocation[ichem][index] += tempPtr[index];

                            }

#ifdef MODEL_3D
                        }
#endif
                    }
                }
            }
#ifdef PROFILE_THREAD_LEVEL_CHEM_DIFF
            end_times[tid] = omp_get_wtime();
            elapsed2[tid] = end_times[tid] - start_times[tid];
#pragma omp barrier
        }

        cout << "Chemical " << ichem << ":" << endl;
        for(int t = 0; t < num_threads; t++) {
            cout << "       thread " << t << " took: [" << elapsed1[t] << "]        [" <<
                    elapsed2[t] << "]" << endl;
        }
#endif

        delete[] tempPtr;

    }
#endif      // OPT_CHEM
}



#endif  // GPU_DIFFUSE
