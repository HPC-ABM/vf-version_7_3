/*
 * convolutionFFT_common.h
 *
 *  Created on: Oct 16, 2016
 *      Author: nseekhao
 */

#ifndef CONVOLUTIONFFT_COMMON_H_
#define CONVOLUTIONFFT_COMMON_H_

#include "../enums.h"
#include <cufft.h>


typedef unsigned int uint;

#ifdef __CUDACC__
typedef float2 fComplex;
#else
typedef struct
{
    float x;
    float y;
} fComplex;
#endif


////////////////////////////////////////////////////////////////////////////////
// Helper functions
////////////////////////////////////////////////////////////////////////////////
//Round a / b to nearest higher integer value
inline int iDivUp(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

//Align a to nearest higher multiple of b
inline int iAlignUp(int a, int b)
{
    return (a % b != 0) ? (a - a % b + b) : a;
}

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
        float  baseChem);

extern "C" void convolutionClampToBorderCPU(
    float *h_Result,
    float *h_Data,
    float *h_Kernel,
    int dataH,
    int dataW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
);

extern "C" void extractCenter(
        float *d_Dst,
        float *d_Src,
        int dataH,
        int dataW,
        int windowH,
        int windowW
);

extern "C" void unpadResult(
        float *d_Dst,
        float *d_Src,
        int dataH,
        int dataW,
        int fftH,
        int fftW
);

extern "C" void extractCenter3D(
            float *d_Dst,
            float *d_Src,
            int dataD,
            int dataH,
            int dataW,
            int windowD,
            int windowH,
            int windowW
);

extern "C" void unpadResult3D(
            float *d_Dst,
            float *d_Src,
            int dataD,
            int dataH,
            int dataW,
            int fftD,
            int fftH,
            int fftW
);



extern "C" void padKernel(
    float *d_PaddedKernel,
    float *d_Kernel,
    int fftH,
    int fftW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
);

extern "C" void padDataConstantVF(
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
    float  baseChem
);

extern "C" void padDataRightWall(
    float *d_PaddedData,
    float *d_Data,
    int fftH,
    int fftW,
    int dataH,
    int dataW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
);

extern "C" void padDataMirror(
    float *d_PaddedData,
    float *d_Data,
    int fftH,
    int fftW,
    int dataH,
    int dataW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
);

extern "C" void padKernel3D(
        float *d_Dst,
        float *d_Src,
        int fftD,
        int fftH,
        int fftW,
        int kernelD,
        int kernelH,
        int kernelW,
        int kernelZ,
        int kernelY,
        int kernelX
);

extern "C" void padDataClampToBorder(
    float *d_PaddedData,
    float *d_Data,
    int fftH,
    int fftW,
    int dataH,
    int dataW,
    int kernelH,
    int kernelW,
    int kernelY,
    int kernelX
);

extern "C" void padDataClampToBorder3D(
    float *d_PaddedData,
    float *d_Data,
    int fftD,
    int fftH,
    int fftW,
    int dataD,
    int dataH,
    int dataW,
    int kernelD,
    int kernelH,
    int kernelW,
    int kernelZ,
    int kernelY,
    int kernelX
);

extern "C" void padDataConstantVF3D(
		float *d_Dst,
		float *d_Src,
		int fftD,
		int fftH,
		int fftW,
		int dataD,
		int dataH,
		int dataW,
		int kernelD,
		int kernelH,
		int kernelW,
		int kernelZ,
		int kernelY,
		int kernelX,
		int epiBoundary,
		float  baseChem
);


extern "C" void modulateAndNormalize(
    fComplex *d_Dst,
    fComplex *d_Src,
    int fftH,
    int fftW,
    int padding
);

extern "C" void complexPower(
    fComplex *d_Dst,
    int fftD,
    int fftH,
    int fftW,
    int padding,
    int n
);

extern "C" void modulate3D(
    fComplex *d_Dst,
    fComplex *d_Src,
    int fftD,
    int fftH,
    int fftW,
    int padding
);

extern "C" void modulateAndNormalize3D(
    fComplex *d_Dst,
    fComplex *d_Src,
    int fftD,
    int fftH,
    int fftW,
    int padding
);

extern "C" void spPostprocess2D(
    void *d_Dst,
    void *d_Src,
    uint DY,
    uint DX,
    uint padding,
    int dir
);

extern "C" void spPreprocess2D(
    void *d_Dst,
    void *d_Src,
    uint DY,
    uint DX,
    uint padding,
    int dir
);

extern "C" void spProcess2D(
    void *d_Data,
    void *d_Data0,
    void *d_Kernel0,
    uint DY,
    uint DX,
    int dir
);


#define KCoeffsH                3
#define KCoeffsW                3

#define N_CHEM                  8

#ifdef MODEL_3D
#define N_GPU                   2
#define KCoeffsD                3
#else
#define N_GPU                   1
#define KCoeffsD                1
#endif

typedef struct CCTX             // convolution context
{
    float dt;
    float dx2;
    int kernelRadius;
    int niter;

    int KD;
    int KH;
    int KW;

    int KZ;
    int KX;
    int KY;

    int DD;
    int DH;
    int DW;

    int FFTD;
    int FFTH;
    int FFTW;

    int windowD;
    int windowH;
    int windowW;

    float *h_ibuffs[N_CHEM];            // Pinned/system buffer
    float *h_obuffs[N_CHEM];            // Pinned/system buffer

    int    gpu_id[N_CHEM];	// index
    int    dev_id[N_CHEM];	// actual ID

    float *d_data[N_CHEM];              // Used in computing kernel to store extracted window
                                        // and fftDiffuse3D to hold input data
    float *h_data[N_CHEM];              // Used in computing kernel to store extracted window

    fComplex *h_kernelspectrum[N_CHEM];
    fComplex *d_kernelspectrum_h[N_CHEM];
} c_ctx;



extern "C" int snapTransformSize(int dataSize);
extern "C" void H2D(int ic, c_ctx* chem_cctx, int np);
extern "C" void D2H(int ic, c_ctx* chem_cctx, int np);

#ifndef MODEL_3D

extern "C" bool computeKernel(
                float           *d_Window,
                int              kernelRadius,
                float            lambda,
                float            gamma,                                 // decay constant
                float            dt,
                c_ctx            cctx);


extern "C" bool computeKernelSpectrum(
                fComplex        *d_KernelSpectrum,
                float           *d_Kernel,
                c_ctx           kernel_cctx,
                c_ctx           chem_cctx
                );

extern "C" bool fftDiffuse2D(
                float           *d_Data,
                fComplex        *d_KernelSpectrum0,
                c_ctx            cctx,
                int              epiBoundary,
                float            baseChem);



#else   // MODEL_3D

//#define N_ITER                  10
#define N_GPU                   2

// Chem GPU ownership;
extern "C" int gpu_id[N_CHEM];

extern "C" void reportMemUsageGPU();


extern "C" bool computeKernel3DBatch(
        int,
        float[N_CHEM],
        float[N_CHEM],
        float,
        c_ctx,
        c_ctx);

extern "C" bool computeKernel3D(
        int,
        float,
        float,
        float,
        c_ctx,
        cufftHandle,
        cufftHandle,
        short int);

extern "C" bool computeKernelSpectrum3D(
        c_ctx       kernel_cctx,
        c_ctx       chem_cctx,
        short int   ic);

extern "C" bool fftDiffuse3D(
        float       *d_Data,
        fComplex    *d_KernelSpectrum,
        cufftHandle fftPlanFwd,
        cufftHandle fftPlanInv,
        c_ctx       cctx,
        int         epiBoundary,
        float       baseChem);


#endif  // MODEL_3D





#endif /* CONVOLUTIONFFT_COMMON_H_ */
