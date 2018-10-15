/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

// Modified by Nuttiiya Seekhao to support volume rendering of float value
// from main memory

// Simple 3D volume renderer

#ifndef _VOLUMERENDER_KERNEL_CU_
#define _VOLUMERENDER_KERNEL_CU_

#include <assert.h>
#include <helper_cuda.h>
#include <helper_math.h>

#include "../../enums.h"

// Note: Originally, including common_vis worked. However, it stopped working due to GLX inclusion.
// TODO: Separate GLX includes and definitions
#include "../../common.h"
//#include "./common_vis.h"
//#include "./VolumeManager.h"


typedef unsigned int  uint;
typedef unsigned char uchar;

//typedef unsigned char VolumeType;
typedef float VolumeType;

VolumeType *d_svBuffer[m_ecmtotal] = {0};
//cudaArray *d_svArray[m_ecmtotal] = {0};
cudaArray *d_volumeArray[m_ecmtotal] = {0};
#ifdef ECV_SAMPLE_CHEM
cudaArray *d_chemsample_h[TOTAL_CHEM];
#endif



cudaArray *d_transferFuncArrayCol = {0};
cudaArray *d_transferFuncArrayEla = {0};
cudaArray *d_transferFuncArrayHya = {0};


#ifdef AVEP
surface<void, cudaSurfaceType3D> srfCol;
surface<void, cudaSurfaceType3D> srfEla;
surface<void, cudaSurfaceType3D> srfHya;
#endif	// AVEP
texture<VolumeType, 3, cudaReadModeElementType> texCol;
texture<VolumeType, 3, cudaReadModeElementType> texEla;
texture<VolumeType, 3, cudaReadModeElementType> texHya;

//texture<VolumeType, 3, cudaReadModeNormalizedFloat> tex;         // 3D texture
texture<float4, 1, cudaReadModeElementType>     transferTexCol; // 1D transfer function texture
texture<float4, 1, cudaReadModeElementType>     transferTexEla;
texture<float4, 1, cudaReadModeElementType>     transferTexHya;

#ifdef ECV_SAMPLE_CHEM
texture<VolumeType, 3, cudaReadModeElementType> texChem0;
texture<VolumeType, 3, cudaReadModeElementType> texChem1;
texture<VolumeType, 3, cudaReadModeElementType> texChem2;
texture<VolumeType, 3, cudaReadModeElementType> texChem3;
texture<VolumeType, 3, cudaReadModeElementType> texChem4;
texture<VolumeType, 3, cudaReadModeElementType> texChem5;
texture<VolumeType, 3, cudaReadModeElementType> texChem6;
texture<VolumeType, 3, cudaReadModeElementType> texChem7;

surface<void, cudaSurfaceType3D>                srfChem0;
surface<void, cudaSurfaceType3D>                srfChem1;
surface<void, cudaSurfaceType3D>                srfChem2;
surface<void, cudaSurfaceType3D>                srfChem3;
surface<void, cudaSurfaceType3D>                srfChem4;
surface<void, cudaSurfaceType3D>                srfChem5;
surface<void, cudaSurfaceType3D>                srfChem6;
surface<void, cudaSurfaceType3D>                srfChem7;

texture<float4, 1, cudaReadModeElementType>	    transferTexChem0;
texture<float4, 1, cudaReadModeElementType>	    transferTexChem1;
texture<float4, 1, cudaReadModeElementType>	    transferTexChem2;
texture<float4, 1, cudaReadModeElementType>	    transferTexChem3;
texture<float4, 1, cudaReadModeElementType>	    transferTexChem4;
texture<float4, 1, cudaReadModeElementType>	    transferTexChem5;
texture<float4, 1, cudaReadModeElementType>	    transferTexChem6;
texture<float4, 1, cudaReadModeElementType>	    transferTexChem7;

cudaArray *d_transferFuncArrayChem0 = {0};
cudaArray *d_transferFuncArrayChem1 = {0};
cudaArray *d_transferFuncArrayChem2 = {0};
cudaArray *d_transferFuncArrayChem3 = {0};
cudaArray *d_transferFuncArrayChem4 = {0};
cudaArray *d_transferFuncArrayChem5 = {0};
cudaArray *d_transferFuncArrayChem6 = {0};
cudaArray *d_transferFuncArrayChem7 = {0};

#endif	// ECV_SAMPLE_CHEM

typedef struct
{
    float4 m[3];
} float3x4;

__constant__ float3x4 c_invViewMatrix;  // inverse view matrix

struct Ray
{
    float3 o;   // origin
    float3 d;   // direction
};

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

//Round a / b to nearest higher integer value
int iDivUp_AVEP(int a, int b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

#ifdef ECV_SAMPLE_CHEM
#ifdef ECV_SAMPLE_CHEM_TEST

__device__ float smult[TOTAL_CHEM] =
	 {
			 50000.0f,
			 10000.0f,
			 1000000.0f,
			 100.0f,
			 10000.0f,
			 100000.0f,
			 100000.0f,
			 10000.0f
	 };	// sample multiplier

__global__ void sampleChem_kernel(
		float *d_Src,
		int dataD,
		int dataH,
		int dataW,
		int chemIndex
		)
{

	const int z = blockDim.z * blockIdx.z + threadIdx.z;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;
	const int x = blockDim.x * blockIdx.x + threadIdx.x;

	const bool validZ = (0 <= z) && (z < dataD);
	const bool validY = (0 <= y) && (y < dataH);
	const bool validX = (0 <= x) && (x < dataW);

	const bool validZ_h = validZ && (z%ECV_SAMPLE_STRIDE_HGH == 0);
	const bool validY_h = validY && (y%ECV_SAMPLE_STRIDE_HGH == 0);
	const bool validX_h = validX && (x%ECV_SAMPLE_STRIDE_HGH == 0);

	const bool validZ_l = validZ && (z%ECV_SAMPLE_STRIDE_LOW == 0);
	const bool validY_l = validY && (y%ECV_SAMPLE_STRIDE_LOW == 0);
	const bool validX_l = validX && (x%ECV_SAMPLE_STRIDE_LOW == 0);


	const int sampleW_l = dataW / ECV_SAMPLE_STRIDE_LOW;
	const int sampleH_l = dataH / ECV_SAMPLE_STRIDE_LOW;
	const int sampleD_l = dataD / ECV_SAMPLE_STRIDE_LOW;

	const int sampleW_h = dataW / ECV_SAMPLE_STRIDE_HGH;
	const int sampleH_h = dataH / ECV_SAMPLE_STRIDE_HGH;
	const int sampleD_h = dataD / ECV_SAMPLE_STRIDE_HGH;

	int dx_l = x/ECV_SAMPLE_STRIDE_LOW;
	int dy_l = y/ECV_SAMPLE_STRIDE_LOW;
	int dz_l = z/ECV_SAMPLE_STRIDE_LOW;

	int dx_h = x/ECV_SAMPLE_STRIDE_HGH;
	int dy_h = y/ECV_SAMPLE_STRIDE_HGH;
	int dz_h = z/ECV_SAMPLE_STRIDE_HGH;

	const bool validDx_l = (dx_l < sampleW_l);
	const bool validDy_l = (dy_l < sampleH_l);
	const bool validDz_l = (dz_l < sampleD_l);

	const bool validDx_h = (dx_h < sampleW_h);
	const bool validDy_h = (dy_h < sampleH_h);
	const bool validDz_h = (dz_h < sampleD_h);

	if (validZ_h && validY_h && validX_h && validDx_h && validDy_h && validDz_h)
	{

		float sample = d_Src[z * dataH * dataW + y * dataW + x] * 50000000.0f;//*smult[chemIndex];
		if (sample > 1.0f) sample = 1.0f;
//		if (sample < 0.01f) sample = 1.0f;

		surf3Dwrite(sample, srfChem0, dx_h * sizeof(float), dy_h, dz_h);
	}

	if (validZ_l && validY_l && validX_l && validDx_l && validDy_l && validDz_l)
	{

		float sample = d_Src[z * dataH * dataW + y * dataW + x] * 50000000.0f;//*smult[chemIndex];
		if (sample > 1.0f) sample = 1.0f;
//		if (sample < 0.001f) sample = 1.0f;

		surf3Dwrite(sample, srfChem2, dx_l * sizeof(float), dy_l, dz_l);
	}

}
#else	// ECV_SAMPLE_CHEM_TEST
__global__ void sampleChem_kernel(
		float *d_Src,
		int dataD,
		int dataH,
		int dataW,
		int chemIndex
		)
{

	const int z = blockDim.z * blockIdx.z + threadIdx.z;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;
	const int x = blockDim.x * blockIdx.x + threadIdx.x;

	const bool validZ = (0 <= z) && (z < dataD) && (z%ECV_SAMPLE_STRIDE == 0);
	const bool validY = (0 <= y) && (y < dataH) && (y%ECV_SAMPLE_STRIDE == 0);
	const bool validX = (0 <= x) && (x < dataW) && (x%ECV_SAMPLE_STRIDE == 0);

//	const bool validz = (z >= 0) && (z < datad);
//	const bool validy = (y >= 0) && (y < datah);
//	const bool validx = (x >= 0) && (x < dataw);

	const int sampleW = dataW / ECV_SAMPLE_STRIDE;
	const int sampleH = dataH / ECV_SAMPLE_STRIDE;
	const int sampleD = dataD / ECV_SAMPLE_STRIDE;

	int dx = x/ECV_SAMPLE_STRIDE;
	int dy = y/ECV_SAMPLE_STRIDE;
	int dz = z/ECV_SAMPLE_STRIDE;

	const bool validDx = (dx < sampleW);
	const bool validDy = (dy < sampleH);
	const bool validDz = (dz < sampleD);

	if (validZ && validY && validX && validDx && validDy && validDz)
	{

		float sample = d_Src[z * dataH * dataW + y * dataW + x];

		switch (chemIndex)
		{
		case 0:
			sample *= 50000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem0, dx * sizeof(float), dy, dz);
			break;
		case 1:
			sample *= 10000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem1, dx * sizeof(float), dy, dz);
			break;
		case 2:
			sample *= 1000000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem2, dx * sizeof(float), dy, dz);
			break;
		case 3:
			sample *= 100.0f;//5000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem3, dx * sizeof(float), dy, dz);
			break;
		case 4:
			sample *= 10000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem4, dx * sizeof(float), dy, dz);
			break;
		case 5:
			sample *= 100000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem5, dx * sizeof(float), dy, dz);
			break;
		case 6:
			sample *= 100000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem6, dx * sizeof(float), dy, dz);
			break;
		case 7:
			sample *= 10000.0f;
			if (sample > 1.0f) sample = 1.0f;
			surf3Dwrite(sample, srfChem7, dx * sizeof(float), dy, dz);
			break;
		}

	}

}
#endif	// ECV_SAMPLE_CHEM_TEST

extern "C" void sampleChem(
	    float *d_Src,
	    int dataD,
	    int dataH,
	    int dataW,
	    int chemIndex
)
{
    dim3 threads(8, 8, 4);
    dim3 grid(iDivUp_AVEP(dataW, threads.x), iDivUp_AVEP(dataH, threads.y), iDivUp_AVEP(dataD, threads.z));

    printf("   sampling chem [%dx%dx%d] ...\n", dataW, dataH, dataD);
    sampleChem_kernel<<<grid, threads>>>(
        d_Src,
        dataD,
        dataH,
        dataW,
        chemIndex
    );
    getLastCudaError("sampleChem_kernel<<<>>> execution failed\n");
}
#endif	// ECV_SAMPLE_CHEM

#ifdef AVEP

#ifdef AVEP_INC

__global__
void bufferToVolumeAVEP_round0_kernel(
		float *d_Dst,
		float *d_Src,
		int svW,
		int svH,
		int svD,
		int volumeW,
		int volumeH,
		int volumeD,
		cudaPos offset,
		ecm_i ecmType)
{
	const int z = blockDim.z * blockIdx.z + threadIdx.z;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;
	const int x = blockDim.x * blockIdx.x + threadIdx.x;

	bool isInBound = (0 <= x && x < svW) &&
			   	 	 	 	   (0 <= y && y < svH) &&
			   	 	 	 	   (0 <= z && z < svD);


	if (isInBound)
	{

		int vx = x + (offset.x)/sizeof(VolumeType);
		int vy = y + offset.y;
		int vz = z + offset.z;

		float sample_dst = 0.0f;


		switch (ecmType)
		{
		case m_col:
			surf3Dread(&sample_dst, srfCol, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		case m_ela:
			surf3Dread(&sample_dst, srfEla, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		case m_hya:
			surf3Dread(&sample_dst, srfHya, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		default:
			surf3Dread(&sample_dst, srfCol, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		}

		// write diff to source
		d_Src[vz*volumeW*volumeH + vy*volumeW + vx] -= sample_dst;
		if (d_Src[vz*volumeW*volumeH + vy*volumeW + vx] == 0.0f) d_Src[vz*volumeW*volumeH + vy*volumeW + vx] = 0.4f;
	}

}

__global__
void bufferToVolumeAVEP_kernel(
		float *d_Dst,
		float *d_Src,
		int svW,
		int svH,
		int svD,
		int volumeW,
		int volumeH,
		int volumeD,
		cudaPos offset,
		ecm_i ecmType,
		int incRound,
		float incFactor)
{
	const int z = blockDim.z * blockIdx.z + threadIdx.z;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;
	const int x = blockDim.x * blockIdx.x + threadIdx.x;

	bool isInBound = (0 <= x && x < svW) &&
			   	 	 	 	   (0 <= y && y < svH) &&
			   	 	 	 	   (0 <= z && z < svD);


	if (isInBound)
	{

		int vx = x + (offset.x)/sizeof(VolumeType);
		int vy = y + offset.y;
		int vz = z + offset.z;

		float multiplier = ((float) (incRound + 1))*incFactor;
		float sample = 0.0f;

		switch (ecmType)
		{
		case m_col:
			surf3Dread(&sample, srfCol, vx * sizeof(VolumeType), vy, vz);
			sample += multiplier*d_Src[vz*volumeW*volumeH + vy*volumeW + vx];
			surf3Dwrite(sample, srfCol, vx * sizeof(VolumeType), vy, vz);
			break;
		case m_ela:
			surf3Dread(&sample, srfEla, vx * sizeof(VolumeType), vy, vz);
			sample += multiplier*d_Src[vz*volumeW*volumeH + vy*volumeW + vx];
			surf3Dwrite(sample, srfEla, vx * sizeof(VolumeType), vy, vz);
			break;
		case m_hya:
			surf3Dread(&sample, srfHya, vx * sizeof(VolumeType), vy, vz);
			sample += multiplier*d_Src[vz*volumeW*volumeH + vy*volumeW + vx];
			surf3Dwrite(sample, srfHya, vx * sizeof(VolumeType), vy, vz);
			break;
		default:
			surf3Dread(&sample, srfCol, vx * sizeof(VolumeType), vy, vz);
			sample += multiplier*d_Src[vz*volumeW*volumeH + vy*volumeW + vx];
			surf3Dwrite(sample, srfCol, vx * sizeof(VolumeType), vy, vz);
			break;
		}
	}

}

extern "C"
void bufferToVolumeAVEP(
		float *d_Dst,
		float *d_Src,
		int svW,
		int svH,
		int svD,
		int volumeW,
		int volumeH,
		int volumeD,
		cudaPos offset,
		ecm_i ecmType,
		int incRound,
		float incFactor)
{
  assert(d_Src != d_Dst);
  assert(svW <= volumeW);
  assert(svH <= volumeH);
  assert(svD <= volumeD);
  dim3 threads(8, 8, 4);
  dim3 grid(iDivUp_AVEP(svW, threads.x), iDivUp_AVEP(svH, threads.y), iDivUp_AVEP(svD, threads.z));

  if (!incRound)	// round 0
  {
    bufferToVolumeAVEP_round0_kernel<<<grid, threads>>>(
        d_Dst,
        d_Src,
        svW,
        svH,
        svD,
        volumeW,
        volumeH,
        volumeD,
        offset,
        ecmType
    );
    bufferToVolumeAVEP_kernel<<<grid, threads>>>(
        d_Dst,
        d_Src,
        svW,
        svH,
        svD,
        volumeW,
        volumeH,
        volumeD,
        offset,
        ecmType,
        incRound,
        incFactor
    );
  } else {
    bufferToVolumeAVEP_kernel<<<grid, threads>>>(
        d_Dst,
        d_Src,
        svW,
        svH,
        svD,
        volumeW,
        volumeH,
        volumeD,
        offset,
        ecmType,
        incRound,
        incFactor
    );
  }

  getLastCudaError("bufferToVolumeAVEP_kernel<<<>>> execution failed\n");
}


#else	// AVEP_INC


__global__
void bufferToVolumeAVEP_kernel(
		float *d_Dst,
		float *d_Src,
		int svW,
		int svH,
		int svD,
		int volumeW,
		int volumeH,
		int volumeD,
		cudaPos offset,
		ecm_i ecmType)
{
	const int z = blockDim.z * blockIdx.z + threadIdx.z;
	const int y = blockDim.y * blockIdx.y + threadIdx.y;
	const int x = blockDim.x * blockIdx.x + threadIdx.x;

	bool isInBound = (0 <= x && x < svW) &&
			   	 	 	 	   (0 <= y && y < svH) &&
			   	 	 	 	   (0 <= z && z < svD);


	if (isInBound)
	{

		int vx = x + (offset.x)/sizeof(VolumeType);
		int vy = y + offset.y;
		int vz = z + offset.z;

		float sample = d_Src[z*SV_W*SV_H + y*SV_W + x];	// use buffer dimensions


		switch (ecmType)
		{
		case m_col:
			surf3Dwrite(sample, srfCol, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		case m_ela:
			surf3Dwrite(sample, srfEla, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		case m_hya:
			surf3Dwrite(sample, srfHya, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		default:
			surf3Dwrite(sample, srfCol, vx * sizeof(VolumeType), vy, vz);//, cudaBoundaryModeZero);
			break;
		}
	}

}

extern "C"
void bufferToVolumeAVEP(
		float *d_Dst,
		float *d_Src,
		int svW,
		int svH,
		int svD,
		int volumeW,
		int volumeH,
		int volumeD,
		cudaPos offset,
		ecm_i ecmType)
{
  assert(d_Src != d_Dst);
  assert(svW <= volumeW);
  assert(svH <= volumeH);
  assert(svD <= volumeD);
  dim3 threads(8, 8, 4);
  dim3 grid(iDivUp_AVEP(svW, threads.x), iDivUp_AVEP(svH, threads.y), iDivUp_AVEP(svD, threads.z));

  bufferToVolumeAVEP_kernel<<<grid, threads>>>(
      d_Dst,
      d_Src,
      svW,
      svH,
      svD,
      volumeW,
      volumeH,
      volumeD,
      offset,
      ecmType
  );
  getLastCudaError("bufferToVolumeAVEP_kernel<<<>>> execution failed\n");
}

#endif	// AVEP_INC
#endif

__device__
int intersectBox(Ray r, float3 boxmin, float3 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float3 invR = make_float3(1.0f) / r.d;
    float3 tbot = invR * (boxmin - r.o);
    float3 ttop = invR * (boxmax - r.o);

    // re-order intersections to find smallest and largest on each axis
    float3 tmin = fminf(ttop, tbot);
    float3 tmax = fmaxf(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = fmaxf(fmaxf(tmin.x, tmin.y), fmaxf(tmin.x, tmin.z));
    float smallest_tmax = fminf(fminf(tmax.x, tmax.y), fminf(tmax.x, tmax.z));

    *tnear = largest_tmin;
    *tfar = smallest_tmax;

    return smallest_tmax > largest_tmin;
}

// transform vector by matrix (no translation)
__device__
float3 mul(const float3x4 &M, const float3 &v)
{
    float3 r;
    r.x = dot(v, make_float3(M.m[0]));
    r.y = dot(v, make_float3(M.m[1]));
    r.z = dot(v, make_float3(M.m[2]));
    return r;
}

// transform vector by matrix with translation
__device__
float4 mul(const float3x4 &M, const float4 &v)
{
    float4 r;
    r.x = dot(v, M.m[0]);
    r.y = dot(v, M.m[1]);
    r.z = dot(v, M.m[2]);
    r.w = 1.0f;
    return r;
}

__device__ uint rgbaFloatToInt(float4 rgba)
{
    rgba.x = __saturatef(rgba.x);   // clamp to [0.0, 1.0]
    rgba.y = __saturatef(rgba.y);
    rgba.z = __saturatef(rgba.z);
    rgba.w = __saturatef(rgba.w);
    return (uint(rgba.w*255)<<24) | (uint(rgba.z*255)<<16) | (uint(rgba.y*255)<<8) | uint(rgba.x*255);
}

__global__ void
d_render(uint *d_output, uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale)
{
    const int maxSteps = 500;
    const float tstep = 0.01f;
    const float opacityThreshold = 0.95f;

    const float3 boxMin = make_float3(-1.0f, -1.0f, -1.0f);
    const float3 boxMax = make_float3(1.0f, 1.0f, 1.0f);

    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;

    if ((x >= imageW) || (y >= imageH)) return;

    float u = (x / (float) imageW)*2.0f-1.0f;
    float v = (y / (float) imageH)*2.0f-1.0f;

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
    eyeRay.d = normalize(make_float3(u, v, -2.0f));
    eyeRay.d = mul(c_invViewMatrix, eyeRay.d);



    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

    if (!hit) return;

    if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating color
    float4 sum = make_float4(0.0f);
    float t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;


    for (int i=0; i<maxSteps; i++)
    {
        // read from 3D texture
        // remap position to [0, 1]
        float sample = tex3D(texCol, pos.x*0.5f+0.5f, pos.y*0.5f+0.5f, pos.z*0.5f+0.5f);

        //sample *= 64.0f;    // scale for 10-bit data

        // lookup in transfer function texture
        float4 col = tex1D(transferTexCol, (sample-transferOffset)*transferScale);

        col.w *= density;

        // "under" operator for back-to-front blending
        //sum = lerp(sum, col, col.w);

        // pre-multiply alpha
        col.x *= col.w;
        col.y *= col.w;
        col.z *= col.w;
        // "over" operator for front-to-back blending
        sum = sum + col*(1.0f - sum.w);

        // exit early if opaque
        if (sum.w > opacityThreshold)
            break;

        t += tstep;

        if (t > tfar) break;

        pos += step;

    }

    sum *= brightness;

    // write output color
    d_output[y*imageW + x] = rgbaFloatToInt(sum);
}

#ifdef ECV_SAMPLE_CHEM_TEST
__global__ void
d_render_test_dim(uint *d_output, uint nx, uint ny, uint nz, uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale,
         int chemType,
         bool isHighRes)
{
    const int maxSteps = 500;
    const float tstep = 0.01f;
    const float opacityThreshold = 0.95f;

    // Calculate box dimensions using largest dimension as reference
    const float a = -1.0f;
    const float b = +1.0f;
    const float ref = (float) max(nx, max(ny, nz));

    const float x_halfwidth = (((float) nx)/(2.0f * ref))*(b-a);
    const float y_halfwidth = (((float) ny)/(2.0f * ref))*(b-a);
    const float z_halfwidth = (((float) nz)/(2.0f * ref))*(b-a);



    const float3 boxMin = make_float3(-1.0f*x_halfwidth, -1.0f*y_halfwidth, -1.0f*z_halfwidth);
    const float3 boxMax = make_float3( 1.0f*x_halfwidth,  1.0f*y_halfwidth,  1.0f*z_halfwidth);


    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;

    if ((x >= imageW) || (y >= imageH)) return;

    // Calculate ray vector direction
    float u = ((float) x / (float) imageW)*2.0f-1.0f;
    float v = ((float) y / (float) imageH)*2.0f-1.0f;

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
    eyeRay.d = normalize(make_float3(u, v, -2.0f));
    eyeRay.d = mul(c_invViewMatrix, eyeRay.d);


    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

    if (!hit) return;

    if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating color
    float4 sum = make_float4(0.0f);
    float t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;


    for (int i=0; i<maxSteps; i++)
    {
        // read from 3D texture
        // remap position to [0, 1]
  	  float posx = (pos.x + x_halfwidth)/(2.0f*x_halfwidth);
  	  float posy = (pos.y + y_halfwidth)/(2.0f*y_halfwidth);
  	  float posz = (pos.z + z_halfwidth)/(2.0f*z_halfwidth);


  	  float sample_chem;
  	  float4 col_chem;

  	  if (isHighRes)
  	  {
  	  	sample_chem = tex3D(texChem0, posx, posy, posz);
  	  } else {	// low resolution
  	  	sample_chem = tex3D(texChem2, posx, posy, posz);
  	  }

  	  // lookup in transfer function texture
 	  	switch (chemType)
  	  	{
  	  	case 0:
  	  		col_chem = tex1D(transferTexChem0, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	case 1:
  	  		col_chem = tex1D(transferTexChem1, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	case 2:
  	  		col_chem = tex1D(transferTexChem2, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	case 3:
  	  		col_chem = tex1D(transferTexChem3, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	case 4:
  	  		col_chem = tex1D(transferTexChem4, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	case 5:
  	  		col_chem = tex1D(transferTexChem5, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	case 6:
  	  		col_chem = tex1D(transferTexChem6, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	case 7:
  	  		col_chem = tex1D(transferTexChem7, (sample_chem-transferOffset)*transferScale);
  	  		break;
  	  	}


      col_chem.w *= density;

      // "under" operator for back-to-front blending
      //sum = lerp(sum, col, col.w);

      // pre-multiply alpha
      col_chem.x *= col_chem.w;
      col_chem.y *= col_chem.w;
      col_chem.z *= col_chem.w;
      // "over" operator for front-to-back blending
      sum = sum + col_chem*(1.0f - sum.w);

      // exit early if opaque
      if (sum.w > opacityThreshold)
          break;

      t += tstep;

      if (t > tfar) break;

      pos += step;
    }

    sum *= brightness;

    // write output color
    d_output[y*imageW + x] = rgbaFloatToInt(sum);
}

#endif	// ECV_SAMPLE_CHEM_TEST

__global__ void
d_render_sp_dim(uint *d_output, uint nx, uint ny, uint nz, uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale, int gpu_id)
{
    const int maxSteps = 500;
    const float tstep = 0.01f;
    const float opacityThreshold = 0.95f;

    // Calculate box dimensions using largest dimension as reference
    const float a = -1.0f;
    const float b = +1.0f;
    const float ref = (float) max(nx, max(ny, nz));

    const float x_halfwidth = (((float) nx)/(2.0f * ref))*(b-a);
    const float y_halfwidth = (((float) ny)/(2.0f * ref))*(b-a);
    const float z_halfwidth = (((float) nz)/(2.0f * ref))*(b-a);

//    const float x_halfwidth_chem = x_halfwidth/ECV_SAMPLE_STRIDE;
//    const float y_halfwidth_chem = y_halfwidth/ECV_SAMPLE_STRIDE;
//    const float z_halfwidth_chem = z_halfwidth/ECV_SAMPLE_STRIDE;


    const float3 boxMin = make_float3(-1.0f*x_halfwidth, -1.0f*y_halfwidth, -1.0f*z_halfwidth);
    const float3 boxMax = make_float3( 1.0f*x_halfwidth,  1.0f*y_halfwidth,  1.0f*z_halfwidth);

//    const float3 boxMin_chem = make_float3(-1.0f*x_halfwidth_chem, -1.0f*y_halfwidth_chem, -1.0f*z_halfwidth_chem);
//    const float3 boxMax_chem = make_float3( 1.0f*x_halfwidth_chem,  1.0f*y_halfwidth_chem,  1.0f*z_halfwidth_chem);

    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;

    if ((x >= imageW) || (y >= imageH)) return;

    // Calculate ray vector direction using largest dimension as reference
#ifdef 	ECV_SEPARATE
    float u = ((float) x / (float) imageW)*2.0f-1.0f;
    float v = ((float) y / (float) imageH)*2.0f-1.0f;
#else	// ECV_SEPARATE
    const float ray_ref = (float) max(imageW, imageH);
    float u = ((float) x / ray_ref)*2.0f - (imageW/ray_ref);//(float) imageW)*2.0f-1.0f;
    float v = ((float) y / ray_ref)*2.0f - (imageH/ray_ref);//(float) imageH)*2.0f-1.0f;
#endif	// ECV_SEPARATE

    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
    eyeRay.d = normalize(make_float3(u, v, -2.0f));
    eyeRay.d = mul(c_invViewMatrix, eyeRay.d);


    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

    if (!hit) return;

    if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating color
    float4 sum = make_float4(0.0f);
    float t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;


    for (int i=0; i<maxSteps; i++)
    {
        // read from 3D texture
        // remap position to [0, 1]
  	  float posx = (pos.x + x_halfwidth)/(2.0f*x_halfwidth);
  	  float posy = (pos.y + y_halfwidth)/(2.0f*y_halfwidth);
  	  float posz = (pos.z + z_halfwidth)/(2.0f*z_halfwidth);

  	  float posx_chem = posx;//(pos.x + x_halfwidth_chem)/(2.0f*x_halfwidth_chem);
  	  float posy_chem = posy;//(pos.y + y_halfwidth_chem)/(2.0f*y_halfwidth_chem);
  	  float posz_chem = posz;//(pos.z + z_halfwidth_chem)/(2.0f*z_halfwidth_chem);

  	  float sample_chem0, sample_chem1;
  	  float sample_chem2, sample_chem3;
  	  float4 col_chem;
  	  float4 col_chem0, col_chem1;
  	  float4 col_chem2, col_chem3;

//  	  float sample;
//  	  float4 col;

    	// Assuming 2 GPUs
    	// TODO: GPUs > 2
  	  float blendFactor = 0.5f;
    	if (gpu_id == 0) {
  	  	sample_chem0 = 10.0f * tex3D(texChem0, posx_chem, posy_chem, posz_chem);
  	  	sample_chem1 = 1.0f * tex3D(texChem2, posx_chem, posy_chem, posz_chem);
  	  	sample_chem2 = 1.0f * tex3D(texChem4, posx_chem, posy_chem, posz_chem);
  	  	sample_chem3 = 1.0f * tex3D(texChem6, posx_chem, posy_chem, posz_chem);
  	  	// lookup in transfer function texture
  	  	col_chem0 = tex1D(transferTexChem0, (sample_chem0-transferOffset)*transferScale);
  	  	col_chem1 = tex1D(transferTexChem2, (sample_chem1-transferOffset)*transferScale);
  	  	col_chem2 = tex1D(transferTexChem4, (sample_chem2-transferOffset)*transferScale);
  	  	col_chem3 = tex1D(transferTexChem6, (sample_chem3-transferOffset)*transferScale);
  	  	// blend
  	  	col_chem.x = 0.25f*col_chem0.x + 0.25f*col_chem1.x + 0.25f*col_chem2.x + 0.25f*col_chem3.x;
  	  	col_chem.y = 0.25f*col_chem0.y + 0.25f*col_chem1.y + 0.25f*col_chem2.y + 0.25f*col_chem3.y;
  	  	col_chem.z = 0.25f*col_chem0.z + 0.25f*col_chem1.z + 0.25f*col_chem2.z + 0.25f*col_chem3.z;
  	  	col_chem.w = 0.25f*col_chem0.w + 0.25f*col_chem1.w + 0.25f*col_chem2.w + 0.25f*col_chem3.w;


    	} else {
  	  	sample_chem0 = 1.0f * tex3D(texChem1, posx_chem, posy_chem, posz_chem);
  	  	sample_chem1 = 1.0f * tex3D(texChem3, posx_chem, posy_chem, posz_chem);
  	  	sample_chem2 = 1.0f * tex3D(texChem5, posx_chem, posy_chem, posz_chem);
  	  	sample_chem3 = 1.0f * tex3D(texChem7, posx_chem, posy_chem, posz_chem);
  	  	// lookup in transfer function texture
  	  	col_chem0 = tex1D(transferTexChem1, (sample_chem0-transferOffset)*transferScale);
  	  	col_chem1 = tex1D(transferTexChem3, (sample_chem1-transferOffset)*transferScale);
  	  	col_chem2 = tex1D(transferTexChem5, (sample_chem2-transferOffset)*transferScale);
  	  	col_chem3 = tex1D(transferTexChem7, (sample_chem3-transferOffset)*transferScale);
  	  	// blend
  	  	col_chem.x = 1.0f*col_chem0.x + 1.0f*col_chem1.x + 1.0f*col_chem2.x + 1.0f*col_chem3.x;
  	  	col_chem.y = 1.0f*col_chem0.y + 1.0f*col_chem1.y + 1.0f*col_chem2.y + 1.0f*col_chem3.y;
  	  	col_chem.z = 1.0f*col_chem0.z + 1.0f*col_chem1.z + 1.0f*col_chem2.z + 1.0f*col_chem3.z;
  	  	col_chem.w = 1.0f*col_chem0.w + 1.0f*col_chem1.w + 1.0f*col_chem2.w + 1.0f*col_chem3.w;

    	}

      col_chem.w *= density;

      // "under" operator for back-to-front blending
      //sum = lerp(sum, col, col.w);

      // pre-multiply alpha
      col_chem.x *= col_chem.w;
      col_chem.y *= col_chem.w;
      col_chem.z *= col_chem.w;
      // "over" operator for front-to-back blending
      sum = sum + col_chem*(1.0f - sum.w);

      // exit early if opaque
      if (sum.w > opacityThreshold)
          break;

      t += tstep;

      if (t > tfar) break;

      pos += step;
    }

    sum *= brightness;

    // write output color
    d_output[y*imageW + x] = rgbaFloatToInt(sum);
}


__global__ void
d_render_dim(uint *d_output, uint nx, uint ny, uint nz, uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale, int ecmChemType, bool isChem)
{
    const int maxSteps = 500;
    const float tstep = 0.01f;
    const float opacityThreshold = 0.95f;

    // Calculate box dimensions using largest dimension as reference
    const float a = -1.0f;
    const float b = +1.0f;
    const float ref = (float) max(nx, max(ny, nz));

    const float x_halfwidth = (((float) nx)/(2.0f * ref))*(b-a);
    const float y_halfwidth = (((float) ny)/(2.0f * ref))*(b-a);
    const float z_halfwidth = (((float) nz)/(2.0f * ref))*(b-a);

    const float3 boxMin = make_float3(-1.0f*x_halfwidth, -1.0f*y_halfwidth, -1.0f*z_halfwidth);
    const float3 boxMax = make_float3( 1.0f*x_halfwidth,  1.0f*y_halfwidth,  1.0f*z_halfwidth);

    uint x = blockIdx.x*blockDim.x + threadIdx.x;
    uint y = blockIdx.y*blockDim.y + threadIdx.y;

    if ((x >= imageW) || (y >= imageH)) return;

    // Calculate ray vector direction using largest dimension as reference
#ifdef 	ECV_SEPARATE
    float u = ((float) x / (float) imageW)*2.0f-1.0f;
    float v = ((float) y / (float) imageH)*2.0f-1.0f;
#else	// ECV_SEPARATE
    const float ray_ref = (float) max(imageW, imageH);
    float u = ((float) x / ray_ref)*2.0f - (imageW/ray_ref);//(float) imageW)*2.0f-1.0f;
    float v = ((float) y / ray_ref)*2.0f - (imageH/ray_ref);//(float) imageH)*2.0f-1.0f;
#endif	// ECV_SEPARATE



    // calculate eye ray in world space
    Ray eyeRay;
    eyeRay.o = make_float3(mul(c_invViewMatrix, make_float4(0.0f, 0.0f, 0.0f, 1.0f)));
    eyeRay.d = normalize(make_float3(u, v, -2.0f));
    eyeRay.d = mul(c_invViewMatrix, eyeRay.d);


    // find intersection with box
    float tnear, tfar;
    int hit = intersectBox(eyeRay, boxMin, boxMax, &tnear, &tfar);

    if (!hit) return;

    if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from front to back, accumulating color
    float4 sum = make_float4(0.0f);
    float t = tnear;
    float3 pos = eyeRay.o + eyeRay.d*tnear;
    float3 step = eyeRay.d*tstep;


    for (int i=0; i<maxSteps; i++)
    {
        // read from 3D texture
        // remap position to [0, 1]

    	  float posx = (pos.x + x_halfwidth)/(2.0f*x_halfwidth);
    	  float posy = (pos.y + y_halfwidth)/(2.0f*y_halfwidth);
    	  float posz = (pos.z + z_halfwidth)/(2.0f*z_halfwidth);

#ifdef ECV_SAMPLE_CHEM
#ifdef ECV_INTERLEAVE
//    	  float posx_chem = (pos.x + x_halfwidth_chem)/(2.0f*x_halfwidth_chem);
//    	  float posy_chem = (pos.y + y_halfwidth_chem)/(2.0f*y_halfwidth_chem);
//    	  float posz_chem = (pos.z + z_halfwidth_chem)/(2.0f*z_halfwidth_chem);

    	  float sample_chem0, sample_chem2;
    	  float sample_chem4, sample_chem6;
    	  float4 col_chem;
    	  float4 col_chem0, col_chem2;
    	  float4 col_chem4, col_chem6;
#endif	// ECV_INTERLEAVE
#endif	// ECV_SAMPLE_CHEM

    	  float sample;
    	  float4 col;

    	  if (isChem)
    	  {
#ifdef ECV_SAMPLE_CHEM
    	  	switch (ecmChemType)
    	  	{
    	  	case 0:
      	  	sample = tex3D(texChem0, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem0, (sample-transferOffset)*transferScale);
    	  		break;
    	  	case 1:
      	  	sample = tex3D(texChem1, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem1, (sample-transferOffset)*transferScale);
    	  		break;
    	  	case 2:
      	  	sample = tex3D(texChem2, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem2, (sample-transferOffset)*transferScale);
    	  		break;
    	  	case 3:
      	  	sample = tex3D(texChem3, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem3, (sample-transferOffset)*transferScale);
    	  		break;
    	  	case 4:
      	  	sample = tex3D(texChem4, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem4, (sample-transferOffset)*transferScale);
    	  		break;
    	  	case 5:
      	  	sample = tex3D(texChem5, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem5, (sample-transferOffset)*transferScale);
    	  		break;
    	  	case 6:
      	  	sample = tex3D(texChem6, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem6, (sample-transferOffset)*transferScale);
    	  		break;
    	  	case 7:
      	  	sample = tex3D(texChem7, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexChem7, (sample-transferOffset)*transferScale);
    	  		break;
    	  	}

#endif	// ECV_SAMPLE_CHEM
    	  } else {

      	  switch (ecmChemType)
      	  {
      	  case m_col:
      	  {
      	  	sample = tex3D(texCol, posx, posy, posz);
      	  	//sample *= 64.0f;    // scale for 10-bit data

      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexCol, (sample-transferOffset)*transferScale);
#ifdef ECV_SAMPLE_CHEM
#ifdef ECV_INTERLEAVE
      	  	float ecmf = 0.6f;
      	  	float chmf = 1.0f-ecmf;
      	  	sample_chem0 = 10.0f * tex3D(texChem0, posx, posy, posz);
//      	  	sample_chem2 = 1.0f * tex3D(texChem2, posx_chem, posy_chem, posz_chem);
//      	  	sample_chem4 = 1.0f * tex3D(texChem4, posx_chem, posy_chem, posz_chem);
//      	  	sample_chem6 = 1.0f * tex3D(texChem6, posx_chem, posy_chem, posz_chem);
      	  	// lookup in transfer function texture
      	  	col_chem0 = tex1D(transferTexChem0, (sample_chem0-transferOffset)*transferScale);
//      	  	col_chem2 = tex1D(transferTexChem2, (sample_chem2-transferOffset)*transferScale);
//      	  	col_chem4 = tex1D(transferTexChem4, (sample_chem4-transferOffset)*transferScale);
//      	  	col_chem6 = tex1D(transferTexChem6, (sample_chem6-transferOffset)*transferScale);
      	  	// blend
//      	  	col_chem.x = 0.25f*col_chem0.x + 0.25f*col_chem2.x + 0.25f*col_chem4.x + 0.25f*col_chem6.x;
//      	  	col_chem.y = 0.25f*col_chem0.y + 0.25f*col_chem2.y + 0.25f*col_chem4.y + 0.25f*col_chem6.y;
//      	  	col_chem.z = 0.25f*col_chem0.z + 0.25f*col_chem2.z + 0.25f*col_chem4.z + 0.25f*col_chem6.z;
//      	  	col_chem.w = 0.25f*col_chem0.w + 0.25f*col_chem2.w + 0.25f*col_chem4.w + 0.25f*col_chem6.w;
      	  	col_chem.x = 1.0f*col_chem0.x;// + 1.0f*col_chem2.x + 1.0f*col_chem4.x + 1.0f*col_chem6.x;
      	  	col_chem.y = 1.0f*col_chem0.y;// + 1.0f*col_chem2.y + 1.0f*col_chem4.y + 1.0f*col_chem6.y;
      	  	col_chem.z = 1.0f*col_chem0.z;// + 1.0f*col_chem2.z + 1.0f*col_chem4.z + 1.0f*col_chem6.z;
      	  	col_chem.w = 1.0f*col_chem0.w;// + 1.0f*col_chem2.w + 1.0f*col_chem4.w + 1.0f*col_chem6.w;

//      	  	col.x = col_chem.x < 0.1f? chmf*col_chem.x + ecmf*col.x : chmf*col_chem.x + 1.0f*col.x;
//      	  	col.y = col_chem.y < 0.1f? chmf*col_chem.y + ecmf*col.y : chmf*col_chem.y + 1.0f*col.y;
//      	  	col.z = col_chem.z < 0.1f? chmf*col_chem.z + ecmf*col.z : chmf*col_chem.z + 1.0f*col.z;
//      	  	col.w = col_chem.w < 0.1f? chmf*col_chem.w + ecmf*col.w : chmf*col_chem.w + 1.0f*col.w;
      	  	col.x = chmf*col_chem.x + 1.0f*col.x;
      	  	col.y = chmf*col_chem.y + 1.0f*col.y;
      	  	col.z = chmf*col_chem.z + 1.0f*col.z;
      	  	col.w = chmf*col_chem.w + 1.0f*col.w;
#endif	// ECV_INTERLEAVE
#endif	// ECV_SAMPLE_CHEM

      	  	break;
      	  }
      	  case m_ela:
      	  	sample = tex3D(texEla, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexEla, (sample-transferOffset)*transferScale);
      	  	break;
      	  case m_hya:
      	  	sample = tex3D(texHya, posx, posy, posz);
      	  	// lookup in transfer function texture
      	  	col = tex1D(transferTexHya, (sample-transferOffset)*transferScale);
      	  	break;
      	  }
    	  }

        col.w *= density;

        // "under" operator for back-to-front blending
        //sum = lerp(sum, col, col.w);

        // pre-multiply alpha
        col.x *= col.w;
        col.y *= col.w;
        col.z *= col.w;
        // "over" operator for front-to-back blending
        sum = sum + col*(1.0f - sum.w);

        // exit early if opaque
        if (sum.w > opacityThreshold)
            break;

        t += tstep;

        if (t > tfar) break;

        pos += step;

    }

    sum *= brightness;

    // write output color
    d_output[y*imageW + x] = rgbaFloatToInt(sum);
}

extern "C"
void setTextureFilterMode(bool bLinearFilter)
{
    texCol.filterMode = bLinearFilter ? cudaFilterModeLinear : cudaFilterModePoint;
    texEla.filterMode = bLinearFilter ? cudaFilterModeLinear : cudaFilterModePoint;
    texHya.filterMode = bLinearFilter ? cudaFilterModeLinear : cudaFilterModePoint;
}

extern "C"
void printCpyParams(cudaMemcpy3DParms cp){
	/**
	 *   struct cudaArray     *srcArray;
  struct cudaPos        srcPos;
  struct cudaPitchedPtr srcPtr;
  struct cudaArray     *dstArray;
  struct cudaPos        dstPos;
  struct cudaPitchedPtr dstPtr;
  struct cudaExtent     extent;
  enum cudaMemcpyKind   kind;
	 */
	printf("copy params:\n");
	printf("\tsrcArray: %p\n", cp.srcArray);
	printf("\tsrcPos: %d, %d, %d\n", cp.srcPos.x, cp.srcPos.y, cp.srcPos.z);
	printf("\tsrcPtr:\n");
//	if(cp.srcPtr != 0)
//	{
		printf("\t\tpitch: %d\n", cp.srcPtr.pitch);
		printf("\t\tptr: %p\n",   cp.srcPtr.ptr);
		printf("\t\txsize: %d\n", cp.srcPtr.xsize);
		printf("\t\tysize: %d\n", cp.srcPtr.ysize);
//	}
	printf("\tdstArray: %p\n", cp.dstArray);
	printf("\tdstPos: %d, %d, %d\n", cp.dstPos.x, cp.dstPos.y, cp.dstPos.z);
	printf("\tdstPtr:\n");
//	if(cp.dstPtr != 0)
//	{
		printf("\t\tpitch: %d\n", cp.dstPtr.pitch);
		printf("\t\tptr: %p\n",   cp.dstPtr.ptr);
		printf("\t\txsize: %d\n", cp.dstPtr.xsize);
		printf("\t\tysize: %d\n", cp.dstPtr.ysize);
//	}
	printf("\textent: %d, %d, %d\n", cp.extent.width, cp.extent.height, cp.extent.depth);
}


#ifdef AVEP
#ifdef AVEP_INC

extern "C"
void bufferECMmapAVEP(
		cudaMemcpy3DParms copyParams,
		cudaMemcpy3DParms svCopyParams,
		ecm_i ecmType,
		int incRound,
		float incFactor)
{
	if (!incRound)	// first round
	{
		// copy a subvolume from host into device buffer
		printf("\t\tcopying...\n");
		checkCudaErrors(cudaMemcpy3D(&svCopyParams));
		printf("\t\tdone copying\n");
	}

	// copy data from device buffer to device volume array
	printf("\t\tbuffering...\n");
	VolumeType *d_Src = (VolumeType *) svCopyParams.dstPtr.ptr;
	VolumeType *d_Dst = (VolumeType *) copyParams.dstArray;
	bufferToVolumeAVEP(
			d_Dst,
			d_Src,
			(svCopyParams.extent.width)/sizeof(VolumeType),
			svCopyParams.extent.height,
			svCopyParams.extent.depth,
			copyParams.extent.width,
			copyParams.extent.height,
			copyParams.extent.depth,
			svCopyParams.srcPos,
			ecmType,
			incRound,
			incFactor);
	printf("\t\tdone buffering\n");

}

#else	// AVEP_INC

extern "C"
void bufferECMmapAVEP(
		cudaMemcpy3DParms copyParams,
		cudaMemcpy3DParms svCopyParams,
		ecm_i ecmType)
{
	// copy a subvolume from host into device buffer
	printf("\t\tcopying...\n");
	checkCudaErrors(cudaMemcpy3D(&svCopyParams));
	printf("\t\tdone copying\n");

	// copy data from device buffer to device volume array
	printf("\t\tbuffering...\n");
	VolumeType *d_Src = (VolumeType *) svCopyParams.dstPtr.ptr;
	VolumeType *d_Dst = (VolumeType *) copyParams.dstArray;
	bufferToVolumeAVEP(
			d_Dst,
			d_Src,
			(svCopyParams.extent.width)/sizeof(VolumeType),
			svCopyParams.extent.height,
			svCopyParams.extent.depth,
			copyParams.extent.width,
			copyParams.extent.height,
			copyParams.extent.depth,
			svCopyParams.srcPos,
			ecmType);
	printf("\t\tdone buffering\n");

}
#endif	// AVEP_INC
#endif	// AVEP

#ifdef ECV_SAMPLE_CHEM
#ifdef ECV_SAMPLE_CHEM_TEST
// gets called in DiffusionHelper.cpp
extern "C"
void initCudaChemSample(cudaExtent volumeSize, int chemIndex)
{
  // create 3D array
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  checkCudaErrors(cudaMalloc3DArray(&(d_chemsample_h[chemIndex]), &channelDesc, volumeSize, cudaArraySurfaceLoadStore));

  cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();

  switch (chemIndex)
  {
  case 0:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem0, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem0.normalized = true;                      // access with normalized texture coordinates
    texChem0.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem0.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem0.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem0, d_chemsample_h[chemIndex], channelDesc));

    // TNF:
  	float4 transferFunc[] =
  	{
  			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
  			{  0.23,  0.11,  0.44, 0.30, },	// purple
  			{  0.84,  0.42,  0.47, 0.40, }, // salmon pink
  			{  1.00,  0.69,  0.48, 0.80, }, // mild orange
  			{  1.00,  0.79,  0.58, 1.00, }, // light mild orange
  	};


    cudaArray *d_transferFuncArrayChem0;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem0, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem0, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem0.filterMode = cudaFilterModeLinear;
    transferTexChem0.normalized = true;    // access with normalized texture coordinates
    transferTexChem0.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem0, d_transferFuncArrayChem0, channelDesc2));

  	break;
  }
  case 2:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem2, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem2.normalized = true;                      // access with normalized texture coordinates
    texChem2.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem2.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem2.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem2, d_chemsample_h[chemIndex], channelDesc));

    // TGF: Purple-Turquoise
  	float4 transferFunc[] =
  	{
  			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
  			{  0.623520,  0.372549,  0.623529, 0.30, },	// 0.10	// purple
//  			{  1.00, 0.32, 0.18, 0.30, }, // bright orange
  			{  0.60,  0.80,  0.196078, 0.50, }, // 0.20	// yellow-green
  			{  1.00,  1.00,  0.00, 0.60, }, // 0.60 // yellow
  			{  0.196078,  0.60,  0.80, 0.80, }, // 0.80 // sky blue
  			{  0.439216,  0.858824,  0.576471, 1.00, }, // Turquoise
  	};


    cudaArray *d_transferFuncArrayChem2;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem2, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem2, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem2.filterMode = cudaFilterModeLinear;
    transferTexChem2.normalized = true;    // access with normalized texture coordinates
    transferTexChem2.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem2, d_transferFuncArrayChem2, channelDesc2));

  	break;
  }
  default:
  {
  	printf("Chem Resolution Comparison: Wrong buffer index %d\n", chemIndex);
  	exit(-1);
  }
  }
}

#else	// ECV_SAMPLE_CHEM_TEST
// gets called in DiffusionHelper.cpp
extern "C"
void initCudaChemSample(cudaExtent volumeSize, int chemIndex)
{
  // create 3D array
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  checkCudaErrors(cudaMalloc3DArray(&(d_chemsample_h[chemIndex]), &channelDesc, volumeSize, cudaArraySurfaceLoadStore));

  cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();

  switch (chemIndex)
  {
  case 0:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem0, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem0.normalized = true;                      // access with normalized texture coordinates
    texChem0.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem0.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem0.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem0, d_chemsample_h[chemIndex], channelDesc));

    // TNF:
  	float4 transferFunc[] =
  	{
  			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
  			{  0.23,  0.11,  0.44, 0.30, },	// purple
  			{  0.84,  0.42,  0.47, 0.40, }, // salmon pink
  			{  1.00,  0.69,  0.48, 0.80, }, // mild orange
  			{  1.00,  0.79,  0.58, 1.00, }, // light mild orange
  	};


    cudaArray *d_transferFuncArrayChem0;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem0, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem0, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem0.filterMode = cudaFilterModeLinear;
    transferTexChem0.normalized = true;    // access with normalized texture coordinates
    transferTexChem0.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem0, d_transferFuncArrayChem0, channelDesc2));

  	break;
  }
  case 1:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem1, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem1.normalized = true;                      // access with normalized texture coordinates
    texChem1.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem1.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem1.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem1, d_chemsample_h[chemIndex], channelDesc));

    // TGF: Purple-Turquoise
  	float4 transferFunc[] =
  	{
  			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
  			{  0.623520,  0.372549,  0.623529, 0.30, },	// 0.10	// purple
//  			{  1.00, 0.32, 0.18, 0.30, }, // bright orange
  			{  0.60,  0.80,  0.196078, 0.50, }, // 0.20	// yellow-green
  			{  1.00,  1.00,  0.00, 0.60, }, // 0.60 // yellow
  			{  0.196078,  0.60,  0.80, 0.80, }, // 0.80 // sky blue
  			{  0.439216,  0.858824,  0.576471, 1.00, }, // Turquoise
  	};


    cudaArray *d_transferFuncArrayChem1;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem1, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem1, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem1.filterMode = cudaFilterModeLinear;
    transferTexChem1.normalized = true;    // access with normalized texture coordinates
    transferTexChem1.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem1, d_transferFuncArrayChem1, channelDesc2));

  	break;
  }
  case 2:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem2, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem2.normalized = true;                      // access with normalized texture coordinates
    texChem2.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem2.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem2.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem2, d_chemsample_h[chemIndex], channelDesc));

    // FGF: Brown

  	float4 transferFunc[] =
  	{
  			{  0.000000,  0.000000,  0.000000, 0.00, },	// 0.00
  			{  0.647059,  0.164706,  0.164706, 0.05, },	// 0.05
  			{  0.647059,  0.164706,  0.164706, 0.10, },	// 0.10
  			{  0.647059,  0.164706,  0.164706, 0.15, },	// 0.15
  			{  0.647059,  0.164706,  0.164706, 0.20, }, // 0.20
  			{  0.647059,  0.164706,  0.164706, 0.25, }, // 0.25
  			{  0.647059,  0.164706,  0.164706, 0.30, }, // 0.30
  			{  0.647059,  0.164706,  0.164706, 0.35, }, // 0.35
  			{  0.647059,  0.164706,  0.164706, 0.40, }, // 0.40
  			{  0.647059,  0.164706,  0.164706, 0.45, }, // 0.45
  			{  0.647059,  0.164706,  0.164706, 0.50, }, // 0.50
  			{  0.647059,  0.164706,  0.164706, 0.55, }, // 0.55
  			{  0.647059,  0.164706,  0.164706, 0.60, }, // 0.60
  			{  0.647059,  0.164706,  0.164706, 0.65, }, // 0.65
  			{  0.647059,  0.164706,  0.164706, 0.70, }, // 0.70
  			{  0.647059,  0.164706,  0.164706, 0.75, }, // 0.75
  			{  0.647059,  0.164706,  0.164706, 0.80, }, // 0.80
  			{  0.647059,  0.164706,  0.164706, 0.85, }, // 0.85
  			{  0.647059,  0.164706,  0.164706, 0.90, }, // 0.90
  			{  0.647059,  0.164706,  0.164706, 0.95, }, // 0.95
  			{  0.647059,  0.164706,  0.164706, 1.00, }, // 1.00
  	};

    cudaArray *d_transferFuncArrayChem2;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem2, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem2, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem2.filterMode = cudaFilterModeLinear;
    transferTexChem2.normalized = true;    // access with normalized texture coordinates
    transferTexChem2.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem2, d_transferFuncArrayChem2, channelDesc2));

  	break;
  }
  case 3:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem3, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem3.normalized = true;                      // access with normalized texture coordinates
    texChem3.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem3.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem3.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem3, d_chemsample_h[chemIndex], channelDesc));

    // MMP8: Sky blue

  	float4 transferFunc[] =
  	{
  			{  0.000000,  0.000000,  0.000000, 0.00, },	// 0.00
  			{  0.196078,  0.600000,  0.800000, 0.05, },	// 0.05
  			{  0.196078,  0.600000,  0.800000, 0.10, },	// 0.10
  			{  0.196078,  0.600000,  0.800000, 0.15, },	// 0.15
  			{  0.196078,  0.600000,  0.800000, 0.20, }, // 0.20
  			{  0.196078,  0.600000,  0.800000, 0.25, }, // 0.25
  			{  0.196078,  0.600000,  0.800000, 0.30, }, // 0.30
  			{  0.196078,  0.600000,  0.800000, 0.35, }, // 0.35
  			{  0.196078,  0.600000,  0.800000, 0.40, }, // 0.40
  			{  0.196078,  0.600000,  0.800000, 0.45, }, // 0.45
  			{  0.196078,  0.600000,  0.800000, 0.50, }, // 0.50
  			{  0.196078,  0.600000,  0.800000, 0.55, }, // 0.55
  			{  0.196078,  0.600000,  0.800000, 0.60, }, // 0.60
  			{  0.196078,  0.600000,  0.800000, 0.65, }, // 0.65
  			{  0.196078,  0.600000,  0.800000, 0.70, }, // 0.70
  			{  0.196078,  0.600000,  0.800000, 0.75, }, // 0.75
  			{  0.196078,  0.600000,  0.800000, 0.80, }, // 0.80
  			{  0.196078,  0.600000,  0.800000, 0.85, }, // 0.85
  			{  0.196078,  0.600000,  0.800000, 0.90, }, // 0.90
  			{  0.196078,  0.600000,  0.800000, 0.95, }, // 0.95
  			{  0.196078,  0.600000,  0.800000, 1.00, }, // 1.00
  	};


    cudaArray *d_transferFuncArrayChem;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem3, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem3, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem3.filterMode = cudaFilterModeLinear;
    transferTexChem3.normalized = true;    // access with normalized texture coordinates
    transferTexChem3.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem3, d_transferFuncArrayChem3, channelDesc2));

  	break;
  }
  case 4:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem4, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem4.normalized = true;                      // access with normalized texture coordinates
    texChem4.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem4.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem4.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem4, d_chemsample_h[chemIndex], channelDesc));

    // IL1: Green Beach from https://digitalsynopsis.com/design/beautiful-color-ui-gradients-backgrounds/
   	float4 transferFunc[] =
    	{
    			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
    			{  0.01,  0.67,  0.69, 0.50, },	// blue-mild green
    			{  0.00,  0.80,  0.67, 1.00, }, // light blue-mild green
    	};


    cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
    cudaArray *d_transferFuncArrayChem4;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem4, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem4, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem4.filterMode = cudaFilterModeLinear;
    transferTexChem4.normalized = true;    // access with normalized texture coordinates
    transferTexChem4.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem4, d_transferFuncArrayChem4, channelDesc2));

  	break;
  }
  case 5:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem5, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem5.normalized = true;                      // access with normalized texture coordinates
    texChem5.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem5.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem5.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem5, d_chemsample_h[chemIndex], channelDesc));

    // IL6: Pink
  	float4 transferFunc[] =
  	{
  			{  0.000000,  0.000000,  0.000000, 0.00, },	// 0.00
  			{  0.737255,  0.560784,  0.560784, 0.05, },	// 0.05
  			{  0.737255,  0.560784,  0.560784, 0.10, },	// 0.10
  			{  0.737255,  0.560784,  0.560784, 0.15, },	// 0.15
  			{  0.737255,  0.560784,  0.560784, 0.20, }, // 0.20
  			{  0.737255,  0.560784,  0.560784, 0.25, }, // 0.25
  			{  0.737255,  0.560784,  0.560784, 0.30, }, // 0.30
  			{  0.737255,  0.560784,  0.560784, 0.35, }, // 0.35
  			{  0.737255,  0.560784,  0.560784, 0.40, }, // 0.40
  			{  0.737255,  0.560784,  0.560784, 0.45, }, // 0.45
  			{  0.737255,  0.560784,  0.560784, 0.50, }, // 0.50
  			{  0.737255,  0.560784,  0.560784, 0.55, }, // 0.55
  			{  0.737255,  0.560784,  0.560784, 0.60, }, // 0.60
  			{  0.737255,  0.560784,  0.560784, 0.65, }, // 0.65
  			{  0.737255,  0.560784,  0.560784, 0.70, }, // 0.70
  			{  0.737255,  0.560784,  0.560784, 0.75, }, // 0.75
  			{  0.737255,  0.560784,  0.560784, 0.80, }, // 0.80
  			{  0.737255,  0.560784,  0.560784, 0.85, }, // 0.85
  			{  0.737255,  0.560784,  0.560784, 0.90, }, // 0.90
  			{  0.737255,  0.560784,  0.560784, 0.95, }, // 0.95
  			{  0.737255,  0.560784,  0.560784, 1.00, }, // 1.00
  	};


    cudaArray *d_transferFuncArrayChem5;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem5, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem5, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem5.filterMode = cudaFilterModeLinear;
    transferTexChem5.normalized = true;    // access with normalized texture coordinates
    transferTexChem5.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem5, d_transferFuncArrayChem5, channelDesc2));

  	break;
  }
  case 6:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem6, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem6.normalized = true;                      // access with normalized texture coordinates
    texChem6.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem6.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem6.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem0, d_chemsample_h[chemIndex], channelDesc));

    // IL8: Orange-Yellow
  	float4 transferFunc[] =
  	{
  			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
  			{  0.23,  0.11,  0.44, 0.30, },	// purple
  			{  0.84,  0.42,  0.47, 0.40, }, // salmon pink
  			{  1.00,  0.69,  0.48, 0.80, }, // mild orange
  			{  1.00,  0.79,  0.58, 1.00, }, // light mild orange
  	};

    cudaArray *d_transferFuncArrayChem6;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem6, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem6, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem6.filterMode = cudaFilterModeLinear;
    transferTexChem6.normalized = true;    // access with normalized texture coordinates
    transferTexChem6.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem6, d_transferFuncArrayChem6, channelDesc2));

  	break;
  }
  case 7:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem7, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem7.normalized = true;                      // access with normalized texture coordinates
    texChem7.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem7.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem7.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem7, d_chemsample_h[chemIndex], channelDesc));

    // IL10: White
  	float4 transferFunc[] =
  	{
  			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
  			{  0.26,  0.81,  0.64, 0.30, },	// 0.10	// green
  			{  0.60,  0.80,  0.196078, 0.50, }, // 0.20	// yellow-green
  			{  1.00,  0.11,  0.68, 0.60, }, // 0.60
  			{  0.678431,  0.917647,  0.917647, 0.80, }, // 0.80
  			{  0.00,  0.00,  1.00, 1.00, }, // 1.00
  	};

    cudaArray *d_transferFuncArrayChem7;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem7, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem7, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem7.filterMode = cudaFilterModeLinear;
    transferTexChem7.normalized = true;    // access with normalized texture coordinates
    transferTexChem7.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem7, d_transferFuncArrayChem7, channelDesc2));

  	break;
  }
  default:
  {
  	// bind array to 3D surface
  	checkCudaErrors(cudaBindSurfaceToArray(srfChem0, d_chemsample_h[chemIndex], channelDesc));

    // set texture parameters
    texChem0.normalized = true;                      // access with normalized texture coordinates
    texChem0.filterMode = cudaFilterModeLinear;      // linear interpolation
    texChem0.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
    texChem0.addressMode[1] = cudaAddressModeClamp;

    // bind array to 3D texture
    checkCudaErrors(cudaBindTextureToArray(texChem0, d_chemsample_h[chemIndex], channelDesc));

  	float4 transferFunc[] =
  	{
  			{  0.00,  0.00,  0.00, 0.00, },	// 0.00
  			{  1.00,  0.00,  1.00, 0.05, },	// 0.05
  			{  1.00,  0.05,  0.90, 0.10, },	// 0.10
  			{  0.80,  0.10,  0.80, 0.15, },	// 0.15
  			{  0.60,  0.15,  0.70, 0.20, }, // 0.20
  			{  0.40,  0.20,  0.60, 0.25, }, // 0.25
  			{  0.20,  0.25,  0.50, 0.30, }, // 0.30
  			{  0.00,  0.30,  0.40, 0.35, }, // 0.35
  			{  0.40,  0.35,  0.30, 0.40, }, // 0.40
  			{  0.60,  0.40,  0.20, 0.45, }, // 0.45
  			{  0.70,  0.45,  0.10, 0.50, }, // 0.50
  			{  0.80,  0.45,  0.00, 0.55, }, // 0.55
  			{  0.90,  0.50,  0.00, 0.60, }, // 0.60
  			{  1.00,  0.50,  0.00, 0.65, }, // 0.65
  			{  1.00,  0.50,  0.00, 0.70, }, // 0.70
  			{  1.00,  0.50,  0.00, 0.75, }, // 0.75
  			{  1.00,  0.50,  0.00, 0.80, }, // 0.80
  			{  1.00,  0.50,  0.00, 0.85, }, // 0.85
  			{  1.00,  0.50,  0.00, 0.90, }, // 0.90
  			{  1.00,  0.50,  0.00, 0.95, }, // 0.95
  			{  1.00,  0.50,  0.00, 1.00, }, // 1.00
  	};

    cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
    cudaArray *d_transferFuncArrayChem0;
    checkCudaErrors(cudaMallocArray(&d_transferFuncArrayChem0, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
    checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayChem0, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

    transferTexChem0.filterMode = cudaFilterModeLinear;
    transferTexChem0.normalized = true;    // access with normalized texture coordinates
    transferTexChem0.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

    // Bind the array to the texture
    checkCudaErrors(cudaBindTextureToArray(transferTexChem0, d_transferFuncArrayChem0, channelDesc2));

  	break;
  }
  }

}
#endif	// ECV_SAMPLE_CHEM_TEST
#endif	// ECV_SAMPLE_CHEM

extern "C"
void bufferECMmap(cudaMemcpy3DParms copyParams)
{
	checkCudaErrors(cudaMemcpy3D(&copyParams));
}

extern "C"
void initCuda(void *h_volume, cudaExtent volumeSize, cudaMemcpy3DParms &copyParams, ecm_i ecmType)
{
    // create 3D array
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<VolumeType>();
    checkCudaErrors(cudaMalloc3DArray(&(d_volumeArray[ecmType]), &channelDesc, volumeSize, cudaArraySurfaceLoadStore));

    // copy data to 3D array
    copyParams.srcPtr   = make_cudaPitchedPtr(h_volume, volumeSize.width*sizeof(VolumeType), volumeSize.width, volumeSize.height);
    copyParams.dstArray = d_volumeArray[ecmType];
    copyParams.extent   = volumeSize;
    copyParams.kind     = cudaMemcpyHostToDevice;

    checkCudaErrors(cudaMemcpy3D(&copyParams));

    //create transfer function texture
    switch(ecmType)
		{
    case m_col:
    {
    	// Collagen
    	float4 transferFunc[] =
    	{
    			{  0.00,  0.00,  0.00, 0.0, },	// 0.00
    			{  1.00,  0.00,  0.00, 0.5, },	// 0.05 - SLP ILP
    			{  1.00,  0.30,  0.30, 0.8, },	// 0.10 - SLP ILP
    			{  0.00,  0.00,  0.00, 0.0, },	// 0.15 - SLP ILP
    			{  0.80,  0.15,  0.10, 0.2, }, // 0.20 - DLP ILP SLP
    			{  0.80,  0.15,  0.10, 0.3, }, // 0.25 - DLP
    			{  0.80,  0.15,  0.10, 0.4, }, // 0.30 - DLP
    			{  0.80,  0.15,  0.10, 0.5, }, // 0.35 - DLP
    			{  0.80,  0.15,  0.10, 0.6, }, // 0.40 - DLP
    			{  0.80,  0.15,  0.10, 0.7, }, // 0.45
    			{  0.80,  0.15,  0.10, 0.8, }, // 0.50
    			{  0.80,  0.15,  0.10, 0.9, }, // 0.55
    			{  0.85,  0.10,  0.15, 0.8, }, // 0.60
    			{  0.90,  0.10,  0.15, 0.7, }, // 0.65
    			{  0.95,  0.10,  0.10, 0.6, }, // 0.70
    			{  1.00,  0.10,  0.10, 0.5, }, // 0.75
    			{  1.00,  0.10,  0.10, 0.6, }, // 0.80
    			{  1.00,  0.10,  0.10, 0.7, }, // 0.85
    			{  1.00,  0.20,  0.20, 0.8, }, // 0.90
    			{  1.00,  0.30,  0.30, 0.9, }, // 0.95
    			{  1.00,  0.60,  0.00, 1.0, }, // 1.00
    			{  0.60,  0.40,  0.32, 1.0, },
    	};


#ifdef AVEP
    	// bind array to 3D surface
    	checkCudaErrors(cudaBindSurfaceToArray(srfCol, d_volumeArray[ecmType], channelDesc));
//    	checkCudaErrors(cudaBindSurfaceToArray(srfCol, d_volumeArray[ecmType]));
#endif
      // set texture parameters
      texCol.normalized = true;                      // access with normalized texture coordinates
      texCol.filterMode = cudaFilterModeLinear;      // linear interpolation
      texCol.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
      texCol.addressMode[1] = cudaAddressModeClamp;

      // bind array to 3D texture
      checkCudaErrors(cudaBindTextureToArray(texCol, d_volumeArray[ecmType], channelDesc));


      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
      cudaArray *d_transferFuncArrayCol;
      checkCudaErrors(cudaMallocArray(&d_transferFuncArrayCol, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
      checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayCol, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

      transferTexCol.filterMode = cudaFilterModeLinear;
      transferTexCol.normalized = true;    // access with normalized texture coordinates
      transferTexCol.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

      // Bind the array to the texture
      checkCudaErrors(cudaBindTextureToArray(transferTexCol, d_transferFuncArrayCol, channelDesc2));
    	break;
    }

    case m_ela:
    {
/      // Elastin
      float4 transferFunc[] =
      {
      		{  0.0,  0.0,  0.0, 0.0, },	// 0.00
      		{  0.0,  1.0,  0.0, 0.5, },	// 0.05 - SLP ILP
      		{  0.0,  0.0,  0.0, 0.0, },	// 0.10 - SLP ILP
      		{  0.63, 0.12,  0.4, 0.3, },	// 0.15 - SLP ILP
      		{  0.63, 0.12,  0.4, 0.4, }, // 0.20 - DLP ILP SLP
      		{  0.63, 0.12,  0.4, 0.5, }, // 0.25 - DLP
      		{  0.63, 0.12,  0.4, 0.6, }, // 0.30 - DLP
      		{  0.63, 0.12,  0.4, 0.7, }, // 0.35 - DLP
      		{  0.63, 0.12,  0.4, 0.8, }, // 0.40 - DLP
      		{  0.63, 0.12,  0.4, 0.9, }, // 0.45
      		{  0.0, 1.0,  0.30, 1.0, }, // 0.50
      		{  0.0, 1.0,  0.30, 1.0, }, // 0.55
      		{  0.0, 1.0,  0.30, 1.0, }, // 0.60
      		{  0.0, 1.0,  0.30, 1.0, }, // 0.65
      		{  0.0, 1.0,  0.30, 1.0, }, // 0.70
      		{  0.0, 1.0,  0.30, 1.0, }, // 0.75
      		{  0.0,  0.0,  0.0, 0.0, }, // 0.80
      		{  0.0, 1.0,  0.30, 1.0, }, // 0.85
      		{  0.0, 1.0,  0.40, 0.5, }, // 0.90
      		{  0.0, 1.0,  0.50, 0.7, }, // 0.95
      		{  0.0, 1.0,  0.60, 1.0, }, // 1.00
      };
      // set texture parameters
      texEla.normalized = true;                      // access with normalized texture coordinates
      texEla.filterMode = cudaFilterModeLinear;      // linear interpolation
      texEla.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
      texEla.addressMode[1] = cudaAddressModeClamp;

      // bind array to 3D texture
#ifdef AVEP
    	checkCudaErrors(cudaBindSurfaceToArray(srfEla, d_volumeArray[ecmType], channelDesc));
#endif
      checkCudaErrors(cudaBindTextureToArray(texEla, d_volumeArray[ecmType], channelDesc));

      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
      cudaArray *d_transferFuncArrayEla;
      checkCudaErrors(cudaMallocArray(&d_transferFuncArrayEla, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
      checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayEla, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

      transferTexEla.filterMode = cudaFilterModeLinear;
      transferTexEla.normalized = true;    // access with normalized texture coordinates
      transferTexEla.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

      // Bind the array to the texture
      checkCudaErrors(cudaBindTextureToArray(transferTexEla, d_transferFuncArrayEla, channelDesc2));

    	break;
    }

    case m_hya:
    {
      // Hyaluronan
#ifdef RAT_VF
      float4 transferFunc[] =
      {
      		{  0.0,  0.00,  0.00, 0.0, },	// 0.00
      		{  0.0,  0.00,  1.00, 0.5, },	// 0.05 - SLP ILP
      		{  0.3,  0.30,  1.00, 0.8, },	// 0.10 - SLP ILP
      		{  0.0,  0.00,  0.00, 0.0, },	// 0.15 - SLP ILP
      		{  1.0,  0.43,  0.78, 0.2, }, // 0.20 - DLP ILP SLP
      		{  1.0,  0.43,  0.78, 0.3, }, // 0.25 - DLP
      		{  1.0,  0.43,  0.78, 0.4, }, // 0.30 - DLP
      		{  1.0,  0.43,  0.78, 0.5, }, // 0.35 - DLP
      		{  1.0,  0.43,  0.78, 0.6, }, // 0.40 - DLP
      		{  1.0,  0.43,  0.78, 0.7, }, // 0.45
      		{  1.0,  0.43,  0.78, 0.8, }, // 0.50
      		{  1.0,  0.43,  0.78, 0.9, }, // 0.55
      		{  0.8,  0.33,  0.85, 0.8, }, // 0.60
      		{  0.5,  0.23,  0.90, 0.7, }, // 0.65
      		{  0.3,  0.13,  0.95, 0.6, }, // 0.70
      		{  0.0,  0.00,  1.00, 0.5, }, // 0.75
      		{  0.1,  0.10,  1.00, 0.6, }, // 0.80
      		{  0.2,  0.20,  1.00, 0.7, }, // 0.85
      		{  0.3,  0.30,  1.00, 0.8, }, // 0.90
      		{  0.4,  0.40,  1.00, 0.9, }, // 0.95
      		{  0.7,  0.70,  1.00, 1.0, }, // 1.00
      };
#else	// RAT_VF
      float4 transferFunc[] =
      {
      		{  0.0,  0.00,  0.00, 0.0, },	// 0.00
      		{  0.0,  0.00,  1.00, 0.5, },	// 0.05 - SLP ILP
      		{  0.3,  0.30,  1.00, 0.8, },	// 0.10 - SLP ILP
      		{  0.0,  0.00,  0.00, 0.0, },	// 0.15 - SLP ILP
      		{  0.1,  0.43,  0.78, 0.2, }, // 0.20 - DLP ILP SLP
      		{  0.2,  0.43,  0.78, 0.3, }, // 0.25 - DLP
      		{  0.3,  0.43,  0.78, 0.4, }, // 0.30 - DLP
      		{  0.4,  0.43,  0.78, 0.5, }, // 0.35 - DLP
      		{  0.5,  0.43,  0.78, 0.6, }, // 0.40 - DLP
      		{  0.6,  0.43,  0.78, 0.7, }, // 0.45
      		{  0.7,  0.43,  0.78, 0.8, }, // 0.50
      		{  0.7,  0.43,  0.78, 0.9, }, // 0.55
      		{  0.8,  0.33,  0.85, 0.8, }, // 0.60
      		{  0.5,  0.23,  0.90, 0.7, }, // 0.65
      		{  0.3,  0.13,  0.95, 0.6, }, // 0.70
      		{  0.0,  0.00,  1.00, 0.5, }, // 0.75
      		{  0.1,  0.10,  1.00, 0.6, }, // 0.80
      		{  0.2,  0.20,  1.00, 0.7, }, // 0.85
      		{  0.3,  0.30,  1.00, 0.8, }, // 0.90
      		{  0.4,  0.40,  1.00, 0.9, }, // 0.95
      		{  0.3,  0.10,  1.00, 1.0, }, // 1.00
      		{  0.0,  0.00,  1.00, 1.0, }, // 1.00
      };
#endif	// RAT_VF
      // set texture parameters
      texHya.normalized = true;                      // access with normalized texture coordinates
      texHya.filterMode = cudaFilterModeLinear;      // linear interpolation
      texHya.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
      texHya.addressMode[1] = cudaAddressModeClamp;

      // bind array to 3D texture
#ifdef AVEP
    	checkCudaErrors(cudaBindSurfaceToArray(srfHya, d_volumeArray[ecmType], channelDesc));
#endif
      checkCudaErrors(cudaBindTextureToArray(texHya, d_volumeArray[ecmType], channelDesc));

      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
      cudaArray *d_transferFuncArrayHya;
      checkCudaErrors(cudaMallocArray(&d_transferFuncArrayHya, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
      checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayHya, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

      transferTexHya.filterMode = cudaFilterModeLinear;
      transferTexHya.normalized = true;    // access with normalized texture coordinates
      transferTexHya.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

      // Bind the array to the texture
      checkCudaErrors(cudaBindTextureToArray(transferTexHya, d_transferFuncArrayHya, channelDesc2));


    	break;
    }

    default:
    {
    	// WEF
    	// White bg
    	float4 transferFunc[] =
    	{
    			{  0.0,  0.0,  0.0, 0.0, },
    			{  1.0,  0.0,  0.0, 1.0, },
    			{  1.0,  0.0,  0.0, 1.0, },
    			{  0.97, 0.8, 0.72, 1.0, },
    			{  0.97, 0.8, 0.72, 0.5, },
    			{  0.80, 0.6, 0.52, 0.7, },
    			{  0.60, 0.4, 0.32, 1.0, },//0.5, },
    	};

//    	// Black bg
//    	float4 transferFunc[] =
//    	{
//    			{  0.0,  0.0,  0.0, 0.0, },
//    			{  1.0,  0.0,  0.0, 1.0, },
//    			{  1.0,  0.0,  0.0, 1.0, },
//    			{  0.97, 0.8, 0.72, 1.0, },
//    			{  0.97, 0.4, 0.30, 1.0, },
//    			{  0.97, 0.6, 0.50, 0.7, },
//    			{  0.97, 0.8, 0.72, 0.8, },//0.5, },
//		  };
      // set texture parameters
      texCol.normalized = true;                      // access with normalized texture coordinates
      texCol.filterMode = cudaFilterModeLinear;      // linear interpolation
      texCol.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
      texCol.addressMode[1] = cudaAddressModeClamp;

      // bind array to 3D texture
#ifdef AVEP
    	checkCudaErrors(cudaBindSurfaceToArray(srfCol, d_volumeArray[ecmType], channelDesc));
#endif
      checkCudaErrors(cudaBindTextureToArray(texCol, d_volumeArray[ecmType], channelDesc));

      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
      cudaArray *d_transferFuncArrayCol;
      checkCudaErrors(cudaMallocArray(&d_transferFuncArrayCol, &channelDesc2, sizeof(transferFunc)/sizeof(float4), 1));
      checkCudaErrors(cudaMemcpyToArray(d_transferFuncArrayCol, 0, 0, transferFunc, sizeof(transferFunc), cudaMemcpyHostToDevice));

      transferTexCol.filterMode = cudaFilterModeLinear;
      transferTexCol.normalized = true;    // access with normalized texture coordinates
      transferTexCol.addressMode[0] = cudaAddressModeClamp;   // wrap texture coordinates

      // Bind the array to the texture
      checkCudaErrors(cudaBindTextureToArray(transferTexCol, d_transferFuncArrayCol, channelDesc2));

    	break;
    }

		}



}

#ifdef AVEP
extern "C"
void initCudaAVEP(
		void *h_volume,
		cudaExtent volumeSize,
		cudaMemcpy3DParms &copyParams,
		cudaMemcpy3DParms &svCopyParams,
		ecm_i ecmType)
{
	initCuda(h_volume, volumeSize, copyParams, ecmType);

#ifdef AVEP_INC
	// Allocate buffer device memory
	checkCudaErrors(cudaMalloc(&(d_svBuffer[ecmType]),
			volumeSize.width*volumeSize.height*volumeSize.depth*sizeof(VolumeType)));
	svCopyParams.dstPtr = make_cudaPitchedPtr(
																			d_svBuffer[ecmType],
																			volumeSize.width*sizeof(VolumeType),
																			volumeSize.width,
																			volumeSize.height);
#else	// AVEP_INC
	// Allocate buffer device memory
	checkCudaErrors(cudaMalloc(&(d_svBuffer[ecmType]), SV_W*SV_H*SV_D*sizeof(VolumeType)));
	svCopyParams.dstPtr = make_cudaPitchedPtr(
																			d_svBuffer[ecmType],
																			SV_W*sizeof(VolumeType),
																			SV_W,
																			SV_H);
#endif	// AVEP_INC

	// initialize copy params for sub-volumes
	svCopyParams.srcPtr  = copyParams.srcPtr;
	svCopyParams.extent = make_cudaExtent(SV_W*sizeof(VolumeType), SV_H, SV_D);
	svCopyParams.kind = cudaMemcpyHostToDevice;

}
#endif

extern "C"
void freeCudaBuffers()
{
	for(int ei = 0; ei < m_ecmtotal; ei++) {
    checkCudaErrors(cudaFreeArray(d_volumeArray[ei]));
	}

	checkCudaErrors(cudaFreeArray(d_transferFuncArrayCol));
	checkCudaErrors(cudaFreeArray(d_transferFuncArrayEla));
	checkCudaErrors(cudaFreeArray(d_transferFuncArrayHya));
}

extern "C"
void render_kernel(dim3 gridSize, dim3 blockSize, uint *d_output, uint imageW, uint imageH,
                   float density, float brightness, float transferOffset, float transferScale)
{
    d_render<<<gridSize, blockSize>>>(d_output, imageW, imageH, density,
                                      brightness, transferOffset, transferScale);
}

extern "C"
void render_kernel_dim(dim3 gridSize, dim3 blockSize, uint *d_output, uint nx, uint ny, uint nz, uint imageW, uint imageH,
                   float density, float brightness, float transferOffset, float transferScale, int ecmChemType, bool isChem)
{
    d_render_dim<<<gridSize, blockSize>>>(d_output, nx, ny, nz, imageW, imageH, density,
                                      brightness, transferOffset, transferScale, ecmChemType, isChem);
}

extern "C"
void render_sp_kernel_dim(dim3 gridSize, dim3 blockSize, uint *d_output, uint nx, uint ny, uint nz, uint imageW, uint imageH,
                   float density, float brightness, float transferOffset, float transferScale, int gpu_id)
{
    d_render_sp_dim<<<gridSize, blockSize>>>(d_output, nx, ny, nz, imageW, imageH, density,
                                      brightness, transferOffset, transferScale, gpu_id);
}

#ifdef ECV_SAMPLE_CHEM_TEST
extern "C"
void render_test_kernel_dim(dim3 gridSize, dim3 blockSize, uint *d_output, uint nx, uint ny, uint nz, uint imageW, uint imageH,
                   float density, float brightness, float transferOffset, float transferScale, int chemType, bool isHigh)
{
    d_render_test_dim<<<gridSize, blockSize>>>(d_output, nx, ny, nz, imageW, imageH, density,
                                      brightness, transferOffset, transferScale, chemType, isHigh);
}
#endif

extern "C"
void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix)
{
    checkCudaErrors(cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeofMatrix));
}


#endif // #ifndef _VOLUMERENDER_KERNEL_CU_
