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

#include <helper_cuda.h>
#include <helper_math.h>

#include "../../enums.h"


typedef unsigned int  uint;
typedef unsigned char uchar;

cudaArray *d_volumeArray[m_ecmtotal] = {0};
cudaArray *d_transferFuncArrayCol = {0};
cudaArray *d_transferFuncArrayEla = {0};
cudaArray *d_transferFuncArrayHya = {0};

//typedef unsigned char VolumeType;
typedef float VolumeType;

texture<VolumeType, 3, cudaReadModeElementType> texCol;
texture<VolumeType, 3, cudaReadModeElementType> texEla;
texture<VolumeType, 3, cudaReadModeElementType> texHya;
//texture<VolumeType, 3, cudaReadModeNormalizedFloat> tex;         // 3D texture
texture<float4, 1, cudaReadModeElementType>     transferTexCol; // 1D transfer function texture
texture<float4, 1, cudaReadModeElementType>     transferTexEla;
texture<float4, 1, cudaReadModeElementType>     transferTexHya;

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

__global__ void
d_render_dim(uint *d_output, uint nx, uint ny, uint nz, uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale, ecm_i ecmType)
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

    	  float posx = (pos.x + x_halfwidth)/(2.0f*x_halfwidth);
    	  float posy = (pos.y + y_halfwidth)/(2.0f*y_halfwidth);
    	  float posz = (pos.z + z_halfwidth)/(2.0f*z_halfwidth);

    	  float sample;
    	  float4 col;

    	  switch (ecmType)
    	  {
    	  case m_col:
    	  	sample = tex3D(texCol, posx, posy, posz);

    	  	//sample *= 64.0f;    // scale for 10-bit data

    	  	// lookup in transfer function texture
    	  	col = tex1D(transferTexCol, (sample-transferOffset)*transferScale);
    	  	break;
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
void bufferECMmap(cudaMemcpy3DParms copyParams)
{
	checkCudaErrors(cudaMemcpy3D(&copyParams));
}

extern "C"
void initCuda(void *h_volume, cudaExtent volumeSize, cudaMemcpy3DParms &copyParams, ecm_i ecmType)
{
    // create 3D array
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<VolumeType>();
    checkCudaErrors(cudaMalloc3DArray(&(d_volumeArray[ecmType]), &channelDesc, volumeSize));

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
//    	float4 transferFunc[] =
//    	{
//    			{  0.0,  0.0,  0.0, 0.0, },	// 0.00
//    			{  1.0,  0.0,  0.0, 0.5, },	// 0.05 - SLP ILP
//    			{  0.0,  0.0,  0.0, 0.0, },	// 0.10 - SLP ILP
//    			{  0.0,  0.0,  0.0, 0.0, },	// 0.15 - SLP ILP
//    			{  0.0,  0.0,  0.0, 0.0, }, // 0.20 - DLP ILP SLP
//    			{  0.0,  0.0,  0.0, 0.0, }, // 0.25 - DLP
//    			{  1.0,  0.0,  0.0, 0.2, }, // 0.30 - DLP
//    			{  1.0,  0.0,  0.0, 0.2, }, // 0.35 - DLP
//    			{  1.0,  0.0,  0.0, 0.2, }, // 0.40 - DLP
//    			{  1.0,  0.0,  0.0, 0.2, }, // 0.45
//    			{  1.0,  0.0,  0.0, 0.2, }, // 0.50
//    			{  1.0,  0.0,  0.0, 0.3, }, // 0.55
//    			{  1.0,  0.0,  0.0, 0.4, }, // 0.60
//    			{  1.0,  0.0,  0.0, 0.5, }, // 0.65
//    			{  1.0,  0.0,  0.0, 0.6, }, // 0.70
//    			{  1.0,  0.0,  0.0, 0.7, }, // 0.75
//    			{  1.0,  0.0,  0.0, 0.8, }, // 0.80
//    			{  1.0,  0.1,  0.1, 1.0, }, // 0.85
//    			{  1.0,  0.2,  0.2, 0.5, }, // 0.90
//    			{  1.0,  0.3,  0.3, 0.7, }, // 0.95
//    			{  1.0,  0.6,  0.0, 0.8, }, // 1.00
//    	};
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
      // Elastin
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
      // set texture parameters
      texHya.normalized = true;                      // access with normalized texture coordinates
      texHya.filterMode = cudaFilterModeLinear;      // linear interpolation
      texHya.addressMode[0] = cudaAddressModeClamp;  // clamp texture coordinates
      texHya.addressMode[1] = cudaAddressModeClamp;

      // bind array to 3D texture
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


//        float4 transferFunc[] =
//        {
//        		{  0.0, 0.0, 0.0, 0.0, },
//        		{  1.0, 0.0, 0.0, 0.1, },
//        		{  0.0, 0.0, 0.0, 0.0, },
//        		{  0.0, 0.0, 0.0, 0.0, },
//        		{  1.0, 0.0, 0.0, 0.00001, },//0.1
//        		{  0.97, 0.8, 0.72, 0.01, },
//        		{  0.0, 0.0, 0.0, 0.001, },
//        		{  0.0, 0.0, 0.0, 0.001, },
//        		{  0.0, 0.0, 0.0, 0.001, },		// red
//        		{  1.0, 0.0, 0.0, 0.001, },		// red
//        		{  0.97, 0.8, 0.72, 0.01, },		// red
//        		{  0.97, 0.8, 0.72, 0.01, },		// red
//        		//    		{  1.0, 0.0, 0.0, 0.1, },		// red
//        		{  0.97, 0.8, 0.72, 0.12, },		// flesh
//        };


//    float4 transferFunc[] =
//    {
//    		{  0.0, 0.0, 0.0, 0.0, },
//    		{  1.0, 0.0, 0.0, 0.1, },
//    		{  0.0, 0.0, 0.0, 0.0, },
//    		{  0.0, 0.0, 0.0, 0.0, },
//    		{  1.0, 0.0, 0.0, 0.00001, },//0.1
//    		{  0.97, 0.8, 0.72, 0.01, },
//    		{  0.0, 0.0, 0.0, 0.001, },
//    		{  0.0, 0.0, 0.0, 0.001, },
//    		{  0.0, 0.0, 0.0, 0.001, },		// red
//    		{  1.0, 0.0, 0.0, 0.001, },		// red
//    		{  0.97, 0.8, 0.72, 0.01, },		// red
//    		{  0.97, 0.8, 0.72, 0.01, },		// red
//    		//    		{  1.0, 0.0, 0.0, 0.1, },		// red
//    		{  0.97, 0.8, 0.72, 0.12, },		// flesh
//    };
//    float4 transferFunc[] =
//    {
//        {  0.5, 0.0, 0.2, 0.0 },      // 0.0
//        {  0.5, 0.0, 0.2, 0.100 },
//        {0.576, 1.439, 1.000, 0.1000 },  // 0.1
//        {1.000, 0.753, 0.796, 0.300 }, // 0.2
//        {0.000, 0.749, 1.000, 0.400 },  // 0.3
//        {0.000, 0.749, 1.000, 0.500 },  // 0.4
//        {0.498, 1.000, 0.831, 0.600 },  // 0.5
//        {1.000, 0.489, 0.314, 0.700 }, // 0.6
//        {1.000, 0.750, 0.700, 0.750 }, // 0.7
//        {1.000, 1.000, 0.000, 0.800 }, // 0.8
//        {0.500, 1.000, 0.200, 0.900 }, // 0.9
//        {0.000, 1.000, 0.000, 1.000 }, // 1.0
//        {0.999, 1.000, 0.999, 1.000 }, // 1.0
//    };



}

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
                   float density, float brightness, float transferOffset, float transferScale, ecm_i ecmType)
{
    d_render_dim<<<gridSize, blockSize>>>(d_output, nx, ny, nz, imageW, imageH, density,
                                      brightness, transferOffset, transferScale, ecmType);
}

extern "C"
void copyInvViewMatrix(float *invViewMatrix, size_t sizeofMatrix)
{
    checkCudaErrors(cudaMemcpyToSymbol(c_invViewMatrix, invViewMatrix, sizeofMatrix));
}


#endif // #ifndef _VOLUMERENDER_KERNEL_CU_
