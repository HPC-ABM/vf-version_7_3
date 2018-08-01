/* 
 * WHChemical.cpp
 * 
 * File Contents: Contains WHChemical class
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *
 * Created on Jun 19, 2013, 9:58 PM
 *
 *****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#include "WHChemical.h"
#include "../../World/Usr_World/woundHealingWorld.h"
#include <iostream>

using namespace std;

WHWorld* WHChemical::chemWorldPtr = NULL;

WHChemical::WHChemical() {

  this->np = 110*1390*1006; // Default world dimensions
  this->nx = 110;
  this->ny = 1390;
  this->nz = 1006;

  this->pbaseline = NULL;

#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
  this->pneugrad  = NULL;
  this->pmacgrad  = NULL;
  this->pfibgrad  = NULL;
  for (int igrad = 0; igrad < 3; igrad++)
      this->grads[igrad]     = NULL;
#else
  this->grads     = NULL;
#endif

#ifdef OPT_CHEM
  this->total = NULL;
#endif

  this->allocateWorldChem();

}

WHChemical::WHChemical(int nx, int ny, int nz) {

  this->np = nx*ny*nz;
  this->nx = nx;
  this->ny = ny;
  this->nz = nz;

  this->pbaseline = NULL;

#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
  this->pneugrad  = NULL;
  this->pmacgrad  = NULL;
  this->pfibgrad  = NULL;
  for (int igrad = 0; igrad < 3; igrad++)
      this->grads[igrad]     = NULL;
#else
  this->grads     = NULL;
#endif

#ifdef OPT_CHEM
  this->total = NULL;
#endif

  this->allocateWorldChem();
}

WHChemical::~WHChemical() {
printf("~WHChemical::WHChemical\n");
#ifdef OPT_CHEM
if (this->np <= 0) return;
#ifdef CHEM_PACK8
  if(this->pChem) delete(this->pChem);
  if(this->dChem) delete(this->dChem);
  if(this->tChem) delete(this->tChem);
#elif defined (CHEM_PACK4)
  if(this->pChem[0]) delete(this->pChem[0]);
  if(this->dChem[0]) delete(this->dChem[0]);
  if(this->tChem[0]) delete(this->tChem[0]);
  if(this->pChem[1]) delete(this->pChem[1]);
  if(this->dChem[1]) delete(this->dChem[1]);
  if(this->tChem[1]) delete(this->tChem[1]);
#else
  for (int ic = 0; ic < 8; ic++) {
    if(this->pChem[ic]) delete(this->pChem[ic]);
    if(this->dChem[ic]) delete(this->dChem[ic]);
    if(this->tChem[ic]) delete(this->tChem[ic]);
  }
#endif
  if (this->total) delete(this->total);

#else   // OPT_CHEM

#ifdef PAGELOCK_DIFFUSE
  printf("WHChemical: PAGELOCK DeAllocation\n");
  for (int ic = 0; ic < 8; ic++)
  {
      if (this->pChem[ic])
      {
        checkCudaErrors(cudaFreeHost(this->pChem[ic]));
        this->pChem[ic] = NULL;
      }
      if (this->dChem[ic]) delete (this->dChem[ic]);
      if (this->tChem[ic])
      {
        checkCudaErrors(cudaFreeHost(this->tChem[ic]));
        this->tChem[ic] = NULL;
      }
  }
  printf("WHChemical: Done deallocation\n");
#else	// PAGELOCK_DIFFUSE
  for (int ic = 0; ic < 8; ic++)
  {
      if (this->pChem[ic]) delete (this->pChem[ic]);
      if (this->dChem[ic]) delete (this->dChem[ic]);
#ifdef GPU_DIFFUSE
      if (this->tChem[ic]) delete (this->tChem[ic]);
#endif  // GPU_DIFFUSE
  }
#endif
#endif	// OPT_CHEM

#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
  if (this->pfibgrad) delete(this->pfibgrad);
  if (this->pmacgrad) delete(this->pmacgrad);
  if (this->pneugrad) delete(this->pneugrad);
#endif
}

void WHChemical::allocateWorldChem()
{
#ifdef OPT_CHEM

if (this->np <= 0)
{
  cerr << "ERROR: WHChemical not allocated" << endl;
  exit(-1);
}
#ifdef CHEM_PACK8
printf("CHEM_PACK8: requesting for %f (x3) GB\n", (sizeof(CS8)*np)/(1024.0*1024.0*1024.0));
  this->pChem = new CS8[np];
  this->dChem = new CS8[np];
  this->tChem = new CS8[np];

  this->grads = new CS3[np];
#elif defined (CHEM_PACK4)
printf("CHEM_PACK4: requesting for %f (x6) GB\n", (sizeof(CS4)*np)/(1024.0*1024.0*1024.0));
  this->pChem[0] = new CS4[np];
  this->dChem[0] = new CS4[np];
  this->tChem[0] = new CS4[np];
  this->pChem[1] = new CS4[np];
  this->dChem[1] = new CS4[np];
  this->tChem[1] = new CS4[np];

  this->grads = new CS3[np];
#else
  printf("CHEM_PACK1: requesting for %f (x3) GB\n", (sizeof(float)*np*8)/(1024.0*1024.0*1024.0));
  for (int ic = 0; ic < 8; ic++) {
    this->pChem[ic] = new float[np];
    this->dChem[ic] = new float[np];
    this->tChem[ic] = new float[np];
  }
  // Allocating gradients
  for (int igrad = 0; igrad < 3; igrad ++)
  {
      this->grads[igrad] = new float[np];
  }
#endif

  this->total = new float[8];

#else   // OPT_CHEM
#ifndef PAGELOCK_DIFFUSE
  for (int ic = 0; ic < 8; ic++)
  {
      this->pChem[ic] = new float[np];
      this->dChem[ic] = new float[np];
#ifdef GPU_DIFFUSE
      this->tChem[ic] = new float[np];
#endif  // GPU_DIFFUSE
  }

#else	// !PAGELOCK_DIFFUSE
  printf("WHChemical: PAGELOCK Allocation\n");
  for (int ic = 0; ic < 8; ic++)
  {
      this->dChem[ic] = new float[np];
      checkCudaErrors(cudaMallocHost((void**)&(this->pChem[ic]), np * sizeof(float)));
      checkCudaErrors(cudaMallocHost((void**)&(this->tChem[ic]), np * sizeof(float)));
  }
#endif	// PAGELOCK_DIFFUSE


  // Allocating gradients
  for (int igrad = 0; igrad < 3; igrad ++)
  {
      this->grads[igrad] = new float[np];
  }

  this->linkAttributes();
#endif	// OPT_CHEM

}

void WHChemical::resetTotals()
{
	if (total)
	{
		for (int ic = 0; ic < 8; ic++)
			this->total[ic] = 0.0f;
	} else {
		printf("WHChemical::resetTotals: this->total not yet allocated\n");
		exit(-1);
	}
}

void WHChemical::linkAttributes()
{
#ifdef OPT_CHEM

#else	// OPT_CHEM
	/* Link World attribute chemAllocation with WHWorldChem
   * (WHWorldChem->is linked to WHChemical) */
	this->pTNF      = this->pChem[TNF];
	this->pTGF      = this->pChem[TGF];
	this->pFGF      = this->pChem[FGF];
	this->pMMP8     = this->pChem[MMP8];
	this->pIL1beta  = this->pChem[IL1beta];
	this->pIL6 	    = this->pChem[IL6];
	this->pIL8 	    = this->pChem[IL8];
	this->pIL10     = this->pChem[IL10];
	this->dTNF 	    = this->dChem[TNF];
	this->dTGF 	    = this->dChem[TGF];
	this->dFGF 	    = this->dChem[FGF];
	this->dMMP8     = this->dChem[MMP8];
	this->dIL1beta  = this->dChem[IL1beta];
	this->dIL6 	    = this->dChem[IL6];
	this->dIL8 	    = this->dChem[IL8];
	this->dIL10     = this->dChem[IL10];
#ifdef GPU_DIFFUSE
    this->tTNF      = this->tChem[TNF];
    this->tTGF      = this->tChem[TGF];
    this->tFGF      = this->tChem[FGF];
    this->tMMP8     = this->tChem[MMP8];
    this->tIL1beta  = this->tChem[IL1beta];
    this->tIL6      = this->tChem[IL6];
    this->tIL8      = this->tChem[IL8];
    this->tIL10     = this->tChem[IL10];
#endif	// GPU_DIFFUSE
    this->pfibgrad  = this->grads[FIBgrad];
    this->pmacgrad  = this->grads[MACgrad];
    this->pneugrad  = this->grads[NEUgrad];
#endif	// OPT_CHEM
}

float WHChemical::getGrad(int ic, int in)
{
#ifdef OPT_CHEM
#ifdef CHEM_PACK8
    return this->grads[in].chem[ic];
#elif defined (CHEM_PACK4)
    return this->grads[in].chem[ic];
#else
    return this->grads[ic][in];
#endif
#else       // OPT_CHEM
    return this->grads[ic][in];
#endif      // OPT_CHEM
}

void  WHChemical::setGrad(int ic, int in, float level)
{
#ifdef OPT_CHEM
#ifdef CHEM_PACK8
    this->grads[in].chem[ic] = level;
#elif defined (CHEM_PACK4)
    this->grads[in].chem[ic] = level;
#else
    this->grads[ic][in] = level;
#endif
#else       // OPT_CHEM
    this->grads[ic][in] = level;
#endif      // OPT_CHEM
}

void  WHChemical::incGrad(int ic, int in, float inc)
{
#ifdef OPT_CHEM
#ifdef CHEM_PACK8
    this->grads[in].chem[ic] += inc;
#elif defined (CHEM_PACK4)
    this->grads[in].chem[ic] += inc;
#else
    this->grads[ic][in] += inc;
#endif
#else       // OPT_CHEM
    this->grads[ic][in] += inc;
#endif      // OPT_CHEM
}

float WHChemical::getLevel(int ic, int in, bool isGrad)
{
    if (isGrad) {
        return this->getGrad(ic, in);
    } else {
#ifdef OPT_CHEM
        return this->getPchem(ic, in);
#else
        return this->pChem[ic][in];
#endif
    }
}

void WHChemical::update(int totaldamage, float* cytokineDecay, chemupdate_t update_type)
{

#ifdef OPT_CHEM

	for (int ic = 0; ic < 8; ic++) {
		this->total[ic] = 0;
	}

#pragma omp parallel
	{
		float sum[8] = {0};
#pragma omp for
		for (int zi = 0; zi < nz; zi++) {
			for (int yi = 0; yi < ny; yi++) {
				for (int xi = 0; xi < nx; xi++) {
					int in = xi + yi*nx + zi*nx*ny;
					float dChem;
					float rChem;
					float plevel;
					// For each chem in this patch, update.
					for (int ic = 0; ic < 8; ic++) {
						dChem = this->getDchem(ic, in);
						this->setDchem(ic, in, 0);

						rChem = (update_type == gpuUpdate?
									this->getTchem(ic, in):
									this->getDchem(ic, in));

						plevel = dChem + rChem * cytokineDecay[ic];
						this->setPchem(ic, in, plevel);
						sum[ic] += plevel;
					}
					// Update gradients
					float patchIL1beta = this->getPchem(IL1beta, in);
					float patchTNF     = this->getPchem(TNF, in);
					float patchTGF     = this->getPchem(TGF, in);
					float patchFGF     = this->getPchem(FGF, in);
					float patchIL6     = this->getPchem(IL6, in);
					float patchIL8     = this->getPchem(IL8, in);

					float grad = patchIL1beta + patchTNF + patchTGF + patchFGF +
					        this->chemWorldPtr->worldECM[in].fragmnt_col[read_t];
                    this->setGrad(FIBgrad, in, patchTGF);
                    this->setGrad(NEUgrad, in, grad + patchIL6 + patchIL8);
                    this->setGrad(MACgrad, in, grad + this->chemWorldPtr->worldECM[in].fragmnt_ela[read_t]);
				}
			}
		}

#pragma omp critical
		{
			//Initialize total chemical concentration
			for (int ic = 0; ic < 8; ic++)
				this->total[ic] += sum[ic];
		}
	}
#else		// OPT_CHEM

    float *rTNF,     *rTGF, *rFGF, *rMMP8;
    float *rIL1beta, *rIL6, *rIL8, *rIL10;

    switch (update_type)
    {
	case cpuUpdate:
		rTNF		= dTNF;
		rTGF		= dTGF;
		rFGF		= dFGF;
		rMMP8		= dMMP8;
		rIL1beta	= dIL1beta;
		rIL6		= dIL6;
		rIL8		= dIL8;
		rIL10		= dIL10;
		break;
	case gpuUpdate:
		rTNF		= tTNF;
		rTGF		= tTGF;
		rFGF		= tFGF;
		rMMP8		= tMMP8;
		rIL1beta	= tIL1beta;
		rIL6		= tIL6;
		rIL8		= tIL8;
		rIL10		= tIL10;
		break;
    }

    for (int ic = 0; ic < 8; ic++) {
        this->total[ic] = 0;

    }

#pragma omp parallel
{
    double sum[8] = {0};
#pragma omp for
    for (int in = 0; in < np; in++)
    {
        for (int ic = 0; ic < 8; ic ++)
        {
            float level = getDchem(ic, in) + getTchem(ic, in) * cytokineDecay[ic]; 
            setDchem(ic, in, 0);
            setPchem(ic, in, level);
            sum[ic] += level;
        }
	// Update gradient
	float patchIL1beta = this->pIL1beta[in];
	float patchTNF = this->pTNF[in];
	float patchTGF = this->pTGF[in];
	float patchFGF = this->pFGF[in];
	float grad = patchIL1beta + patchTNF + patchTGF + patchFGF + this->chemWorldPtr->worldECM[in].fragmnt_col[read_t];
	this->pfibgrad[in] = patchTGF;
	this->pneugrad[in] = grad + this->pIL6[in] + this->pIL8[in];
	this->pmacgrad[in] = grad + this->chemWorldPtr->worldECM[in].fragmnt_ela[read_t];
    }
#pragma omp critical
    {
        for (int ic = 0; ic < 8; ic++)
        {
            this->total[ic] += sum[ic];
        }
    }
}   
#endif		// OPT_CHEM

}

void WHChemical::updateChemGPU(int totaldamage, float* cytokineDecay) {
	this->update(totaldamage, cytokineDecay, gpuUpdate);
}

// Always called in non-GPU_DIFFUSE
// Called only at beginning otherwise
void WHChemical::updateChemCPU(int totaldamage, float* cytokineDecay) {
	this->update(totaldamage, cytokineDecay, cpuUpdate);
}


/****************************************************************
 * GETTERS                                                      *
 ****************************************************************/

float* WHChemical::getPchemPtr(int ic)
{
#ifdef OPT_CHEM

// return selected column
#ifdef CHEM_PACK8
  return &(this->pChem[0].chem[ic]);
#elif defined (CHEM_PACK4)
  return &(this->pChem[ic/4][0].chem[ic%4]);
#else
  return this->pChem[ic];
#endif

#else   // OPT_CHEM

  return this->pChem[ic];

#endif  // OPT_CHEM
}

float* WHChemical::getTchemPtr(int ic)
{
#ifdef OPT_CHEM

// return selected column
#ifdef CHEM_PACK8
  return &(this->tChem[0].chem[ic]);
#elif defined (CHEM_PACK4)
  return &(this->tChem[ic/4][0].chem[ic%4]);
#else
  return this->tChem[ic];
#endif

#else   // OPT_CHEM

  return this->tChem[ic];

#endif  // OPT_CHEM
}

// TODO: make inline
float WHChemical::getPchem(int ic, int in)
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
  return this->pChem[ic][in];
#else

#ifdef CHEM_PACK8
  return this->pChem[in].chem[ic];
#elif defined (CHEM_PACK4)
  return this->pChem[ic/4][in].chem[ic%4];
#endif

#endif
}

// TODO: make inline
float WHChemical::getDchem(int ic, int in)
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
  return this->dChem[ic][in];
#else

#ifdef CHEM_PACK8
  return this->dChem[in].chem[ic];
#elif defined (CHEM_PACK4)
  return this->dChem[ic/4][in].chem[ic%4];
#endif

#endif
}

// TODO: make inline
float WHChemical::getTchem(int ic, int in)
{
#ifndef GPU_DIFFUSE
  return this->pChem[ic][in];
#else	// GPU_DIFFUSE

#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
  return this->tChem[ic][in];
#else

#ifdef CHEM_PACK8
  return this->tChem[in].chem[ic];
#elif defined (CHEM_PACK4)
  return this->tChem[ic/4][in].chem[ic%4];
#endif

#endif

#endif	// GPU_DIFFUSE
}





/****************************************************************
 * SETTERS/MODIFIERS                                            *
 ****************************************************************/

/*inline*/ void WHChemical::setPchem(int ic, int in, float clevel) 
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
    this->pChem[ic][in] = clevel;
#else

#ifdef CHEM_PACK8
  this->pChem[in].chem[ic] = clevel;
#elif defined (CHEM_PACK4)
  this->pChem[ic/4][in].chem[ic%4] = clevel;
#endif

#endif
}

/*inline*/ void WHChemical::setDchem(int ic, int in, float clevel) 
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
    this->dChem[ic][in] = clevel;
#else

#ifdef CHEM_PACK8
  this->dChem[in].chem[ic] = clevel;
#elif defined (CHEM_PACK4)
  this->dChem[ic/4][in].chem[ic%4] = clevel;
#endif

#endif
}

/*inline*/ void WHChemical::setTchem(int ic, int in, float clevel) 
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
    this->tChem[ic][in] = clevel;
#else

#ifdef CHEM_PACK8
  this->tChem[in].chem[ic] = clevel;
#elif defined (CHEM_PACK4)
  this->tChem[ic/4][in].chem[ic%4] = clevel;
#endif

#endif
}


/*inline*/ void WHChemical::incPchem(int ic, int in, float clevel) 
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
    this->pChem[ic][in] += clevel;
#else

#ifdef CHEM_PACK8
  this->pChem[in].chem[ic] += clevel;
#elif defined (CHEM_PACK4)
  this->pChem[ic/4][in].chem[ic%4] += clevel;
#endif

#endif
}

/*inline*/ void WHChemical::incDchem(int ic, int in, float clevel) 
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
    this->dChem[ic][in] += clevel;
#else

#ifdef CHEM_PACK8
  this->dChem[in].chem[ic] += clevel;
#elif defined (CHEM_PACK4)
  this->dChem[ic/4][in].chem[ic%4] += clevel;
#endif

#endif
}

/*inline*/ void WHChemical::incTchem(int ic, int in, float clevel) 
{
#if !defined (CHEM_PACK4) && !defined (CHEM_PACK8)
    this->tChem[ic][in] += clevel;
#else

#ifdef CHEM_PACK8
  this->tChem[in].chem[ic] += clevel;
#elif defined (CHEM_PACK4)
  this->tChem[ic/4][in].chem[ic%4] += clevel;
#endif

#endif
}


#ifdef OPT_CHEM

void WHChemical::unpackChem(c_ctx* chem_cctx, int ic)
{
	int np = this->np;
	float **ibuffs = chem_cctx->h_ibuffs;
#pragma omp parallel for num_threads(PACKTH)
	for (int in = 0; in < np; in++) {
		float temp = this->getPchem(ic, in);
		ibuffs[ic][in] = temp;
	}
}

void WHChemical::packChem(c_ctx* chem_cctx, int ic)
{
	int np = this->np;
	float **obuffs = chem_cctx->h_obuffs;
#pragma omp parallel for num_threads(PACKTH)
	for (int in = 0; in < np; in++) {
		float level = obuffs[ic][in];
		this->setTchem(ic, in , level);
	}
}

void WHChemical::unpackAllChem(c_ctx* chem_cctx)
{
	int np = this->np;
	float **ibuffs = chem_cctx->h_ibuffs;
#pragma omp parallel for num_threads(PACKTH)
//    for (int ic = 0; ic < 8; ic++) {
	for (int in = 0; in < np; in++) {
		for (int ic = 0; ic < 8; ic++) {
			float temp = this->getPchem(ic, in);
			ibuffs[ic][in] = temp;
		}
	}
}

void WHChemical::packAllChem(c_ctx* chem_cctx)
{
	int np = this->np;
	float **obuffs = chem_cctx->h_obuffs;
#pragma omp parallel for num_threads(PACKTH)
	for (int in = 0; in < np; in++) {
		for (int ic = 0; ic < 8; ic++) {
			float level = obuffs[ic][in];
			this->setTchem(ic, in , level);
		}
	}
}

#endif		// OPT_CHEM



