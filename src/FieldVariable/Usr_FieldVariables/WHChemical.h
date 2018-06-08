/* 
 * WHChemical.h
 * 
 * File Contents: Contains declarations for WHChemical class
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

#ifndef WHChemical_H
#define	WHChemical_H

#include "../FieldVariable.h"
#include "../../Diffusion/convolutionFFT_common.h"

#include "../../common.h"
#include <omp.h>
#include <cuda_runtime.h>

#include <vector>
using namespace std;
class World;
class WHWorld;

/* 
 * WHCHEMICAL CLASS DESCRIPTION:
 * WHChemical is a derived class of the parent class FieldVariable. It contains
 * all data members related to chemical concentration and gradients. Can be
 * used to access chemical concentrations by name.
 */
class WHChemical: public FieldVariable {
 public:

    /*
     * Description:	Default WHChemical constructor
     *
     * Return: void
     *
     * Parameters: void
     */
    WHChemical();

    /*
     * Description:	WHChemical constructor
     *
     * Return: void
     *
     * Parameters: nx  -- x-coordinate of the patch the WHChemical is on
     *             ny  -- y-coordinate of the patch the WHChemical is on
     *             nz  -- z-coordinate of the patch the WHChemical is on
     */
    WHChemical(int nx, int ny, int nz);

    /*
     * Description: WHChemical destructor
     *
     * Return: void
     *
     * Parameters: void
     */
    ~WHChemical();

  void allocateWorldChem();
  void resetTotals();
  void linkAttributes();
  void update(int totaldamage, float* cytokineDecay, chemupdate_t update_type);
  void updateChemGPU(int totaldamage, float* cytokineDecay);

  float* getPchemPtr(int ic);
  float* getTchemPtr(int ic);

  float getPchem(int ic, int in);
  float getDchem(int ic, int in);
  float getTchem(int ic, int in);
  void setPchem(int ic, int in, float clevel);
  void setDchem(int ic, int in, float clevel);
  void setTchem(int ic, int in, float clevel);
  void incPchem(int ic, int in, float clevel);
  void incDchem(int ic, int in, float clevel);
  void incTchem(int ic, int in, float clevel);

  float getGrad(int ic, int in);
  void  setGrad(int ic, int in, float level);
  void  incGrad(int ic, int in, float inc);
  float getLevel(int ic, int in, bool isGrad);

#ifdef OPT_CHEM
  void unpackChem(c_ctx* chem_cctx, int ic);
  void packChem(c_ctx* chem_cctx, int ic);
  void unpackAllChem(c_ctx* chem_cctx);
  void packAllChem(c_ctx* chem_cctx);
#endif

  /*
   * Description:	Update chemicals to reflect next tick's states.
   * 			This is called once at the initialization state for GPU_DIFFUSE.
     * 			Update in the following manner:
     * 				p<chem> = d<chem> + t<chem>*(1-gamma)
     * 			where gamma is a cytokine specific constant derived from the cytokine's
     * 			halflife.
     * 				gamma = 1 - 2^(-1/halflife)
   *
   * Return: void
   *
   * Parameters: void
   */
  void updateChemCPU(int totaldamage, float* cytokineDecay);

  typedef struct ChemStruct3 {      // for gradient packing
   float chem[3];
  } CS3;

  typedef struct ChemStruct4 {
   float chem[4];
  } CS4;

  typedef struct ChemStruct8 {
   float chem[8];
  } CS8;


  // TODO: Make private

#ifdef OPT_CHEM

#ifdef CHEM_PACK8
  CS8 *pChem;
  CS8 *dChem;
  CS8 *tChem;

  CS3 *grads;       // gradients
#elif defined (CHEM_PACK4)
  CS4 *pChem[2];
  CS4 *dChem[2];
  CS4 *tChem[2];

  CS3 *grads;       // gradients
#else
  float *pChem[8];
  float *dChem[8];
  float *tChem[8];

  float *grads[3];
  float *pfibgrad, *pmacgrad, *pneugrad;       // gradients
#endif

  float *total;



#else		// OPT_CHEM
#ifdef GPU_DIFFUSE

    /* The following variables keep track of a patch's cytokine level immediately
     * after the GPU diffusion is performed */
    float *tTNF,     *tTGF, *tFGF, *tMMP8;
    float *tIL1beta, *tIL6, *tIL8, *tIL10;

    float *tChem[8];

#endif

    float *pChem[8];
    float *dChem[8];
    float *grads[3];

    /* The following variables keep track of a patch's cytokine levels (pTNF,
     * pTGF, pFGF,...) and the change in a patch's cytokine levels throughout
     * the current tick (dTNF, dTGF, dFGF,...) */
    float* pTNF, *dTNF;
    float* pTGF, *dTGF;
    float* pFGF, *dFGF;
    float* pMMP8, *dMMP8;
    float* pIL1beta, *dIL1beta;
    float* pIL6, *dIL6;
    float* pIL8, *dIL8;
    float* pIL10, *dIL10;

    /* Keep track of the strength of the gradients that attract fibroblasts,
     * macrophages, and neutrophils respectively. */
    float *pfibgrad, *pmacgrad, *pneugrad;

    // Keep track of the total cytokine levels in the world
    float total[8];
//    float totalTNF, totalTGF, totalFGF, totalMMP8, totalIL1beta, totalIL6, totalIL8,totalIL10; 
#endif		// OPT_CHEM

    float *pbaseline;

    int np;     // number of patches
    int nx;
    int ny;
    int nz;

    // Pointer from an chemical data object to a WHWorld
    static WHWorld* chemWorldPtr;

};

#endif	/* WHChemical_H */

