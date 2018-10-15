/* 
 * File: Fibroblast.cpp
 *
 * File Contents: Contains the Fibroblast class.
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

#include "Fibroblast.h"
#include "../../enums.h"
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;
int Fibroblast::numOfFibroblasts = 0;


float Fibroblast::cytokineSynthesis[32] = {0.0f};
float Fibroblast::activation[5] = {10.0, 50.0, 0, 25.0, 2.5};
float Fibroblast::ECMsynthesis[19] = {1, 0, 0.01, 50, 25, 2, 10, 5, 1, 0, 25, 2, 1, 0, 50, 5, 10, 12, 1};
float Fibroblast::proliferation[6] = {24, 10, 1, 0, 25, 3};



Fibroblast::Fibroblast() {
}

Fibroblast::Fibroblast(Patch* patchPtr) {
  int tid = 0;
#ifdef _OMP
  // Get thread id in order to access the seed that belongs to this thread
  tid = omp_get_thread_num();
#endif
  this->ix[write_t] = patchPtr->indice[0];
  this->iy[write_t] = patchPtr->indice[1];
  this->iz[write_t] = patchPtr->indice[2];
  this->index[write_t] = patchPtr->index;
  this->alive[write_t] = true;
#ifdef DISABLE_RAND
  this->life[write_t] = WHWorld::reportTick(0, 103);  // 0 corresponds to hours. TODO(Kim:) INSERT REF
#else
  this->life[write_t] = WHWorld::reportTick(0, 100 + rand_r(&(agentWorldPtr->seeds[tid]))%7);  // Unactivated fibroblasts live for 100 to 107 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
#endif
  this->activate[write_t] = false;
  this->color[write_t] = cfibroblast;
  this->size[write_t] = 2;
  this->type[write_t] = fibroblast;
  this->ix[read_t] = patchPtr->indice[0];
  this->iy[read_t] = patchPtr->indice[1];
  this->iz[read_t] = patchPtr->indice[2];
  this->index[read_t] = patchPtr->index;

  /* In OMP version, we wait to add cells at the end of the tick,
   * whereas in serial version, cells are always added right away.
   * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
  this->alive[read_t] = true;
  this->life[read_t] = this->life[write_t];
#else
  this->alive[read_t] = false;
  this->life[read_t] = 0;
#endif
  this->activate[read_t] = false;
  this->color[read_t] = cfibroblast;
  this->size[read_t] = 2;
  this->type[read_t] = fibroblast;

  // TODO(NS): INSERT REF
  this->maxProlif[0] = 1;
  this->maxProlif[1] = 2;
  this->maxProlif[2] = 1;

  Fibroblast::numOfFibroblasts++;
}

Fibroblast::Fibroblast(int x, int y, int z) {
  // DEBUG
  /*
#pragma omp critical
{
   (Agent::agentWorldPtr->newfibs)++;
}
   */
  int tid = 0;

#ifdef _OMP
  // Get thread id in order to access the seed that belongs to this thread
  tid = omp_get_thread_num();
#endif
  this->ix[write_t] = x;
  this->iy[write_t] = y;
  this->iz[write_t] = z;
  this->index[write_t] = x + y*nx + z*nx*ny;
  this->alive[write_t] = true;
#ifdef DISABLE_RAND
  this->life[write_t] = WHWorld::reportTick(0, 103);  // 0 corresponds to hours. TODO(Kim:) INSERT REF
#else
  this->life[write_t] = WHWorld::reportTick(0, 100 + rand_r(&(agentWorldPtr->seeds[tid]))%7);  // Unactivated fibroblasts live for 100 to 107 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
#endif
  this->activate[write_t] = false;
  this->color[write_t] = cfibroblast;
  this->size[write_t] = 2;
  this->type[write_t] = fibroblast;

  /* In OMP version, we wait to add cells at the end of the tick,
   * whereas in serial version, cells are always added right away.
   * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
  this->ix[read_t] = x;
  this->iy[read_t] = y;
  this->iz[read_t] = z;
  this->index[read_t] = x + y*nx + z*nx*ny;
  this->alive[read_t] = true;
  this->life[read_t] = this->life[write_t];
  this->activate[read_t] = false;
  this->color[read_t] = cfibroblast;
  this->size[read_t] = 2;
  this->type[read_t] = fibroblast;
#else
  this->ix[read_t] = x;
  this->iy[read_t] = y;
  this->iz[read_t] = z;
  this->index[read_t] = x + y*nx + z*nx*ny;
  this->alive[read_t] = false;
  this->life[read_t] = 0;
  this->activate[read_t] = false;
  this->color[read_t] = cfibroblast;
  this->size[read_t] = 2;
  this->type[read_t] = fibroblast;
#endif
  // TODO(NS): INSERT REF
  this->maxProlif[0] = 1;
  this->maxProlif[1] = 2;
  this->maxProlif[2] = 1;

  Fibroblast::numOfFibroblasts++;
}

Fibroblast::~Fibroblast() {
}

void Fibroblast::cellFunction() {
  if (this->alive[read_t] == false) return;

  if (this->activate[read_t] == false) {
      this->fib_cellFunction();
  } else {
      this->afib_cellFunction();
  }
}

void Fibroblast::fib_cellFunction() {

  /*************************************************************************
   * MOVEMENT                                                              *
   *************************************************************************/
  this->fib_migrate();

  /*************************************************************************
   * ACTIVATION                                                            *
   *************************************************************************/
  this->fib_activate();

  /*************************************************************************
   * DEGRADE                                                               *
   *************************************************************************/
  this->fib_degrade();

}

void Fibroblast::afib_cellFunction() {

  /*************************************************************************
   * MOVEMENT                                                              *
   *************************************************************************/
  bool isInDam = this->afib_migrate();

  /*************************************************************************
   * PROLIFERATION                                                         *
   *************************************************************************/
  this->afib_proliferate();

  /*************************************************************************
   * CYTOKINE SYNTHESIS                                                    *
   *************************************************************************/
// DEBUG
  float countnHA = this->countNeighborECM(nha);
  this->afib_produce_cytokines(countnHA);

  /*************************************************************************
   * ECM PROTEIN SYNTHESIS                                                 *
   *************************************************************************/
  if (isInDam) this->afib_produce_ecms();


  /*************************************************************************
   * DEACTIVATION                                                          *
   *************************************************************************/
  this->afib_deactivate();


  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  this->afib_degrade();

}


void Fibroblast::fib_migrate(){
  this->wiggle();
}

void Fibroblast::fib_activate(){
  /* An unactivated fibroblast can be activated anywhere in the world.
   * TODO(NS): INSERT REF? */

#ifndef CALIBRATION
  int chance_low =   1;//25;
  int chance_med =  50;
  int chance_hig = 100;
#else   // CALIBRATION
  int chance_low = Fibroblast::activation[0];
  int chance_med = Fibroblast::activation[1];
  int chance_hig = Fibroblast::activation[2];

  assert(chance_low < chance_med);
  assert(chance_med < chance_hig);
#endif

  int in = this->index[read_t];
  int patchTGF = agentWorldPtr->WHWorldChem->pTGF[in];

  bool isHig = (patchTGF > FIB_ACT_TGF_UB);
  bool isMed = (FIB_ACT_TGF_LB < patchTGF && patchTGF <= FIB_ACT_TGF_UB);
  bool isLow = (patchTGF <= FIB_ACT_TGF_LB);

  // Larger than upper_bound       -> High activation chance
  // Between lower_ and upper_bound -> Medium activation chance
  // Lower than lower_bound        -> Low activation chance
  bool isAct = (isHig && rollDice(chance_hig)) ||
      (isMed && rollDice(chance_med)) ||
      (isLow && rollDice(chance_low)) ;

  if (isAct) this->fibActivation();

}

void Fibroblast::fib_degrade(){
  // Decrement life by 1 tick
  // Unactivated fibroblasts can die naturally
  this->life[write_t] = this->life[read_t] - 1;
  if (this->life[read_t] <= 0) {
      this->die();
  }
}

void Fibroblast::afib_degrade(){
  // Number of proliferation left in the last interval
  int prolifLeft = this->maxProlif[2];

  // Dead if maximum number of proliferation reached
  if (prolifLeft <= 0) this->die();
}

void Fibroblast::afib_proliferate(){

  double hr = Agent::agentWorldPtr->reportHour();
  int TGFrelated = 0;

  // Update quota 2 to include 1
  if (hr == FIB_PROLIF_I2) {
      this->maxProlif[2] += this->maxProlif[1];
      this->maxProlif[1] = 0;
  }

  bool isInI0 = (FIB_PROLIF_I0 < hr) && (hr < FIB_PROLIF_I1);     // day 4-7
  bool isInI1 = (FIB_PROLIF_I1 < hr) && (hr < FIB_PROLIF_I2);     // day 8-14
  bool isInI2 = (FIB_PROLIF_I2 < hr);                             // day 15+

  bool isInQ0 = isInI0 && (this->maxProlif[0] > 0);
  bool isInQ1 = isInI1 && (this->maxProlif[1] > 0);
  bool isInQ2 = isInI2 && (this->maxProlif[2] > 0);

  if (isInQ0 || isInQ1 || isInQ2) {
#ifndef CALIBRATION
      float a  = 1.0/5.0;//1.0;
      float b  = 0;//25;
#else
      float a  = 1.0/Fibroblast::proliferation[5];
      float b  = Fibroblast::proliferation[3];
#endif

      float meanTGF = this->meanNeighborChem(TGF);
      if      (meanTGF < FIB_PRO_TGF_MX) TGFrelated = +1;
      else if (meanTGF > FIB_PRO_TGF_MN) TGFrelated = -1;

      float meanTNF = this->meanNeighborChem(TNF);
      float meanFGF = this->meanNeighborChem(FGF);
      float meanIL1 = this->meanNeighborChem(IL1beta);
      int  countfHA = this->countNeighborECM(fha);

      float fibProlif = log10(1 + meanTNF + meanFGF + meanIL1 +
          TGFrelated*meanTGF + countfHA);
      // TODO(NS):INSERT REF? and Update Equation
      if (rollDice(a*fibProlif + b)) {
          this->hatchnewfibroblast(2);
          this->die();
          return;
      }
  }
}


float Fibroblast::calc_stim(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset,
    float max_stim)
{
  float stim = this->calc_chem(uregs, dregs, coefs, offset);

  if (stim <= 1.0f) return 0.0f;
  if (stim > max_stim) return 1.0f;	// clamp stim

  float min_stim = 1.0;

  return (stim - min_stim)/(max_stim - min_stim);
}

float Fibroblast::produce_tgf(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset){

  int in = this->getIndex();
  float TGFinc;
  if (Agent::agentPatchPtr[in].isInDamZone()) {
    TGFinc = this->calc_chem(uregs, dregs, coefs, offset);
  } else {
    // TODO(NS): TGF secretion rate outside of wound (initial a-fibs)
    TGFinc = 0.0f;
  }

  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dTGF[in]) += TGFinc;
  return TGFinc;
}

void Fibroblast::afib_produce_cytokines(float neighborHAs) {
  // Calculates chemical gradients and patch chemical concentrations
  int in = this->index[read_t];
  unsigned *seed_arr = (Agent::agentWorldPtr)->seeds;

  float patchTNF  = this->agentWorldPtr->WHWorldChem->pTNF[in];
  float patchTGF  = this->agentWorldPtr->WHWorldChem->pTGF[in];
  float patchIL1  = this->agentWorldPtr->WHWorldChem->pIL1beta[in];
  float patchIL6  = this->agentWorldPtr->WHWorldChem->pIL6[in];
  float patchIL10 = this->agentWorldPtr->WHWorldChem->pIL10[in];
  float patchfHA  = Agent::agentECMPtr[in].getnfHA();

  // Regulator coefficients and offsets
#ifndef CALIBRATION
  float tgf_cTNF  = 0.0000001f;//1.0f;
  float tgf_cIL10 = 100.0f;
  float tgf_ofst  = Agent::baseline[TGF];//0.0f;

  float fgf_cTGF = 5.2e-3;
  float fgf_ofst = util::randFloatRange(1.0f, 1.6f, seed_arr)*(1e-7);

  float tnf_cTGF  = 0.00000005f;//0.00000002f;//0.00001f;
  float tnf_cIL10 = 10000000.0f;
  float tnf_cnHA  = 1000.0f;
  float tnf_ofst  = Agent::baseline[TNF];//0.0f;

  float il6_cIL6  = 0.0000001f;//1.0f;
  float il6_cTGF  = 0.00000001f;//1.0f;
  float il6_cTNF  = 0.000001f;//1.0f;
  float il6_cIL1  = 0.0000001f;//1.0f;
  float il6_cnHA  = 100000.0f;
  float il6_cIL10 = 100000.0f;
  float il6_ofst  = Agent::baseline[IL6];//0.0f;

  float il8_cTGF  = 0.00000001f;//1.0f;
  float il8_cnHA  = 10000.0f;
  float il8_cIL10 = 10000.0f;
  float il8_ofst  = Agent::baseline[IL8];//0.0f;

	float il1_cTGF  = 0.0000000f;	// up
	float il1_cTNF  = 0.00000f;	// up
	float il1_cfHA  = 0.0000000f; // up
	float il1_cIL10 = 0.0f;//1000000.0f;	// down
	float il1_ofst  = Agent::baseline[IL1beta];

#else
  // TODO: Reorder parameters
  float tgf_cTNF  = Fibroblast::cytokineSynthesis[0];
  float tgf_cIL10 = Fibroblast::cytokineSynthesis[1];
  float tgf_ofst  = Agent::baseline[TGF];

  float fgf_cTGF = 5.2e-3;
  float fgf_ofst = util::randFloatRange(1.0f, 1.6f, seed_arr)*(1e-7);

  float tnf_cTGF  = Fibroblast::cytokineSynthesis[2];
  float tnf_cIL10 = Fibroblast::cytokineSynthesis[3];
  float tnf_cnHA  = Fibroblast::cytokineSynthesis[4];
  float tnf_ofst  = Agent::baseline[TNF];

  float il6_cIL6  = Fibroblast::cytokineSynthesis[5];
  float il6_cTGF  = Fibroblast::cytokineSynthesis[6];
  float il6_cTNF  = Fibroblast::cytokineSynthesis[7];
  float il6_cIL1  = Fibroblast::cytokineSynthesis[8];
  float il6_cnHA  = Fibroblast::cytokineSynthesis[9];
  float il6_cIL10 = Fibroblast::cytokineSynthesis[10];
  float il6_ofst  = Agent::baseline[IL6];

  float il8_cTGF  = Fibroblast::cytokineSynthesis[11];
  float il8_cnHA  = Fibroblast::cytokineSynthesis[12];
  float il8_cIL10 = Fibroblast::cytokineSynthesis[13];
  float il8_ofst  = Agent::baseline[IL8];

	float il1_cTGF  = Fibroblast::cytokineSynthesis[14];	// up
	float il1_cTNF  = Fibroblast::cytokineSynthesis[15];	// up
	float il1_cfHA  = Fibroblast::cytokineSynthesis[16];  // up
	float il1_cIL10 = Fibroblast::cytokineSynthesis[17];	// down
	float il1_ofst  = Agent::baseline[IL1beta];
#endif

  REGULATORS_T tgf_coefs{tgf_cTNF, tgf_cIL10};
  REGULATORS_T fgf_coefs{fgf_cTGF};
  REGULATORS_T tnf_coefs{tnf_cTGF, tnf_cIL10, tnf_cnHA};
  REGULATORS_T il6_coefs{il6_cIL6, il6_cTGF,  il6_cTNF, il6_cIL1, il6_cnHA, il6_cIL10};
  REGULATORS_T il8_coefs{il8_cTGF, il8_cnHA,  il8_cIL10};

  REGULATORS_T il1_ucoefs{il1_cTGF, il1_cTNF, il1_cfHA};
  REGULATORS_T il1_dcoefs{il1_cIL10};

  // Up- and down-regulators for each cytokine produced by myofibroblasts
  REGULATORS_T tgf_uregs{patchTNF};
  REGULATORS_T tgf_dregs{patchIL10};

  REGULATORS_T fgf_uregs{patchTGF};
  REGULATORS_T fgf_dregs{};        // none

  REGULATORS_T tnf_uregs{patchTGF};
  REGULATORS_T tnf_dregs{patchIL10, neighborHAs};

  REGULATORS_T il6_uregs{patchIL6, patchTGF, patchTNF, patchIL1};
  REGULATORS_T il6_dregs{patchIL10, neighborHAs};

  REGULATORS_T il8_uregs{patchTGF};
  REGULATORS_T il8_dregs{patchIL10, neighborHAs};

  REGULATORS_T il1_uregs{patchTGF, patchTNF, patchfHA};
  REGULATORS_T il1_dregs{patchIL10};


  // Call cytokine production functions
  float TGFinc = produce_tgf(tgf_uregs, tgf_dregs,  tgf_coefs, tgf_ofst);
  float FGFinc = produce_fgf(fgf_uregs, fgf_dregs,  fgf_coefs, fgf_ofst);
  float TNFinc = produce_tnf(tnf_uregs, tnf_dregs,  tnf_coefs, tnf_ofst);
  float IL6inc = produce_il6(il6_uregs, il6_dregs,  il6_coefs, il6_ofst);
  float IL8inc = produce_il8(il8_uregs, il8_dregs,  il8_coefs, il8_ofst);
//	float IL1inc = produce_il1(il1_uregs, il1_dregs, il1_ucoefs, il1_dcoefs, il1_ofst);


#ifdef PRINT_SECRETION
  int x = this->ix[read_t];
  int y = this->iy[read_t];
  int z = this->iz[read_t];
  printCytRelease(3, TGF, x, y, z, TGFinc);
  printCytRelease(3, FGF, x, y, z, FGFinc);
  printCytRelease(3, TNF, x, y, z, TNFinc);
  printCytRelease(3, IL6, x, y, z, IL6inc);
  printCytRelease(3, IL8, x, y, z, IL8inc);
  printCytRelease(3, IL1beta, x, y, z, IL1inc);
#endif  // PRINT_SECRETION

}

bool Fibroblast::afib_migrate(){
  int in = this->index[read_t];
  bool inDamZone = (Agent::agentPatchPtr[in]).inDamzone;
  if (!inDamZone) {
      if (this->moveToHighestChem(FIBgrad) != true) this->wiggle();
  } else {
      //this->wiggle();
  	this->wiggle_wound();
  }

  this->isMigrated = true;

  return inDamZone;
}

void Fibroblast::afib_produce_ecms(){
  // This function should only be called if this fib has migrated
//  assert(this->isMigrated);
  this->isMigrated = false;     // reset migration flag

  /*************************************************************************
   * Make damaged neighbor list                                            *
   *************************************************************************/
  int read_index;
  // Check if the location has been modified in this tick
  if (isModified(this->index)) {
      // If it has, work off of the intermediate value
      read_index = write_t;
  } else {
      // If it has NOT, work off of the original value
      read_index = read_t;
  }

  // DEBUG
  read_index = read_t;

  int dx, dy, dz;
  // Location of fibroblast in x,y,z dimensions of world.
  int x = this->ix[read_index];
  int y = this->iy[read_index];
  int z = this->iz[read_index];
  // Number of patches in x,y,z dimensions of world
  int nx = Agent::nx;
  int ny = Agent::ny;
  int nz = Agent::nz;
  int randInt, target, in;
  vector <int> damagedneighbors;

  // Make a list of damaged neighboring patches
#ifndef MODEL_3D
  int i_bgn =  9;
  int i_end = 18;
#else
  int i_bgn =  0;
  int i_end = 27;
#endif

  for (int i = i_bgn; i < i_end; i++) {
      dx = Agent::dX[i];
      dy = Agent::dY[i];
      dz = Agent::dZ[i];
      in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
      // Try a new neighboring patch if this one is outside the world dimensions.
      if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
      // Add the valid damaged neighboring patch to the list.
      if (Agent::agentPatchPtr[in].damage[read_t] != 0) {
          damagedneighbors.push_back(i);
      }
  }


  // Target a random damaged neighboring patch, if there are any.
  if (damagedneighbors.size() > 0) {
      int dX_col, dY_col, dZ_col;
      int dX_ela, dY_ela, dZ_ela;
      int dX_hya, dY_hya, dZ_hya;
      // get index changes for damaged neighbor for ecms
      this->get_rn_from_list(damagedneighbors, dX_col, dY_col, dZ_col);
      this->get_rn_from_list(damagedneighbors, dX_ela, dY_ela, dZ_ela);
      this->get_rn_from_list(damagedneighbors, dX_hya, dY_hya, dZ_hya);


      int dI_targets[3][3] = {
          dX_col, dY_col, dZ_col,
          dX_ela, dY_ela, dZ_ela,
          dX_hya, dY_hya, dZ_hya};

      // make monomers
      this->make_monomers(dI_targets);
  }

}

void Fibroblast::afib_deactivate() {
  // TODO(NS): Change to local stimulus sensing
  /* Activated fibroblasts might be deactivated once the damage is cleared. TODO(Kim): INSERT REF? */
  int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[pdamage];
	int initialdam = (Agent::agentWorldPtr)->getInitialDam();
	int threshold_dam = (initialdam*5)/100;	// 5% if the initial damage
#ifndef CALIBRATION
  if (totaldamage <= threshold_dam && rollDice(50)) this->fibDeactivation();
#else  // CALIBRATION
  if (totaldamage <= threshold_dam && rollDice(Fibroblast::activation[4])) this->fibDeactivation();
#endif  // CALIBRATION
}

void Fibroblast::get_rn_from_list(
    vector <int> damagedneighbors,
    int &dX,
    int &dY,
    int &dZ)
{
  int tid = 0;
#ifdef _OMP
  // Get thread id in order to access the seed that belongs to this thread
  tid = omp_get_thread_num();
#endif

  int target, randInt;
//  int n = damagedneighbors.size();
//  int i8 = -2;	// index of target = 8
//  int i = 0;
//  for (i = 0; i < n; i++) {
//  	int target = damagedneighbors[i];
//  	if (target > 8) {
//  		i8 = i-1;
//  		break;
//  	}
//  }
//
//  if (i8 == -2) i8 = i-1;
//
//
//  if (i8 == -1) {
//  	randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
//  } else if (i8 == 0) {
//  	randInt = 0;
//  } else {
//  	randInt = rand_r(&(agentWorldPtr->seeds[tid])) % i8;
//  }


  randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();

  target = damagedneighbors[randInt];
	dX = Agent::dX[target];
	dY = Agent::dY[target];
	dZ = Agent::dZ[target];
}

void Fibroblast::make_monomers(int dI_targets[3][3]) {

  const int icol = 0, iela = 1, ihya = 2;
  const int   ix = 0,   iy = 1,   iz = 2;

  int x, y, z, read_index;

  // Regulator coefficients and offsets
#ifndef CALIBRATION
  float col_cTGF  = 1.0f;
  float col_cTNF  = 1.0f;
  float col_cIL1  = 1.0f;
  float col_cIL6  = 1.0f;
  float col_cfHA  = 1.0f;
  float col_cFGF  = 1.0f;
  float col_cIL8  = 1.0f;
  float col_cIL10 = 1.0f;
  float col_ofst  = 0.0f;

  float ela_cTGF = 10.0f;
  float ela_cFGF = 1.0f;
  float ela_cIL1 = 1.0f;
  float ela_cTNF = 1.0f;
  float ela_ofst = 0.0f;

#ifdef RAT_VF
  float hya_cTGF  = 1.0f;
#else
  float hya_cTGF  = 20.0f;//10.0f;
#endif
  float hya_cIL1  = 1.0f;
  float hya_cIL6  = 1.0f;
  float hya_cIL8  = 1.0f;
  float hya_cIL10 = 1.0f;
  float hya_cFGF  = 1.0f;
  float hya_ofst  = 0.0f;
#else
  // TODO: Reorder parameters
  float col_cTGF  = Fibroblast::ECMsynthesis[0];
  float col_cTNF  = Fibroblast::ECMsynthesis[1];
  float col_cIL1  = Fibroblast::ECMsynthesis[2];
  float col_cIL6  = Fibroblast::ECMsynthesis[3];
  float col_cfHA  = Fibroblast::ECMsynthesis[4];
  float col_cFGF  = Fibroblast::ECMsynthesis[5];
  float col_cIL8  = Fibroblast::ECMsynthesis[6];
  float col_cIL10 = Fibroblast::ECMsynthesis[7];
  float col_ofst  = Fibroblast::ECMsynthesis[8];

  float ela_cTGF = Fibroblast::ECMsynthesis[9];
  float ela_cFGF = Fibroblast::ECMsynthesis[10];
  float ela_cIL1 = Fibroblast::ECMsynthesis[11];
  float ela_cTNF = Fibroblast::ECMsynthesis[12];
  float ela_ofst = Fibroblast::ECMsynthesis[13];

  float hya_cTGF  = Fibroblast::ECMsynthesis[14];
  float hya_cIL1  = Fibroblast::ECMsynthesis[15];
  float hya_cIL6  = Fibroblast::ECMsynthesis[16];
  float hya_cIL8  = Fibroblast::ECMsynthesis[17];
  float hya_cIL10 = Fibroblast::ECMsynthesis[18];
  float hya_cFGF  = Fibroblast::ECMsynthesis[19];
  float hya_ofst  = Fibroblast::ECMsynthesis[20];
#endif

  float meanTNF  = Agent::agentWorldPtr->WHWorldChem->pChem[TNF][this->getIndex()];//this->meanNeighborChem(TNF);
  float meanTGF  = Agent::agentWorldPtr->WHWorldChem->pChem[TGF][this->getIndex()];//this->meanNeighborChem(TGF);
  float meanFGF  = Agent::agentWorldPtr->WHWorldChem->pChem[FGF][this->getIndex()];//this->meanNeighborChem(FGF);
  float meanIL1  = Agent::agentWorldPtr->WHWorldChem->pChem[IL1beta][this->getIndex()];//this->meanNeighborChem(IL1beta);
  float meanIL6  = Agent::agentWorldPtr->WHWorldChem->pChem[IL6][this->getIndex()];//this->meanNeighborChem(IL6);
  float meanIL8  = Agent::agentWorldPtr->WHWorldChem->pChem[IL8][this->getIndex()];//this->meanNeighborChem(IL8);
  float meanIL10 = Agent::agentWorldPtr->WHWorldChem->pChem[IL10][this->getIndex()];//this->meanNeighborChem(IL10);

  float countfHA = this->countNeighborECM(fha);

  REGULATORS_T uregs_col{ meanTGF,  meanTNF,  meanIL1,  meanIL6, countfHA};
  REGULATORS_T dregs_col{ meanFGF,  meanIL8,  meanIL10};
  REGULATORS_T coefs_col{col_cTGF, col_cTNF, col_cIL1, col_cIL6, col_cfHA,
    col_cFGF, col_cIL8, col_cIL10};

  REGULATORS_T uregs_ela{ meanTGF};
  REGULATORS_T dregs_ela{ meanFGF,  meanIL1,  meanTNF};
  REGULATORS_T coefs_ela{ela_cTGF,
    ela_cFGF, ela_cIL1, ela_cTNF};

  REGULATORS_T uregs_hya{ meanTGF,  meanIL1,   meanIL6};
  REGULATORS_T dregs_hya{ meanIL8,  meanIL10,  meanFGF};
  REGULATORS_T coefs_hya{hya_cTGF, hya_cIL1,  hya_cIL6,
    hya_cIL8, hya_cIL10, hya_cFGF};

  // Make collagen monomers
  float newcol = this->make_coll_monomers(
      dI_targets[icol],
      uregs_col,
      dregs_col,
      coefs_col,
      col_ofst);

  // Make elastin monomers
  float newela = this->make_elas_monomers(
      dI_targets[iela],
      uregs_ela,
      dregs_ela,
      coefs_ela,
      ela_ofst);

  // Make hyaluronans
  float newhya = this->make_hyaluronans(
      dI_targets[ihya],
      uregs_hya,
      dregs_hya,
      coefs_hya,
      hya_ofst);

  if (newcol > 0.0f) {
  	this->move(dI_targets[icol][ix], dI_targets[icol][iy], dI_targets[icol][iz], read_t);
  } else if (newela > 0.0f) {
  	this->move(dI_targets[iela][ix], dI_targets[iela][iy], dI_targets[iela][iz], read_t);
  } else if (newhya > 0.0f) {
  	this->move(dI_targets[ihya][ix], dI_targets[ihya][iy], dI_targets[ihya][iz], read_t);
  } else {
  	wiggle_wound();
  }


  // TODO(NS): add print debuging section
  //  float newcoll = this->make_coll_monomers();
  //  float newelas = this->make_elas_monomers();
  //  float newha   = this->make_hyaluronans();

}


int Fibroblast::get_current_index(int &x, int &y, int &z, int &read_index)
{
  // Check if the location has been modified in this tick
  if (isModified(this->index)) {
      // If it has, work off of the intermediate value
      read_index = write_t;
  } else {
      // If it has NOT, work off of the original value
      read_index = read_t;
  }

  // Location of fibroblast in x,y,z dimensions of world.
  x = this->ix[read_index];
  y = this->iy[read_index];
  z = this->iz[read_index];
}

float Fibroblast::make_coll_monomers(
    int dI_target[3],
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset)
{
  int x, y, z, read_index;
  this->get_current_index(x, y, z, read_index);

  // Number of patches in x,y,z dimensions of world
  int nx = Agent::nx;
  int ny = Agent::ny;
  int nz = Agent::nz;

  // calculate stimulation (normalized)
  float max_stim = 2.0f;
  float norm_stim = this->calc_stim(uregs, dregs, coefs, offset, max_stim);
  float newMonomers = norm_stim * FIB_COL_PROD_RATE;

//  int dx = dI_target[0],
//      dy = dI_target[1],
//      dz = dI_target[2];

  int dx = 0,
      dy = 0,
      dz = 0;

  int target = (z+dz) * nx * ny + (y+dy) * nx + (x+dx);


	// DEBUG vis
//	Agent::agentWorldPtr->incECM(target, m_col, newMonomers/CONV_RATE_COL);
//	return 0;

  if (norm_stim == 0.0f) return 0.0f;


  Agent::agentECMPtr[target].addMonomerCol(newMonomers);


  return newMonomers;
}


float Fibroblast::make_elas_monomers(
    int dI_target[3],
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset)
{
  int x, y, z, read_index;
  this->get_current_index(x, y, z, read_index);

  // Number of patches in x,y,z dimensions of world
  int nx = Agent::nx;
  int ny = Agent::ny;
  int nz = Agent::nz;

  // calculate stimulation (normalized)
  float max_stim = 2.0f;
  float norm_stim = this->calc_stim(uregs, dregs, coefs, offset, max_stim);
  float newMonomers = norm_stim * FIB_ELA_PROD_RATE;

  int dx = dI_target[0],
      dy = dI_target[1],
      dz = dI_target[2];

//  this->move(dx, dy, dz, read_index);
  int target = (z+dz) * nx * ny + (y+dy) * nx + (x+dx);

	// DEBUG vis
//	Agent::agentWorldPtr->incECM(target, m_col, 1.0f);
  if (norm_stim == 0.0f) return 0.0f;



  Agent::agentECMPtr[target].addMonomerEla(newMonomers);
  // debug vis
#ifdef HUMAN_VF
  Agent::agentWorldPtr->incECM(target, m_ela, 1.0f);
#endif
  return newMonomers;
}


float Fibroblast::make_hyaluronans(
    int dI_target[3],
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset)
{
  int x, y, z, read_index;
  this->get_current_index(x, y, z, read_index);

  // Number of patches in x,y,z dimensions of world
  int nx = Agent::nx;
  int ny = Agent::ny;
  int nz = Agent::nz;

  // calculate stimulation (normalized)
  float max_stim = 2.0f;
  float norm_stim = this->calc_stim(uregs, dregs, coefs, offset, max_stim);
  float newhya = norm_stim * (FIB_HYA_PROD_RATE/CONV_RATE_HYA);

  int dx = dI_target[0],
      dy = dI_target[1],
      dz = dI_target[2];

//  this->move(dx, dy, dz, read_index);
  int target = (z+dz) * nx * ny + (y+dy) * nx + (x+dx);

	// DEBUG vis
//	Agent::agentWorldPtr->incECM(target, m_col, 1.0f);
  if (norm_stim == 0.0f) return 0.0f;


  // Default HA lifespan = 100 ticks
  // TODO(NS): ADD REF
  Agent::agentECMPtr[target].addHAs(Agent::agentWorldPtr->clock + 100, newhya);
  // debug vis
//  Agent::agentWorldPtr->incECM(target, m_hya, norm_stim*3.75f);

  return newhya;

}

void Fibroblast::hatchnewfibroblast(int number) {
  int newfibs = 0;
  // Location of fibroblast in x,y,z dimensions of the world
  int x = this->ix[read_t];
  int y = this->iy[read_t];
  int z = this->iz[read_t];
  // Number of patches in x,y,z dimensions of world
  int nx = Agent::nx;
  int ny = Agent::ny;
  int nz = Agent::nz;

  // Shuffle neighboring patches and go through them in a random order
#ifdef MODEL_3D
  const int nn = 27;  // number of neighbors + self
#else
  const int nn =  8;  // number of neighbors + self
#endif
  int nb[nn];
  for (int i = 0; i < nn; i++) {
      nb[i] = Agent::neighbor[i];
  }
  random_shuffle(&nb[0], &nb[nn-1]);
  for (int i = 0; i < nn && newfibs < number; i++) {
      // Distance away from target neighboring patch in x,y,z dimensions
      int dx = Agent::dX[nb[i]];
      int dy = Agent::dY[nb[i]];
      int dz = Agent::dZ[nb[i]];
      // Patch row major index of target neighboring patch
      int in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

      /* Try a new target neighboring patch if this one is not inside the world
       * dimensions, or is occupied, or is a capillary patch, or is an epithelial
       * patch. TODO(Kim): INSERT REF? (fibroblast can't hatch in epithelium or capillary) */
      if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
      int targetType = agentPatchPtr[in].type[read_t];
      //		if (agentPatchPtr[in].isOccupied() || targetType == capillary || targetType == epithelium) continue;

      if (targetType == capillary || targetType == epithelium) continue;

      /* If patch is unoccupied, craete new instance on new patch at
       * (x + dX, y + dY, z + dZ). setOccupied() sets the new patch as occupied and
       * returns true if the patch was already occupied */
      if (!Agent::agentPatchPtr[in].setOccupied())
        {
          // Create a new fibroblast at the valid target neighboring patch
          Fibroblast* newfibroblast = new Fibroblast(x + dx, y + dy, z + dz);
          newfibroblast->activate[write_t] = true;
          newfibs++;
          // Update target neighboring patch as occupied
          Agent::agentPatchPtr[in].occupiedby[write_t] = afibroblast;
#ifdef _OMP
          /* If executing OMP version, add the pointer to this new fibroblast to the
           * thread-local list first. WHWorld::UpdateCells() will take care of
           * putting it in the global list at the end */
          int tid = omp_get_thread_num();
          // DEBUG

          Agent::agentWorldPtr->localNewFibs[tid]->push_back(newfibroblast);
#else
          // If executing serial version, add the pointer to this new fibroblast to the global list right away
          Agent::agentWorldPtr->fibs.addData(newfibroblast, DEFAULT_TID);
#endif
        }
  }

}


void Fibroblast::die() {
  // DEBUG
  /*
#pragma omp critical
{
  if (this->activate[read_t] == false) (Agent::agentWorldPtr->deadfibs)++;
  else (Agent::agentWorldPtr->dead_afibs)++;
}
   */
  int in = this->index[read_t];
  Agent::agentPatchPtr[in].clearOccupied();
  Agent::agentPatchPtr[in].occupiedby[write_t] = nothing;
  this->alive[write_t] = false;
  this->life[write_t] = 0;
}

void Fibroblast::fibActivation() {

  int in = this->index[read_t];
  int target = this->index[write_t]; /* This assumes that after this function
  is called, no more move() would be called in the same tick -- should be taken cared by multiple moves using Agent::isModify() (? and the cell will not naturally die ?)*/
  if (this->activate[read_t] == false){ // && this->life[read_t] > 1) {
      // DEBUG
      /*
#pragma omp critical
{
   (Agent::agentWorldPtr->actfibs)++;
}
       */
      this->type[write_t] = afibroblast;
      this->activate[write_t] = true;
      this->color[write_t] = cafibroblast;
      Agent::agentPatchPtr[target].setOccupied();
      Agent::agentPatchPtr[target].occupiedby[write_t] = afibroblast;
  }
}

void Fibroblast::fibDeactivation() {
  int in = this->index[read_t];
  int target = this->index[write_t]; /* This assumes that after this function
  is called, no more move() would be called in the same tick */
  if (this->activate[read_t] == true&& this->life[read_t] > 1) {
      // DEBUG
      /*
#pragma omp critical
{
   (Agent::agentWorldPtr->deactfibs)++;
}
       */
      this->type[write_t] = fibroblast;
      this->activate[write_t] = false;
      this->color[write_t] = cfibroblast;
      Agent::agentPatchPtr[in].setOccupied();
      Agent::agentPatchPtr[in].occupiedby[write_t] = fibroblast;
  }
}

void Fibroblast::copyAndInitialize(Agent* original, int dx, int dy, int dz) {

  int tid = 0;
#ifdef _OMP
  // Get thread id in order to access the seed that belongs to this thread
  tid = omp_get_thread_num();
#endif
  int in = this->index[read_t];
  // Initializes location of new Fibroblast relative to original agent
  this->ix[write_t] = original->getX() + dx;
  this->iy[write_t] = original->getY() + dy;
  this->iz[write_t] = original->getZ() + dz;
  this->index[write_t] = this->ix[write_t] + this->iy[write_t]*Agent::nx + this->iz[write_t]*Agent::nx*Agent::ny;
  // Initializes new Fibroblast
  this->alive[read_t] = true;
#ifdef DISABLE_RAND
  this->life[write_t] = WHWorld::reportTick(0, 103);  // 0 corresponds to hours. TODO(Kim:) INSERT REF
#else
  this->life[write_t] = WHWorld::reportTick(0, 100 + rand_r(&(agentWorldPtr->seeds[tid]))%7);  // Unactivated fibroblasts live for 100 to 107 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
#endif
  this->activate[read_t] = false;
  this->color[read_t]= cfibroblast;
  this->size[read_t] = 2;
  this->type[read_t] = fibroblast;
  this->alive[write_t] = true;
  this->life[write_t] = this->life[read_t];
  this->activate[write_t] = false;
  this->color[write_t]= cfibroblast;
  this->size[write_t] = 2;
  this->type[write_t] = fibroblast;

  // TODO(NS): INSERT REF
  this->maxProlif[0] = 1;
  this->maxProlif[1] = 2;
  this->maxProlif[2] = 1;

  Fibroblast::numOfFibroblasts++;

  // Assigns new Fibroblast to this patch if it is unoccupied
  if (Agent::agentPatchPtr[in].isOccupied() == false) {
      Agent::agentPatchPtr[in].setOccupied();
      Agent::agentPatchPtr[in].occupiedby[write_t] = this->type[read_t];
  } else {
      cout << "error in hatching and initialization!!!" << dx << " " << dy << endl;
  }
}


