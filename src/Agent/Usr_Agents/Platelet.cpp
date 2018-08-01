/* 
 * File: Platelet.cpp
 *
 * File Contents: Contains the Platelet class.
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

#include "Platelet.h"
#include "../../enums.h"
#include <iostream>

using namespace std;
int Platelet::numOfPlatelets = 0;

float Platelet::cytokineSynthesis[3] = {0.1, 0.1, 0.5};


Platelet::Platelet() {
	this->color[read_t] = cplatelet;
	this->color[write_t] = cplatelet;
	cout << "Platelet is generated" << endl;
}

Platelet::Platelet(int x, int y, int z) {
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
//	this->life[write_t] = WHWorld::reportTick(rand_r(&(agentWorldPtr->seeds[tid]))%12 + 12, 0); // Platelets live for 12-23 hours in tissue. 0 corresponds to days. TODO(Kim): INSERT REF
	this->life[write_t] = WHWorld::reportTick(rand_r(&(agentWorldPtr->seeds[tid]))%25, 1); // Platelets live for 24-48 hours in tissue. 0 corresponds to days. TODO(Kim): INSERT REF
	this->color[write_t] = cplatelet;
	this->size[write_t] = 2;
	this->type[write_t] = platelet;

	this->ix[read_t] = x;
	this->iy[read_t] = y;
	this->iz[read_t] = z;
	this->index[read_t] = x + y*nx + z*nx*ny;

  /* In OMP version, we wait to add cells at the end of the tick,
   * whereas in serial version, cells are always added right away.
   * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
	this->alive[read_t] = true;
	this->life[read_t] = this->life[write_t];
#else
	this->alive[read_t] = false;
	this->life[read_t]  = 0;

#endif
	this->color[read_t] = cplatelet;
	this->size[read_t] = 2;
	this->type[read_t] = platelet;

	Platelet::numOfPlatelets++;
}

Platelet::Platelet(Patch* patchPtr) {
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
//	this->life[write_t] = WHWorld::reportTick(rand_r(&(agentWorldPtr->seeds[tid]))%12 + 12, 0); // Platelets live for 12-23 hours in tissue. 0 corresponds to days. TODO(Kim): INSERT REF
	this->life[write_t] = WHWorld::reportTick(rand_r(&(agentWorldPtr->seeds[tid]))%25, 1); // Platelets live for 24-48 hours in tissue. 0 corresponds to days. TODO(Kim): INSERT REF
	this->color[write_t] = cplatelet;
	this->size[write_t] = 2;
	this->type[write_t] = platelet;
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
	this->life[read_t]  = 0;

#endif
	this->color[read_t] = cplatelet;
	this->size[read_t] = 2;
	this->type[read_t] = platelet;
	Platelet::numOfPlatelets++;
}

Platelet::~Platelet() {
}

void Platelet::cellFunction() {

  int in = this->index[read_t];
  if ((this->life[read_t] > 0) && (this->alive[read_t] == true)) {
    /*************************************************************************
     * CHEMICAL SYNTHESIS                                                    *
     *************************************************************************/
  // Platelet makes TGF, IL1-beta, MMP8 TODO(Kim): INSERT REF?
#ifndef CALIBRATION
  float factor  = 1.0;//0.000001;
#ifdef RAT_VF
  float factorIL1 = 11.5;//115.0;//100.0;
#else
  float factorIL1 = 0.00115;
#endif
  float factorTGF = 0.000001;

  float TGFinc  = 0.1;
  float IL1inc  = 0.5*Agent::agentWorldPtr->baselineChem[MMP8];
  float MMP8inc = 0.1*Agent::agentWorldPtr->baselineChem[IL1beta];
  TGFinc  *= factorTGF;
  IL1inc  *= factorIL1;
  MMP8inc *= factor;
#else  // CALIBRATION
  float TGFinc  = Platelet::cytokineSynthesis[0];
  float IL1inc  = Platelet::cytokineSynthesis[1]*Agent::agentWorldPtr->baselineChem[IL1beta];
  float MMP8inc = Platelet::cytokineSynthesis[2]*Agent::agentWorldPtr->baselineChem[MMP8];

#endif  // CALIBRATION
#ifdef OPT_CHEM
  Agent::agentWorldPtr->WHWorldChem->incDchem(TGF,     in, TGFinc);
  Agent::agentWorldPtr->WHWorldChem->incDchem(MMP8,    in, MMP8inc);
  Agent::agentWorldPtr->WHWorldChem->incDchem(IL1beta, in, IL1inc);
#else		// OPT_CHEM
  Agent::agentWorldPtr->WHWorldChem->dTGF[in]     += TGFinc;
  Agent::agentWorldPtr->WHWorldChem->dMMP8[in]    += MMP8inc;
  Agent::agentWorldPtr->WHWorldChem->dIL1beta[in] += IL1inc;
#endif		// OPT_CHEM
#ifdef PRINT_SECRETION
  int x = this->ix[read_t];
  int y = this->iy[read_t];
  int z = this->iz[read_t];
  printCytRelease(0, TGF, x, y, z, Platelet::cytokineSynthesis[0]);
  printCytRelease(0, IL1beta, x, y, z, Platelet::cytokineSynthesis[1]*Agent::agentWorldPtr->baselineChem[IL1beta]);
  printCytRelease(0, MMP8, x, y, z, Platelet::cytokineSynthesis[2]*Agent::agentWorldPtr->baselineChem[MMP8]);
#endif  // PRINT_SECRETION
    /*************************************************************************
     * DEATH                                                                 *
     *************************************************************************/
    // Platelets can die naturally
		this->life[write_t] = this->life[read_t] - 1;
		if (this->life[read_t] <= 0) {
			this->die();
		}

	} else if (this->life[read_t] == 0 || this->alive[read_t] == false) {
    // Already dead, do nothing
	} else {
		cout << "an error in dying is made from platelet!!!" << endl;
	}
}

void Platelet::die() {
	int in = this->index[read_t];
	Agent::agentPatchPtr[in].clearOccupied();
	Agent::agentPatchPtr[in].occupiedby[write_t] = nothing;
	this->alive[write_t] = false;
	this->life[write_t] = 0;
}
