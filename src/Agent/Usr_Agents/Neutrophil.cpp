/* 
 * File: Neutrophil.cpp
 *
 * File Contents: Contains the Neutrophil class.
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
 *****************************************************************************/ // TODO(Kim): Update the file comment once we figure out the copyright issues

#include "Neutrophil.h"
#include "../../enums.h"
#include <iostream>

using namespace std;
int Neutrophil::numOfNeutrophil = 0; 
float Neutrophil::cytokineSynthesis[21] = {1, 1, 20, 1, 1, 1, 1, 250, 1, 1, 1, 1, 10, 2, 1, 0.5, 15, 100, 3, 1, 1};
float Neutrophil::activation[4] = {0.1, 0, 25, 10};
float Neutrophil::death[2] = {10, 0.01};


Neutrophil::Neutrophil() {
}

Neutrophil::Neutrophil(Patch* patchPtr) {
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
	this->life[write_t] = WHWorld::reportTick(9 + rand_r(&(agentWorldPtr->seeds[tid]))%55, 0);
 /* Unactivated neutrophils live for 12-20 hours 
  * (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1748424/). 
  * Second source: 5-90 hours
  * (Tak T, Tesselaar K, Pillay J, Borghans JA, Koenderman L (2013).
  *   "What's your age again? Determination of human neutrophil half-lives revisited".
  *   Journal of Leukocyte Biology. 94 (4): 595–601.)
  * 0 corresponds to days. */
	this->activate[write_t] = false;
	this->color[write_t] = cneutrophil;
	this->size[write_t] = 2;
	this->type[write_t] = neutrophil;
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
	this->color[read_t] = cneutrophil;
	this->size[read_t] = 2;
	this->type[read_t] = neutrophil;

	Neutrophil::numOfNeutrophil++;
}

Neutrophil::~Neutrophil() {
}

void Neutrophil::cellFunction() {
	if (this->alive[read_t] == false) return;
	if (this->activate[read_t] == false) {
		neu_cellFunction();
	} else {
		aneu_cellFunction();
	}
}

void Neutrophil::neuActivation() {
	int in = this->index[read_t];
	int target = this->index[write_t]; /* This assumes that after this function
  is called, no more move() would be called in the same tick and will not die naturally*/

	if (this->activate[read_t] == false && this->life[read_t] > 1) {
		this->type[write_t] = aneutrophil;
		this->color[write_t] = caneutrophil;
		this->size[write_t] = 2;
		this->activate[write_t] = true;
		int tid = 0;
#ifdef _OMP
    // Get thread id in order to access the seed that belongs to this thread
		tid = omp_get_thread_num();
#endif
		this->life[write_t] = WHWorld::reportTick(0, 1 + (rand_r(&(agentWorldPtr->seeds[tid]))%2)); 
    /* Activated neutrophils live for 1-2 days
     * (http://en.wikipedia.org/wiki/Neutrophil#Lifespan). 
     * 0 corresponds to hours. */

		Agent::agentPatchPtr[target].setOccupied();
		Agent::agentPatchPtr[target].occupiedby[write_t] = aneutrophil;
	}
}

void Neutrophil::neuSniff() {
	if (this->moveToHighestChem(NEUgrad)) {
		return;
	} else {
		this->wiggle();
	}
}

void Neutrophil::die() {
// DEBUG
   /*
 #pragma omp critical
 {
   (Agent::agentWorldPtr->deadneus)++;
 }
*/
	int in = this->index[read_t];
	Agent::agentPatchPtr[in].clearOccupied();
	Agent::agentPatchPtr[in].occupiedby[write_t] = nothing;
	this->alive[write_t] = false;
	this->life[write_t] = 0;
}

void Neutrophil::neu_cellFunction() {
  int in = this->index[read_t];
  int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[pdamage];
#ifdef OPT_CHEM
  float patchTNF = this->agentWorldPtr->WHWorldChem->getPchem(TNF, in);
  float patchIL10 = this->agentWorldPtr->WHWorldChem->getPchem(IL10, in);
#else		// OPTC_CHEM
  float patchTNF = this->agentWorldPtr->WHWorldChem->pTNF[in];
  float patchIL10 = this->agentWorldPtr->WHWorldChem->pIL10[in];
#endif	 	// OPT_CHEM
  /* Unactivated neutrophils only move along their preferred gradient if there 
   * is damage. TODO(Kim): Insert ref? */
  if (totaldamage == 0 ) {
//		cout << "	Neu wiggle() totaldam = 0" << endl;
    this->wiggle();
  } else {
    /*************************************************************************
     * MOVEMENT                                                              *
     *************************************************************************/
    this->neuSniff();

    /*************************************************************************
     * ACTIVATION                                                            *
     *************************************************************************/
		/* An unactivated neutrophil can be activated if it is in the damage zone.
     * TODO(Kim): Insert ref? */
    if (Agent::agentPatchPtr[in].inDamzone == 1) {
#ifndef CALIBRATION
      float activationFactor = 5.0f;//0.1;
      int chance1 = 15;
      int chance2 = 5;
      if (patchTNF >= patchIL10*activationFactor || patchTNF > 0 && rollDice(chance1) || rollDice(chance2))  // TODO(Kim): INSERT REF?
#else  // CALIBRATION
      if (patchTNF >= patchIL10*Neutrophil::activation[0] || patchTNF > Neutrophil::activation[1] && rollDice(Neutrophil::activation[2]) || rollDice(Neutrophil::activation[3]))  // TODO(Kim): INSERT REF?
#endif  // CALIBRATION
      {
        this->neuActivation();
      }
    }
	
  }

  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  // Unactivated neutrophils can die naturally
  this->life[write_t] = this->life[read_t] - 1;
  if (this->life[read_t] <= 0) {
    this->die();
  }
}

void Neutrophil::aneu_cellFunction() {
  /*************************************************************************
   * MOVEMENT                                                              *
   *************************************************************************/
  // Activated neutrophils always move along their preferred gradient. TODO(Kim): INSERT REF?
  this->neuSniff();

  /*************************************************************************
   * CHEMICAL SYNTHESIS                                                    *
   *************************************************************************/
  int in = this->index[read_t];

#ifdef OPT_CHEM 
  float patchTNF = this->agentWorldPtr->WHWorldChem->getPchem(TNF, in);
  float patchTGF = this->agentWorldPtr->WHWorldChem->getPchem(TGF, in);
  float patchIL10 = this->agentWorldPtr->WHWorldChem->getPchem(IL10, in);
#else		// OPT_CHEM
  float patchTNF = this->agentWorldPtr->WHWorldChem->pTNF[in];
  float patchTGF = this->agentWorldPtr->WHWorldChem->pTGF[in];
  float patchIL10 = this->agentWorldPtr->WHWorldChem->pIL10[in];
#endif		// OPT_CHEM
  /* Activated neutrophils synthesize new cytokines in quantities dependent on
   * the vocal treatment type. TODO(Kim): INSERT REFS? */
  float factor = 0.0001;

  float TNFinc = 0.0;
  float MMP8inc = 0.0;
#ifndef CALIBRATION
  if (this->agentWorldPtr->treatmentOption == voicerest) {
    TNFinc  = 1/(1 + patchTGF + patchIL10);
    MMP8inc = (250 + patchTNF)/(1 + patchTGF);
  } else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
    TNFinc  = 20/(1 + patchTGF + patchIL10);
    MMP8inc = (10 + patchTNF*2)/(1 + patchTGF*0.5);
  } else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
    TNFinc  = 1/(1 + patchTGF + patchIL10);
    MMP8inc = 15*(100 + patchTNF*3)/(1 + patchTGF);
  }

  TNFinc  *= factor;
  MMP8inc *= factor;
#else  // CALIBRATION
  if (this->agentWorldPtr->treatmentOption == voicerest) {
    TNFinc  = Neutrophil::cytokineSynthesis[0]/
             (Neutrophil::cytokineSynthesis[1] + patchTGF + patchIL10);
    MMP8inc = Neutrophil::cytokineSynthesis[6]*
             (Neutrophil::cytokineSynthesis[7] + patchTNF*Neutrophil::cytokineSynthesis[8])/
             (Neutrophil::cytokineSynthesis[9] + patchTGF*Neutrophil::cytokineSynthesis[10]);

  } else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
    TNFinc  = Neutrophil::cytokineSynthesis[2]/
             (Neutrophil::cytokineSynthesis[3] + patchTGF + patchIL10);
    MMP8inc = Neutrophil::cytokineSynthesis[11]*
             (Neutrophil::cytokineSynthesis[12] + patchTNF*Neutrophil::cytokineSynthesis[13])/
             (Neutrophil::cytokineSynthesis[14] + patchTGF*Neutrophil::cytokineSynthesis[15]);

  } else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
    TNFinc  = Neutrophil::cytokineSynthesis[4]/
             (Neutrophil::cytokineSynthesis[5] + patchTGF + patchIL10);
    MMP8inc = Neutrophil::cytokineSynthesis[16]*
             (Neutrophil::cytokineSynthesis[17] + patchTNF*Neutrophil::cytokineSynthesis[18])/
             (Neutrophil::cytokineSynthesis[19] + patchTGF*Neutrophil::cytokineSynthesis[20]);
  }
#endif  // CALIBRATION
#ifdef OPT_CHEM
  this->agentWorldPtr->WHWorldChem->incDchem(TNF,  in, TNFinc);
  this->agentWorldPtr->WHWorldChem->incDchem(MMP8, in, MMP8inc);
#else		// OPT_CHEM
  (this->agentWorldPtr->WHWorldChem->dTNF[in])  += TNFinc;
  (this->agentWorldPtr->WHWorldChem->dMMP8[in]) += MMP8inc;
#endif		// OPT_CHEM

#ifdef PRINT_SECRETION
  int x = this->ix[read_t];
  int y = this->iy[read_t];
  int z = this->iz[read_t];
  printCytRelease(1, TNF,  x, y, z, TNFinc);
  printCytRelease(1, MMP8, x, y, z, MMP8inc);
#endif  // PRINT_SECRETION
  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  // Activated neutrophils might die once the damage is cleared
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[pdamage];
#ifndef CALIBRATION
	if (totaldamage == 0 && rollDice(10)) {  // TODO(Kim): INSERT REFS?
#else  // CALIBRATION
	if (totaldamage == 0 && rollDice(Neutrophil::death[0])) {  // TODO(Kim): INSERT REFS?
#endif  // CALIBRATION
		this->die(); 
		return;  // Return immediately after cell's death
	}
  // Activated neutrophils can die naturally
#ifndef CALIBRATION
	this->life[write_t] = this->life[read_t] - 1 - 0.01*patchIL10;  // TODO(Kim): INSERT REFS?
#else  // CALIBRATION
	this->life[write_t] = this->life[read_t] - 1 - Neutrophil::death[1]*patchIL10;  // TODO(Kim): INSERT REFS?
#endif  // CALIBRATION
	if (this->life[read_t] <= 0) {
    this->die();
  }

  /*************************************************************************
   * SIGNAL TRANSDUCTION                                                   *
   *************************************************************************/                                           
	Agent::agentWorldPtr->highTNFdamage = true;
}
