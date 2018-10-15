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
float Neutrophil::cytokineSynthesis[21] = {0.0f};
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
  this->aneu_produce_cytokines();

  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  // Activated neutrophils might die once the damage is cleared
  int in = this->index[read_t];
  float patchIL10 = this->agentWorldPtr->WHWorldChem->pIL10[in];
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[pdamage];
#ifndef CALIBRATION
	if (totaldamage == 0 && rollDice(10))  // TODO(Kim): INSERT REFS?
#else  // CALIBRATION
	if (totaldamage == 0 && rollDice(Neutrophil::death[0]))  // TODO(Kim): INSERT REFS?
#endif  // CALIBRATION
	{
		this->die(); 
		return;  // Return immediately after cell's death
	}
  // Activated neutrophils can die naturally
#ifndef CALIBRATION
	this->life[write_t] = this->life[read_t] - 1 - 0.01*patchIL10;  // TODO(Kim): INSERT REFS?
#else  // CALIBRATION
	this->life[write_t] = this->life[read_t] - 1 - Neutrophil::death[1]*patchIL10;  // TODO(Kim): INSERT REFS?
#endif  // CALIBRATION
	if (this->life[read_t] <= 0)
	{
    this->die();
  }

  /*************************************************************************
   * SIGNAL TRANSDUCTION                                                   *
   *************************************************************************/                                           
	Agent::agentWorldPtr->highTNFdamage = true;
}


void Neutrophil::aneu_produce_cytokines() {
		// Calculates chemical gradients and patch chemical concentrations
		int in = this->index[read_t];

		float patchTNF  = this->agentWorldPtr->WHWorldChem->pTNF[in];
		float patchTGF  = this->agentWorldPtr->WHWorldChem->pTGF[in];
		float patchIL8  = this->agentWorldPtr->WHWorldChem->pIL8[in];
		float patchIL10 = this->agentWorldPtr->WHWorldChem->pIL10[in];

		/* Activated neutrophils synthesize new cytokines in quantities dependent on
		 * the vocal treatment type. TODO(Kim): INSERT REFS? */

		float TNFinc  = 0.0f;
		float MMP8inc = 0.0f;
		float IL1inc  = 0.0f;

		float tnf_cNum  = 0.0f;
		float tnf_cDen  = 0.0f;
		float tnf_cTGF  = 0.0f;
		float tnf_cIL10 = 0.0f;
		float tnf_ofst  = 0.0f;

		float il1_cTNF  = 0.0f;	// up
		float il1_cIL8  = 0.0f;	// up
		float il1_cIL10 = 0.0f;	// down
		float il1_ofst  = 0.0f;

		float mmp8_cNum = 0.0f;
		float mmp8_cTNF = 0.0f;
		float mmp8_cDen = 0.0f;
		float mmp8_cTGF = 0.0f;
		float mmp8_ofst = 0.0f;

		tnf_ofst  = Agent::baseline[TNF];//0.0f;
		il1_ofst  = Agent::baseline[IL1beta];
		mmp8_ofst = Agent::baseline[MMP8];//0.0f;

		// Setting up, regulators, coefficients and offsets
#ifndef CALIBRATION
		float factor = 0.0001;

		if (this->agentWorldPtr->treatmentOption == voicerest) {
			tnf_cNum  = factor*1.0f;
			tnf_cDen  = 1.0f;
			tnf_cTGF  = 1.0f;
			tnf_cIL10 = 1.0f;

			mmp8_cNum  = factor*250.0f;
			mmp8_cTNF  = factor*1.0f;
			mmp8_cDen  = 1.0f;
			mmp8_cTGF  = 1.0f;
		} else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
			tnf_cNum  = factor*20.0f;
			tnf_cDen  = 1.0f;
			tnf_cTGF  = 1.0f;
			tnf_cIL10 = 1.0f;

			mmp8_cNum  = factor*10.0f;
			mmp8_cTNF  = factor*2.0f;
			mmp8_cDen  = 1.0f;
			mmp8_cTGF  = 0.5f;
		} else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
			tnf_cNum  = factor*1.0f;
			tnf_cDen  = 1.0f;
			tnf_cTGF  = 1.0f;
			tnf_cIL10 = 1.0f;

			mmp8_cNum  = factor*1500.0f;
			mmp8_cTNF  = factor*45.0f;
			mmp8_cDen  = 1.0f;
			mmp8_cTGF  = 1.0f;
		}

		il1_cTNF  = 0.1f;
		il1_cIL8  = 0.001f;
		il1_cIL10 = 10000.0f;

		REGULATORS_T tnf_ucoefs{tnf_cNum};
		REGULATORS_T tnf_dcoefs{tnf_cDen, tnf_cTGF, tnf_cIL10};
		REGULATORS_T tnf_uregs{1.0f};
		REGULATORS_T tnf_dregs{1.0f, patchTGF, patchIL10};

		REGULATORS_T mmp8_ucoefs{mmp8_cNum, mmp8_cTNF};
		REGULATORS_T mmp8_dcoefs{mmp8_cDen, mmp8_cTGF};
		REGULATORS_T mmp8_uregs{1.0f, patchTNF};
		REGULATORS_T mmp8_dregs{1.0f, patchTGF};

#else
		// TODO: Reorder parameters

		if (this->agentWorldPtr->treatmentOption == voicerest) {
			tnf_cNum  = 1.0f;
			tnf_cTGF  = Neutrophil::cytokineSynthesis[0];
			tnf_cIL10 = Neutrophil::cytokineSynthesis[1];

			mmp8_cTNF  = Neutrophil::cytokineSynthesis[2];
			mmp8_cTGF  = Neutrophil::cytokineSynthesis[3];
		} else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
			tnf_cNum  = 1.0f;
			tnf_cTGF  = Neutrophil::cytokineSynthesis[4];
			tnf_cIL10 = Neutrophil::cytokineSynthesis[5];

			mmp8_cTNF  = Neutrophil::cytokineSynthesis[6];
			mmp8_cTGF  = Neutrophil::cytokineSynthesis[7];
		} else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
			tnf_cNum  = 1.0f;
			tnf_cTGF  = Neutrophil::cytokineSynthesis[8];
			tnf_cIL10 = Neutrophil::cytokineSynthesis[9];

			mmp8_cTNF  = Neutrophil::cytokineSynthesis[10];
			mmp8_cTGF  = Neutrophil::cytokineSynthesis[11];
		}

		il1_cTNF  = Neutrophil::cytokineSynthesis[12];
		il1_cIL8  = Neutrophil::cytokineSynthesis[13];
		il1_cIL10 = Neutrophil::cytokineSynthesis[14];

		REGULATORS_T tnf_ucoefs{tnf_cNum};
		REGULATORS_T tnf_dcoefs{tnf_cTGF, tnf_cIL10};
		REGULATORS_T tnf_uregs{1.0f};
		REGULATORS_T tnf_dregs{patchTGF, patchIL10};

		REGULATORS_T mmp8_ucoefs{mmp8_cTNF};
		REGULATORS_T mmp8_dcoefs{mmp8_cTGF};
		REGULATORS_T mmp8_uregs{patchTNF};
		REGULATORS_T mmp8_dregs{patchTGF};

#endif

		REGULATORS_T il1_ucoefs{il1_cTNF,il1_cIL8 };
		REGULATORS_T il1_dcoefs{il1_cIL10};
		REGULATORS_T il1_uregs{patchTNF, patchIL8};
		REGULATORS_T il1_dregs{patchIL10};

		// Call cytokine production functions
		TNFinc  = produce_tnf ( tnf_uregs,  tnf_dregs,  tnf_ucoefs,  tnf_dcoefs,  tnf_ofst);
		MMP8inc = produce_mmp8(mmp8_uregs, mmp8_dregs, mmp8_ucoefs, mmp8_dcoefs, mmp8_ofst);
//		IL1inc  = produce_il1 ( il1_uregs,  il1_dregs,  il1_ucoefs,  il1_dcoefs,  il1_ofst);

#ifdef PRINT_SECRETION
		int x = this->ix[read_t];
		int y = this->iy[read_t];
		int z = this->iz[read_t];
		printCytRelease(3, TNF,     x, y, z, TNFinc);
		printCytRelease(3, MMP8, 	  x, y, z, MMP8inc);
		printCytRelease(3, IL1beta, x, y, z, IL1inc);
#endif  // PRINT_SECRETION

	}

