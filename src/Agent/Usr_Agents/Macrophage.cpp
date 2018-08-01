/* 
 * File: Macrophage.cpp
 *
 * File Contents: Contains the Macrophage class.
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

#include "Macrophage.h"
#include "../../enums.h"
#include <iostream>

using namespace std;
int Macrophage::numOfMacrophage = 0;

float Macrophage::cytokineSynthesis[82] = {0.5, 1, 1, 1, 0.5, 2, 1, 0.1, 1, 0.1, 1, 1, 5, 0.1, 1, 1, 1, 1, 2, 1, 15, 15, 1, 0.5, 1, 1, 0.5, 1, 13, 1, 1, 1, 0.01, 1, 1, 1, 0.5, 1, 10, 0.5, 1, 1, 1, 1, 10, 1, 1, 4, 1, 0.5, 2, 1, 10, 5, 5, 1, 1, 100, 100, 0.1, 1, 0.5, 10, 1, 7, 0.1, 1, 0.5, 0.01, 1, 1, 0.1, 0.01, 1, 4, 0.001, 0.001, 0.05, 1, 1, 0.0005, 0.0005};
float Macrophage::activation[5] = {0.1, 0, 25, 10, 3};


Macrophage::Macrophage() {
}

Macrophage::Macrophage(int x, int y, int z, int bloodORtiss) {
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
	this->activate[write_t] = false;
	this->color[write_t] = cmacrophage;
	this->size[write_t] = 2;
	this->type[write_t] = macrophag;
	this->bloodORtissue = bloodORtiss;

  /* In OMP version, we wait to add cells at the end of the tick,
   * whereas in serial version, cells are always added right away.
   * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
	this->ix[read_t] = x;
	this->iy[read_t] = y;
	this->iz[read_t] = z;
	this->index[read_t] = x + y*nx + z*nx*ny;
	this->alive[read_t] = true;
	this->activate[read_t] = false;
	this->color[read_t] = cmacrophage;
	this->size[read_t] = 2;
	this->type[read_t] = macrophag;
#else
	this->ix[read_t] = x;
	this->iy[read_t] = y;
	this->iz[read_t] = z;
	this->index[read_t] = x + y*nx + z*nx*ny;
	this->alive[read_t] = false;
	this->activate[read_t] = false;
	this->color[read_t] = cmacrophage;
	this->size[read_t] = 2;
	this->type[read_t] = macrophag;
#endif

	if (bloodORtissue == blood) {

#ifdef DISABLE_RAND
		this->life[write_t] = WHWorld::reportTick(39, 0);  // Unactivated macrophages live for 8 to 69 hours in blood. 0 corresponds to days. TODO(Kim): INSERT REF
#else
		this->life[write_t] = WHWorld::reportTick(8 + rand_r(&(agentWorldPtr->seeds[tid]))%62, 0);  // Unactivated macrophages live for 8 to 69 hours in blood. 0 corresponds to days. TODO(Kim): INSERT REF
#endif
    /* In OMP version, we wait to add cells at the end of the tick,
     * whereas in serial version, cells are always added right away.
     * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
		this->life[read_t] = this->life[write_t];
#else
		this->life[read_t] = 0;
#endif

	} else{
#ifdef DISABLE_RAND
		this->life[write_t] = WHWorld::reportTick(0, 90);
#else
		this->life[write_t] = WHWorld::reportTick(0, 60 + rand_r(&(agentWorldPtr->seeds[tid]))%60);
#endif
    /* Unactivated macrophages live for 60 to 119 days in tissue 
     * (http://www.nanomedicine.com/NMIIA/15.4.3.1.htm). 
     * 0 corresponds to hours.*/

    /* In OMP version, we wait to add cells at the end of the tick,
     * whereas in serial version, cells are always added right away.
     * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
        	this->life[read_t] = this->life[write_t];
#else
		this->life[read_t] = 0;
#endif
	}
	Macrophage::numOfMacrophage++;
}

Macrophage::Macrophage(Patch* patchPtr, int bloodORtiss) {
  int tid = 0;
#ifdef _OMP
// Get thread id in order to access the seed that belongs to this thread
  tid = omp_get_thread_num();
#endif
        this->ix[write_t] = patchPtr->indice[0];
	this->iy[write_t] = patchPtr->indice[1];
	this->iz[write_t] = patchPtr->indice[2];
	this->index[write_t] = this->ix[write_t] + this->iy[write_t]*nx + this->iz[write_t]*nx*ny;
	this->alive[write_t] = true;
	this->activate[write_t] = false;
	this->color[write_t] = cmacrophage;
	this->size[write_t] = 2;
	this->type[write_t] = macrophag;
	this->bloodORtissue = bloodORtiss;
	this->ix[read_t] = patchPtr->indice[0];
	this->iy[read_t] = patchPtr->indice[1];
	this->iz[read_t] = patchPtr->indice[2];
	this->index[read_t] = this->ix[read_t] + this->iy[read_t]*nx + this->iz[read_t]*nx*ny;

  /* In OMP version, we wait to add cells at the end of the tick,
   * whereas in serial version, cells are always added right away.
   * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
	this->alive[read_t] = true;
#else
	this->alive[read_t] = false;
#endif
	this->activate[read_t] = false;
	this->color[read_t] = cmacrophage;
	this->size[read_t] = 2;
	this->type[read_t] = macrophag;


	if (bloodORtissue == blood) {
#ifdef DISABLE_RAND
		this->life[write_t] = WHWorld::reportTick(39, 0);  // Unactivated macrophages live for 8 to 69 hours in blood. 0 corresponds to days. TODO(Kim): INSERT REF
#else
		this->life[write_t] = WHWorld::reportTick(8 + rand_r(&(agentWorldPtr->seeds[tid]))%62, 0);  // Unactivated macrophages live for 8 to 69 hours in blood. 0 corresponds to days. TODO(Kim): INSERT REF
#endif
    /* In OMP version, we wait to add cells at the end of the tick,
     * whereas in serial version, cells are always added right away.
     * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
		this->life[read_t] = this->life[write_t];
#else
		this->life[read_t] = 0;
#endif

	} else {
#ifdef DISABLE_RAND
		this->life[write_t] = WHWorld::reportTick(0, 90);
#else
		this->life[write_t] = WHWorld::reportTick(0, 60 + rand_r(&(agentWorldPtr->seeds[tid]))%60);
#endif
    /* Unactivated macrophages live for 60 to 119 days in tissue 
     * (http://www.nanomedicine.com/NMIIA/15.4.3.1.htm). 
     * 0 corresponds to hours. */

    /* In OMP version, we wait to add cells at the end of the tick,
     * whereas in serial version, cells are always added right away.
     * Thus, cells, when added in OMP, should be alive right away. */
#ifdef _OMP
		this->life[read_t] = this->life[write_t];
#else
		this->life[read_t] = 0;
#endif
	}
	Macrophage::numOfMacrophage++;
}

Macrophage::~Macrophage() {
}

void Macrophage::cellFunction() {
	if (this->alive[read_t] == false) return;
	if (this->activate[read_t] == false) {
		this->mac_cellFunction();
	} else {
		this->activatedmac_cellFunction();
	}
}

void Macrophage::mac_cellFunction() {
	int in = this->index[read_t];
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[pdamage];
#ifdef OPT_CHEM
	float patchTNF     = this->agentWorldPtr-> WHWorldChem->getPchem(TNF, in);
	float patchIL1beta = this->agentWorldPtr-> WHWorldChem->getPchem(IL1beta, in);
	float patchIL10    = this->agentWorldPtr-> WHWorldChem->getPchem(IL10, in);
#else		// OPT_CHEM
	float patchTNF     = this->agentWorldPtr->WHWorldChem->pTNF[in];
	float patchIL1beta = this->agentWorldPtr->WHWorldChem->pIL1beta[in];
	float patchIL10    = this->agentWorldPtr->WHWorldChem->pIL10[in];
#endif		// OPT_CHEM


    /*************************************************************************
     * MOVEMENT                                                              *
     *************************************************************************/

		this->macSniff();
//		this->wiggle();

    /*************************************************************************
     * ACTIVATION                                                            *
     *************************************************************************/
		/* An unactivated macrophage can be activated if it is in the damage zone.
     * TODO(Kim): Insert ref? */
#ifndef CALIBRATION
		float threshold_fact = 0.1;
		int chance2 = 25;//25;//25;
		int chance3 = 10;//10;
#else
		float threshold_fact = Macrophage::activation[0];
		int chance2 = Macrophage::activation[2];
		int chance3 = Macrophage::activation[3];
#endif

		if (this->agentPatchPtr[in].inDamzone) {
			if ((patchTNF + patchIL1beta > patchIL10*threshold_fact)) {  // TODO(Kim): Insert ref?
				this->macActivation();
			} else if (patchTNF + patchIL1beta > 0 && rollDice(chance2)) {  // TODO(Kim): Insert ref?
				this->macActivation();
			} else if (rollDice(chance3)) {
				this->macActivation();
			}
		}


  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  // Unactivated macrophages can die naturally
	if (this->life[read_t] != this->life[write_t])      // been modified in the current tick
	{                                                   // for example, activation life
	    this->life[write_t] = this->life[write_t] - 1;
	} else {
	    this->life[write_t] = this->life[read_t] - 1;
	}
	if (this->life[read_t] <=0 ) { 
    		this->die();
  	}
}

void Macrophage::macSniff() {
	if (this->moveToHighestChem(MACgrad) == true) {
		return;
	} else {
		this->wiggle();
	}
}

void Macrophage::macActivation() {
	int in = this->index[read_t];
	int target = this->index[write_t]; /* This assumes that after this function
  is called, no more move() would be called in the same tick and cell will not die naturally*/

	if (this->activate[read_t] == false && this->life[read_t] > 1) {
		this->type[write_t] = amacrophag;
		this->color[write_t] = camacrophage;
		this->activate[write_t] = true;
		int tid = 0;
#ifdef _OMP
    // Get thread id in order to access the seed that belongs to this thread
		tid = omp_get_thread_num();
#endif
	this->life[write_t] = WHWorld::reportTick(0, 2) +
	        (rand_r(&(agentWorldPtr->seeds[tid])) % WHWorld::reportTick(0, 3));
//		this->life[write_t] = WHWorld::reportTick(0, 2 + rand_r(&(agentWorldPtr->seeds[tid]))%3);
    /* Activated macrophages live for 2-4 days 
     * (http://rumi.gdcb.iastate.edu/wiki/ThackerMeeting20060124).
     * 0 corresponds to hours. */

		Agent::agentPatchPtr[target].setOccupied();
		Agent::agentPatchPtr[target].occupiedby[write_t] = amacrophag;
	}
}

void Macrophage::macDeactivation() {
	int in = this->index[read_t];
	int target = this->index[write_t]; /* This assumes that after this function
  is called, no more move() would be called in the same tick and will not die naturally*/
	if (this->activate[read_t] == true && this->life[read_t] > 1) {
		this->type[write_t] = macrophag;
		this->color[write_t] = camacrophage;
		this->activate[write_t] = false;
		Agent::agentPatchPtr[in].setOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = macrophag;
	}
}

void Macrophage::die() {
	int in = this->index[read_t];
	Agent::agentPatchPtr[in].clearOccupied();
	Agent::agentPatchPtr[in].occupiedby[write_t] = nothing;
	this->alive[write_t] = false;
	this->life[write_t] = 0;
}

void Macrophage::activatedmac_cellFunction() {
  /*************************************************************************
   * MOVEMENT                                                              *
   *************************************************************************/
	// Activated macrophages always move along their preferred gradient. TODO(Kim): Insert ref?
	this->macSniff();


  /*************************************************************************
   * CHEMICAL SYNTHESIS                                                    *
   *************************************************************************/
	// Calculates neighboring ECM and patch chemical concentrations
	int in = this->index[read_t];
#ifndef CALIBRATION
	int rvis = Agent::agentWorldPtr->RVIS;
	int rvvs = Agent::agentWorldPtr->RVVS;
	int ssis = Agent::agentWorldPtr->SSIS;
	int ssvs = Agent::agentWorldPtr->SSVS;
#else  // CALIBRATION
	float rvis = Agent::agentWorldPtr->RVIS;
	float rvvs = Agent::agentWorldPtr->RVVS;
  float ssis = Agent::agentWorldPtr->SSIS;
	float ssvs = Agent::agentWorldPtr->SSVS;
#endif  // CALIBRATION
	int countfHA = this->countNeighborECM(fha);
	float bTNF = (Agent::agentWorldPtr)->baselineChem[TNF];
	float bIL1beta = (Agent::agentWorldPtr)->baselineChem[IL1beta];
	float bIL6 = (Agent::agentWorldPtr)->baselineChem[IL6];
	float bIL8 = (Agent::agentWorldPtr)->baselineChem[IL8];
	float bIL10 = (Agent::agentWorldPtr)->baselineChem[IL10];

#ifdef OPT_CHEM
	float patchTNF     = this->agentWorldPtr-> WHWorldChem->getPchem(TNF,in);
	float patchTGF     = this->agentWorldPtr-> WHWorldChem->getPchem(TGF, in);
	float patchIL1beta = this->agentWorldPtr-> WHWorldChem->getPchem(IL1beta, in);
	float patchIL6     = this->agentWorldPtr-> WHWorldChem->getPchem(IL6, in);
	float patchIL8     = this->agentWorldPtr-> WHWorldChem->getPchem(IL8, in);
	float patchIL10    = this->agentWorldPtr-> WHWorldChem->getPchem(IL10, in);
#else		// OPT_CHEM
	float patchTNF     = this->agentWorldPtr->WHWorldChem->pTNF[in];
	float patchTGF     = this->agentWorldPtr->WHWorldChem->pTGF[in];
	float patchIL1beta = this->agentWorldPtr->WHWorldChem->pIL1beta[in];
	float patchIL6     = this->agentWorldPtr->WHWorldChem->pIL6[in];
	float patchIL8     = this->agentWorldPtr->WHWorldChem->pIL8[in];
	float patchIL10    = this->agentWorldPtr->WHWorldChem->pIL10[in];
#endif		// OPT_CHEM
 	//Change in chemicals due to cells. TODO(Kim): INSERT REF?
  float factor = 0.0005;//0.1;//0.0001;
  float factorIL1 = 0.0001;//0.01;//0.097;
  float factorTNF = 0.000001;//0.85;//0.4;//0.0001;
  float factorTNF_fHA = 0.1;//0.1;//100.0;
  float factorIL1_fHA = 1.0;
  float factorIL8_fHA = 1.0;
  float factorIL10 = 0.1;//0.1;
  float factorFGF = 0.0001;
  float factorTGF = 0.1;//1.0;//0.0001;
  float factorIL6 = 0.01;//0.0001
#ifndef CALIBRATION
  float TGFinc = 1.0 + patchTNF + patchIL10;
  float FGFinc = 1.0;

  TGFinc *= factorTGF;
  FGFinc *= factorFGF;
#else  // CALIBRATION
  float TGFinc = Macrophage::cytokineSynthesis[15]*
                (Macrophage::cytokineSynthesis[16] + patchTNF + patchIL10);
  float FGFinc = Macrophage::cytokineSynthesis[17];
#endif  // CALIBRATION

#ifdef OPT_CHEM
  this->agentWorldPtr->WHWorldChem->incDchem(TGF, in, TGFinc);
  this->agentWorldPtr->WHWorldChem->incDchem(FGF, in, FGFinc);
#else
  (this->agentWorldPtr->WHWorldChem->dTGF[in]) += TGFinc;
  (this->agentWorldPtr->WHWorldChem->dFGF[in]) += FGFinc;
#endif

  /* Activated macrophages synthesize new cytokines in quantities dependent on
   * the vocal treatment type. TODO(Kim): INSERT REFS? */
  float TNFinc  = 0;
  float IL1inc  = 0;
  float IL6inc  = 0;
  float IL8inc  = 0;
  float IL10inc = 0;
#ifndef CALIBRATION
  if (this->agentWorldPtr->treatmentOption == voicerest) {
    TNFinc  = bTNF*0.05  + (1 + patchTNF + patchIL1beta + factorTNF_fHA*countfHA)/
		(1 + patchTGF + patchIL10);
    IL1inc  = bIL1beta*2 + (15 + patchTNF + patchIL1beta*15 + factorIL1_fHA*countfHA)/(1 + patchTGF + patchIL10);
    IL6inc  = bIL6*0.01  + ((1 + patchTNF +patchIL1beta)/(1 + patchIL10));
    IL8inc  = bIL8*2     + ((10 + patchTNF*5 + patchIL1beta*5 + factorIL8_fHA*countfHA)/(1 + patchIL10));
    IL10inc = bIL10*0.01 + (1 + patchIL10*0.1);
  } else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
    TNFinc  = bTNF*0.5 + 2*(1 + patchTNF + patchIL1beta + countfHA + rvis*0.1)/
             (1 + patchTGF + patchIL6*0.1 + patchIL10);
    IL1inc  = bIL1beta*0.5 + (1 + patchTNF*0.5 + patchIL1beta + countfHA + rvis)/
             (1 + patchTGF + patchIL6 + patchIL10);
    IL8inc  = bIL8 + (100 + patchTNF + patchIL1beta*100 + countfHA + rvis*0.1)/(1 + patchIL10*0.5);
    IL10inc  = bIL10*0.01 + (4 + patchIL6*0.001 + patchIL10*0.001);
    if (((Agent::agentWorldPtr)->reportDay()) <= 7)
    {
      IL6inc = 0.5*(1 + rvvs*10);
    } else {
      IL6inc = bIL6 + 0.5*(1 + patchTNF + patchIL1beta)/(1 + patchIL10);
    } 
  } else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
    TNFinc  = bTNF + (1 + patchTNF + patchIL1beta*5 + countfHA + ssis*0.1)/
             (1 + patchTGF + patchIL6 + patchIL10);
    IL1inc  = bIL1beta*13 + 1*(1 + patchTNF + patchIL1beta + countfHA + ssis)/
             (1 + patchTGF + patchIL6 + patchIL10);
    IL8inc  = bIL8 + 10*(1 + patchTNF + patchIL1beta*7 + countfHA + ssis*0.1)/(1 + patchIL10*0.5);
    IL10inc = bIL10 *0.05 + (1 + patchIL6*0.0005 + patchIL10*0.0005);
    if ((Agent::agentWorldPtr)->reportDay() <= 7)
    {
      IL6inc = (1 + ssvs*10);
    } else {
      IL6inc = bIL6 + (1 + patchTNF + patchIL1beta*4)/(1 + patchIL10*0.5);
    }
  }

  TNFinc *= factorTNF;
  IL1inc *= factorIL1;
  IL6inc *= factor;
  IL8inc *= factor;
  IL10inc *= factorIL10;

#else		// CALIBRATION
  if (this->agentWorldPtr->treatmentOption == voicerest) {
    TNFinc  = bTNF*Macrophage::cytokineSynthesis[0] + Macrophage::cytokineSynthesis[1]*
              (Macrophage::cytokineSynthesis[2] + patchTNF + patchIL1beta + countfHA)/
              (Macrophage::cytokineSynthesis[3] + patchTGF + patchIL10);
    IL1inc  = bIL1beta*Macrophage::cytokineSynthesis[18] + Macrophage::cytokineSynthesis[19]*
              (Macrophage::cytokineSynthesis[20] + patchTNF +
               patchIL1beta*Macrophage::cytokineSynthesis[21] + countfHA)/
              (Macrophage::cytokineSynthesis[22] + patchTGF + patchIL10);
    IL6inc  = bIL6*Macrophage::cytokineSynthesis[32] + Macrophage::cytokineSynthesis[33]*
              (Macrophage::cytokineSynthesis[34] + patchTNF + patchIL1beta)/
              (Macrophage::cytokineSynthesis[35] + patchIL10);
    IL8inc  = bIL8*Macrophage::cytokineSynthesis[50] + Macrophage::cytokineSynthesis[51]*
              (Macrophage::cytokineSynthesis[52] + patchTNF*Macrophage::cytokineSynthesis[53] +
               patchIL1beta*Macrophage::cytokineSynthesis[54] + countfHA)/
              (Macrophage::cytokineSynthesis[55] + patchIL10);
    IL10inc = bIL10*Macrophage::cytokineSynthesis[68] + Macrophage::cytokineSynthesis[69]*
              (Macrophage::cytokineSynthesis[70] + patchIL10*Macrophage::cytokineSynthesis[71]);
  } else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
    TNFinc  = bTNF*Macrophage::cytokineSynthesis[4] + Macrophage::cytokineSynthesis[5]*
              (Macrophage::cytokineSynthesis[6] + patchTNF + patchIL1beta + countfHA +
               rvis*Macrophage::cytokineSynthesis[7])/
              (Macrophage::cytokineSynthesis[8] + patchTGF +
               patchIL6*Macrophage::cytokineSynthesis[9] + patchIL10); 
    IL1inc  = bIL1beta*Macrophage::cytokineSynthesis[23] + Macrophage::cytokineSynthesis[24]*
             (Macrophage::cytokineSynthesis[25] +
              patchTNF*Macrophage::cytokineSynthesis[26] + patchIL1beta + countfHA + rvis)/
             (Macrophage::cytokineSynthesis[27] + patchTGF + patchIL6 + patchIL10); 
    IL8inc  = bIL8 + Macrophage::cytokineSynthesis[56]*
              (Macrophage::cytokineSynthesis[57] + patchTNF +
               patchIL1beta*Macrophage::cytokineSynthesis[58] + countfHA +
               rvis*Macrophage::cytokineSynthesis[59])/
              (Macrophage::cytokineSynthesis[60] + patchIL10*Macrophage::cytokineSynthesis[61]);
    IL10inc = bIL10*Macrophage::cytokineSynthesis[72] + Macrophage::cytokineSynthesis[73]*
              (Macrophage::cytokineSynthesis[74] + patchIL6*Macrophage::cytokineSynthesis[75] +
               patchIL10*Macrophage::cytokineSynthesis[76]);
    if (((Agent::agentWorldPtr)->reportDay()) <= 7)
    {
      IL6inc = Macrophage::cytokineSynthesis[36]*
               (Macrophage::cytokineSynthesis[37] + rvvs*Macrophage::cytokineSynthesis[38]);
    } else {
      IL6inc = bIL6 + Macrophage::cytokineSynthesis[39]*(Macrophage::cytokineSynthesis[40] +
               patchTNF + patchIL1beta)/(Macrophage::cytokineSynthesis[41] + patchIL10);
    } 
  } else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
    TNFinc  = bTNF + Macrophage::cytokineSynthesis[10]*
              (Macrophage::cytokineSynthesis[11] + patchTNF +
               patchIL1beta*Macrophage::cytokineSynthesis[12] + countfHA +
               ssis*Macrophage::cytokineSynthesis[13])/
              (Macrophage::cytokineSynthesis[14] + patchTGF + patchIL6 + patchIL10); 
    IL1inc  = bIL1beta*Macrophage::cytokineSynthesis[28] + Macrophage::cytokineSynthesis[29]*
              (Macrophage::cytokineSynthesis[30] + patchTNF + patchIL1beta + countfHA + ssis)/
              (Macrophage::cytokineSynthesis[31] + patchTGF + patchIL6 + patchIL10); 
    IL8inc  = bIL8 + Macrophage::cytokineSynthesis[62]*
              (Macrophage::cytokineSynthesis[63] + patchTNF + 
               patchIL1beta*Macrophage::cytokineSynthesis[64] + countfHA + 
               ssis*Macrophage::cytokineSynthesis[65])/
              (Macrophage::cytokineSynthesis[66] + patchIL10*Macrophage::cytokineSynthesis[67]); 
    IL10inc = bIL10*Macrophage::cytokineSynthesis[77] + Macrophage::cytokineSynthesis[78]*
              (Macrophage::cytokineSynthesis[79] + patchIL6*
               Macrophage::cytokineSynthesis[80] + patchIL10*Macrophage::cytokineSynthesis[81]); 
    if ((Agent::agentWorldPtr)->reportDay() <= 7)
    {
      IL6inc = Macrophage::cytokineSynthesis[42]*
               (Macrophage::cytokineSynthesis[43] + ssvs*Macrophage::cytokineSynthesis[44]); 
    } else {
      IL6inc = bIL6 + Macrophage::cytokineSynthesis[45]*
               (Macrophage::cytokineSynthesis[46] + patchTNF + 
                patchIL1beta*Macrophage::cytokineSynthesis[47])/
               (Macrophage::cytokineSynthesis[48] + patchIL10*Macrophage::cytokineSynthesis[49]);
    }
  }

#endif		// CALIBRATION

#ifdef OPT_CHEM
  this->agentWorldPtr->WHWorldChem->incDchem(TNF,     in, TNFinc);
  this->agentWorldPtr->WHWorldChem->incDchem(IL1beta, in, IL1inc);
  this->agentWorldPtr->WHWorldChem->incDchem(IL6,     in, IL6inc);
  this->agentWorldPtr->WHWorldChem->incDchem(IL8,     in, IL8inc);
  this->agentWorldPtr->WHWorldChem->incDchem(IL10,    in, IL10inc);
#else		// OPT_CHEM
  (this->agentWorldPtr->WHWorldChem->dTNF[in])     += TNFinc; 
  (this->agentWorldPtr->WHWorldChem->dIL1beta[in]) += IL1inc;
  (this->agentWorldPtr->WHWorldChem->dIL6[in])     += IL6inc;
  (this->agentWorldPtr->WHWorldChem->dIL8[in])     += IL8inc;
  (this->agentWorldPtr->WHWorldChem->dIL10[in])    += IL10inc;
#endif		// OPT_CHEM


#ifdef PRINT_SECRETION
  int x = this->ix[read_t];
  int y = this->iy[read_t];
  int z = this->iz[read_t];
  printCytRelease(2, TGF,     x, y, z, TGFinc);
  printCytRelease(2, FGF,     x, y, z, FGFinc);
  printCytRelease(2, TNF,     x, y, z, TNFinc);
  printCytRelease(2, IL1beta, x, y, z, IL1inc);
  printCytRelease(2, IL6,     x, y, z, IL6inc);
  printCytRelease(2, IL8,     x, y, z, IL8inc);
  printCytRelease(2, IL10,    x, y, z, IL10inc);
#endif  // PRINT_SECRETION


  /*************************************************************************
   * DEACTIVATION                                                          *
   *************************************************************************/
  /* Activated macrophages might be deactivated once the damage is cleared. TODO(Kim): INSERT REF? */
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[pdamage];
	int initialdam = (Agent::agentWorldPtr)->getInitialDam();
	int threshold_dam = (initialdam*0)/10;	// 0% if the initial damage
#ifndef CALIBRATION
	if (totaldamage <= threshold_dam && rollDice(50)) this->macDeactivation();
#else  // CALIBRATION
	if (totaldamage <= threshold_dam && rollDice(Macrophage::activation[4])) this->macDeactivation();
#endif  // CALIBRATION

  /*************************************************************************
   * SIGNAL TRANSDUCTION                                                   *
   *************************************************************************/
	Agent::agentWorldPtr->highTNFdamage = true;

  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  // Activated macrophages can die naturally
	this->life[write_t] = this->life[read_t] - 1;
	if (this->life[read_t] <= 0) { 
    this->die();
  }
}

