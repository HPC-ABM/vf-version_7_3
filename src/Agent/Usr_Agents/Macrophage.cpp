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

float Macrophage::cytokineSynthesis[82] = {0.0};
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

void Macrophage::amac_produce_cytokines()
{
	// Calculates neighboring ECM and patch chemical concentrations
	const int currentDay = (Agent::agentWorldPtr)->reportDay();
	int treatment = this->agentWorldPtr->treatmentOption;
	int in = this->index[read_t];
	float rvis = Agent::agentWorldPtr->RVIS;
	float rvvs = Agent::agentWorldPtr->RVVS;
  float ssis = Agent::agentWorldPtr->SSIS;
	float ssvs = Agent::agentWorldPtr->SSVS;

	float TGFinc  = 0.0f;
	float FGFinc  = 0.0f;
	float TNFinc  = 0.0f;
	float IL1inc  = 0.0f;
	float IL6inc  = 0.0f;
	float IL8inc  = 0.0f;
	float IL10inc = 0.0f;

	float countfHA = this->countNeighborECM(fha);
//	float bTNF     = (Agent::agentWorldPtr)->baselineChem[TNF];
//	float bIL1beta = (Agent::agentWorldPtr)->baselineChem[IL1beta];
//	float bIL6     = (Agent::agentWorldPtr)->baselineChem[IL6];
//	float bIL8     = (Agent::agentWorldPtr)->baselineChem[IL8];
//	float bIL10    = (Agent::agentWorldPtr)->baselineChem[IL10];

	float patchTNF  = this->agentWorldPtr->WHWorldChem->pTNF[in];
	float patchTGF  = this->agentWorldPtr->WHWorldChem->pTGF[in];
	float patchIL1  = this->agentWorldPtr->WHWorldChem->pIL1beta[in];
	float patchIL6  = this->agentWorldPtr->WHWorldChem->pIL6[in];
	float patchIL8  = this->agentWorldPtr->WHWorldChem->pIL8[in];
	float patchIL10 = this->agentWorldPtr->WHWorldChem->pIL10[in];

	// TGF
	float tgf_cNum  = 0.0f;
	float tgf_cDen  = 0.0f;
	float tgf_cTNF  = 0.0f;
	float tgf_cIL10 = 0.0f;

	// FGF
	float fgf_cNum = 0.0f;
	float fgf_cDen = 0.0f;

	// TNF
	float tnf_cNum  = 0.0f;	 float tnf_cDen  = 0.0f;
	float tnf_cTNF  = 0.0f;  float tnf_cTGF  = 0.0f;
	float tnf_cIL1  = 0.0f;  float tnf_cIL10 = 0.0f;
	float tnf_cFHA  = 0.0f;  float tnf_cIL6  = 0.0f;
	float tnf_cRVIS = 0.0f;

	// IL1
	float il1_cNum  = 0.0f;   float il1_cDen  = 0.0f;
	float il1_cTNF  = 0.0f;   float il1_cTGF  = 0.0f;
	float il1_cIL1  = 0.0f;   float il1_cIL10 = 0.0f;
	float il1_cFHA  = 0.0f;   float il1_cIL6  = 0.0f;
	float il1_cRVIS = 0.0f;
	float il1_cSSIS = 0.0f;

	// IL6
	float il6_cNum = 0.0f;   float il6_cDen  = 0.0f;
	float il6_cTNF = 0.0f;   float il6_cIL10 = 0.0f;
	float il6_cIL1 = 0.0f;

	REGULATORS_T il6_ucoefs;
	REGULATORS_T il6_dcoefs;
	REGULATORS_T il6_uregs;
	REGULATORS_T il6_dregs;

	// IL8
	float il8_cNum  = 0.0f;   float il8_cDen  = 0.0f;
	float il8_cTNF  = 0.0f;   float il8_cIL10 = 0.0f;
	float il8_cIL1  = 0.0f;
	float il8_cFHA  = 0.0f;
	float il8_cRVIS = 0.0f;
	float il8_cSSIS = 0.0f;

	// IL10
  float il10_cNum  = 0.0f; float il10_cDen  = 0.0f;
  float il10_cIL10 = 0.0f;
  float il10_cIL6  = 0.0f;
  float il6_cRVVS  = 0.0f;
  float il6_cSSVS  = 0.0f;

	float tgf_ofst  = Agent::baseline[TGF];
	float fgf_ofst  = Agent::baseline[FGF];
	float tnf_ofst  = Agent::baseline[TNF];
	float il1_ofst  = Agent::baseline[IL1beta];
	float il6_ofst  = Agent::baseline[IL6];
	float il8_ofst  = Agent::baseline[IL8];
	float il10_ofst = Agent::baseline[IL10];


#ifndef CALIBRATION

  float factorIL8  = 0.0005;//0.1;//0.0001;
  float factorIL1  = 0.0001;//0.01;//0.097;
  float factorTNF  = 0.000001;//0.85;//0.4;//0.0001;
  float factorIL10 = 0.1;//0.1;
  float factorFGF  = 0.0001;
  float factorTGF  = 0.1;//1.0;//0.0001;
  float factorIL6  = 0.0005;//0.01;//0.0001
  float factorTNF_fHA = 0.1;//0.1;//100.0;
  float factorIL1_fHA = 1.0;
  float factorIL8_fHA = 1.0;

  tgf_cNum  = 1.0f*factorTGF;
  tgf_cTNF  = 1.0f*factorTGF;
  tgf_cIL10 = 1.0f*factorTGF;
  tgf_cDen  = 1.0f;

  fgf_cNum = 1.0f*factorFGF;
  fgf_cDen = 1.0f;

	REGULATORS_T tgf_ucoefs{tgf_cNum, tgf_cTNF, tgf_cIL10};
	REGULATORS_T tgf_dcoefs{tgf_cDen};
	REGULATORS_T tgf_uregs{1.0f, patchTNF, patchIL10};
	REGULATORS_T tgf_dregs{1.0f};

	REGULATORS_T fgf_ucoefs{fgf_cNum};
	REGULATORS_T fgf_dcoefs{fgf_cDen};
	REGULATORS_T fgf_uregs{1.0f};
	REGULATORS_T fgf_dregs{1.0f};

	switch (treatment)
	{
	case voicerest:
	{
		// TNF up-regulator           | down-regulator coeffs
		tnf_cNum = 1.0f*factorTNF;			tnf_cDen  = 1.0f;
		tnf_cTNF = 1.0f*factorTNF; 			tnf_cTGF  = 1.0f;
		tnf_cIL1 = 1.0f*factorTNF; 			tnf_cIL10 = 1.0f;
		tnf_cFHA = factorTNF_fHA*factorTNF;

		REGULATORS_T tnf_ucoefs{tnf_cNum, tnf_cTNF, tnf_cIL1, tnf_cFHA};
		REGULATORS_T tnf_dcoefs{tnf_cDen, tnf_cTGF, tnf_cIL10};
		REGULATORS_T tnf_uregs{1.0f, patchTNF, patchIL1, countfHA};
		REGULATORS_T tnf_dregs{1.0f, patchTGF, patchIL10};

		// IL1 up-regulator           | down-regulator coeffs
		il1_cNum = 15.0f*factorIL1;     il1_cDen  = 1.0f;
		il1_cTNF =  1.0f*factorIL1;     il1_cTGF  = 1.0f;
		il1_cIL1 = 15.0f*factorIL1;     il1_cIL10 = 1.0f;
		il1_cFHA = factorIL1_fHA*factorIL1;

		REGULATORS_T il1_ucoefs{il1_cNum, il1_cTNF, il1_cIL1, il1_cFHA};
		REGULATORS_T il1_dcoefs{il1_cDen, il1_cTGF, il1_cIL10};
		REGULATORS_T il1_uregs{1.0f, patchTNF, patchIL1, countfHA};
		REGULATORS_T il1_dregs{1.0f, patchTGF, patchIL10};

		// IL6 up-regulator | down-regulator coeffs
		il6_cNum = 1.0f*factorIL6;      il6_cDen  = 1.0f;
		il6_cTNF = 1.0f*factorIL6;      il6_cIL10 = 1.0f;
		il6_cIL1 = 1.0f*factorIL6;

		il6_ucoefs = {il6_cNum, il6_cTNF, il6_cIL1};
		il6_dcoefs = {il6_cDen, il6_cIL10};
		il6_uregs  = {1.0f, patchTNF, patchIL1};
		il6_dregs  = {1.0f, patchIL10};

		// IL8 up-regulator | down-regulator coeffs
		il8_cNum = 10.0f*factorIL8;     il8_cDen  = 1.0f;
		il8_cTNF =  5.0f*factorIL8;     il8_cIL10 = 1.0f;
		il8_cIL1 =  5.0f*factorIL8;
		il8_cFHA = factorIL8_fHA*factorIL8;

		REGULATORS_T il8_ucoefs{il8_cNum, il8_cTNF, il8_cIL1, il8_cFHA};
		REGULATORS_T il8_dcoefs{il8_cDen, il8_cIL10};
		REGULATORS_T il8_uregs{1.0f, patchTNF, patchIL1, countfHA};
		REGULATORS_T il8_dregs{1.0f, patchIL10};

		// IL10 up-regulator        | down-regulator coeffs
    il10_cNum  = 1.0f*factorIL10; il10_cDen  = 1.0f;
    il10_cIL10 = 0.1f*factorIL10;

		REGULATORS_T il10_ucoefs{il10_cNum, il10_cIL10};
		REGULATORS_T il10_dcoefs{il10_cDen};
		REGULATORS_T il10_uregs{1.0f, patchIL10};
		REGULATORS_T il10_dregs{1.0f};

		// Call cytokine production functions
		TGFinc  = produce_tgf ( tgf_uregs,  tgf_dregs,  tgf_ucoefs,  tgf_dcoefs,  tgf_ofst);
		FGFinc  = produce_fgf ( fgf_uregs,  fgf_dregs,  fgf_ucoefs,  fgf_dcoefs,  fgf_ofst);
		TNFinc  = produce_tnf ( tnf_uregs,  tnf_dregs,  tnf_ucoefs,  tnf_dcoefs,  tnf_ofst);
		IL1inc  = produce_il1 ( il1_uregs,  il1_dregs,  il1_ucoefs,  il1_dcoefs,  il1_ofst);
	  IL6inc  = produce_il6 ( il6_uregs,  il6_dregs,  il6_ucoefs,  il6_dcoefs,  il6_ofst);
	  IL8inc  = produce_il8 ( il8_uregs,  il8_dregs,  il8_ucoefs,  il8_dcoefs,  il8_ofst);
	  IL10inc = produce_il10(il10_uregs, il10_dregs, il10_ucoefs, il10_dcoefs, il10_ofst);

		break;
	}
	case resonantvoice:
	{
		// TNF up-regulator | down-regulator coeffs
		tnf_cNum  = 2.0f;     tnf_cDen  = 1.0f;
		tnf_cTNF  = 2.0f;     tnf_cTGF  = 1.0f;
		tnf_cIL1  = 2.0f;     tnf_cIL6  = 0.1f;
		tnf_cFHA  = 2.0f;     tnf_cIL10 = 1.0f;
		tnf_cRVIS = 0.2f;

		REGULATORS_T tnf_ucoefs{tnf_cNum, tnf_cTNF, tnf_cIL1, tnf_cFHA, tnf_cRVIS};
		REGULATORS_T tnf_dcoefs{tnf_cDen, tnf_cTGF, tnf_cIL6, tnf_cIL10};
		REGULATORS_T tnf_uregs{1.0f, patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T tnf_dregs{1.0f, patchTGF, patchIL6, patchIL10};

		// IL1 up-regulator | down-regulator coeffs
		il1_cNum  = 1.0f;     il1_cDen  = 1.0f;
		il1_cTNF  = 0.5f;     il1_cTGF  = 1.0f;
		il1_cIL1  = 1.0f;     il1_cIL6  = 1.0f;
		il1_cFHA  = 1.0f;     il1_cIL10 = 1.0f;
		il1_cRVIS = 1.0f;

		REGULATORS_T il1_ucoefs{il1_cNum, il1_cTNF, il1_cIL1, il1_cFHA, il1_cRVIS};
		REGULATORS_T il1_dcoefs{il1_cDen, il1_cTGF, il1_cIL10};
		REGULATORS_T il1_uregs{1.0f, patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T il1_dregs{1.0f, patchTGF, patchIL6, patchIL10};

		// IL8 up-regulator | down-regulator coeffs
		il8_cNum  = 100.0f;     il8_cDen  = 1.0f;
		il8_cTNF  =   1.0f;     il8_cIL10 = 0.5f;
		il8_cIL1  = 100.0f;
		il8_cFHA  =   1.0f;
		il8_cRVIS =   0.1f;

		REGULATORS_T il8_ucoefs{il8_cNum, il8_cTNF, il8_cIL1, il8_cFHA, il8_cRVIS};
		REGULATORS_T il8_dcoefs{il8_cDen, il8_cIL10};
		REGULATORS_T il8_uregs{1.0f, patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T il8_dregs{1.0f, patchIL10};

		// IL10 up-regulator | down-regulator coeffs
    il10_cNum  = 4.000f; il10_cDen  = 1.0f;
    il10_cIL6  = 0.001f;
    il10_cIL10 = 0.001f;

		REGULATORS_T il10_ucoefs{il10_cNum, il10_cIL6, il10_cIL10};
		REGULATORS_T il10_dcoefs{il10_cDen};
		REGULATORS_T il10_uregs{1.0f, patchIL6, patchIL10};
		REGULATORS_T il10_dregs{1.0f};

		// IL6 up-regulator | down-regulator coeffs
    if (currentDay <= 7)
    {
  		il6_cNum  = 0.5f;      il6_cDen  = 1.0f;
  		il6_cRVVS = 5.0f;

  		il6_ucoefs = {il6_cNum, il6_cRVVS};
  		il6_dcoefs = {il6_cDen};
  		il6_uregs  = {1.0f, rvvs};
  		il6_dregs  = {1.0f};

  		il6_ofst = 0.0f;	// TODO: Verify --- Old code has no baseline
    } else {
  		il6_cNum = 0.5f;      il6_cDen  = 1.0f;
  		il6_cTNF = 0.5f;      il6_cIL10 = 1.0f;
  		il6_cIL1 = 0.5f;

  		il6_ucoefs = {il6_cNum, il6_cTNF, il6_cIL1};
  		il6_dcoefs = {il6_cDen, il6_cIL10};
  		il6_uregs  = {1.0f, patchTNF, patchIL1};
  		il6_dregs  = {1.0f, patchIL10};
    }

  	// Call cytokine production functions
  	TGFinc  = produce_tgf ( tgf_uregs,  tgf_dregs,  tgf_ucoefs,  tgf_dcoefs,  tgf_ofst);
  	FGFinc  = produce_fgf ( fgf_uregs,  fgf_dregs,  fgf_ucoefs,  fgf_dcoefs,  fgf_ofst);
  	TNFinc  = produce_tnf ( tnf_uregs,  tnf_dregs,  tnf_ucoefs,  tnf_dcoefs,  tnf_ofst);
  	IL1inc  = produce_il1 ( il1_uregs,  il1_dregs,  il1_ucoefs,  il1_dcoefs,  il1_ofst);
    IL6inc  = produce_il6 ( il6_uregs,  il6_dregs,  il6_ucoefs,  il6_dcoefs,  il6_ofst);
    IL8inc  = produce_il8 ( il8_uregs,  il8_dregs,  il8_ucoefs,  il8_dcoefs,  il8_ofst);
    IL10inc = produce_il10(il10_uregs, il10_dregs, il10_ucoefs, il10_dcoefs, il10_ofst);

		break;
	}

	case spontaneousspeech:
	{

		// TNF up-regulator | down-regulator coeffs
		tnf_cNum  = 1.0f;     tnf_cDen  = 1.0f;
		tnf_cTNF  = 1.0f;     tnf_cTGF  = 1.0f;
		tnf_cIL1  = 5.0f;     tnf_cIL6  = 1.0f;
		tnf_cFHA  = 1.0f;     tnf_cIL10 = 1.0f;
		tnf_cRVIS = 0.1f;

		REGULATORS_T tnf_ucoefs{tnf_cNum, tnf_cTNF, tnf_cIL1, tnf_cFHA};
		REGULATORS_T tnf_dcoefs{tnf_cDen, tnf_cTGF, tnf_cIL6, tnf_cIL10};
		REGULATORS_T tnf_uregs{1.0f, patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T tnf_dregs{1.0f, patchTGF, patchIL6, patchIL10};

		// IL1 up-regulator | down-regulator coeffs
		il1_cNum  = 1.0f;     il1_cDen  = 1.0f;
		il1_cTNF  = 1.0f;     il1_cTGF  = 1.0f;
		il1_cIL1  = 1.0f;     il1_cIL6  = 1.0f;
		il1_cFHA  = 1.0f;     il1_cIL10 = 1.0f;
		il1_cSSIS = 1.0f;

		REGULATORS_T il1_ucoefs{il1_cNum, il1_cTNF, il1_cIL1, il1_cFHA, il1_cSSIS};
		REGULATORS_T il1_dcoefs{il1_cDen, il1_cTGF, il1_cIL10};
		REGULATORS_T il1_uregs{1.0f, patchTNF, patchIL1, countfHA, ssis};
		REGULATORS_T il1_dregs{1.0f, patchTGF, patchIL6, patchIL10};

		// IL8 up-regulator | down-regulator coeffs
		il8_cNum  = 10.0f;     il8_cDen  = 1.0f;
		il8_cTNF  = 10.0f;     il8_cIL10 = 0.5f;
		il8_cIL1  = 70.0f;
		il8_cFHA  = 10.0f;
		il8_cSSIS =  0.1f;

		REGULATORS_T il8_ucoefs{il8_cNum, il8_cTNF, il8_cIL1, il8_cFHA, il8_cSSIS};
		REGULATORS_T il8_dcoefs{il8_cDen, il8_cIL10};
		REGULATORS_T il8_uregs{1.0f, patchTNF, patchIL1, countfHA, ssis};
		REGULATORS_T il8_dregs{1.0f, patchIL10};

		// IL10 up-regulator | down-regulator coeffs
    il10_cNum  = 1.0000f; il10_cDen  = 1.0f;
    il10_cIL6  = 0.0005f;
    il10_cIL10 = 0.0005f;

		REGULATORS_T il10_ucoefs{il10_cNum, il10_cIL6, il10_cIL10};
		REGULATORS_T il10_dcoefs{il10_cDen};
		REGULATORS_T il10_uregs{1.0f, patchIL6, patchIL10};
		REGULATORS_T il10_dregs{1.0f};

		// IL6 up-regulator | down-regulator coeffs
    if (currentDay <= 7)
    {
  		il6_cNum  =  1.0f;  il6_cDen  = 1.0f;
  		il6_cSSVS = 10.0f;

  		il6_ucoefs = {il6_cNum, il6_cSSVS};
  		il6_dcoefs = {il6_cDen};
  		il6_uregs = {1.0f, ssvs};
  		il6_dregs = {1.0f};

  		il6_ofst = 0.0f;	// TODO: Verify --- Old code has no baseline
    } else {
  		il6_cNum = 1.0f;    il6_cDen  = 1.0f;
  		il6_cTNF = 1.0f;    il6_cIL10 = 0.5f;
  		il6_cIL1 = 4.0f;

  		il6_ucoefs = {il6_cNum, il6_cTNF, il6_cIL1};
  		il6_dcoefs = {il6_cDen, il6_cIL10};
  		il6_uregs = {1.0f, patchTNF, patchIL1};
  		il6_dregs = {1.0f, patchIL10};
    }

  	// Call cytokine production functions
  	TGFinc  = produce_tgf ( tgf_uregs,  tgf_dregs,  tgf_ucoefs,  tgf_dcoefs,  tgf_ofst);
  	FGFinc  = produce_fgf ( fgf_uregs,  fgf_dregs,  fgf_ucoefs,  fgf_dcoefs,  fgf_ofst);
  	TNFinc  = produce_tnf ( tnf_uregs,  tnf_dregs,  tnf_ucoefs,  tnf_dcoefs,  tnf_ofst);
  	IL1inc  = produce_il1 ( il1_uregs,  il1_dregs,  il1_ucoefs,  il1_dcoefs,  il1_ofst);
    IL6inc  = produce_il6 ( il6_uregs,  il6_dregs,  il6_ucoefs,  il6_dcoefs,  il6_ofst);
    IL8inc  = produce_il8 ( il8_uregs,  il8_dregs,  il8_ucoefs,  il8_dcoefs,  il8_ofst);
    IL10inc = produce_il10(il10_uregs, il10_dregs, il10_ucoefs, il10_dcoefs, il10_ofst);
	}

	default:
		break;
	}


#else  // CALIBRATION
  tgf_cTNF  = Macrophage::cytokineSynthesis[0];
  tgf_cIL10 = Macrophage::cytokineSynthesis[1];
  tgf_cDen  = 1.0f;

  fgf_cNum = Macrophage::cytokineSynthesis[2];
  fgf_cDen = 1.0f;

	REGULATORS_T tgf_ucoefs{tgf_cTNF, tgf_cIL10};
	REGULATORS_T tgf_dcoefs{tgf_cDen};
	REGULATORS_T tgf_uregs{patchTNF, patchIL10};
	REGULATORS_T tgf_dregs{1.0f};

	REGULATORS_T fgf_ucoefs{fgf_cNum};
	REGULATORS_T fgf_dcoefs{fgf_cDen};
	REGULATORS_T fgf_uregs{1.0f};
	REGULATORS_T fgf_dregs{1.0f};

	switch (this->agentWorldPtr->treatmentOption)
	{
	case voicerest:
	{
		// TNF up-regulator                             | down-regulator coeffs
		tnf_cTNF = Macrophage::cytokineSynthesis[3]; 			tnf_cTGF  = Macrophage::cytokineSynthesis[6];
		tnf_cIL1 = Macrophage::cytokineSynthesis[4]; 			tnf_cIL10 = Macrophage::cytokineSynthesis[7];
		tnf_cFHA = Macrophage::cytokineSynthesis[5];

		REGULATORS_T tnf_ucoefs{tnf_cTNF, tnf_cIL1, tnf_cFHA};
		REGULATORS_T tnf_dcoefs{tnf_cTGF, tnf_cIL10};
		REGULATORS_T tnf_uregs{patchTNF, patchIL1, countfHA};
		REGULATORS_T tnf_dregs{patchTGF, patchIL10};

		// IL1 up-regulator                             | down-regulator coeffs
		il1_cTNF = Macrophage::cytokineSynthesis[8];     il1_cTGF  = Macrophage::cytokineSynthesis[11];
		il1_cIL1 = Macrophage::cytokineSynthesis[9];     il1_cIL10 = Macrophage::cytokineSynthesis[12];
		il1_cFHA = Macrophage::cytokineSynthesis[10];

		REGULATORS_T il1_ucoefs{il1_cTNF, il1_cIL1, il1_cFHA};
		REGULATORS_T il1_dcoefs{il1_cTGF, il1_cIL10};
		REGULATORS_T il1_uregs{patchTNF, patchIL1, countfHA};
		REGULATORS_T il1_dregs{patchTGF, patchIL10};

		// IL6 up-regulator                             | down-regulator coeffs
		il6_cTNF = Macrophage::cytokineSynthesis[13];      il6_cIL10 = Macrophage::cytokineSynthesis[15];
		il6_cIL1 = Macrophage::cytokineSynthesis[14];

		il6_ucoefs = {il6_cTNF, il6_cIL1};
		il6_dcoefs = {il6_cIL10};
		il6_uregs  = {patchTNF, patchIL1};
		il6_dregs  = {patchIL10};

		// IL8 up-regulator                             | down-regulator coeffs
		il8_cTNF = Macrophage::cytokineSynthesis[16];     il8_cIL10 = Macrophage::cytokineSynthesis[19];
		il8_cIL1 = Macrophage::cytokineSynthesis[17];
		il8_cFHA = Macrophage::cytokineSynthesis[18];

		REGULATORS_T il8_ucoefs{il8_cTNF, il8_cIL1, il8_cFHA};
		REGULATORS_T il8_dcoefs{il8_cIL10};
		REGULATORS_T il8_uregs{patchTNF, patchIL1, countfHA};
		REGULATORS_T il8_dregs{patchIL10};

		// IL10 up-regulator                             | down-regulator coeffs
    il10_cIL10 = Macrophage::cytokineSynthesis[20];

		REGULATORS_T il10_ucoefs{il10_cIL10};
		REGULATORS_T il10_dcoefs{1.0f};
		REGULATORS_T il10_uregs{patchIL10};
		REGULATORS_T il10_dregs{1.0f};

		// Call cytokine production functions
		TGFinc  = produce_tgf ( tgf_uregs,  tgf_dregs,  tgf_ucoefs,  tgf_dcoefs,  tgf_ofst);
		FGFinc  = produce_fgf ( fgf_uregs,  fgf_dregs,  fgf_ucoefs,  fgf_dcoefs,  fgf_ofst);
		TNFinc  = produce_tnf ( tnf_uregs,  tnf_dregs,  tnf_ucoefs,  tnf_dcoefs,  tnf_ofst);
		IL1inc  = produce_il1 ( il1_uregs,  il1_dregs,  il1_ucoefs,  il1_dcoefs,  il1_ofst);
	  IL6inc  = produce_il6 ( il6_uregs,  il6_dregs,  il6_ucoefs,  il6_dcoefs,  il6_ofst);
	  IL8inc  = produce_il8 ( il8_uregs,  il8_dregs,  il8_ucoefs,  il8_dcoefs,  il8_ofst);
	  IL10inc = produce_il10(il10_uregs, il10_dregs, il10_ucoefs, il10_dcoefs, il10_ofst);

		break;
	}
	case resonantvoice:
	{
		// TNF up-regulator | down-regulator coeffs
		tnf_cTNF  = Macrophage::cytokineSynthesis[21];     tnf_cTGF  = Macrophage::cytokineSynthesis[28];
		tnf_cIL1  = Macrophage::cytokineSynthesis[22];     tnf_cIL6  = Macrophage::cytokineSynthesis[29];
		tnf_cFHA  = Macrophage::cytokineSynthesis[23];     tnf_cIL10 = Macrophage::cytokineSynthesis[30];
		tnf_cRVIS = Macrophage::cytokineSynthesis[24];

		REGULATORS_T tnf_ucoefs{tnf_cTNF, tnf_cIL1, tnf_cFHA, tnf_cRVIS};
		REGULATORS_T tnf_dcoefs{tnf_cTGF, tnf_cIL6, tnf_cIL10};
		REGULATORS_T tnf_uregs{patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T tnf_dregs{patchTGF, patchIL6, patchIL10};

		// IL1 up-regulator | down-regulator coeffs
		il1_cTNF  = Macrophage::cytokineSynthesis[31];     il1_cTGF  = Macrophage::cytokineSynthesis[35];
		il1_cIL1  = Macrophage::cytokineSynthesis[32];     il1_cIL6  = Macrophage::cytokineSynthesis[36];
		il1_cFHA  = Macrophage::cytokineSynthesis[33];     il1_cIL10 = Macrophage::cytokineSynthesis[37];
		il1_cRVIS = Macrophage::cytokineSynthesis[34];

		REGULATORS_T il1_ucoefs{il1_cTNF, il1_cIL1, il1_cFHA, il1_cRVIS};
		REGULATORS_T il1_dcoefs{il1_cTGF, il1_cIL10};
		REGULATORS_T il1_uregs{patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T il1_dregs{patchTGF, patchIL6, patchIL10};

		// IL8 up-regulator | down-regulator coeffs
		il8_cTNF  = Macrophage::cytokineSynthesis[38];     il8_cIL10 = Macrophage::cytokineSynthesis[42];
		il8_cIL1  = Macrophage::cytokineSynthesis[39];
		il8_cFHA  = Macrophage::cytokineSynthesis[40];
		il8_cRVIS = Macrophage::cytokineSynthesis[41];

		REGULATORS_T il8_ucoefs{il8_cTNF, il8_cIL1, il8_cFHA, il8_cRVIS};
		REGULATORS_T il8_dcoefs{il8_cIL10};
		REGULATORS_T il8_uregs{patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T il8_dregs{patchIL10};

		// IL10 up-regulator | down-regulator coeffs
    il10_cIL6  = Macrophage::cytokineSynthesis[43];
    il10_cIL10 = Macrophage::cytokineSynthesis[44];

		REGULATORS_T il10_ucoefs{il10_cIL6, il10_cIL10};
		REGULATORS_T il10_dcoefs{1.0f};
		REGULATORS_T il10_uregs{patchIL6, patchIL10};
		REGULATORS_T il10_dregs{1.0f};

		// IL6 up-regulator | down-regulator coeffs
    if (currentDay <= 7)
    {
  		il6_cNum  = Macrophage::cytokineSynthesis[45];      il6_cDen  = Macrophage::cytokineSynthesis[47];
  		il6_cRVVS = Macrophage::cytokineSynthesis[46];

  		il6_ucoefs = {il6_cRVVS};
  		il6_dcoefs = {1.0f};
  		il6_uregs  = {rvvs};
  		il6_dregs  = {1.0f};

  		il6_ofst = 0.0f;	// TODO: Verify --- Old code has no baseline
    } else {
  		il6_cTNF = Macrophage::cytokineSynthesis[47];      il6_cIL10 = Macrophage::cytokineSynthesis[49];
  		il6_cIL1 = Macrophage::cytokineSynthesis[48];

  		il6_ucoefs = {il6_cTNF, il6_cIL1};
  		il6_dcoefs = {il6_cIL10};
  		il6_uregs  = {patchTNF, patchIL1};
  		il6_dregs  = {patchIL10};
    }

  	// Call cytokine production functions
  	TGFinc  = produce_tgf ( tgf_uregs,  tgf_dregs,  tgf_ucoefs,  tgf_dcoefs,  tgf_ofst);
  	FGFinc  = produce_fgf ( fgf_uregs,  fgf_dregs,  fgf_ucoefs,  fgf_dcoefs,  fgf_ofst);
  	TNFinc  = produce_tnf ( tnf_uregs,  tnf_dregs,  tnf_ucoefs,  tnf_dcoefs,  tnf_ofst);
  	IL1inc  = produce_il1 ( il1_uregs,  il1_dregs,  il1_ucoefs,  il1_dcoefs,  il1_ofst);
    IL6inc  = produce_il6 ( il6_uregs,  il6_dregs,  il6_ucoefs,  il6_dcoefs,  il6_ofst);
    IL8inc  = produce_il8 ( il8_uregs,  il8_dregs,  il8_ucoefs,  il8_dcoefs,  il8_ofst);
    IL10inc = produce_il10(il10_uregs, il10_dregs, il10_ucoefs, il10_dcoefs, il10_ofst);

		break;
	}

	case spontaneousspeech:
	{
		// TNF up-regulator | down-regulator coeffs
		tnf_cTNF  = Macrophage::cytokineSynthesis[50];     tnf_cTGF  = Macrophage::cytokineSynthesis[54];
		tnf_cIL1  = Macrophage::cytokineSynthesis[51];     tnf_cIL6  = Macrophage::cytokineSynthesis[55];
		tnf_cFHA  = Macrophage::cytokineSynthesis[52];     tnf_cIL10 = Macrophage::cytokineSynthesis[56];
		tnf_cRVIS = Macrophage::cytokineSynthesis[53];

		REGULATORS_T tnf_ucoefs{tnf_cTNF, tnf_cIL1, tnf_cFHA};
		REGULATORS_T tnf_dcoefs{tnf_cTGF, tnf_cIL6, tnf_cIL10};
		REGULATORS_T tnf_uregs{patchTNF, patchIL1, countfHA, rvis};
		REGULATORS_T tnf_dregs{patchTGF, patchIL6, patchIL10};

		// IL1 up-regulator | down-regulator coeffs
		il1_cTNF  = Macrophage::cytokineSynthesis[57];     il1_cTGF  = Macrophage::cytokineSynthesis[61];
		il1_cIL1  = Macrophage::cytokineSynthesis[58];     il1_cIL6  = Macrophage::cytokineSynthesis[62];
		il1_cFHA  = Macrophage::cytokineSynthesis[59];     il1_cIL10 = Macrophage::cytokineSynthesis[63];
		il1_cSSIS = Macrophage::cytokineSynthesis[60];

		REGULATORS_T il1_ucoefs{il1_cTNF, il1_cIL1, il1_cFHA, il1_cSSIS};
		REGULATORS_T il1_dcoefs{il1_cTGF, il1_cIL10};
		REGULATORS_T il1_uregs{patchTNF, patchIL1, countfHA, ssis};
		REGULATORS_T il1_dregs{patchTGF, patchIL6, patchIL10};

		// IL8 up-regulator | down-regulator coeffs
		il8_cTNF  = Macrophage::cytokineSynthesis[64];     il8_cIL10 = Macrophage::cytokineSynthesis[68];
		il8_cIL1  = Macrophage::cytokineSynthesis[65];
		il8_cFHA  = Macrophage::cytokineSynthesis[66];
		il8_cSSIS = Macrophage::cytokineSynthesis[67];

		REGULATORS_T il8_ucoefs{il8_cTNF, il8_cIL1, il8_cFHA, il8_cSSIS};
		REGULATORS_T il8_dcoefs{il8_cIL10};
		REGULATORS_T il8_uregs{patchTNF, patchIL1, countfHA, ssis};
		REGULATORS_T il8_dregs{patchIL10};

		// IL10 up-regulator | down-regulator coeffs
    il10_cNum  = Macrophage::cytokineSynthesis[69]; il10_cDen  = Macrophage::cytokineSynthesis[72];
    il10_cIL6  = Macrophage::cytokineSynthesis[70];
    il10_cIL10 = Macrophage::cytokineSynthesis[71];

		REGULATORS_T il10_ucoefs{il10_cIL6, il10_cIL10};
		REGULATORS_T il10_dcoefs{1.0f};
		REGULATORS_T il10_uregs{patchIL6, patchIL10};
		REGULATORS_T il10_dregs{1.0f};

		// IL6 up-regulator | down-regulator coeffs
    if (currentDay <= 7)
    {
  		il6_cSSVS = Macrophage::cytokineSynthesis[73];

  		il6_ucoefs = {il6_cSSVS};
  		il6_dcoefs = {1.0f};
  		il6_uregs  = {ssvs};
  		il6_dregs  = {1.0f};

  		il6_ofst = 0.0f;	// TODO: Verify --- Old code has no baseline
    } else {
  		il6_cNum = Macrophage::cytokineSynthesis[74];    il6_cDen  = Macrophage::cytokineSynthesis[77];
  		il6_cTNF = Macrophage::cytokineSynthesis[75];    il6_cIL10 = Macrophage::cytokineSynthesis[78];
  		il6_cIL1 = Macrophage::cytokineSynthesis[76];

  		il6_ucoefs = {il6_cTNF, il6_cIL1};
  		il6_dcoefs = {il6_cIL10};
  		il6_uregs  = {patchTNF, patchIL1};
  		il6_dregs  = {patchIL10};
    }

  	// Call cytokine production functions
  	TGFinc  = produce_tgf ( tgf_uregs,  tgf_dregs,  tgf_ucoefs,  tgf_dcoefs,  tgf_ofst);
  	FGFinc  = produce_fgf ( fgf_uregs,  fgf_dregs,  fgf_ucoefs,  fgf_dcoefs,  fgf_ofst);
  	TNFinc  = produce_tnf ( tnf_uregs,  tnf_dregs,  tnf_ucoefs,  tnf_dcoefs,  tnf_ofst);
  	IL1inc  = produce_il1 ( il1_uregs,  il1_dregs,  il1_ucoefs,  il1_dcoefs,  il1_ofst);
    IL6inc  = produce_il6 ( il6_uregs,  il6_dregs,  il6_ucoefs,  il6_dcoefs,  il6_ofst);
    IL8inc  = produce_il8 ( il8_uregs,  il8_dregs,  il8_ucoefs,  il8_dcoefs,  il8_ofst);
    IL10inc = produce_il10(il10_uregs, il10_dregs, il10_ucoefs, il10_dcoefs, il10_ofst);

    break;
	}

	default:
		break;
	}

#endif  // CALIBRATION



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

  /* Activated macrophages synthesize new cytokines in quantities dependent on
   * the vocal treatment type. TODO(Kim): INSERT REFS? */
	this->amac_produce_cytokines();

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

