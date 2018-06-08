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
#include <array>

using namespace std;
int Fibroblast::numOfFibroblasts = 0;


float Fibroblast::cytokineSynthesis[32] = {10, 1, 2, 20, 1, 1, 1, 10, 0.5, 1, 5, 0.01, 1, 1, 1, 1, 10, 0.5, 1, 1, 1, 1, 10, 10, 1, 1, 1, 1, 0.05, 1, 5, 1};
float Fibroblast::activation[5] = {10.0, 50.0, 0, 25.0, 2.5};
float Fibroblast::ECMsynthesis[19] = {1, 0, 0.01, 50, 25, 2, 10, 5, 1, 0, 25, 2, 1, 0, 50, 5, 10, 12, 1};
float Fibroblast::proliferation[6] = {24, 10, 1, 0, 25, 3};

int threadSafeRand();

Fibroblast::Fibroblast() {
// DEBUG
/*
#pragma omp critical
{
   (Agent::agentWorldPtr->newfibs)++;
}
*/
	cout << "default fib alloc" << endl;
}

Fibroblast::Fibroblast(Patch* patchPtr) {
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
	this->ix[write_t] = patchPtr->indice[0];
	this->iy[write_t] = patchPtr->indice[1];
	this->iz[write_t] = patchPtr->indice[2];
	this->index[write_t] = patchPtr->index;
	this->alive[write_t] = true;
#ifdef DISABLE_RAND
	this->life[write_t] = WHWorld::reportTick(0, 8);  // Unactivated fibroblasts live for 5 to 11 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
#else
	this->life[write_t] = WHWorld::reportTick(rand_r(&(agentWorldPtr->seeds[tid]))%24, 5 + rand_r(&(agentWorldPtr->seeds[tid]))%7);  // Unactivated fibroblasts live for 5 to 11 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
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
	this->life[write_t] = WHWorld::reportTick(0, 8);  // Unactivated fibroblasts live for 5 to 11 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
#else
	this->life[write_t] = WHWorld::reportTick(0, 5 + rand_r(&(agentWorldPtr->seeds[tid]))%7);  // Unactivated fibroblasts live for 5 to 11 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
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
	int in = this->index[read_t];
	double hours = Agent::agentWorldPtr->reportHour();
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

  /* Unactivated fibroblasts only move along their preferred gradient and 
   * perform biological functions if there is damage. TODO(Kim): INSERT REF? */
	if (totaldamage == 0) {
		this->wiggle();
	} else {
    /*************************************************************************
     * PROLIFERATION                                                         *
     *************************************************************************/

		// Unactivated fibroblasts proliferate every 24 hours. TODO(Kim):INSERT REF?
#ifndef CALIBRATION
		if (fmod(hours, 24) == 0) {
#else  // CALIBRATION
		if (fmod((float) hours, Fibroblast::proliferation[0]) == 0) {
#endif  // CALIBRATION
			float meanTNF = this->meanNeighborChem(TNF);
			float meanTGF = this->meanNeighborChem(TGF);
			float meanFGF = this->meanNeighborChem(FGF);
			float meanIL1 = this->meanNeighborChem(IL1beta);
			int countfHA = this->countNeighborECM(fha);
			int TGFrelated = 0;
#ifndef CALIBRATION
			if (meanTGF <= 10) {
#else  // CALIBRATION
			if (meanTGF <= Fibroblast::proliferation[1]) {
#endif  // CALIBRATION
        // Low TGF (0.1-1nm) stimulate fib proliferation and attraction. TODO(Kim):INSERT REF?
				TGFrelated = 1;
			} else {
        // High TGF (1-10nm) inhibits proliferation. TODO(Kim):INSERT REF?
				TGFrelated = -1;
			}
#ifndef CALIBRATION
			float fibProlif = countfHA + log10(1 + meanTNF + TGFrelated*meanTGF + meanFGF + meanIL1 + countfHA);  // TODO(Kim):INSERT REF?
			if (rollDice(25 + fibProlif/2)) {  // TODO(Kim):INSERT REF?
#else  // CALIBRATION
			float fibProlif = countfHA + Fibroblast::proliferation[2]*(log10(1 + meanTNF + TGFrelated*meanTGF + meanFGF + meanIL1 + countfHA)) + Fibroblast::proliferation[3];  // TODO(Kim):INSERT REF?
			if (rollDice(Fibroblast::proliferation[4] + fibProlif/Fibroblast::proliferation[5])) {  // TODO(Kim):INSERT REF?
#endif  // CALIBRATION
				this->hatchnewfibroblast(2);
				this->die();
				return;
			}
		}

    /*************************************************************************
     * MOVEMENT                                                              *
     *************************************************************************/
		this->fibSniff();

    /*************************************************************************
     * ACTIVATION                                                            *
     *************************************************************************/
		/* An unactivated fibroblast can be activated if it is in the damage zone.
     * TODO(Kim): INSERT REF? */
//		if (Agent::agentPatchPtr[in].inDamzone == true) {
			// Low TGF promote and high TGF inhibit chances of fibroblast activation
#ifdef OPT_CHEM
			float patchTGF = agentWorldPtr->WHWorldChem->getPchem(TGF, in);
#else		// OPT_CHEM
			float patchTGF = agentWorldPtr->WHWorldChem->pTGF[in];
#endif		// OPT_CHEM
#ifndef CALIBRATION
			int chance1 =  50;
			int chance2 =   5;
			int chance3 =   1;	// work with chance4 to reduce chance to < 1
			int chance4 =  50;
			float TGFthresholdHigh = 10.0;//10;
			float TGFthresholdLow  =  0.001;
			float TGFbaseline = 0.0f;//agentWorldPtr->WHWorldChem->pbaseline[TGF];
			if ((patchTGF > TGFthresholdHigh && rollDice(chance1)) ||
				((patchTGF > TGFthresholdLow && patchTGF <= TGFthresholdHigh) && rollDice(chance2)) ||
//					(rollDice(chance3) && rollDice(chance4))) { // TODO(Kim): INSERT REF?
				((patchTGF > TGFbaseline && patchTGF <= TGFthresholdLow) && (rollDice(chance3) && rollDice(chance4))) ) {
#else  // CALIBRATION
//			if ((patchTGF > Fibroblast::activation[0] && rollDice(Fibroblast::activation[1])) || (patchTGF > Fibroblast::activation[2]) || (rollDice(Fibroblast::activation[3]))) { // TODO(Kim): INSERT REF?
			if ((patchTGF > Fibroblast::activation[0] && rollDice(Fibroblast::activation[4]))   // high TGF, low activation chance
			        || ((patchTGF > Fibroblast::activation[2]) &&
			                (patchTGF <= Fibroblast::activation[0]) &&      // Medium TGF
			                    rollDice(Fibroblast::activation[3])) ||     //  , medium activation chance
			                (patchTGF == 0 && rollDice(Fibroblast::activation[1]))) // Zero TGF, high activation chance
			{
#endif  // CALIBRATION
				this->fibActivation();
			}
//		}
	}

  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  // Unactivated fibroblasts can die naturally
	this->life[write_t] = this->life[read_t] - 1;
	if (this->life[read_t] <= 0) {
		this->die();
	}
}

void Fibroblast::hatchnewAfibroblast(int number) {
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
		int nb[27];
		for (int i = 0; i < 27; i++) {
			nb[i] = Agent::neighbor[i];
		}
		random_shuffle(&nb[0], &nb[26]);
		for (int i = 0; i < 27 && newfibs < number; i++) {
	#else
		int nb[8];
		for (int i = 0; i < 8; i++) {
			nb[i] = Agent::neighbor[i];
		}
		random_shuffle(&nb[0], &nb[7]);
		for (int i = 0; i < 8 && newfibs < number; i++) {
	#endif
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
				newfibs++;
				// Update target neighboring patch as occupied
				Agent::agentPatchPtr[in].occupiedby[write_t] = afibroblast;
				newfibroblast->fibActivation();
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
	int nb[27];
	for (int i = 0; i < 27; i++) {
		nb[i] = Agent::neighbor[i];
	}
	random_shuffle(&nb[0], &nb[26]);
	for (int i = 0; i < 27 && newfibs < number; i++) {
#else
	int nb[8];
	for (int i = 0; i < 8; i++) {
		nb[i] = Agent::neighbor[i];
	}
	random_shuffle(&nb[0], &nb[7]);
	for (int i = 0; i < 8 && newfibs < number; i++) {
#endif
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
			newfibs++;
			// Update target neighboring patch as occupied
			Agent::agentPatchPtr[in].occupiedby[write_t] = fibroblast;
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

void Fibroblast::fibSniff() {
  //TODO(Nuttiya): Why do we check for damage in this function? We always check before the function call. And Macrophage::macSniff & Neutrophil::neuSniff don't check for damage.
	int in = this->index[read_t];
	if ((Agent::agentPatchPtr[in]).inDamzone == true) {
		if (this->moveToHighestChem(FIBgrad) != true) this->wiggle();
	} else {
		this->wiggle();
	}

	// FGF can excite fibroblast and overcome gradient TODO(Kim): INSERT REF?
	if ((this->meanNeighborChem(FGF)) > 0) {
		this->wiggle();
	}
}

void Fibroblast::fibSniffOutDam() {
        if (this->moveToHighestChem(FIBgrad) != true) this->wiggle();

//    // FGF can excite fibroblast and overcome gradient TODO(Kim): INSERT REF?
//    if ((this->meanNeighborChem(FGF)) > 0) {
//        this->wiggle();
//    }
}

// DEBUG coll syn
void Fibroblast::fibSniff2() {

    // Location of fibroblast in x,y,z dimensions of the world
      int x = this->ix[read_t];
      int y = this->iy[read_t];
      int z = this->iz[read_t];
    // Number of patches in x,y,z dimensions of world
      int nx = Agent::nx;
      int ny = Agent::ny;
      int nz = Agent::nz;

    vector<array<int, 4>> damagedneighbors;
    vector<array<int, 4>> healthyneighbors;

    vector<int> damagedindex;
    vector<int> healthyindex;

    int dx, dy, dz, in;

    vector<int>::iterator it;

//    bool isInDam = Agent::agentPatchPtr[in].isInDamZone();

    // Make a list of damaged neighboring patches
#ifndef MODEL_3D
    for (int i = 9; i < 18; i++) {
#else
    for (int i = 0; i < 27; i++) {
#endif
        dx = Agent::dX[i];
        dy = Agent::dY[i];
        dz = Agent::dZ[i];
        in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

        // DEBUG coll syn
        if (in == this->index[read_t]) continue;
    // Try a new neighboring patch if this one is outside the world dimensions.
        if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
    // Add the valid damaged neighboring patch to the list.
        if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back({dx, dy, dz, in});
        else healthyneighbors.push_back({dx, dy, dz, in});
    }

    // Move to any of the damaged neighbors
    const int ntargets1 = damagedneighbors.size();
    const int ntargets2 = healthyneighbors.size();
    for (int i = 0; i < ntargets1; i++) damagedindex.push_back(i);
    for (int i = 0; i < ntargets2; i++) healthyindex.push_back(i);

    // shuffle damaged neighbors
    // TODO: Use my random

//    if (!isInDam)
//        std::random_shuffle(damagedindex.begin(), damagedindex.end(), &threadSafeRand);
    for (it = damagedindex.begin(); it != damagedindex.end(); it++)
    {
        if (this->moveToTarget(damagedneighbors[*it][0], damagedneighbors[*it][1],
                damagedneighbors[*it][2], damagedneighbors[*it][3])) return;
    }

    // shuffle healthy neighbors
    // TODO: Use my random
//    if (!isInDam)
//        std::random_shuffle(healthyindex.begin(), healthyindex.end(), &threadSafeRand);
    for (it = healthyindex.begin(); it != healthyindex.end(); it++)
    {
        if (this->moveToTarget(healthyneighbors[*it][0], healthyneighbors[*it][1],
                healthyneighbors[*it][2], healthyneighbors[*it][3])) return;
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
	this->life[read_t] = WHWorld::reportTick(0, 5 + rand_r(&(agentWorldPtr->seeds[tid]))%7);  // Unactivated fibroblasts live for 5 to 11 days. 0 corresponds to hours. TODO(Kim:) INSERT REF
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

	Fibroblast::numOfFibroblasts++;

  // Assigns new Fibroblast to this patch if it is unoccupied
	if (Agent::agentPatchPtr[in].isOccupied() == false) {
		Agent::agentPatchPtr[in].setOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = this->type[read_t];
	} else {
		cout << "error in hatching and initialization!!!" << dx << " " << dy << endl;
	}
}

void Fibroblast::afib_cellFunction() {
    int in = this->index[read_t];

    double hours = Agent::agentWorldPtr->reportHour();

  /*************************************************************************
   * MOVEMENT                                                              *
   *************************************************************************/
  /* Activated fibroblasts only move along their preferred gradient if there
   * is damage. TODO(Kim): Insert ref? */
  int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];
  // DEBUG coll syn
  bool isInDam = Agent::agentPatchPtr[in].isInDamZone();
  if (totaldamage != 0) {
//      if (isInDam) {
          this->fibSniff();//this->fibSniff2();//
//      } else {
//          this->fibSniffOutDam();   // follow gradient
//      }
  }

  /*************************************************************************
   * PROLIFERATION                                                         *
   *************************************************************************/

//	// activated fibroblasts proliferate every 24 hours. TODO(Kim):INSERT REF?
//#ifndef CALIBRATION
//	if (fmod(hours, 24) == 0) {
//#else  // CALIBRATION
//	if (fmod((float) hours, Fibroblast::proliferation[0]) == 0) {
//#endif  // CALIBRATION
//		float meanTNF = this->meanNeighborChem(TNF);
//		float meanTGF = this->meanNeighborChem(TGF);
//		float meanFGF = this->meanNeighborChem(FGF);
//		float meanIL1 = this->meanNeighborChem(IL1beta);
//		int countfHA = this->countNeighborECM(fha);
//		int TGFrelated = 0;
//#ifndef CALIBRATION
//		if (meanTGF <= 10) {
//#else  // CALIBRATION
//		if (meanTGF <= Fibroblast::proliferation[1]) {
//#endif  // CALIBRATION
//      // Low TGF (0.1-1nm) stimulate fib proliferation and attraction. TODO(Kim):INSERT REF?
//			TGFrelated = 1;
//		} else {
//      // High TGF (1-10nm) inhibits proliferation. TODO(Kim):INSERT REF?
//			TGFrelated = -1;
//		}
//#ifndef CALIBRATION
//		float fibProlif = countfHA + log10(1 + meanTNF + TGFrelated*meanTGF + meanFGF + meanIL1 + countfHA);  // TODO(Kim):INSERT REF?
//		if (rollDice(25 + fibProlif/2)) {  // TODO(Kim):INSERT REF?
//#else  // CALIBRATION
//		float fibProlif = countfHA + Fibroblast::proliferation[2]*(log10(1 + meanTNF + TGFrelated*meanTGF + meanFGF + meanIL1 + countfHA)) + Fibroblast::proliferation[3];  // TODO(Kim):INSERT REF?
//		if (rollDice(Fibroblast::proliferation[4] + fibProlif/Fibroblast::proliferation[5])) {  // TODO(Kim):INSERT REF?
//#endif  // CALIBRATION
//			this->hatchnewAfibroblast(2);
//			this->die();
//			return;
//		}
//	}


  /*************************************************************************
   * ECM PROTEIN & CHEMICAL SYNTHESIS                                      *
   *************************************************************************/
	// Calculates chemical gradients and patch chemical concentrations
  float meanTNF = this->meanNeighborChem(TNF);
  float meanTGF = this->meanNeighborChem(TGF);
  float meanFGF = this->meanNeighborChem(FGF);
  float meanIL1 = this->meanNeighborChem(IL1beta);
  float meanIL6 = this->meanNeighborChem(IL6);
  float meanIL8 = this->meanNeighborChem(IL8);
  int countnHA = this->countNeighborECM(nha);
  int countfHA = this->countNeighborECM(fha);
#ifdef OPT_CHEM
  float patchTNF  = this->agentWorldPtr->WHWorldChem->getPchem(TNF,in) ;
  float patchTGF  = this->agentWorldPtr->WHWorldChem->getPchem(TGF,in) ;
  float patchIL6  = this->agentWorldPtr->WHWorldChem->getPchem(IL6,in) ;
  float patchIL10 = this->agentWorldPtr->WHWorldChem->getPchem(IL10,in );
#else		// OPT_CHEM
  float patchTNF  = this->agentWorldPtr->WHWorldChem->pTNF[in];
  float patchTGF  = this->agentWorldPtr->WHWorldChem->pTGF[in];
  float patchIL6  = this->agentWorldPtr->WHWorldChem->pIL6[in];
  float patchIL10 = (this->agentWorldPtr->WHWorldChem->pIL10[in]);
#endif		// OPT_CHEM

  // Makes collagen and elastin every 12 hours. TODO(Kim): INSERT REF?

  // Every 30 minute
//#ifndef CALIBRATION
  int makeCollHours = 24;//12;//1;
  if (fmod(((Agent::agentWorldPtr)->reportHour()), makeCollHours) == 0)
//#else  // CALIBRATION
//  if (fmod((float)((Agent::agentWorldPtr)->reportHour()), Fibroblast::ECMsynthesis[17]) == 0)
//#endif  // CALIBRATION
  {
  // DEBUG coll syn
//    this->makeOCollagen2();
    this->makeOCollagen(meanTGF, meanFGF, meanIL1, meanIL6, meanIL8, countnHA, countfHA);
    this->makeOElastin(meanTNF, meanTGF, meanFGF, meanIL1);
  }

  // Makes hyaluronan every hour. TODO(Kim): INSERT REF?
#ifndef CALIBRATION
  if (fmod((Agent::agentWorldPtr)->reportHour(), 1) == 0)
#else  // CALIBRATION
  if (fmod((float)((Agent::agentWorldPtr)->reportHour()), Fibroblast::ECMsynthesis[18]) == 0)
#endif  // CALIBRATION
  {
    this->makeHyaluronan(meanTNF, meanTGF, meanFGF, meanIL1);
  }

  //Change in chemicals due to cells. TODO(Kim): INSERT REF?
#ifndef CALIBRATION
  float factor = 0.0001;
  float factorFGF = 0.0001;
  float factorTGF = 0.05;
  float TGFinc = 10 + 0.5*(1 + patchTNF + patchIL10);
  float FGFinc = 5;

  TGFinc *= factorTGF; // = 0.00000000f; //
  FGFinc *= factorFGF;

#else  // CALIBRATION
  float TGFinc = Fibroblast::cytokineSynthesis[7] +
                 Fibroblast::cytokineSynthesis[8]*
                  (Fibroblast::cytokineSynthesis[9] + patchTNF + patchIL10);
  float FGFinc = Fibroblast::cytokineSynthesis[10];
#endif  // CALIBRATION

#ifdef OPT_CHEM
  this->agentWorldPtr->WHWorldChem->incDchem(TGF, in, TGFinc);
  this->agentWorldPtr->WHWorldChem->incDchem(FGF, in, FGFinc);
#else		// OPT_CHEM
  (this->agentWorldPtr->WHWorldChem->dTGF[in]) += TGFinc;
  (this->agentWorldPtr->WHWorldChem->dFGF[in]) += FGFinc;
#endif		// OPT_CHEM 

  /* Activated fibroblasts synthesize new cytokines in quantities dependent on
   * the vocal treatment type. TODO(Kim): INSERT REFS? */
  float TNFinc = 0.0;
  float IL6inc = 0.0;
  float IL8inc = 0.0;
#ifndef CALIBRATION
  if (this->agentWorldPtr->treatmentOption == voicerest) {
    TNFinc = 10/(1 + patchTGF + patchIL10*2 + countnHA);
    IL6inc = 0.01*((1 + patchTNF)/(1 + patchIL10));
    IL8inc = 1/(1 + patchIL10 + countnHA);
  } else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
    TNFinc = 20/(1 + patchTGF + patchIL6 + patchIL10 + countnHA);
    IL8inc = 0.05/(1 + patchIL10 + countnHA);
    if (((Agent::agentWorldPtr)->reportDay()) <= 7) {
      IL6inc = (1 + ((Agent::agentWorldPtr)->RVVS)*10);
    } else {
      IL6inc = 0.5*((1 + patchTNF)/(1 + patchIL10));
    }
  } else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
    TNFinc = 1/(1 + patchTGF + patchIL6 + patchIL10 + countnHA);
    IL8inc = 5/(1 + patchIL10 + countnHA);
    if (((Agent::agentWorldPtr)->reportDay()) <= 7) {
      IL6inc = 1 + ((Agent::agentWorldPtr)->SSVS)*10;
    } else { 
      IL6inc = 10*((1 + patchTNF)/(1 + patchIL10));
    } 
  }

  TNFinc *= factor;
  IL6inc *= factor;
  IL8inc *= factor;

#else	// CALIBRATION
  if (this->agentWorldPtr->treatmentOption == voicerest) {
    TNFinc = Fibroblast::cytokineSynthesis[0]/(Fibroblast::cytokineSynthesis[1] +
             patchTGF + patchIL10*Fibroblast::cytokineSynthesis[2] + countnHA);
    IL6inc = Fibroblast::cytokineSynthesis[11]*(Fibroblast::cytokineSynthesis[12] + patchTNF)/
            (Fibroblast::cytokineSynthesis[13] + patchIL10);
    IL8inc = Fibroblast::cytokineSynthesis[26]/
            (Fibroblast::cytokineSynthesis[27] + patchIL10 + countnHA);
  } else if (this->agentWorldPtr->treatmentOption == resonantvoice) {
    TNFinc = Fibroblast::cytokineSynthesis[3]/
            (Fibroblast::cytokineSynthesis[4] + patchTGF + patchIL6 + patchIL10 + countnHA); 
    IL8inc = Fibroblast::cytokineSynthesis[28]/
            (Fibroblast::cytokineSynthesis[29] + patchIL10 + countnHA);
    if (((Agent::agentWorldPtr)->reportDay()) <= 7) {
      IL6inc = Fibroblast::cytokineSynthesis[14]*(Fibroblast::cytokineSynthesis[15] +
             ((Agent::agentWorldPtr)->RVVS)*Fibroblast::cytokineSynthesis[16]);
    } else {
      IL6inc = Fibroblast::cytokineSynthesis[17]*(Fibroblast::cytokineSynthesis[18] + patchTNF)/
              (Fibroblast::cytokineSynthesis[19] + patchIL10);
    }
  } else if (this->agentWorldPtr->treatmentOption == spontaneousspeech) {
    TNFinc = Fibroblast::cytokineSynthesis[5]/
            (Fibroblast::cytokineSynthesis[6] + patchTGF + patchIL6 + patchIL10 + countnHA);
    IL8inc = Fibroblast::cytokineSynthesis[30]/
            (Fibroblast::cytokineSynthesis[31] + patchIL10 + countnHA);
    if (((Agent::agentWorldPtr)->reportDay()) <= 7) {
      IL6inc = Fibroblast::cytokineSynthesis[20]*(Fibroblast::cytokineSynthesis[21] +
             ((Agent::agentWorldPtr)->SSVS)*Fibroblast::cytokineSynthesis[22]);
    } else {
      IL6inc = Fibroblast::cytokineSynthesis[23]*(Fibroblast::cytokineSynthesis[24] + patchTNF)/
              (Fibroblast::cytokineSynthesis[25] + patchIL10);
    }
  }
#endif

#ifdef OPT_CHEM
  this->agentWorldPtr->WHWorldChem->incDchem(TNF, in, TNFinc);
  this->agentWorldPtr->WHWorldChem->incDchem(IL8, in, IL8inc);
  this->agentWorldPtr->WHWorldChem->incDchem(IL6, in, IL6inc);
#else		// OPT_CHEM    
  (this->agentWorldPtr->WHWorldChem->dTNF[in]) += TNFinc;
  (this->agentWorldPtr->WHWorldChem->dIL6[in]) += IL6inc;
  (this->agentWorldPtr->WHWorldChem->dIL8[in]) += IL8inc;
#endif		// OPT_CHEM

#ifdef PRINT_SECRETION
  int x = this->ix[read_t];
  int y = this->iy[read_t];
  int z = this->iz[read_t];
  printCytRelease(3, TGF, x, y, z, TGFinc);
  printCytRelease(3, FGF, x, y, z, FGFinc);
  printCytRelease(3, TNF, x, y, z, TNFinc);
  printCytRelease(3, IL6, x, y, z, IL6inc);
  printCytRelease(3, IL8, x, y, z, IL8inc);
#endif  // PRINT_SECRETION

  /*************************************************************************
   * DEACTIVATION                                                          *
   *************************************************************************/
  /* Activated fibroblasts might be deactivated once the damage is cleared. TODO(Kim): INSERT REF? */
	totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	// DEBUG
	// Change to percentage of initial damage
	int damageThreshold = 0;//5000;

#ifndef CALIBRATION
	// DEBUG
	if (totaldamage == 0 && rollDice(2.5)) this->fibDeactivation();
//	if (totaldamage < damageThreshold && rollDice(2.5)) this->fibDeactivation();
//	if (totaldamage == damageThreshold && rollDice(2.5)) this->fibDeactivation();
#else  // CALIBRATION
	if (totaldamage == 0 && rollDice(Fibroblast::activation[4])) this->fibDeactivation();
#endif  // CALIBRATION

  /*************************************************************************
   * DEATH                                                                 *
   *************************************************************************/
  // Activated fibroblasts can die naturally
	this->life[write_t] = this->life[read_t] - 1;
	if (this->life[read_t] <= 0) {
		this->die();
	}
}

// TODO: Move to util
int threadSafeRand()
{
    int tid = 0;
#ifdef _OMP
  // Get thread id in order to access the seed that belongs to this thread
    tid = omp_get_thread_num();
#endif


    return rand_r(&((Agent::agentWorldPtr)->seeds[tid]));

}



int randProductionColl()
{
	int range = fibProductionTCollUB - fibProductionTCollLB;
	int ret = fibProductionTCollLB + (threadSafeRand()%range);

	return ret;
}


void Fibroblast::makeOCollagen2() {


    // DEBUG coll syn
    if (Agent::agentPatchPtr[this->index[write_t]].damage[read_t] != 0)
    {
        Agent::agentECMPtr[this->index[write_t]].addColls(randProductionColl());
    }

}



void Fibroblast::makeOCollagen(float meanTGF, float meanFGF, float meanIL1, float meanIL6, float meanIL8, int countnHA, int countfHA) {
	int read_index;
	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		// If it has, work off of the intermediate value
		read_index = write_t;
	} else {
		// If it has NOT, work off of the original value
		read_index = read_t;
	}
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

	// DEBUG coll syn
	if (Agent::agentPatchPtr[this->index[read_t]].damage[read_t] != 0)
	{
	    Agent::agentECMPtr[this->index[read_t]].addColls(randProductionColl());
	}


	// Make a list of damaged neighboring patches
#ifndef MODEL_3D
	for (int i = 9; i < 18; i++) {
#else
	for (int i = 0; i < 27; i++) {
#endif
		dx = Agent::dX[i];
		dy = Agent::dY[i];
		dz = Agent::dZ[i];
		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// DEBUG coll syn
		if (in == this->index[read_t]) continue;
    // Try a new neighboring patch if this one is outside the world dimensions.
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
    // Add the valid damaged neighboring patch to the list.
		if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back(i);
	}

	// Target a random damaged neighboring patch, if there are any.
	if (damagedneighbors.size() > 0) {
		int tid = 0;
#ifdef _OMP
    // Get thread id in order to access the seed that belongs to this thread
		tid = omp_get_thread_num();
#endif
		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
		target = damagedneighbors[randInt];
		dx = Agent::dX[target];
		dy = Agent::dY[target];
		dz = Agent::dZ[target];

		/* Based on number of fHA, nHA, chance, chemical concentrations, move to 
     * new patch and sprout ocollagen */
#ifndef CALIBRATION
		int stimulation = log10(1 + (meanTGF + meanIL6)/(1 + meanFGF + meanIL1 + meanIL8*0.01)); //TODO(Kim): INSERT REFS?
		if ((countfHA > countnHA && rollDice(50 + stimulation)) ||
				(countfHA == countnHA && rollDice(25 + stimulation/2)) ||
				(rollDice(10+stimulation/5))) {
#else  // CALIBRATION
		int stimulation = (Fibroblast::ECMsynthesis[0]*(log10(1 + (meanTGF + meanIL6)) + Fibroblast::ECMsynthesis[1])/(1 + meanFGF + meanIL1 + meanIL8*Fibroblast::ECMsynthesis[2])); //TODO(Kim): INSERT REFS?
		if ((countfHA > countnHA && rollDice(Fibroblast::ECMsynthesis[3] + stimulation)) ||
				(countfHA == countnHA && rollDice(Fibroblast::ECMsynthesis[4] + stimulation/Fibroblast::ECMsynthesis[5])) ||
				//if (rollDice(100)){
				(rollDice(Fibroblast::ECMsynthesis[6] + stimulation/Fibroblast::ECMsynthesis[7]))) {
#endif  // CALIBRATION
			in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
			this->move(dx, dy, dz, read_index);
			// TODO: accelerate with TGF
			// DEBUG - verification
			Agent::agentECMPtr[in].addColls(10 + rand_r(&(agentWorldPtr->seeds[tid]))%20);
//			Agent::agentECMPtr[in].addColls(1 + rand_r(&(agentWorldPtr->seeds[tid]))%2);
//			Agent::agentECMPtr[in].addColls(randProductionColl());
			// debug vis
//			if (z == 72)
//				Agent::agentWorldPtr->incECM(this->index[read_t], m_col, 1000);//Agent::agentWorldPtr->maxFCol);


/*
#ifdef _OMP
                         Agent::agentECMPtr[in].lock();     // collagen ecm mutex --------- begin
#endif
			Agent::agentECMPtr[in].ocollagen[write_t] = Agent::agentECMPtr[in].ocollagen[read_t] + 1 + rand_r(&(agentWorldPtr->seeds[tid]))%2;
#ifdef OPT_ECM
                        Agent::agentECMPtr[in].set_dirty(); 
#endif

#ifdef _OMP
                         Agent::agentECMPtr[in].unlock();   // collagen ecm mutex --------- end
#endif
*/
		}
	}
}

void Fibroblast::makeOElastin(float meanTNF, float meanTGF, float meanFGF, float meanIL1) {
	int read_index;
	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		// If it has, work off of the intermediate value
		read_index = write_t;
	} else {
		// If it has NOT, work off of the original value
		read_index = read_t;
	}
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

	// Make list of damaged neighboring patches
#ifndef MODEL_3D
	for (int i = 9; i < 18; i++) {
#else
	for (int i = 0; i < 27; i++) {
#endif
		dx = Agent::dX[i];
		dy = Agent::dY[i];
		dz = Agent::dZ[i];
		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
    // Try a new neighboring patch if this one is outside the world dimensions.
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
    // Add the valid damaged neighboring patch to the list.
		if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back(i);
	}

	// Target a random damaged neighboring patch, if there are any.
	if (damagedneighbors.size() > 0) {
		int tid = 0;
#ifdef _OMP
    // Get thread id in order to access the seed that belongs to this thread
		tid = omp_get_thread_num();
#endif
		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
		target = damagedneighbors[randInt];
		dx = Agent::dX[target];
		dy = Agent::dY[target];
		dz = Agent::dZ[target];
		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// Based on mean TGf, IL1, TNF, chance, move to new patch and sprout oelastin. TODO(Kim): INSERT REFS?
#ifndef CALIBRATION
		int stimulation = log10((1 + meanTGF)/(1 + meanFGF + meanIL1 + meanTNF));
		if (rollDice(25 + stimulation/2)){
#else  // CALIBRATION
		int stimulation = (Fibroblast::ECMsynthesis[8]*(log10((1 + meanTGF)) + Fibroblast::ECMsynthesis[9])/(1 + meanFGF + meanIL1 + meanTNF));
		if (rollDice(Fibroblast::ECMsynthesis[10] + stimulation/Fibroblast::ECMsynthesis[11])) {
#endif  // CALIBRATION
			this->move(dx, dy, dz, read_index);
#ifdef _OMP
			 Agent::agentECMPtr[in].lock();    // elastin ecm mutex --------- begin
#endif
			 // DEBUG - verification
			 Agent::agentECMPtr[in].oelastin[write_t] = Agent::agentECMPtr[in].oelastin[read_t] + 5 + rand_r(&(agentWorldPtr->seeds[tid]))%20;
			 //			Agent::agentECMPtr[in].oelastin[write_t] = Agent::agentECMPtr[in].oelastin[read_t] + 1 + rand_r(&(agentWorldPtr->seeds[tid]))%2;
#ifdef OPT_ECM
                        Agent::agentECMPtr[in].set_dirty(); 
#endif
#ifdef _OMP
                        Agent::agentECMPtr[in].unlock();   // elastin ecm mutex --------- end
#endif

//                        cout << " sprout new elastin at " << in << endl; 
		}
	}
}

void Fibroblast::makeHyaluronan(float meanTNF, float meanTGF, float meanFGF, float meanIL1) {

	int read_index;
	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		// If it has, work off of the intermediate value
		read_index = write_t;
	} else {
		// If it has NOT, work off of the original value
		read_index = read_t;
	}
  // Location of fibroblast in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;
	vector <int> damagedneighbors;
#ifndef CALIBRATION
	int stimFactor = 1.0;
	int stimulation = stimFactor*(log10(1 + meanTGF + meanFGF + meanTNF + meanIL1)); // TODO(Kim): INSERT REF?
#else  // CALIBRATION
	int stimulation = Fibroblast::ECMsynthesis[12]*(log10(1 + meanTGF + meanFGF + meanTNF + meanIL1)) + Fibroblast::ECMsynthesis[13]; // TODO(Kim): INSERT REF?
#endif  // CALIBRATION
	int dx, dy, dz, in, randInt, target;

	// Make list of damaged neighboring patches
#ifndef MODEL_3D
	for (int i = 9; i < 18; i++) {
#else
	for (int i = 0; i < 27; i++) {
#endif
		dx = Agent::dX[i];
		dy = Agent::dY[i];
		dz = Agent::dZ[i];
		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
    // Try a new neighboring patch if this one is outside the world dimensions.
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
    // Add the valid damaged neighboring patch to the list.
		if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back(i);
	}

	// Target a random damaged neighboring patch, if there are any.
	if (damagedneighbors.size() > 0) {
		int tid = 0;
#ifdef _OMP
    // Get thread id in order to access the seed that belongs to this thread
		tid = omp_get_thread_num();
#endif
		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
		target = damagedneighbors[randInt];
		dx = Agent::dX[target];
		dy = Agent::dY[target];
		dz = Agent::dZ[target];
		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// Based on mean TGF, FGF, TNF, IL1, chance, move to new patch and sprout HA
#ifndef CALIBRATION
		int activationEnhancement = 0; // 3;
		if (rollDice(activationEnhancement + stimulation)) { // TODO(Kim): INSERT REF?
//		if (rollDice(50 + stimulation)) { // TODO(Kim): INSERT REF?
#else  // CALIBRATION
		if (rollDice(Fibroblast::ECMsynthesis[14] + stimulation)) { // TODO(Kim): INSERT REF?
#endif  // CALIBRATION
			this->move(dx, dy, dz, read_index);
			Agent::agentECMPtr[in].addHAs(1);
		}
#ifndef CALIBRATION
	} else if (rollDice(5 + stimulation/10)) {
#else  // CALIBRATION
//	}else if (rollDice(Fibroblast::ECMsynthesis[15] + stimulation/Fibroblast::ECMsynthesis[16])) {
	}else if (rollDice(5 + stimulation/10)) {
#endif  // CALIBRATION
	        in = this->index[read_t];
		Agent::agentECMPtr[in].addHAs(1);
	}
}

