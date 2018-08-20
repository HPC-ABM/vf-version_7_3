/* 
 * File: Agent.cpp
 * 
 * File Contents: Contains the Agent class.
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 * 
 * Created on May 19, 2013, 7:37 AM
 */
/*****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#include "Agent.h"
#include "../World/Usr_World/woundHealingWorld.h"
#include "../enums.h"
#include <iostream>
#include <vector>

WHWorld* Agent::agentWorldPtr = NULL;
Patch* Agent::agentPatchPtr = NULL;
ECM* Agent::agentECMPtr = NULL;
float Agent::baseline[TOTAL_CHEM] = {0.0f};
int Agent::nx = 0;
int Agent::ny = 0;
int Agent::nz = 0;
int Agent::dZ[27] = {-1,  0,  1, -1,  0,  1, -1,  0,  1, -1,  0,  1, -1, 0, 1, -1, 0, 1, -1,  0,  1, -1, 0, 1, -1, 0, 1};    //was x
int Agent::dY[27] = {-1, -1, -1,  0,  0,  0,  1,  1,  1, -1, -1, -1,  0, 0, 0,  1, 1, 1, -1, -1, -1,  0, 0, 0,  1, 1, 1};     //was y
int Agent::dX[27] = {-1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0, 0, 0,  0, 0, 0,  1,  1,  1,  1, 1, 1,  1, 1, 1};     //was z
#ifdef MODEL_3D         //3D case
int Agent::neighbor[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
#else	                // 2D case
// We do not include neighbor 13 (no movement to self (0,0,0))
int Agent::neighbor[8] = {9, 10, 11, 12, 14, 15, 16, 17};
#endif

using namespace std;

Agent::Agent() {
}

Agent::~Agent() {
}

bool Agent::isAlive() {
	if (this->life[write_t] <= 0) this->alive[write_t] = false;
	else this->alive[write_t] = true;
	return this->alive[read_t];
} 

void Agent::cellFunction() {
}

void Agent::copyAndInitialize(Agent* original, int dx, int dy, int dz) {
}

bool Agent::isActivated()
{
	return this->activate[read_t];
}

void Agent::printCytRelease(int celltype, int chemtype, int x, int y, int z, float level) {
  //printf("\t[%d,%d,%d], %s, %s, %f\n", x, y, z, celltype, chemtype, level);
  // SUM
//  this->agentWorldPtr->chemSecreted[celltype][chemtype] += level;

  // MAX
  if (this->agentWorldPtr->chemSecreted[celltype][chemtype] < level)
    {
      this->agentWorldPtr->chemSecreted[celltype][chemtype] = level;
      this->agentWorldPtr->chemSecretedCoord[celltype][chemtype][0] = x;
      this->agentWorldPtr->chemSecretedCoord[celltype][chemtype][1] = y;
      this->agentWorldPtr->chemSecretedCoord[celltype][chemtype][2] = z;
    }
}

void die() {
}

int Agent::getX() {
	return this->ix[read_t];
}

int Agent::getY() {
	return this->iy[read_t];
}

int Agent::getZ() {
	return this->iz[read_t];
}

int Agent::getIndex() {
	return this->index[read_t];
}

//bool Agent::rollDice(float percent) {
bool Agent::rollDice(int percent) {
	int tid = 0;
#ifdef _OMP
        // Get thread id in order to access the seed that belongs to this thread
	tid = omp_get_thread_num();
#endif
	int randNum = rand_r(&(agentWorldPtr->seeds[tid]))%100;
	if (randNum < percent) {
		return 1;
	} else {
		return 0;
	}
}

float Agent::calc_chem(
    REGULATORS_T  uregs,
    REGULATORS_T  dregs,
    REGULATORS_T ucoefs,
    REGULATORS_T dcoefs,
    float offset) {
  assert(uregs.size() == ucoefs.size());
  assert(dregs.size() == dcoefs.size());

  float inc_num = 0; // numerator
  float inc_den = 0; // denominator
  int   uci = 0;
  int   dci = 0;

  for (REGULATORS_T::iterator it = uregs.begin(); it != uregs.end(); it++) {
      // numerator = numerator + coeff[i]*uregs[i]
      inc_num += ucoefs[uci++]*(*it);
  }

  for (REGULATORS_T::iterator it = dregs.begin(); it != dregs.end(); it++) {
      // denominator = denominator + coeff[i]*dregs[i]
      inc_den += dcoefs[dci++]*(*it);
  }

  // Check for zero denom
  if (inc_den == 0.0f) return 0.0f;

  return (inc_num/inc_den) + offset;
}

float Agent::calc_chem(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset) {
  assert(uregs.size() + dregs.size() == coefs.size());

  float inc_num = 0; // numerator
  float inc_den = 0; // denominator
  int   ci = 0;

  for (REGULATORS_T::iterator it = uregs.begin(); it != uregs.end(); it++) {
      // numerator = numerator + coeff[i]*uregs[i]
      inc_num += coefs[ci++]*(*it);
  }

  for (REGULATORS_T::iterator it = dregs.begin(); it != dregs.end(); it++) {
      // denominator = denominator + coeff[i]*dregs[i]
      inc_den += coefs[ci++]*(*it);
  }

  // Check for zero denom
  if (inc_den == 0.0f) return 0.0f;

  return (inc_num/inc_den) + offset;
}

float Agent::produce_tnf(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T coefs,
		float offset){

	int in = this->getIndex();
	float TNFinc = this->calc_chem(uregs, dregs, coefs, offset);
	// Update chem change
	(this->agentWorldPtr->WHWorldChem->dTNF[in]) += TNFinc;
	return TNFinc;
}

float Agent::produce_tnf(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
		float offset){

	int in = this->getIndex();
	float TNFinc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
	// Update chem change
	(this->agentWorldPtr->WHWorldChem->dTNF[in]) += TNFinc;
	return TNFinc;
}

float Agent::produce_tgf(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset){

  int in = this->getIndex();
  float TGFinc = this->calc_chem(uregs, dregs, coefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dTGF[in]) += TGFinc;
  return TGFinc;
}

float Agent::produce_tgf(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
		float offset){

	int in = this->getIndex();
	float TGFinc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
	// Update chem change
	(this->agentWorldPtr->WHWorldChem->dTGF[in]) += TGFinc;
	return TGFinc;
}


float Agent::produce_fgf(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset){

  int in = this->getIndex();
  float FGFinc = this->calc_chem(uregs, dregs, coefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dFGF[in]) += FGFinc;
  return FGFinc;
}

float Agent::produce_fgf(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
		float offset){

	int in = this->getIndex();
	float FGFinc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
	// Update chem change
	(this->agentWorldPtr->WHWorldChem->dFGF[in]) += FGFinc;
	return FGFinc;
}


float Agent::produce_mmp8(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T coefs,
		float offset){

	int in = this->getIndex();
	float MMP8inc = this->calc_chem(uregs, dregs, coefs, offset);
	// Update chem change
	(this->agentWorldPtr->WHWorldChem->dMMP8[in]) += MMP8inc;
	return MMP8inc;
}

float Agent::produce_mmp8(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
		float offset){

	int in = this->getIndex();
	float MMP8inc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
	// Update chem change
	(this->agentWorldPtr->WHWorldChem->dMMP8[in]) += MMP8inc;
	return MMP8inc;
}


float Agent::produce_il1(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset){

  int in = this->getIndex();
  float IL1inc = this->calc_chem(uregs, dregs, coefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dIL1beta[in]) += IL1inc;
  return IL1inc;
}

float Agent::produce_il1(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
		float offset){

	int in = this->getIndex();
	float IL1inc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
	// Update chem change
	(this->agentWorldPtr->WHWorldChem->dIL1beta[in]) += IL1inc;
	return IL1inc;
}

float Agent::produce_il6(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset){

  int in = this->getIndex();
  float IL6inc = this->calc_chem(uregs, dregs, coefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dIL6[in]) += IL6inc;
  return IL6inc;
}

float Agent::produce_il6(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
		float offset){

  int in = this->getIndex();
  float IL6inc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dIL6[in]) += IL6inc;
  return IL6inc;
}


float Agent::produce_il8(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset){

  int in = this->getIndex();
  float IL8inc = this->calc_chem(uregs, dregs, coefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dIL8[in]) += IL8inc;
  return IL8inc;
}

float Agent::produce_il8(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
    float offset){

  int in = this->getIndex();
  float IL8inc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dIL8[in]) += IL8inc;
  return IL8inc;
}


float Agent::produce_il10(
    REGULATORS_T uregs,
    REGULATORS_T dregs,
    REGULATORS_T coefs,
    float offset){

  int in = this->getIndex();
  float IL10inc = this->calc_chem(uregs, dregs, coefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dIL10[in]) += IL10inc;
  return IL10inc;
}

float Agent::produce_il10(
		REGULATORS_T uregs,
		REGULATORS_T dregs,
		REGULATORS_T ucoefs,
		REGULATORS_T dcoefs,
    float offset){

  int in = this->getIndex();
  float IL10inc = this->calc_chem(uregs, dregs, ucoefs, dcoefs, offset);
  // Update chem change
  (this->agentWorldPtr->WHWorldChem->dIL10[in]) += IL10inc;
  return IL10inc;
}


bool Agent::move(int dX, int dY, int dZ, int read_index) {
  // Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;
  // Location of target patch in x,y,z dimensions of world.
	int targetX = x + dX;
	int targetY = y + dY;
	int targetZ = z + dZ;

  /* Abort movement if trying to move to own patch, or to a patch outside world
   * dimensions, or to the patch the agent was at at the beginning of the tick,
   * or to an occupied patch. */
	if (dX == 0 && dY== 0 && dZ == 0) return false;
	if (targetX < 0 || targetX >= nx || targetY < 0 || targetY >= ny || targetZ < 0 || targetZ >= nz) return false;
	int newIndex = targetX + targetY*nx + targetZ*nx*ny;
	if (newIndex == this->index[read_t]) return false;
  /* If patch is unoccupied, move current instance to new patch at
   * (x + dX, y + dY, z + dZ). setOccupied() sets the new patch as occupied and
   * returns true if the patch was already occupied */
	if (!Agent::agentPatchPtr[newIndex].setOccupied()) { 
		// Get the current location
		int in = this->index[read_index];
		// Update residing patch to 'available'
		Agent::agentPatchPtr[in].clearOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = unoccupied;
		// Update location coordinates
		this->ix[write_t] = targetX;
		this->iy[write_t] = targetY;
		this->iz[write_t] = targetZ;
		this->index[write_t] = newIndex;
		Agent::agentPatchPtr[newIndex].occupiedby[write_t] = this->type[read_t];
		return true;
	} else {
		//cout << "ERROR error in cell moving" << endl;
		return false;
	}
}

bool Agent::moveToTarget(int dx, int dy, int dz, int index)
{
    int read_index;
    // Check if the location has been modified in this tick
    if (isModified(this->index)) {
        // If it has, work off of the intermediate value
        read_index = write_t;
    } else {
        // If it has NOT, work off of the original value
        read_index = read_t;
    }

    return this->move(dx, dy, dz, read_index);
}

void Agent::wiggle() {                                                      
	int read_index;
	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		// If it has, work off of the intermediate value
		read_index = write_t;
	} else {
		// If it has NOT, work off of the original value
		read_index = read_t;
	}
  // Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];
	int currentindex = this->index[read_index];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	int tid = 0;
#ifdef _OMP
	// Get thread id in order to access the seed that belongs to this thread
	tid = omp_get_thread_num();
#endif
	int trial = 0;
	do {
		// Pick a neighbor to move to at random
#ifndef MODEL_3D
		int i = rand_r(&(agentWorldPtr->seeds[tid])) % 8;
#else
		int i = rand_r(&(agentWorldPtr->seeds[tid])) % 27;
#endif
		int dx = Agent::dX[neighbor[i]];
		int dy = Agent::dY[neighbor[i]];
		int dz = Agent::dZ[neighbor[i]];

		// Calculate target index
		int newindex = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// If trying to move off the side boundaries, die
		if (z + dz < 0 || z + dz >= nz || y + dy < 0 || y + dy >= ny) {
			this->life[write_t] = 0;
			this->die();
			return;
		}

		// Check if target patch is valid
		bool isValid = (x + dx >= 0 || x + dx < nx) &&  // within x-boundaries
				!(Agent::agentPatchPtr[newindex].isOccupiedWrite()) && // not occupied
				(Agent::agentPatchPtr[newindex].isHealthy()); // healthy

		// If the target patch is not valid, pick a new neighbor
		if (!isValid) {
			continue;
		}


    // Patch type constraints on movement: TODO(Kim): INSERT REFS?
		if (Agent::agentPatchPtr[newindex].type[read_t] == tissue) {
			// If target patch is tissue, ok to move
		} else if (Agent::agentPatchPtr[newindex].type[read_t] == epithelium) {
			// If target patch is epithelium, move in the opposite direction
			dx = -dx;
			dy = -dy;
		} else if ((Agent::agentPatchPtr[newindex].type[read_t] == capillary) && (Agent::agentPatchPtr[currentindex].type[read_t] == capillary)) {
			// If target patch is capillary and current patch is capillary, move
		} else if ((Agent::agentPatchPtr[newindex].type[read_t] == capillary) && (Agent::agentPatchPtr[currentindex].type[read_t] != capillary)) {
			// If target patch is capillary but current patch is not
			vector<int> xtarget;
			vector<int> ytarget;
			vector<int> ztarget;

			/* Look for neighbors that are not capillary that are inside world dimensions */
			for (int dxx = -1; dxx <= 1; dxx++) {
				for (int dyy = -1; dyy <= 1; dyy++) {
					for (int dzz = -1; dzz <= 1; dzz++) {
						if (x + dxx < 0 || x + dxx >= nx || y + dyy < 0 || y + dyy >= ny || z + dzz < 0 || z + dzz >= nz) continue;
						int in = (x + dxx) + (y + dyy)*nx + (z + dzz)*nx*ny;
						if(Agent::agentPatchPtr[in].type[read_t] != capillary) {
							xtarget.push_back(dxx);
							ytarget.push_back(dyy);
							ztarget.push_back(dzz);
						}
					}
				}
			}

			/* Move to a random neighbor that is not capillary that is inside world dimensions */
			int randInt = rand_r(&(agentWorldPtr->seeds[tid]))%(xtarget.size());
			dx = xtarget[randInt];
			dy = ytarget[randInt];
			dz = ztarget[randInt];
		} else {
			cout << "exception! encountered by Agent on " << this->index[read_t] << " " << Agent::agentPatchPtr[newindex].type[read_t];
			cout <<  " " << Agent::agentPatchPtr[this->index[read_t]].type[read_t] << endl;
		}

		if (this->move(dx, dy, dz, read_index) == true) {
			// If move() was successful, get out of the while loop
			break;
		} else {
			// If move() was NOT successful, pick a new neighbor
			continue;
		}
	} while(++trial < 27); //TODO(Nuttiiya): Why is this here?

}

void Agent::wiggle_wound(){
	int read_index;
	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		// If it has, work off of the intermediate value
		read_index = write_t;
	} else {
		// If it has NOT, work off of the original value
		read_index = read_t;
	}
  // Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	// Wiggle with <- x bias
	for (int dx = -1; dx <= 1; dx++) {
//	int dx = -1;
		for (int dy = -1; dy <= 1; dy++) {
			for (int dz = -1; dz <= 1; dz++) {
				// Calculate target index
				int newindex = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

				// If trying to move off the side boundaries, die
				if (z + dz < 0 || z + dz >= nz || y + dy < 0 || y + dy >= ny) {
					this->life[write_t] = 0;
					this->die();
					return;
				}

				// Check if target patch is valid
				bool isValid = (x + dx >= 0 || x + dx < nx) &&  // within x-boundaries
						!(Agent::agentPatchPtr[newindex].isOccupiedWrite()) && // not occupied
						(Agent::agentPatchPtr[newindex].isHealthy()); // healthy

				// If the target patch is not valid, pick a new neighbor
				if (!isValid) {
					continue;
				}

				if (this->move(dx, dy, dz, read_index) == true) {
					// If move() was successful, get out of the for loop
					break;
				}

			}
		}
	}


}

float Agent::meanNeighborChem(int chemIndex) {
	int totalchemical = 0, numberofpatches = 0;
  // Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

  /* Count the number of chemicals of type chemIndex in all neighbors that are
   * inside world dimensions */
  for (int dZ = -1; dZ <= 1; dZ++) {
    for (int dY = -1; dY <= 1; dY++) {
      for (int dX = -1; dX <= 1; dX++) {
        if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz)
        { continue; }
          int in = (x + dX) + (y + dY)*nx + (z + dZ)*nx*ny;
#ifdef OPT_CHEM
          totalchemical += Agent::agentWorldPtr->WHWorldChem->getPchem(chemIndex, in);
#else
          totalchemical += Agent::agentWorldPtr->WHWorldChem->pChem[chemIndex][in];
#endif
          numberofpatches++;
      }
    }
  }
	return totalchemical/numberofpatches;
}

float Agent::countNeighborECM(int ECMIndex) {
	float numberofecm = 0.0f;
  // Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

  /* Count the number of ECM proteins of type ECMIndex in all neighbors that
   * are inside world dimensions */
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz) continue;
				int in = (x + dX) + (y + dY)*nx + (z + dZ)*nx*ny;
				switch (ECMIndex) {
          case oc:
            numberofecm += Agent::agentECMPtr[in].monomer_col[read_t];
            break;
          case nc:
            numberofecm += Agent::agentECMPtr[in].polymer_col;
            break;
          case fc:
            numberofecm += Agent::agentECMPtr[in].fragmnt_col[read_t];
            break;
          case oe:
            numberofecm += Agent::agentECMPtr[in].monomer_ela[read_t];
            break;
          case ne:
            numberofecm += Agent::agentECMPtr[in].polymer_ela;
            break;
          case fe:
            numberofecm += Agent::agentECMPtr[in].fragmnt_ela[read_t];
            break;
          case oha:
            numberofecm += Agent::agentECMPtr[in].getnHA();//HA[read_t];
            break;
          case nha:
            numberofecm += Agent::agentECMPtr[in].getnHA();//HA[read_t];
            break;
          case fha:
            numberofecm += Agent::agentECMPtr[in].getnfHA();//fHA[read_t];
            break;
				}
			}
		}
	}
	return numberofecm;
}

int Agent::countNeighborCells(int cellIndex) {
	int numberofcells = 0;
  // Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

  /* Count the number of cells of type cellIndex in all neighbors that are
   * inside world dimensions */
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz) continue;
				int in = (x + dX) + (y + dY)*nx + (z + dZ)*nx*ny;
				if (Agent::agentPatchPtr[in].isOccupied() == false) continue;
				if (Agent::agentPatchPtr[in].occupiedby[read_t] == cellIndex) numberofcells++;
			}
		}
	}
	return numberofcells;
}


bool Agent::moveToHighestChem(int chemIndex) {
    bool isGrad = (chemIndex == NEUgrad) || (chemIndex == MACgrad) || (chemIndex == FIBgrad);
    return this->moveToHighestChem(chemIndex, isGrad);
}

bool Agent::moveToHighestChem(int chemIndex, bool isGrad) {
	int read_index;
	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		// If it has, work off of the intermediate value
		read_index = write_t;
	} else {
		// If it has NOT, work off of the original value
		read_index = read_t;
	}

  // Location of agent in x,y,z dimensions of world.
	int ix = this->ix[read_index];
	int iy = this->iy[read_index];
	int iz = this->iz[read_index];
	int index = this->index[read_index];
  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	double highestchem = Agent::agentWorldPtr->WHWorldChem->getLevel(chemIndex, index, isGrad);
	int dx = 0, dy = 0, dz = 0;

  /* Find the neighbor inside world dimensions with the highest concentration
   * of chemical of type chemIndex */
  for (int dzz = -1; dzz <= 1; dzz++) {
    for (int dyy = -1; dyy <= 1; dyy++) {
      for (int dxx = -1; dxx <= 1; dxx++) {
        if (ix + dxx < 0 || ix + dxx >= nx || iy + dyy < 0 || iy + dyy >= ny || iz + dzz < 0 || iz + dzz >= nz) continue;
        int in = (ix + dxx) + (iy + dyy)*nx + (iz + dzz)*nx*ny;
        float currentChem = Agent::agentWorldPtr->WHWorldChem->getLevel(chemIndex, in, isGrad);
//        if (!Agent::agentPatchPtr[in].isHealthy) continue;
        if (currentChem > highestchem) {
          highestchem = currentChem;
          dx = dxx;
          dy = dyy;
          dz = dzz;
        }
      }
    }
  }

  /* If highestchem is zero, no movement */
  if (fabs(highestchem - 0.0f) < 0.00000000001) return false;
  /* Move to the neighbor with the highest concentration of chemical of type
   * chemIndex if it is not the agent's current patch */
	if (dx == 0 && dy == 0 && dz == 0) return false;
	int newIndex = (ix + dx) + (iy + dy)*nx + (iz + dz)*nx*ny;
	return this->move(dx, dy, dz, read_index);

}

void Agent::updateAgent() {
  //DEBUG
//  if (this->type[read_t] == 1)
//    printf("-- Fib update agent\n");
	this->ix[read_t] = this->ix[write_t];
	this->iy[read_t] = this->iy[write_t];
	this->iz[read_t] = this->iz[write_t];
	this->index[read_t] = this->index[write_t];
	this->alive[read_t] = this->alive[write_t];
	this->life[read_t] =  this->life[write_t];

	//DEBUG
//	if (this->activate[read_t] != this->activate[write_t] && this->activate[write_t] == true)
//	  printf("   activating %d\n", this->type[read_t]);
	this->activate[read_t] = this->activate[write_t];
	this->color[read_t] = this->color[write_t];
	this->size[read_t] = this->size[write_t];
	this->type[read_t] = this->type[write_t];
}

