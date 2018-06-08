/* 
 * File: ECM.cpp
 *
 * File Contents: Contains ECM class
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *
 * Created on Oct 20, 2014, 7:37 AM
 */
/*****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#include "ECM.h"
#include "../World/Usr_World/woundHealingWorld.h"
#include "../enums.h"
#include <iostream>
#include <vector>
#include <cstdlib>                                      
#include <stdio.h>                                     
#include <string.h>                                     
#include <algorithm>

#ifdef PROFILE_ECM
#include <time.h>
#include <sys/time.h>
#endif

using namespace std;

//FIXME: Update max num of ECM on each patch

Patch* ECM::ECMPatchPtr = NULL; 
WHWorld* ECM::ECMWorldPtr = NULL;

int ECM::dx[27] = {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1};
int ECM::dy[27] = {-1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1};
int ECM::dz[27] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
int ECM::d[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};    

// First 9 elements are neighbors and self in the same Z-plane (2D)
int ECM::nid[27] = {9, 10, 11, 12, 13, 14, 15, 16, 17,
		0,  1,  2,  3,  4,  5,  6,  7,  8,
		18, 19, 20, 21, 22, 23, 24, 25, 26};

ECM::ECM() {
}

ECM::ECM(int x, int y, int z, int index) {

#ifdef _OMP
	omp_init_lock(&ECMmutex);
#endif

#ifdef OPT_ECM
	this->HAlife = new std::vector<int>;
#else   //  OPT_ECM
	this->HAlife[write_t] = new std::vector<int>;
	this->HAlife[read_t]  = new std::vector<int>;
#endif

	this->dirty = true;
	this->request_dirty = true;
	this->dirty_from_neighbors = true;
	this->indice[0] = x;
	this->indice[1] = y;
	this->indice[2] = z;
	this->index = index;
	this->empty[read_t] = true;
	this->ocollagen[read_t] = 0;
	this->ncollagen[read_t] = 0;
	this->fcollagen[read_t] = 0;
	this->oelastin[read_t] = 0;
	this->nelastin[read_t] = 0;
	this->felastin[read_t] = 0;
	this->HA[read_t] = 0;
	this->fHA[read_t] = 0;
	this->fcDangerSignal[read_t] = false;
	this->feDangerSignal[read_t] = false;
	this->fHADangerSignal[read_t] = false;
	this->scarIndex[read_t] = false;
	this->empty[write_t] = true;
	this->ocollagen[write_t] = 0;
	this->ncollagen[write_t] = 0;
	this->fcollagen[write_t] = 0;
	this->oelastin[write_t] = 0;
	this->nelastin[write_t] = 0;
	this->felastin[write_t] = 0;
	this->HA[write_t] = 0;
	this->fHA[write_t] = 0;
	this->fcDangerSignal[write_t] = false;
	this->feDangerSignal[write_t] = false;
	this->fHADangerSignal[write_t] = false;
	this->scarIndex[write_t] = false;
#ifdef OPT_FECM_REQUEST
	fcollagenrequest = 0;
	felastinrequest  = 0;
	fHArequest       = 0;
#else
	memset(requestfcollagen, 0, 27*sizeof(int));
	memset(requestfelastin, 0, 27*sizeof(int));
	memset(requestfHA, 0, 27*sizeof(int));
#endif	// OPT_FECM_REQUEST

}

ECM::~ECM() {
}

void ECM::set_dirty() {
#ifdef OPT_ECM
	this->dirty = true;
#endif
}

void ECM::set_request_dirty() {
#ifdef OPT_ECM
	this->request_dirty = true;
#endif
}

void ECM::set_dirty_from_neighbors() {
#ifdef OPT_ECM
	this->dirty_from_neighbors = true;
#endif
}

void ECM::reset_dirty() {
#ifdef OPT_ECM
	this->dirty = false;
#endif
}

void ECM::reset_request_dirty() {
#ifdef OPT_ECM
	this->request_dirty = false;
#endif
}

void ECM::reset_dirty_from_neighbors() {
#ifdef OPT_ECM
	this->dirty_from_neighbors = false;
#endif
}

void ECM::decrement(int n) {
	--n;
}


void ECM::lock() {
#ifdef _OMP
	omp_set_lock(&ECMmutex);
#endif
}

void ECM::unlock() {
#ifdef _OMP
	omp_unset_lock(&ECMmutex);
#endif
}


void ECM::cleanFragmentedECM(){
	if (ECMPatchPtr[this->index].occupiedby[read_t] == amacrophag || ECMPatchPtr[this->index].occupiedby[read_t] == aneutrophil) {
		if (this->fcollagen[read_t] > 0) {
			int wFC = this->fcollagen[write_t];
			int rFC = this->fcollagen[read_t];
			if (wFC == rFC)
				this->fcollagen[write_t] = rFC - 1;
			else
				this->fcollagen[write_t] = wFC - 1;
		}

		if (this->felastin[read_t] > 0) {
			int wFE = this->felastin[write_t];
			int rFE = this->felastin[read_t];
			if (wFE == rFE)
				this->felastin[write_t] = rFE - 1;
			else
				this->felastin[write_t] = wFE - 1;
		}

		if (this->fHA[read_t] > 0) {
			int wFHA = this->fHA[write_t];
			int rFHA = this->fHA[read_t];
			int nHArm = 1; // 6
			if (wFHA == rFHA)
				this->fHA[write_t] = rFHA - nHArm;
			else
				this->fHA[write_t] = wFHA - nHArm;
		}
	}
}


#ifdef OPT_FECM_REQUEST

void ECM::addCollagenFragment()
{
	set_request_dirty();
#pragma omp atomic
	(this->fcollagenrequest)++;
}

void ECM::addElastinFragment()
{
	set_request_dirty();
#pragma omp atomic
	(this->felastinrequest)++;
}

void ECM::addHAFragment()
{
	set_request_dirty();
#pragma omp atomic
	(this->fHArequest)++;
}

void ECM::addHAFragments(int nfr)
{
	set_request_dirty();
#pragma omp atomic
	(this->fHArequest) += nfr;
}

#endif	// OPT_FECM_REQUEST

void ECM::ECMFunction() {



	struct timeval t0, t1;
	/*************************************************************************
	 * DAMAGE REPAIR                                                         *
	 *************************************************************************/
	// New ECM proteins can repair damage on their patch  // TODO(Kim): INSERT REF?

	if (this->ncollagen[read_t] > 0 || this->nelastin[read_t] > 0 || this->HA[read_t] > 0) {
#ifdef PROFILE_ECM
		gettimeofday(&t0, NULL);
#endif
		this->repairDamage();

#ifdef PROFILE_ECM
		gettimeofday(&t1, NULL);
		ECM::ECMWorldPtr->ECMrepairTime += (t1.tv_sec-t0.tv_sec)*1000000 + (t1.tv_usec-t0.tv_usec);
#endif
	}

	/*************************************************************************
	 * DEATH                                                                 *
	 *************************************************************************/
	/* Each hyaluronan life decreases at each time step and they die naturally
	 * Collagen & elastin have 'infinite life' because their half lives are on 
	 * the order of years. */
#ifdef OPT_ECM
	/*
// DEBUG 
      int HAsizet = HAlife->size();
      if (HAsizet > ECM::ECMWorldPtr->maxHAsize) 
      {
        ECM::ECMWorldPtr->maxHAsize = HAsizet;
        //printf("\tNew HA size max: %d\n", HAsize);
      }
	 */
#ifdef PROFILE_ECM
	gettimeofday(&t0, NULL);
#endif
	//	for_each(HAlife->begin(), HAlife->end(), decrement);
	int HAsize = HAlife->size();
	for (int i = 0; i < HAsize; i++) {
		decHAlife(i);
	}
#ifdef PROFILE_ECM
	gettimeofday(&t1, NULL);
	ECM::ECMWorldPtr->HAlifeTime += (t1.tv_sec-t0.tv_sec)*1000000 + (t1.tv_usec-t0.tv_usec);
#endif
#else  // OPT_ECM
	for (int i = 0; i < HAlife[read_t]->size(); i++) {
		decHAlife(i);
	}
#endif  // OPT_ECM

	/*************************************************************************
	 * DANGER SIGNALLING                                                     *
	 *************************************************************************/
	// Fragmented ECM proteins can signal danger one time  //TODO(Kim): INSERT REF?
	if (fcDangerSignal[read_t] == true || feDangerSignal[read_t] == true || fHADangerSignal[read_t] == true) {
#ifdef PROFILE_ECM
		gettimeofday(&t0, NULL);
#endif
		this->set_dirty();
		this->dangerSignal();
		fcDangerSignal[write_t] = false;
		feDangerSignal[write_t] = false;
		fHADangerSignal[write_t] = false;
#ifdef PROFILE_ECM
		gettimeofday(&t1, NULL);
		ECM::ECMWorldPtr->ECMdangerTime += (t1.tv_sec-t0.tv_sec)*1000000 + (t1.tv_usec-t0.tv_usec);
#endif
	}


	/*************************************************************************
	 * FRAGMENTS CLEAN UP                                                    *
	 *************************************************************************/
	this->cleanFragmentedECM();

	/*************************************************************************
	 * SCAR FORMATION                                                        *
	 *************************************************************************/
	// Original collagen can create a scar if above threshold of 100 // TODO(Kim): INSERT REF?
#ifdef PROFILE_ECM
	gettimeofday(&t0, NULL);
#endif
	if (ncollagen[read_t] >= 10) {  // TODO after sensitivity
		this->set_dirty();
		scarIndex[write_t] = true;                                 // FIXME
	}
#ifdef PROFILE_ECM
	gettimeofday(&t1, NULL);
	ECM::ECMWorldPtr->ECMscarTime += (t1.tv_sec-t0.tv_sec)*1000000 + (t1.tv_usec-t0.tv_usec);
#endif
}

void ECM::repairDamage() {
	// Location of neighbor in x,y,z dimensions of the world:
	int tempX, tempY, tempZ, tempIndex;
	int nx, ny, nz;
	Patch* tempPatchPtr;
	// Number of patches in x,y,z dimensions of the world:
	nx = ECM::ECMWorldPtr->nx;
	ny = ECM::ECMWorldPtr->ny;
	nz = ECM::ECMWorldPtr->nz;
	//Repair damage on current and Neighbor Patches
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				// Location of patch to be repaired in x,y,z dimensions of world.
				tempX = this->indice[0] + dX;
				tempY = this->indice[1] + dY;
				tempZ = this->indice[2] + dZ;

				// Try a new patch if this one is outside the world dimensions.
				if (tempX < 0 || tempX >= nx || tempY < 0 || tempY >= ny || tempZ < 0 || tempZ >= nz) continue;

				// Get access to the patch on which this ECM manager resides
				tempIndex = tempX + tempY*nx + tempZ*(nx)*(ny);
				tempPatchPtr = &(ECM::ECMPatchPtr[tempIndex]);

				if (tempPatchPtr->damage[write_t] > 0) {        // NS: Why not read from read_t?
					// Data race allowed, since we're just overwriting values
					tempPatchPtr->dirty = true;
					tempPatchPtr->damage[write_t] = 0;
					tempPatchPtr->type[write_t] = tissue;
					tempPatchPtr->color[write_t] = ctissue;
					tempPatchPtr->health[write_t] = 100;
					//                                        cout << " repair damage at " <<tempIndex << endl;
				}
			}
		}
	}
}

void ECM::dangerSignal() {
	this->ECMPatchPtr[this->index].dirty = true;
	this->ECMPatchPtr[this->index].damage[write_t]++;
	this->ECMPatchPtr[this->index].health[write_t] = 0;
	this->ECMPatchPtr[this->index].color[write_t] = cdamage;

	/* Activated macrophages and activated neutrophils can remove newly created
	 * damage if they are present. */  // TODO(Kim): INSERT REF?
	if (ECMPatchPtr[this->index].isOccupied() == false) return;
	if (ECMPatchPtr[this->index].occupiedby[read_t] == amacrophag || ECMPatchPtr[this->index].occupiedby[read_t] == aneutrophil) {
		ECMPatchPtr[this->index].damage[write_t] = 0;
		ECMPatchPtr[this->index].type[write_t] = tissue;
		ECMPatchPtr[this->index].color[write_t] = ctissue;
		ECMPatchPtr[this->index].health[write_t] = 100;
	}
}

void ECM::fragmentNCollagen() {
	// Distance to neighbor in x,y,z dimensions of the world:
	int dX, dY, dZ;
	int newfragments = 0;
	int dn[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = indice[0];
	int iy = indice[1];
	int iz = indice[2];
	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

	// Alert change in status of original collagen
	if (this->ncollagen[write_t] > 0) {
		this->set_dirty();
		//		this->set_request_dirty();
	}

	/* Request hatching two fragmented collagens on neighboring patches for each 
	 * original collagen */
	while (this->ncollagen[write_t] > 0) {
		this->ncollagen[write_t]--;
		newfragments = 0;

		// Request one fragmented collagen at a random inbounds neighboring patch.


		// TODO(Caroline) Might want to make the radius = 2
		std::random_shuffle(&dn[0], &dn[27]);
		for (int i = 0; i < 27; i++) {
			dX = dx[dn[i]];
			dY = dy[dn[i]];
			dZ = dz[dn[i]];
			if (newfragments >= 2) break;
			if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue;
			int in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
#ifdef OPT_FECM_REQUEST
			this->ECMWorldPtr->worldECM[in].addCollagenFragment();
#else	// OPT_FECM_REQUEST
			//'dX + dY*3 + dZ*3*3 + 13' determines which neighbor
			this->requestfcollagen[dX + dY*3 + dZ*3*3 + 13]++;
			// Alert change in status of collagen on this patch
			this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
#endif	// OPT_FECM_REQUEST
			newfragments++;
		}
	}
//#ifdef VISUALIZATION
	// Update ECM polymer map
	this->ECMWorldPtr->setECM(this->index, m_col, 0);
//#endif	// VISUALIZATION
	this->isEmpty();
}

void ECM::fragmentNElastin() {
	// Distance to neighbor in x,y,z dimensions of the world:
	int dX, dY, dZ;
	int dn[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
	int newfragments = 0;
	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = indice[0];
	int iy = indice[1];
	int iz = indice[2];
	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

	// Alert change in status of original elastin
	if (this->nelastin[write_t] > 0) {
		this->set_dirty();
		//		this->set_request_dirty();
	}

	/* Request hatching two fragmented elastins on neighboring patches for each
	 * original elastin */
	while (this->nelastin[write_t] > 0) {
		this->nelastin[write_t]--;
		newfragments = 0;

		// Request one fragmented elastin at a random inbounds neighboring patch.


		// TODO(Caroline) Might want to make the radius = 2
		std::random_shuffle(&dn[0], &dn[27]);
		for (int i = 0; i < 27; i++) {
			dX = dx[dn[i]];
			dY = dy[dn[i]];
			dZ = dz[dn[i]];
			if (newfragments >= 2) break;
			if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue;
			int in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
#ifdef OPT_FECM_REQUEST
			this->ECMWorldPtr->worldECM[in].addElastinFragment();
#else   // OPT_FECM_REQUEST
			//'dX + dY*3 + dZ*3*3 + 13' determines which neighbor
			this->requestfelastin[dX + dY*3 + dZ*3*3 + 13]++;
			// Alert change in status of elastin on this patch
			this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
#endif	// OPT_FECM_REQUEST
			newfragments++;
		}
	}
//#ifdef VISUALIZATION
	// Update ECM polymer map
	this->ECMWorldPtr->setECM(this->index, m_ela, 0);
//#endif	// VISUALIZATION

	this->isEmpty();
}

void ECM::fragmentHA() {
	int tid = 0;
#ifdef _OMP
	// Get thread id in order to access the seed that belongs to this thread
	tid = omp_get_thread_num();
#endif
	// Distance to neighbor in x,y,z dimensions of the world:
	int dX, dY, dZ;
	int dn[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
	int newfragments = 0;
	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = indice[0];
	int iy = indice[1];
	int iz = indice[2];
	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

	// Alert change in status of hyaluronan
	if (HA[write_t] > 0) {
		this->set_dirty();
		//		this->set_request_dirty();
	}

	/* Request hatching two fragmented hyaluronans on neighboring patches for 
	 * each hyaluronan */


#ifndef OPT_ECM
	vector<int> *HAlife = this->HAlife[write_t];
#endif
	int nHA = HAlife->size();

	// Kills all alive hyaluronans and fragment
	for (int iHA = 0; iHA < nHA; iHA++)
	{

		newfragments = 0;

		if (HAlife->at(iHA) == 0) continue;
#ifdef OPT_ECM
		rmvHA(iHA);
#else
		HAlife[iHA] = 0;
#endif

		// Request one fragmented hyaluronan at a random inbounds neighboring patch

		// TODO(Caroline) Might want to make the radius = 2
		std::random_shuffle(&dn[0], &dn[27]);
		for (int i = 0; i < 27; i++) {
			dX = dx[dn[i]];
			dY = dy[dn[i]];
			dZ = dz[dn[i]];
			if (newfragments >= 2) break;
			if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue;
			int in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
#ifdef OPT_FECM_REQUEST
			this->ECMWorldPtr->worldECM[in].addHAFragment();
#else   // OPT_FECM_REQUEST
			//'dX + dY*3 + dZ*3*3 + 13' determines which neighbor in radius 1
			this->requestfHA[dX + dY*3 + dZ*3*3 + 13]++;
			// Alert change in status of hyaluronan on this patch
			this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
#endif	// OPT_FECM_REQUEST
			newfragments++;
		}

	}

//#ifdef VISUALIZATION
			// Update ECM polymer map
			this->ECMWorldPtr->setECM(this->index, m_hya, 0);
//#endif	// VISUALIZATION



	this->isEmpty();
}

void ECM::updateECM() {
	// Patch row major index of neighbor:
	int in;
#ifndef OPT_FECM_REQUEST
	// Amount of requested fragmented ECM proteins:
	int fcollagenrequest = 0, felastinrequest = 0, fHArequest = 0;
#endif
	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = this->indice[0];
	int iy = this->indice[1];
	int iz = this->indice[2];
	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

	/*************************************************************************
	 * FRAGMENTED ECM REQUESTS                                               *
	 *************************************************************************/
	// Iterate through neighboring patches, count any fcollage/elastin/HA requests for ECM manager

#ifdef OPT_FECM_REQUEST
	if (this->request_dirty)
	{
		if (ocollagen[read_t] + ncollagen[read_t] + fcollagen[read_t] + fcollagenrequest > MAX_COL) {
			cout << "  Error fcollagen request" << endl;
		} else if (oelastin[read_t] + nelastin[read_t] + felastin[read_t] + felastinrequest > MAX_ELA){
			cout << "  Error felastin request" << endl;
		} else if (HA[read_t] + fHA[read_t] + fHArequest > MAX_HYA) {
			cout << "  Error fcollagen request" << endl;
		} else {
			// Fragmented ECM proteins serve as danger signals once
			this->fcollagen[write_t] += fcollagenrequest;
			this->fcDangerSignal[write_t] += fcollagenrequest;
			this->felastin[write_t] += felastinrequest;
			this->feDangerSignal[write_t] += felastinrequest;
			this->fHA[write_t] += fHArequest;
			this->fHADangerSignal[write_t] += fHArequest;
			this->isEmpty();
		}
		// NS: Check semantics
		this->fcollagen[read_t] = this->fcollagen[write_t];
		this->felastin[read_t] = this->felastin[write_t];
		this->fHA[read_t] = this->fHA[write_t];

		this->fcollagenrequest = 0;
		this->felastinrequest  = 0;
		this->fHArequest       = 0;
	}
#else	// OPT_FECM_REQUEST

#ifdef OPT_ECM
	// Only process requests if a neighbor has indicated that it's made a fragment request to this ECM manager
	if (this->dirty_from_neighbors) {
#endif	// OPT_ECM
#ifdef ECM_UNROLL_LOOP
		// Distance to neighbors in x,y,z dimensions of the world:
		int dX[8] = {-1, 0, 1, -1, 1, -1, 0, 1};
		int dY[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
		int dZ = 0;  // 2D
		int targetZ = iz + dZ;
		for (int i = 0; i < 8; i++) {
			int targetX = ix + dX[i];
			int targetY = iy + dY[i];
			// Only consider requests from neighbors that are inside the world.
			if (!(targetX < 0 || targetX >= nx || targetY < 0 || targetY >= ny || targetZ < 0 || targetZ >= nz)) {
				in = targetX + targetY*nx + targetZ*nx*ny;
				ECM* neighborECMPtr = &(ECMWorldPtr->worldECM[in]);
				/* self_neighbor_in is the index of this ECM manager in the list of
				 * its neighbor's list of neighbors. */
				int self_neighbor_in = (-dX[i]) + (-dY[i])*3 + (-dZ)*3*3 + 13;
				fcollagenrequest += neighborECMPtr->requestfcollagen[self_neighbor_in];
				felastinrequest += neighborECMPtr->requestfelastin[self_neighbor_in];
				fHArequest += neighborECMPtr->requestfHA[self_neighbor_in];
			}
		}
#else
		// dX, dY, dZ are distance to neighbors in x,y,z dimensions of the world
		for (int dX = -1; dX < 2; dX++) {
			for (int dY = -1; dY < 2; dY++) {
				for (int dZ = -1; dZ < 2; dZ++) {
					// Try another neighbor if this one is out of bounds
					if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue;
					in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
					int a = ECMWorldPtr->worldECM[in].requestfcollagen[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
					int b = ECMWorldPtr->worldECM[in].requestfelastin[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
					int c = ECMWorldPtr->worldECM[in].requestfHA[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
					if (a != 0) {
						//cout << "  fcollagen requested at " << ix + iy*nx + iz*nx*ny << " by " << in << endl;
					}
					if (b != 0) {
						//cout << "  felastin requested at " << ix + iy*nx + iz*nx*ny << " by " << in << endl;
					}
					if (c != 0) {
						//cout << "  fHA requested at " << ix + iy*nx + iz*nx*ny << " by " << in << endl;
					}
					fcollagenrequest += a;
					felastinrequest += b;
					fHArequest += c;
				}
			}
		}
#endif

		// Fragmented ECM requests can only be accepted if there is enough space. // TODO(Kim): INSERT REF?
		if (ocollagen[read_t] + ncollagen[read_t] + fcollagen[read_t] + fcollagenrequest > MAX_COL) {
			cout << "  Error fcollagen request" << endl;
		} else if (oelastin[read_t] + nelastin[read_t] + felastin[read_t] + felastinrequest > MAX_ELA) {
			cout << "  Error felastin request" << endl;
		} else if (HA[read_t] + fHA[read_t] + fHArequest > MAX_HYA) {
			cout << "  Error fcollagen request" << endl;
		} else {
			// Fragmented ECM proteins serve as danger signals once
			this->fcollagen[write_t] += fcollagenrequest;
			this->fcDangerSignal[write_t] += fcollagenrequest;
			this->felastin[write_t] += felastinrequest;
			this->feDangerSignal[write_t] += felastinrequest;
			this->fHA[write_t] += fHArequest;
			this->fHADangerSignal[write_t] += fHArequest;
			this->isEmpty();
		}
		this->fcollagen[read_t] = this->fcollagen[write_t];
		this->felastin[read_t] = this->felastin[write_t];
		this->fHA[read_t] = this->fHA[write_t];

#ifdef OPT_ECM
	}	// if (this->dirty_from_neighbors)
#endif

#endif	// OPT_FECM_REQUEST
	/*************************************************************************
	 * READ/WRITE SYNCHRONIZATION                                            *
	 *************************************************************************/
#ifdef OPT_ECM
	// Only synchronize the read and write entries if there's a change in value
	if (this->dirty) {
#endif
		/*
		this->empty[read_t] = this->empty[write_t];
		this->ocollagen[read_t] = this->ocollagen[write_t];
		this->ncollagen[read_t] = this->ncollagen[write_t];
		this->oelastin[read_t] = this->oelastin[write_t];
		this->nelastin[read_t] = this->nelastin[write_t];
		this->HA[read_t] = this->HA[write_t];
		 */
		// Convert ocollagen (tropocollagen monomer) to ncollagen (polymer)
		int o_nCollRatio = 2;
		int o_nElasRatio = 2;

		int ocoll = this->ocollagen[read_t];
		int oelas = this->oelastin[read_t];

		int leftoverOcoll = ocoll%o_nCollRatio;
		int leftoverOelas = oelas%o_nElasRatio;

		// Turn o<ECM> to n<ECM> (monomer to polymer)
		this->ncollagen[write_t] += ocoll/o_nCollRatio;
		this->ncollagen[read_t]  =  this->ncollagen[write_t];
		this->nelastin[write_t]  += oelas/o_nElasRatio;
		this->nelastin[read_t]   =  this->nelastin[write_t];
//#ifdef VISUALIZATION
		// Update ECM polymer map
		// DEBUG vis
		if (this->indice[2] == 14)
			this->ECMWorldPtr->incECM(this->index, m_col, (ocoll/o_nCollRatio));
		this->ECMWorldPtr->incECM(this->index, m_ela, (oelas/o_nElasRatio));
//#endif	// VISUALIZATION

		// Combine leftover monomer with newly deposited monomer
		this->ocollagen[write_t] += leftoverOcoll;
		this->oelastin[write_t]  += leftoverOelas;

		// Synchronize
		this->empty[read_t] = this->empty[write_t];
		this->ocollagen[read_t] = this->ocollagen[write_t];
		this->ncollagen[read_t] = this->ncollagen[write_t];
		this->oelastin[read_t] = this->oelastin[write_t];
		this->nelastin[read_t] = this->nelastin[write_t];
		this->HA[read_t] = this->HA[write_t];

		// Reset o<ECM>[write_t]
		                // Note: [write_t] field is use as buffer for amount added earlier in this tick
		//       It does NOT accumulate values from ealier ticks
		this->ocollagen[write_t] = 0;
		this->oelastin[write_t]  = 0;

		// Get rid of all dead hyaluronans
#ifdef OPT_ECM
		vector<int> *vec = this->HAlife;
#else
		vector<int> *vec = this->HAlife[write_t];
#endif
		vector<int>::iterator first	= vec->begin();
		for (vector<int>::iterator it = first; it != vec->end();) {
			int life = *it;
			if (life <= 0) {
				it = vec->erase(it);
			} else {
				++it;
			}
		}

#ifndef OPT_ECM
		this->HAlife->at(read_t) = this->HAlife->at(write_t);
#endif
		this->fcDangerSignal[read_t] = this->fcDangerSignal[write_t];
		this->feDangerSignal[read_t] = this->feDangerSignal[write_t];
		this->fHADangerSignal[read_t] = this->fHADangerSignal[write_t];
		this->scarIndex[read_t] = this->scarIndex[write_t];
#ifdef OPT_ECM
	}	// if (this->dirty)
#endif

	// If number of oECM exceed threshold, replace with nECM 
	//  define ocollagen to ncollagen (collagen molecule) conversion threshold, after sensitivity analysis


	// Remove all dirty flags
	this->reset_dirty();
#ifndef OPT_FECM_REQUEST
	this->reset_dirty_from_neighbors();
#endif	// ! OPT_FECM_REQUEST



}

void ECM::isEmpty() {
	int totalcollagen = ocollagen[write_t] + ncollagen[write_t] + fcollagen[write_t] ;
	int totalelastin = oelastin[write_t] + nelastin[write_t] + felastin[write_t];
	int totalHA = HA[write_t] + fHA[write_t];
	if (totalcollagen + totalelastin + totalHA == 0) {
		this->empty[write_t] = true;
	} else {
		this->empty[write_t] = false;
	}
}


void ECM::addHAs(int nHA)
{

	this->lock();		// prevent data race

	this->HA[write_t] += nHA;
	for (int i = 0; i < nHA; i++)
	{
#ifdef OPT_ECM
		this->HAlife->push_back(100);
#else
		this->HAlife[write_t]->push_back(100);
		this->HAlife[read_t]->push_back(100);
#endif
	}

//#ifdef VISUALIZATION
	// Update ECM polymer map
	this->ECMWorldPtr->incECM(this->index, m_hya, nHA);
//#endif	// VISUALIZATION

	this->unlock();	// prevent data race

#ifdef OPT_ECM
	this->set_dirty();
#endif

#ifdef DEBUG_HA_SYNC
	if (!checknHAbool()) {
		printf("\tadd HA\n");
		exit(-1);
	}
#endif	// DEBUG_HA_SYNC
}

bool ECM::isModifiedHA()
{
	return (this->HA[read_t]) != (this->HA[write_t]);
}

int ECM::getnHA()
{
	return this->HA[read_t];
}

int ECM::getnfHA()
{
	return this->fHA[read_t];
}

int ECM::getnHAlife()
{
	int nHA = 0;
	for (std::vector<int>::iterator it = HAlife->begin(); it != HAlife->end(); ++it)
	{
		if (*it != 0) nHA++;
	}
	return nHA;
}

bool ECM::rmvHA(int iHA)
{
#ifdef OPT_ECM
	if (this->HAlife->size() <= iHA) return false;
	if (this->HAlife->at(iHA) == 0)     return false;
	this->HAlife->at(iHA) = 0;
	this->set_dirty();
#else
	if (this->HAlife[write_t]->size() <= iHA) return false;
	if (this->HAlife[write_t].at(iHA) == 0)     return false;
	this->HAlife[write_t]->at(iHA) = 0;
#endif
	this->HA[write_t]--;

#ifdef DEBUG_HA_SYNC
	if (!checknHAbool()) {
		printf("\tremove HA\n");
		exit(-1);
	}
#endif	// DEBUG_HA_SYNC

	return true;
}

void ECM::checkHAdiff()
{
	int diff = this->HA[write_t] - this->HA[read_t];
	if ((diff != 1) & (diff != 0)) {printf("DIFF: %d\n", diff); exit(-1);}
}

void ECM::checknHA()
{
	int nHA1 = this->HA[write_t];
	int nHA2 = this->getnHAlife();
	if (nHA1 != nHA2)
	{
		printf("ECM[%d, %d, %d] -- nHA1: %d\tnHA2: %d\n", indice[0], indice[1], indice[2],
				nHA1, nHA2);
		exit(-1);
	}
}

bool ECM::checknHAbool()
{
	int nHA1 = this->HA[write_t];
	int nHA2 = this->getnHAlife();
	if (nHA1 != nHA2)
	{
		printf("ECM[%d, %d, %d] -- nHA1: %d\tnHA2: %d\n", indice[0], indice[1], indice[2],
				nHA1, nHA2);
		return false;
	}
	return true;
}

bool ECM::killHA(int iHA)
{
#ifdef OPT_ECM
	if (this->HAlife->at(iHA) <= 0) return false;
	this->HAlife->at(iHA) = 0;
#else
	if (this->HAlife[write_t]->at(iHA) <= 0) return false;
	this->HAlife[write_t]->at(iHA);
#endif
	this->HA[write_t]--;
	return true;
}

bool ECM::killHA(std::vector<int>::iterator itHA)
{
	if (*itHA <= 0) return false;
	*itHA = 0;
	this->HA[write_t]--;
	return true;
}

bool ECM::decHAlife(int iHA)
{

#ifdef DEBUG_HA_SYNC
	if (!checknHAbool()) {
		printf("\t(1) dec HA\n");
		exit(-1);
	}
#endif	// DEBUG_HA_SYNC

	int oldlife = this->HAlife->at(iHA);
	int oldHA = this->HA[write_t];

#ifdef OPT_ECM
	if (this->HAlife->at(iHA) <= 0) return false;
	int life = --(this->HAlife->at(iHA));
#else
	if (this->HAlife[write_t]->at(iHA) <= 0) return false;
	int life = --(this->HAlife[write_t]->at(iHA));
#endif
	if (life == 0) this->HA[write_t]--;

	return true;

}

void ECM::addColls(int nColls)
{

	this->set_dirty();
#pragma omp atomic
	this->ocollagen[write_t] += nColls;

}

int ECM::getnColl()
{
	return this->ocollagen[read_t];
}

void ECM::resetrequests() {
#ifdef OPT_ECM
	// Only clear the requests if there are any in this tick
	if (this->request_dirty) {
#endif
#ifdef OPT_FECM_REQUEST
		// Already reset in updateECM()
#else	// OPT_FECM_REQUEST
		memset(this->requestfcollagen, 0, 27*sizeof(int));
		memset(this->requestfelastin, 0, 27*sizeof(int));
		memset(this->requestfHA, 0, 27*sizeof(int));
#endif	// OPT_FECM_REQUEST
		this->reset_request_dirty();
#ifdef OPT_ECM
	}
#endif
	//	this->reset_request_dirty();
}
