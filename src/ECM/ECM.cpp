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

	//#ifdef _OMP
	omp_init_lock(&ECMlock);
	//#endif

	this->dirty = true;
	this->request_dirty = true;
	this->dirty_from_neighbors = true;
	this->indice[0] = x;
	this->indice[1] = y;
	this->indice[2] = z;
	this->index = index;
	this->polymer_col = 0;
	this->monomer_col[read_t] = 0;
	this->fragmnt_col[read_t] = 0;
	this->polymer_ela = 0;
	this->monomer_ela[read_t] = 0;
	this->fragmnt_ela[read_t] = 0;
	this->totalHA[read_t] = 0.0f;
	this->fHA[read_t] = 0;
	this->fcDangerSignal[read_t] = false;
	this->feDangerSignal[read_t] = false;
	this->fHADangerSignal[read_t] = false;
	this->scarIndex = false;
	this->polymer_col = 0;
	this->monomer_col[write_t] = 0;
	this->fragmnt_col[write_t] = 0;
	this->polymer_ela = 0;
	this->monomer_ela[write_t] = 0;
	this->fragmnt_ela[write_t] = 0;
	this->totalHA[write_t] = 0.0f;
	this->fHA[write_t] = 0;
	this->fcDangerSignal[write_t] = false;
	this->feDangerSignal[write_t] = false;
	this->fHADangerSignal[write_t] = false;
	this->scarIndex = false;
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
#ifdef _OMP
	omp_destroy_lock(&ECMlock);
#endif
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


inline void ECM::lock() {
#ifdef _OMP
	omp_set_lock(&ECMlock);
#endif
}

inline void ECM::unlock() {
#ifdef _OMP
	omp_unset_lock(&ECMlock);
#endif
}


void ECM::cleanFragmentedECM(){
	if (ECMPatchPtr[this->index].occupiedby[read_t] == amacrophag || ECMPatchPtr[this->index].occupiedby[read_t] == aneutrophil) {
		if (this->fragmnt_col[read_t] > 0) {
			int wFC = this->fragmnt_col[write_t];
			int rFC = this->fragmnt_col[read_t];
			if (wFC == rFC)
				this->fragmnt_col[write_t] = rFC - 1;
			else
				this->fragmnt_col[write_t] = wFC - 1;
		}

		if (this->fragmnt_ela[read_t] > 0) {
			int wFE = this->fragmnt_ela[write_t];
			int rFE = this->fragmnt_ela[read_t];
			if (wFE == rFE)
				this->fragmnt_ela[write_t] = rFE - 1;
			else
				this->fragmnt_ela[write_t] = wFE - 1;
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

	int repairAmount = (int) (this->polymer_col + this->polymer_ela + this->totalHA[read_t]);
#ifdef PROFILE_ECM
	gettimeofday(&t0, NULL);
#endif
	if (repairAmount) this->repairDamage(repairAmount);

#ifdef PROFILE_ECM
	gettimeofday(&t1, NULL);
	ECM::ECMWorldPtr->ECMrepairTime += (t1.tv_sec-t0.tv_sec)*1000000 + (t1.tv_usec-t0.tv_usec);
#endif



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
	if (monomer_col[read_t] >= 10) {  // TODO after sensitivity
		this->set_dirty();
		scarIndex = true;                                 // FIXME
	}
#ifdef PROFILE_ECM
	gettimeofday(&t1, NULL);
	ECM::ECMWorldPtr->ECMscarTime += (t1.tv_sec-t0.tv_sec)*1000000 + (t1.tv_usec-t0.tv_usec);
#endif
}

void ECM::repairDamage(int repairAmount) {

	int in = this->index;
	Patch *tempPatchPtr = &(ECM::ECMPatchPtr[in]);
	// if no damage on patch, do nothing
	if (!tempPatchPtr->damage[read_t]) return;

	if (tempPatchPtr->damage[read_t] > repairAmount) {
		// Data race allowed, since we're just overwriting values
		tempPatchPtr->damage[write_t] -= repairAmount;
	} else {
		tempPatchPtr->damage[write_t] = 0;
		tempPatchPtr->color[write_t] = ctissue;
	}
	tempPatchPtr->dirty = true;
}

void ECM::dangerSignal() {
	this->ECMPatchPtr[this->index].dirty = true;
	this->ECMPatchPtr[this->index].damage[write_t]++;
	this->ECMPatchPtr[this->index].color[write_t] = cdamage;

	/* Activated macrophages and activated neutrophils can remove newly created
	 * damage if they are present. */  // TODO(Kim): INSERT REF?
	if (ECMPatchPtr[this->index].isOccupied() == false) return;
	if (ECMPatchPtr[this->index].occupiedby[read_t] == amacrophag || ECMPatchPtr[this->index].occupiedby[read_t] == aneutrophil) {
		ECMPatchPtr[this->index].damage[write_t]--;
		ECMPatchPtr[this->index].color[write_t] = ctissue;
	}
}

void ECM::fragmentCol() {
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

	int cols = (int) this->polymer_col;

	// Alert change in status of original collagen
	if (this->polymer_col >= 1.0f) {
		this->set_dirty();


		/* Request hatching two fragmented collagens on neighboring patches for each
		 * original collagen */
		while (this->polymer_col >= 1.0f) {
			this->polymer_col--;
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
#else   // OPT_FECM_REQUEST
				//'dX + dY*3 + dZ*3*3 + 13' determines which neighbor
				this->requestfragmnt_col[dX + dY*3 + dZ*3*3 + 13]++;
				// Alert change in status of collagen on this patch
				this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
#endif  // OPT_FECM_REQUEST
				newfragments++;
			}
		}

		this->polymer_col = 0.0f;

		//#ifdef VISUALIZATION
		// Update ECM polymer map
		this->ECMWorldPtr->resetECM(this->index, m_col);
		//#endif    // VISUALIZATION
	}



}

void ECM::fragmentEla() {
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
	if (this->polymer_ela >= 1.0f) {
		this->set_dirty();
		//		this->set_request_dirty();


		/* Request hatching two fragmented elastins on neighboring patches for each
		 * original elastin */
		while (this->polymer_ela >= 1.0f) {
			this->polymer_ela--;
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
				this->requestfragmnt_ela[dX + dY*3 + dZ*3*3 + 13]++;
				// Alert change in status of elastin on this patch
				this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
#endif	// OPT_FECM_REQUEST
				newfragments++;
			}
		}

		this->polymer_ela = 0.0f;
		//#ifdef VISUALIZATION
		// Update ECM polymer map
		this->ECMWorldPtr->resetECM(this->index, m_ela);
		//#endif	// VISUALIZATION
	}

}

void ECM::fragmentHA(bool isInit) {

	if (isInit) {
		// each molecule becomes 1 fragment
		this->fHA[write_t] += this->totalHA[write_t];	// during initialization
	} else {
		// each molecule becomes 1 fragment
		this->fHA[write_t] += this->totalHA[read_t];
	}

	// remove all original HA molecules
	this->rmvAllHAs();


	//#ifdef VISUALIZATION
	// Update ECM polymer map
	this->ECMWorldPtr->resetECM(this->index, m_hya);
	//#endif	// VISUALIZATION


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
		//		if (polymer_col + monomer_col[read_t] + fragmnt_col[read_t] + fcollagenrequest > MAX_COL) {
		//			cout << "  Error fcollagen request" << endl;
		//		} else if (polymer_ela + monomer_ela[read_t] + fragmnt_ela[read_t] + felastinrequest > MAX_ELA){
		//			cout << "  Error felastin request" << endl;
		//		} else if (totalHA[read_t] + fHA[read_t] + fHArequest > MAX_HYA) {
		//			cout << "  Error fcollagen request" << endl;
		//		} else {
		// Fragmented ECM proteins serve as danger signals once
		this->fragmnt_col[write_t] += fcollagenrequest;
		this->fcDangerSignal[write_t] += fcollagenrequest;
		this->fragmnt_ela[write_t] += felastinrequest;
		this->feDangerSignal[write_t] += felastinrequest;
		this->fHA[write_t] += fHArequest;
		this->fHADangerSignal[write_t] += fHArequest;
		//		}
		// NS: Check semantics
		this->fragmnt_col[read_t] = this->fragmnt_col[write_t];
		this->fragmnt_ela[read_t] = this->fragmnt_ela[write_t];
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
				fcollagenrequest += neighborECMPtr->requestfragmnt_col[self_neighbor_in];
				felastinrequest += neighborECMPtr->requestfragmnt_ela[self_neighbor_in];
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
					int a = ECMWorldPtr->worldECM[in].requestfragmnt_col[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
					int b = ECMWorldPtr->worldECM[in].requestfragmnt_ela[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
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
		//		if (polymer_col + monomer_col[read_t] + fragmnt_col[read_t] + fcollagenrequest > MAX_COL) {
		//			cout << "  Error fcollagen request" << endl;
		//		} else if (polymer_ela + monomer_ela[read_t] + fragmnt_ela[read_t] + felastinrequest > MAX_ELA) {
		//			cout << "  Error felastin request" << endl;
		//		} else if (HA[read_t] + fHA[read_t] + fHArequest > MAX_HYA) {
		//			cout << "  Error fcollagen request" << endl;
		//		} else {
		// Fragmented ECM proteins serve as danger signals once
		this->fragmnt_col[write_t] += fcollagenrequest;
		this->fcDangerSignal[write_t] += fcollagenrequest;
		this->fragmnt_ela[write_t] += felastinrequest;
		this->feDangerSignal[write_t] += felastinrequest;
		this->fHA[write_t] += fHArequest;
		this->fHADangerSignal[write_t] += fHArequest;
		//		}
		this->fragmnt_col[read_t] = this->fragmnt_col[write_t];
		this->fragmnt_ela[read_t] = this->fragmnt_ela[write_t];
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
		this->monomer_col[read_t] = this->monomer_col[write_t];
		this->monomer_ela[read_t] = this->monomer_ela[write_t];
		this->HA[read_t] = this->HA[write_t];
		 */
		// Convert polymer_col (troppolymer_col monomer) to monomer_col (polymer)
		int mp_CollRatio = CONV_RATE_COL;//2;
		int mp_ElasRatio = CONV_RATE_ELA;//2;

		int mcoll = this->monomer_col[read_t];
		int melas = this->monomer_ela[read_t];

		float npc_floor;
		float npe_floor;
		float nha_new = this->totalHA[write_t] - this->totalHA[read_t];

		//		if (nha_new > 0.0f) this->ECMWorldPtr->incECM(this->index, m_hya, nha_new);
		// Added in addHA

		if (mcoll >= mp_CollRatio) {      // enough monomers to convert
			float np = mcoll/mp_CollRatio;
			float npc_floor = std::floor(np);
			this->polymer_col += npc_floor;
			this->monomer_col[write_t] = np-npc_floor;     // monomer after conversion
			this->ECMWorldPtr->incECM(this->index, m_col, npc_floor);
		}

		if (melas >= mp_ElasRatio) {      // enough monomers to convert
			float np = melas/mp_ElasRatio;
			float npe_floor = std::floor(np);
			this->polymer_ela += npe_floor;
			this->monomer_ela[write_t] = np-npe_floor;     // monomer after conversion
			this->ECMWorldPtr->incECM(this->index, m_ela, npe_floor);
		}

		// Synchronize
		this->monomer_col[read_t] = this->monomer_col[write_t];
		this->monomer_ela[read_t] = this->monomer_ela[write_t];
		this->totalHA[read_t] = this->totalHA[write_t];

		//#ifdef VISUALIZATION
		// Update ECM polymer map
		// DEBUG vis
		//if (this->indice[2] == 14)




		//#endif        // VISUALIZATION

		// Get rid of all dead hyaluronans
		this->rmvDeadHAs();

		this->fcDangerSignal[read_t] = this->fcDangerSignal[write_t];
		this->feDangerSignal[read_t] = this->feDangerSignal[write_t];
		this->fHADangerSignal[read_t] = this->fHADangerSignal[write_t];
#ifdef OPT_ECM
	}	// if (this->dirty)
#endif


	// Remove all dirty flags
	this->reset_dirty();
#ifndef OPT_FECM_REQUEST
	this->reset_dirty_from_neighbors();
#endif	// ! OPT_FECM_REQUEST

}

void ECM::addMonomerEla(float nEla)
{
#pragma omp atomic
	this->monomer_ela[write_t] += nEla;

#ifdef OPT_ECM
	this->set_dirty();
#endif
}

void ECM::addPolymerEla(float nEla)
{
#pragma omp atomic
	this->polymer_ela += nEla;

#ifdef OPT_ECM
	this->set_dirty();
#endif
}

void ECM::addHAs(int expTick, float nHA)
{
	// nHA already converted to relative unit in Fibroblast::make_hyaluronans()

	float currHA = 0.0f;
	int tick = expTick;
	// default life = 100 ticks
	if (expTick == -1) tick = this->ECMWorldPtr->getTick() + 100;


	//	this->lock(); 		// prevent data race
#ifdef _OMP
#pragma omp critical
#endif
	{
		// find other HAs with the same expiration tick
		ITERATOR_HL currHA_it = HAlife.find(tick);

		if (currHA_it != HAlife.end())
		{
			// get current amount a key 'tick'
			currHA = HAlife.at(tick);
			// add new amount to original amount
			currHA += nHA;

			// add total amount to the map
			currHA_it->second = currHA;

		}

		// update total HA with NEW amount
		this->totalHA[write_t] += nHA;

		//#ifdef VISUALIZATION
		// Update ECM polymer map
		this->ECMWorldPtr->incECM(this->index, m_hya, nHA);
		//#endif	// VISUALIZATION

	}
	//	this->unlock();	// prevent data race

#ifdef OPT_ECM
	this->set_dirty();
#endif

}

bool ECM::isModifiedHA()
{
	return (this->totalHA[read_t]) != (this->totalHA[write_t]);
}

float ECM::getnHA()
{
	return this->totalHA[read_t];
}

int ECM::getnfHA()
{
	return this->fHA[read_t];
}


void ECM::rmvHA(int tick)
{

	float removedHA = 0.0f;

	// find the amount of HAs with the specified expiration tick
	ITERATOR_HL rmvHA_it = HAlife.find(tick);
	if (rmvHA_it != HAlife.end()) removedHA = HAlife.at(tick);
	else return;


	// remove total amount to the map
	this->HAlife.erase(rmvHA_it);
	// update total HA with the amount removed
	this->totalHA[write_t] -= removedHA;


	//#ifdef VISUALIZATION
	// Update ECM polymer map
	this->ECMWorldPtr->decECM(this->index, m_hya, removedHA/MAX_HYA);
	//#endif    // VISUALIZATION


#ifdef OPT_ECM
	this->set_dirty();
#endif

}

void ECM::rmvDeadHAs()
{
	int tick = this->ECMWorldPtr->getTick();
	this->rmvHA(tick);
}

void ECM::rmvAllHAs()
{
	// reset HAlife map
	this->HAlife.clear();

	// reset total HA
	this->totalHA[write_t] = 0;

	// update ECM polymer map
	//#ifdef VISUALIZATION
	// Update ECM polymer map
	this->ECMWorldPtr->resetECM(this->index, m_hya);
	//#endif    // VISUALIZATION
#ifdef OPT_ECM
	this->set_dirty();
#endif
}


void ECM::addMonomerCol(float nColls)
{

#pragma omp atomic
	this->monomer_col[write_t] += nColls;

#ifdef OPT_ECM
	this->set_dirty();
#endif

}

void ECM::addPolymerCol(float nColls)
{

#pragma omp atomic
	this->polymer_col += nColls;

#ifdef OPT_ECM
	this->set_dirty();
#endif

}

int ECM::getnColl()
{
	return this->polymer_col;
}

bool ECM::isEmpty()
{
	float total = monomer_col[read_t] + monomer_ela[read_t] +
			fragmnt_col[read_t] + fragmnt_ela[read_t] +
			polymer_col + polymer_ela +
			totalHA[read_t] + fHA[read_t];

	return (total <= 0.0f);
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
