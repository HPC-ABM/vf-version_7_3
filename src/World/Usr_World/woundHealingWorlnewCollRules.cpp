/*
 * woundHealingWorld.cpp
 *
 * File contents: Contains the WHWorld class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *
 * Created on Jun 17, 2013, 6:59 PM
 *****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <cstring>
#include <algorithm>
#include <cstdlib>  
#define PI 3.14159

#include <string>
#include <sstream>
#include <omp.h>

#include "woundHealingWorld.h"
#include "../../enums.h"
#include "../../Utilities/timer.h"
#include "../../Utilities/error_utils.h"
#include "../../Utilities/input_utils.h"
#include "../../Utilities/mem_utils.h"

#include "../../Diffusion/DiffusionHelper.h"
#include "../../Diffusion/convolutionFFT_common.h"

using namespace std;


/********************************************
 * STATIC VARS INITIALIZATIONS              *
 ********************************************/
double WHWorld::clock = 0;
unsigned WHWorld::seed = 2000;
int WHWorld::bmcx = 0;
int WHWorld::bmcy = 0;
float WHWorld::epithickness = 0;
bool WHWorld::highTNFdamage = false;
float WHWorld::patchpermm = 0;
int WHWorld::initialTissue = 0; 

#ifdef MODEL_VOCALFOLD
int WHWorld::SLPxmax = 0;
int WHWorld::SLPxmin = 0;
int WHWorld::ILPxmax = 0;
int WHWorld::ILPxmin = 0;
int WHWorld::DLPxmax = 0;
int WHWorld::DLPxmin = 0;
float WHWorld::VFvolumefraction=0;
#endif

#ifndef CALIBRATION
int WHWorld::RVIS = 0;
int WHWorld::RVVS = 0;
int WHWorld::SSIS = 0;
int WHWorld::SSVS = 0;
#else
int WHWorld::RVIS = 5;
int WHWorld::RVVS = 10;
int WHWorld::SSIS = 10;
int WHWorld::SSVS = 10;
#endif
float WHWorld::thresholdTNFdamage = 10.0;
float WHWorld::thresholdMMP8damage = 10.0;
float WHWorld::sproutingFrequency[6] = {2.0, 4.0, 2.0, 4.0, 6.0, 12.0};
float WHWorld::sproutingAmount[14] = {20, 1, 8, 1, 8, 0.01, 1, 0.01, 1, 0.01, 1, 1, 0.01, 1};
float WHWorld::cytokineDecay[8] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.5};
float WHWorld::halfLifes_static[8] = {13, 13, 13, 13, 13, 13, 13, 30};


//WHWorld::WHWorld(double length, double width, double height, double plength) {
  WHWorld::WHWorld(double width, double length, double height, double plength) {
#ifdef DEBUG_WORLD
	// DEBUG
	maxHAsize = 0;
	HAfrags = 0;
	newOcolls = 0;
	deadneus = 0;
	newfibs = 0;
	actfibs = 0;
	deadfibs = 0;
	dead_afibs = 0;
	deactfibs = 0;
	deadfibs = 0;
	for (int icell = 0; icell < 4; icell++) {
		for (int ichem =0; ichem < 8; ichem++) {
			chemSecreted[icell][ichem] = 0.0;
			chemSecretedCoord[icell][ichem][0] = -1;
			chemSecretedCoord[icell][ichem][1] = -1;
			chemSecretedCoord[icell][ichem][2] = -1;
		}
	}
	for (int ichem = 0; ichem < 8; ichem++) {
		maxPatchChem[ichem][0] = 0.0;
		maxPatchChem[ichem][1] = 0.0;
		maxPatchChem[ichem][2] = 0.0;

		maxPatchChemCoord[ichem][0] = -1;
		maxPatchChemCoord[ichem][1] = -1;
		maxPatchChemCoord[ichem][2] = -1;

		maxOldPatchChem[ichem] = 0.0;

		maxOldPatchChemCoord[ichem][0] = -1;
		maxOldPatchChemCoord[ichem][1] = -1;
		maxOldPatchChemCoord[ichem][2] = -1;

		minOldPatchChem[ichem] = 90000.0;

		minOldPatchChemCoord[ichem][0] = -1;
		minOldPatchChemCoord[ichem][1] = -1;
		minOldPatchChemCoord[ichem][2] = -1;
	}
#endif

#ifdef PROFILE_ECM
this->ECMrepairTime = 0;
this->HAlifeTime    = 0;
this->ECMdangerTime = 0;
this->ECMscarTime   = 0;
#endif


  // Generate random seeds
  for(int i = 0; i < NUM_THREAD; i++) {
  	seeds[i] = seed + 17*i;//25234 + 17*i;
  }

  // Allocate memory for local lists of cell pointers to add
  for (int i = 0; i < MAX_NUM_THREADS; i++) {
  	localNewPlats[i] = new vector<Platelet*>;
  	localNewNeus[i] = new vector<Neutrophil*>;
  	localNewMacs[i] = new vector<Macrophage*>;
  	localNewFibs[i] = new vector<Fibroblast*>;
  }



  /********************************************
   * AGGREGRATED STATS INITIALIZATIONS        *
   ********************************************/
  totalOC   = 0, totalNC   = 0, totalFC   = 0;
  totalOE   = 0, totalNE   = 0, totalFE   = 0;
  totalHA   = 0, totalFHA  = 0;
  for(int i = 0; i < p_celltotal; i++)
	  this->totalCell[i] = 0;

  /********************************************
   * GRID SETUP                               *
   ********************************************/
  this->patchlength = plength;	// mm

  // Number of patches in x,y,z dimensions:
  int nx = width/patchlength;
  int ny = length/patchlength;
  int nz = height/patchlength;
  World::setupGrid(
  		nx,             // number of grid points (patches) in x dimension
  		ny,             // number of grid points (patches) in y dimension
  		nz,             // number of grid points (patches) in z dimension
  		0.0,            // min coordinates in x
  		width,          // max coordinates in x
  		0.0,            // min coordinates in y
  		length,         // max coordinates in y
  		0.0,            // min coordinates in z
  		height          // max coordinates in z
  );

  // Read input parameters (chem baseline, wound dimensions, initial cells) from config file
  int temp = WHWorld::userInput();
  cout << "length, width, height: " << length << " " << width << " " << height << " " << endl;
  cout << "nx, ny, nz: " << nx << " " << ny << " " << nz << " " << endl;

  cout << "   allocating ECM Managers (also Patches) with best-effort first touch" << endl;
  // Allocate and initialize Patches/ECM
  util::allocate<Patch>(&(this->worldPatch), nx*ny*nz);
  util::allocate<ECM>  (&(this->worldECM),   nx*ny*nz);

  cout << "worldPatch size: " << (nx)*(ny)*(nz) << endl;
  cout << "worldECM size: " << (nx)*(ny)*(nz) << endl;

  /* Try initializing Patches and ECMs with the threads that will access them
   * later since the default allocation policy on Linux platforms is
   * first-touch. This is a best-effort implementation, since we cannot
   * guarantee size of data accessed per thread to be an integer multiple of
   * page size. */

  // Calculate maximum number of collagen fibers allowed per patch
  // (w_patch mm/0.000086 mm)^2 - 1
  float numPerRow = patchlength/0.000086;
  float numPerPlane = numPerRow*numPerRow;
  this->maxFCol     = numPerPlane  - 1.0f;
  this->maxFCol /= COLL_UNIT;

  cout << "maximum number of collagen fibers per patch: " << this->maxFCol * COLL_UNIT << endl;

  this->maxChangeFCol = fibProductionTCollUB/collagenConversion;

#ifdef MODEL_3D
#pragma omp parallel for
#endif
  for (int iz = 0; iz < nz; iz++) {
#ifndef MODEL_3D
  #pragma omp parallel for
#endif
  	for (int iy = 0; iy < ny; iy++) {
  		for (int ix = 0; ix < nx; ix++) {
  			int in = ix + iy*nx + iz*nx*ny;
  			this->worldECM[in]   = ECM(ix, iy, iz, in);
  		}
  	}
  }

//#ifdef VISUALIZATION
  // Allocate memory for ECM map
  for (int ie = 0; ie < m_ecmtotal; ie++)
  {
  	this->ecmMap[ie]        = new float[np]();
  	this->ecmGradMap[ie]    = new float[np]();
  	this->ecmPreProcMap[ie] = new float[np]();
  }

//#endif

#ifdef MODEL_3D
#pragma omp parallel for
#endif
  for (int iz = 0; iz < nz; iz++) {
#ifndef MODEL_3D
  #pragma omp parallel for
#endif
  	for (int iy = 0; iy < ny; iy++) {
  		for (int ix = 0; ix < nx; ix++) {
  			int in = ix + iy*nx + iz*nx*ny;
  			this->worldPatch[in] = Patch(ix, iy, iz, in);
  		}
  	}
  }

  // Allocate memory for WHWorldChem
  this->WHWorldChem = new WHChemical(this->nx, this->ny, this->nz);


  // Define Class static variables and pointers
  WHWorld::patchpermm = nx/width;
  Agent::nx = this->nx;
  Agent::ny = this->ny;
  Agent::nz = this->nz;
  Agent::agentWorldPtr = this;
  Agent::agentPatchPtr = this->worldPatch;
  Agent::agentECMPtr = this->worldECM;
  ECM::ECMWorldPtr = this;
  ECM::ECMPatchPtr = this->worldPatch;
  WHChemical::chemWorldPtr = this;

  /********************************************
   * INITIALIZATION SUBROUTINES               *
   ********************************************/
  /* Define initial attributes of patches, damage, ECM, chem and cells
   *  based on user defined values (in config file) and traits of native tissue */
  cerr << "(i0/10) initializePatches ..." << endl;
  this->initializePatches();
  cerr << "(i1/10) initializeECM ..." << endl;
  this->initializeECM();
  cerr << "(i2/10) initializeChem ..." << endl;
  this->initializeChem();
  cerr << "(i3/10) initializeFibroblasts ..." << endl;
  this->initializeFibroblasts();
  cerr << "(i4/10) initializeMacrophages ..." << endl;
  this->initializeMacrophages();
  cerr << "(i5/10) initializeNeutrophils ..." << endl;
  this->initializeNeutrophils();
  cerr << "(i6/10) initializeDamage ..." << endl;
  this->initializeDamage();

  /* Calling update functions to synchronize read and write portion of the
   * attributes */
  //WHWorldChe->updateChemCPU();
  cerr << "(i7/10) updateCellsInitial ..." << endl;
  this->updateCellsInitial();  // Add cells to list before removal and updates
  cerr << "(i8/10) updateECMManagers ..." << endl;
  this->updateECMManagers();
#ifdef DEBUG_WORLD
    int nColl = 0;
    int nElas = 0;
    int nHA = 0;
    int nwColl = 0, nwElas = 0;
    for (int i = 0; i < (nx*ny*nz - 1); i++)
    {
        nColl += this->worldECM[i].ncollagen[read_t];
        nElas += this->worldECM[i].nelastin[read_t];
        nwColl += this->worldECM[i].ncollagen[write_t];
        nwElas += this->worldECM[i].nelastin[write_t];
        nHA += this->worldECM[i].getnHA();
    }
    cout << "   NEW collagens: " << nColl << endl;
    cout << "   NEW elastins: " << nElas << endl;
    cout << "   NEW wcollagens: " << nwColl << endl;
    cout << "   NEW welastins: " << nwElas << endl;
    cout << "   NEW HAs: " << nHA << endl;
#endif	// DEBUG_WORLD
  cerr << "(i9/10) updatePatches ..." << endl;
  this->updatePatches();
#ifdef DEBUG_WORLD
     nColl = 0;
     nElas = 0;
     nHA = 0;
     nwColl = 0, nwElas = 0;
    for (int i = 0; i < nx*ny*nz; i++)
    {
        nColl += this->worldECM[i].ncollagen[read_t];
        nElas += this->worldECM[i].nelastin[read_t];
        nwColl += this->worldECM[i].ncollagen[write_t];
        nwElas += this->worldECM[i].nelastin[write_t];
        nHA += this->worldECM[i].getnHA();
    }
    cout << "-----------------------" << endl;
    cout << "   NEW collagens: " << nColl << endl;
    cout << "   NEW elastins: " << nElas << endl;
    cout << "   NEW wcollagens: " << nwColl << endl;
    cout << "   NEW welastins: " << nwElas << endl;
    cout << "   NEW HAs: " << nHA << endl;
#endif	// DEBUG_WORLD

  this->updateCellStats();
  cout << "setupGrid complete" << endl;
}

WHWorld:: ~WHWorld(){
	cout << "WHWorld:: ~WHWorld()" << endl;

	free(this->D);
	free(this->HalfLifes);

#ifdef GPU_DIFFUSE

	cout << "	diffusion_helper::deallocConvCtxBuffers(this->chem_cctx)" << endl;
	diffusion_helper::deallocConvCtxBuffers(this->chem_cctx);
	free(this->chem_cctx);

       /**********************
	* Clear GPUs         *
	**********************/
/* 	for (int ig = 0; ig < N_GPU; ig++)
 	{
   		checkCudaErrors(cudaSetDevice(ig));
   		checkCudaErrors(cudaGetLastError());
   		checkCudaErrors(cudaDeviceReset());
 	}
*/
#endif

	for (int i = 0; i < MAX_NUM_THREADS; i++) {

		delete localNewPlats[i];
		delete localNewNeus[i];
		delete localNewMacs[i];
		delete localNewFibs[i];

	}

	cerr << " removing dead plats" << endl;
	int platsSize = plats.size();

	for (int i = 0; i < platsSize; i++) {

		Platelet* plat = plats.getDataAt(i);
		if (!plat) continue;
		plats.deleteData(i, DEFAULT_TID);
		delete plat;

	}

	cerr << " removing dead neus" << endl;
	int neusSize = neus.size();
#pragma omp parallel for
	for (int i = 0; i < neusSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		Neutrophil* neu = neus.getDataAt(i);
		if (!neu) continue;
		neus.deleteData(i, tid);
		delete neu;
	}

	cerr << " removing dead macs" << endl;
	int macsSize = macs.size();
#pragma omp parallel for
	for (int i = 0; i < macsSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		Macrophage* mac = macs.getDataAt(i);
		if (!mac) continue;
		macs.deleteData(i, tid);
		delete mac;
	}

	cerr << " removing dead fibs" << endl;
	int fibsSize = fibs.size();
#pragma omp parallel for
	for (int i = 0; i < fibsSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		Fibroblast* fib = fibs.getDataAt(i);
		if (!fib) continue;
		fibs.deleteData(i, tid);
		delete fib;
	}


	if (worldPatch != NULL)	delete [] worldPatch;
	if (worldECM != NULL) delete [] worldECM;

	cout << "WHWorld has been successfully destructed" << endl;
}

void destroyPlat(Platelet* &agent) {

	if (agent) {

		delete agent;
		//		agent = NULL;
	}
}

void destroyNeu(Neutrophil* &agent) {

	if (agent) {

		delete agent;
		//		agent = NULL;
	}
}

void destroyMac(Macrophage* &agent) {

	if (agent) {

		delete agent;
		//		agent = NULL;
	}
}

void destroyFib(Fibroblast* &agent) {

	if (agent) {

		delete agent;
		//		agent = NULL;
	}
}

void WHWorld::initializePatches() {

	// Define static world variables
	this->epithickness = epitheliumthickness;
	this->capRadius = capillaryradius;
	float bmcx_mm = x_max - epithickness;   // Basement membrane center x-coordinate (in mm)
	this->bmcx = nx - mmToPatch(epithickness); // basement membrane center (in patches)
	this->bmcy = 0;                         // Basement membrane center y-coordinate (in patches)

	// DEBUG


#ifdef MODEL_VOCALFOLD  

	// Define boundaries of lamina propria tissue layers (superficial, intermediate, deep)
	this->SLPxmax = (this->bmcx);
	this->SLPxmin = (this->bmcx - mmToPatch(fractionSLP*LPx)) < 0? 0:(this->bmcx - mmToPatch(fractionSLP*LPx));
	this->ILPxmax = (this->bmcx - mmToPatch(fractionSLP*LPx)) < 0? 0:(this->bmcx - mmToPatch(fractionSLP*LPx));
	this->ILPxmin = (this->bmcx - mmToPatch((fractionSLP+fractionILP)*LPx)) < 0? 0:(this->bmcx - mmToPatch((fractionSLP+fractionILP)*LPx));
	this->DLPxmax = (this->bmcx - mmToPatch((fractionSLP+fractionILP)*LPx)) < 0? 0:(this->bmcx - mmToPatch((fractionSLP+fractionILP)*LPx));
	this->DLPxmin = (this->bmcx - mmToPatch(LPx)) < 0? 0:(this->bmcx - mmToPatch(LPx));
	// Volume Fraction of ABM World to vocal fold (24.9mm x 1.6mm x 17.4mm); used to scale initial cell population
	this->VFvolumefraction = (x_max - x_min)*(y_max - y_min)*(z_max - z_min)/(LPx*LPy*LPz);
#endif
#ifndef CALIBRATION
	this->RVIS = 5;                         // Resonant voice impact stress
	this->SSIS = 10;                        // Spontaneous speech impact stress
	this->RVVS = 10;                        // Resonant voice vibratory stress
	this->SSVS = 10;                        // Spontaneous speech vibratory stress
#endif
	this->highTNFdamage = false;

	// Calculate number of capillaries
	double intercapXDistance = capillaryXdistance;	// Distance between capillaries along x (in mm)
	double intercapYDistance = capillaryYdistance;	// Distance between capillaries along y (in mm)
	this->capX.resize(ceil((bmcx_mm - x_min)/intercapXDistance));
	cout << "capX #: " << capX.size() << endl;
	this->capY.resize(ceil((y_max - y_min)/(intercapYDistance)));
	cout << "capY #: " << capY.size() << endl;

//	printf("x_max: %f\tepithickness: %f\tbmcx_mm: %f\n", x_max, epithickness, bmcx_mm);
//	printf("this->bmcx: %d\tmmToPatch(epithickness): %d\n", this->bmcx, mmToPatch(epithickness));
//	printf("bmcx: %d\tmmToPatch(capRadius): %d\tmmToPatch(intercapXDistance): %d\n", bmcx, mmToPatch(capRadius),
//			mmToPatch(intercapXDistance));


#ifdef RAT_VF
	int capRadiusP = 1;
#else
	int capRadiusP = mmToPatch(capRadius);
#endif

	// Calculate capillary center coordinates
	int maxCapX = 0;
	// x-coordinate:
	for (int icap = 0; icap < capX.size(); icap++) {
		int xTemp = (bmcx - capRadiusP) - icap*mmToPatch(intercapXDistance);
		this->capX[icap] = xTemp;
		maxCapX = maxCapX < xTemp? xTemp:maxCapX;
	}
	// y-coordinate:
	for (int icap = 0; icap < capY.size(); icap++) {
		this->capY[icap] = (capRadiusP) + icap*mmToPatch(intercapYDistance);
	}

	// Assign patches as either tissue or epithelium
	for (int iz = 0; iz < nz; iz++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {

				int in = ix + iy*nx + iz*nx*ny;
				// Epithelium
				if ((ix >= this->bmcx) && (ix < this->bmcx + mmToPatch(epithickness))) {
					this->worldPatch[in].type[read_t] = epithelium;
					this->worldPatch[in].type[write_t] = epithelium;
					this->worldPatch[in].color[read_t] = cepithelium;
					this->worldPatch[in].color[write_t] = cepithelium;
				}
				// Tissue
				if (ix < this->bmcx) {
					this->worldPatch[in].type[read_t] = tissue;
					this->worldPatch[in].type[write_t] = tissue;
					this->worldPatch[in].color[read_t] = ctissue;
					this->worldPatch[in].color[write_t] = ctissue;

#ifdef MODEL_VOCALFOLD
					// Define Vocal Fold tissue layers
					if (ix > this->SLPxmin) {       // Superficial LP (SLP)
						this->worldPatch[in].LP = SLP;

					} else if (ix > this->ILPxmin){	// Intermediate LP (ILP)
						this->worldPatch[in].LP = ILP;

					} else if (ix > this->DLPxmin){	// Deep LP (ILP)
						this->worldPatch[in].LP = DLP;
					} else { 			// Muscle
						this->worldPatch[in].LP = muscle;
					}
#endif
				}
			}
		}
	}

	cout << "finished assigning tissue" << endl;

	/* Assign patches that are within the capillary radius (capradius) from the
	 * capillary center coordinate (capX[ixCap], capY[iyCap]) as capillaries. */
#ifdef RAT_VF
//	// DEBUG vis
//	for (int iz = 0; iz < nz; iz++) {
//		for (int iy = 0; iy < ny; iy++) {
//			for (int ix = 0; ix < nx; ix++) {
//
//				int in = ix + iy*nx + iz*nx*ny;
//
//				if (worldPatch[in].type[read_t] == tissue || worldPatch[in].type[read_t] == capillary)
//					this->setECM(in, m_col, 0.1);
//
////				if ((iz == 0) || (iz == (this->nz - 1)))
////				{
////					worldPatch[in].type[read_t]   = capillary;
////					worldPatch[in].type[write_t]  = capillary;
////					worldPatch[in].color[read_t]  = ccapillary;
////					worldPatch[in].color[write_t] = ccapillary;
////				}
//			}
//		}
//	}


	for (int iyCap = 0; iyCap < capY.size(); iyCap ++) {
		for (int ixCap = 0; ixCap < capX.size(); ixCap ++) {
			for (int iz = 0; iz < nz; iz++) {
				int ix = capX[ixCap];
				int iy = capY[iyCap];

				int in = ix + iy*nx + iz*nx*ny;
				worldPatch[in].type[read_t]   = capillary;
				worldPatch[in].type[write_t]  = capillary;
				worldPatch[in].color[read_t]  = ccapillary;
				worldPatch[in].color[write_t] = ccapillary;
			}
		}
	}
#else
	for (int ixCap = 0; ixCap < capX.size(); ixCap ++) {

		int iyCapStart;
		if (ixCap % 2) iyCapStart = 0;
		else iyCapStart = 1;

		for (int iyCap = iyCapStart; iyCap < capY.size(); iyCap += 1) {// 2) {

			for (int iz = 0; iz < nz; iz++) {
				for (int iy = capY[iyCap]; iy <= capY[iyCap] + 2*mmToPatch(capRadius); iy++) {
					for (int ix = capX[ixCap] ; ix >= capX[ixCap] - 2*mmToPatch(capRadius); ix--) {
						int in = ix + iy*nx + iz*nx*ny;
						// Try another patch if this one is outside the bound of the world
						if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) continue;
						int a = (ix - (capX[ixCap] - mmToPatch(capRadius)))*(ix - (capX[ixCap] - mmToPatch(capRadius)));
						int b = (iy - (capY[iyCap] + mmToPatch(capRadius)))*(iy - (capY[iyCap] + mmToPatch(capRadius)));

						if (a + b <= mmToPatch(capRadius)*mmToPatch(capRadius)) {

							worldPatch[in].type[read_t]   = capillary;
							worldPatch[in].type[write_t]  = capillary;
							worldPatch[in].color[read_t]  = ccapillary;
							worldPatch[in].color[write_t] = ccapillary;

						}
					}
				}
			}
		}
	}
#endif

//	for (int iz = 0; iz < nz; iz++)
//	{
//		for (int iy = 0; iy < ny; iy++){
//			for (int ix = 0; ix < nx; ix++)
//			{
//				int in = ix + iy*nx + iz*nx*ny;
//				printf("%d, ", worldPatch[in].type[read_t]);
//			}
//			printf("\n");
//		}
//		printf("\n\n\n");
//	}


}


int WHWorld::getInitialDam()
{
	return initialDam;
}

int WHWorld::countDamage()
{
	return totaldamage;
}


void WHWorld::initializeChem(){
#ifdef GPU_DIFFUSE
	this->initializeChemGPU();
#else
	this->initializeChemCPU();
#endif
}



void WHWorld::initPatchChem()
{
#ifdef OPT_CHEM
	// Initialize chemical concentrations

	int isTissue = 0;

	for (int ic = 0; ic < this->typesOfBaseChem; ic++)
	{
		this->WHWorldChem->total[ic] = 0;
	}

#pragma omp parallel
	{
		float sum[8] = {0};
		float plevel = 0;
		/* Try initializing chemicals with the threads that will access them
		 * later since the default allocation policy on Linux platforms is
		 * first-touch. This is a best-effort implementation, since we cannot
		 * guarantee size of data accessed per thread to be an integer multiple
		 * of page size. */
#ifndef MODEL_3D
#pragma omp for
#endif  // !MODEL_3D
		for (int iz = 0; iz < this->nz; iz++) {
#ifdef MODEL_3D
#pragma omp for
#endif  // MODEL_3D
			for (int iy = 0; iy < this->ny; iy++) {
				for (int ix = 0; ix < this->nx; ix++) {
					int in = ix + iy*nx + iz*nx*ny;
					isTissue = (this->worldPatch[in].type[read_t] == tissue);

					for (int ic = 0; ic < typesOfBaseChem; ic++)
					{
						WHWorldChem->setDchem(ic, in, 0);
						/* Baseline chemical concentrations are initialized in tissue
						 * (not epithelium or capillaries) */ //TODO(Kim): INSERT REF?
						plevel = (this->baselineChem[ic]/WHWorld::initialTissue);//*isTissue;
						WHWorldChem->setPchem(ic, in, plevel);
						sum[ic] += plevel;
					}

					// Initialize chemical gradient levels that agents are attracted by
					float patchIL1 = this->WHWorldChem->getPchem(IL1beta, in);
					float patchIL6 = this->WHWorldChem->getPchem(IL6, in);
					float patchIL8 = this->WHWorldChem->getPchem(IL8, in);
					float patchTNF = this->WHWorldChem->getPchem(TNF, in);
					float patchTGF = this->WHWorldChem->getPchem(TGF, in);
					float patchFGF = this->WHWorldChem->getPchem(FGF, in);
					float patchcollagen = this->worldECM[in].fcollagen[read_t];
					float grad = patchIL1 + patchTNF + patchTGF + patchFGF + patchcollagen;
					this->WHWorldChem->setGrad(FIBgrad, in, patchTGF);  // TODO(Kim): INSERT REF?
					this->WHWorldChem->setGrad(NEUgrad, in, grad + patchIL6 + patchIL8);  // TODO(Kim): INSERT REF?
					this->WHWorldChem->setGrad(MACgrad, in, grad + this->worldECM[in].felastin[read_t]);  // TODO(Kim): INSERT REF?

				}
			}
		}
#pragma omp critical
		{
			//Initialize total chemical concentration
			for (int ic = 0; ic < typesOfBaseChem; ic++)
				this->WHWorldChem->total[ic] += sum[ic];
		}
	}
	cout << "results from inside initialization: " << endl;
	for (int ic = 0; ic< typesOfBaseChem; ic++)
	{
		cout << "	" << this->WHWorldChem->total[ic] << endl;
	}
#else

	int countTissue = WHWorld::initialTissue;
	WHWorldChem->resetTotals();



#pragma omp parallel
{
	// local sum
	double *sum = new double[typesOfBaseChem];
	for (int ic = 0; ic < typesOfBaseChem; ic++)
	{
		sum[ic] = 0.0f;
	}

#ifdef MODEL_3D
#pragma omp for
#endif
	for (int iz = 0; iz < this->nz; iz++) {
		/* Try initializing chemicals with the threads that will access them
		 * later since the default allocation policy on Linux platforms is
		 * first-touch. This is a best-effort implementation, since we cannot
		 * guarantee size of data accessed per thread to be an integer multiple
		 * of page size. */
#ifndef MODEL_3D
#pragma omp for
#endif
		for (int iy = 0; iy < this->ny; iy++) {
			for (int ix = 0; ix < this->nx; ix++) {
				int in = ix + iy*nx + iz*nx*ny;
				this->WHWorldChem->dTNF[in] = 0;
				this->WHWorldChem->dTGF[in] = 0;
				this->WHWorldChem->dFGF[in] = 0;
				this->WHWorldChem->dMMP8[in] = 0;
				this->WHWorldChem->dIL1beta[in] = 0;
				this->WHWorldChem->dIL6[in] = 0;
				this->WHWorldChem->dIL8[in] = 0;
				this->WHWorldChem->dIL10[in] = 0;

				/* Baseline chemical concentrations are initialized in tissue
				 * epithelium and capillaries */ //TODO(Kim): INSERT REF?
				// DEBUG rat
				this->WHWorldChem->pTNF[in] = this->baselineChem[TNF]/countTissue;
				this->WHWorldChem->pTGF[in] = this->baselineChem[TGF]/countTissue;
				this->WHWorldChem->pFGF[in] = this->baselineChem[FGF]/countTissue;
				this->WHWorldChem->pMMP8[in] = this->baselineChem[MMP8]/countTissue;
				this->WHWorldChem->pIL1beta[in] = this->baselineChem[IL1beta]/countTissue;
				this->WHWorldChem->pIL6[in] = this->baselineChem[IL6]/countTissue;
				this->WHWorldChem->pIL8[in] = this->baselineChem[IL8]/countTissue;
				this->WHWorldChem->pIL10[in] = this->baselineChem[IL10]/countTissue;

				// Initialize chemical gradient levels that agents are attracted by
				float patchIL1 = this->WHWorldChem->pIL1beta[in];
				float patchIL6 =  this->WHWorldChem->pIL6[in];
				float patchIL8 = this->WHWorldChem->pIL8[in];
				float patchTNF = this->WHWorldChem->pTNF[in];
				float patchTGF = this->WHWorldChem->pTGF[in];
				float patchFGF = this->WHWorldChem->pFGF[in];
				float patchcollagen = this->worldECM[in].fcollagen[read_t];
				float grad = patchIL1 + patchTNF + patchTGF + patchFGF + patchcollagen;


				this->WHWorldChem->pfibgrad[in] = patchTGF;  // TODO(Kim): INSERT REF?
				this->WHWorldChem->pneugrad[in] = grad + patchIL6 + patchIL8;  // TODO(Kim): INSERT REF?
				this->WHWorldChem->pmacgrad[in] = grad + this->worldECM[in].felastin[read_t];  // TODO(Kim): INSERT REF?

				sum[TNF]     += this->WHWorldChem->pTNF[in];
				sum[TGF]     += this->WHWorldChem->pTGF[in];
				sum[FGF]     += this->WHWorldChem->pFGF[in];
				sum[MMP8]    += this->WHWorldChem->pMMP8[in];
				sum[IL1beta] += this->WHWorldChem->pIL1beta[in];
				sum[IL6]     += this->WHWorldChem->pIL6[in];
				sum[IL8]     += this->WHWorldChem->pIL8[in];
				sum[IL10]    += this->WHWorldChem->pIL10[in];
				
			}
		}
	}
#pragma omp critical
	{
					//Initialize total chemical concentration
					this->WHWorldChem->total[TNF]  += sum[TNF];
					this->WHWorldChem->total[TGF]  += sum[TGF];
					this->WHWorldChem->total[FGF]  += sum[FGF];
					this->WHWorldChem->total[MMP8] += sum[MMP8];
					this->WHWorldChem->total[IL1beta] += sum[IL1beta];
					this->WHWorldChem->total[IL6] += sum[IL6];
					this->WHWorldChem->total[IL8] += sum[IL8];
					this->WHWorldChem->total[IL10] += sum[IL10];
	}
}
	cout << "results from inside initialization: " << this->WHWorldChem->total[TNF] << " ";
	cout << this->WHWorldChem->total[TGF] << " " << this->WHWorldChem->total[FGF];
	cout << " " << this->WHWorldChem->total[MMP8]
		<< " " << this->WHWorldChem->total[IL1beta];
	cout << " " << this->WHWorldChem->total[IL6] << " ";
	cout << this->WHWorldChem->total[IL8] << " " << this->WHWorldChem->total[IL10] << endl;
#endif
}

void WHWorld::initializeChemCPU() {

	this->typesOfBaseChem = this->baselineChem.size();
	this->typesOfChem   	= typesOfBaseChem*2 + 3;
// DEBUG(*)
	int countTissue = this->np;//this->countPatchType(tissue);
	WHWorld::initialTissue = countTissue;

#ifdef OPT_CHEM

#else	// OPT_CHEM
	WHWorldChem->linkAttributes();
#endif	// OPT_CHEM

	/*
	 * Initialize chemical concentrations on each patch w.r.t. user baseline input.
	 * Initialize total chemical accordingly.
	 */
	if (typesOfBaseChem == 8)
	{
		initPatchChem();
	} else {
		if (util::ABMerror(
				1,
				"Error initializing chemicals!!",
				__FILE__,
				__LINE__))
			exit(1);
	}
	cout << "Finished initializing chem" << endl;

}

#ifdef GPU_DIFFUSE

void WHWorld::initKernelConstants(float* lambda,
		float *gamma,
		int    numChem,
		c_ctx *kernel_cctx)
{


#ifdef PRINT_KERNEL
	cerr << "Allocated " << numChem << " elements for lambda (" << lambda << ")  and gamma ("
			<< gamma << ")" << endl;
	cerr << "D: " << D << "		halflife: " << HalfLifes << endl;
#endif
	float maxD = *std::max_element(D, D+numChem);
	printf("maxD: %f\n", maxD);

	float dx	= (this->dx) * 1000.0f;	// um
	float dx2	= dx * dx;
#ifdef RAT_VF
	float dt	= 0.5;		// TODO: remove hard code
#else
	float dt	= floor(dx2/(6*maxD));	// 2.0;
#endif

	kernel_cctx->dx2 = dx2;
	kernel_cctx->dt  = dt ;

	for (int ic = 0; ic < numChem; ic++) {
		lambda[ic] = (D[ic]*dt)/dx2;
		// Overwrite halflifes from config file
		// Use instead values from sensitivity analysis input file (Sample.txt)

		HalfLifes[ic] = halfLifes_static[ic] * 60;		// minutes * 60
		// TODO: get rid of gamma?
		gamma[ic] = 1 - pow(2, -(1/HalfLifes[ic]));//0.00;
		// Initizalize baseline per patch
		// DEBUG rat
		this->pbaseline[ic] = this->baselineChem[ic] / WHWorld::initialTissue;
		printf("D: %f\tdt: %f\tdx: %f\n", D[ic], dt, dx);
		printf("pbaseline[%d]:\t\t%.12f\tlambda: %f\n", ic, pbaseline[ic], lambda[ic]);
	}

	this->WHWorldChem->pbaseline = this->pbaseline;

#ifdef PRINT_KERNEL
	cerr << "Done--" << endl;

	cerr << "Initializing kernel and chemical data context" << endl;
#endif
}



void WHWorld::initializeChemGPU() {

#ifdef MODEL_3D
    diffusion_helper::findGPUs();
    //DEBUG:
    printf("++++ in initChemGPU: \n");
    reportMemUsageGPU();
#endif  // MODEL_3D



    this->typesOfBaseChem = this->baselineChem.size();
    this->typesOfChem     = typesOfBaseChem*2 + 3;

	const int numChem = this->typesOfBaseChem;
// DEBUG(*)
	int countTissue = this->np;//this->countPatchType(tissue);
	WHWorld::initialTissue = countTissue;

	/********************************************
	 * Kernel Computation	                    *
	 ********************************************/

	// Kernel convolution context
	c_ctx kernel_cctx;

	// Allocate and set dimensions in chemical convolution context
	this->chem_cctx         = new c_ctx();

	float *lambda = new float[numChem];
	float *gamma  = new float[numChem];
	pbaseline     = new float[numChem];

	initKernelConstants(lambda, gamma, numChem, &kernel_cctx);
	diffusion_helper::initKernelDimensions(&kernel_cctx);

	diffusion_helper::printContext( kernel_cctx);

	diffusion_helper::initChemDimensions(this->chem_cctx,
	            kernel_cctx, this->nx, this->ny, this->nz);
  diffusion_helper::allocKernelBuffers(&kernel_cctx);
  diffusion_helper::allocChemBuffers(this->chem_cctx, this->WHWorldChem, this->nx, this->ny, this->nz);

	diffusion_helper::prepKernel(this->chem_cctx, kernel_cctx, lambda, gamma, numChem);

	// After context allocation and initialization completes, print values
//	diffusion_helper::printContext( kernel_cctx);
	diffusion_helper::printContext(*chem_cctx);

	// Kernel buffer pointers deallocations
	diffusion_helper::deallocConvCtxBuffers(&kernel_cctx);


	/********************************************
	 * GPU Diffusion Preparation	              *
	 ********************************************/

	cout << "number of chem allocated is "<< this->typesOfChem << endl;
	cout << "number of chem allocated for GPU diffusion is " << typesOfBaseChem << endl;


	// Instantiate diffusion manager
	this->diffuserPtr = new Diffuser(this->typesOfBaseChem,
	                                    this->chem_cctx,
	                                    &kernel_cctx,
	                                    this->WHWorldChem);

#ifdef OPT_CHEM

#else	// OPT_CHEM
	WHWorldChem->linkAttributes();
#endif	// OPT_CHEM

	/*
	 * Initialize chemical concentrations on each patch w.r.t. user baseline input.
	 * Initialize total chemical accordingly.
	 */
	if (typesOfBaseChem == 8)
	{
		initPatchChem();
	} else {
		if (util::ABMerror(
				1,
				"Error initializing chemicals!!",
				__FILE__,
				__LINE__))
			exit(1);
	}
	cout << "Finished initializing chem" << endl;

	//	sdkDeleteTimer(&hTimer);
#ifdef MODEL_3D

#else

	delete lambda;
	delete gamma;

#endif

}
#endif // GPU_DIFFUSE

void WHWorld::initializeFibroblasts() {

	// Instantiate Fibroblast list
#ifdef RAT_VF
	fibs = ArrayChain<Fibroblast *>(DEFAULT_DATA_SMALL, 4, NULL, NULL);//&destroyFib);
#else
	fibs = ArrayChain<Fibroblast *>(DEFAULT_DATA_XXLARGE, 4, NULL, NULL);//&destroyFib);
#endif

#ifdef MODEL_VOCALFOLD

	int ymin = 0;
	int ymax = ny;
	int zmin = 0;
	int zmax = nz;
	// Vocal Fold celluarlity vary along five sections (each representing 20% of LP)
	// [Ref] Catten, Michael, et al. "Analysis of cellular location and concentration in vocal fold
	// lamina propria" Otolaryngology--Head and Neck Surgery 118.5 (1998):663-667.
	float twentypercent = bmcx*0.2;
	// Initial Fibroblast population scaled to volume fraction of ABM world to vocal fold
	int initialFib = this->initialCells[0]*VFvolumefraction;
	int initialAFib = this->initialCells[3]*VFvolumefraction;
	cout << " VF volume fraction:" << VFvolumefraction << endl;
	cout << " initialfib:" << initialFib << endl;
	cout << " initialAfib:" << initialAFib << endl;

	/*
	 * Non-uniform cellularity in vocal fold lamina propria with depth
	 * Sprout fibroblast in variable concentrations at 20% intervals along LP depth(x-direction)
	 * [Ref] Catten, Michael, et al. "Analysis of cellular location and concentration in vocal fold
	 * lamina propria." Otolaryngology--Head and Neck Surgery 118.5 (1998): 663-667.
	 */

	// Sprout initial fibroblast in vocal fold lamina propria layers (SLP, ILP, DLP)


	sproutAgent(                    // Sprout fibroblast in superficial 0-20% of LP
			initialFib*fibroblastOne,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			fibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - twentypercent,   //  -- left
			bmcx,                   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 20-40% LP
			initialFib*fibroblastTwo,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			fibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 2*twentypercent, //  -- left
			bmcx - twentypercent,   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 40-60% LP
			initialFib*fibroblastThree,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			fibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 3*twentypercent, //  -- left
			bmcx - 2*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 60-80% LP
			initialFib*fibroblastFour,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			fibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 4*twentypercent, //  -- left
			bmcx - 3*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 80-100% LP
			initialFib*fibroblastFive,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			fibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 5*twentypercent, //  -- left
			bmcx - 4*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear


	// Activated Fibroblasts

	sproutAgent(                    // Sprout fibroblast in superficial 0-20% of LP
			initialAFib*afibroblastOne,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			afibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - twentypercent,   //  -- left
			bmcx,                   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 20-40% LP
			initialAFib*afibroblastTwo,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			afibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 2*twentypercent, //  -- left
			bmcx - twentypercent,   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 40-60% LP
			initialAFib*afibroblastThree,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			afibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 3*twentypercent, //  -- left
			bmcx - 2*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 60-80% LP
			initialAFib*afibroblastFour,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			afibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 4*twentypercent, //  -- left
			bmcx - 3*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear
	sproutAgent(                    // Sprout fibroblast in superficial 80-100% LP
			initialAFib*afibroblastFive,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			afibroblast,             // Type of agent to sprout
			// Physical Boundary
			bmcx - 5*twentypercent, //  -- left
			bmcx - 4*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax);                  //  -- rear

	cout << "Finished initializing fibs" << endl;
#else
	/* Sprout initial fibroblast randomly on tissue within bounds using reservoir
	 * sampling */ // TODO(Kim): INSERT REF?
	int initialFib = this->initialCells[0];
	//	int xmin = (this->nx - 1)/2 - mmToPatch(0.6);
	int xmin = 0;
	int xmax = nx;
	int ymin = 0;
	int ymax = ny;
	int zmin = 0;
	int zmax = nz;
	sproutAgent(initialFib,         // Number of cells to sprout
			tissue,         // Type of patch to sprout on
			fibroblast,     // Type of agent to sprout
			// Physical Boundary
			xmin,           //  -- left
			xmax,           //  -- right
			ymin,           //  -- top
			ymax,           //  -- bottom
			zmin,           //  -- front
			zmax);          //  -- rear

	cout << "Finished initializing fibs" << endl;
#endif
}

void WHWorld::initializeNeutrophils() {

	// Instantiate Neutrophil list
#ifdef RAT_VF
	neus = ArrayChain<Neutrophil*>(DEFAULT_DATA_SMALL, 4, NULL, NULL);  //WHWorld::destroyNeu);
#else
	neus = ArrayChain<Neutrophil*>(DEFAULT_DATA_XLARGE, 4, NULL, NULL);  //WHWorld::destroyNeu);
#endif
	int ymin = 0;
	int ymax = ny;
	int zmin = 0;
	int zmax = nz;
	int xmin = 0;
	int xmax = nx;

#ifdef MODEL_VOCALFOLD

	// Initial Neutrophil population scaled to volume fraction of ABM world to vocal fold
	size_t initNeutrophils = this->initialCells[1]*VFvolumefraction;
	cout << " initialneu:" << initNeutrophils << endl;

	// Neutrophil in capillary
	sproutAgent(initNeutrophils,    // Number of cells to sprout
			capillary,      // Type of patch to sprout on
			neutrophil,     // type of agent to sprout
			// Physical Boundary
			xmin,           //  -- left
			xmax,           //  -- right
			ymin,           //  -- top
			ymax,           //  -- bottom
			zmin,           //  -- front
			zmax);          //  -- rear


	cout << "Finished initializing neus" << endl;
#else
	/* Sprout initial neutrophil randomly on capillary within bounds using
	 * reservoir sampling */ // TODO(Kim): INSERT REF?
	size_t initNeutrophils = this->initialCells[1];
	sproutAgent(initNeutrophils,    // Number of cells to sprout
			capillary,      // Type of patch to sprout on
			neutrophil,     // type of agent to sprout
			// Physical Boundary
			xmin,           //  -- left
			xmax,           //  -- right
			ymin,           //  -- top
			ymax,           //  -- bottom
			zmin,           //  -- front
			zmax);          //  -- rear

	cout << "Finished initializing neus" << endl;
#endif
}

void WHWorld::initializeMacrophages() {

	// Instantiate Macrophage list
#ifdef RAT_VF
	macs = ArrayChain<Macrophage*>(DEFAULT_DATA_SMALL, 4, NULL, NULL);  //WHWorld::destroyMac);
#else
	macs = ArrayChain<Macrophage*>(DEFAULT_DATA_XXLARGE, 4, NULL, NULL);  //WHWorld::destroyMac);
#endif
	int ymin = 0;
	int ymax = ny;
	int zmin = 0;
	int zmax = nz;

#ifdef MODEL_VOCALFOLD

	// Vocal Fold celluarlity vary along five sections (each representing 20% of LP)
	// [Ref] Catten, Michael, et al. "Analysis of cellular location and concentration in vocal fold
	// lamina propria" Otolaryngology--Head and Neck Surgery 118.5 (1998):663-667.
	float twentypercent = bmcx*0.2;
	// Initial Neutrophil population scaled to volume fraction of ABM world to vocal fold
	int initMacrophages = this->initialCells[2]*VFvolumefraction;
	cout << " initialmac:" << initMacrophages << endl;

	// Sprout half initial Macrophages in tissue
	cout << " initialmac: 1" << endl;
	sproutAgent(
			0.5*initMacrophages*macrophageOne,// Initial cells in section 1 (superficial 0-20% of LP)
			tissue,                 // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - twentypercent,   //  -- left, superficial 0% LP
			bmcx,                   //  -- right, superficial 20% of LP
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			tissue);                // Sprout in tissue
	cout << " initialmac: 2" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 20-40% LP
			0.5*initMacrophages*macrophageTwo,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 2*twentypercent, //  -- left
			bmcx - twentypercent,   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			tissue);                // Sprout in tissue
	cout << " initialmac: 3" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 40-60% LP
			0.5*initMacrophages*macrophageThree,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 3*twentypercent, //  -- left
			bmcx - 2*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			tissue);                // Sprout in tissue
	cout << " initialmac: 4" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 60-80% LP
			0.5*initMacrophages*macrophageFour,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 4*twentypercent, //  -- left
			bmcx - 3*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			tissue);                // Sprout in tissue
	cout << " initialmac: 5" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 80-100% LP
			0.5*initMacrophages*macrophageFive,// Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 5*twentypercent, //  -- left
			bmcx - 4*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			tissue);                // Sprout in tissue

	// Sprout half initial Macrophages in capillary
	cout << " initialmac: 6" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 0-20% of LP
			0.5*initMacrophages*macrophageOne,// Number of cells to sprout
			capillary,              // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - twentypercent,   //  -- left
			bmcx,                   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			blood);                 // Sprout in blood
	cout << " initialmac: 7" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 20-40% LP
			0.5*initMacrophages*macrophageTwo,// Number of cells to sprout
			capillary,              // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 2*twentypercent, //  -- left
			bmcx - twentypercent,   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			blood);                 // Sprout in blood
	cout << " initialmac: 8" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 40-60% LP
			0.5*initMacrophages*macrophageThree,// Number of cells to sprout
			capillary,              // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 3*twentypercent, //  -- left
			bmcx - 2*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			blood);                 // Sprout in blood
	cout << " initialmac: 9" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 60-80% LP
			0.5*initMacrophages*macrophageFour,// Number of cells to sprout
			capillary,              // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 4*twentypercent, //  -- left
			bmcx - 3*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			blood);                 // Sprout in blood
	cout << " initialmac: 10" << endl;
	sproutAgent(                    // Sprout macrophage in superficial 80-100% LP
			0.5*initMacrophages*macrophageFive,// Number of cells to sprout
			capillary,              // Type of patch to sprout on
			macrophag,              // Type of agent to sprout
			// Physical Boundary
			bmcx - 5*twentypercent, //  -- left
			bmcx - 4*twentypercent, //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			blood);                 // Sprout in blood

	cout << "Finished initializing macs" << endl;

#else
	/* Sprout initial Macrophage half in capillary and half in tissue using
	 * reservoir sampling */ // TODO(Kim): INSERT REF?
	int initMacrophages = this->initialCells[2];
	//	int xmin = (nx - 1)/2 + mmToPatch(0.8);
	int xmin = 0;
	int xmax = nx;
	sproutAgent(0.5*initMacrophages,        // Number of cells to sprout
			capillary,              // Type of patch to sprout on
			macrophag,              // type of agent to sprout
			// Physical Boundary
			xmin,                   //  -- left
			xmax,                   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			blood);                 // Sprout in blood

	sproutAgent(0.5*initMacrophages,        // Number of cells to sprout
			tissue,                 // Type of patch to sprout on
			macrophag,              // type of agent to sprout
			// Physical Boundary
			xmin,                   //  -- left
			xmax,                   //  -- right
			ymin,                   //  -- top
			ymax,                   //  -- bottom
			zmin,                   //  -- front
			zmax,                   //  -- rear
			tissue);                // Sprout in tissue

	cout << "Finished initializing macs" << endl;
#endif
}

void WHWorld::initializeECM() {

	/*******************************************
	 * COLLAGEN                                *
	 *******************************************/
	double collagenGap;
	int ix;

#ifdef MODEL_VOCALFOLD
	int ymin = 0;
	int ymax = ny;
	int zmin = 0;
	int zmax = nz;
	cout << " SLPxmin, SLPxmax: " << this->SLPxmin << " , " << this->SLPxmax << endl;
	cout << " ILPxmin, ILPxmax: " << this->ILPxmin << " , " << this->ILPxmax << endl;
	cout << " DLPxmin, DLPxmax: " << this->DLPxmin << " , " << this->DLPxmax << endl;

	/* ECM (collagen, elastin, HA) fractional area in vocal fold layers (superficial LP (SLP),
	 * intermediate LP (ILP), deep LP (DLP))
	 * [1] Hahn, Mariah S., et al. "Quantitative and comparative studies of the vocal fold extracellular
	 * matrix II: collagen." Annals of Otology, Rhinology & Laryngology 115.3 (2006): 225-232.
	 * [2] Hahn, Mariah S., et al. "Quantitative and comparative studies of the vocal fold extracellular
	 * matrix I: elastic fibers and hyaluronic acid." Annals of Otology, Rhinology & Laryngology 115.2 (2006): 156-164.
	 * NOTE: assume ECM fractional area is proportional to fractional volume for 3D case*/
	double collagenFracAreaSLP = 0.71;//0.307;
	double collagenFracAreaILP = 0.32;//0.360;
	double collagenFracAreaDLP = 0.27;//0.333;

	// Volume (number of patches) in each SLP, ILP, DLP
	double volumeSLP = (SLPxmax - SLPxmin)*(ymax - ymin)*(zmax - zmin);
	double volumeILP = (ILPxmax - ILPxmin)*(ymax - ymin)*(zmax - zmin);
	double volumeDLP = (DLPxmax - DLPxmin)*(ymax - ymin)*(zmax - zmin);
	cout << " volumeSLP, ILP, DLP: " << volumeSLP << " , " << volumeILP << " , " << volumeDLP <<endl;

	// Volume (number of patches) of ECM (collagen, elastin, HA) in each SLP, ILP, DLP
	double collagenSLP = collagenFracAreaSLP*volumeSLP;
	double collagenILP = collagenFracAreaILP*volumeILP;
	double collagenDLP = collagenFracAreaDLP*volumeDLP;
	cout << " collagenSLP, ILP, DLP: " << collagenSLP << " , " << collagenILP << " , " << collagenDLP <<endl;

	sproutAgentInArea(collagenSLP,                  // Number of patches to sprout agents
			tissue,                       // patch type
			nc,                           // original collagen agent type
			SLPxmin,                      // lowest x-coordinate agent could be sprouted at
			SLPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // collagen in SLP
	sproutAgentInArea(collagenILP,                  // Number of patches to sprout agents
			tissue,                       // patch type
			nc,                           // original collagen agent type
			ILPxmin,                      // lowest x-coordinate agent could be sprouted at
			ILPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // collagen in ILP
	sproutAgentInArea(collagenDLP,                  // Number of patches to sprout agents
			tissue,                       // patch type
			nc,                           // original collagen agent type
			DLPxmin,                      // lowest x-coordinate agent could be sprouted at
			DLPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // collagen in DLP
	cout << "Finished initializing collagens" << endl;

#else
	// Sprout original collagen on tissue every 0.32mm  // TODO(Kim): INSERT REF?
	collagenGap = 0.32;             // Spacing between collagen (in mm)
	ix = this->bmcx - epithickness;
	for (int iz = 0; iz < nz ; iz++) {
		for (int iy = 0; iy < ny ; iy++) {

			int in = ix + iy*nx + iz*nx*ny;
			if ((this->worldPatch[in].type[read_t] == tissue) && (iy%mmToPatch(collagenGap) == 0)) {

				this->worldPatch[in].initcollagen = true;
				this->worldECM[in].ncollagen[write_t] = this->worldECM[in].ncollagen[read_t] + 1;
				this->worldECM[in].isEmpty();
			}
		}
	}
	cout << "Finished initializing collagens" << endl;
#endif

	/*******************************************
	 * ELASTIN                                 *
	 *******************************************/
#ifdef MODEL_VOCALFOLD

	/* ECM (collagen, elastin, HA) fractional area in vocal fold layers (superficial LP (SLP),
	 * intermediate LP (ILP), deep LP (DLP))
	 * [1] Hahn, Mariah S., et al. "Quantitative and comparative studies of the vocal fold extracellular
	 * matrix II: collagen." Annals of Otology, Rhinology & Laryngology 115.3 (2006): 225-232.
	 * [2] Hahn, Mariah S., et al. "Quantitative and comparative studies of the vocal fold extracellular
	 * matrix I: elastic fibers and hyaluronic acid." Annals of Otology, Rhinology & Laryngology 115.2 (2006): 156-164.
	 * Note: assume ECM fractional area is proportional to fractional volume and same for 2 and 3D case*/
	double elastinFracAreaSLP = 0.312;
	double elastinFracAreaILP = 0.354;
	double elastinFracAreaDLP = 0.335;

	// Volume (number of patches) in each SLP, ILP, DLP
	volumeSLP = (SLPxmax - SLPxmin)*(ymax - ymin)*(zmax - zmin);
	volumeILP = (ILPxmax - ILPxmin)*(ymax - ymin)*(zmax - zmin);
	volumeDLP = (DLPxmax - DLPxmin)*(ymax - ymin)*(zmax - zmin);

	// Volume (number of patches) of ECM (collagen, elastin, HA) in each SLP, ILP, DLP
	double elastinSLP = elastinFracAreaSLP*volumeSLP;
	double elastinILP = elastinFracAreaILP*volumeILP;
	double elastinDLP = elastinFracAreaDLP*volumeDLP;
	cout << " elastinSLP, ILP, DLP: " << elastinSLP << " , " << elastinILP << " , " << elastinDLP <<endl;

	sproutAgentInArea(elastinSLP,                  // Number of agents to sprout
			tissue,                       // patch type
			ne,                           // original elastin agent type
			SLPxmin,                      // lowest x-coordinate agent could be sprouted at
			SLPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // elastin in SLP
	sproutAgentInArea(elastinILP,                  // Number of agents to sprout
			tissue,                       // patch type
			ne,                           // original elastin agent type
			ILPxmin,                      // lowest x-coordinate agent could be sprouted at
			ILPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // elastin in ILP
	sproutAgentInArea(elastinDLP,                  // Number of agents to sprout
			tissue,                       // patch type
			ne,                           // original elastin agent type
			DLPxmin,                      // lowest x-coordinate agent could be sprouted at
			DLPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // elastin in DLP
	cout << "Finished initializing elastin" << endl;
#else
	// Sprout original elastin at basement membrane on tissue, every 0.24mm
	double elastinGap = 0.24;              // Spacing between elastin (in mm)
	ix = this->bmcx - epithickness;

	for (int iy = 0; iy < ny; iy++) {
		for (int iz = 0; iz < nz; iz++) {

			int in = ix + iy*nx + iz*nx*ny;
			if ((this->worldPatch[in].type[read_t] == tissue) && (iy%mmToPatch(elastinGap) == 0)) {

				this->worldPatch[in].initelastin = true;
				this->worldECM[in].nelastin[write_t] = this->worldECM[in].nelastin[read_t] + 1;
				this->worldECM[in].isEmpty();
			}
		}
	}
	// Store indices of tissue patches not at basement membrane
	vector <int> patchlist;
	for (int ix = mmToPatch(0.2); ix < this->bmcx-mmToPatch(0.16); ix++) {  // TODO(Kim): INSERT REF? (Why 0.02 and 0.16?)
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {

				int in = ix + iy*nx + iz*nx*ny;
				if ((this->worldPatch[in].type[read_t] == tissue)) {
					patchlist.push_back(in);
				}
			}
		}
	}
	// Sprout oelastin on 0.04% of tissue patches  // TODO(Kim): INSERT REF?
	for (int i = 0; i < 0.0004*nx*ny*nz; i++) {

		if (patchlist.size() <=0) continue;
		int randInt = rand() % patchlist.size();
		int patchIndex = patchlist[randInt];
		this->worldPatch[patchIndex].initelastin = true;
		this->worldECM[patchIndex].nelastin[write_t] = this->worldECM[patchIndex].nelastin[read_t] + 1;
		this->worldECM[patchIndex].isEmpty();
	}
	cout << "Finished initializing elastins" << endl;
#endif
	/*******************************************
	 * HYALURONAN                              *
	 *******************************************/

#ifdef MODEL_VOCALFOLD
	/*
	 * ECM (collagen, elastin, HA) fractional area in vocal fold layers (superficial LP (SLP),
	 * intermediate LP (ILP), deep LP (DLP))
	 * [1] Hahn, Mariah S., et al. "Quantitative and comparative studies of the vocal fold extracellular
	 * matrix II: collagen." Annals of Otology, Rhinology & Laryngology 115.3 (2006): 225-232.
	 * [2] Hahn, Mariah S., et al. "Quantitative and comparative studies of the vocal fold extracellular
	 * matrix I: elastic fibers and hyaluronic acid." Annals of Otology, Rhinology & Laryngology 115.2 (2006): 156-164.
	 * Note: assume ECM fractional area is proportional to fractional volume and same for 2 and 3D case
	 */
	double HAFracAreaSLP = 0.312;
	double HAFracAreaILP = 0.354;
	double HAFracAreaDLP = 0.335;

	// Volume (number of patches) in each SLP, ILP, DLP
	volumeSLP = (SLPxmax - SLPxmin)*(ymax - ymin)*(zmax - zmin);
	volumeILP = (ILPxmax - ILPxmin)*(ymax - ymin)*(zmax - zmin);
	volumeDLP = (DLPxmax - DLPxmin)*(ymax - ymin)*(zmax - zmin);

	// Volume (number of patches) of ECM (collagen, elastin, HA) in each SLP, ILP, DLP
	double HASLP = HAFracAreaSLP*volumeSLP;
	double HAILP = HAFracAreaILP*volumeILP;
	double HADLP = HAFracAreaDLP*volumeDLP;
	cout << " HASLP, ILP, DLP: " << HASLP << " , " << HAILP << " , " << HADLP <<endl;

	sproutAgentInArea(HASLP,                        // Number of agents to sprout
			tissue,                       // patch type
			oha,                          // original collagen agent type
			SLPxmin,                      // lowest x-coordinate agent could be sprouted at
			SLPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // HA in SLP
	sproutAgentInArea(HAILP,                        // Number of agents to sprout
			tissue,                       // patch type
			oha,                          // original collagen agent type
			ILPxmin,                      // lowest x-coordinate agent could be sprouted at
			ILPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // HA in ILP
	sproutAgentInArea(HADLP,                        // Number of agents to sprout
			tissue,                       // patch type
			oha,                          // original collagen agent type
			DLPxmin,                      // lowest x-coordinate agent could be sprouted at
			DLPxmax,                      // highest x-coordinate agent
			ymin,                         // lowest y-coordinate
			ymax,                         // highest y-coordinate
			zmin,                         // lowest z-coordinate
			zmax                          // highest z-coordinate
	);    // HA in DLP
	cout << "Finished initializing HA" << endl;
#ifdef DEBUG_WORLD
	// DEBUG
	int nColl = 0;
	int nElas = 0;
	int nHA = 0;
	int nwColl = 0, nwElas = 0;
	cout << "	checking nHA ..." << endl;
	for (int i = 0; i < (nx*ny*nz - 1); i++)
	{
		this->worldECM[i].checknHA();
		nColl += this->worldECM[i].ncollagen[read_t];
		nElas += this->worldECM[i].nelastin[read_t];
        nwColl += this->worldECM[i].ncollagen[write_t];
        nwElas += this->worldECM[i].nelastin[write_t];
		nHA += this->worldECM[i].getnHA();
	}
	cout << "	finished checking nHA" << endl;
	cout << "NEW collagens: " << nColl << endl;
	cout << "NEW elastins: " << nElas << endl;
    cout << "NEW wcollagens: " << nwColl << endl;
    cout << "NEW welastins: " << nwElas << endl;
	cout << "NEW HAs: " << nHA << endl;
#endif // DEBUG_WORLD
#else
	// DEBUG
	// NOTE: Same as NetLogo initialization of HA
	int HAradius = 5;
	int HAgap = 1;
	//    int HAradius = mmToPatch(0.22); 	// 14.667
	//    int HAgap = mmToPatch(0.04); 		// 2.66
	int iy, ymin, ymax;

	/* Find patches which can be a center for sprouting original hyaluronan
	 * (initHAcenters) using a sweep from right to left (x = nx to x=0). */
	for (int iz = 0; iz < nz; iz++) {
		for (int ix = nx - 1; ix >= 0; ix--) {

			ymin = 0;
			ymax = ny - 1;
			while (ymin <= ymax ) {

				int randInt = rand()%2;

				if (randInt == 0) {
					iy = ymin++;
				} else if (randInt == 1) {
					iy = ymax--;
				}

				// If patch is tissue and qualifies as initHAcenter, sprout 6 HA on patch.  // TODO(Kim): INSERT REF?
				int in = ix + iy*nx + iz*nx*ny; // Patch row major index

				if (this->worldPatch[in].type[read_t] != tissue) continue;

				if (initHARadius(ix, iy, iz, HAradius, HAgap) == true) {

					this->initHAcenters.push_back(in);
					this->worldPatch[in].initHA = true;
					this->worldECM[in].addHAs(6);
					/*
					this->worldECM[in].HA[write_t] = this->worldECM[in].HA[read_t] + 6;

					for (int i = 0; i < 6;i++) {
#ifdef OPT_ECM
						this->worldECM[in].HAlife.push_back(100);
#else
						this->worldECM[in].HAlife[write_t].push_back(100);
						this->worldECM[in].HAlife[read_t].push_back(100);
#endif
					}
					 */
					this->worldECM[in].empty[write_t] = false;
				}
			}
		}
	}
	cout << "initHA size: " << initHAcenters.size() << endl;
#endif
}


/* Patches that can sprout original hyaluronan must	have only tissue 
 * surrounding them within a specified radius and must not have another
 * center within 2*HAradius + HAgap from them. */ // TODO(Kim): INSERT REF?
bool WHWorld::initHARadius(int a, int b, int c, double radius, double gap) {

	for (int iz = c - 2*radius - gap; iz <= c + 2*radius + gap; iz++) {
		for (int iy = b - 2*radius - gap; iy <= b + 2*radius + gap; iy++) {
			for (int ix = a - 2*radius - gap; ix <= a + 2*radius + gap; ix++) {

				// Try another neighbbor if this one is outside world boundaries
				if (iz < 0 || iz >= nz || iy < 0 || iy >= ny || ix < 0 || ix >= nx) continue;
				int in = ix + iy*nx + iz*nx*ny; // Patch row major index
				float distance = sqrt((ix - a)*(ix - a) + (iy - b)*(iy - b) + (iz -c)*(iz -c));
				float radiusgap = (2*radius + gap);

				// All neighboring patches within radius must be tissue  // TODO(Kim): INSERT REF?
				if (distance < radius && (this->worldPatch[in].type[read_t] != tissue)) {
					return false;
				}

				// There must not be any other sprouting centers within 2*radius+gap  // TODO(Kim): INSERT REF?
				if (distance < radiusgap && (this->worldPatch[in].initHA == true)) {
					return false;
				}
			}
		}
	}
	return true;
}

void WHWorld::initializeDamage() {
	// DEBUG vis
//	return;
	/*************************************
	 * BUILD WOUND                       *
	 *************************************/
	// Caclulate center, length and depth of elliptical wound
	double woundX = this->bmcx;
	double woundY = ny/2;
	double woundZ = nz/2;
	double woundDepth	= mmToPatch(this->wound[0]);
	double woundRadiusY	= mmToPatch(this->wound[1]/2);
	double woundRadiusZ	= mmToPatch(this->wound[2]/2);
	double woundSeverity    = this->wound[3];
	double randInt=0;

	double damageFraction  = woundSeverity/100.0;
#ifdef MODEL_3D
	double damZoneArea     = (4.0/3.0)*PI*woundDepth*woundRadiusY*woundRadiusZ; // #patches in damage zone
#else	// MODEL_3D
        double damZoneArea     = 0.5*PI*woundDepth*woundRadiusY;
#endif	// MODEL_3D
	double nDamagePatchesD = damZoneArea * damageFraction;
	int    nDamagePatches  = (int) nDamagePatchesD;

	int initialDamage = 0;
	this->totaldamage = 0;
	this->woundX = bmcx - 1;
	this->woundY = ny/2;
	vector <int> damzonepatches;

	/* Build reservoir so that it holds the patches in the damaged zone. Some
	 * patches will appear more than once (and some possibly not at all) in
	 * reservoir; this represents the varying level of damage at each damaged
	 * patch. */
	// Old code. No reference.
//	for (int i = 0; i < nDamagePatches; i++) {
//		initialDamage += (rand()%5 + 2);
//	}
	// New code. 'woundSeverity'% of patches in damage zone are damaged
	initialDamage = nDamagePatches;


	int *reservoir= new int [initialDamage];
	std::fill_n(reservoir, initialDamage, 0);
	cout << "initial damage: " << initialDamage << endl;
	initialDam = initialDamage;

	// Set patches within ellipse/ellipsoid centered at woundX, woundY (, woundZ) as damageZone
#ifdef MODEL_3D
	for (int iz = woundZ - woundRadiusZ; iz <= woundZ + woundRadiusZ; iz++)
#else	// MODEL_3D
        int iz = 0;
#endif	// MODEL_3D
		for (int iy= woundY - woundRadiusY; iy <= woundY + woundRadiusY; iy++){
			for (int ix= woundX - woundDepth; ix <= woundX; ix++){
				int in = ix + iy*nx + iz*nx*ny;

				// Try another neighbbor if this one is non-tissue or outside world boundaries
				if ((ix < 0) || (ix >= nx) || (iy < 0) || (iy >= ny) || (iz < 0) || (iz >= nz)) continue;
				if (worldPatch[in].type[read_t] != tissue) continue;

				float a = pow((float) (ix - woundX)/woundDepth, 2.0);
				float b = pow((float) (iy - woundY)/woundRadiusY, 2.0);
#ifdef MODEL_3D
				float c = pow((float) (iz - woundZ)/woundRadiusZ, 2.0);
#else	// MODEL_3D
                                float c = 0;
#endif	// MODEL_3D
				if (a + b + c <= 1) {
					worldPatch[in].inDamzone = true;
					damzonepatches.push_back(in);
				}
			}
		}

	for (int i = 0; i < initialDamage; i++) {
		int randompatch = rand()%damzonepatches.size();
		reservoir[i] = damzonepatches[randompatch];
	}
	printf("Wound Spec in cubes:\n\trx: %f\try:%f\trz:%f\tarea:%d\n"
		,woundDepth, woundRadiusY, woundRadiusZ, damzonepatches.size());
	cout << "Finished building wound" << endl;
	/*************************************
	 * GENERATE PLATELETS & FRAGMENT ECM *
	 *************************************/
	// Instantiate Platelet list
#ifdef MODEL_3D
	// or LARGE for acute?
#ifdef RAT_VF
	    plats = ArrayChain<Platelet*>(DEFAULT_DATA_SMALL, 4, NULL, NULL);
#else
        plats = ArrayChain<Platelet*>(DEFAULT_DATA_XLARGE, 4, NULL, NULL);
#endif
        #else	// MODEL_3D
	plats = ArrayChain<Platelet*>(DEFAULT_DATA_SMALL, 4, NULL, NULL);  //WHWorld::destroyPlat);
#endif	// MODEL_3D
	Patch* tempPatchPtr;
	Agent* tempAgentPtr;
	int plateletcount = 0;

	// Assign damage and sprout platelet at each damaged patch  // TODO(Kim): INSERT REF?
	for (int i = 0; i < initialDamage; i++) {

		int in = reservoir[i];  // Patch row major index
		int ix = worldPatch[in].indice[0];
		int iy = worldPatch[in].indice[1];
		int iz = worldPatch[in].indice[2];
		worldPatch[in].color[write_t] = cdamage;
		worldPatch[in].health[write_t] = 0;
		worldPatch[in].damage[write_t] = worldPatch[in].damage[read_t]++;

		// Sprout platelet
		if (this->worldPatch[in].isOccupied() == false) {

			tempPatchPtr = &(this->worldPatch[in]);
			Platelet* newPlatelet = new Platelet(tempPatchPtr);
			if (!this->plats.addData(newPlatelet, DEFAULT_TID)) {
				cerr << "Error: Could not add platelet in init damage()" << endl;
				exit(-1);
			}
			this->worldPatch[in].setOccupied();
			this->worldPatch[in].occupiedby[write_t] = platelet;
			plateletcount++;
		}
		// ECM fragmentation due to initial damage, within 2 unit radius of damage  // TODO(Kim): INSERT REF?

		for (int dX = -2; dX <= 2; dX++) {
			for (int dY = -2; dY <= 2; dY++) {
				for (int dZ = -2; dZ <= 2; dZ++) {
					int i = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
					// Try another patch if this one in outside world boundaries or has no ECM
					if ((ix + dX < 0) || (ix + dX >= nx) || (iy + dY < 0) || (iy + dY >= ny) || (iz + dZ < 0) || (iz + dZ >= nz)) continue;
					// DEBUG vis
					//setECM(in, m_col, 0);

					if (worldECM[i].empty[write_t] == true) continue;
					this->worldECM[i].fragmentNCollagen();
					this->worldECM[i].fragmentNElastin();
					if(!(this->worldECM[i].checknHAbool())){
						printf("before fragment\n");
						exit(-1);
					}
					this->worldECM[i].fragmentHA();
				}
			}
		}

	}
	this->totaldamage = this->countPatchType(damage);
	cout <<  this->totaldamage << " damage created" << endl;
	cout << plateletcount << " platelets sprouted " << endl;
	cout << "Finished initializing plats" << endl;
	delete[] reservoir;
}

void WHWorld::updateStats()
{
    this->totalOC = 0, this->totalNC = 0, this->totalFC = 0;
    this->totalOE = 0, this->totalNE = 0, this->totalFE = 0;
    this->totalHA = 0, this->totalFHA = 0;

    for (int i = 0; i < this->np; i++)
    {
        this->totalNC += this->worldECM[i].ncollagen[read_t];
        this->totalOC += this->worldECM[i].ocollagen[read_t];
        this->totalFC += this->worldECM[i].fcollagen[read_t];

        this->totalNE += this->worldECM[i].nelastin[read_t];
        this->totalOE += this->worldECM[i].oelastin[read_t];
        this->totalFE += this->worldECM[i].felastin[read_t];

        this->totalHA += this->worldECM[i].getnHA();
        this->totalFHA += this->worldECM[i].getnfHA();
    }
}

void WHWorld::updateCellStats()
{
	int f = 0, af = 0, m = 0, am = 0, n = 0, an = 0;

	int fibsSize = this->fibs.size();
	for (int i = 0; i < fibsSize; i++) {
		Fibroblast* fib = this->fibs.getDataAt(i);
		if (!fib) continue;
		if (fib->isAlive() == false) continue;
		if (fib->activate[read_t] == false) f++;
		else af++;
	}

	int macsSize = this->macs.size();
	for (int i = 0; i<macsSize; i++) {
		Macrophage* mac = this->macs.getDataAt(i);
		if (!mac) continue;
		if (mac->isAlive() == false) continue;
		if (mac->activate[read_t] == false) m++;
		else am++;
	}
	int neusSize = this->neus.size();
	for (int i = 0; i<neusSize; i++){
		Neutrophil* neu = this->neus.getDataAt(i);
		if (!neu) continue;
		if (neu->isAlive() == false) continue;
		if (neu->activate[read_t] == false) n++;
		else an++;
	}

	// TODO: parallelize
	this->totalCell[p_plt] = this->plats.actualSize();
	this->totalCell[p_neu]  = n;
	this->totalCell[p_mac]  = m;
	this->totalCell[p_fib]  = f;
	this->totalCell[p_anu] = an;
	this->totalCell[p_amc] = am;
	this->totalCell[p_afb] = af;
	this->totalCell[p_dam] = this->totaldamage;
}

/*
 * Each call to WHWorld::go() performs the following major steps:
 * 	(0) Cell seedings
 * 	(1) Chemical diffusion
 * 	(2) Cell function
 * 	(3) ECM function
 * 	(4) Attributes synchronization
 * 			a) Update chemicals
 * 			b) Update cells
 * 			c) Update ECM managers
 * 			d) Update patches
 */
int WHWorld::go() {

	Patch* tempPatchPtr;
	Agent* tempAgentPtr;
	double hours = this->reportHour();
	double days = this->reportDay();

//	for (int i = 0; i < this->np; i++)
//	{
//		// DEBUG vis
//		if(this->worldPatch[i].getType() == capillary) this->setECM(i, m_col, 0.5);//(ocoll/o_nCollRatio)*1000.0);
////		if(this->worldPatch[i].isInDamZone()) this->setECM(i, m_col, 0.5);//(ocoll/o_nCollRatio)*1000.0);
//	}

	// Profiling options defined in common.h
#ifdef PROFILE_MAJOR_STEPS
	struct timeval start, end;
	long elapsed_time;  // in milliseconds
#endif

	// Increment Clock in ticks (1 tick = 30 min)
	WHWorld::clock++;
	cout << "tick " << clock << " , hour " << hours << endl;


#ifdef PROFILE_MAJOR_STEPS
#if defined(GPU_DIFFUSE) && defined(_OMP)

	/* TIME_STAGE is a macro for timing a command/function and printing the
	 * timing info (See Utilities/time.h) */
#pragma omp parallel num_threads(2)
	{
		int tid = omp_get_thread_num();

		if (tid == 1) {
#ifdef OVERLAP_VIS
			cout << "go() waiting for visualization before diffuse" << endl;
			visualization::waitForVis();	
#endif	// OVERLAP_VIS
			cout << "go() task 1" << endl;
			/********************************************
			 * CHEMICAL DIFFUSION                       *
			 ********************************************/
			TIME_STAGE(this->diffuseCytokines(), "Chemical diffusion", "1");
		} else if (tid == 0) {
			cout << "go() task 2" << endl;
			/********************************************
			 * CELL SEEDING                             *
			 ********************************************/
			TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");

			/********************************************
			 * CELL FUNCTION                            *
			 ********************************************/
			TIME_STAGE(this->executeCells(), "Cells function", "2");

			/********************************************
			 * ECM FUNCTION                             *
			 ********************************************/
			TIME_STAGE(this->executeECMs(), "ECM function", "3(a)");
			// Request fragment<ECM>
			TIME_STAGE(this->requestECMfragments(), "ECM fragmentation", "3(b)");

			/********************************************
			 * ATTRIBUTES SYNCHRONIZATION               *
			 ********************************************/
			cerr << " begin update " << endl;
			TIME_STAGE(this->updateCells(), "Update cells", "4(b)");
			TIME_STAGE(this->updateECMManagers(), "Update ECM Managers", "4(c)");
			TIME_STAGE(this->updatePatches(), "Update Patches", "4(d)");
			cerr << " finished Cell, ECM and Patch updates" << endl;
		}
	}

	cerr << "begin chem update" << endl;
	TIME_STAGE(this->updateChem(), "Update chem", "4(a)");

#ifdef OVERLAP_VIS
	visualization::notifyComDone();	
#endif	// OVERLAP_VIS

#else	// GPU_DIFFUSE && _OMP
	/* TIME_STAGE is a macro for timing a command/function and printing the
	 * timing info (See Utilities/time.h) */

	/********************************************
	 * CHEMICAL DIFFUSION                       *
	 ********************************************/
	TIME_STAGE(this->diffuseCytokines(), "Chemical diffusion", "1");

	/********************************************
	 * CELL SEEDING                             *
	 ********************************************/
	TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");


	/********************************************
	 * CELL FUNCTION                            *
	 ********************************************/
	TIME_STAGE(this->executeCells(), "Cells function", "2");

	/********************************************
	 * ECM FUNCTION                             *
	 ********************************************/
	TIME_STAGE(this->executeECMs(), "ECM function", "3(a)");
	// Request fragment<ECM>
	TIME_STAGE(this->requestECMfragments(), "ECM fragmentation", "3(b)");

	/********************************************
	 * ATTRIBUTES SYNCHRONIZATION               *
	 ********************************************/
	cerr << " begin update " << endl;

	TIME_STAGE(this->updateCells(), "Update cells", "4(b)");
	TIME_STAGE(this->updateECMManagers(), "Update ECM Managers", "4(c)");
	TIME_STAGE(this->updatePatches(), "Update Patches", "4(d)");
	TIME_STAGE(this->updateChem(), "Update chem", "4(a)");
#endif	// GPU_DIFFUSE && _OMP
#else	// PROFILE_MAJOR_STEPS


	/********************************************
	 * CHEMICAL DIFFUSION                       *
	 ********************************************/
#ifdef GPU_DIFFUSE
	this->diffuseCytokines();
#endif

	/********************************************
	 * CELL SEEDING                             *
	 ********************************************/
	this->seedCells(hours);

	/********************************************
	 * CELL FUNCTION                            *
	 ********************************************/
	this->executeCells();

	/********************************************
	 * ECM FUNCTION                             *
	 ********************************************/
	this->executeECMs();
	this->requestECMfragments();

	/********************************************
	 * ATTRIBUTES SYNCHRONIZATION               *
	 ********************************************/
	cerr << " begin update " << endl;

#ifndef GPU_DIFFUSE
	this->diffuseCytokines();
#endif

	this->updateCells();
	this->updateECMManagers();
	this->updatePatches();
	cerr << " finished Cell, ECM and Patch updates" << endl;
	this->updateChem();
	cerr << "begin chem update" << endl;

#endif

	// Update totaldamage to be used in the next iteration
	this->totaldamage = this->countPatchType(damage);

#ifdef PROFILE_ECM
	cout << "ECM Repair took: 		" << ECMrepairTime << " us" << endl;
	cout << "ECM HA life decrements took:	" << HAlifeTime << " us" << endl;
	cout << "ECM danger signal took:	" << ECMdangerTime << " us" << endl;
	cout << "ECM scar formation took:	" << ECMscarTime << " us" << endl;

	this->ECMrepairTime = 0;
	this->HAlifeTime    = 0;
	this->ECMdangerTime = 0;
	this->ECMscarTime   = 0;
#endif
#ifdef DEBUG_WORLD
	cout << "Max HA size: " << maxHAsize << endl;
	maxHAsize = 0;

	cout << "tick " << clock << " dead neus " << deadneus << endl;
	cout << "tick " << clock << " dead afibs " << dead_afibs << endl;
	cout << "tick " << clock << " new fibs " << newfibs << endl;
	cout << "tick " << clock << " deac fibs " << deactfibs << endl;
	cout << "tick " << clock << " acti fibs " << actfibs << endl;
	cout << "tick " << clock << " dead fibs " << deadfibs << endl;
	deadfibs = 0;
	deadneus = 0;
	newfibs = 0;
	actfibs = 0;
	dead_afibs = 0;
	deactfibs = 0;

	cout << "HA fragments: +" << HAfrags << endl;
	HAfrags = 0;

	cout << "New collagen monomer: +" << newOcolls << endl;
	newOcolls = 0;
#endif


	this->updateCellStats();


	return 0;
}

/****************************************************************
 * MAJOR SECTION SUBROUTINES - begin                            *
 ****************************************************************/

void WHWorld::seedCells(float hours) {

	// DEBUG rat
//	return;
#ifdef NO_SECRETION
  return;
#endif

	int impactStress = -1;
	/**********************************************
	 * Sprout Agents based on impact stress	*
	 **********************************************/
	if ((hours <= 6.5) && (hours >= 2.5)) {  //TODO(Kim): INSERT REF? (Why 2.5-6.5?)

		if (this->treatmentOption == resonantvoice) impactStress = RVIS;
		if (this->treatmentOption == spontaneousspeech) impactStress = SSIS;

		if (impactStress >= 0) {  //TODO(Kim): INSERT REF? (Why not for voice rest?)
#ifndef CALIBRATION
			int platSproutingFactor = 20; // TODO(Kim): INSERT REF? (Why *20?)
			int neuSproutingFactor  = 1;
#else	// !CALIBRATION
			int platSproutingFactor = WHWorld::sproutingAmount[0];
			int neuSproutingFactor  = WHWorld::sproutingAmount[1];
#endif  // !CALIBRATION
			cout << " sprout plat and neut due to impact stress " << endl;
			sproutAgentInArea(impactStress*platSproutingFactor,
					tissue,  // TODO(Kim): INSERT REF? (Why tissue?)
					platelet,
					bmcx - 1,  // TODO(Kim): INSERT REF? (Why bmcx?)
					bmcx,  // TODO(Kim): INSERT REF? (Why bmcx?)
					bmcy - 1,  // TODO(Kim): INSERT REF? (Why bmcy?)
					bmcy,  // TODO(Kim): INSERT REF? (Why bmcy?)
					0,  // lowest z-coordinate agent could be sprouted at
					nz);
			sproutAgent(impactStress*neuSproutingFactor,  // TODO(Kim): INSERT REF? (Why not *20?)
					capillary,    // TODO(Kim): INSERT REF? (Why capillary?)
					neutrophil,
					0,  // lowest x-coordinate agent could be sprouted at
					nx,
					0,  // lowest y-coordinate agent could be sprouted at
					ny,
					0,  // lowest z-coordinate agent could be sprouted at
					nz);
		}
	}

	/**********************************************
	 * Sprout Agents based on total damage	*
	 **********************************************/
	int totaldamage = this->totaldamage;//this->countPatchType(damage);
	if (hours == 0) {
		// NS: Do something?
	} else if (totaldamage == 0) {  // TODO(Kim): INSERT REF? (Why only if no damage & hours==0?)
		cout << " no damage, sprout neutrophil and macrophage " << endl;
#ifndef CALIBRATION
		float neuSproutingFreq = 2.0f;
		float macSproutingFreq = 4.0f;

		int neuSproutingAmount = 8; // TODO(Kim): INSERT REF? (Why 8?)
		int macSproutingAmount = 1; // TODO(Kim): INSERT REF? (Why 1?)
#else	// !CALIBRATION
		float neuSproutingFreq = WHWorld::sproutingFrequency[0];
		float macSproutingFreq = WHWorld::sproutingFrequency[1];

		int neuSproutingAmount = WHWorld::sproutingAmount[2];
		int macSproutingAmount = WHWorld::sproutingAmount[3];
#endif	// !CALIBRATION
		if (fmod (hours, neuSproutingFreq) == 0) {  // Sprout neutrophil every 2 hours
			sproutAgent (neuSproutingAmount,  // Number of agents to sprout
					capillary,  // TODO(Kim): INSERT REF? (Why capillary?)
					neutrophil,
					0,  // lowest x-coordinate agent could be sprouted at
					nx,
					0,  // lowest y-coordinate agent could be sprouted at
					ny,
					0,  // lowest z-coordinate agent could be sprouted at
					nz);
		}
		if (fmod(hours, macSproutingFreq) == 0) {  // Sprout macrophage every 4 hours
			sproutAgent (macSproutingAmount,  // Number of agents to sprout
					capillary,  // TODO(Kim): INSERT REF? (Why capillary?)
					macrophag,
					0,  // lowest x-coordinate agent could be sprouted at
					nx,
					0,  // lowest y-coordinate agent could be sprouted at
					ny,
					0,  // lowest z-coordinate agent could be sprouted at
					nz,
					blood);
		}
	} else {	// hours != 0 && totaldamage != 0
		cout << " damage " << endl;
#ifndef CALIBRATION
		cout << "no calib" << endl;
		float neuSproutingFreq  =  2.0f;
		float bmacSproutingFreq =  4.0f;
		float tmacSproutingFreq =  6.0f;
		float fibSproutingFreq  = 12.0f;

		int neuSproutingAmount   = 8 + totaldamage*0.00001;//0.01; //TODO(Kim): INSERT REF?(Why 8+totaldamage*0.01?)
		int bmacSproutingAmount  = 1 + totaldamage*0.00001;//0.01; //TODO(Kim): INSERT REF?(Why 1+totaldamage*0.01?)
		int tmacSproutingAmount1 = 1 + totaldamage*0.0001;//0.1; //TODO(Kim): INSERT REF?(Why 1+totaldamage*0.01?)
		int tmacSproutingAmount2 = 1; 		    //TODO(Kim): INSERT REF?(Why 1?)
		int fibSproutingAmount1  = 1 + totaldamage*0.0001;//0.01;
		int fibSproutingAmount2  = 1;
#else	// !CALIBRATION
		float neuSproutingFreq  = WHWorld::sproutingFrequency[2];
		float bmacSproutingFreq = WHWorld::sproutingFrequency[3];
		float tmacSproutingFreq = WHWorld::sproutingFrequency[4];
		float fibSproutingFreq  = WHWorld::sproutingFrequency[5];

		int neuSproutingAmount   = WHWorld::sproutingAmount[4] + totaldamage*WHWorld::sproutingAmount[5];
		int bmacSproutingAmount  = WHWorld::sproutingAmount[6] + totaldamage*WHWorld::sproutingAmount[7];
		int tmacSproutingAmount1 = WHWorld::sproutingAmount[8] + totaldamage*WHWorld::sproutingAmount[9];
		int tmacSproutingAmount2 = WHWorld::sproutingAmount[10];
		int fibSproutingAmount1 = WHWorld::sproutingAmount[11] + totaldamage*WHWorld::sproutingAmount[12];
		int fibSproutingAmount2 = WHWorld::sproutingAmount[13];
#endif	// !CALIBRATION
		if (fmod(hours, neuSproutingFreq) == 0) {  // Sprout neutrophil every 2 hours
			cout << " every " << neuSproutingFreq   <<" hours, sprout "
					<< neuSproutingAmount << " neutrophil" << endl;
			sproutAgent(neuSproutingAmount,  // Number of agents to sprout
					capillary,  // TODO(Kim): INSERT REF? (Why capillary?)
					neutrophil,
					0,  // lowest x-coordinate agent could be sprouted at
					nx,
					0,  // lowest y-coordinate agent could be sprouted at
					ny,
					0,  // lowest z-coordinate agent could be sprouted at
					nz);
		}
		if (fmod(hours, bmacSproutingFreq) == 0) {  // Sprout capillary macrophage every 4 hours
			cout << " every " << bmacSproutingFreq << " hours, sprout "
					<< bmacSproutingAmount << " bmacrophages" << endl;
			sproutAgent(bmacSproutingAmount,  // Number of agents to sprout
					capillary,  // TODO(Kim): INSERT REF? (Why capillary?)
					macrophag,
					0,  // lowest x-coordinate agent could be sprouted at
					nx,
					0,  // lowest y-coordinate agent could be sprouted at
					ny,
					0,  // lowest z-coordinate agent could be sprouted at
					nz,
					blood);
		}
		if (fmod(hours, tmacSproutingFreq) == 0) {  // Sprout tissue macrophages every 6 hours
			cout << " every " << tmacSproutingFreq << " hours, sprout "
					<< tmacSproutingAmount1 * tmacSproutingAmount2 * 2 << " macrophages" << endl;
//			for (int i = 0; i < tmacSproutingAmount1; i++) {
//				int iz = rand() % nz;
				sproutAgentInArea(tmacSproutingAmount1 * tmacSproutingAmount2,
						tissue,  // TODO(Kim): INSERT REF? (Why tissue?)
						macrophag,
						0,  // lowest x-coordinate agent could be sprouted at
						bmcx,  // TODO(Kim): INSERT REF? (Why bmcx?)
						0,  // lowest y-coordinate agent could be sprouted at
						1,  // highest y-coordinate agent could be sprouted at, top of the world
						0,  // lowest z-coordinate agent could be sprouted at
						nz - 1,
						tissue);
				sproutAgentInArea(tmacSproutingAmount1 * tmacSproutingAmount2,
						tissue,  // TODO(Kim): INSERT REF? (Why tissue?)
						macrophag,
						0,  // lowest x-coordinate agent could be sprouted at
						bmcx,  // TODO(Kim): INSERT REF? (Why bmcx?)
						ny - 1,  // TODO(Kim): INSERT REF? (Why ny-1?)
						ny,
						0,  // lowest z-coordinate agent could be sprouted at
						nz - 1,
						tissue);
//			}
		}
		if (fmod(hours, fibSproutingFreq) == 0) {  // Sprout fibroblasts every 12 hours
			cout << " every " << fibSproutingFreq << " hours, sprout "
					<< fibSproutingAmount1 * fibSproutingAmount2 * 2 << " fibroblasts" << endl;
//			for (int i = 0; i < fibSproutingAmount1; i++) {
//				int iz = rand() % nz;
				sproutAgentInArea(fibSproutingAmount1 * fibSproutingAmount2,
						tissue,  // TODO(Kim): INSERT REF? (Why tissue?)
						fibroblast,
						0,  // lowest x-coordinate agent could be sprouted at
						bmcx,  // TODO(Kim): INSERT REF? (Why bmcx?)
						0,  // lowest y-coordinate agent could be sprouted at
						1,  // highest y-coordinate agent could be sprouted at, top of the world
						0,  // lowest z-coordinate agent could be sprouted at
						nz - 1);
				sproutAgentInArea(fibSproutingAmount1 * fibSproutingAmount2,
						tissue,  // TODO(Kim): INSERT REF? (Why tissue?)
						fibroblast,
						0,  // lowest x-coordinate agent could be sprouted at
						bmcx,  // TODO(Kim): INSERT REF? (Why bmcx?)
						ny - 1,  // TODO(Kim): INSERT REF? (Why ny-1?), bottom of the world
						ny,
						0,  // lowest z-coordinate agent could be sprouted at
						nz - 1);
//			}
		}
	}
}

void WHWorld::diffuseCytokines(float dt) {
#ifdef PDE_DIFFUSE

#ifdef GPU_DIFFUSE
	this->diffuserPtr->diffuseChemGPU(this->bmcx-1);
#else
	/* To satisfy stability conditions, at given dx, dt repeat
	 * central approximation finite difference diffusion until reach 30 min tick*/
	for (int i = 0; i < 30/dt; i++) this->diffuserPtr->diffuseChemCPU(dt);
#endif

#else
	this->NetlogoDiffuse();
#endif
}

void WHWorld::executePlats() {
	int platsSize = plats.size(); /* This is only an upper bound on cell list
	 * size. It is NOT an actual count of cells (some entries are NULL) */
	//#pragma omp parallel for
	for (int i = 0; i < platsSize; i++) {
		Platelet* plat = plats.getDataAt(i);
		if (!plat) continue;
		plat->cellFunction();
	}
}

void WHWorld::executeNeus() {
	int neusSize = neus.size(); /* This is only an upper bound on cell list
	 * size. It is NOT an actual count of cells (some entries are NULL) */
#pragma omp parallel for schedule(dynamic, 1000)
	for (int i = 0; i < neusSize; i++) {
		Neutrophil* neu = neus.getDataAt(i);
		if (!neu) continue;
		neu->cellFunction();
	}
}

void WHWorld::executeMacs() {
	int macsSize = macs.size(); /* This is only an upper bound on cell list
	 * size. It is NOT an actual count of cells (some entries are NULL) */
#pragma omp parallel for schedule(dynamic, 1000)
	for (int i = 0; i < macsSize; i++) {
		Macrophage* mac = macs.getDataAt(i);
		if (!mac) continue;
		mac->cellFunction();
	}
}

void WHWorld::executeFibs() {
	int fibsSize = fibs.size(); /* This is only an upper bound on cell list
	 * size. It is NOT an actual count of cells (some entries are NULL) */
        
#pragma omp parallel for schedule(dynamic, 1000)
	for (int i = 0; i < fibsSize; i++) {
		Fibroblast* fib = fibs.getDataAt(i);
		if (!fib) continue;
		fib->cellFunction();
	}
}

void WHWorld::executeCells() {
	// DEBUG rat
//	return;

#ifdef PROFILE_CELL_FUNC
	TIME_STAGE(this->executePlats(),	"	Cell Function: Platelets",		"	");
	TIME_STAGE(this->executeNeus(),		"	Cell Function: Neutrophils",	"	");
	TIME_STAGE(this->executeMacs(),		"	Cell Function: Macrophages",	"	");
	TIME_STAGE(this->executeFibs(),		"	Cell Function: Fibroblasts",	"	");
#else
	cerr << "executing plats ..." << endl;
	this->executePlats();
	cerr << "executing neus ..." << endl;
	this->executeNeus();
	cerr << "executing macs ..." << endl;
	this->executeMacs();
	cerr << "executing fibs ..." << endl;
	this->executeFibs();
	cerr << "cell execution completed ..." << endl;
#endif
}

void WHWorld::executeECMs(){
	cerr << " ECM function  " << endl;
	int numPatches = (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny;
	int numNonEmpty[MAX_NUM_THREADS] = {0};
#pragma omp parallel
	{
		int sum = 0;
#pragma omp for schedule(dynamic, 1000)
		for (int in = 0; in < numPatches; in++) {
			if (worldECM[in].empty[read_t] == false) {
				this->worldECM[in].ECMFunction();
				sum++;
			}
		}
#ifdef _OMP
		numNonEmpty[omp_get_thread_num()] = sum;
#else
		numNonEmpty[0] = sum;
#endif
	}


}

void WHWorld::requestECMfragments() {

#ifndef CALIBRATION
	float TNFthreshold  = 10.0;//50.0;//10.0;
	float MMP8threshold = 10.0;//50.0;//10.0;
#else	// CALIBRATION
	float TNFthreshold  = WHWorld::thresholdTNFdamage;
	float MMP8threshold = WHWorld::thresholdMMP8damage;
#endif	// CALIBRATION

	if (WHWorld::highTNFdamage == true) {  //TODO(Kim): INSERT REF?
		cout << " high TNF damage " << endl;
		WHWorld::highTNFdamage = false;
#pragma omp parallel firstprivate(TNFthreshold, MMP8threshold)
	{
// DEBUG
//  printf("[%d] -- 1\n", omp_get_thread_num());
#pragma omp for
		for (int in = 0; in < np; in++) {
			float patchTNF  = this->WHWorldChem->getPchem(TNF,  in);
			float patchMMP8 = this->WHWorldChem->getPchem(MMP8, in);	// if OPT_CHEM we want to reuse
			// multiple chems on each patch
			// as much as possible
			if (patchTNF > TNFthreshold)
			{
				this->worldECM[in].fragmentNCollagen();
				this->worldECM[in].fragmentNElastin();
				this->worldECM[in].fragmentHA();
			}
		}
// DEBUG
//  printf("[%d] -- 2\n", omp_get_thread_num());
	}
	}


	// This always get executed if OPT_CHEM is NOT turned on.
	// Otherwise, only executed if highTNFdamage is not set.
	cout << " high MMP8 damage " << endl;
#pragma omp parallel for firstprivate(MMP8threshold)
	for (int in = 0; in < np; in++) {
		float patchMMP8 = this->WHWorldChem->getPchem(MMP8, in);
		if (patchMMP8 > MMP8threshold)
		{
			this->worldECM[in].fragmentNCollagen();
		}
	}
	/*
#ifdef OPT_CHEM
  }
#endif	// OPT_CHEM
	 */
}

/* Each patch diffuses 50% of its chemical equally to its 8 neighboring
 * patches. (Each neighbor receives 1/8 of 50% of the patch's original amount
 * of chemical neighboring patch. */
void WHWorld::NetlogoDiffuse() {
	// NetLogoDiffuse is not obsolete	-- 10/10/2016
	return;
}



/*  TNF = 0,
  TGF = 1,
  FGF = 2,
  MMP8 = 3,
  IL1beta = 4,
  IL6 = 5,
  IL8 = 6,
  IL10 = 7,*/
void WHWorld::printChemInfo(){
	char cellname[4][6]  = {"plat", "neu", "mac", "fib"};
	char chemname[8][10] = {"TNF", "TGF", "FGF", "MMP8", "IL1beta", "IL6", "IL8", "IL10"};

	for (int icell = 0; icell < 4; icell++) {
		for (int ichem = 0; ichem < 8; ichem++) {
			printf("\t%s,\t%s,\t%f\t(%d,%d,%d)\n", cellname[icell], chemname[ichem],
					chemSecreted[icell][ichem],
					chemSecretedCoord[icell][ichem][0],
					chemSecretedCoord[icell][ichem][1],
					chemSecretedCoord[icell][ichem][2]);
		}
	}
	printf("\n");
	//
	//  for (int ichem = 0; ichem < 8; ichem++) {
	//      printf("\t%s,\tp:\t%f,\td:\t%f,\tt:\t%f\t(%d,%d,%d)\n", chemname[ichem],
	//          maxPatchChem[ichem][0],
	//          maxPatchChem[ichem][1],
	//          maxPatchChem[ichem][2],
	//          maxPatchChemCoord[ichem][0],
	//          maxPatchChemCoord[ichem][1],
	//          maxPatchChemCoord[ichem][2]);
	//  }
	//
	//  printf("\n");
	//  for (int ichem = 0; ichem < 8; ichem++) {
	//      printf("\t%s,\toldpmax:\t%f\t(%d,%d,%d)\n", chemname[ichem],
	//          maxOldPatchChem[ichem],
	//          maxOldPatchChemCoord[ichem][0],
	//          maxOldPatchChemCoord[ichem][1],
	//          maxOldPatchChemCoord[ichem][2]);
	//  }
	//
	//  printf("\n");
	//  for (int ichem = 0; ichem < 8; ichem++) {
	//      printf("\t%s,\toldpmin:\t%f\t(%d,%d,%d)\n", chemname[ichem],
	//          minOldPatchChem[ichem],
	//          minOldPatchChemCoord[ichem][0],
	//          minOldPatchChemCoord[ichem][1],
	//          minOldPatchChemCoord[ichem][2]);
	//  }

}
void WHWorld::updateChem() {

// HERE (2) MOVE BELOW -- WHWorldChem->update(totaldamage);

	int totaldam = this->totaldamage;//countPatchType(damage);

	// DEBUG rat
#ifdef CALIBRATION
	float cytokineDecay[8];
	for (int ic = 0; ic < 8; ic++)
	{
		cytokineDecay[ic] = WHWorld::cytokineDecay[ic];
	}
#else		// CALIBRATION

	float cytokineDecay[8] = {	0.319,		// TNF
					0.00003052,	// TGF
					0.67,		// FGF
					0.2,
					0.2,
					0.7937,		// IL6
					0.9284,		// IL8
					0.2};

#endif		// CALIBRATION

#ifdef GPU_DIFFUSE
	WHWorldChem->updateChemGPU(totaldamage, cytokineDecay);
#else	// GPU_DIFFUSE
	WHWorldChem->updateChemCPU(totaldamage, cytokineDecay);
#endif	// GPU_DIFFUSE

#ifdef DEBUG_WORLD
	// DEBUG
	printChemInfo();

	// DEBUG
	for (int icell = 0; icell < 4; icell++) {
		for (int ichem =0; ichem < 8; ichem++) {
			chemSecreted[icell][ichem] = 0.0;
			chemSecretedCoord[icell][ichem][0] = -1;
			chemSecretedCoord[icell][ichem][1] = -1;
			chemSecretedCoord[icell][ichem][2] = -1;
		}
	}

	//	for (int ichem = 0; ichem < 8; ichem++) {
	//	    maxPatchChem[ichem][0] = 0.0;
	//	    maxPatchChem[ichem][1] = 0.0;
	//	    maxPatchChem[ichem][2] = 0.0;
	//
	//	    maxPatchChemCoord[ichem][0] = -1;
	//	    maxPatchChemCoord[ichem][1] = -1;
	//	    maxPatchChemCoord[ichem][2] = -1;
	//
	//	    maxOldPatchChem[ichem] = 0.0;
	//
	//	    maxOldPatchChemCoord[ichem][0] = -1;
	//	    maxOldPatchChemCoord[ichem][1] = -1;
	//	    maxOldPatchChemCoord[ichem][2] = -1;
	//
	//	    minOldPatchChem[ichem] = 90000.0;
	//
	//	    minOldPatchChemCoord[ichem][0] = -1;
	//	    minOldPatchChemCoord[ichem][1] = -1;
	//	    minOldPatchChemCoord[ichem][2] = -1;
	//	}
#endif	// DEBUG_WORLD

}

void WHWorld::executeAllECMUpdates() {

    int np = this->np;
#pragma omp parallel
{

//#pragma omp for nowait
#pragma omp for schedule(dynamic, 1000)
    for (int in = 0; in < np; in++)
    {
        this->worldECM[in].updateECM();
    }

//#pragma omp for schedule(guided, 1000)
/*#pragma omp for schedule(dynamic, 1000)
    for (int in = 0; in < np; in++)
    {
        this->worldPatch[in].updatePatch();
    }
*/
}
/*
#pragma omp parallel for  
	for (int iz = 0; iz < nz; iz++) {
		//#pragma omp parallel for
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {

				int in = ix + iy*nx + iz*nx*ny;
				this->worldECM[in].updateECM();
			}
		}       
	}               
*/
}                       

void WHWorld::executeAllECMResetRequests() {
#ifndef OPT_FECM_REQUEST
#pragma omp parallel for
	for (int iz = 0; iz < nz; iz++) {
		//#pragma omp parallel for
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				int in = ix + iy*nx + iz*nx*ny;
				this->worldECM[in].resetrequests();
			}
		}       
	}               
#endif
}                       

void WHWorld::updateECMManagers() {
#ifdef PROFILE_ECM_UPDATE
	TIME_STAGE(this->executeAllECMUpdates(),		"	updateECM()",		"	");
	TIME_STAGE(this->executeAllECMResetRequests(),	"	resetrequests()",	"	");
#else                   
	this->executeAllECMUpdates();
	this->executeAllECMResetRequests();
#endif                  
}                       

void WHWorld::updatePatches() {
    int np = this->np;
#pragma omp parallel for schedule(dynamic, 1000)
    for (int in = 0; in < np; in++)
    {
        this->worldPatch[in].updatePatch();
    }
}                       

void WHWorld::updatePlats() {
	cerr << "	removing dead plats" << endl;
	int platsSize = plats.size();
	// Currently number of platelets are not big enough to benefit from parallelization
	//#pragma omp parallel for
	for (int i = 0; i < platsSize; i++) {
		// Get pointer of platelet i from the array chain
		Platelet* plat = plats.getDataAt(i);
		if (!plat) continue;  // Platelet was deleted
		plat->updateAgent();
		// Remove dead platelets
		if (plat->isAlive() == false) {
			// Get residing patch index and update its occupancy
			int in = plat->getIndex();
			this->worldPatch[in].clearOccupied();
			this->worldPatch[in].occupiedby[write_t] = nothing;
			plats.deleteData(i, DEFAULT_TID);
			delete plat;
		}
	}
	Platelet::numOfPlatelets= plats.actualSize();
}

void WHWorld::updateNeus() {
	cerr << "	removing dead neus" << endl;
	int neusSize = neus.size();
#pragma omp parallel for
	for (int i = 0; i < neusSize; i++) {

#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		// Get pointer of neutrophil i from the array chain
		Neutrophil* neu = neus.getDataAt(i);
		if (!neu) continue;  // Neutrophil was deleted
		neu->updateAgent();
		// Remove dead neutrophils
		if (neu->isAlive() == false) {
			// Get residing patch index and update its occupancy
			int in = neu->getIndex();
			this->worldPatch[in].clearOccupied();
			this->worldPatch[in].occupiedby[write_t] = nothing;
			neus.deleteData(i, tid);
			delete neu;
		}
	}
	Neutrophil::numOfNeutrophil = neus.actualSize();
}

void WHWorld::updateMacs() {
	cerr << "	removing dead macs" << endl;
	int macsSize = macs.size();
#pragma omp parallel for
	for (int i = 0; i < macsSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		// Get pointer of macrophage i from the array chain
		Macrophage* mac = macs.getDataAt(i);
		if (!mac) continue;  // Macrophage was deleted
		mac->updateAgent();
		// Remove dead macrophages
		if (mac->isAlive() == false) {
			// Get residing patch index and update its occupancy
			int in = mac->getIndex();
			this->worldPatch[in].clearOccupied();
			this->worldPatch[in].occupiedby[write_t] = nothing;
			macs.deleteData(i, tid);
			delete mac;
		}
	}
	Macrophage::numOfMacrophage = macs.actualSize();
}

void WHWorld::updateFibs() {
	cerr << "	removing dead fibs" << endl;
	int fibsSize = fibs.size();
#pragma omp parallel for
	for (int i = 0 ; i < fibsSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		// Get pointer of fibroblast i from the array chain
		Fibroblast* fib = fibs.getDataAt(i);
		if (!fib) continue;  // Fibroblast was deleted
		fib->updateAgent();
		// Remove dead fibroblasts
		if (fib->isAlive() == false) {
			// Get residing patch index and update its occupancy
			int in = fib->getIndex();
			this->worldPatch[in].clearOccupied();
			this->worldPatch[in].occupiedby[write_t] = nothing;
			fibs.deleteData(i, tid);
			delete fib;
		}
	}
	Fibroblast::numOfFibroblasts = fibs.actualSize();
}

/*
 * Steps:
 * 1. Perform updates
 * 2. Remove all dead cells
 * 3. If OMP, add cells from thread-local lists to corresponding global lists
 */
void WHWorld::updateCells() {
	//	this->updateFibs();
	//	this->updateMacs();
	//	this->updateNeus();
	//	this->updatePlats();

	// Add new cells
	// NS: Only fibroblast part is necessary?
#ifdef _OMP
	/* In OMP, cells were only added to each thread's local list when
	 * sproutAgent() was called. Thus, this step is needed to add those cells
	 * onto the global lists.	 */
	cerr << "	updateCell() _OMP" << endl;
	// TODO: parallelize
	//	int numThreads = omp_get_num_threads();
	int numThreads = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
	cout << "		numThreads = " << numThreads << endl;
	for (int tid = 0; tid < numThreads; tid++) {
		// Platelets
		vector<Platelet*>* pvec_ptr = localNewPlats[tid];
		for (vector<Platelet*>::iterator plat_it = pvec_ptr->begin(); plat_it != pvec_ptr->end(); plat_it++) {
			Platelet* newPlat = *plat_it;
			if (!plats.addData(newPlat, tid)) {
				cerr << "Error: Could not add platelet" << endl;
				exit(-1);
			}
		}
		pvec_ptr->clear();

		// Neutrophils
		vector<Neutrophil*>* nvec_ptr = localNewNeus[tid];
		for (vector<Neutrophil*>::iterator neu_it = nvec_ptr->begin(); neu_it != nvec_ptr->end(); neu_it++) {
			Neutrophil* newNeu = *neu_it;
			if (!neus.addData(newNeu, tid)) {
				cerr << "Error: Could not add neutrophil" << endl;
				exit(-1);
			}
		}
		nvec_ptr->clear();

		// Macrophages
		vector<Macrophage*>* mvec_ptr = localNewMacs[tid];
		for (vector<Macrophage*>::iterator mac_it = mvec_ptr->begin(); mac_it != mvec_ptr->end(); mac_it++) {
			Macrophage* newMac = *mac_it;
			if (!macs.addData(newMac, tid)) {
				cerr << "Error: Could not add macrophage" << endl;
				exit(-1);
			}
		}
		mvec_ptr->clear();

		// Fibroblasts
		vector<Fibroblast*>* fvec_ptr = localNewFibs[tid];
		for (vector<Fibroblast*>::iterator fib_it = fvec_ptr->begin(); fib_it != fvec_ptr->end(); fib_it++) {
			Fibroblast* newFib = *fib_it;
			if (!fibs.addData(newFib, tid)) {
				cerr << "Error: Could not add fibroblast" << endl;
				exit(-1);
			}
		}
		fvec_ptr->clear();
	}
#endif
	this->updateFibs();
	this->updateMacs();
	this->updateNeus();
	this->updatePlats();

}

/****************************************************************
 * MAJOR SECTION SUBROUTINES - end                              *
 ****************************************************************/

//#ifdef VISUALIZATION


void WHWorld::incECM(int index, ecm_i ecmType, float count)
{
	(this->ecmMap[ecmType][index])          += count/((float) (this->maxFCol));
	(this->ecmGradMap[ecmType][index])      += count/((float) (this->maxChangeFCol));
	(this->ecmPreProcMap[ecmType][index])    = RAW_ECM_FACTOR*ecmMap[ecmType][index] + RAW_dECM_FACTOR*ecmGradMap[ecmType][index];
}

void WHWorld::decECM(int index, ecm_i ecmType, float count)
{
	(this->ecmMap[ecmType][index])        -= count/((float) (this->maxFCol));
	(this->ecmGradMap[ecmType][index])    -= count/((float) (this->maxChangeFCol));
	(this->ecmPreProcMap[ecmType][index])  = RAW_ECM_FACTOR*ecmMap[ecmType][index] + RAW_dECM_FACTOR*ecmGradMap[ecmType][index];
}

void WHWorld::setECM(int index, ecm_i ecmType, float count)
{
	(this->ecmMap[ecmType][index]) = count/((float) (this->maxFCol));
}

void WHWorld::resetECMmap()
{
	std::memset(this->ecmMap[m_col],        0, np * sizeof(float));
	std::memset(this->ecmGradMap[m_col],    0, np * sizeof(float));
	std::memset(this->ecmPreProcMap[m_col], 0, np * sizeof(float));
}

//#endif		// VISUALIZATION


void WHWorld::getWndPos(
               int &wnd_xb,
               int &wnd_xe,
               int &wnd_yb,
               int &wnd_ye,
	       int &wnd_zb,
               int &wnd_ze)
{
  wnd_xb = this->wnd_xb;
  wnd_xe = this->wnd_xe;
  wnd_yb = this->wnd_yb;
  wnd_ye = this->wnd_ye;
  wnd_zb = this->wnd_zb;
  wnd_ze = this->wnd_ze;
}

//NOTE: only use this function to sprout nc, ne in initialization.
void WHWorld::sproutAgent(int num, int patchType, agent_t agentType,
		int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int bloodORtiss) {

#ifdef OPT_CELL_SEEDING
	if (xmin != 0 || xmax != nx || ymin != 0 || ymax != ny || zmin != 0 || zmax != nz){
		sproutAgentInArea (num, patchType, agentType, xmin, xmax, ymin, ymax, zmin, zmax, bloodORtiss);
	} else {
		sproutAgentInWorld (num, patchType, agentType);
	}

#else
	// Target a specific area of the world
	sproutAgentInArea (num, patchType, agentType, xmin, xmax, ymin, ymax, zmin, zmax, bloodORtiss);
#endif
}


void WHWorld::createCell(int cellType, bool bloodORtiss, Patch* tempPatchPtr)
{
	// Create cell and bind to patch
	switch (cellType) {

	case platelet: {
		Platelet* newPlatelet = new Platelet(tempPatchPtr);
#ifdef _OMP
		int tid = omp_get_thread_num();
		this->localNewPlats[tid]->push_back(newPlatelet);
#else
		if (!this->plats.addData(newPlatelet, DEFAULT_TID)) {
			cerr << "Error: Could not add platelet in sproutAgent()" << endl;
			exit(-1);
		}
#endif
		break;
	}
	case fibroblast: {
		Fibroblast* newFib = new Fibroblast(tempPatchPtr);
#ifdef _OMP
		int tid = omp_get_thread_num();
		this->localNewFibs[tid]->push_back(newFib);
#else
		if (!this->fibs.addData(newFib, DEFAULT_TID)) {
			cerr << "Error: Could not add fibroblast in sproutAgent()" << endl;
			exit(-1);
		}
#endif
		break;
	}
	case afibroblast: {
		Fibroblast* newFib = new Fibroblast(tempPatchPtr);
		newFib->fibActivation();
#ifdef _OMP
		int tid = omp_get_thread_num();
		this->localNewFibs[tid]->push_back(newFib);
#else
		if (!this->fibs.addData(newFib, DEFAULT_TID)) {
			cerr << "Error: Could not add fibroblast in sproutAgent()" << endl;
			exit(-1);
		}
#endif
		//        newFib->fibActivation();
		break;
	}
	case macrophag: {
		Macrophage* newMac = new Macrophage(tempPatchPtr, bloodORtiss);
#ifdef _OMP
		int tid = omp_get_thread_num();
		this->localNewMacs[tid]->push_back(newMac);
#else
		if (!this->macs.addData(newMac, DEFAULT_TID)) {
			cerr << "Error: Could not add macrophage in sproutAgent()" << endl;
			exit(-1);
		}
#endif
		break;
	}
	case neutrophil: {
		Neutrophil* newNeu = new Neutrophil(tempPatchPtr);
#ifdef _OMP
		int tid = omp_get_thread_num();
		this->localNewNeus[tid]->push_back(newNeu);
#else
		if (!this->neus.addData(newNeu, DEFAULT_TID)) {
			cerr << "Error: Could not add neutrophil in sproutAgent()" << endl;
			exit(-1);
		}
#endif
		break;
	}
	}
}

void WHWorld::sproutCellsInArea(int num, int patchType, int cellType,
		int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, bool bloodORtiss)
{
	cout << "sprout cells in area" << endl;
#pragma omp parallel
	{
		Patch* tempPatchPtr;
		int tid = 0;
		int nth = 1;
		int nsprouted = 0;
		int zrange = zmax - zmin;
		int yrange = ymax - ymin;
		int xrange = xmax - xmin;

#ifdef _OMP
		// Get thread id in order to access the seed that belongs to this thread
		tid = omp_get_thread_num();
		nth = omp_get_num_threads();
#endif
		int local_num = num/nth;
		if (tid < (num%nth)) ++local_num;

		int itercount = 0;
		int iterlimit = local_num * 3;
		int *reservoir = new int[local_num];

		while ((nsprouted <= local_num) && (itercount < iterlimit))
		{
			int iz = (rand_r(&(this->seeds[tid]))%zrange) + zmin;
			int iy = (rand_r(&(this->seeds[tid]))%yrange) + ymin;
			int ix = (rand_r(&(this->seeds[tid]))%xrange) + xmin;

			int in = ix + iy*nx + iz*nx*ny;
			if (in < 0 || in >= this->np)
			{
				cout << "ERROR: in sproutAgentnArea: [x, y, z]: " << ix << ", " << iy
						<< ", " << iz << endl;
				exit(-1);
			}
			// try to sprout on this patch for cells
			tempPatchPtr = &(this->worldPatch[in]);
			if (tempPatchPtr->type[read_t] != patchType)
			{
				continue;
			}
			if (!tempPatchPtr->setOccupied()) {
				tempPatchPtr->occupiedby[write_t] = cellType;
				reservoir[nsprouted] = in;
				createCell(cellType, bloodORtiss, tempPatchPtr);
				++nsprouted;
			}
			itercount++;
		}
		// DEBUG
		if (nsprouted < local_num)
		{
			cout << "[" << tid << "]" << "Warning: only " << nsprouted <<
					" out of " << local_num << "cells sprouted" << endl;
			exit(-1);
		}
	}

}

bool WHWorld::searchIntArray(int *arr, int size, int v)
{
	for (int i = 0; i < size; i++)
	{
		if (arr[i] == v)
		{
			return true;
		}
	}
	return false;
}

// return random float in [fractionProteinLB, fractionProteinUB]
float getProteinPatchFraction()
{
	int LB 	= static_cast<int>(fractionProteinLB*100.0f);
	int UB 	= static_cast<int>(fractionProteinUB*100.0f);
	int range = UB - LB;
	int r		= LB + (rand()%range);
	return static_cast<float>(r)/100.0f;
}

void WHWorld::sproutAgentInArea(int num, int patchType, agent_t agentType,
		int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, int bloodORtiss) {

	// DEBUG vis
//	return;

    bool isECM = (oc <= agentType) && (agentType <= fha);
//	int count = 0, temp = 0;
	vector <int> patchlist;
	int *reservoir = new int [num];
	int *patchesToSprout = NULL;
	for (int i = 0; i < num; i++) reservoir[i] = -1;
	Patch* tempPatchPtr;
	Agent* tempAgentPtr;
	int in, agentIndex, max;
#ifdef SPARSE_WORLD_OPT
	if (agentType == platelet || agentType == neutrophil
			|| agentType == macrophag || agentType == fibroblast) {
		sproutCellsInArea(num, patchType, agentType,
				xmin, xmax, ymin, ymax, zmin, zmax, bloodORtiss);
		return;
	} 
#endif	// SPARSE_WORLD_OPT
	for (int izz = zmin; izz < zmax + 1; izz++) {
		for (int iyy = ymin; iyy < ymax + 1; iyy++) {
			for (int ixx = xmin; ixx < xmax + 1; ixx++) {
				in = ixx + iyy*nx + izz*nx*ny;
				// Try another patch if this one is out of bounds or the wrong type or occupied
				if (ixx < 0 || ixx >= nx || iyy < 0 || iyy >= ny || izz < 0 || izz >= nz) continue;
				if (WHWorld::worldPatch[in].type[read_t] != patchType) continue;
				// DEBUG
//				temp++;
				if (this->worldPatch[in].isOccupied() == false) {
					patchlist.push_back(in);
				}
			}
		}
	}

//	printf("possible patches: %d\n", temp);

	int nAvailPatches = patchlist.size();

    if (isECM) {
        if (patchlist.size() < 0) {  // No available patches
            cerr << "Error: SproutAgent, no available patches within bounds! " << endl;
            delete[] reservoir;
            exit(-1);
        }

        for (int i = 0; i < num; i++) {
            int randnumber = rand() % nAvailPatches;
            reservoir[i] = patchlist[randnumber];  // Prepare 'num' random patches
        }

        patchesToSprout = reservoir;

    } else {
        if (patchlist.size() < num) {  // Not enough available patches
            cerr << "Error: SproutAgent, not enough available patches within bounds! " << endl;
            cerr << "   patchlist.size(): " << nAvailPatches << "  num: " << num << endl;
            cerr << "   bounds x: [" << xmin << ", " << xmax << "]";
            cerr <<        "   y: [" << ymin << ", " << ymax << "]";
            cerr <<        "   z: [" << zmin << ", " << zmax << "]" << endl;
            delete[] reservoir;
            exit(-1);
        }

        // Shuffle available patch list
        std::random_shuffle(patchlist.begin(), patchlist.end());

        patchesToSprout = &(patchlist[0]);
    }


	// Sprout agent on each patch in patchesToSprout
//#pragma omp parallel for
	for (int i = 0; i < num; i++) {
		int in = patchesToSprout[i];
		if (in < 0 || in > (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny) continue;
#ifdef NO_CELLS_IN_DAMZONE
		if (this->worldPatch[in].isInDamZone()) continue; 
#endif
//	// NEW DEBUG
//		if (!isECM && this->worldPatch[in].indice[0] > 100) continue;

		switch (agentType) {

		case platelet: {
			tempPatchPtr = &(this->worldPatch[in]);
			Platelet* newPlatelet = new Platelet(tempPatchPtr);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewPlats[tid]->push_back(newPlatelet);
#else
			if (!this->plats.addData(newPlatelet, DEFAULT_TID)) {
				cerr << "Error: Could not add platelet in sproutAgent()" << endl;
				exit(-1);
			}
#endif
			this->worldPatch[in].setOccupied();
			this->worldPatch[in].occupiedby[write_t] = platelet;
			break;
		}
		case fibroblast: {
			tempPatchPtr = &(this->worldPatch[in]);
			Fibroblast* newFib = new Fibroblast(tempPatchPtr);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewFibs[tid]->push_back(newFib);
#else
			if (!this->fibs.addData(newFib, DEFAULT_TID)) {
				cerr << "Error: Could not add fibroblast in sproutAgent()" << endl;
				exit(-1);
			}
#endif
			this->worldPatch[in].setOccupied();
			this->worldPatch[in].occupiedby[write_t] = fibroblast;
			break;
		}
		case afibroblast: {
			tempPatchPtr = &(this->worldPatch[in]);
			Fibroblast* newFib = new Fibroblast(tempPatchPtr);
			newFib->fibActivation();
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewFibs[tid]->push_back(newFib);
#else
			if (!this->fibs.addData(newFib, DEFAULT_TID)) {
				cerr << "Error: Could not add fibroblast in sproutAgent()" << endl;
				exit(-1);
			}
#endif
			//        newFib->fibActivation();
			this->worldPatch[in].setOccupied();
			this->worldPatch[in].occupiedby[write_t] = fibroblast;
			break;
		}
		case macrophag: {
			tempPatchPtr = &(this->worldPatch[in]);
			Macrophage* newMac = new Macrophage(tempPatchPtr, bloodORtiss);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewMacs[tid]->push_back(newMac);
#else
			if (!this->macs.addData(newMac, DEFAULT_TID)) {
				cerr << "Error: Could not add macrophage in sproutAgent()" << endl;
				exit(-1);
			}
#endif
			this->worldPatch[in].setOccupied();
			this->worldPatch[in].occupiedby[write_t] = macrophag;
			break;
		}
		case neutrophil: {
			tempPatchPtr = &(this->worldPatch[in]);
			Neutrophil* newNeu = new Neutrophil(tempPatchPtr);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewNeus[tid]->push_back(newNeu);
#else
			if (!this->neus.addData(newNeu, DEFAULT_TID)) {
				cerr << "Error: Could not add neutrophil in sproutAgent()" << endl;
				exit(-1);
			}
#endif
			this->worldPatch[in].setOccupied();
			this->worldPatch[in].occupiedby[write_t] = neutrophil;
			break;
		}
		case oc: {
			this->worldECM[in].ocollagen[write_t] = this->worldECM[in].ocollagen[read_t] + 1;
			this->worldECM[in].isEmpty();
			break;
		}
		case oe: {
			this->worldECM[in].oelastin[write_t] = this->worldECM[in].oelastin[read_t] + 1;
			this->worldECM[in].isEmpty();
			break;
		}
		case nc: {
			float pf      = getProteinPatchFraction();
			float numcollf = pf*this->maxFCol; //0.9*this->maxFCol;
			int   numcoll = static_cast<int>(numcollf);
			this->worldPatch[in].initcollagen = true;
//			this->worldECM[in].ncollagen[write_t] = this->worldECM[in].ncollagen[read_t] + 1;
			this->worldECM[in].ncollagen[write_t] = this->worldECM[in].ncollagen[read_t] + numcoll;
			this->worldECM[in].isEmpty();

			// Update ECM polymer map
			// DEBUG Vis
			int ix = worldECM[in].indice[0];
			int iy = worldECM[in].indice[1];
			int iz = worldECM[in].indice[2];
//			if(iz == 72)
				setECM(in, m_col, numcollf);

			break;
		}
		case ne: {
			this->worldPatch[in].initelastin = true;
			this->worldECM[in].nelastin[write_t] = this->worldECM[in].nelastin[read_t] + 1;
			this->worldECM[in].isEmpty();
#ifdef VISUALIZATION
			// Update ECM polymer map
			incECM(in, m_ela, 2);
#endif	// VISUALIZATION
			break;
		}
		case oha: {
			this->initHAcenters.push_back(in);
			this->worldPatch[in].initHA = true;
			this->worldECM[in].addHAs(6);
			this->worldECM[in].empty[write_t] = false;
#ifdef VISUALIZATION
			// Update ECM polymer map
			incECM(in, m_hya, 6);
#endif	// VISUALIZATION
			break;
		}

		}
	}
	delete[] reservoir;
	return;
}

#ifdef OPT_CELL_SEEDING
/*
 * Optimized by:
 *  - If sprout in tissue:
 *     (*) Randomly choosing a target patch:
 *          - If is tissue and unoccupied, sprout
 * 			    - Else repeat (*)
 *  - Else (sprout in blood):
 *     (**)	Look at the list of capillary patches initialized in the setup stage
 *          - Create a list of unoccupied capillary patches
 *          - Pick randomly and sprout
 *          - Repeat until 'num' cells are sprouted
 */
void WHWorld::sproutAgentInWorld(int num, int patchType, int agentType, bool bloodORtiss) {
	int count = 0;
	vector <int> patchlist;
	int* reservoir = new int [num];
	for (int i = 0; i < num; i++) reservoir[i] = -1;
	Patch* tempPatchPtr;
	Agent* tempAgentPtr;
	int in, agentIndex, max;
	int totalNumPatches = this->(n-1)x*this->(ny-1)*this->(nz-1);	 // Is this accurate?
	int numfound;
	int counter;
	int threshold;
	vector<int> unoccupiedCaps;
	switch (patchType) {
	case tissue:
		numfound = 0;
		counter = 0;  // counter is guard for infinite loop
		threshold = totalNumPatches/100;
		while (num < numfound && counter < threshold) {
			int i = rand() % totalNumPatches;
			if (this->worldPatch[i].type == tissue && !this->worldPatch[i].isOccupied()) {
				// Found unoccupied tissue patches
				reservoir[numfound++] = i;
				// Ensure we do not get index duplicates from later iterations
				this->worldPatch[i].setOccupiedLight();
			}
			counter++;
		}
		if (numfound < num) {
			cout << " sprout agent error, only found " << numfound << " out of " << num << endl;
			exit(-1);
		}
		break;
	case capillary:
		/* Since we do not have that many capillary patches, and the list of
		 * capillary patches do not change, we keep the list so we can pick one off
		 * of the list */
		for (vector<int>::iterator cpi = this->CapPatchIndices.begin();
				cpi != this->CapPatchIndices.end(); cpi++) {
			int index = *cpi;
			if (this->worldPatch[index].isOccupied()) {
				// Found an unoccupied capillary
				unoccupiedCaps.push_back(index);
			}
		}
		int numCandidates = unoccupiedCaps.size();
		if (numCandidates < num) {
			cout << " sprout agent error, no available cap patch within bounds! (" << numCandidates
					<< "/" << num << ")" << endl;
			exit(-1);
		}
		numfound = 0;
		counter = 0;
		threshold = numCandidates*4;
		while (num < numfound && counter < threshold) {
			int i = rand() % numCandidates;
			if (!this->worldPatch[i].isOccupied()) {		// TODO: Fix, since this check is almost always redundant
				reservoir[numfound++] = i;
				// Ensure we do not get index duplicates from later iterations
				this->worldPatch[i].setOccupiedLight();
			}
			counter++;
		}
		break;
	case epithelium:
		cerr << "Error: spoutAgentInWorld() EPITHELIUM not yet supported" << endl;
		exit(-1);
		break;
	case damage:
		cerr << "Error: spoutAgentInWorld() DAMAGE not yet supported" << endl;
		exit(-1);
		break;
	default:
		cerr << "Error: spoutAgentInWorld() Unrecognized patch type" << endl;
		exit(-1);
		break;
	}

	// Sprout agent on each patch in reservoir



	for (int i = 0; i < num; i++) {

		if (reservoir[i] < 0 || reservoir[i] > (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny) continue;
		switch (agentType) {
		case platelet: {
			tempPatchPtr = &(this->worldPatch[reservoir[i]]);
			Platelet* newPlatelet = new Platelet(tempPatchPtr);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewPlats[tid]->push_back(newPlatelet);
#else
			if (!this->plats.addData(newPlatelet, DEFAULT_TID)) {
				cerr << "Error: Could not add platelet in sproutAgentInWorld()" << endl;
				exit(-1);
			}
#endif
			//			this->worldPatch[reservoir[i]].occupied = true;
			this->worldPatch[reservoir[i]].setOccupied();
			this->worldPatch[reservoir[i]].occupiedby = platelet;
			break;
		}
		case fibroblast: {
			tempPatchPtr = &(this->worldPatch[reservoir[i]]);
			Fibroblast* newFib = new Fibroblast(tempPatchPtr);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewFibs[tid]->push_back(newFib);
#else
			if (!this->fibs.addData(newFib, DEFAULT_TID)) {
				cerr << "Error: Could not add fibroblast in sproutAgentInWorld()" << endl;
				exit(-1);
			}
#endif
			//			this->worldPatch[reservoir[i]].occupied = true;
			this->worldPatch[reservoir[i]].setOccupied();
			this->worldPatch[reservoir[i]].occupiedby = fibroblast;
			break;
		}
		case macrophag: {
			tempPatchPtr = &(this->worldPatch[reservoir[i]]);
			Macrophage* newMac = new Macrophage(tempPatchPtr, bloodORtiss);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewMacs[tid]->push_back(newMac);
#else
			if (!this->macs.addData(newMac, DEFAULT_TID)) {
				cerr << "Error: Could not add macrophage in sproutAgentInWorld()" << endl;
				exit(-1);
			}
#endif
			//			this->worldPatch[reservoir[i]].occupied = true;
			this->worldPatch[reservoir[i]].setOccupied();
			this->worldPatch[reservoir[i]].occupiedby = macrophag;
			break;
		}
		case neutrophil: {
			tempPatchPtr = &(this->worldPatch[reservoir[i]]);
			Neutrophil* newNeu = new Neutrophil(tempPatchPtr);
#ifdef _OMP
			int tid = omp_get_thread_num();
			this->localNewNeus[tid]->push_back(newNeu)
#else
						if (!this->neus.addData(newNeu, DEFAULT_TID)) {
							cerr << "Error: Could not add neutrophil in sproutAgentInWorld()" << endl;
							exit(-1);
						}
#endif
			//			this->worldPatch[reservoir[i]].occupied = true;
			this->worldPatch[reservoir[i]].setOccupied();
			this->worldPatch[reservoir[i]].occupiedby = neutrophil;
			break;
		}
		}
	}
	delete[] reservoir;
	return;
}
#endif  //OPT_CELL_SEEDING

int WHWorld::countPatchType(int whichType) {
	if (whichType == tissue || whichType == epithelium || whichType == capillary) {
		Patch::numOfEachTypes[whichType] = 0;
		for (int iz = 0; iz < this->nz; iz++) {
			int currCount = 0;
#pragma omp parallel for reduction(+:currCount)
			for (int iy = 0; iy < this->ny; iy++) {
				for (int ix = 0; ix < this->nx; ix++) {
					int in = ix + iy*nx + iz*nx*ny;
					if (this->worldPatch[in].type[read_t] == whichType) {
						currCount++;
						//						Patch::numOfEachTypes [whichType]++;
					}
				}
			}
			Patch::numOfEachTypes [whichType] += currCount;
		}
	} else if (whichType == damage) {
		Patch::numOfEachTypes[whichType] = 0;
		for (int iz = 0; iz < this->nz; iz++) {
			int currCount = 0;
#pragma omp parallel for reduction(+:currCount)
			for (int iy = 0; iy < this->ny; iy++) {
				for (int ix = 0; ix < this->nx; ix++) {
					int in = ix + iy*nx + iz*nx*ny;
					currCount += this->worldPatch[in].damage[read_t];
				}
			}
			Patch::numOfEachTypes [whichType] += currCount;
		}
	} else
		cout << "type must be 1, 2 , 3 or 4!" << endl;
	return Patch::numOfEachTypes[whichType];
}

int WHWorld::mmToPatch(double mm) {
	return mm*(this->patchpermm);
}

int WHWorld::reportTick(int hour, int day) {
	return (hour*2 + day*48);
}


void WHWorld::setSeed(unsigned s) {
	WHWorld::seed = s;
}


double WHWorld::reportHour() {
	return (WHWorld::clock)/2;
}

double WHWorld::reportDay() {
	return (WHWorld::clock)/48;
}

int WHWorld::countNeighborPatchType(int ix, int iy, int iz,  int patchType) {
	int neighborcount = 0;
	for (int dx = -1; dx < 2; dx++) {
		for (int dy = -1; dy < 2; dy++) {
			for (int dz = -1; dz < 2; dz++) {
				int neighborindex = (ix + dx) + (iy + dy)*nx + (iz + dz)*ny*nx;
				if (ix + dx < 0 || ix + dx >= nx || iy + dy < 0 || iy + dy >= ny || iz + dz < 0 || iz + dz >= nz) continue;
				if (dx == 0 && dy ==0 && dz ==0 ) continue;
				if (Agent::agentPatchPtr[neighborindex].type[read_t] == patchType) neighborcount++;
			}
		}
	}
	return neighborcount;
}

/*
 * Steps:
 *   1. If OMP, add cells from thread-local lists to corresponding global lists
 *   2. Perform updates
 */
void WHWorld::updateCellsInitial() {
	// Cell lists should be empty
	// Add new cells
	// NS: Only fibroblast part is necessary?
#ifdef _OMP
	cerr << "	updateCell() _OMP" << endl;
	// TODO: parallelize
	int numThreads = omp_get_num_threads();
	for (int tid = 0; tid < numThreads; tid++) {
		/********************************************
		 * PLATELETS                                *
		 ********************************************/
		vector<Platelet*>* pvec_ptr = localNewPlats[tid];
		for (vector<Platelet*>::iterator plat_it = pvec_ptr->begin(); plat_it != pvec_ptr->end(); plat_it++) {
			Platelet* newPlat = *plat_it;
			if (!plats.addData(newPlat, tid)) {
				cerr << "Error: Could not add platelet" << endl;
				exit(-1);
			}
		}
		pvec_ptr->clear();
		/********************************************
		 * NEUTROPHILS                              *
		 ********************************************/
		vector<Neutrophil*>* nvec_ptr = localNewNeus[tid];
		for (vector<Neutrophil*>::iterator neu_it = nvec_ptr->begin(); neu_it != nvec_ptr->end(); neu_it++) {
			Neutrophil* newNeu = *neu_it;
			if (!neus.addData(newNeu, tid)) {
				cerr << "Error: Could not add neutrophil" << endl;
				exit(-1);
			}
		}
		nvec_ptr->clear();
		/********************************************
		 * MACROPHAGES                              *
		 ********************************************/
		vector<Macrophage*>* mvec_ptr = localNewMacs[tid];
		for (vector<Macrophage*>::iterator mac_it = mvec_ptr->begin(); mac_it != mvec_ptr->end(); mac_it++) {
			Macrophage* newMac = *mac_it;
			if (!macs.addData(newMac, tid)) {
				cerr << "Error: Could not add macrophage" << endl;
				exit(-1);
			}
		}
		mvec_ptr->clear();
		/********************************************
		 * FIBROBLASTS                              *
		 ********************************************/
		vector<Fibroblast*>* fvec_ptr = localNewFibs[tid];
		for (vector<Fibroblast*>::iterator fib_it = fvec_ptr->begin(); fib_it != fvec_ptr->end(); fib_it++) {
			Fibroblast* newFib = *fib_it;
			if(!fibs.addData(newFib, tid)) {
				cerr << "Error: Could not add fibroblast" << endl;
				exit(-1);
			}
		}
		fvec_ptr->clear();
	}
#endif
	/********************************************
	 * FIBROBLASTS                              *
	 ********************************************/
	// No need for deletion since these are new cells


	int fibsSize = fibs.size();

#pragma omp parallel for
	for (int i = 0 ; i < fibsSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		Fibroblast* fib = fibs.getDataAt(i);
		if (!fib) continue;

		fib->updateAgent();

	}
	Fibroblast::numOfFibroblasts = fibs.actualSize();


	/********************************************
	 * MACROPHAGES                              *
	 ********************************************/
	int macsSize = macs.size();
#pragma omp parallel for
	for (int i = 0; i < macsSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		Macrophage* mac = macs.getDataAt(i);
		if (!mac) continue;
		mac->updateAgent();
	}
	Macrophage::numOfMacrophage = macs.actualSize();
	/********************************************
	 * NEUTROPHILS                              *
	 ********************************************/
	int neusSize = neus.size();
#pragma omp parallel for
	for (int i = 0; i < neusSize; i++) {
#ifdef _OMP
		int tid = omp_get_thread_num();
#else
		int tid = DEFAULT_TID;
#endif
		Neutrophil* neu = neus.getDataAt(i);
		if (!neu) continue;
		neu->updateAgent();
	}
	Neutrophil::numOfNeutrophil = neus.actualSize();
	/********************************************
	 * PLATELETS                                *
	 ********************************************/
	int platsSize = plats.size();
	//#pragma omp parallel for
	for (int i = 0; i < platsSize; i++) {
		Platelet* plat = plats.getDataAt(i);
		plat->updateAgent();
	}
	Platelet::numOfPlatelets = plats.actualSize();
}


// TODO: Relocate?

int WHWorld::userInput() {

	// Read input parameters from user-specified file
	ifstream infile(util::getInputFileName());

	int numChem = -1;
	// TODO: Make this check for specific tag (field name)

	if (infile.is_open()) {
		char garbage[100];
		/********************************************
		 * CHEMICALS                                *
		 ********************************************/
		cout << "Reading the number of baseline chemicals" << endl;
		float temp;
		infile >> garbage;
		infile >> temp;
		//cout << garbage;
		//cin >> garbage;

		numChem = temp;
		this->baselineChem.resize(temp);
		cout << "The number of baseline chemicals are " << this->baselineChem.size() << endl;
		for (int ichem = 0; ichem < baselineChem.size(); ichem++) {
			cout << "Reading the baselineChemical " << ichem << endl;
			infile >> garbage;
			infile >> this->baselineChem[ichem];
			//cout << garbage;
			//cin >> garbage;
			cout << "baselineChem " << ichem << " is " << this->baselineChem[ichem] << endl;
		}

		// convert average chemical concentration (pg/mL) to total chem (pg)
		float pw = this->patchlength;
		int   np = nx*ny*nz;
		float worldVolume = pw*pw*pw * (float) np;		// mm^3
		this->worldVmL  = worldVolume * 0.001;		// mL
		printf("pw: %f, np: %d\n", pw, np);
		printf("World volume: %f\tWorld Volume (mL): %f\n", worldVolume, worldVmL);
		for (int ichem = 0; ichem < baselineChem.size(); ichem++) {
			this->baselineChem[ichem] *= worldVmL;		// pg
		}


		/********************************************
		 * WOUND                                    *
		 ********************************************/
		cout << "Reading the wound depth (x-radius)" << endl;
		infile >> garbage;
		infile >> this->wound[0];
		//cout << garbage;
		//cin >> garbage;
		cout << "wound depth (x) is " << this->wound[0] << endl;

		cout << "Reading the wound y-radius" << endl;
		infile >> garbage;
		infile >> this->wound[1];
		//cout << garbage;
		//cin >> garbage;
		cout << "wound y-radius " << this->wound[1] << endl;

		cout << "Reading the wound z-radius" << endl;
		infile >> garbage;
		infile >> this->wound[2];
		//cout << garbage;
		//cin >> garbage;
		cout << "wound z-radius is " << this->wound[2] << endl;

		cout << "Reading the wound severity" << endl;
		infile >> garbage;
		infile >> this->wound[3];
		//cout << garbage;
		//cin >> garbage;
		cout << "wound severity is " << this->wound[3] << endl;

		/********************************************
		 * CELLS                                    *
		 ********************************************/
		cout << "Reading the initial number of types of cells" << endl;
		int tempCells;
		infile >> garbage;
		infile >> tempCells;
		//cout << garbage;
		//cin >> garbage;
		this->initialCells.resize(tempCells);
		cout << "The number of types of cells are " << this->initialCells.size() << endl;
		for (int icell = 0; icell< initialCells.size(); icell++) {
			cout << "Reading the initial # of cell type " << icell + 1 << endl;
			infile >> garbage;
			infile >> this->initialCells[icell];
			//cout << garbage;
			//cin >> garbage;
			cout<<"# of cell " << icell << " is " << this->initialCells[icell] << endl;
		}

		/********************************************
		 * TREATMENT                                *
		 ********************************************/
		cout << "Reading the treatment option" << endl;
		int treatment;
		infile >> garbage;
		infile >> treatment;
		//cout << garbage;
		//cin >> garbage;
		this->treatmentOption = treatment;
		if (treatment == 0) {
			cout << "The treatment option is voice rest" << endl;
		}
		else if (treatment == 1) {
			cout << "The treatment option is resonant voice" << endl;
		}
		else if (treatment == 2) {
			cout << "The treatment option is spontaneous voice" << endl;
		}
		else {
			cout << "Wrong option! Enter either 0(voice rest), 1(resonant voice), or 2(spontaneous) please!" << endl;
			cin >> treatment;
			this->treatmentOption = treatment;
		}



		/********************************************
		 * CYTOKINE PROPERTIES                      *
		 ********************************************/
		cout << "Reading Cytokine Properties" << endl;
		float D;
		int HL_s;
		char tag[100];
		do {
			if (numChem < 1) {
				cerr << "Warning: No chem allocated!!!" << endl;
				break;
			}
			// Allocate memory for diffusion coefficients and half-lifes
			this->D = (float*) malloc(sizeof(float) * numChem);
			this->HalfLifes = (int *) malloc(sizeof(int) * numChem);
			while (infile >> tag) {
				//infile >> tag;
				if (!strcmp(tag, "D:")) {
					infile >> D;
					for (int ic = 0; ic < numChem; ic++) {
						this->D[ic] = D;
					}
					cout << "	setting all D to: " << D << endl;
				} else if (!strcmp(tag, "HL:")) {
					infile >> HL_s;
					for (int ic = 0; ic < numChem; ic++) {
						this->HalfLifes[ic] = HL_s;
					}
					cout << "	seeting all HLs to: " << HL_s << endl;
				} else if (!strcmp(tag, "D_TNF:")) {
					infile >> D;
					this->D[TNF] = D;
					cout << "	D_TNF: " << D << endl;
				} else if (!strcmp(tag, "HL_TNF:")) {
					infile >> HL_s;
					this->HalfLifes[TNF] = HL_s;
					cout << "	HL_TNF: " << HL_s << endl;
				} else if (!strcmp(tag, "D_TGF:")) {
					infile >> D;
					this->D[TGF] = D;
					cout << "	D_TGF: " << D << endl;
				} else if (!strcmp(tag, "HL_TGF:")) {
					infile >> HL_s;
					this->HalfLifes[TGF] = HL_s;
					cout << "	HL_TGF: " << HL_s << endl;
				} else if (!strcmp(tag, "D_FGF:")) {
					infile >> D;
					this->D[FGF] = D;
					cout << "	D_FGF: " << D << endl;
				} else if (!strcmp(tag, "HL_FGF:")) {
					infile >> HL_s;
					this->HalfLifes[FGF] = HL_s;
					cout << "	HL_FGF: " << HL_s << endl;
				} else if (!strcmp(tag, "D_MMP8:")) {
					infile >> D;
					this->D[MMP8] = D;
					cout << "	D_MMP8: " << D << endl;
				} else if (!strcmp(tag, "HL_MMP8:")) {
					infile >> HL_s;
					this->HalfLifes[MMP8] = HL_s;
					cout << "	HL_MMP8: " << HL_s << endl;
				} else if (!strcmp(tag, "D_IL1beta:")) {
					infile >> D;
					this->D[IL1beta] = D;
					cout << "	D_IL1beta: " << D << endl;
				} else if (!strcmp(tag, "HL_IL1beta:")) {
					infile >> HL_s;
					this->HalfLifes[IL1beta] = HL_s;
					cout << "	HL_IL1beta: " << HL_s << endl;
				} else if (!strcmp(tag, "D_IL6:")) {
					infile >> D;
					this->D[IL6] = D;
					cout << "	D_IL6: " << D << endl;
				} else if (!strcmp(tag, "HL_IL6:")) {
					infile >> HL_s;
					this->HalfLifes[IL6] = HL_s;
					cout << "	HL_IL6: " << HL_s << endl;
				} else if (!strcmp(tag, "D_IL8:")) {
					infile >> D;
					this->D[IL8] = D;
					cout << "	D_IL8: " << D << endl;
				} else if (!strcmp(tag, "HL_IL8:")) {
					infile >> HL_s;
					this->HalfLifes[IL8] = HL_s;
					cout << "	HL_IL8: " << HL_s << endl;
				} else if (!strcmp(tag, "D_IL10:")) {
					infile >> D;
					this->D[IL10] = D;
					cout << "	D_IL10: " << D << endl;
				} else if (!strcmp(tag, "HL_IL10:")) {
					infile >> HL_s;
					this->HalfLifes[IL10] = HL_s;
					cout << "	HL_IL10: " << HL_s << endl;
				} else {
					cout << "	invalid tag: " << tag << endl;
				}
			}
		} while (0);

		infile.close();
	}  // end of if file opens properly
	else {
		cout << "Cannot open file!" << endl;
	}

	return 0;
}

