/* 
 * File: enums.h
 * 
 * File Contents: Contains enum declarations
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *
 * Copyright (c) 2013 by A. Najafi-Yazdi
 *
 * Created on May 19, 2013, 7:55 AM
 */
/*****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues
#ifndef ENUMS_H
#define	ENUMS_H

#define REAL double

// Chemical update types
enum chemupdate_t {
	cpuUpdate,
	gpuUpdate
};

// Diffusion padding schemes
enum pad_t {
	pClampToBorder,
	pRightWall,
	pMirror,
	pConstantVF
};

// Types of agents & ECM managers
enum agent_t { 
  unoccupied = -1,
  platelet = 0,
  fibroblast = 1,
  neutrophil = 3,
  macrophag = 5,
  afibroblast = 6,
  aneutrophil = 7,
  amacrophag = 8,
  oc = 9,  // Original collagen
  nc = 10,  // New collagen
  fc = 11,  // Fragmented collagen
  oe = 12,  // Original elastin
  ne = 13,  // New elastin
  fe = 14,  // Fragmented elastin
  oha = 15,  // Original hyaluronan
  nha = 16,  // New hyaluronan
  fha = 17  // Fragmented hyaluronan
};

// Cell plots index
enum cell_i{
	p_neu,
	p_anu,
	p_mac,
	p_amc,
	p_fib,
	p_afb,
	p_plt,
	p_dam,
	p_celltotal
};

// ECM map index
enum ecm_i{
	m_col,
	m_ela,
	m_hya,
	m_ecmtotal
};

// LP section index
enum lp_i{
    lp_di,
    lp_ii,
    lp_si,
    lp_total
};

/* Colors of patches and agents for visualization with Paraview 3.0
 * (Kitware(Clifton Park, New York), Sandia National Labs(Livermore, CA),
 * CSimSoft(American Fork, Utah)). */
enum color_t {
  cplatelet = 25,
  cafibroblast = 105,
  cfibroblast = 100,
  caneutrophil = 85,
  cneutrophil = 80,
  camacrophage = 65,
  cmacrophage = 60,

  ccollagen = 139,
  celastin = 129,
  cHA = 119,
  cfcollagen = 0,
  cfelastin = 0,
  cfHA = 0, 

  cnothing = 0,
  ctissue = 9,
  cepithelium = 87,
  ccapillary = 16,
  cdamage = 0,
  cunidentifiable = 0,
  cSLP =10,
  cILP =11,
  cDLP =12,
  cmuscle = 13
};

// Types of chemicals
enum chemical_t {
  // p: patch, d: delta (change during the tick)
  TNF = 0,
  TGF = 1,
  FGF = 2,
  MMP8 = 3,
  IL1beta = 4,
  IL6 = 5,
  IL8 = 6,
  IL10 = 7,
  NEUgrad = 0,
  MACgrad = 1,
  FIBgrad = 2
};

// Types of patches
enum patches_t {
  nothing = 0,
  blood = 1,
  tissue = 2,
  epithelium = 3,
  capillary = 4,
  pdamage = 5,
  unidentifiable = 6,
//  SLP = 7,
//  ILP = 8,
//  DLP = 9,
//  muscle = 10
};


// Types of treatment
enum treatment_t {
  voicerest = 0,
  resonantvoice = 1,
  spontaneousspeech = 2
};

// Time points within a tick
enum readwrite_t {
  read_t = 0,  // Start of a tick
  write_t = 1  // End of a tick
}; 

#endif	/* ENUMS_H */
