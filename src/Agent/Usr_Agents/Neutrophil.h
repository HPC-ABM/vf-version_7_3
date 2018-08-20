/* 
 * File: Neutrophil.h
 *
 * File Contents: Contains declarations for the Neutrophil class.
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

#ifndef NEUTROPHIL_H
#define	NEUTROPHIL_H

#include "../Agent.h"
#include "../../World/Usr_World/woundHealingWorld.h"
#include <stdlib.h>
#include <vector>

using namespace std;

/*
 * NEUTROPHIL CLASS DESCRIPTION:
 * Neutrophil is a derived class of the parent class Agent. 
 * It manages all neutrophil agents. 
 * It is used to initialize a neutrophil, carry out its biological function,
 * activate it, make it move, and make it die.
 */
class Neutrophil: public Agent {
 public:
    /*
     * Description:	Default neutrophil constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    Neutrophil();

    /*
     * Description:	Neutrophil constructor. Initializes neutrophil attributes.
     *
     * Return: void
     *
     * Parameters: patchPtr  -- Pointer to patch on which neutrophil will
     *                          reside. NOTE: The pointer cannot be NULL.
     */ 
    Neutrophil(Patch* patchPtr);

    /*
     * Description:	Neutrophil destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    ~Neutrophil(); 

    /*
     * Description:	Performs biological function of a neutrophil.
     *
     * Return: void
     *
     * Parameters: void
     */
    void cellFunction();

    /*
     * Description:	Performs biological function of an unactivated neutrophil.
     *
     * Return: void
     *
     * Parameters: void
     */
    void neu_cellFunction(); 

    /*
     * Description:	Performs biological function of an activated neutrophil.
     *
     * Return: void
     *
     * Parameters: void
     */
    void aneu_cellFunction();

    /*
     * Description:     Call cytokine production functions for activated neutrophils
     *
     *                  TODO(NS): INSERT REF
     *
     *
     * Return: void
     *
     * Parameters: void
     */
    void aneu_produce_cytokines();

    /*
     * Description:	Activates an unactivated neutrophil. 
     *              Updates the neutrophil class members.
     *
     * Return: void
     *
     * Parameters: void
     */
    void neuActivation();

    /*
     * Description:	Moves a neutrophil along its preferred chemical gradient.
     *
     * Return: void
     *
     * Parameters: void
     */
    void neuSniff();

    /*
     * Description:	Performs neutrophil death. Updates the 
     *              neutrophil class members. Does not update numOfNeutrophil;
     *              this must be done elsewhere.
     *
     * Return: void
     *
     * Parameters: void
     */
    void die();
  
    // Keeps track of the quantitiy of living neutrophils.
    static int numOfNeutrophil;


    // Parameters involved in synthesis of TNF, MMP8 by activated neutrophils:
    static float cytokineSynthesis[21];
    // Parameters involved in neutrophil activation
    static float activation[4];
    // Parameters involved in neutrophil death
    static float death[2];


};

#endif	/* Neutrophil_H */
