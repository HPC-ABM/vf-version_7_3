/* 
 * File: Macrophage.h
 *
 * File Contents: Contains declarations for the Macrophage class.
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

#ifndef MACROPHAGE_H
#define	MACROPHAGE_H

#include "../Agent.h"
#include "../../World/Usr_World/woundHealingWorld.h"
#include <stdlib.h>
#include <vector>

using namespace std;

/*
 * MACROPHAGE CLASS DESCRIPTION:
 * Macrophage is a derived class of the parent class Agent. 
 * It manages all macrophage agents. 
 * It is used to initialize a macrophage, carry out its biological function,
 * activate it, deactivate it, make it move, and make it die.
 */
class Macrophage: public Agent {
 public:

    /*
     * Description:	Default macrophage constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    Macrophage();

    /*
     * Description:	Macrophage constructor. 
     *              Initializes macrophage class members.
     *
     * Return: void
     * 
     * Parameters: x            -- Position of macrophage in x dimension
     *             y            -- Position of macrophage in y dimension
     *             z            -- Position of macrophage in z dimension
     *             bloodOrTiss  -- Enumic type for whether macrophage is 
     *                             located in blood or in tissue.
     */
    Macrophage(int x, int y, int z, int bloodOrTiss);

    /*
     * Description:	Macrophage constructor. Initializes macrophage attributes.
     *
     * Return: void
     *
     * Parameters: patchPtr     -- Pointer to patch on which macrophage will
     *                             reside. NOTE: The pointer cannot be NULL.
     *             bloodOrTiss  -- Enumic type for whether macrophage is 
     *                             located in blood or in tissue.
     */ 
    Macrophage(Patch* patchPtr, int bloodOrTiss);

    /*
     * Description:	Macrophage destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    ~Macrophage();

    /*
     * Description:	Performs biological function of a macrophage.
     *
     * Return: void
     *
     * Parameters: void
     */
    void cellFunction();

    /*
     * Description:	Performs biological function of an unactivated macrophage.
     *
     * Return: void
     *
     * Parameters: void
     */
    void mac_cellFunction();

    /*
     * Description:	Performs biological function of an activated macrophage.
     *
     * Return: void
     *
     * Parameters: void
     */
    void activatedmac_cellFunction();

    /*
     * Description:     Call cytokine production functions for myofibroblasts
     *
     *                  TODO(NS): INSERT REF
     *
     *
     * Return: void
     *
     * Parameters: void
     */
    void amac_produce_cytokines();


    /*
     * Description:	Moves a macrophage along its preferred chemical gradient.
     *
     * Return: void
     *
     * Parameters: void
     */
    void macSniff();

    /*
     * Description:	Activates an unactivated macrophage. 
     *              Updates the macrophage class members.
     *
     * Return: void
     *
     * Parameters: void
     */
    void macActivation();

    /*
     * Description:	Deactivates an activated macrophage. 
     *              Updates the macrophage class members.
     *
     * Return: void
     *
     * Parameters: void
     */
    void macDeactivation(); 

    /*
     * Description:	Performs macrophage death. Updates the 
     *              macrophage class members. Does not update numOfMacrophage;
     *              this must be done elsewhere.
     *
     * Return: void
     *
     * Parameters: void
     */
    void die();

    // Keeps track of the quantitiy of living macrophages.
    static int numOfMacrophage;
    /* Parameters involved in synthesis of TNF, TGF, FGF, IL1, IL6, IL8, IL10
     * by activated macrophages: */
    static float cytokineSynthesis[82];
    // Parameters involved in macrophage activation and deactivation
    static float activation[5];

    // Keeps track of whether the macrophage is in blood or in tissue.
    int bloodORtissue;

};

#endif	/* MACROPHAGE_H */
