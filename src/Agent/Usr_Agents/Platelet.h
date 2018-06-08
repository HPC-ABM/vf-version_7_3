/* 
 * File: Platelet.h
 *
 * File Contents: Contains declarations for the Platelet class.
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

#ifndef PLATELET_H
#define	PLATELET_H

#include "../Agent.h"
#include "../../World/Usr_World/woundHealingWorld.h"
#include <stdlib.h>
#include <vector>

using namespace std;

/*
 * PLATELET CLASS DESCRIPTION:
 * Platelet is a derived class of the parent class Agent. 
 * It manages all platelet agents. 
 * It is used to initialize a platelet, carry out its biological function,
 * and make it die.
 */
class Platelet: public Agent {
 public:
    /*
     * Description:	Default platelet constructor. 
     *              Initializes the color class member only.
     *
     * Return: void
     *
     * Parameters: void
     */
    Platelet();

    /*
     * Description:	Platelet constructor. Initializes platelet class members.
     *
     * Return: void
     * 
     * Parameters: x  -- Position of platelet in x dimension
     *             y  -- Position of platelet in y dimension
     *             z  -- Position of platelet in z dimension
     */
    Platelet(int x, int y, int z);

    /*
     * Description:	Platelet constructor. Initializes platelet attributes.
     *
     * Return: void
     *
     * Parameters: patchPtr  -- Pointer to patch on which platelet will reside
     *                          NOTE: The pointer cannot be NULL. 
     */
    Platelet(Patch* patchPtr);

    /*
     * Description:	Platelet destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    ~Platelet(); 

    /*
     * Description:	Performs biological function of a platelet.
     *
     * Return: void
     *
     * Parameters: void
     */
    void cellFunction();

    /*
     * Description:	Performs platelet death. Updates the platelet class members. 
     *              Does not update numOfPlatelets; this must be done elsewhere.
     *
     * Return: void
     *
     * Parameters: void
     */
    void die();	    

    // Keeps track of the quantitiy of living platelets.
    static int numOfPlatelets;


    // Parameters involved in synthesis of TGF, IL1, MMP8 by platelets:
    static float cytokineSynthesis[3];
    

};

#endif	/* PLATELET_H */
