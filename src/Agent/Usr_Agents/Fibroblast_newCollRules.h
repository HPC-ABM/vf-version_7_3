/* 
 * File: Fibroblast.h
 *
 * File Contents: Contains declarations for the Fibroblast class.
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

#ifndef FIBROBLAST_H
#define	FIBROBLAST_H

#include "../Agent.h"
#include "../../Patch/Patch.h"
#include "../../World/Usr_World/woundHealingWorld.h"

class ECM; 

#include <stdlib.h>
#include <vector>
#include <cmath>

#include <omp.h>

using namespace std;

/*
 * FIBROBLAST CLASS DESCRIPTION:
 * Fibroblast is a derived class of the parent class Agent. 
 * It manages all fibroblast agents. 
 * It is used to initialize a fibroblast, carry out its biological function,
 * activate it, deactivate it, make it move, and make it die.
 */
class Fibroblast: public Agent {
 public:
    /*
     * Description:	Default fibroblast constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    Fibroblast();

    /*
     * Description:	Fibroblast constructor. Initializes fibroblast attributes.
     *
     * Return: void
     *
     * Parameters: patchPtr  -- Pointer to patch on which fibroblast will
     *                          reside. NOTE: The pointer cannot be NULL.
     */ 
    Fibroblast(Patch* patchPtr);

    /*
     * Description:	Fibroblast constructor. 
     *              Initializes fibroblast class members.
     *
     * Return: void
     * 
     * Parameters: x  -- Position of fibroblast in x dimension
     *             y  -- Position of fibroblast in y dimension
     *             z  -- Position of fibroblast in z dimension
     */
    Fibroblast(int x, int y, int z); 

    /*
     * Description:	Fibroblast destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    ~Fibroblast(); 

    /*
     * Description:	Performs biological function of a fibroblast.
     *
     * Return: void
     *
     * Parameters: void
     */
    void cellFunction();

    /*
     * Description:	Performs biological function of an unactivated fibroblast.
     *
     * Return: void
     *
     * Parameters: void
     */                                                                                                
    void fib_cellFunction();                                                                                                                        

    /*
     * Description:	Performs biological function of an activated fibroblast.
     *
     * Return: void
     *
     * Parameters: void
     */
    void afib_cellFunction();                         

    /*
     * Description:	Moves a fibroblast along its preferred chemical gradient.
     *
     * Return: void
     *
     * Parameters: void
     */
    void fibSniff();
    void fibSniff2();
    void fibSniffOutDam();

    /*
     * Description:	Performs fibroblast death. Updates the 
     *              fibroblast class members. Does not update numOfFibroblasts;
     *              this must be done elsewhere.
     *
     * Return: void
     *
     * Parameters: void
     */
    void die();						

    /*
     * Description:	Activates an unactivated fibroblast. 
     *              Updates the fibroblast class members.
     *
     * Return: void
     *
     * Parameters: void
     */
    void fibActivation();

    /*
     * Description:	Deactivates an activated fibroblast. 
     *              Updates the fibroblast class members.
     *
     * Return: void
     *
     * Parameters: void
     */
    void fibDeactivation();

    /* 
     * Description: Copies the location of 'original' agent and initializes a 
     *              new fibroblast at a distance away determined by dx, dy, dz.
     *              NOTE: Target patch at distance of dx,dy,dz must be 
     *              unoccupied for proper functionality.
     *
     * Return: void
     *
     * Parameters: original  -- Agent to be copied
     *             dx        -- Difference in x-coordinate of the fibroblast's 
     *                          location relative to original's location.
     *             dy        -- Difference in y-coordinate of the fibroblast's 
     *                          location relative to original's location.
     *             dz        -- Difference in z-coordinate of the fibroblast's
     *                          location relative to original's location.
     *                          NOTE: dz = 0 because it is only 2D for now. TODO(Nuttiiya): I'm guessing this needs to change when you implement 3D?
    */
    void copyAndInitialize(Agent* original, int dx, int dy, int dz = 0);  

    /*
     * Description:	Sprouts original collagen on one of the activated fibroblast's 
     *              damaged neighbor patches.
     *
     * Return: void
     *
     * Parameters: meanTGF   -- Average TGF concentration of the activated 
     *                          fibroblast's neighbors
     * 				     meanFGF   -- Average FGF concentration of the activated
     *                          fibroblast's neighbors
     * 				     meanIL1   -- Average IL1 concentration of the activated
     *                          fibroblast's neighbors
     * 				     meanIL6   -- Average IL6 concentration of the activated
     *                          fibroblast's neighbors
     * 				     meanIL8   -- Average IL8 concentration of the activated
     *                          fibroblast's neighbors
     * 				     countnHA	 -- Number of new hyaluronan on the activated
     *                          fibroblast's neighbors
     * 				     countfHA	 -- Number of fragment hyaluronan on the activated
     *                          fibroblast's neighbors
     */
    void makeOCollagen(float meanTGF, float meanFGF, float meanIL1, float meanIL6, float meanIL8, int countnHA, int countfHA);
    void makeOCollagen2();

    /*
     * Description:	Sprouts original elastin on one of the activated fibroblast's 
     *              damaged neighbor patches.
     *
     * Return: void
     *
     * Parameters: meanTNF  -- Average TNF concentration of the activated
     *                         fibroblast's neighbors
     * 				     meanTGF  -- Average TGF concentration of the activated
     *                         fibroblast's neighbors
     * 				     meanFGF  -- Average FGF concentration of the activated
     *                         fibroblast's neighbors
     * 				     meanIL1  -- Average IL1 concentration of the activated
     *                         fibroblast's neighbors
     */
    void makeOElastin(float meanTNF, float meanTGF, float meanFGF, float meanIL1);

    /*
     * Description:	Sprouts new hyaluronan on one of the activated fibroblast's
     *              damaged neighbor patches.
     *
     * Return: void
     *
     * Parameters: meanTNF  -- Average TNF concentration of the activated
     *                         fibroblast's neighbors
     * 				     meanTGF  -- Average TGF concentration of the activated
     *                         fibroblast's neighbors
     * 				     meanFGF  -- Average FGF concentration of the activated
     *                         fibroblast's neighbors
     * 				     meanIL1  -- Average IL1 concentration of the activated
     *                         fibroblast's neighbors
     */
    void makeHyaluronan(float meanTNF, float meanTGF, float meanFGF, float meanIL1);

    /*
     * Description:	Hatches a new fibroblast on 'number' unoccupied neighbors.
     *              Does not update numOfFibroblasts; this must be done
     *              elsewhere.
     *
     * Return: void
     *
     * Parameters: number  -- Number of new fibroblasts to hatch
     */
    void hatchnewfibroblast(int number);

    void hatchnewAfibroblast(int number);

    // Keeps track of the quantitiy of living neutrophils.
    static int numOfFibroblasts;
    /* Parameters involved in synthesis of TNF, TGF, FGF, IL6, IL8
     * by activated fibroblasts: */
    static float cytokineSynthesis[32];
    // Parameters involved in fibroblast activation and deactivation
    static float activation[5];
    // Parameters involved in ECM synthesis
    static float ECMsynthesis[19];
    // Parameters invloved in FIbroblast proliferation
    static float proliferation[6];

};

#endif	/* FIBROBLAST_H */
