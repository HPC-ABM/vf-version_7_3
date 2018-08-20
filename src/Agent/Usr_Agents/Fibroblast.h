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
#include "../../Utilities/math_utils.h"

class ECM; 

#include <stdlib.h>
#include <vector>
#include <cmath>
#include <assert.h>

#include <omp.h>

// TODO (NS): INSERT REF
#define FIB_ACT_TGF_UB  3.375e-5
#define FIB_ACT_TGF_LB  6.750e-6
#define FIB_PRO_TGF_MX  3.375e-6        // Max stim level
#define FIB_PRO_TGF_MN  3.375e-5        // Min inhibitorial level

#define FIB_PROLIF_I0    96
#define FIB_PROLIF_I1   192
#define FIB_PROLIF_I2   360


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
     * Description:     Performs fibroblast death. Updates the
     *              fibroblast class members. Does not update numOfFibroblasts;
     *              this must be done elsewhere.
     *
     * Return: void
     *
     * Parameters: void
     */
    void die();


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

    /*************************************************************************
     * INACTIVATED FIBROBLAST RELATED FUNCTIONS                              *
     *************************************************************************/
    /*
     * Description:	Performs biological function of an unactivated fibroblast.
     *
     * Return: void
     *
     * Parameters: void
     */                                                                                                
    void fib_cellFunction();

    /*
     * Description:     wiggle()
     *
     * Return: void
     *
     * Parameters: void
     */
    void fib_migrate();

    /*
     * Description:     If enough stimulation, activate with some chance
     *                  Patch TGF:
     *                  (-inf,     3.375e-5] pg -> low chance    (25%)
     *                  (3.375e-5, 6.750e-6] pg -> medium chance (50%)
     *                  (6.750e-6,     +inf) pg -> high chance  (100%)
     *
     * Return: void
     *
     * Parameters: void
     */
    void fib_activate();

    /*
     * Description:     Performs biological aging by 1 tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void fib_degrade();

    /*
     * Description:     Activates an unactivated fibroblast.
     *              Updates the fibroblast class members.
     *
     * Return: void
     *
     * Parameters: void
     */
    void fibActivation();


    /*************************************************************************
     * ACTIVATED FIBROBLAST RELATED FUNCTIONS                                *
     *************************************************************************/
    /*
     * Description:	Performs biological function of an activated fibroblast.
     *
     * Return: void
     *
     * Parameters: void
     */
    void afib_cellFunction();

    /*
     * Description:     Proliferate at variable rates.
     *                  After injury:
     *                  Day 4-7  --   3.34 days
     *                  Day 8-14 --  10.10 days
     *                  Day 15+  -- 103.00 days
     *
     *                  Day  4 = Hour  96 -> FIB_PROLIF_I0
     *                  Day  8 = Hour 192 -> FIB_PROLIF_I1
     *                  Day 15 = Hour 360 -> FIB_PROLIF_I2
     *
     *                  TODO(NS): INSERT REF
     *
     *
     * Return: void
     *
     * Parameters: void
     */
    void afib_proliferate();

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
    void afib_produce_cytokines(float neighborHAs);

    float produce_tgf(REGULATORS_T uregs, REGULATORS_T dregs, REGULATORS_T coefs, float offset);


    /*
     * Description:     Case:   - Outside Wound -- Follow gradient to a healthy patch
     *                          - Inside Wound  -- wiggle()
     *
     * Return:          True,  if inside wound
     *                  False, otherwise
     *
     * Parameters: void
     */

    bool afib_migrate();

    /*
     * Description:     Make a list of damaged neighbors
     *                  If the list is NOT empty,
     *                          call make_<ecm>_monomers() functions
     *
     *
     *                  TODO(NS): INSERT REF
     *
     * Return: void
     *
     * Parameters: void
     */
    void afib_produce_ecms();

    void make_monomers(int dI_targets[3][3]);

    float  make_coll_monomers(
        int dI_target[3],
        REGULATORS_T uregs,
        REGULATORS_T dregs,
        REGULATORS_T coefs,
        float offset);
    float  make_elas_monomers(
        int dI_target[3],
        REGULATORS_T uregs,
        REGULATORS_T dregs,
        REGULATORS_T coefs,
        float offset);
    float  make_hyaluronans(
        int dI_target[3],
        REGULATORS_T uregs,
        REGULATORS_T dregs,
        REGULATORS_T coefs,
        float offset);
    float calc_stim(
        REGULATORS_T uregs,
        REGULATORS_T dregs,
        REGULATORS_T coefs,
        float offset,
        float max_stim = 2.0f);

    // Get random neighbor from the list
    void get_rn_from_list(
        vector <int> damagedneighbors,
        int &dX,
        int &dY,
        int &dZ);

    /*
     * Description:     Call die() function if max #proliferation reached
     *
     * Return: void
     *
     * Parameters: void
     */
    void afib_degrade();

    /*
     * Description:     Call de-activation function if damage is cleared
     *
     * Return: void
     *
     * Parameters: void
     */
    void afib_deactivate();

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
     * Description:	Hatches a new fibroblast on 'number' unoccupied neighbors.
     *              Does not update numOfFibroblasts; this must be done
     *              elsewhere.
     *
     * Return: void
     *
     * Parameters: number  -- Number of new fibroblasts to hatch
     */
    void hatchnewfibroblast(int number);


    int get_current_index(int &x, int &y, int &z, int &read_index);

    /*************************************************************************
     * ACTIVATED and INACTIVATED FIBROBLAST RELATED FUNCTIONS                *
     *************************************************************************/

    // Number of maximum proliferation per interval
    int maxProlif[3];
    // For make ecm function to use to check for migration
    bool isMigrated;


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
