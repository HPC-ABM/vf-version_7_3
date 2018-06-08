/*
 * File: ECM.h
 *
 * File Contents: Contains declarations for ECM class
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */
/*****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#ifndef ECM_H
#define	ECM_H

#include <stdlib.h>
#include <vector>
#include "../Agent/Agent.h"
#include "../FieldVariable/FieldVariable.h"
#include "../enums.h"

#include <omp.h>

class World;
class WHWorld;
class Patch;

/*
 * ECM CLASS DESCRIPTION:
 * The ECM class manages ECM proteins. 
 * It is used to fragment ECM proteins, 
 * and signal and repair damage. 
 */
class ECM {
 typedef std::vector<int>::iterator HL_ITERATOR;
 public:
    /*
     * Description:	Default ECM constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    ECM();

    /*
     * Description:	ECM constructor. 
     *              Initializes ECM class members.
     *
     * Return: void
     * 
     * Parameters: x      -- Position of ECM manager in x dimension
     *             y      -- Position of ECM manager in y dimension
     *             z      -- Position of ECM manager in z dimension
     *             index  -- Patch row major index of ECM manager in world
     */
    ECM(int x, int y, int z, int index); 

    /*
     * Description:	ECM destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    ~ECM();

    /*
     * Description:	Performs biological function of ECM proteins.
     *
     * Return: void
     *
     * Parameters: void
     */
    void ECMFunction();

    /*
     * Description:	Repairs the damage on the patch where the ECM manager is.
     *
     * Return: void
     *
     * Parameters: void
     */
    void repairDamage();

    /*
     * Description:	Creates damage on the patch where the ECM manager is. If
     *              there are activated macrophages or activated neutrophils 
     *              already on the patch, the damage is cleared.
     *
     * Return: void
     *
     * Parameters: void
     */
    void dangerSignal();

    /*
     * Description:	Clears all flags present in the current tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void resetrequests();

    /*
     * Description:	Breaks down each new collagen ECM protein on this ECM 
     *              manager's patch into two fragmented collagen ECM proteins
     *
     * Return: void
     *
     * Parameters: void
     */
    void fragmentNCollagen();

    /*
     * Description:	Breaks down each new elastin ECM protein on this ECM 
     *              manager's patch into two fragmented elastin ECM proteins
     *
     * Return: void
     *
     * Parameters: void
     */
    void fragmentNElastin();

    /*
     * Description:	Breaks down each hyaluronan ECM protein on this ECM 
     *              manager's patch into two fragmented hyaluronan ECM proteins
     *
     * Return: void
     *
     * Parameters: void
     */
    void fragmentHA();

    /*
     * Description:	Update given ECM manager with all requests for fragment ECM
     *              from neighboring patches
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateECM();

    /* Description: Assigns true to the ECM manager's empty[] array at the end 
     *              of the tick if there are no ECM proteins for this manager.
     * 
     * Return: void
     * 
     * Parameters: void
     */
    void isEmpty();

    /* Description: 	HA related fields are set private to put more control on read/write
     *              synchronization of HA and HAlife. HA[write_t] should always be <=
     *              than HAlife.size(), reflecting current number of HA, whereas
     *              HA[read_t] should reflect the number of HAs that were alive at the
     *              beginning of the current tick.
     * 
     * Return: void
     * 
     * Parameters: nHA	--	Number of new HAs to add
     */
    void addHAs(int nHA);

    /* Description: 	HA related fields are set private to put more control on read/write
     *              synchronization of HA and HAlife. HA[write_t] should always be <=
     *              than HAlife.size(), reflecting current number of HA, whereas
     *              HA[read_t] should reflect the number of HAs that were alive at the
     *              beginning of the current tick. This function removes 1 HA at index
     *              iHA from HAlife and decrement HA[write_t].
     * 
     * Return: true if removal is successful, false otherwise
     * 
     * Parameters: iHA	--	Index into HAlife of HA to be removed
     */
    bool rmvHA(int iHA);

    bool decHAlife(int iHA);

    int getnHA();

    int getnHAlife();

    int getnfHA();

    void addColls(int nColl);
   
    bool isModifiedHA(); 
    int getnColl();

    void cleanFragmentedECM();


    // DEBUG
    void checknHA();
    bool checknHAbool();

    void checkHAdiff();

#ifdef OPT_ECM
    /*
     * Description:	Set dirty flag to indicate change in this ECM's attributes.
     *              Should be called whenever this ECM manager modifies its 
     *              own attributes. Dirty flag is checked when ECMupdate() is 
     *              called to see if update is necessary.
     *
     * Return: void
     *
     * Parameters: void
     */
    void set_dirty();

    /*
     * Description:	Set request-dirty flag to indicate that this ECM manager 
     *              has requested ECM fragments on a neighbor.
     *              Should be called whenever this ECM manager submits a
     *              request for fragments on a neighboring patch.
     * 				      Request-dirty flag is checked when reset_request() is 
     *              called to see if memset() on requests is necessary.
     *
     * Return: void
     *
     * Parameters: void
     */
    void set_request_dirty();

    /*
     * Description:	Set dirty-from-neighbor flag to indicate that a neighbor
     *              has requested ECM fragments from this ECM manager.
     *              Should be called whenever this ECM manager submits a 
     *              request for fragments to a neighboring patch as follows:
     * 					
     *  this->ECMWorldPtr->worldECM[neighbor_index].set_dirty_from_neighbors();
     *
     *              Dirty-from-neighbor flag is checked when updateECM() is 
     *              called to see if there's a fragment request to be processed
     *
     * Return: void
     *
     * Parameters: void
     */
    void set_dirty_from_neighbors();

    /*
     * Description:	Reset dirty flag. Should be called after done using this 
     *              flag in this tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void reset_dirty();

    /*
     * Description:	Reset request-dirty flag. Should be called after done using
     *              this flag in this tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void reset_request_dirty();

    /*
     * Description:	Reset dirty-from-neighbor flag. Should be called after done
     *              using this flag in this tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void reset_dirty_from_neighbors();
#endif

#ifdef OPT_FECM_REQUEST

    void addCollagenFragment();
    void addElastinFragment();
    void addHAFragment();
    void addHAFragments(int nfr);
#endif

    bool killHA(int iHA);
    bool killHA(HL_ITERATOR itHA);

    void lock();
    void unlock();

    /*
     * Description:	Decreases the input parameter by 1. 
     *              Used to decrease the life of ECM hyaluronan proteins.
     *
     * Return: void
     *
     * Parameters: n  -- An integer for the number of lives for hyaluronan.
     */
    static void decrement(int n);


    /*************************************************************************
     * CONSTANT VARIABLES                                                    *
     *************************************************************************/
    // Used to store the ECM's location in x,y,z dimensions of world
    int indice[3];
    // Used to store the patch row major index of the ECM manager.
    int index;

    /*************************************************************************
     * NON-CONSTANT VARIABLES                                                *
     *************************************************************************/
    // Whether there is ECM or not at the beginning and end of each tick
    bool empty[2];
    /* The number of original collagen (tropocollagen monomer), new collagen(collagen) and fragmented collagen at
     * the beginning and end of each tick
     */
    int ocollagen[2], ncollagen[2], fcollagen[2];
    /* The number of original elastin (monomer), new elastin and fragmented elastin at
     * the beginning and end of each tick
     */
    int oelastin[2], nelastin[2], felastin[2];
    /* The number of hyaluronan and fragmented hyaluronan at the beginning and
     * end of each tick
     */ 
//    int HA[2];
//    int fHA[2]; 
#ifdef OPT_FECM_REQUEST
    int fcollagenrequest;
    int felastinrequest;
    int fHArequest;
#else
    // Keeps track of all neighbors' requests for fragmented collagen.
    int requestfcollagen[27];
    // Keeps track of all neighbors' requests for fragmented elastin.
    int requestfelastin[27];
    // Keeps track of all neighbors' requests for fragmented hyaluronan.
    int requestfHA[27];
#endif	// OPT_FECM_REQUEST

#ifdef _OMP
    // Lock to this ECM object
    omp_lock_t ECMmutex;
#endif
    /* Whether or not there is a fragmented collagen signalling danger at the
     * beginning and end of each tick 
     */
    bool fcDangerSignal[2];
    /* Whether or not there is a fragmented elastin signalling danger at the
     * beginning and end of each tick 
     */
    bool feDangerSignal[2];
    /* Whether or not there is a fragmented hyaluronan signalling danger at the
     * beginning and end of each tick 
     */
    bool fHADangerSignal[2]; 
    /* Whether there is enough collagen to form a scar at the beginning and end
     * of each tick
     */
    bool scarIndex[2];

#ifdef OPT_ECM
    // Whether there has been a change in this ECM's attributes    
    bool dirty;
    // Whether this ECM manager has requested ECM fragments on a neighbor
    bool request_dirty;
    // Whether a neighbor has requested ECM fragments from this ECM manager
    bool dirty_from_neighbors;
#endif

    /*************************************************************************
     * STATIC VARIABLES                                                      *
     *************************************************************************/
    // Pointer to array of all patches in the world
    static Patch* ECMPatchPtr; 
    // Pointer to the wound healing world
    static WHWorld* ECMWorldPtr;
    // The maximum amount of collagen (of all types) allowed
    static int maxcollagen;
    // The maximum amount of elastin (of all types) allowed
    static int maxelastin;
    // The maximum amount of hyaluronan (of all types) allowed
    static int maxHA; 
    // x dimension displacement to each neighbor
    static int dx[27];
    // y dimension displacement to each neighbor
    static int dy[27];
    // z dimension displacement to each neighbor
    static int dz[27];
    // Index of each neighbor in dx, dy, dz arrays
    static int d[27];
    // Index of each neighbor in dx, dy, dz, with first 9 entry being middle Z plane
    static int nid[27];
private:
    int HA[2];
    int fHA[2];
#ifdef OPT_ECM
    /* Number of lives for each hyaluronan (HA) remaining in this tick.
     * With OPT_ECM defined, we assume that we access HAlife in this order:
     * 1. Decrement life by calling decrement(int n) on HAlife in ECMFunction()
     * 2. Determine number of HA by calling HAlife.size() in fragmentHA()
     * 3. Remove dead HAs in HAlife in updateECM()
     * NOTE: Steps 1 & 2 can occur in either order but must preceed step 3.
     */
    std::vector<int> *HAlife;
#else
    // Number of lives for hyaluronan at the beginning and end of each tick.
    std::vector<int> *HAlife[2];
#endif
};

#endif	/* ECM_H */

