/* 
 * File: Agent.h
 *
 * File Contents: Contains declarations for the Agent class.
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 * 
 * Created on May 19, 2013, 7:37 AM
 */
/*****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#ifndef AGENT_H
#define	AGENT_H

#include <stdlib.h>
#include <algorithm>
#include <vector>

#include "../common_inc.h"

#include "../common.h"
#include "../enums.h"

//#include "../World/World.h"

class World;
class WHWorld;
class Patch;
class ECM; 

using std::vector;

/*
 * AGENT CLASS DESCRIPTION:
 * Agent is a parent class of all cell "agents". 
 * It is used to move agents on 2D patches, to make agents die, 
 * and to make them move to highest chemical concentration.
 */
class Agent {
 public:
    /*
     * Description:	Default agent constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    Agent();

    /*
     * Description:	Virtual agent destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    virtual ~Agent();

    // Enumic type to keep track of the type of agent.
    enum agenttype_t {platelet, fibroblast, macrophage, neutrophil}; 
    // Enumic type to keep track of the type of chemical.
    enum chemtype_t {TNF, TGF, FGF, MMP8, IL1beta, IL6, IL8, IL10}; 

    /*
     * Description:	Determines whether an agent is alive or not and updates the
     *              agent's properties for the end of the current tick.
     *
     * Return: True if the agent is alive, false otherwise.
     *
     * Parameters: void
     */
    bool isAlive();

    /*
     * Description:	Determines the location of the agent in the x dimension of 
     *              the world.
     *
     * Return: x-coordinate of agent
     *
     * Parameters: void
     */
    int getX();

    /*
     * Description:	Determines the location of the agent in the y dimension of 
     *              the world.
     *
     * Return: y-coordinate of agent
     *
     * Parameters: void
     */
    int getY();

    /*
     * Description:	Determines the location of the agent in the z dimension of 
     *              the world.
     *
     * Return: z-coordinate of agent
     *
     * Parameters: void
     */
    int getZ();

    /*
     * Description:	Determines the patch row major index of the agent.
     *
     * Return: patch row major index of agent
     *
     * Parameters: void
     */
    int getIndex();

    /*
     * Description:	Rolls a hypothetical dice with 'percent' chance of 
     *              successful roll
     *
     * Return: True if a successful roll, false otherwise
     *
     * Parameters: percent  -- percentage that dictates the chance of a 
     *                         successful roll
     */
    static bool rollDice(int percent);
//    static bool rollDice(float percent);

    /*
     * Description:	Moves the cell to the target patch if it is valid and 
     *              available. Updates cell location and patch availability.
     *
     * Returns: True if the move was successful, false otherwise.
     * 
     * Parameters: dX         -- distance to move in the x direction
     *                           (-: left, +: right)
     *				     dY         -- distance to move in the y direction
     *                           (-: up, +: down)
     *				     dZ         -- distance to move in the z direction
     *                           (-: inwards, +: outwards)
     *				     readIndex  -- index into cell's current location array
     *	          						   NOTE: If the location has been modified in
     *                           this tick, readIndex should be write_t
     *                           Else, readIndex should be read_t
     *
     * Usage: This function should be called after a call to check whether the
     *        location is dirty (modified in this tick)
     *        Example call sequence:
     *				    if (isModified(this->index)) {
     *					    move (dX, dY, dZ, write_t);
     *				    } else {
     *					    move (dX, dY, dZ, read_t);
     *				    }
     *        The check is to make sure that we are reading from the most 
     *        up-to-date information of the cell. This check is necessary,
     *        since move() can be called multiple times in a tick, thus could
     * 		    potentially modify each attribute more than once. The 
     *        intermediate values need to be kept track of.
     */
    bool move(int dX, int dY, int dZ, int readIndex);

    /*
     * Description:	Determines the mean concentration of chemical of type 
     *              chemIndex from 27 neighboring patches
     *
     * Return: mean concentration of chemical of type chemIndex from 27 
     *         neighboring patches
     *
     * Parameters: chemIndex  -- enumic value of chemical_t for chemical type
     */
    float meanNeighborChem(int chemIndex);

    /*
     * Description:	Determines the number of cells of type cellIndex from 27
     *              neighboring patches.
     *
     * Return: number of cells of type cellIndex from 27 neighboring patches.
     *
     * Parameters: cellIndex  -- enumic value of agent_t for cell type
     */
    int countNeighborCells(int cellIndex);

    /*
     * Description:	Determines the number of ECM proteins of type ECMIndex from
     *              27 neighboring patches.
     *
     * Return: number of ECM proteins of type ECMIndex from 27 neighboring
     *         patches.
     *
     * Parameters: ECMIndex  -- enumic value of agent_t for ECM protein type
     */
    float countNeighborECM(int ECMIndex);

    /*
     * Description:	Moves to neighboring patch with highest concentration of
     *              chemical of type chemIndex.
     *
     * Return: True if moved successfully, false otherwise.
     *         neighboring patches
     *
     * Parameters: chemIndex    -- enumic value of chemical_t for chemical type
     *             isGrad       -- flag indicating whether or not to move to a gradient
     */
    bool moveToHighestChem(int chemIndex, bool isGrad);

    /*
     * Description: Moves to neighboring patch with highest concentration of
     *              chemical of type chemIndex.
     *
     * Return: True if moved successfully, false otherwise.
     *         neighboring patches
     *
     * Parameters: chemIndex    -- enumic value of chemical_t for chemical type
     */
    bool moveToHighestChem(int chemIndex);

    /*
     * Description:	Move to a random neighboring patch while respecting rules
     *                  dictating movement between tissue types.
     *                  wiggle() ONLY moves onto a healhty neighbor
     *
     * Return: True if moved successfully, false otherwise.
     *         neighboring patches
     *
     * Parameters: chemIndex  -- enumic value of chemical_t for chemical type
     */
    void wiggle();

    void wiggle_wound();

    bool moveToTarget(int dX, int dY, int dZ, int index);

    /*
     * Description:	Updates class members for the end of the tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateAgent();

    /*
     * Description:	Virtual function for cell function
     *
     * Return: void
     *
     * Parameters: void
     */
    virtual void cellFunction();

    /*
     * Description:	Virtual function for agent death
     *
     * Return: void
     *
     * Parameters: void
     */
    virtual void die() = 0;

    /* 
     * Description: Virtual function for copying and initializing a new agent
     *
     * Return: void
     *
     * Parameters: original  -- Agent to be copied
     *             dx        -- Difference in x-coordinate between current and
     *                          new agents.
     *             dy        -- Difference in y-coordinate between current and 
     *                          new agents.
     *             dz        -- Difference in z-coordinate between current and
     *                          new agents.
     *                          NOTE: dz = 0 because it is only 2D for now. TODO(Nuttiiya): I'm guessing this needs to change when you implement 3D?
    */
    virtual void copyAndInitialize(Agent* original, int dx, int dy, int dz = 0);

    bool isActivated();
    void printCytRelease(int celltype, int chemtype, int x, int y, int z, float level);

    // Number of lives remaining at the beginning and end of each tick
    int life[2];
    // Whether agent is activated or not at the beginning and end of each tick
    bool activate[2];
    // Agent's color at the beginning and end of each tick
    int color[2];
    // Agent's size at the beginning and end of each tick
    float size[2];
    // Agent's type at the beginning and end of each tick
    int type[2];

    // Pointer from an agent to a WHWorld    
    static WHWorld* agentWorldPtr;
    // Pointer from an agent to a Patch
    static Patch* agentPatchPtr;
    // Pointer from an agent to an ECM
    static ECM* agentECMPtr;
    // Number of patches in x,y,z dimensions of the world
    static int nx, ny, nz;
    // Difference in x-,y-,z-coordinates to neighboring patches
    static int dX[27], dY[27], dZ[27];
#ifdef MODEL_3D
    // Array of neighboring patches in 3D model
    static int neighbor[27];
#else
    // Array of neighboring patches in 2D model
    static int neighbor[8];
#endif

 protected:
    // Agent position in x,y,z dimensions at the beginning and end of each tick
    int ix[2],iy[2],iz[2];
    // Patch row major index for agent at the beginning and end of each tick
    int index[2];
    // Life status of agent at the beginning and end of each tick
    bool alive[2];
};

#endif	/* AGENT_H */

