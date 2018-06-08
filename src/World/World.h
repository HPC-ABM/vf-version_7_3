/*
 * File: World.h 
 *
 * File Contents: Contains declarations for the World class
 *   
 * Author: alireza
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *
 * Created on May 19, 2013, 7:41 AM
 *****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#ifndef WORLD_H
#define	WORLD_H

#include "../Agent/Agent.h"
#include "../FieldVariable/FieldVariable.h"
#include "../Patch/Patch.h"

#include "../common_inc.h"

#include <stdlib.h>
#include <vector>

using namespace std;
//extern class dev_Patch;

/*
 * WORLD CLASS DESCRIPTION:
 * The World class manages all worlds in the model. 
 * It is used to set up the dimensions of the world.
 */
class World {
 public:
    /*
     * Description:	Default World constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    World();

    /*
     * Description:	World constructor. 
     *
     * Return: void
     *
     * Parameters: orig  -- Pointer to an original World
     */
    World(const World& orig);

    /*
     * Description:	Virtual World destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    virtual ~World();

    /*
     * Description:	Initializes the grid size and dimensions
     *
     * Return: void
     *
     * Parameters: nx, ny, nz    -- Number of grid points (patches) 
     *                               in x,y,z dimensions
     *             x_min, x_max  -- Min and max coordinates in x
     *             y_min, y_max  -- Min and max coordinates in y
     *             z_min, z_max  -- Min and max coordinates in z
     */
    void setupGrid(int nx, int ny,int nz, REAL x_min, REAL x_max, REAL y_min, REAL y_max, REAL z_min, REAL z_max); //!< set up the dimensions

    /*
     * Description:	Write vtk output file "filename" for animation
     *
     * Return: void
     *
     * Parameters: filename  -- Path to output file
     *             t         -- time
     */
    void outputWorld_VTK_binary(const char* filename, REAL t);

    // Number of grid points (lattices) (patches) on this world in x,y,z
    int nx, ny, nz;
    // Number of patches
    int np;
    // World volume in mL
    float worldVmL;
    // Mesh (patch) size in x,y,z dimenstions
    REAL dx, dy, dz;
    // Max and min coordinate in x 
    REAL x_min, x_max;
    // Max and min coordinate in y 
    REAL y_min, y_max;
    // Max and min coordinate in z 
    REAL z_min, z_max;
#ifndef MIN_MEM
    // Spatial coordinates of the grid
    vector<REAL> x,y,z;
#endif
    // Two dimensional array of field variables
    vector<vector<REAL> > field_var;            

    // For generating random numbers
    unsigned seed;
};

#endif	/* WORLD_H */
