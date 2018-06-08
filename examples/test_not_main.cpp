/* 
 * File: test_not_main.cpp
 *
 * File Contents: Contains main method of model and various output functions.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

// Include C/C++ libraries
#include <cstdlib>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <fstream>
#include <cstring>
#include <sstream>
#include <time.h>
#include <sys/time.h>

//Include local libraries
#include "../src/Driver/Driver.h"
#include "../src/Utilities/output_utils.h"
#include "../src/World/World.h"
#include "../src/World/Usr_World/woundHealingWorld.h"
#include "../src/Agent/Agent.h"
#include "../src/Agent/Usr_Agents/Platelet.h"
#include "../src/Agent/Usr_Agents/Fibroblast.h"
#include "../src/FieldVariable/FieldVariable.h"
#include "../src/Patch/Patch.h"
#include "../src/ECM/ECM.h"


#include "../src/enums.h"
#include "../src/FieldVariable/Usr_FieldVariables/WHChemical.h"


#ifdef VISUALIZATION

#ifdef MODEL_3D
#include "../src/Visualization/3D/Visualization.h"
#endif	// MODEL_3D

#endif	// VISUALIZATION

#include "../src/Utilities/parameters.h"

#include "../src/Utilities/input_utils.h"
using namespace std;

/*
 * Description:	Outputs the color of each patch to the given file.
 *              Outputs the color of the agent that is one the patch and
 *              if there is no agent, then outputs the color of the patch type.
 *              Also prints the number of patches that are tissue, epithelium, 
 *              blood, fibroblast, new fibroblast, platelet, macrophage,
 *              neutrophil, black, and the total number of patches.
 *
 * Return: 0 if succesful
 *
 * Parameters: WHWorld*  -- Pointer to the wound healing world whose patches'
 *                          colors should be outputted.
 *             char*     -- Output file name
 */
//int outputColor(WHWorld*,char*);

/*
 * Description:	Outputs the color of each patch to the given file.
 *              Outputs the color of the ECM that is one the patch and
 *              if there is no ECM, then outputs the color of the patch type.
 *
 * Return: 0 if succesful
 *
 * Parameters: myWorld   -- Pointer to the wound healing world whose patches'
 *                          colors should be outputted.
 *             fileName  -- Output file name
 */
//int outputECM(WHWorld* myWorld, char* fileName);

/*
 * Description:	Outputs the chemical concentration of the given chemical on
 *              each patch to the given file.
 *
 * Return: 0 if succesful
 *
 * Parameters: WHWorld*  -- Pointer to the wound healing world whose patches'
 *                          colors should be outputted.
 *             char*     -- Output file name
 *             int       -- Enumic value of the chemical type to output
 */
//int outputChem(WHWorld*, char*, int);

/*
 * Description:	Main method for model simulation. Sets up the world and
 *              executes each tick of the simulation. Calculates cell insertion
 *              and deletion statistics as well as execution time information.
 *
 * Return: 0 if succesful
 *
 * Parameters: argc  -- Number of command-line arguments
 *             argv  -- Command-line arguments
 */
int main(int argc, char** argv) {

  /********************************
   * SET UP A WOUND HEALING WORLD *
   ********************************/
  // Get baseline cell and chemical values that user specified
	util::processOptions(argc, argv);
	util::printOptions();
	util::processParameters("Sample.txt");
	clock_t tStart = clock();
	WHWorld::setSeed(util::getSeed());
	WHWorld myWorld = WHWorld(util::getWorldXWidth(), util::getWorldYWidth(), util::getWorldZWidth(), util::getPatchWidth());
	printf("Setup Execution time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);


  /****************************************************
   * EXECUTE EACH TICK (30 MINUTES EACH) OF THE MODEL *
   ****************************************************/
	int numTicks = util::getNumTicks();

	Driver modelDriver = Driver(&myWorld);
	modelDriver.init(argc, argv);
	modelDriver.run(numTicks);

	return 0;
}

