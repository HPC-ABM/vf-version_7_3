
#include "Driver.h"

Driver::Driver(WHWorld *worldPtr)
{
  this->worldPtr = worldPtr;
}

Driver::~Driver(){}


void Driver::init(int argc, char ** argv)
{
#ifdef VISUALIZATION
  visualization::init(argc, argv, worldPtr);
#endif
}

void Driver::run(int numTicks)
{

  if (!worldPtr)
  {
    printf("Error: Driver.run() null worldPtr\n");
    exit(-1);
  }

  struct timeval start, end;      // Timing information
  long elapsed_times[numTicks];   // in milliseconds 

#ifdef VISUALIZATION

/*****************************************/
#ifdef OVERLAP_VIS

#pragma omp parallel num_threads(2)
{
  int tid = omp_get_thread_num();
  if (tid == 0) {
    // Render
    cout << "Driver(): start vis" << endl;
    visualization::start();
  } else {
    // Compute
    cout << "Driver(): start go()" << endl;
    for (int tick = 0; tick < numTicks; tick++) {
     /**********************
      * TIMING INFORMATION *
      **********************/
      gettimeofday(&start, NULL);
      worldPtr->go();
      gettimeofday(&end, NULL);

      // Store execution time for this tick
      elapsed_times[tick] = (end.tv_sec*1000 + end.tv_usec/1000) -
                              (start.tv_sec*1000 + start.tv_usec/1000);
      cout << "       this go() execution took " << elapsed_times[tick] << " ms" << endl;
    }

   /**********************
    * TIMING INFORMATION *
    **********************/
    int acc = 0;
    for (int i = 0; i < numTicks; i++) {
            cout << "Tick " << i << " took " << elapsed_times[i] << " ms" << endl;
            acc += elapsed_times[i];
    }
    cout << "Average time (excluding iter 0):  "
         << (acc - elapsed_times[0])/(numTicks - 1)
         << " ms" << endl;
    exit(-1);
  }

}
#elif defined(INTERACTIVE_VIS)
#pragma omp parallel num_threads(2)
{
  int tid = omp_get_thread_num();
  if (tid == 0) {
    // Render
    cout << "Driver(): start vis" << endl;
    visualization::start();
  } else {
    // Compute
    cout << "Driver(): start go()" << endl;
    int tick = 0;
    while (!visualization::isFinished()) {
    	if (visualization::isPaused()) continue;

     /**********************
      * TIMING INFORMATION *
      **********************/
      gettimeofday(&start, NULL);
      worldPtr->go();
      gettimeofday(&end, NULL);

      // Store execution time for this tick
      elapsed_times[tick%numTicks] = (end.tv_sec*1000 + end.tv_usec/1000) -
                              (start.tv_sec*1000 + start.tv_usec/1000);
      cout << "       this go() execution took " << elapsed_times[tick%numTicks] << " ms" << endl;
      tick++;
    }

   /**********************
    * TIMING INFORMATION *
    **********************/
    int acc = 0;
    for (int i = 0; i < numTicks; i++) {
            cout << "Tick " << i << " took " << elapsed_times[i] << " ms" << endl;
            acc += elapsed_times[i];
    }
    cout << "Average time (excluding iter 0):  "
         << (acc - elapsed_times[0])/(numTicks - 1)
         << " ms" << endl;
    exit(-1);
  }

}
#else 	// OVERLAP_VIS
  cout << "+++++++VIS+++++++" << endl;
  visualization::start();
#endif	// OVERLAP_VIS
/****************************************/

#else	// VISUALIZATION
 
  for (int tick = 0; tick < numTicks; tick++)
  {
    cout << "+++++++NO  VIS+++++++" << endl;
#ifdef PARAVIEW_RENDERING
    // Prepare output filenames and output current state of wound healing world
    char simulation[50] = "output/Simulation/5daysimulation_";
    char ECMsim[50] = "output/Simulation/ECM_";
    char chem0Sim[50] = "output/Simulation/chem0_";
    char extension[10] = ".vtk";
    char tempNumber[20] = "";
    sprintf(tempNumber, "%d", tick); // %d makes the result be a signed decimal integer
    strcat(simulation, tempNumber);
    strcat(simulation, extension);
    strcat(ECMsim, tempNumber);
    strcat(ECMsim, extension);

    // Ouputs the cytokine values:
    char ptnf[50] = "output/Simulation/tnf_";
    char ptgf[50] = "output/Simulation/tgf_";
    char pfgf[50] = "output/Simulation/fgf_";
    char pil1[50] = "output/Simulation/il1_";
    char pil6[50] = "output/Simulation/il6_";
    char pil8[50] = "output/Simulation/il8_";
    char pil10[50] = "output/Simulation/il10_";
    char pmmp[50] = "output/Simulation/mmp_";
    char pmac[50] = "output/Simulation/pmacgrad_";
    char pneu[50] = "output/Simulation/pneugrad_";
    char pfib[50] = "output/Simulation/pfibgrad_";
    sprintf (tempNumber, "%d", tick); // %d makes the result be a decimal integer
    strcat(ptnf, tempNumber); strcat(ptgf, tempNumber); strcat(pfgf, tempNumber);
    strcat(pil1, tempNumber); strcat(pil6, tempNumber); strcat(pil8, tempNumber);
    strcat(pil10, tempNumber); strcat(pmmp, tempNumber); strcat(pmac, tempNumber);
    strcat(pneu, tempNumber); strcat(pfib, tempNumber);
    strcat(ptnf, extension); strcat(ptgf, extension); strcat(pfgf, extension);
    strcat(pil1, extension); strcat(pil6, extension); strcat(pil8, extension);
    strcat(pil10, extension);strcat(pmmp, extension); strcat(pmac, extension);
    strcat(pneu, extension); strcat(pfib, extension);
    util::outputColor(worldPtr, simulation);
    util::outputECM(worldPtr, ECMsim);
    util::outputChem(worldPtr, ptnf, TNF);
    util::outputChem(worldPtr, ptgf, TGF);
    util::outputChem(worldPtr, pfgf, FGF);
    util::outputChem(worldPtr, pil1, IL1beta);
    util::outputChem(worldPtr, pil6, IL6);
    util::outputChem(worldPtr, pil8, IL8);
    util::outputChem(worldPtr, pil10, IL10);
    util::outputChem(worldPtr, pmmp, MMP8);
//    util::outputChem(worldPtr, pmac, pmacgrad);
//    util::outputChem(worldPtr, pneu, pneugrad);
//    util::outputChem(worldPtr, pfib, pfibgrad);
#endif

#ifdef BIOMARKER_OUTPUT
    util::outputWorld_csv(worldPtr);
#endif

#ifdef WRITE_RAW_ECM
    if (tick % RAW_TIME_STRIDE == 0)
    {
    	char collfile[100];
    	sprintf(collfile, "output/Volume/Collagen_%d.raw", tick);
    	util::writeNormalizedFloatRaw(
    	        collfile,
    	        worldPtr->ecmMap[m_col],
    	        worldPtr->nx,
    	        worldPtr->ny,
    	        worldPtr->nz);

//    	util::writeNormalizedFloatRawFrames("output/Volume/Collagen_10frames.raw",
//                worldPtr->ecmMap[m_col],
//                worldPtr->nx,
//                worldPtr->ny,
//                worldPtr->nz);

    	char dcollfile[100];
    	sprintf(dcollfile, "output/Volume/dCollagen_%d.raw", tick);
    	util::writeSuperPositionedECMgrad(
                dcollfile,
                worldPtr->ecmMap[m_col],
                worldPtr->ecmGradMap[m_col],
                worldPtr->nx,
                worldPtr->ny,
                worldPtr->nz);

    	// DEBUG collagen
        char fcollfile[100];
        sprintf(fcollfile, "output/Volume/Collagen_%d.csv", tick);
        util::writeFloatCsv(fcollfile,
                worldPtr->ecmMap[m_col],
                worldPtr->nx,
                worldPtr->ny,
                worldPtr->nz);
    }
#endif

#ifdef WRITE_RAW_CELL
    if (tick % RAW_TIME_STRIDE == 0)
    {
        util::writeCellRaw(tick,
                worldPtr,
                worldPtr->nx,
                worldPtr->ny,
                worldPtr->nz);
    }
#endif

    // Run the simulation for 1 tick (30 min)
    clock_t t1 = clock();
    cerr << "executing go() at tick " << tick << " ..." << endl;
    cout << "entering go() at tick " << tick << endl;
    gettimeofday(&start, NULL);
    worldPtr->go();
    gettimeofday(&end, NULL);

    // Store execution time for this tick
    elapsed_times[tick] = (end.tv_sec*1000 + end.tv_usec/1000) - (start.tv_sec*1000 + start.tv_usec/1000);
    cout << "       this go() execution took " << elapsed_times[tick] << " ms" << endl;
  }

#ifdef CALIBRATION
    // Used for sensitivity analysis
    if (worldPtr->treatmentOption == voicerest) {

            util::outputTotalChem(worldPtr, "output/SensitivityAnalysis/FinalTotalChemVR.dat");

    } else if (worldPtr->treatmentOption == resonantvoice){

            util::outputTotalChem(worldPtr, "output/SensitivityAnalysis/FinalTotalChemRV.dat");

    } else if (worldPtr->treatmentOption == spontaneousspeech){

            util::outputTotalChem(worldPtr, "output/SensitivityAnalysis/FinalTotalChemSS.dat");
    }
#endif

#ifdef COLLECT_CELL_INS_DEL_STATS
   int addedSum = 0;
   int erasedSum = 0;
   int addedMax = worldPtr->numAddedCells[0];
   int erasedMax = worldPtr->numErasedCells[0];
   int addedMin = worldPtr->numAddedCells[0];
   int erasedMin = worldPtr->numErasedCells[0];
   int addedMin_NoZero = worldPtr->numAddedCells[0];
   int erasedMin_NoZero = worldPtr->numErasedCells[0];

  /*****************************
   * CELL INSERTION STATISTICS *
   *****************************/
   cout << "       Cell Insertion Stats: " << endl;
   for (int i = 0; i < numTicks; i++) {
           cout << worldPtr->numAddedCells[i] << endl;
           addedSum += worldPtr->numAddedCells[i];
           addedMax = worldPtr->numAddedCells[i] > addedMax? worldPtr->numAddedCells[i] : addedMax;
           addedMin = worldPtr->numAddedCells[i] < addedMin? worldPtr->numAddedCells[i] : addedMin;
           addedMin_NoZero = (worldPtr->numAddedCells[i] < addedMin) && (worldPtr->numAddedCells[i] > 0)?
                           worldPtr->numAddedCells[i] : addedMin;
   }
   cout << "               Min:            " << addedMin << endl;
   cout << "               Min (no 0):     " << erasedMin << endl;
   cout << "               Max:            " << addedMax << endl;
   cout << "               Average:        " << addedSum/numTicks << endl;

  /****************************
   * CELL DELETION STATISTICS *
   ****************************/
   cout << endl;
   cout << "       Cell Deletion Stats: " << endl;
   for (int i = 0; i < numTicks; i++) {
           cout << worldPtr->numErasedCells[i] << endl;
           erasedSum += worldPtr->numErasedCells[i];
           erasedMax = worldPtr->numErasedCells[i] > addedMax? worldPtr->numErasedCells[i] : addedMax;
           erasedMin = worldPtr->numErasedCells[i] < addedMin? worldPtr->numErasedCells[i] : addedMin;
           erasedMin_NoZero = (worldPtr->numErasedCells[i] < addedMin) && (worldPtr->numErasedCells[i] > 0)? worldPtr->numErasedCells[i] : addedMin;
   }
   cout << "               Min:            " << erasedMin << endl;
   cout << "               Min (no 0):     " << erasedMin << endl;
   cout << "               Max:            " << erasedMax << endl;
   cout << "               Average:        " << erasedSum/numTicks << endl;
#endif

  
  /**********************
   * TIMING INFORMATION *
   **********************/
   int acc = 0;
   for (int i = 0; i < numTicks; i++) {
           cout << "Tick " << i << " took " << elapsed_times[i] << " ms" << endl;
           acc += elapsed_times[i];
   }
   cout << "Average time:  " << acc/numTicks << " ms" << endl;

#endif	// VISUALIZATION

}



