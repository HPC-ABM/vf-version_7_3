
#ifndef DRIVER_H_

#define DRIVER_H_


#include "../Utilities/output_utils.h"
#include "../Utilities/parameters.h"
#include "../common.h"
#include "../World/Usr_World/woundHealingWorld.h"

#ifdef VISUALIZATION

#ifdef MODEL_3D
#include "../Visualization/3D/Visualization.h"
#endif  // MODEL_3D

#endif  // VISUALIZATION

class World;
class WHWorld;

using namespace std;

class Driver
{
public:
  Driver(WHWorld *worldPtr);
  virtual ~Driver();

  void run(int ticks);
  void init(int argc, char ** argv);
private:
  WHWorld *worldPtr;

};

#endif	// DRIVER_H_
