#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#pragma once

#include "../World/Usr_World/woundHealingWorld.h"

#include <string>

class World;
class WHWorld;

using namespace std;
namespace util {

int threadSafeRand(unsigned *seed_arr);


float randFloatRange(float a, float b, unsigned *seed_arr);

bool divisible(long double a, long double b);

bool divisible(float a, float b);

}

#endif
