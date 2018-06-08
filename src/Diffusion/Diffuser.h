/*
 * Diffuser.h
 *
 *  Created on: Oct 16, 2016
 *      Author: nseekhao
 */

#ifndef DIFFUSER_H_
#define DIFFUSER_H_

#include "../common_inc.h"
#include "../common.h"

#include "convolutionFFT_common.h"
#include "DiffusionHelper.h"

#include <unistd.h>
#include <stdlib.h>
#include <vector>

using namespace std;

class Diffuser {
public:
    /*
     * Description: Default Diffuser constructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    Diffuser();

    /*
     * Description: Diffuser constructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    Diffuser(int nBaseChem, c_ctx *cc, c_ctx *kc, WHChemical *WHWorldChem);

    /*
     * Description: Virtual Diffuser destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    virtual ~Diffuser();

#ifdef GPU_DIFFUSE

void diffuseChemGPU(int rightWallIndex);

#else   // GPU_DIFFUSE

void diffuseChemCPU(float dt);

#endif  // GPU_DIFFUSE

private:
    int nBaseChem;
    c_ctx *chem_cctx;
    c_ctx *kernel_cctx;
    WHChemical *WHWorldChem;

};


#endif /* DIFFUSER_H_ */
