
#include "DiffusionHelper.h"


#ifdef GPU_DIFFUSE

void diffusion_helper::findGPUs()
{
#if (N_GPU > 1)
    int totalGPUs;
    checkCudaErrors(cudaGetDeviceCount(&totalGPUs));

    if (N_GPU > totalGPUs){
        printf("Error: Only have %d GPUs. Require %d.\n", totalGPUs, N_GPU);
        exit(-1);
    }

    //Print the device information to run the code
    for (int i = 0 ; i < N_GPU ; i++)
    {
        cudaDeviceProp deviceProp;
        checkCudaErrors(cudaGetDeviceProperties(&deviceProp, i));
        printf("GPU Device %d: \"%s\" with compute capability %d.%d\n", i, deviceProp.name, deviceProp.major, deviceProp.minor);

    }
#endif  // N_GPU
}

void diffusion_helper::initKernelDimensions(c_ctx *kernel_cctx)
{
    /********************************************
     * GPU IDs Initialization                   *
     ********************************************/
    for (int ic = 0; ic < N_CHEM; ic++)
    {
        kernel_cctx->gpu_id[ic] = ic % N_GPU;
        // DEBUG
        if (kernel_cctx->gpu_id[ic] == 1) kernel_cctx->dev_id[ic] = GPU2;
        else kernel_cctx->dev_id[ic] = GPU1;
        printf("chemical %d is assigned to GPU %d\n", ic, kernel_cctx->gpu_id[ic]);
    }

    /********************************************
     * Dimension calculation                    *
     ********************************************/

#ifdef MODEL_3D	//------
    int   totalSec  = 1800;//1800;//900;//1800;

#ifdef RAT_VF
    int   nSec      =  200;//150;			// Does not work for > 200
#else	// RAT_VF
    int   nSec      =  300;//300;//16;//180;//900;
#endif	// RAT_VF

    int   k_niter   = totalSec/nSec;
    float dt		= kernel_cctx->dt;	// 2.0;

    if (totalSec % nSec)
    {
        printf("ERROR: totalSec (%d) NOT integer multiple of nSec (%d)\n", totalSec, nSec);
    }
#else   //------ MODEL_3D
    int nSec    = 1800;//90;//180;//600;
    int k_niter = -1;
    // TODO: calculate dt w.r.t. max D (for now assume max D = 20.0)
    float dt = 2.4;//2.5;
#endif  //------ MODEL_3D

    int kernelRadius            = nSec/dt;

    // DEBUG rat
#ifdef RAT_VF
    kernelRadius = nSec;
    const int outKernelH        = (kernelRadius-1) * 2 + 1;
    const int outKernelW        = outKernelH;
          int outKernelD        = outKernelH;

//    kernelRadius = 204;
//    const int outKernelH        = 509;
//    const int outKernelW        = 509;
//          int outKernelD        = 509;
#else	// RAT_VF
    const int outKernelH        = (kernelRadius-1) * 2 + 1;
    const int outKernelW        = outKernelH;
          int outKernelD        = outKernelH;

//          printf(":::::: DEBUG: nsec: %d, kernelradius: %d\n", nSec, kernelRadius);
#endif	// RAT_VF

#ifdef MODEL_3D
#ifdef HUMAN_VF

    const int    kfftD = 616;
    const int    kfftH = 616;
    const int    kfftW = 616;

    const int windowD = 147;
    const int windowH = 147;
    const int windowW = 147;
#elif defined(RAT_VF)

    const int    kfftD = snapTransformSize(outKernelD + KCoeffsD - 1);
    const int    kfftH = snapTransformSize(outKernelH + KCoeffsH - 1);
    const int    kfftW = snapTransformSize(outKernelW + KCoeffsW - 1);
//    const int    kfftD = 512;	// 512 is the closest power of 2
//    const int    kfftH = 512;	// 	to KCoeffsD + nSec (base kernel size) - 1
//    const int    kfftW = 512;	// 	or 3 + 300 - 1

    // DEBUG rat
    const int windowD = outKernelD / 2;//(kernelRadius/2) + 1;//283;//285;	// (1.0 mm/0.007 mm)*2 + 1;
    const int windowH = outKernelH / 2;//(kernelRadius/2) + 1;//283;//285;	// (1.0 mm/0.007 mm)*2 + 1;
    const int windowW = outKernelW / 2;//(kernelRadius/2) + 1;//283;// 57;	// (0.2 mm/0.007 mm)*2 + 1;

#else

    const int    kfftD = snapTransformSize(outKernelD + KCoeffsD - 1);
    const int    kfftH = snapTransformSize(outKernelH + KCoeffsH - 1);
    const int    kfftW = snapTransformSize(outKernelW + KCoeffsW - 1);

    const int windowD = (kernelRadius/2) + 1;
    const int windowH = (kernelRadius/2) + 1;
    const int windowW = (kernelRadius/2) + 1;
#endif  // HUMAN_VF
#else  // MODEL_3D
    outKernelD        = 1;

    const int    kfftH = snapTransformSize(outKernelH + KCoeffsH - 1);
    const int    kfftW = snapTransformSize(outKernelW + KCoeffsW - 1);
    const int    kfftD = 1;

    const int windowH = (kernelRadius/2) * 2 + 1; //outKernelH/2;
    const int windowW = (kernelRadius/2) * 2 + 1; //outKernelW/2;
    const int windowD = 1;
#endif  // MODEL_3D

    /********************************************
     * Set dimension                            *
     ********************************************/
    kernel_cctx->niter          = k_niter;
    kernel_cctx->dt             = dt;
    kernel_cctx->kernelRadius   = kernelRadius;

    kernel_cctx->DD     = outKernelD;
    kernel_cctx->DH     = outKernelH;
    kernel_cctx->DW     = outKernelW;

    kernel_cctx->KD     = KCoeffsD;
    kernel_cctx->KH     = KCoeffsH;
    kernel_cctx->KW     = KCoeffsW;

    kernel_cctx->KZ     = 1;//0;
    kernel_cctx->KX     = 1;//0;
    kernel_cctx->KY     = 1;//0;

    kernel_cctx->FFTD   = kfftD;
    kernel_cctx->FFTH   = kfftH;
    kernel_cctx->FFTW   = kfftW;

    kernel_cctx->windowD     = windowD;
    kernel_cctx->windowH     = windowH;
    kernel_cctx->windowW     = windowW;
}


void diffusion_helper::allocKernelBuffers(c_ctx *kernel_cctx)
{
    int windowD = kernel_cctx->windowD;
    int windowH = kernel_cctx->windowH;
    int windowW = kernel_cctx->windowW;

    int kfftD = kernel_cctx->FFTD;
    int kfftH = kernel_cctx->FFTH;
    int kfftW = kernel_cctx->FFTW;

    int *gpu_id = kernel_cctx->gpu_id;

    /********************************************
     * Kernel Buffer Allocations                *
     ********************************************/

    /*
     *      h_ibuffs            -- NULL (only used for diffusion NOT kernel computation)
     *      h_obuffs            -- NULL (only used for diffusion NOT kernel computation)
     *      h_kernelspectrum    -- NULL (only used for diffusion NOT kernel computation)
     *
     * 3D:
     *      d_data              -- Point to one buffer per device
     *      h_data              -- Point to one buffer per device
     *      d_kernelspectrum_h   -- Point to one buffer per device
     *
     * 2D:
     *      d_data              -- Point to one buffer per chemical type
     *      h_data              -- Point to one buffer per chemical type
     *      d_kernelspectrum_h    -- Point to one buffer per chemical type
     *
     */

    // Input, output buffers, and h_kernelspectrum are
    // for chemical diffusion NOT kernel computation

    for (int ic = 0; ic < N_CHEM; ic++)
    {
        kernel_cctx->h_ibuffs[ic] = NULL;
        kernel_cctx->h_obuffs[ic] = NULL;

        kernel_cctx->h_kernelspectrum[ic] = NULL;
    }

    int k_wsize_b = windowD * windowH * windowW * sizeof(float);
    int fftsize_b = kfftD * kfftH * (kfftW / 2 + 1) * sizeof(fComplex);

#ifdef MODEL_3D
    // Allocate one buffer for each of the following per device:
    //  - h_data
    //  - d_data
    //  - d_kernelspectrum
    for (int ig = 0; ig < N_GPU; ig++) {
        kernel_cctx->h_data[ig]	= (float *) malloc(k_wsize_b);

	// DEBUG
	int devID;
	if (ig == 0) devID = GPU1;
	else devID = GPU2;
        checkCudaErrors(cudaSetDevice(devID));//ig));
        //DEBUG
        printf("------ before alloc kernel buffers gpu %d ------\n", ig);
        reportMemUsageGPU();

        checkCudaErrors(cudaMalloc((void **)&kernel_cctx->d_data[ig],
                k_wsize_b));
        checkCudaErrors(cudaMalloc((void **)&kernel_cctx->d_kernelspectrum_h[ig], fftsize_b));
        //DEBUG
        printf("------ 1 after alloc kernel buffers gpu %d ------\n", ig);
        reportMemUsageGPU();
    }

    for (int ic = N_GPU; ic < N_CHEM; ic++){
        int ig = gpu_id[ic];
	// DEBUG
	int devID;
	if (ig == 0) devID = GPU1;
	else devID = GPU2;
        checkCudaErrors(cudaSetDevice(devID));//ig));

        kernel_cctx->h_data[ic]              = kernel_cctx->h_data[ig];
        kernel_cctx->d_data[ic]              = kernel_cctx->d_data[ig];
        kernel_cctx->d_kernelspectrum_h[ic]  = kernel_cctx->d_kernelspectrum_h[ig];

        //DEBUG
        printf("------ 2 after alloc kernel buffers gpu %d ------\n", ig);
        reportMemUsageGPU();
    }

#else	// MODEL_3D
    for (int ic = 0; ic < N_CHEM; ic++){
        // HOST
        kernel_cctx->h_data[ic]  = (float *) malloc(k_wsize_b);

        // DEVICE
        checkCudaErrors(cudaMalloc((void **)&kernel_cctx->d_data[ic],
                k_wsize_b));
        checkCudaErrors(cudaMalloc((void **)&kernel_cctx->d_kernelspectrum_h[ic], fftsize_b));
    }
#endif	// MODEL_3D

}

void diffusion_helper::deallocConvCtxBuffers(c_ctx *ctx)
{
    for (int ic = 0; ic < 8; ic++)
    {
#if defined(OPT_CHEM) || defined(PAGELOCK_DIFFUSE)
        if (ctx->h_ibuffs[ic]) checkCudaErrors(cudaFreeHost(ctx->h_ibuffs[ic]));
#else   // OPT_CHEM
        if (ctx->h_ibuffs[ic]) free(ctx->h_ibuffs[ic]);
#endif// OPT_CHEM
        // reset pointer
        ctx->h_ibuffs[ic] = NULL;
    }


    for (int ic = 0; ic < 8; ic++)
    {
#if defined(OPT_CHEM) || defined(PAGELOCK_DIFFUSE)
        if (ctx->h_obuffs[ic]) checkCudaErrors(cudaFreeHost(ctx->h_obuffs[ic]));
#else   // OPT_CHEM
        if (ctx->h_obuffs[ic]) free(ctx->h_obuffs[ic]);
#endif// OPT_CHEM
        // reset pointer
        ctx->h_obuffs[ic] = NULL;
    }


    for (int ic = 0; ic < 8; ic++)
    {
#ifdef  PAGELOCK_DIFFUSE
        if (ctx->h_kernelspectrum[ic]) checkCudaErrors(cudaFreeHost(ctx->h_kernelspectrum[ic]));
#else	//  PAGELOCK_DIFFUSE 
        if (ctx->h_kernelspectrum[ic]) free(ctx->h_kernelspectrum[ic]);
#endif	//  PAGELOCK_DIFFUSE
        // reset pointer
        ctx->h_kernelspectrum[ic] = NULL;
    }


    int *gpu_id = ctx->gpu_id;
    for (int ic = 0; ic < 8; ic++)
    {
        int ig = gpu_id[ic];
        if (ctx->d_kernelspectrum_h[ig]) checkCudaErrors(cudaFree(ctx->d_kernelspectrum_h[ig]));
        // reset pointer
        ctx->d_kernelspectrum_h[ic] = NULL;
    }

}


void diffusion_helper::printContext(c_ctx c)
{
	cout << "dt"  << "	" << c.dt << endl;
	cout << "dx2" << "	" << c.dx2 << endl;
	cout << "kernelRadius" << "	" << c.kernelRadius << endl;
	cout << "niter" << "	" << c.niter << endl << endl;

	cout << "KD" << "	" << c.KD << endl;
	cout << "KH" << "	" << c.KH << endl;
	cout << "KW" << "	" << c.KW << endl << endl;

	cout << "KZ" << "	" << c.KZ << endl;
	cout << "KX" << "	" << c.KX << endl;
	cout << "KY" << "	" << c.KY << endl << endl;

	cout << "DD" << "	" << c.DD << endl;
	cout << "DH" << "	" << c.DH << endl;
	cout << "DW" << "	" << c.DW << endl << endl;

	cout << "FFTD" << "	" << c.FFTD << endl;
	cout << "FFTH" << "	" << c.FFTH << endl;
	cout << "FFTW" << "	" << c.FFTW << endl << endl;

	cout << "windowD" << "	" << c.windowD << endl;
	cout << "windowH" << "	" << c.windowH << endl;
	cout << "windowW" << "	" << c.windowW << endl << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "h_ibuffs[" << i << "]" << c.h_ibuffs[i] << endl;
	cout << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "h_obuffs[" << i << "]" << c.h_obuffs[i] << endl;
	cout << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "gpu_id[" << i << "]" << c.gpu_id[i] << endl;
	cout << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "dev_id[" << i << "]" << c.dev_id[i] << endl;
	cout << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "d_data[" << i << "]" << c.d_data[i] << endl;
	cout << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "h_data[" << i << "]" << c.h_data[i] << endl;
	cout << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "h_kernelspectrum[" << i << "]" << c.h_kernelspectrum[i] << endl;
	cout << endl;

	for (int i = 0; i < N_CHEM; i++)
		cout << "d_kernelspectrum_h[" << i << "]" << c.d_kernelspectrum_h[i] << endl;
	cout << endl;

}


void diffusion_helper::initChemDimensions(c_ctx *chem_cctx, c_ctx kernel_cctx, int nx, int ny, int nz)
{
    /********************************************
     * GPU IDs Initialization                   *
     ********************************************/
    for (int ic = 0; ic < N_CHEM; ic++)
    {
        chem_cctx->gpu_id[ic] = ic % N_GPU;
	// DEBUG
	if (chem_cctx->gpu_id[ic] == 1) chem_cctx->dev_id[ic] = GPU2;
	else chem_cctx->dev_id[ic] = GPU1;
    }
    /********************************************
     * Dimension calculation                    *
     ********************************************/
    const int    windowW = kernel_cctx.windowW;
    const int    windowH = kernel_cctx.windowH;
    const int    windowD = kernel_cctx.windowD;

    const int    chemH = ny;
    const int    chemW = nx;
    const int    chemD = nz;

#ifdef MODEL_3D
#ifdef HUMAN_VF
    const int    cfftD = 1152;
    const int    cfftW =  256;
    const int    cfftH = 1536;
#else   // HUMAN_VF
    const int    cfftD = snapTransformSize(chemD + kernel_cctx.windowD - 1);
    const int    cfftW = snapTransformSize(chemW + kernel_cctx.windowW - 1);
    const int    cfftH = snapTransformSize(chemH + kernel_cctx.windowH - 1);
#endif  // HUMAN_VF
#else   // MODEL_3D
    const int    cfftD = 1;
    const int    cfftH = snapTransformSize(chemH + kernel_cctx.windowH - 1);
    const int    cfftW = snapTransformSize(chemW + kernel_cctx.windowW - 1);
#endif  // MODEL_3D

    int rc = 0;

    if (cfftD - windowD + 1 < chemD) {
        printf("InitializeChemGPU Error: bad dimensions\n\tcfftD: %d\twindowD: %d\tchemD: %d\n"
                , cfftD, windowD, chemD);
        rc = -1;
    }
    if (cfftH - windowH + 1 < chemH) {
        printf("InitializeChemGPU Error: bad dimensions\n\tcfftH: %d\twindowH: %d\tchemH: %d\n"
                , cfftH, windowH, chemH);
        rc = -1;
    }
    if (cfftW - windowW + 1 < chemW) {
        printf("InitializeChemGPU Error: bad dimensions\n\tcfftW: %d\twindowW: %d\tchemW: %d\n"
                , cfftW, windowW, chemW);
        rc = -1;
    }
    if (rc == -1) exit(rc);
    /********************************************
     * Set dimension                            *
     ********************************************/
    chem_cctx->DD   = chemD;
    chem_cctx->DH   = chemH;
    chem_cctx->DW   = chemW;

    chem_cctx->KD   = kernel_cctx.windowD;
    chem_cctx->KH   = kernel_cctx.windowH;
    chem_cctx->KW   = kernel_cctx.windowW;

    int isEvenX			= (kernel_cctx.windowW % 2)? 0 : 1;
    int isEvenY			= (kernel_cctx.windowH % 2)? 0 : 1;
    int isEvenZ			= (kernel_cctx.windowD % 2)? 0 : 1;

    chem_cctx->KZ   = (kernel_cctx.windowD / 2) - isEvenX;
    chem_cctx->KX   = (kernel_cctx.windowW / 2) - isEvenY;
    chem_cctx->KY   = (kernel_cctx.windowH / 2) - isEvenZ;

    chem_cctx->FFTD = cfftD;
    chem_cctx->FFTH = cfftH;
    chem_cctx->FFTW = cfftW;

}


void diffusion_helper::allocChemBuffers(c_ctx *chem_cctx, WHChemical *WHWorldChem, int nx, int ny, int nz)
{
    int cfftD = chem_cctx->FFTD;
    int cfftH = chem_cctx->FFTH;
    int cfftW = chem_cctx->FFTW;

    int fftsizeb  = cfftD * cfftH * (cfftW / 2 + 1) * sizeof(fComplex);
    int datasizeb = nx*ny*nz * sizeof(float);

    /********************************************
     * Diffusion Buffer Allocations             *
     ********************************************/

    /*
     *      h_ibuffs            -- Point to one buffer per chemical type
     *      h_obuffs            -- Point to one buffer per chemical type
     *      h_kernelspectrum    -- Point to one buffer per chemical type
     *      h_data              -- NULL  (only used for kernel computation NOT diffusion )
     *
     * 3D:
     *      d_data              -- Point to one buffer per device
     *      d_kernelspectrum_h   -- Point to one buffer per device
     *
     * 2D:
     *      d_data              -- Point to one buffer per chemical type
     *      d_kernelspectrum_h   -- Point to one buffer per chemical type
     *
     */

    // System host buffer h_data is for chemical kernel computation
    for (int ic = 0; ic < N_CHEM; ic++)
        chem_cctx->h_data[ic] = NULL;

    // Pinned host buffer h_ibuffs is used here for chemical data input
    //  (which gets copied to d_data).
    // However, if OPT_CHEM is not defined, then we set h_ibuffs to pointer to patch chem.
#ifdef OPT_PINNED_MEM
    // Pinned host memory
    checkCudaErrors(cudaMallocHost((void**)&(chem_cctx->h_obuffs[0]), datasizeb));
    checkCudaErrors(cudaMallocHost((void**)&(chem_cctx->h_ibuffs[0]), datasizeb));
#endif	// OPT_PINNED_MEM
    for (int ic = 0; ic < N_CHEM; ic++)
    {
#ifdef OPT_CHEM
#ifdef OPT_PINNED_MEM
        if (ic)
        {
            chem_cctx->h_obuffs[ic] = chem_cctx->h_obuffs[0];
            chem_cctx->h_ibuffs[ic] = chem_cctx->h_ibuffs[0];
        }
#else	// OPT_PINNED_MEM
        // Pinned host memory
        checkCudaErrors(cudaMallocHost((void**)&(chem_cctx->h_obuffs[ic]), datasizeb));
        checkCudaErrors(cudaMallocHost((void**)&(chem_cctx->h_ibuffs[ic]), datasizeb));
                //, cudaHostAllocWriteCombined));
#endif	// OPT_PINNED_MEM
#else   // OPT_CHEM
        // System host memory
        chem_cctx->h_obuffs[ic] = WHWorldChem->getTchemPtr(ic);
        chem_cctx->h_ibuffs[ic] = WHWorldChem->getPchemPtr(ic);
#endif  // OPT_CHEM
    }


#ifdef MODEL_3D

    int *gpu_id = chem_cctx->gpu_id;

    /**
     * For index 0 to N_GPU-1, Allocate:
     *  - d_data
     *  - d_kernelspectrum_h
     *  - h_kernelspectrum
     */
    for (int ig = 0; ig < N_GPU; ig++) {
        checkCudaErrors(cudaSetDevice(chem_cctx->dev_id[ig]));//ig));
        //DEBUG
        printf("------ gpu %d ------\n", ig);
        reportMemUsageGPU();

        printf("Allocating chem data buffer on GPU device %d\n", ig);
        // DEVICE
        checkCudaErrors(cudaMalloc((void **)&(chem_cctx->d_data[ig]),           datasizeb));
        checkCudaErrors(cudaMalloc((void **)&(chem_cctx->d_kernelspectrum_h[ig]), fftsizeb));

        // HOST
#ifdef PAGELOCK_DIFFUSE
        checkCudaErrors(cudaHostAlloc(
                (void**)&(chem_cctx->h_kernelspectrum[ig]), fftsizeb, cudaHostAllocWriteCombined));
#else   // PAGELOCK_DIFFUSE
        chem_cctx->h_kernelspectrum[ig] = (fComplex *) malloc(fftsizeb);
#endif  // PAGELOCK_DIFFUSE
    }

    /**
     * For index N_GPU to N_CHEM-1, Allocate:
     *  - h_kernelspectrum
     *
     * For index N_GPU to N_CHEM-1, copy allocated pointer values from index 0 to N_GPU:
     *  - d_data
     *  - d_data
     */
    for (int ic = N_GPU; ic < N_CHEM; ic++){
        int ig = gpu_id[ic];
        checkCudaErrors(cudaSetDevice(chem_cctx->dev_id[ic]));//ig));

        printf("\tCopying pointer from chem_cctx[%d] to chem_cctx[%d]\t%p\t%p\n",
		ig, ic, chem_cctx->d_data[ig], chem_cctx->d_kernelspectrum_h[ig]);

        // DEVICE
        chem_cctx->d_data[ic]             = chem_cctx->d_data[ig];
#ifndef M40
        chem_cctx->d_kernelspectrum_h[ic] = chem_cctx->d_kernelspectrum_h[ig];
#else	// !M40
        checkCudaErrors(cudaMalloc((void **)&(chem_cctx->d_kernelspectrum_h[ic]), fftsizeb));	
#endif	// !M40

        // HOST
#ifdef PAGELOCK_DIFFUSE
        checkCudaErrors(cudaHostAlloc(
                (void**)&(chem_cctx->h_kernelspectrum[ic]),fftsizeb, cudaHostAllocWriteCombined));
#else	// PAGELOCK_DIFFUSE
        chem_cctx->h_kernelspectrum[ic] = (fComplex *) malloc(fftsizeb);
#endif	// PAGELOCK_DIFFUSE
    }


#else	//MODEL_3D

    for (int ic = 0; ic < N_CHEM; ic++){
        // DEVICE
        checkCudaErrors(cudaMalloc((void **)&(chem_cctx->d_data[ic]),           datasizeb));
        checkCudaErrors(cudaMalloc((void **)&(chem_cctx->d_kernelspectrum_h[ic]), fftsizeb));

        // HOST
#ifdef PAGELOCK_DIFFUSE
        checkCudaErrors(cudaHostAlloc(
                (void**)&(chem_cctx->h_kernelspectrum[ic]), fftsizeb, cudaHostAllocWriteCombined));
#else   // PAGELOCK_DIFFUSE
        chem_cctx->h_kernelspectrum[ic] = (fComplex *) malloc(fftsizeb);
#endif  // PAGELOCK_DIFFUSE
    }

#endif	//MODEL_3D
}


void diffusion_helper::prepAndComputeKernelSpectrum(
        c_ctx* chem_cctx, c_ctx kernel_cctx, float ** h_dWindow)
{
#ifdef MODEL_3D
    /*
     * Computing kernel spectrum is performed in computeKernel3DBatch(), which
     * is called in prepAndComputeKernel()
     */
#else   // MODEL_3D
    /********************************************
     * Kernel Spectrum Computation              *
     ********************************************/
    const int    chemH = chem_cctx->DH;
    const int    chemW = chem_cctx->DW;

    const int    cfftH = snapTransformSize(chemH + kernel_cctx.windowH - 1);
    const int    cfftW = snapTransformSize(chemW + kernel_cctx.windowW - 1);

    printf("Chem world: %d x %d\n", chemH, chemW);

    StopWatchInterface *hTimer = NULL;
    sdkCreateTimer(&hTimer);

    sdkResetTimer(&hTimer);
    sdkStartTimer(&hTimer);

    // COMPUTE KERNEL SPECTRUM
    for (int ic = 0; ic < N_CHEM; ic++) {
        computeKernelSpectrum(
                chem_cctx->d_kernelspectrum_h[ic],    // output
                h_dWindow[ic],                      // input
                kernel_cctx,
                *chem_cctx
        );
    }

    sdkStopTimer(&hTimer);
    double kernelSpectrumComputationTime = sdkGetTimerValue(&hTimer);
    printf("\tTotal kernel computation: %f MPix/s (%f ms)\n",
            (double)chemH * (double)chemW * 1e-6 / (kernelSpectrumComputationTime * 0.001),
            kernelSpectrumComputationTime);

#endif  // MODEL_3D
}


void diffusion_helper::prepKernel(c_ctx* chem_cctx, c_ctx kernel_cctx, float *lambda, float *gamma, int numChem)
{
    const int chemW = chem_cctx->DW;
    const int chemH = chem_cctx->DH;
    const int chemD = chem_cctx->DD;

    const int windowW = kernel_cctx.windowW;
    const int windowH = kernel_cctx.windowH;
    const int windowD = kernel_cctx.windowD;

    const int kernelRadius = kernel_cctx.kernelRadius;

    float dt = kernel_cctx.dt;

#ifdef MODEL_3D

    int c_dsize_b   = chemD     * chemH    * chemW * sizeof(float);
    int *gpu_id = kernel_cctx.gpu_id;

    // Compute 3D kernel + spectrum in batch
    computeKernel3DBatch(
            kernelRadius,
            lambda,
            gamma,
            dt,
            kernel_cctx,
            *chem_cctx
    );



#else   // MODEL_3D

    /*
     * |----------- SYSTEM HOST MEMORY -----------|
     * |                                          |
     * |---------- |                              |  |---------- DEVICE MEMORY ----------|
     * | h_dWindow | ---------> | dev pointer w0 ||--|-------->  |f|f|...|f|f|           |
     * |---------- |            | dev pointer w1 ||--|------     |...|...|...|           |
     * |                        |      ...       ||  |     |     |f|f|...|f|f|           |
     * |                        | dev pointer w7 ||--|---- |                             |
     * |                                          |  |   | --->  |f|f|...|f|f|           |
     * |                                          |  |   |       |...|...|...|           |
     * |                                          |  |   |       |f|f|...|f|f|           |
     * |                                          |  |   |                               |
     * |                                          |  |   ----->  |f|f|...|f|f|           |
     * |                                          |  |           |...|...|...|           |
     * |                                          |  |           |f|f|...|f|f|           |
     */
    float
    **h_dWindow;

    // Allocate array of 'numChem' pointers for device kernel center windows
    h_dWindow = (float **) malloc (numChem * sizeof(float*));

    // Allocate 'numChem' arrays for device kernels
    for (int ic = 0; ic < numChem; ic++) {
#ifdef PRINT_KERNEL
        cout << "Allocating kernel arrays" << endl;
#endif
        checkCudaErrors(cudaMalloc((void **)&h_dWindow[ic], windowH * windowW * sizeof(float)));
#ifdef PRINT_KERNEL
        cout << "Computing kernel arrays " << ic << endl;
#endif
        // KERNEL COMPUTATION
        if (!computeKernel(
                h_dWindow[ic],
                kernelRadius,
                lambda[ic],
                gamma[ic],
                dt,
                kernel_cctx))
        {
            // TODO: Error Handling
        }
    }

    /********************************************
     * Kernel Spectrum Computation              *
     ********************************************/
    prepAndComputeKernelSpectrum(chem_cctx, kernel_cctx, h_dWindow);

#ifdef PRINT_KERNEL
        cout << "Freeing device pointers ..." << endl;
#endif
    // Free the rest of the pointers
    if (h_dWindow)
    {
        for (int ic = 0; ic < numChem; ic++) {
            checkCudaErrors(cudaFree(h_dWindow[ic]));
        }
        free(h_dWindow);   h_dWindow = NULL;
    }
#ifdef PRINT_KERNEL
        cout << "Finished freeing device pointers" << endl;
#endif

#endif  // MODEL_3D
}


#endif	// GPU_DIFFUSE






