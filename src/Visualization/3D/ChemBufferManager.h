#ifndef CHEM_BUFFER_MANAGER
#define CHEM_BUFFER_MANAGER

#ifdef ECV_SAMPLE_CHEM


extern "C" void initCudaChemSample(cudaExtent volumeSize, int chemIndex);

extern "C" void sampleChem(
		  float *d_Src,
	    int dataD,
	    int dataH,
	    int dataW,
	    int chemIndex
);

#endif	// ECV_SAMPLE_CHEM

#endif	// CHEM_BUFFER_MANAGER
