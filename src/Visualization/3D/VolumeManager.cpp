
#include "VolumeManager.h"

#ifdef VISUALIZATION

template<typename T>
VolumeManager::VolumeManager(int x, int y, int z, T* volumeDataPtr): x(x), y(y), z(z)
{
    if (!volumeDataPtr)
        this->volumeBuffer = new T[x*y*z];
}

VolumeManager::~VolumeManager()
{
    if (this->volumeBuffer) delete [] this->volumeBuffer;
}


// TODO: Move this to child class of ECM volume
// Read content change, amplify and superimpose on scaled initial value
void VolumeManager::preprocessECMdata(const float *ecm, const float *dEcm)
{

//    if (util::ABMerror(
//            !this->volumeBuffer,
//            "Error pre-processing ECM data!!",
//            __FILE__,
//            __LINE__))
//        exit(2);
    if (!this->volumeBuffer)
    {
        std::cerr << "VolumeManager::preprocessECMdata: NULL buffer" << std::endl;
        exit(2);
    }

    const int dataCount = this->x * this->y * this->z;
    unsigned char buff[dataCount];
    for (int iz = 0; iz < z; iz++)
        for (int iy = 0; iy < y; iy++)
            for (int ix = 0; ix < x; ix++)
            {
                int index = iz * this->y * this->x + iy * this->x + ix;
                float f = RAW_ECM_FACTOR*ecm[index] + RAW_dECM_FACTOR*dEcm[index];
                this->volumeBuffer[index] = (f >= 1.0 ? 255 : (f <= 0.0 ? 0 : (int)floor(f * 256.0)));
            }
}

#endif  // VISUALIZATION
