
#ifndef OUTPUT_UTILS_H_
#define OUTPUT_UTILS_H_

#pragma once

#include "../World/Usr_World/woundHealingWorld.h"
#include "../Patch/Patch.h"

#include <string>

class World;
class WHWorld;
class Patch;

using namespace std;
namespace util {


//writes raw file
bool writeRawFile(const char *fileName, const char *data, const int dataCount){
    using namespace std;
    ofstream out(fileName,ios::binary);
    if(!out)
        return false;

    out.write((char *)data, dataCount);
    out.close();
    return true;
}

bool appendRawFile(const char *fileName, const char *data, const int dataCount){
    using namespace std;
    ofstream out(fileName,ios::binary|ios::app);
    if(!out)
        return false;

    out.write((char *)data, dataCount);
    out.close();
    return true;
}



bool writeFloatCsv(const char *fileName, const float *data,
        int width, int height, int depth)
{
    using namespace std;
    ofstream out(fileName);
    if(!out)
        return false;

    for (int z = 0; z < depth; z++)
    {
        for (int y = 0; y < height; y++)
        {
            for (int x = 0; x < width; x++)
            {
                int index = z * height * width + y * width + x;
                out << data[index] << ", ";
            }
            out << endl;
        }
        out << endl;
    }

    out.close();
    return true;

}

bool writeSuperPositionedECMgrad(const char *fileName, const float *ecm, const float *dEcm,
        int width, int height, int depth)
{

    const int dataCount = width * height * depth;
    unsigned char buff[dataCount];
    for (int z = 0; z < depth; z++)
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                int index = z * height * width + y * width + x;
                float f = RAW_ECM_FACTOR*ecm[index] + RAW_dECM_FACTOR*dEcm[index];
                buff[index] = (f >= 1.0 ? 255 : (f <= 0.0 ? 0 : (int)floor(f * 256.0)));
            }

    return writeRawFile(fileName, (char *)buff, dataCount);

}

bool writeNormalizedFloatRaw(const char *fileName, const float *data,
		int width, int height, int depth)
{

	const int dataCount = width * height * depth;
	unsigned char buff[dataCount];
	for (int z = 0; z < depth; z++)
		for (int y = 0; y < height; y++)
			for (int x = 0; x < width; x++)
			{
				int index = z * height * width + y * width + x;
				float f = data[index];
				buff[index] = (f >= 1.0 ? 255 : (f <= 0.0 ? 0 : (int)floor(f * 256.0)));
			}

	return writeRawFile(fileName, (char *)buff, dataCount);

}

bool writeNormalizedFloatRawFrames(const char *fileName, const float *data,
        int width, int height, int depth)
{

    const int dataCount = width * height * depth;
    unsigned char buff[dataCount];
    for (int z = 0; z < depth; z++)
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                int index = z * height * width + y * width + x;
                float f = data[index];
                buff[index] = (f >= 1.0 ? 255 : (f <= 0.0 ? 0 : (int)floor(f * 256.0)));
            }

    return appendRawFile(fileName, (char *)buff, dataCount);

}

bool writeCellRaw(int tick, WHWorld *worldPtr,
        int width, int height, int depth)
{

    bool rc = true;
    char activated = 10;
    char inactivated = 100;
    short numcelltypes = 3;
    const int dataCount = width * height * depth;
    unsigned char buff[numcelltypes][dataCount];
    for (int z = 0; z < depth; z++)
        for (int y = 0; y < height; y++)
            for (int x = 0; x < width; x++)
            {
                int index = z * height * width + y * width + x;
                int ctype = worldPtr->worldPatch[index].occupiedby[read_t];
                switch (ctype)
                {
                case neutrophil:
                    buff[0][index] = inactivated;
                    break;
                case aneutrophil:
                    buff[0][index] = activated;
                    break;
                case macrophag:
                    buff[1][index] = inactivated;
                    break;
                case amacrophag:
                    buff[1][index] = activated;
                    break;
                case fibroblast:
                    buff[2][index] = inactivated;
                    break;
                case afibroblast:
                    buff[2][index] = activated;
                    break;
                default:
                    break;
                }

            }

    char filenames[numcelltypes][100];

    for (int i = 0; i < numcelltypes; i++)
    {
        sprintf(filenames[i], "output/Volume/cell%d_%d.raw", i, tick);
        rc &= writeRawFile(filenames[i], (char *)buff, dataCount);
    }

    return rc;

}


void outputWorld_csv(WHWorld *world) {

    int nx = world->nx;
    int ny = world->ny;
    int nz = world->nz;

	char outfilename[100];
#ifdef MULT_SEEDS
	sprintf(outfilename, "output/random/Output_Biomarkers_s%d.csv", seed);
#else
	sprintf(outfilename, "output/Output_Biomarkers.csv");
#endif
	if (world->clock == 0) {
#ifndef MULT_SEEDS
		remove(outfilename);
#endif
		ofstream output_file(outfilename, ios::app);
		output_file << "clock" << "," << "TNF" << "," << "TGF" << "," << "FGF" << "," << "Il1beta" << "," << "IL6" << "," << "IL8" << ",";
		output_file << "IL10" << "," << "MMP8" << "," << "Tropocollagen" << "," << "Collagen" << "," << "FragentedCollagen" << "," << "Tropoelastin" << "," << "Elastin" << ",";
		output_file << "FragmentedElastin" << "," << "HA" << "," << "FragmentedHA" << "," << "Damage" << "," << "ActivatedFibroblast" << ",";
		output_file << "ActivatedMacrophage" << "," << "ActivatedNeutrophil" << "," << "Fibroblast" << "," << "Macrophage" << "," << "Neutrophil" << "," << "Platelet" << endl;
		output_file.close();
	}

	ofstream output_file(outfilename, ios::app);
	int f = 0; int af = 0; int m = 0; int am = 0; int n = 0; int an = 0; int plt = 0;
	float oc = 0.0f; float fc = 0.0f; float nc = 0.0f;
	float oe = 0.0f; float fe = 0.0f; float ne = 0.0f;
	float HA = 0.0f; float fHA = 0.0f;
	cout << " fibs " <<  world->fibs.actualSize() << endl;
	cout << " macs " <<  world->macs.actualSize() << endl;
	cout << " neus " <<  world->neus.actualSize() << endl;
	cout << " plats " << world->plats.actualSize() << endl;

	int fibsSize = world->fibs.size();
	for (int i = 0; i < fibsSize; i++) {
		Fibroblast* fib = world->fibs.getDataAt(i);
		if (!fib) continue;
		if (fib->isAlive() == false) continue;
		if (fib->activate[read_t] == false) f++;
		else af++;
	}


	int macsSize = world->macs.size();
	for (int i = 0; i<macsSize; i++) {
		Macrophage* mac = world->macs.getDataAt(i);
		if (!mac) continue;
		if (mac->isAlive() == false) continue;
		if (mac->activate[read_t] == false) m++;
		else am++;
	}
	int neusSize = world->neus.size();
	for (int i = 0; i<neusSize; i++){
		Neutrophil* neu = world->neus.getDataAt(i);
		if (!neu) continue;
		if (neu->isAlive() == false) continue;
		if (neu->activate[read_t] == false) n++;
		else an++;
	}
	plt = world->plats.actualSize(); //Agent::agentWorldPtr->plats.actualSize();

//	for (int in = 0; in < (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny; in++) {
//	for (int in = 0; in < nx*ny*nz; in++) {
//
//		oc += world->worldECM[in].ocollagen[read_t];
//		nc += world->worldECM[in].ncollagen[read_t];
//		fc += world->worldECM[in].fcollagen[read_t];
//		oe += world->worldECM[in].oelastin[read_t];
//		ne += world->worldECM[in].nelastin[read_t];
//		fe += world->worldECM[in].felastin[read_t];
//		HA += world->worldECM[in].getnHA();//HA[read_t];
//		fHA += world->worldECM[in].getnfHA();//fHA[read_t];
//	}

	world->updateStats();
	nc = world->totalNC;
	oc = world->totalOC;
	fc = world->totalFC;

    ne = world->totalNE;
    oe = world->totalOE;
    fe = world->totalFE;

    HA  = world->totalHA;
    fHA = world->totalFHA;

	//	world->countPatchType(damage);


	cout << "  nc  " << nc << endl;
	cout << "  ne  " << ne << endl;
	cout << "  HA  " << HA << endl;
	cout << "  fHA " << fHA << endl;
	output_file << world->clock << ",";
	float worldVmL = world->worldVmL;
	output_file << world->WHWorldChem->total[TNF]/worldVmL << ",";
	output_file << world->WHWorldChem->total[TGF]/worldVmL << ",";
	output_file << world->WHWorldChem->total[FGF]/worldVmL << ",";
	output_file << world->WHWorldChem->total[IL1beta]/worldVmL << ",";
	output_file << world->WHWorldChem->total[IL6]/worldVmL << ",";
	output_file << world->WHWorldChem->total[IL8]/worldVmL << ",";
	output_file << world->WHWorldChem->total[IL10]/worldVmL << ",";
	output_file << world->WHWorldChem->total[MMP8]/worldVmL << ",";
	output_file << oc << "," << nc << "," << fc << "," << oe << "," ;
	output_file << ne << "," << fe << "," << HA << "," << fHA << "," << Patch::numOfEachTypes[pdamage] << "," ;
	output_file << af << "," << am << "," << an << "," << f << "," << m  << ","  << n << "," << plt << endl;

	output_file.close();

}



void patchassign_csv(WHWorld *world) {

    int nx = world->nx;
    int ny = world->ny;
    int nz = world->nz;

	int in = 0;
	int Number = 0;
	for (int iz = 0; iz < nz; iz++) {
		char patchassign[50] = "output/patchassign";
		char cells[50] = "output/cells_read";
		char cells_w[50] = "output/cells_write";
		char initHA[50] = "output/initHA";
		char initcollagen[50] = "output/initcollagen";
		char initelastin[50] = "output/initelastin";
		char damagezone[50] = "output/damagezone";
		char initialdamage[50] = "output/initialdamage";
		char extension[10] = ".csv";
		char tempNumber[20] = "";
		sprintf (tempNumber, "_t%3.0f_z%d", world->clock, Number);

		strcat(patchassign, tempNumber);
		strcat(patchassign, extension);
		strcat(cells, tempNumber);
		strcat(cells, extension);
		strcat(cells_w, tempNumber);
		strcat(cells_w, extension);
		strcat(initHA, tempNumber);
		strcat(initHA, extension);
		strcat(initcollagen, tempNumber);
		strcat(initcollagen, extension);
		strcat(initelastin, tempNumber);
		strcat(initelastin, extension);
		strcat(damagezone, tempNumber);
		strcat(damagezone, extension);
		strcat(initialdamage, tempNumber);
		strcat(initialdamage, extension);

		// Patch Assign
		ofstream output_file(patchassign, ios::app);

		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (world->worldPatch[in].type[read_t] == pdamage ||
				        world->worldPatch[in].damage[read_t] != 0) {
					output_file << "x";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == nothing) {
					output_file << "-";
				}
				if (world->worldPatch[in].type[read_t] == tissue) {
					output_file << "t";
				}
				if (world->worldPatch[in].type[read_t] == epithelium) {
					output_file << "e";
				}
				if (world->worldPatch[in].type[read_t] == capillary) {
					output_file << "o";
				}
				if (world->worldPatch[in].type[read_t] == unidentifiable) {
					output_file << "?";
				}
			}
			output_file << endl;
		}
		output_file.close();

		// initHA
		ofstream output_file1(initHA, ios::app);

		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (world->worldPatch[in].initHA == true) {
					output_file1 << "u";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == pdamage
				        || world->worldPatch[in].damage[read_t] != 0) {
					output_file1 << "x";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == nothing) {
					output_file1 << "-";
				}
				if (world->worldPatch[in].type[read_t] == tissue) {
					output_file1 << "t";
				}
				if (world->worldPatch[in].type[read_t] == epithelium) {
					output_file1 << "e";
				}
				if (world->worldPatch[in].type[read_t] == capillary) {
					output_file1 << "o";
				}
				if (world->worldPatch[in].type[read_t] == unidentifiable) {
					output_file1 << "?";
				}
			}
			output_file1 << endl;
		}

		output_file1.close();

		// Damage Zone
		ofstream output_file2(damagezone, ios::app);

		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;

				if (world->worldPatch[in].inDamzone == true) {
					output_file2 << "z";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == pdamage
				        || world->worldPatch[in].damage[read_t] != 0) {
					output_file2 << "x";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == nothing) {
					output_file2 << "-";
				}
				if (world->worldPatch[in].type[read_t] == tissue) {
					output_file2 << "t";
				}
				if (world->worldPatch[in].type[read_t] == epithelium) {
					output_file2 << "e";
				}
				if (world->worldPatch[in].type[read_t] == capillary) {
					output_file2 << "o";
				}
				if (world->worldPatch[in].type[read_t] == unidentifiable) {
					output_file2 << "?";
				}
			}
			output_file2 << endl;
		}
		output_file2.close();

		// Initial Damage
		ofstream output_file3(initialdamage, ios::app);

		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (world->worldECM[in].monomer_ela[read_t] !=0 &&
				        world->worldPatch[in].damage[read_t] != 0) {
					output_file3 << "g";
					continue;
				}
				if (world->worldECM[in].monomer_ela[read_t] !=0) {
					output_file3 << "m";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == pdamage ||
				        world->worldPatch[in].damage[read_t] != 0) {
					output_file3 << "x";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == nothing) {
					output_file3 << "-";
				}
				if (world->worldPatch[in].type[read_t] == tissue) {
					output_file3 << "t";
				}
				if (world->worldPatch[in].type[read_t] == epithelium) {
					output_file3<<"e";
				}
				if (world->worldPatch[in].type[read_t] == capillary) {
					output_file3 << "o";
				}
				if (world->worldPatch[in].type[read_t] == unidentifiable) {
					output_file3 << "?";
				}
			}
			output_file3 << endl;
		}
		output_file3.close();

		// initCollagen
		ofstream output_file4(initcollagen, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (world->worldECM[in].monomer_col[read_t] != 0) {
					output_file4 << "k";
					continue;
				}
				if (world->worldECM[in].fragmnt_col[read_t] != 0) {
					output_file4 << "f";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == nothing) {
					output_file4 << "-";
				}
				if (world->worldPatch[in].type[read_t] == tissue) {
					output_file4 << "t";
				}
				if (world->worldPatch[in].type[read_t] == epithelium) {
					output_file4 << "e";
				}
				if (world->worldPatch[in].type[read_t] == capillary) {
					output_file4 << "o";
				}
				if (world->worldPatch[in].type[read_t] == unidentifiable) {
					output_file4 << "?";
				}
			}
			output_file4 << endl;
		}
		output_file4.close();

		// initElastin
		ofstream output_file5(initelastin, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (world->worldECM[in].monomer_ela[read_t] !=0) {
					output_file5 << "m";
					continue;
				}
				if (world->worldECM[in].fragmnt_ela[read_t] !=0) {
					output_file5 << "f";
					continue;
				}
				if (world->worldPatch[in].type[read_t] == nothing) {
					output_file5 << "-";
				}
				if (world->worldPatch[in].type[read_t] == tissue) {
					output_file5 << "t";
				}
				if (world->worldPatch[in].type[read_t] == epithelium) {
					output_file5 << "e";
				}
				if (world->worldPatch[in].type[read_t] == capillary) {
					output_file5 << "o";
				}
				if (world->worldPatch[in].type[read_t] == unidentifiable) {
					output_file5 << "?";
				}
			}
			output_file5 << endl;
		}
		output_file5.close();

		// Cells
		ofstream output_file6(cells, ios::app);

		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (world->worldPatch[in].isOccupied()) {
					if (world->worldPatch[in].occupiedby[read_t] == fibroblast) {
						output_file6 << "f";
						continue;
					}
					if (world->worldPatch[in].occupiedby[read_t] == platelet) {
						output_file6 << "p";
						continue;
					}
					if (world->worldPatch[in].occupiedby[read_t] == macrophag) {
						output_file6 << "m";
						continue;
					}
					if (world->worldPatch[in].occupiedby[read_t] == neutrophil) {
						output_file6 << "n";
						continue;
					}
				}
				if (world->worldPatch[in].type[read_t] == nothing) {
					output_file6 << "-";
				}
				if (world->worldPatch[in].type[read_t] == tissue) {
					output_file6 << "t";
				}
				if (world->worldPatch[in].type[read_t] == epithelium) {
					output_file6 << "e";
				}
				if (world->worldPatch[in].type[read_t] == capillary) {
					output_file6 << "o";
				}
				if (world->worldPatch[in].type[read_t] == unidentifiable) {
					output_file6 << "?";
				}
			}
			output_file6 << endl;
		}
		output_file6.close();
		ofstream output_file7(cells_w, ios::app);
		for (int iy = 0; iy < ny; iy++) {
			for (int ix = 0; ix < nx; ix++) {
				in = ix + iy*nx + iz*nx*ny;
				if (world->worldPatch[in].isOccupiedWrite()) {
					if (world->worldPatch[in].occupiedby[write_t] == fibroblast) {
						output_file7 << "f";
						continue;
					}
					if (world->worldPatch[in].occupiedby[write_t] == platelet) {
						output_file7 << "p";
						continue;
					}
					if (world->worldPatch[in].occupiedby[write_t] == macrophag) {
						output_file7 << "m";
						continue;
					}
					if (world->worldPatch[in].occupiedby[write_t] == neutrophil) {
						output_file7 << "n";
						continue;
					}
				}
				if (world->worldPatch[in].type[write_t] == nothing) {
					output_file7 << "-";
				}
				if (world->worldPatch[in].type[write_t] == tissue) {
					output_file7 << "t";
				}
				if (world->worldPatch[in].type[write_t] == epithelium) {
					output_file7 << "e";
				}
				if (world->worldPatch[in].type[write_t] == capillary) {
					output_file7 << "o";
				}
				if (world->worldPatch[in].type[write_t] == unidentifiable) {
					output_file7 << "?";
				}
			}
			output_file7 << endl;
		}
		output_file7.close();
		Number++;
	}
}

int outputColor(WHWorld* myWorld, char* fileName) {

	int dam = 0, tissue = 0, epi = 0, blood = 0, fib = 0, total = 0, plat = 0, black = 0, actDam = 0, newFib = 0, mac = 0, neu = 0;
	ofstream outfile(fileName);

  /* Prepare legacy VTK file format for for visualization with Paraview 3.0
   * (Kitware(Clifton Park, New York), Sandia National Labs(Livermore, CA),
   * CSimSoft(American Fork, Utah)). */
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Really cool data " << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_POINTS " <<endl;
	outfile << "DIMENSIONS " << myWorld->nx << " " << myWorld->ny << " " << myWorld->nz << endl;
	outfile << "ORIGIN 0 0 0 " << endl;
	outfile << "SPACING 1 1 1 " << endl;
	outfile << "POINT_DATA " << myWorld->nx*myWorld->ny*myWorld->nz << endl;
	outfile << "SCALARS Color int 1 " << endl;
	outfile << "LOOKUP_TABLE default " << endl;

	int in = 0;
	for (int iz = 0; iz < myWorld->nz; iz++) {

		for (int iy = 0; iy < myWorld->ny; iy++) {

			for (int ix = 0; ix < myWorld->nx; ix++) {

				if (ix == (myWorld->nx - 1) && iy == (myWorld->ny - 1)) {

					outfile << "139";  // Visualization color legend upper bound

				} else if (ix == (myWorld->nx - 2) && iy == (myWorld->ny - 1)) {

					outfile << "0 ";  // Visualization color legend lower bound

				} else {

					in = ix + iy*myWorld->nx + iz*myWorld->nx*myWorld->ny;

                                        // Output the color on each patch
					outfile << (myWorld->worldPatch[in].color[read_t]) << " ";

                                        // Count the number of cells and patch types
					if (myWorld->worldPatch[in].type[read_t] == tissue)
						tissue++;
					else if (myWorld->worldPatch[in].type[read_t] == epithelium)
						epi++;
					else if (myWorld->worldPatch[in].type[read_t] == capillary)
						blood++;
					else if (myWorld->worldPatch[in].type[read_t] == pdamage)
						black++;
					else if (myWorld->worldPatch[in].occupiedby[read_t] == afibroblast || myWorld->worldPatch[in].occupiedby[read_t] == fibroblast)
						fib++;
					else if (myWorld->worldPatch[in].occupiedby[read_t] == cplatelet)
						plat++;
					else if (myWorld->worldPatch[in].occupiedby[read_t] == amacrophag || myWorld->worldPatch[in].occupiedby[read_t] == macrophag)
						mac++;
					else if (myWorld->worldPatch[in].occupiedby[read_t] == aneutrophil || myWorld->worldPatch[in].occupiedby[read_t] == neutrophil)
						neu++;
					else {
						//cout << "the color is " << myWorld->worldPatch[in].color << " and location is ";
						//cout << ix << " " << iy << endl;
					}
					total++;
				}
			}
			outfile << endl;
		}
	}
  // Output the cell & patch type counts:
	//cout << "file name is " << fileName << endl;
	cout << " the counts are " << tissue << " " << epi << " " << blood << " " << fib << " " ;
        cout << newFib << " " << plat << " " << mac << " " << neu << " " << black << " " << total << endl;
	cout << "  to check, the num of fib and plat are " << Fibroblast::numOfFibroblasts << " ";
        cout << Platelet::numOfPlatelets << " " << Macrophage::numOfMacrophage << " " << Neutrophil::numOfNeutrophil << endl;
	return 0;
}

int outputECM(WHWorld* myWorld, char* fileName) {

  /* Prepare legacy VTK file format for for visualization with Paraview 3.0
   * (Kitware(Clifton Park, New York), Sandia National Labs(Livermore, CA),
   * CSimSoft(American Fork, Utah)). */
	ofstream outfile(fileName);
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Really cool data " << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_POINTS " << endl;
	outfile << "DIMENSIONS " << myWorld->nx << " " << myWorld->ny << " " << myWorld->nz << endl;
	outfile << "ORIGIN 0 0 0 "<< endl;
	outfile << "SPACING 1 1 1 "<< endl;
	outfile << "POINT_DATA " << myWorld->nx*myWorld->ny*myWorld->nz << endl;
	outfile << "SCALARS Color int 1 " << endl;
	outfile << "LOOKUP_TABLE default " << endl;

  // Assign the appropriate color to each patch:
  // Output the color on each patch:
        int in = 0; 
	for (int iz = 0; iz < myWorld->nz; iz++) {

		for (int iy = 0; iy < myWorld->ny; iy++) {

			for (int ix = 0; ix < myWorld->nx; ix++) {

                                in = ix + iy*myWorld->nx + iz*myWorld->nx*myWorld->ny;

                                if (myWorld->worldECM[in].isEmpty() == false) {

					if (myWorld->worldECM[in].monomer_ela[read_t] != 0 || myWorld->worldECM[in].polymer_ela != 0) {

						outfile << celastin << " ";

					} else if ( myWorld->worldECM[in].fragmnt_ela[read_t] != 0) {

						outfile << cfelastin << " ";

					} else if (myWorld->worldECM[in].monomer_col[read_t] != 0 || myWorld->worldECM[in].polymer_col != 0) {

						outfile << ccollagen << " ";

					} else if ( myWorld->worldECM[in].fragmnt_col[read_t] != 0) {

						outfile << cfcollagen << " ";

					//} else if (myWorld->worldECM[in].HA[read_t] != 0 ) {
					} else if (myWorld->worldECM[in].getnHA() != 0 ) {

						outfile << cHA << " ";

					//} else if (myWorld->worldECM[in].fHA[read_t]!= 0) {
					} else if (myWorld->worldECM[in].getnfHA() != 0) {
						outfile << cfHA << " ";
					}


				} else {

					outfile << myWorld->worldPatch[in].getColorfromType() << " ";

				}

			}
			outfile << endl;
		}
	}
	return 0;
}

void outputDiffusionKernel(float *kernel, int kw, int kh, int kd, char* fileName) {

	ofstream outfile(fileName);

  /* Prepare legacy VTK file format for for visualization with Paraview 3.0
   * (Kitware(Clifton Park, New York), Sandia National Labs(Livermore, CA),
   * CSimSoft(American Fork, Utah)). */
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Really cool data " << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_POINTS " << endl;
	outfile << "DIMENSIONS " << kw << " " << kh << " " << kd << endl;
	outfile << "ORIGIN 0 0 0 " << endl;
	outfile << "SPACING 1 1 1 " << endl;
	outfile << "POINT_DATA " << kw*kh*kd << endl;
	outfile << "SCALARS Color float 1 " << endl;
	outfile << "LOOKUP_TABLE default " << endl;

  // Output the chemical concentration on each patch:
	for (int iz = 0; iz < kd; iz++) {
		for (int iy = 0; iy < kh; iy++) {
			for (int ix = 0; ix < kw; ix++) {

				int in = ix + iy*kw + iz*kw*kh;

				outfile << (kernel[in]) << " ";

			}
			outfile << endl;
		}
	}
}

int outputChem(WHWorld* myWorld, char* fileName, int chemIndex) {
	ofstream outfile(fileName);
  /* Prepare legacy VTK file format for for visualization with Paraview 3.0
   * (Kitware(Clifton Park, New York), Sandia National Labs(Livermore, CA),
   * CSimSoft(American Fork, Utah)). */
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Really cool data " << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_POINTS " << endl;
	outfile << "DIMENSIONS " << myWorld->nx << " " << myWorld->ny << " " << myWorld->nz << endl;
	outfile << "ORIGIN 0 0 0 " << endl;
	outfile << "SPACING 1 1 1 " << endl;
	outfile << "POINT_DATA " << myWorld->nx*myWorld->ny*myWorld->nz << endl;
	outfile << "SCALARS Color float 1 " << endl;
	outfile << "LOOKUP_TABLE default " << endl;

  // Output the chemical concentration on each patch:
	for (int iz = 0; iz < myWorld->nz; iz++) {

		for (int iy = 0; iy < myWorld->ny; iy++) {

			for (int ix = 0; ix < myWorld->nx; ix++) {

				int in = ix + iy*myWorld->nx + iz*myWorld->nx*myWorld->ny;
#ifdef OPT_CHEM
// TODO: convert to dChem
					outfile << (myWorld->WHWorldChem->getPchem(chemIndex, in)) << " ";
#else		// OPT_CHEM
					outfile << (myWorld->WHWorldChem->pChem[chemIndex][in]) << " ";
#endif		// OPT_CHEM

			}
			outfile << endl;
		}
	}
	return 0;
}

}  // namespace util

#endif /* OUTPUT_UTILS_H_ */

