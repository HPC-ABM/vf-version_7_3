/*
 * parameters.h
 *
 * File Contents: TODO
 *
 * Created on: July 20, 2015
 * Author: Kimberley Trickey
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

using namespace std;
namespace util {

/*
 * Description:	Initializes parameters from the given input file.
 *
 * Return: void
 *
 * Parameters: filename  -- Name of the file containing the parameter values
 */
void processParameters(string filename) {
#ifdef CALIBRATION
  cout << "Processing Parameters..." << endl;
//  ifstream input_file("Sample.txt", ios::in);
  ifstream input_file;
  input_file.open(filename.c_str(), ios::in);
//  ifstream input_file(filename, ios::in);
//  ifstream input_file("NoExists.txt", ios::in);
  if (!input_file) {
    cerr << "Could not open input file for processing parameters" << endl;
    return;
  }

  string line;
  getline(input_file, line);
  stringstream lineStream(line);
  string value;

  // Set platelet cytokine synthesis parameters
  for (int i = 0; i < 3 ; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Platelet::cytokineSynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Platelet::cytokineSynthesis[" << i << "] = " << setprecision(14) << Platelet::cytokineSynthesis[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to platelet cytokine synthesis parameter #"<< i << endl;
    }
  }

  // Set neutrophil cytokine synthesis parameters
  for (int i = 0; i < 21; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Neutrophil::cytokineSynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Neutrophil::cytokineSynthesis[" << i << "] = " << Neutrophil::cytokineSynthesis[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to neutrophil cytokine synthesis parameter #"<< i << endl;
    }
  }

  // Set macrophage cytokine synthesis parameters
  for (int i = 0; i < 82; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Macrophage::cytokineSynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Macrophage::cytokineSynthesis[" << i << "] = " << Macrophage::cytokineSynthesis[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to macrophage cytokine synthesis parameter #"<< i << endl;
    }
  }

  // Set fibroblast cytokine synthesis parameters
  for (int i = 0; i < 32; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Fibroblast::cytokineSynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Fibroblast::cytokineSynthesis[" << i << "] = " << Fibroblast::cytokineSynthesis[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to fibroblast cytokine synthesis parameter #"<< i << endl;
    }
  }

  // Set resonant voice impact stress parameter
  if (getline(lineStream, value, '\t')) {
    float value_as_float = atof(value.c_str());
    WHWorld::RVIS = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
    cout << "WHWorld::RVIS" << " = " << WHWorld::RVIS << endl;
#endif
  } else {
    cerr << "Error in assigning value to resonant voice impact stress parameter" << endl;
  }

  // Set resonant voice vibratory stress parameter
  if (getline(lineStream, value, '\t')) {
    float value_as_float = atof(value.c_str());
    WHWorld::RVVS = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
    cout << "WHWorld::RVVS" << " = " << WHWorld::RVVS << endl;
#endif
  } else {
    cerr << "Error in assigning value to resonant voice vibratory stress parameter" << endl;
  }

  // Set spontaneous speech impact stress parameter
  if (getline(lineStream, value, '\t')) {
    float value_as_float = atof(value.c_str());
    WHWorld::SSIS = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
    cout << "WHWorld::SSIS" << " = " << WHWorld::SSIS << endl;
#endif
  } else {
    cerr << "Error in assigning value to spontaneous speech impact stress parameter" << endl;
  }

  // Set spontaneous speech vibratory stress parameter
  if (getline(lineStream, value, '\t')) {
    float value_as_float = atof(value.c_str());
    WHWorld::SSVS = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
    cout << "WHWorld::SSVS" << " = " << WHWorld::SSVS << endl;
#endif
  } else {
    cerr << "Error in assigning value to spontaneous speech vibratory stress parameter" << endl;
  }

  // Set TNF damage threshold parameter
  if (getline(lineStream, value, '\t')) {
    float value_as_float = atof(value.c_str());
    WHWorld::thresholdTNFdamage = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
    cout << "WHWorld::thresholdTNFdamage" << " = " << WHWorld::thresholdTNFdamage << endl;
#endif
  } else {
    cerr << "Error in assigning value to TNF damage threshold parameter" << endl;
  }

  // Set MMP8 damage threshold parameter
  if (getline(lineStream, value, '\t')) {
    float value_as_float = atof(value.c_str());
    WHWorld::thresholdMMP8damage = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
    cout << "WHWorld::thresholdMMP8damage" << " = " << WHWorld::thresholdMMP8damage << endl;
#endif
  } else {
    cerr << "Error in assigning value to MMP8 damage threshold parameter" << endl;
  }

  // Set WHWorld sprouting frequency parameters
  for (int i = 0; i < 6; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      WHWorld::sproutingFrequency[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "WHWorld::sproutingFrequency[" << i << "] = " << WHWorld::sproutingFrequency[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to WHWorld sprouting frequency parameter #"<< i << endl;
    }
  }

  // Set WHWorld sprouting amount parameters
  for (int i = 0; i < 14; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      WHWorld::sproutingAmount[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "WHWorld::sproutingAmount[" << i << "] = " << WHWorld::sproutingAmount[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to WHWorld sprouting amount parameter #"<< i << endl;
    }
  }

  // Set fibroblast activation parameters
  for (int i = 0; i < 5; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Fibroblast::activation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Fibroblast::activation[" << i << "] = " << Fibroblast::activation[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to fibroblast activation parameter #"<< i << endl;
    }
  }

  // Set macrophage activation parameters
  for (int i = 0; i < 5; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Macrophage::activation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Macrophage::activation[" << i << "] = " << Macrophage::activation[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to macrophage activation parameter #"<< i << endl;
    }
  }

  // Set neutrophil activation parameters
  for (int i = 0; i < 4; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Neutrophil::activation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Neutrophil::activation[" << i << "] = " << Neutrophil::activation[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to neutrophil activation parameter #"<< i << endl;
    }
  }

  // Set neutrophil death parameters
  for (int i = 0; i < 2; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Neutrophil::death[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Neutrophil::death[" << i << "] = " << Neutrophil::death[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to neutrophil death parameter #"<< i << endl;
    }
  }

  // Set fibroblast ECM synthesis parameters
  for (int i = 0; i < 19; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Fibroblast::ECMsynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Fibroblast::ECMsynthesis[" << i << "] = " << Fibroblast::ECMsynthesis[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to fibroblast ECM synthesis parameter #"<< i << endl;
    }
  }

  // Set fibroblast proliferation parameters
  for (int i = 0; i < 6; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Fibroblast::proliferation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
      cout << "Fibroblast::proliferation[" << i << "] = " << Fibroblast::proliferation[i] << endl;
#endif
    } else {
      cerr << "Error in assigning value to fibroblast proliferation parameter #"<< i << endl;
    }
  }

  // Set WHWorld cytokine decay parameters
  for (int i = 0; i < 8; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      // Convert half-life to decay rate
      if (value_as_float > 0) {
        float decay_rate = pow(0.5, 30/value_as_float);
        WHWorld::cytokineDecay[i] = decay_rate;
	WHWorld::halfLifes_static[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
        cout << "WHWorld::halfLifes_static[" << i << "] = " << value_as_float << 
		", WHWorld::cytokineDecay["<< i << "] = " << WHWorld::cytokineDecay[i] << endl;
#endif
      } else {
        cerr << "Error in assigning value to WHWorld half life parameter and cytokine decay parameter #" << i << " (Half-life <= 0)" << endl;
      }
    } else {
      cerr << "Error in assigning value to WHWorld half life parameter and cytokine decay parameter #" << i << endl;
    }
  }

  input_file.close();
#endif  // ifdef CALIBRATION

  return;
}

void outputTotalChem(WHWorld* myWorld, string filename) {
#ifdef CALIBRATION
  cout << "Outputting total chem to " << filename << endl;
	//ofstream output_file(filename, ios::app);
  ofstream output_file;
  output_file.open(filename.c_str(), ios::app);

  output_file << fixed << myWorld->WHWorldChem->total[TNF] << "\t";
  output_file << myWorld->WHWorldChem->total[TGF] << "\t";
  output_file << myWorld->WHWorldChem->total[FGF] << "\t";
  output_file << myWorld->WHWorldChem->total[IL1beta] << "\t";
  output_file << myWorld->WHWorldChem->total[IL6] << "\t";
  output_file << myWorld->WHWorldChem->total[IL8] << "\t";
  output_file << myWorld->WHWorldChem->total[IL10] << "\t";
  output_file << myWorld->WHWorldChem->total[MMP8] << "\t";
  output_file << endl;
#endif  // ifdef CALIBRATION
	return;
}

}  // namespace util

#endif /* PARAMETERS_H_ */
