/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"

// Define the non-singluar method used for calculations
#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method

#include <DefectiveCrystal.h>

using namespace model;

int main (int argc, char* argv[])
{

#ifdef _MODEL_PYBIND11_ // COMPILED WITH PYBIND11
    pybind11::scoped_interpreter guard{};
#endif

    const std::string folderName(argc>1? std::string(argv[1]) : "./");
    
    
    DislocationDynamicsBase<3> ddBase(folderName);
    DefectiveCrystal<3,0> DC(ddBase);
    DC.runGlideSteps();
    
    return 0;
}
#pragma GCC diagnostic pop
