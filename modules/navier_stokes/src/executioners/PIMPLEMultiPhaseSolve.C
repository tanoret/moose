//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PIMPLEMultiPhaseSolve.h"
#include "FEProblem.h"
#include "SegregatedSolverUtils.h"
#include "LinearSystem.h"

using namespace libMesh;

InputParameters
PIMPLEMultiPhaseSolve::validParams()
{
  InputParameters params = SolveBaseMultiPhase::validParams();
  params.addParam<unsigned int>(
      "num_piso_iterations",
      0,
      "The number of PISO iterations without recomputing the momentum matrix.");

  return params;
}

PIMPLEMultiPhaseSolve::PIMPLEMultiPhaseSolve(Executioner & ex)
  : SolveBaseMultiPhase(ex),
    _num_piso_iterations(getParam<unsigned int>("num_piso_iterations"))
{
}

std::pair<unsigned int, Real>
PIMPLEMultiPhaseSolve::correctVelocities(const bool /*subtract_updated_pressure*/,
                                 const bool /*recompute_face_mass_flux*/,
                                 const SolverParams & solver_params)
{
  std::pair<unsigned int, Real> residual;
  unsigned int piso_iteration_counter = 0;
  while (piso_iteration_counter <= _num_piso_iterations)
  {
    residual = SolveBaseMultiPhase::correctVelocities(
        piso_iteration_counter == 0, piso_iteration_counter == _num_piso_iterations, solver_params);
    piso_iteration_counter++;
  }

  return residual;
}
