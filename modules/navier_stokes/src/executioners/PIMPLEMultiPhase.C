//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "PIMPLEMultiPhase.h"
#include "FEProblem.h"
#include "AuxiliarySystem.h"
#include "LinearSystem.h"

using namespace libMesh;

registerMooseObject("NavierStokesApp", PIMPLEMultiPhase);

InputParameters
PIMPLEMultiPhase::validParams()
{
  InputParameters params = TransientBase::validParams();
  params.addClassDescription(
      "Solves the transient Navier-Stokes equations using the PIMPLEMultiPhase algorithm and "
      "linear finite volume variables.");
  params += PIMPLEMultiPhaseSolve::validParams();

  return params;
}

PIMPLEMultiPhase::PIMPLEMultiPhase(const InputParameters & parameters) : TransientBase(parameters), _PIMPLEMultiPhase_solve(*this)
{
  _fixed_point_solve->setInnerSolve(_PIMPLEMultiPhase_solve);
}

void
PIMPLEMultiPhase::init()
{
  TransientBase::init();
  _PIMPLEMultiPhase_solve.linkRhieChowUserObjects();
  _PIMPLEMultiPhase_solve.setupPressurePin();
}

Real
PIMPLEMultiPhase::relativeSolutionDifferenceNorm(bool check_aux) const
{
  if (check_aux)
    return _aux.solution().l2_norm_diff(_aux.solutionOld()) / _aux.solution().l2_norm();
  else
  {
    // Default criterion for now until we add a "steady-state-convergence-object" option
    Real residual = 0;
    for (const auto sys : _PIMPLEMultiPhase_solve.systemsToSolve())
      residual +=
          std::pow(sys->solution().l2_norm_diff(sys->solutionOld()) / sys->solution().l2_norm(), 2);
    return std::sqrt(residual);
  }
}

std::set<TimeIntegrator *>
PIMPLEMultiPhase::getTimeIntegrators() const
{
  // We use a set because time integrators were added to every system, and we want a unique
  std::set<TimeIntegrator *> tis;
  // Get all time integrators from the systems in the FEProblemSolve
  for (const auto sys : _PIMPLEMultiPhase_solve.systemsToSolve())
    for (const auto & ti : sys->getTimeIntegrators())
      tis.insert(ti.get());
  return tis;
}
