//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TransientBase.h"
#include "PIMPLEMultiPhaseSolve.h"

/**
 * Executioner set up to solve a transient thermal-hydraulics problem using the PIMPLEMultiPhase algorithm.
 * It utilizes segregated linear systems which are solved using a fixed-point iteration within time
 * steps.
 */
class PIMPLEMultiPhase : public TransientBase
{
public:
  static InputParameters validParams();

  PIMPLEMultiPhase(const InputParameters & parameters);

  virtual void init() override;

  virtual Real relativeSolutionDifferenceNorm(bool check_aux) const override;

protected:
  virtual std::set<TimeIntegrator *> getTimeIntegrators() const override;

  /// The solve object performing the PIMPLEMultiPhase iteration
  PIMPLEMultiPhaseSolve _PIMPLEMultiPhase_solve;
};
