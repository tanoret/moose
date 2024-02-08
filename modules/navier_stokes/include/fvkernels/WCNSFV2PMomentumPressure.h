//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVMomentumPressure.h"

/**
 * Adds the coupled pressure term into the two-phase Navier Stokes equations
 */
class WCNSFV2PMomentumPressure : public INSFVMomentumPressure
{
public:
  static InputParameters validParams();
  WCNSFV2PMomentumPressure(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  // Phase fraction
  const Moose::Functor<ADReal> & _alpha;
};
