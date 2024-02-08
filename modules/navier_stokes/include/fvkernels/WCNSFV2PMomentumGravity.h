//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVMomentumGravity.h"

/**
 * Imposes a gravitational force in two-phase Navier Stokes
 */
class WCNSFV2PMomentumGravity : public INSFVMomentumGravity
{
public:
  static InputParameters validParams();
  WCNSFV2PMomentumGravity(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  // Phase fraction
  const Moose::Functor<ADReal> & _alpha;

  // Mixture Density
  const Moose::Functor<ADReal> & _rho_mixture;
};
