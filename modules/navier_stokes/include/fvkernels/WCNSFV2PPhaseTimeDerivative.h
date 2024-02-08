//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVElementalKernel.h"

/**
 * Computes the phase time derivative for the two=phase weakly compressible formulation of the
 * momentum equation, using functor material properties.
 */
class WCNSFV2PPhaseTimeDerivative : public FVElementalKernel
{
public:
  static InputParameters validParams();
  WCNSFV2PPhaseTimeDerivative(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  /// The density
  const Moose::Functor<ADReal> & _rho;
  /// The time derivative of density
  const Moose::Functor<ADReal> & _rho_dot;
};
