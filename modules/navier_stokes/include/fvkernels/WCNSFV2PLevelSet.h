//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFluxKernel.h"
#include "INSFVVelocityVariable.h"

/**
 * Computes the turbulent diffusion of energy term in the weakly compressible formulation
 * of the energy equation, using functor material properties
 */
class WCNSFV2PLevelSet : public FVFluxKernel
{
public:
  static InputParameters validParams();

  WCNSFV2PLevelSet(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  /// the reinitialization parameter - gamma
  const Moose::Functor<ADReal> & _reinitialization_parameter;

  /// the reinitialization parameter - gamma
  const Moose::Functor<ADReal> & _interface_thickness;

  /// Regularizer control
  const Moose::Functor<ADReal> & _regularizer_control;
};
