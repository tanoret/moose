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
 * Computes the mass time derivative for the weakly compressible formulation of the mass
 * equation, using functor material properties
 */
class WCNSFV2PMassTimeDerivative : public FVElementalKernel
{
public:
  static InputParameters validParams();
  WCNSFV2PMassTimeDerivative(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  // The densities
  std::vector<MooseFunctorName> _rho_names;
  std::vector<const Moose::Functor<ADReal> *> _rho;

  // Phase fractions
  std::vector<MooseFunctorName> _alpha_names;
  std::vector<const Moose::Functor<ADReal> *> _alpha;
};
