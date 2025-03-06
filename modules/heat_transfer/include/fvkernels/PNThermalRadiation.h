//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVDiffusion.h"

class PNThermalRadiation : public FVDiffusion
{
public:

  static InputParameters validParams();
  PNThermalRadiation(const InputParameters & params);

protected:

  // Overwriting current residual
  virtual ADReal computeQpResidual() override final;

  // Order of the PN equation
  const unsigned int _n;

  /// Lagged flux order
  const Moose::Functor<ADReal> & _phi_n_minus_2;

  /// Leading flux order
  const Moose::Functor<ADReal> & _phi_n_plus_2;

  /// Lagged scattering cross section
  const Moose::Functor<ADReal> & _sigma_n_minus_1;

  /// Leading scattering cross section
  const Moose::Functor<ADReal> & _sigma_n_plus_1;
};
