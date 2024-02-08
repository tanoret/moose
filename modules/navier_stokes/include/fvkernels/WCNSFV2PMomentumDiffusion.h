//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MathFVUtils.h"
#include "INSFVMomentumDiffusion.h"
#include "INSFVMomentumResidualObject.h"
#include "INSFVVelocityVariable.h"

class WCNSFV2PMomentumDiffusion : public INSFVMomentumDiffusion
{
public:
  static InputParameters validParams();
  WCNSFV2PMomentumDiffusion(const InputParameters & params);

protected:
  virtual ADReal computeStrongResidual(const bool populate_a_coeffs) override;

  // Phase fraction
  const Moose::Functor<ADReal> & _alpha;
};
