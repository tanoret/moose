//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVMomentumAdvection.h"
#include "PiecewiseByBlockLambdaFunctor.h"

/**
 * An advection kernel that implements interpolation schemes specific to Navier-Stokes flow
 * physics
 */
class WCNSFV2PMomentumAdvection : public INSFVMomentumAdvection
{
public:
  static InputParameters validParams();
  WCNSFV2PMomentumAdvection(const InputParameters & params);
  void initialSetup() override;

protected:
  void computeResidualsAndAData(const FaceInfo & fi) override;

  // Phase fraction
  const Moose::Functor<ADReal> & _alpha;
};
