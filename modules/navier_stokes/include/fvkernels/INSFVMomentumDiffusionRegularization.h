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

/// INSFVMomentumDiffusionRegularization implements a standard diffusion term:
///
///     - strong form: \nabla \cdot k \nabla u
///
///     - weak form: \int_{A} k \nabla u \cdot \vec{n} dA
///
/// It uses/requests a material property named "alpha" for k. An average of
/// the elem and neighbor k-values (which should be face-values) is used to
/// compute k on the face. Cross-diffusion correction factors are currently not
/// implemented for the "grad_u*n" term.
class INSFVMomentumDiffusionRegularization : public FVFluxKernel
{
public:
  static InputParameters validParams();
  INSFVMomentumDiffusionRegularization(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;
  const Moose::Functor<ADReal> & _alpha;
};
