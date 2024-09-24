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
/**
 * Computes source the sink terms for the interface area in the mixture model of two-phase flows.
 */
class FVElectrophoresisSource : public FVFluxKernel
{
public:
  static InputParameters validParams();

  FVElectrophoresisSource(const InputParameters & parameters);

  void computeResidual(const FaceInfo & fi) override;
  void computeJacobian(const FaceInfo & fi) override;

protected:
  ADReal computeQpResidual() override;

  /// The dimension of the domain
  const unsigned int _dim;

  /// Model constants
  const Real _F;
  const Real _z;
  const Real _R;

  /// Temperature field
  const Moose::Functor<ADReal> & _T;

  /// Diffusion coefficient
  const Moose::Functor<ADReal> & _D;

  /// Electric potential
  const Moose::Functor<ADReal> & _phi;
};
