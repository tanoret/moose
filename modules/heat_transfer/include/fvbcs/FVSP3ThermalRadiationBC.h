//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFluxBC.h"

/**
 * Robin boundary condition (temperatures) for finite volume scheme between
 * a solid and fluid where the temperatures and heat transfer coefficient
 * are given as a functors
 */
class FVSP3ThermalRadiationBC : public FVFluxBC
{
public:
  FVSP3ThermalRadiationBC(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual ADReal computeQpResidual() override;

  /// Temperature
  const Moose::Functor<ADReal> & _T;

  /// Frquency
  const Moose::Functor<ADReal> & _nu;

  /// Refraction index
  const Moose::Functor<ADReal> & _n1;

  /// kappa Thickness
  const Moose::Functor<ADReal> & _absorptivity;

  /// Optical Thickness
  const Moose::Functor<ADReal> & _optical_thickness;

  /// Radiation from the other moment
  const Moose::Functor<ADReal> & _psi;

  /// Order
  const MooseEnum & _order;

  /// Governing parameters for SP3 moments
  Real _alpha_order;
  Real _beta_order;
  Real _eta_order;

  /// Closure parameters
  const Real _alpha_1 = 5./96. * (34. + 11. * std::sqrt(6./5.));
  const Real _alpha_2 = 5./96. * (34. - 11. * std::sqrt(6./5.));
  const Real _beta_1 = 5./96. * (2. - std::sqrt(6./5.));
  const Real _beta_2 = 5./96. * (2. + std::sqrt(6./5.));
  const Real _eta_1 = 5.*libMesh::pi/2. * (3. + std::sqrt(6./5.));
  const Real _eta_2 = 5.*libMesh::pi/2. * (3. - std::sqrt(6./5.));
};
