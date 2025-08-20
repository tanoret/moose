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

  /// Temperature at the boundary
  const Moose::Functor<ADReal> & _Tb;

  /// Frquency
  const Moose::Functor<ADReal> & _nu;

  /// Frequency bounds for numerical integration
  const Moose::Functor<ADReal> * _nu_low;
  const Moose::Functor<ADReal> * _nu_high;

  /// Refraction index
  const Moose::Functor<ADReal> & _n1;

  /// kappa Thickness
  const Moose::Functor<ADReal> & _absorptivity;

  /// Optical Thickness
  const Moose::Functor<ADReal> & _optical_thickness;

  /// Radiation from the other moment
  const Moose::Functor<ADReal> & _psi;

  // /// Order
  const MooseEnum & _order;

  /// Coefficients
  const Moose::Functor<ADReal> & _alpha;
  const Moose::Functor<ADReal> & _beta;
  const Moose::Functor<ADReal> & _eta;

  /// Governing parameters for SP3 moments
  // Real _alpha_order;
  // Real _beta_order;
  // Real _eta_order;
  Real _squared_mu_order;

  /// Closure parameters
  const Real _squared_mu_1 = 3./7. - 2./7.*std::sqrt(6./5.);
  const Real _squared_mu_2 = 3./7. + 2./7.*std::sqrt(6./5.);


  /// Coefficients
  // const Real _alpha_1 = 0.2115609606560856;
  // const Real _alpha_2 = 0.5970451183268032;
  // const Real _beta_1 = 0.19452369229011088; 
  // const Real _beta_2 = -0.4203259177124563;
  // const Real _eta_1 = -2.623417821661131; 
  // const Real _eta_2 = 9.947147040979628; 

  // const Real _alpha_1 = 5./96. * (34. + 11. * std::sqrt(6./5.));
  // const Real _alpha_2 = 5./96. * (34. - 11. * std::sqrt(6./5.));
  // const Real _beta_1 = 5./96. * (2. - std::sqrt(6./5.));
  // const Real _beta_2 = 5./96. * (2. + std::sqrt(6./5.));
  // const Real _eta_1 = 5.*libMesh::pi/2. * (3. + std::sqrt(6./5.));
  // const Real _eta_2 = 5.*libMesh::pi/2. * (3. - std::sqrt(6./5.));
};