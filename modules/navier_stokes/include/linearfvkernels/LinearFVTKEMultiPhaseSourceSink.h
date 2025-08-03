//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"
#include "NS.h"

/**
 * Kernel that adds contributions to the source and the sink of the turbulent kinetic energy
 * discretized using the finite volume method to a linear system.
 */
class LinearFVTKEMultiPhaseSourceSink : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVTKEMultiPhaseSourceSink(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:
  /// The dimension of the domain
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<Real> & _u_var;
  /// y-velocity
  const Moose::Functor<Real> * _v_var;
  /// z-velocity
  const Moose::Functor<Real> * _w_var;

  /// epsilon - dissipation rate of TKE
  const Moose::Functor<Real> & _epsilon;

  /// Density
  const Moose::Functor<Real> & _rho;

  /// Dynamic viscosity
  const Moose::Functor<Real> & _mu;

  /// Wall distance
  const Moose::Functor<Real> & _d;

  /// The phase fraction
  const Moose::Functor<Real> & _alpha;

  /// C_mu constant
  const Real _C_mu;

  /// Production Limiter Constant
  const Real _C_pl;

  /// Two-layer closure parameters
  const Real _Re_y_star;
  const Real _delta_Re_y;
};
