//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "INSFVVelocityVariable.h"

/**
 * Computes the turbuent viscosity for the k-Epsilon model.
 * Implements two near-wall treatments: equilibrium and non-equilibrium wall functions.
 */
class kOmegaSSTViscosityAux : public AuxKernel
{
public:
  static InputParameters validParams();

  virtual void initialSetup() override;

  kOmegaSSTViscosityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// The dimension of the domain
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<ADReal> & _u_var;
  /// y-velocity
  const Moose::Functor<ADReal> * _v_var;
  /// z-velocity
  const Moose::Functor<ADReal> * _w_var;

  /// Turbulent kinetic energy
  const Moose::Functor<ADReal> & _k;
  /// Turbulent kinetic energy specific dissipation rate
  const Moose::Functor<ADReal> & _omega;

  /// Density
  const Moose::Functor<ADReal> & _rho;
  /// Dynamic viscosity
  const Moose::Functor<ADReal> & _mu;

  /// F2 blending function
  const Moose::Functor<ADReal> & _F2;

  /// Wall boundaries
  const std::vector<BoundaryName> & _wall_boundary_names;

  /// If the user wants to enable bulk wall treatment
  const bool _bulk_wall_treatment;

  /// Method used for wall treatment
  const MooseEnum _wall_treatment;

  /// Method used to limit the k-e time scale
  const MooseEnum _scale_limiter;

  ///@{
  /// Maps for wall bounded elements
  std::map<const Elem *, bool> _wall_bounded;
  std::map<const Elem *, std::vector<Real>> _dist;
  std::map<const Elem *, std::vector<const FaceInfo *>> _face_infos;
  ///@}

  /// Model constants
  static constexpr Real _C_mu = 0.09;
  static constexpr Real _R_k = 6.0;
  static constexpr Real _alpha_0_star = 0.072 / 3.0;
  static constexpr Real _alpha_infty_star = 1.0;
  static constexpr Real _a_1 = 0.31;
};
