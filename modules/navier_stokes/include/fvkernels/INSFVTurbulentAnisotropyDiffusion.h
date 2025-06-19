//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVDiffusion.h"

/// INSFVTurbulentAnisotropyDiffusion implements a standard diffusion term for a turbulent problem:
///
///     - strong form: \nabla \cdot k \nabla u / coef
///
///     - weak form: \int_{A} k \nabla u / coef \cdot \vec{n} dA
///
/// It uses/requests a material property named "coeff" for k. An average of
/// the elem and neighbor k-values (which should be face-values) is used to
/// compute k on the face. Cross-diffusion correction factors are currently not
/// implemented for the "grad_u*n" term.
/// The specialty of this kernel is that it takes into account the wall treatment of the variable with respect to turbulence.
class INSFVTurbulentAnisotropyDiffusion : public FVDiffusion
{
public:
  static InputParameters validParams();
  virtual void initialSetup() override;
  INSFVTurbulentAnisotropyDiffusion(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override final;
  using FVDiffusion::computeResidual;
  void computeResidual(const FaceInfo & fi) override;
  using FVDiffusion::computeJacobian;
  void computeJacobian(const FaceInfo & fi) override;

  const unsigned int _mesh_dim;

  /// Wall boundaries
  const std::vector<BoundaryName> & _wall_boundary_names;

  /// Maps for wall treatment
  std::map<const Elem *, bool> _wall_bounded;

  /// x-velocity
  const Moose::Functor<Real> & _bc0;
  /// y-velocity
  const Moose::Functor<Real> * _bc1;
  /// z-velocity
  const Moose::Functor<Real> * _bc2;
  /// k for viscosity
  const Moose::Functor<ADReal> & _k;
  /// rho for density
  const Moose::Functor<ADReal> & _rho; 
};
