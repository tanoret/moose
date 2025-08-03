//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVDiffusion.h"

/**
 * Kernel that adds contributions from a diffusion term of the turbulent variables,
 * limited in the near-wall regions,
 * discretized using the finite volume method to a linear system.
 */
class LinearFVTurbulentMultiPhaseDiffusion : public LinearFVDiffusion
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVTurbulentMultiPhaseDiffusion(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;

  virtual Real computeNeighborMatrixContribution() override;

  virtual Real computeElemRightHandSideContribution() override;

  virtual Real computeNeighborRightHandSideContribution() override;

  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;

  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;

  virtual void addMatrixContribution() override;

  virtual void addRightHandSideContribution() override;

protected:

  /// Turbulent Kinetic Energy
  const Moose::Functor<Real> * _k;

  /// Density
  const Moose::Functor<Real> * _rho;

  /// Dynamic viscosity
  const Moose::Functor<Real> * _mu;

  /// Wall distance
  const Moose::Functor<Real> * _d;

  /// The phase fraction
  const Moose::Functor<Real> & _alpha;

  /// The functor for the scaling coefficient for the diffusion term
  const Moose::Functor<Real> & _scaling_coeff;

  /// Two-layer closure parameters
  const Real _Re_y_star;
};
