//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

class INSFVVelocityVariable;

/**
 * Computes the interphase drag coefficient for the two-phase mixture model
 */
class WCNSFV2PInterphaseDragFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  WCNSFV2PInterphaseDragFunctorMaterial(const InputParameters & parameters);

protected:
  /// the dimension of the simulation
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<Real> * const _u_var;
  /// y-velocity
  const Moose::Functor<Real> * const _v_var;
  /// z-velocity
  const Moose::Functor<Real> * const _w_var;

  /// x-velocity coupled phase
  const Moose::Functor<Real> * const _u_var_coupled;
  /// y-velocity coupled phase
  const Moose::Functor<Real> * const _v_var_coupled;
  /// z-velocity coupled phase
  const Moose::Functor<Real> * const _w_var_coupled;

  /// Main Phase Void Fraction
  const Moose::Functor<Real> * const _alpha_main;

  /// Coupled Density Void Fraction
  const Moose::Functor<Real> & _alpha;

  /// Density
  const Moose::Functor<Real> & _rho;

  /// Density coupled phase
  const Moose::Functor<Real> & _rho_coupled;

  /// Viscosity
  const Moose::Functor<Real> & _mu;

  /// Viscosity Coupled Phase
  const Moose::Functor<Real> & _mu_coupled;

  /// Particle diameter in the dispersed phase
  const Moose::Functor<Real> & _particle_diameter;

  /// Drag formulation
  const MooseEnum & _drag_formulation_type;

  /// Drag model
  const MooseEnum & _drag_model;
};