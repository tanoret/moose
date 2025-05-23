//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVElementalKernel.h"
#include "INSFVMomentumResidualObject.h"

/**
 * Computes the interphase drag force.
 */
class WCNSFV2PMomentumInterphaseForce : public FVElementalKernel, public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  WCNSFV2PMomentumInterphaseForce(const InputParameters & params);

  // Do not contribute to RC
  void gatherRCData(const Elem &) override {}
  void gatherRCData(const FaceInfo &) override {}

protected:
  ADReal computeQpResidual() override;

  /// the dimension of the simulation
  const unsigned int _dim;

  /// velocity coupled phase
  const Moose::Functor<ADReal> * _u_var_coupled;
  const Moose::Functor<ADReal> * _v_var_coupled;
  const Moose::Functor<ADReal> * _w_var_coupled;

  /// velocity phase
  const Moose::Functor<ADReal> * _u_var;
  const Moose::Functor<ADReal> * _v_var;
  const Moose::Functor<ADReal> * _w_var;

  /// drag coefficient
  const Moose::Functor<ADReal> & _drag_coef;

  /// Coupled Density Void Fraction
  const Moose::Functor<ADReal> & _alpha;

  /// Density
  const Moose::Functor<ADReal> & _rho;

  /// Booleans to activate body force
  const bool & _bool_activate_drag;
  const bool & _bool_activate_lift;
  const bool & _bool_activate_virtual_mass;
};
