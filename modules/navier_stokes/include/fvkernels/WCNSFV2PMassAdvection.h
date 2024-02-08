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
#include "INSFVBCInterface.h"

/**
 * A flux kernel transporting mass across cell faces in weakly compressible simulations
 */
class WCNSFV2PMassAdvection : public FVFluxKernel, public INSFVBCInterface
{
public:
  static InputParameters validParams();
  WCNSFV2PMassAdvection(const InputParameters & params);
  void initialSetup() override;

protected:
  virtual ADReal computeQpResidual() override;
  // virtual bool hasMaterialTimeDerivative() const override { return false; }
  virtual bool hasMaterialTimeDerivative() { return false; };

  bool skipForBoundary(const FaceInfo & fi) const override;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;

  /// The interpolation method to use for the velocity
  Moose::FV::InterpMethod _velocity_interp_method;

  // The densities
  std::vector<MooseFunctorName> _rho_names;
  std::vector<const Moose::Functor<ADReal> *> _rho;

  // Phase fractions
  std::vector<MooseFunctorName> _alpha_names;
  std::vector<const Moose::Functor<ADReal> *> _alpha;

  // Rhie Chow Interpolators
  std::vector<UserObjectName> _rc_names;
  std::vector<const RhieChowInterpolatorBase *> _rc_velocity_providers;
};
