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
 * Simple class to demonstrate off diagonal Jacobian contributions.
 */
class INSFVMomentumSourceRegularization : public FVElementalKernel, public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();

  INSFVMomentumSourceRegularization(const InputParameters & parameters);

  // This object neither contributes to the A coefficients nor to the B (source) coefficients
  void gatherRCData(const Elem &) override {}
  void gatherRCData(const FaceInfo &) override {}

protected:

  ADReal computeQpResidual() override;

  /// Variable to couple with
  const ADVariableValue & _u_reg;

  /// Coupled scaling coefficient
  const Real _scaling_coef;
};
