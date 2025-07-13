//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVMomentumPressure.h"

/**
 * Kernel that adds the component of the pressure gradient in the momentum
 * equations to the right hand side.
 */
class LinearFVMultiPhaseMomentumPressure : public LinearFVMomentumPressure
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVMultiPhaseMomentumPressure(const InputParameters & params);

  virtual Real computeRightHandSideContribution() override;

protected:

  const Moose::Functor<Real> & _alpha;
};
