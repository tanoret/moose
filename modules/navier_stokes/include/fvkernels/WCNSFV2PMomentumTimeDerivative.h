//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "WCNSFVMomentumTimeDerivative.h"

/**
 * Computes the momentum time derivative for the two-phase weakly compressible formulation of the
 * momentum equation, using functor material properties. Only one spatial component is included.
 */
class WCNSFV2PMomentumTimeDerivative : public WCNSFVMomentumTimeDerivative
{
public:
  static InputParameters validParams();
  WCNSFV2PMomentumTimeDerivative(const InputParameters & params);

  using WCNSFVMomentumTimeDerivative::gatherRCData;
  void gatherRCData(const Elem &) override;

protected:
  /// The phase fraction
  const Moose::Functor<ADReal> & _alpha;
};
