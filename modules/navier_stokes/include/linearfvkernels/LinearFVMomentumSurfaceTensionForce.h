//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVElementalKernel.h"
#include "FaceArgInterface.h"

/**
 * Kernel that adds the component of the pressure gradient in the momentum
 * equations to the right hand side.
 */
class LinearFVMomentumSurfaceTensionForce : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVMomentumSurfaceTensionForce(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:

  MooseLinearVariableFV<Real> & getAlphaVariable(const std::string & vname);

  /// The dimension of the problem
  const unsigned int _dim;

  /// Index x|y|z of the momentum equation component
  const unsigned int _index;

  /// The surface tension value
  const Moose::Functor<Real> & _sigma;

  /// The phase fraction
  MooseLinearVariableFV<Real> & _alpha;

  /// Reconstruction approach
};
