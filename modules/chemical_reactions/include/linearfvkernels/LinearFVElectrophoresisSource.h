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

/**
 * Kernel that adds the component of the electrophoresis term
 */
class LinearFVElectrophoresisSource : public LinearFVElementalKernel
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param params The InputParameters for the kernel.
   */
  LinearFVElectrophoresisSource(const InputParameters & params);

  virtual Real computeMatrixContribution() override;

  virtual Real computeRightHandSideContribution() override;

protected:

  MooseLinearVariableFV<Real> & getelectric_potentialVariable(const std::string & vname);

  /// Model constants
  const Real _F;
  const Real _z;
  const Real _R;

  /// Temperature field
  const Moose::Functor<Real> & _T;

  /// Diffusion coefficient
  const Moose::Functor<Real> & _D;

  /// Pointer to the linear finite volume electric potential variable
  MooseLinearVariableFV<Real> & _electric_potential_var;

  /// The volume electric potential variable
  const std::vector<std::unique_ptr<NumericVector<Number>>> & _electric_potential_gradient;

  /// Cache for the volume electric potential variable number
  const unsigned int _electric_potential_var_num;

  /// Cache for the volume electric potential system number
  const unsigned int _electric_potential_sys_num;
};
