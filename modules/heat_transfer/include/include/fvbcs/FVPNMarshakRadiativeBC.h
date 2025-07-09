//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVDirichletBCBase.h"
#include "DenseMatrix.h"

/**
 * A class for turbulent kinetic energy dissipation rate inlet boundary conditions
 */
class FVPNMarshakRadiativeBC : public FVDirichletBCBase
{
public:
  static InputParameters validParams();
  FVPNMarshakRadiativeBC(const InputParameters & params);
  ADReal boundaryValue(const FaceInfo & fi, const Moose::StateArg & state) const override;

protected:

  // Order of the PN equation
  const unsigned int _n;

  /// Marshak matrix of coefficient
  DenseMatrix<Real> * _marshak_inlet_matrix;

  /// Coupled PN fluxes
  const Moose::Functor<ADReal> * _phi_0;
  const Moose::Functor<ADReal> * _phi_2;
  const Moose::Functor<ADReal> * _phi_4;
  const Moose::Functor<ADReal> * _phi_6;

  /// Coupled source
  const Moose::Functor<ADReal> * _boundary_source;

  /// Vector of fluxes
  std::vector<const Moose::Functor<ADReal> *> _phi_vector;

  // Orfer of the current PN equations
  unsigned int _N;
};
