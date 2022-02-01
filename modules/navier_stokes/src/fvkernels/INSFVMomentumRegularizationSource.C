//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumRegularizationSource.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", INSFVMomentumRegularizationSource);

InputParameters
INSFVMomentumRegularizationSource::validParams()
{
  InputParameters params = INSFVElementalKernel::validParams();
  params.addClassDescription(
      "Source part of the regularization kernel.");
  params.addRequiredCoupledVar("var_reg", "The coupled regularized variable");
  params.addParam<Real>("coef", 1.0, "Coefficent multiplier for the coupled force term.");
  return params;
}

INSFVMomentumRegularizationSource::INSFVMomentumRegularizationSource(const InputParameters & params)
  : FVElementalKernel(params),
    INSFVMomentumResidualObject(*this),
    _u_reg(adCoupledValue("var_reg")),
    _coef(getParam<Real>("coef"))
{
}

ADReal
INSFVMomentumRegularizationSource::computeQpResidual()
{
  return _coef * (_u[_qp] - _u_reg[_qp]);
}
