//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMomentumTimeDerivative.h"
#include "SystemBase.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMomentumTimeDerivative);

InputParameters
WCNSFV2PMomentumTimeDerivative::validParams()
{
  InputParameters params = WCNSFVMomentumTimeDerivative::validParams();
  params.addClassDescription(
      "Adds the time derivative term to the two-phase Navier-Stokes momentum equation.");
  params.addRequiredParam<MooseFunctorName>("fd", "The phase fraction functor.");
  return params;
}

WCNSFV2PMomentumTimeDerivative::WCNSFV2PMomentumTimeDerivative(const InputParameters & params)
  : WCNSFVMomentumTimeDerivative(params), _alpha(getFunctor<ADReal>("fd"))
{
}

void
WCNSFV2PMomentumTimeDerivative::gatherRCData(const Elem & elem)
{
  const auto elem_arg = makeElemArg(&elem);
  const auto state = determineState();
  const auto rho = _rho(elem_arg, state);
  const auto rho_dot = _rho_dot(elem_arg, state);
  const auto var = _var(elem_arg, state);
  const auto var_dot = _var.dot(elem_arg, state);
  const auto alpha = _alpha(elem_arg, state);
  const auto alpha_dot = _alpha.dot(elem_arg, state);

  const auto dof_number = elem.dof_number(_sys.number(), _var.number(), 0);
  mooseAssert(var.derivatives()[dof_number] == 1.,
              "This is an implicit assumption in our coefficient calculation.");

  const auto d_alpha_rho_dt = rho * alpha_dot + alpha * rho_dot;
  const auto strong_resid = d_alpha_rho_dt * var + rho * alpha * var_dot;

  // RC diagonal coefficients
  ADReal a = d_alpha_rho_dt;
  a += rho * alpha * var_dot.derivatives()[dof_number];

  const auto volume = _assembly.elementVolume(&elem);
  _rc_uo.addToA(&elem, _index, a * volume);
  addResidualAndJacobian(strong_resid * volume, dof_number);
}
