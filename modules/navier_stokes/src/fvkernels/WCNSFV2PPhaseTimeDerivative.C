//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PPhaseTimeDerivative.h"
#include "SystemBase.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFV2PPhaseTimeDerivative);

InputParameters
WCNSFV2PPhaseTimeDerivative::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "Adds the phase time derivative term to the two-phase, weakly-compressible "
      "Navier-Stokes momentum equation.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "The density material property");
  params.addRequiredParam<MooseFunctorName>(NS::time_deriv(NS::density),
                                            "The time derivative of the density material property");
  return params;
}

WCNSFV2PPhaseTimeDerivative::WCNSFV2PPhaseTimeDerivative(const InputParameters & params)
  : FVElementalKernel(params),
    _rho(getFunctor<ADReal>(NS::density)),
    _rho_dot(getFunctor<ADReal>(NS::time_deriv(NS::density)))
{
}

ADReal
WCNSFV2PPhaseTimeDerivative::computeQpResidual()
{
  const auto elem_arg = makeElemArg(_current_elem);
  const auto state = determineState();
  const auto rho = _rho(elem_arg, state);
  const auto rho_dot = _rho_dot(elem_arg, state);
  const auto alpha_var = _var(elem_arg, state);
  const auto alpha_dot = _var.dot(elem_arg, state);

  return rho_dot * alpha_var + rho * alpha_dot;
}
