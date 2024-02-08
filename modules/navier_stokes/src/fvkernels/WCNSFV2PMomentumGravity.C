//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMomentumGravity.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMomentumGravity);

InputParameters
WCNSFV2PMomentumGravity::validParams()
{
  InputParameters params = INSFVMomentumGravity::validParams();
  params.addClassDescription(
      "Computes a body force due to gravity in two-phase Navier Stokes based simulations.");
  params.addRequiredParam<MooseFunctorName>("fd", "The phase fraction functor.");
  params.addRequiredParam<MooseFunctorName>(NS::density + "_mixture",
                                            "The phase fraction functor.");
  return params;
}

WCNSFV2PMomentumGravity::WCNSFV2PMomentumGravity(const InputParameters & params)
  : INSFVMomentumGravity(params),
    _alpha(getFunctor<ADReal>("fd")),
    _rho_mixture(getFunctor<ADReal>(NS::density + "_mixture"))
{
}

ADReal
WCNSFV2PMomentumGravity::computeQpResidual()
{
  const auto elem = makeElemArg(_current_elem);
  const auto state = determineState();
  return -_alpha(elem, state) * (_rho(elem, state) - _rho_mixture(elem, state)) * _gravity(_index);
}
