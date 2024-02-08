//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMomentumPressure.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMomentumPressure);

InputParameters
WCNSFV2PMomentumPressure::validParams()
{
  InputParameters params = INSFVMomentumPressure::validParams();
  params.addClassDescription(
      "Introduces the coupled pressure term into the two-phase Navier-Stokes momentum equation.");
  params.addRequiredParam<MooseFunctorName>("fd", "The phase fraction functor.");
  return params;
}

WCNSFV2PMomentumPressure::WCNSFV2PMomentumPressure(const InputParameters & params)
  : INSFVMomentumPressure(params), _alpha(getFunctor<ADReal>("fd"))
{
}

ADReal
WCNSFV2PMomentumPressure::computeQpResidual()
{
  const auto elem = Moose::ElemArg{_current_elem, _correct_skewness};
  const auto state = determineState();
  return _alpha(elem, state) * _p.gradient(elem, state)(_index);
}
