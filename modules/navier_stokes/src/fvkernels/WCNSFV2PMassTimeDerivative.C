//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMassTimeDerivative.h"

#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMassTimeDerivative);

InputParameters
WCNSFV2PMassTimeDerivative::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription(
      "Adds the mixture time derivative term to the two-phase weakly-compressible "
      "Navier-Stokes continuity equation.");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      NS::density + "_phases", {}, "The time derivative of the mixture density material property.");
  params.addRequiredParam<std::vector<MooseFunctorName>>(
      "fds", {}, "Volume fractions of the phases.");
  return params;
}

WCNSFV2PMassTimeDerivative::WCNSFV2PMassTimeDerivative(const InputParameters & params)
  : FVElementalKernel(params),
    _rho_names(getParam<std::vector<MooseFunctorName>>(NS::density + "_phases")),
    _alpha_names(getParam<std::vector<MooseFunctorName>>("fds"))
{
  _rho.reserve(_rho_names.size());
  for (const auto & name : _rho_names)
    _rho.push_back(&getFunctor<ADReal>(name));

  _alpha.reserve(_alpha_names.size());
  for (const auto & name : _alpha_names)
    _alpha.push_back(&getFunctor<ADReal>(name));

  if (_alpha.size() != _rho.size())
    paramError("fds",
               "The numebr of specied phase fractions is not the same that the number of provided "
               "densities.");
}

ADReal
WCNSFV2PMassTimeDerivative::computeQpResidual()
{
  const auto elem_arg = makeElemArg(_current_elem);
  const auto state = determineState();

  // Correcting for accumulation error on phase fraction
  Real phase_cell_value = 0.0;
  for (unsigned int i = 0; i < _alpha.size(); ++i)
    phase_cell_value += (*_alpha[i])(elem_arg, state).value();

  // Computing time derivative
  ADReal d_alpha_rho_dt = 0.0;
  for (unsigned int i = 0; i < _alpha.size(); ++i)
    d_alpha_rho_dt += (*_rho[i])(elem_arg, state) * _alpha[i]->dot(elem_arg, state) +
                      (*_alpha[i])(elem_arg, state) * _rho[i]->dot(elem_arg, state);

  return d_alpha_rho_dt / phase_cell_value;
}
