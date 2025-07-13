//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVMultiPhaseTimeDerivative.h"
#include "NS.h"
#include "TimeIntegrator.h"

registerMooseObject("NavierStokesApp", LinearFVMultiPhaseTimeDerivative);

InputParameters
LinearFVMultiPhaseTimeDerivative::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Represents the matrix and right hand side contributions of a "
                             "time derivative term for two-phase flows.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "The density of the phase.");
  params.addRequiredParam<MooseFunctorName>("alpha", "The phase fraction.");
  return params;
}

LinearFVMultiPhaseTimeDerivative::LinearFVMultiPhaseTimeDerivative(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _rho(getFunctor<Real>(NS::density)),
    _alpha(getFunctor<Real>("alpha")),
    _time_integrator(_sys.getTimeIntegrator(_var_num)),
    _rho_alpha_history(_time_integrator.numStatesRequired(), 0.0),
    _state_args(_time_integrator.numStatesRequired(), determineState())
{
  // In case we need older states
  for (const auto i : index_range(_state_args))
    _state_args[i] = Moose::StateArg(i, Moose::SolutionIterationType::Time);
}

Real
LinearFVMultiPhaseTimeDerivative::computeMatrixContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto rho_alpha = _rho(elem_arg, _state_args[0]) * _alpha(elem_arg, _state_args[0]);
  return _time_integrator.timeDerivativeMatrixContribution(rho_alpha) *
         _current_elem_volume;
}

Real
LinearFVMultiPhaseTimeDerivative::computeRightHandSideContribution()
{
  return _time_integrator.timeDerivativeRHSContribution(_dof_id, _rho_alpha_history) *
         _current_elem_volume;
}

void
LinearFVMultiPhaseTimeDerivative::setCurrentElemInfo(const ElemInfo * elem_info)
{
  LinearFVElementalKernel::setCurrentElemInfo(elem_info);

  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  for (const auto i : index_range(_rho_alpha_history))
    _rho_alpha_history[i] = _rho(elem_arg, _state_args[i]) * _alpha(elem_arg, _state_args[i]);
}
