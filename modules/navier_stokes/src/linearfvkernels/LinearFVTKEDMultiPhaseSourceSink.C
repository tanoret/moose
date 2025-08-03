//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVTKEDMultiPhaseSourceSink.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "NavierStokesMethods.h"

registerMooseObject("NavierStokesApp", LinearFVTKEDMultiPhaseSourceSink);

InputParameters
LinearFVTKEDMultiPhaseSourceSink::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Elemental kernel to compute the production and destruction "
                             " terms of turbulent kinetic energy dissipation (TKED).");
  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::TKE, "Coupled turbulent kinetic energy.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "fluid density");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("wall_distance", "Distance to the nearest wall.");

  params.addRequiredParam<MooseFunctorName>("alpha", "The phase fraction.");

  params.addParam<Real>("C1_eps", 1.44, "First epsilon coefficient");
  params.addParam<Real>("C2_eps", 1.92, "Second epsilon coefficient");
  params.addParam<Real>("C_mu", 0.09, "Coupled turbulent kinetic energy closure.");
  params.addParam<Real>("C_pl", 10.0, "Production limiter constant multiplier.");
  params.addParam<Real>("Re_y_star", 60.0, "Cutoff Reynolds number for blending");
  params.addParam<Real>("delta_Re_y", 10.0, "Step in Reynolds number for blending");

  return params;
}

LinearFVTKEDMultiPhaseSourceSink::LinearFVTKEDMultiPhaseSourceSink(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(getFunctor<Real>("u")),
    _v_var(params.isParamValid("v") ? &(getFunctor<Real>("v")) : nullptr),
    _w_var(params.isParamValid("w") ? &(getFunctor<Real>("w")) : nullptr),
    _k(getFunctor<Real>(NS::TKE)),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _d(getFunctor<Real>("wall_distance")),
    _alpha(getFunctor<Real>("alpha")),
    _C1_eps(getParam<Real>("C1_eps")),
    _C2_eps(getParam<Real>("C2_eps")),
    _C_mu(getParam<Real>("C_mu")),
    _C_pl(getParam<Real>("C_pl")),
    _Re_y_star(getParam<Real>("Re_y_star")),
    _delta_Re_y(getParam<Real>("delta_Re_y"))
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied!");

  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied!");

  // Strain tensor term requires velocity gradients;
  // if (auto u = dynamic_cast<const MooseLinearVariableFV<Real> *>(&_u_var))
  //   requestVariableCellGradient(getParam<MooseFunctorName>("u"));
  // if (auto v = dynamic_cast<const MooseLinearVariableFV<Real> *>(_v_var))
  //   requestVariableCellGradient(getParam<MooseFunctorName>("v"));
  // if (auto w = dynamic_cast<const MooseLinearVariableFV<Real> *>(_w_var))
  //   requestVariableCellGradient(getParam<MooseFunctorName>("w"));
}

Real
LinearFVTKEDMultiPhaseSourceSink::computeMatrixContribution()
{

  // Useful variables
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const Real dist = _d(elem_arg, state);
  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  const Real k = _k(elem_arg, state);

  const Real Re_d = rho * std::sqrt(k) * dist / mu;

  if (Re_d <= _Re_y_star)
    // TKED value for near wall element will be directly assigned for this cell
    return 1.0;
  else
  {
    // Convenient definitions
    const auto epsilon = _var.getElemValue(*_current_elem_info, state);
    const Real alpha = _alpha(elem_arg, state);

    // Compute destruction
    const auto destruction = _C2_eps * rho * epsilon / k;

    // Assign to matrix (term gets multiplied by TKED)
    return alpha * destruction * _current_elem_volume;
  }
}

Real
LinearFVTKEDMultiPhaseSourceSink::computeRightHandSideContribution()
{
  // Useful variables
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const Real dist = _d(elem_arg, state);
  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  const Real k = _k(elem_arg, state);
  const Real alpha = _alpha(elem_arg, state);

  const Real Re_d = rho * std::sqrt(k) * dist / mu;

  if (Re_d <= _Re_y_star)
  {
    const Real cl = 0.42 * std::pow(_C_mu, -0.75);
    const Real l_eps = cl * dist * (1.0 - std::exp(-Re_d/2/cl));
    const Real destruction = std::pow(k, 1.5) / l_eps;

    // Assign the computed value of TKED for element near the wall
    return alpha * destruction;
  }
  else
  {
    // Convenient definitions
    const auto epsilon = _var.getElemValue(*_current_elem_info, state);

    const Real Re_d = rho * std::sqrt(k) * dist / mu;
    const Real A = _delta_Re_y / std::atanh(0.98);
    const Real lambda_blend = 0.5*(1.0 + std::tanh((Re_d-_Re_y_star)/A));

    const Real Te = k / epsilon;
    const Real T = std::max(Te, std::sqrt(mu/epsilon/rho));
    const Real mu_t_keps = rho * _C_mu * k * T;

    const Real mu_2_layer = mu * 0.42 * Re_d * std::pow(_C_mu, 0.25) * (1. - std::exp(-Re_d/70.0));

    const Real mu_t = lambda_blend * mu_t_keps + (1.0 - lambda_blend) * mu_2_layer;

    // Compute production of TKE
    const auto symmetric_strain_tensor_sq_norm =
        NS::computeShearStrainRateNormSquared<Real>(_u_var, _v_var, _w_var, elem_arg, state);
    Real production = mu_t* symmetric_strain_tensor_sq_norm;

    // Limit TKE production (needed for flows with stagnation zones)
    const Real production_limit = _C_pl * rho * epsilon;
    production = std::min(production, production_limit);

    // Compute production - recasted with mu_t definition to avoid division by epsilon
    const auto production_epsilon = _C1_eps * production * epsilon / k;

    // Assign to matrix (term gets multiplied by TKED)
    return alpha * production_epsilon * _current_elem_volume;
  }
}
