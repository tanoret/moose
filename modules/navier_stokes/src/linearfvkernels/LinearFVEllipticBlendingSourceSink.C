//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVEllipticBlendingSourceSink.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "NavierStokesMethods.h"

registerMooseObject("MooseApp", LinearFVEllipticBlendingSourceSink);

InputParameters
LinearFVEllipticBlendingSourceSink::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Elemental kernel to compute the production and destruction for the elliptic blending function in the v2f mdoel");

  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::TKE, "Coupled turbulent kinetic energy.");
  params.addRequiredParam<MooseFunctorName>(NS::TKED, "Coupled turbulent kinetic energy dissipation rate.");
  params.addRequiredParam<MooseFunctorName>("theta_squared", "Wall-normal velocity fluctuations.");

  params.addRequiredParam<MooseFunctorName>(NS::density, "fluid density");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity.");

  params.addParam<Real>("C1", 1.4, "First blending production coefficient - production of wall-normal stesses.");
  params.addParam<Real>("C2", 0.3, "Second blending production coefficient - bulk scale.");
  params.addParam<Real>("C3", 6.0, "Third blending coefficient - destruction scaling.");
  params.addParam<Real>("C_mu_theta_2", 0.22, "C_mu_theta_2 closure coefficient.");

  return params;
}

LinearFVEllipticBlendingSourceSink::LinearFVEllipticBlendingSourceSink(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(getFunctor<Real>("u")),
    _v_var(params.isParamValid("v") ? &(getFunctor<Real>("v")) : nullptr),
    _w_var(params.isParamValid("w") ? &(getFunctor<Real>("w")) : nullptr),
    _k(getFunctor<Real>(NS::TKE)),
    _epsilon(getFunctor<Real>(NS::TKED)),
    _theta_squared(getFunctor<Real>("theta_squared")),
    _rho(getFunctor<Real>(NS::density)),
    _mu(getFunctor<Real>(NS::mu)),
    _C1(getParam<Real>("C1")),
    _C2(getParam<Real>("C2")),
    _C3(getParam<Real>("C3")),
    _C_mu_theta_2(getParam<Real>("C_mu_theta_2"))
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied!");

  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied!");
}

Real
LinearFVEllipticBlendingSourceSink::computeMatrixContribution()
{
  // Assign to matrix (term gets multiplied by f)
  return - 1.0 * _current_elem_volume;
}

Real
LinearFVEllipticBlendingSourceSink::computeRightHandSideContribution()
{
  // Convenient definitions
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const Real rho = _rho(elem_arg, state);
  const Real mu = _mu(elem_arg, state);
  const Real TKE = _k(elem_arg, state);
  const Real TKED = _epsilon(elem_arg, state);
  const Real theta_squared = _theta_squared(elem_arg, state);

  // Compute production of TKE
  const auto symmetric_strain_tensor_sq_norm =
      NS::computeShearStrainRateNormSquared<Real>(_u_var, _v_var, _w_var, elem_arg, state);

  // Bulk production
  const auto nu = mu / rho;
  const auto time_scale = std::max(TKE/TKED, _C3*std::sqrt(nu/TKED));
  const auto T1 = 1.0/time_scale*(_C1 - 1.0)*(2./3. - theta_squared/TKE);
  const auto T2 = _C2 * _C_mu_theta_2 * theta_squared * symmetric_strain_tensor_sq_norm / TKED;
  const auto T3 = (_C3 - 1.0) * theta_squared / TKE / time_scale;

  // Compute production
  const auto production = T1 + T2 + T3;

  // Assign to matrix (term gets multiplied by TKED)
  return production * _current_elem_volume;
}
