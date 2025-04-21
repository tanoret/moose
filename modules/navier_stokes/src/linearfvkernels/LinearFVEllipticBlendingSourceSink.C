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
  params.addRequiredParam<MooseFunctorName>(NS::mu_t, "Turbulent viscosity.");

  params.addParam<std::vector<BoundaryName>>("walls", {}, "Boundaries that correspond to solid walls.");
  params.addParam<bool>(
      "linearized_model",
      true,
      "Boolean to determine if the problem should be used in a linear or nonlinear solve");
  MooseEnum wall_treatment("eq_newton eq_incremental eq_linearized neq", "neq");
  params.addParam<MooseEnum>("wall_treatment",
                            wall_treatment,
                            "The method used for computing the wall functions "
                            "'eq_newton', 'eq_incremental', 'eq_linearized', 'neq'");

  params.addParam<Real>("C1", 1.4, "First blending production coefficient - production of wall-normal stesses.");
  params.addParam<Real>("C2", 0.3, "Second blending production coefficient - bulk scale.");
  params.addParam<Real>("C3", 6.0, "Third blending coefficient - destruction scaling.");
  params.addParam<Real>("C_mu", 0.09, "Coupled turbulent kinetic energy closure.");
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
    _mu_t(getFunctor<Real>(NS::mu_t)),
    _wall_boundary_names(getParam<std::vector<BoundaryName>>("walls")),
    _linearized_model(getParam<bool>("linearized_model")),
    _wall_treatment(getParam<MooseEnum>("wall_treatment").getEnum<NS::WallTreatmentEnum>()),
    _C1(getParam<Real>("C1")),
    _C2(getParam<Real>("C2")),
    _C3(getParam<Real>("C3")),
    _C_mu(getParam<Real>("C_mu")),
    _C_mu_theta_2(getParam<Real>("C_mu_theta_2"))
{
  if (_dim >= 2 && !_v_var)
    paramError("v", "In two or more dimensions, the v velocity must be supplied!");

  if (_dim >= 3 && !_w_var)
    paramError("w", "In three or more dimensions, the w velocity must be supplied!");
}

void
LinearFVEllipticBlendingSourceSink::initialSetup()
{
  LinearFVElementalKernel::initialSetup();
  NS::getWallBoundedElements(
      _wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _wall_bounded);
  NS::getWallDistance(_wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _dist);
  NS::getElementFaceArgs(_wall_boundary_names, _fe_problem, _subproblem, blockIDs(), _face_infos);
}

Real
LinearFVEllipticBlendingSourceSink::computeMatrixContribution()
{
  if (_wall_bounded.find(_current_elem_info->elem()) != _wall_bounded.end())
    // TKED value for near wall element will be directly assigned for this cell
    return 1.0;
  else
  {
    // Assign to matrix (term gets multiplied by f)
    return - 1.0 * _current_elem_volume;
  }
}

Real
LinearFVEllipticBlendingSourceSink::computeRightHandSideContribution()
{
  if (_wall_bounded.find(_current_elem_info->elem()) != _wall_bounded.end())
  {
    // Convenient definitions
    const auto state = determineState();
    const auto elem_arg = makeElemArg(_current_elem_info->elem());
    const Real rho = _rho(elem_arg, state);
    const Real mu = _mu(elem_arg, state);
    const Real TKE = _k(elem_arg, state);

    // Convenient variables
    Real f_wall = 0.0;
    std::vector<Real> y_plus_vec, velocity_grad_norm_vec;
    Real tot_weight = 0.0;

    // Get velocity vector
    RealVectorValue velocity(_u_var(elem_arg, state));
    if (_v_var)
      velocity(1) = (*_v_var)(elem_arg, state);
    if (_w_var)
      velocity(2) = (*_w_var)(elem_arg, state);

    // Get near wall faceInfo and distances from cell center to every wall
    const auto & face_info_vec = libmesh_map_find(_face_infos, _current_elem_info->elem());
    const auto & distance_vec = libmesh_map_find(_dist, _current_elem_info->elem());
    mooseAssert(distance_vec.size(), "Should have found a distance vector");
    mooseAssert(distance_vec.size() == face_info_vec.size(),
                "Should be as many distance vectors as face info vectors");

    // Update y+ and wall face cell
    for (unsigned int i = 0; i < distance_vec.size(); i++)
    {
      const auto distance = distance_vec[i];
      mooseAssert(distance > 0, "Should be at a non-zero distance");

      Real y_plus;
      if (_wall_treatment == NS::WallTreatmentEnum::NEQ) // Non-equilibrium / Non-iterative
        y_plus = distance * std::sqrt(std::sqrt(_C_mu) * TKE) * rho / mu;
      else // Equilibrium / Iterative
      {
        const auto parallel_speed = NS::computeSpeed<Real>(
            velocity - velocity * face_info_vec[i]->normal() * face_info_vec[i]->normal());
        y_plus = NS::findyPlus<Real>(mu, rho, std::max(parallel_speed, 1e-10), distance);
      }

      y_plus_vec.push_back(y_plus);
      tot_weight += 1.0;
    }

    // Compute near wall epsilon value
    for (const auto i : index_range(y_plus_vec))
    {
      const auto y_plus = y_plus_vec[i];
      if (y_plus > 11.25)
      {
        const Real TKED = _k(elem_arg, state);
        const Real theta_squared = _theta_squared(elem_arg, state);
        f_wall += _C3 * theta_squared * TKED / (Utility::pow<2>(TKE)) / (std::sqrt(_C_mu) * TKE) / tot_weight;
      }
    }

    // Assign the computed value of TKED for element near the wall
    return f_wall;
  }
  else
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
    // const auto v2fAlpha = 1.0/time_scale*((_C1 - _C3)*theta_squared - 2.0/3.0*TKE*(_C1 - 1.0));

    // Compute production
    const auto production = T1 + T2 + T3;
    // const auto production = v2fAlpha - _C2*(_mu_t(elem_arg, state) * symmetric_strain_tensor_sq_norm / rho);

    // Assign to matrix (term gets multiplied by TKED)
    return - production * _current_elem_volume;
  }
}
