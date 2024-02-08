//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMassAdvection.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMassAdvection);

InputParameters
WCNSFV2PMassAdvection::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params += Moose::FV::interpolationParameters();

  params.addRequiredParam<std::vector<MooseFunctorName>>(
      NS::density + "_phases", "The time derivative of the mixture density material property.");
  params.addRequiredParam<std::vector<MooseFunctorName>>("fds", "Volume fractions of the phases.");
  params.addRequiredParam<std::vector<UserObjectName>>(
      "rhie_chow_user_objects",
      "A vector with the Rhie-Chow user objects for each transported phase.");

  params.set<unsigned short>("ghost_layers") = 2;
  params.suppressParameter<bool>("force_boundary_execution");

  return params;
}

WCNSFV2PMassAdvection::WCNSFV2PMassAdvection(const InputParameters & params)
  : FVFluxKernel(params),
    _rho_names(getParam<std::vector<MooseFunctorName>>(NS::density + "_phases")),
    _alpha_names(getParam<std::vector<MooseFunctorName>>("fds")),
    _rc_names(getParam<std::vector<UserObjectName>>("rhie_chow_user_objects"))
{
  const bool need_more_ghosting =
      Moose::FV::setInterpolationMethods(*this, _advected_interp_method, _velocity_interp_method);
  if (need_more_ghosting && _tid == 0)
  {
    adjustRMGhostLayers(std::max((unsigned short)(3), _pars.get<unsigned short>("ghost_layers")));
    getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")
        ->setErrorOnJacobianNonzeroReallocation(false);
  }

  auto param_check = [&params, this](const auto & param_name)
  {
    if (params.isParamSetByUser(param_name))
      paramError(param_name,
                 "This parameter is not honored by WCNSFV2PMassAdvections like '",
                 name(),
                 "'");
  };

  param_check("force_boundary_execution");

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

  _rc_velocity_providers.reserve(_rc_names.size());
  for (const auto & name : _rc_names)
    _rc_velocity_providers.push_back(&getUserObjectByName<RhieChowInterpolatorBase>(name));

  if (_rc_velocity_providers.size() != _rho.size())
    paramError(
        "rhie_chow_user_objects",
        "The numebr of Rhie Chow user objects provided is not the same that the number of provided "
        "densities.");
}

void
WCNSFV2PMassAdvection::initialSetup()
{
  INSFVBCInterface::initialSetup(*this);
}

bool
WCNSFV2PMassAdvection::skipForBoundary(const FaceInfo & fi) const
{
  // Boundaries to avoid come first, since they are always obeyed
  if (avoidBoundary(fi))
    return true;

  // We're not on a boundary, so technically we're not skipping a boundary
  if (!onBoundary(fi))
    return false;

  // Selected boundaries to force
  for (const auto bnd_to_force : _boundaries_to_force)
    if (fi.boundaryIDs().count(bnd_to_force))
      return false;

  // If we have flux bcs then we do skip
  const auto & [have_flux_bcs, flux_bcs] = _var.getFluxBCs(fi);
  libmesh_ignore(have_flux_bcs);
  for (const auto * const flux_bc : flux_bcs)
    // If we have something like an average-value pressure constraint on a flow boundary, then we
    // still want to execute this advection kernel on the boundary to ensure we're enforcing local
    // conservation (mass in this example)
    if (!dynamic_cast<const FVBoundaryScalarLagrangeMultiplierConstraint *>(flux_bc))
      return true;

  // If we have a flow boundary without a replacement flux BC, then we must not skip. Mass and
  // momentum are transported via advection across boundaries
  for (const auto bc_id : fi.boundaryIDs())
    if (_flow_boundaries.find(bc_id) != _flow_boundaries.end())
      return false;

  // If not a flow boundary, then there should be no advection/flow in the normal direction, e.g. we
  // should not contribute any advective flux
  return true;
}

ADReal
WCNSFV2PMassAdvection::computeQpResidual()
{
  const auto state = determineState();

  std::vector<VectorValue<ADReal>> rc_velocities;
  rc_velocities.reserve(_rc_velocity_providers.size());
  for (unsigned int i = 0; i < _rc_velocity_providers.size(); ++i)
    rc_velocities.push_back(_rc_velocity_providers[i]->getVelocity(
        _velocity_interp_method, *_face_info, state, _tid, hasMaterialTimeDerivative()));

  // Correcting for accumulation error on phase fraction
  Real phase_cell_value = 0.0;
  for (unsigned int i = 0; i < _alpha.size(); ++i)
  {
    const auto rc_v = rc_velocities[i];
    phase_cell_value += (*_alpha[i])(makeFace(*_face_info,
                                              limiterType(_advected_interp_method),
                                              MetaPhysicL::raw_value(rc_v) * _normal > 0),
                                     state)
                            .value();
  }

  ADReal face_two_face_flux = 0.0;
  for (unsigned int i = 0; i < _alpha.size(); ++i)
  {
    const auto rc_v = rc_velocities[i];
    const auto alpha_face = (*_alpha[i])(makeFace(*_face_info,
                                                  limiterType(_advected_interp_method),
                                                  MetaPhysicL::raw_value(rc_v) * _normal > 0),
                                         state);
    const auto rho_face = (*_rho[i])(makeFace(*_face_info,
                                              limiterType(_advected_interp_method),
                                              MetaPhysicL::raw_value(rc_v) * _normal > 0),
                                     state);
    face_two_face_flux += (_normal * rc_v) * alpha_face * rho_face;
  }

  return face_two_face_flux / phase_cell_value;
}
