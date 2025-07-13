//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVMultiPhaseFractionAdvection.h"
#include "MooseLinearVariableFV.h"
#include "NSFVUtils.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", LinearFVMultiPhaseFractionAdvection);

InputParameters
LinearFVMultiPhaseFractionAdvection::validParams()
{
  InputParameters params = LinearFVFluxKernel::validParams();
  params.addClassDescription("Represents the matrix and right hand side contributions of an "
                             "advection term for a the mass-weighted phase fraction.");
  params.addRequiredParam<UserObjectName>(
      "rhie_chow_user_object",
      "The rhie-chow user-object which is used to determine the face velocity.");
  params.addRangeCheckedParam<Real>("c_alpha", 0.0, "0.0<=c_alpha", "The compression velocity scaling constant.");
  params.addParam<MooseFunctorName>(NS::density, "The density.");
  params.addParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addParam<MooseFunctorName>("u_mixture", "The mixture velocity in the x direction.");
  params.addParam<MooseFunctorName>("v_mixture", "The mixture velocity in the y direction.");
  params.addParam<MooseFunctorName>("w_mixture", "The mixture velocity in the z direction.");
  params.addParam<bool>("use_nonorthogonal_correction", false,
                        "If the nonorthogonal correction should be used when computing the normal gradient.");

  params += Moose::FV::advectedInterpolationParameter();

  MooseEnum limiterEnum("average upwind sou min_mod vanLeer quick venkatakrishnan skewness-corrected", "vanLeer");

  params.addParam<MooseEnum>("limiter_method",
                             limiterEnum,
                             "The limiter to use for the advected quantity. Options are "
                             "'upwind', 'average', 'sou' (for second-order upwind), 'min_mod', "
                             "'vanLeer', 'quick', 'venkatakrishnan', and "
                             "'skewness-corrected' with the default being 'upwind'.");
  return params;
}

LinearFVMultiPhaseFractionAdvection::LinearFVMultiPhaseFractionAdvection(const InputParameters & params)
  : LinearFVFluxKernel(params),
    _mass_flux_provider(getUserObject<RhieChowMassFluxMultiPhase>("rhie_chow_user_object")),
    _dim(_subproblem.mesh().dimension()),
    _c_alpha(getParam<Real>("c_alpha")),
    _rho(params.isParamValid(NS::density) ? &(getFunctor<Real>(NS::density)) : nullptr),
    _u_var(params.isParamValid("u") ? &(getFunctor<Real>("u")) : nullptr),
    _v_var(params.isParamValid("v") ? &(getFunctor<Real>("v")) : nullptr),
    _w_var(params.isParamValid("w") ? &(getFunctor<Real>("w")) : nullptr),
    _u_var_mixture(params.isParamValid("u_mixture") ? &(getFunctor<Real>("u_mixture")) : nullptr),
    _v_var_mixture(params.isParamValid("v_mixture") ? &(getFunctor<Real>("v_mixture")) : nullptr),
    _w_var_mixture(params.isParamValid("w_mixture") ? &(getFunctor<Real>("w_mixture")) : nullptr),
    _use_nonorthogonal_correction(getParam<bool>("use_nonorthogonal_correction")),
    _advected_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _total_adv_mass_face_flux(0.0)
{

  if (_use_nonorthogonal_correction)
    _var.computeCellGradients();

  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");
  Moose::FV::setInterpolationMethod(*this, _limiter_method, "limiter_method");

  if(_c_alpha > 1e-42)
  {

    if(!_rho)
      paramError(NS::density,
                "The density must be provided when compression velocity is activated by setting c_alpha>1e-42");

    if(!_u_var)
      paramError("u",
                "The velocity 'u' must be provided when compression velocity is activated by setting c_alpha>1e-42");

    if (_dim >= 2 && !_v_var)
      paramError("v",
                "In two or more dimensions, "
                "the velocity 'v' must be provided when compression velocity is activated by setting c_alpha>1e-42!");

    if (_dim >= 3 && !_w_var)
      paramError("w",
                "In three dimensions, "
                "the velocity 'w' must be provided when compression velocity is activated by setting c_alpha>1e-42!");

    if(!_u_var_mixture)
      paramError("u_mixture",
                 "The velocity mixture 'u' must be provided when compression velocity is activated by setting c_alpha>1e-42");

    if (_dim >= 2 && !_v_var_mixture)
      paramError("v_mixture",
                 "In two or more dimensions, "
                 "the velocity mixture 'v' must be provided when compression velocity is activated by setting c_alpha>1e-42!");

    if (_dim >= 3 && !_w_var_mixture)
      paramError("w_mixture",
                 "In three dimensions, "
                 "the velocity mixture 'w' must be provided when compression velocity is activated by setting c_alpha>1e-42!");
  }
}

Real
LinearFVMultiPhaseFractionAdvection::computeElemMatrixContribution()
{
  Real comp_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
    comp_mass_flux = _lambda_f * computeCompressionVelocityMassFluxMatrixContribution();

  return (_advected_interp_coeffs.first * _total_adv_mass_face_flux + comp_mass_flux) * _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeNeighborMatrixContribution()
{
  Real comp_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
    comp_mass_flux = _lambda_f * computeCompressionVelocityMassFluxMatrixContribution();

  return (_advected_interp_coeffs.second * _total_adv_mass_face_flux + comp_mass_flux) * _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeElemRightHandSideContribution()
{
  Real rhs = 0;
  if (_dim > 1 && _use_nonorthogonal_correction && _c_alpha > 1e-42)
    rhs += _lambda_f * computeCompressionVelocityMassFluxRHSContribution();

  Real comp_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
    comp_mass_flux = _lambda_f * computeCompressionVelocityMassFluxMatrixContribution();

  Real tol_mass_flux = _advected_interp_coeffs.first * _total_adv_mass_face_flux + comp_mass_flux;
  Real alpha_holo = MetaPhysicL::raw_value(_var(_high_order_face, determineState()))
                    - MetaPhysicL::raw_value(_var(_low_order_face, determineState()));
  rhs -= tol_mass_flux * alpha_holo * _current_face_area;

  return rhs;

}

Real
LinearFVMultiPhaseFractionAdvection::computeNeighborRightHandSideContribution()
{
  Real rhs = 0;
  if (_dim > 1 && _use_nonorthogonal_correction && _c_alpha > 1e-42)
    rhs += _lambda_f * computeCompressionVelocityMassFluxRHSContribution();

  Real comp_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
    comp_mass_flux = _lambda_f * computeCompressionVelocityMassFluxMatrixContribution();

  Real tol_mass_flux = _advected_interp_coeffs.second * _total_adv_mass_face_flux + comp_mass_flux;
  Real alpha_holo = MetaPhysicL::raw_value(_var(_high_order_face, determineState()))
                    - MetaPhysicL::raw_value(_var(_low_order_face, determineState()));
  rhs -= tol_mass_flux * alpha_holo * _current_face_area;

  return rhs;
}

Real
LinearFVMultiPhaseFractionAdvection::computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto boundary_value_matrix_contrib = adv_bc->computeBoundaryValueMatrixContribution();

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  // Adding compression velocity contribution
  Real compression_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
  {
    const auto face_arg = singleSidedFaceArg(_current_face_info);
    auto grad_alpha = adv_bc->computeBoundaryGradientMatrixContribution();
    compression_mass_flux = computeCompressionVelocityMassFlux(face_arg, grad_alpha);
  }

  return boundary_value_matrix_contrib * factor *
         (_total_adv_mass_face_flux + compression_mass_flux) * _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0);

  const auto boundary_value_rhs_contrib = adv_bc->computeBoundaryValueRHSContribution();

  // Adding compression velocity contribution
  Real compression_mass_flux = 0.0;
  if(_c_alpha > 1e-42)
  {
    const auto face_arg = singleSidedFaceArg(_current_face_info);
    auto grad_alpha = adv_bc->computeBoundaryGradientRHSContribution();
    compression_mass_flux = computeCompressionVelocityMassFlux(face_arg, grad_alpha);
  }

  return -boundary_value_rhs_contrib * factor *
         (_total_adv_mass_face_flux + compression_mass_flux) * _current_face_area;
}

void
LinearFVMultiPhaseFractionAdvection::setupFaceData(const FaceInfo * face_info)
{
  LinearFVFluxKernel::setupFaceData(face_info);

  // Caching the velocity on the face which will be reused in the advection term's matrix and right
  // hand side contributions
  _total_adv_mass_face_flux = _mass_flux_provider.getUnweightedMassFlux(*face_info);

  // Caching the interpolation coefficients so they will be reused for the matrix and right hand
  // side terms
  _advected_interp_coeffs =
      interpCoeffs(_advected_interp_method, *_current_face_info, true, _total_adv_mass_face_flux);

  // Store low-order face 
  _low_order_face = makeFace(*_current_face_info,
                             limiterType(_advected_interp_method),
                             _total_adv_mass_face_flux);

  // Store higher-order face 
  _high_order_face = makeFace(*_current_face_info,
                              limiterType(_limiter_method),
                              _total_adv_mass_face_flux);

  // CMULES
  const auto total_adv_volume_flux = _mass_flux_provider.getVolumetricFaceFlux(*face_info) * _current_face_area * _dt;

  auto donor    = _total_adv_mass_face_flux > 0 ?
                        _low_order_face.makeElem() :
                        _low_order_face.makeNeighbor();
  auto acceptor = _total_adv_mass_face_flux > 0 ?
                        _low_order_face.makeNeighbor() : 
                        _low_order_face.makeElem();

  auto donor_info    = _total_adv_mass_face_flux > 0 ?
                             _current_face_info->elemInfo() :
                             _current_face_info->neighborInfo();
  auto acceptor_info = _total_adv_mass_face_flux > 0 ?
                             _current_face_info->neighborInfo() :
                             _current_face_info->elemInfo();
  
  if (!donor_info)
  {
    donor_info = acceptor_info;
    donor = acceptor;
  }

  if (!acceptor_info)
  {
    acceptor_info = donor_info;
    acceptor = donor;
  }

  const auto donnor_capacity = MetaPhysicL::raw_value(_var(donor, determineState())) * donor_info->volume();
  const auto acceptor_capacity = (1.0 - MetaPhysicL::raw_value(_var(acceptor, determineState()))) * acceptor_info->volume();

  _lambda_f = std::max(
              std::min(
              std::min(
                  1.0, donnor_capacity/std::abs(total_adv_volume_flux)
                ),
                  acceptor_capacity/std::abs(total_adv_volume_flux)
                ),
                  1e-10
              );
}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFluxMatrixContribution()
{
  const auto face_arg = makeCDFace(*_current_face_info);

  // If we requested nonorthogonal correction, we use the normal component of the
  // cell to face vector.
  const auto d = _use_nonorthogonal_correction
                      ? std::abs(_current_face_info->dCN() * _current_face_info->normal())
                      : _current_face_info->dCNMag();

  return computeCompressionVelocityMassFlux(face_arg, 1.0/d);
}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFluxRHSContribution()
{
  const auto face_arg = makeCDFace(*_current_face_info);

  // Get the gradients from the adjacent cells
  const auto grad_elem = _var.gradSln(*_current_face_info->elemInfo());
  const auto & grad_neighbor = _var.gradSln(*_current_face_info->neighborInfo());

  // Interpolate the two gradients to the face
  const auto interp_coeffs =
      interpCoeffs(Moose::FV::InterpMethod::Average, *_current_face_info, true);

  const auto correction_vector =
      _current_face_info->normal() -
      1 / (_current_face_info->normal() * _current_face_info->eCN()) *
          _current_face_info->eCN();

  //  Compute compression velocity
  Real compression_velocity = computeCompressionVelocityMassFlux(face_arg, 1.0);

  return compression_velocity * (interp_coeffs.first * grad_elem + interp_coeffs.second * grad_neighbor) * correction_vector;
}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFlux(const Moose::FaceArg & face_arg,
                                                                   const Real & grad_alpha)
{
  const auto state = determineState();

  const auto alpha = MetaPhysicL::raw_value(_var(_high_order_face, state));
  const auto grad_alpha_norm = MetaPhysicL::raw_value(_var.gradSln(*_current_face_info->elemInfo()).norm());

  RealVectorValue velocity((*_u_var)(face_arg, state));
  if (_v_var)
    velocity(1) = (*_v_var)(face_arg, state);
  if (_w_var)
    velocity(2) = (*_w_var)(face_arg, state);

  RealVectorValue velocity_mixture((*_u_var_mixture)(face_arg, state));
  if (_v_var)
    velocity_mixture(1) = (*_v_var_mixture)(face_arg, state);
  if (_w_var)
    velocity_mixture(2) = (*_w_var_mixture)(face_arg, state);

  const auto rel_velocity = velocity - velocity_mixture;

  const auto compression_vel = _c_alpha * rel_velocity.norm() * grad_alpha / (grad_alpha_norm + 0.1);

  const auto rho = (*_rho)(face_arg, state);

  const auto compression_mass_flux = rho * alpha * (1. - alpha) * compression_vel;

  return compression_mass_flux;

}

