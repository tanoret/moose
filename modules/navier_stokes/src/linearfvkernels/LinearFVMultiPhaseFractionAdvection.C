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
  params.addRangeCheckedParam<Real>(
      "c_alpha", 0.0, "0.0<=c_alpha", "The compression velocity scaling constant.");
  params.addParam<MooseFunctorName>(NS::density, "The density.");
  params.addParam<bool>(
      "use_nonorthogonal_correction",
      false,
      "If the nonorthogonal correction should be used when computing the normal gradient.");
  params.addParam<bool>(
      "activate_mules",
      true,
      "Flag to aactivate CMULES limiting.");

  params += Moose::FV::advectedInterpolationParameter();

  MooseEnum limiterEnum(
      "min_mod vanLeer vanAlbada sou venkatakrishnan quick average upwind", "vanLeer");

  params.addParam<MooseEnum>("limiter_method",
                             limiterEnum,
                             "The limiter to use for the advected quantity. Options are "
                             "'min_mod', 'vanLeer', 'vanAlbada', 'sou', "
                             "'venkatakrishnan', 'quick', 'average', and "
                             "'upwind' with the default being 'upwind'.");
  return params;
}

LinearFVMultiPhaseFractionAdvection::LinearFVMultiPhaseFractionAdvection(
    const InputParameters & params)
  : LinearFVFluxKernel(params),
    _mass_flux_provider(getUserObject<RhieChowMassFluxMultiPhase>("rhie_chow_user_object")),
    _dim(_subproblem.mesh().dimension()),
    _c_alpha(getParam<Real>("c_alpha")),
    _rho(params.isParamValid(NS::density) ? &(getFunctor<Real>(NS::density)) : nullptr),
    _use_nonorthogonal_correction(getParam<bool>("use_nonorthogonal_correction")),
    _use_mules(getParam<bool>("activate_mules")),
    _advected_interp_coeffs(std::make_pair<Real, Real>(0, 0)),
    _total_adv_mass_face_flux(0.0),
    _limiter_method(getParam<MooseEnum>("limiter_method"))
{

  // if (_use_nonorthogonal_correction)
  //   _var.computeCellGradients();

  Moose::FV::setInterpolationMethod(*this, _advected_interp_method, "advected_interp_method");

  if (_c_alpha > 1e-42)
  {

    if (!_rho)
      paramError(NS::density,
                 "The density must be provided when compression velocity is activated by setting "
                 "c_alpha>1e-42");

    // Gradients are needed for compression velocity
    _var.computeCellGradients();
  }
}

Real
LinearFVMultiPhaseFractionAdvection::computeElemMatrixContribution()
{
  Real comp_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    comp_mass_flux = _lambda_f * computeCompressionVelocityMassFluxMatrixContribution();

  return (_advected_interp_coeffs.first * _total_adv_mass_face_flux + comp_mass_flux) *
         _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeNeighborMatrixContribution()
{
  Real comp_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    comp_mass_flux = _lambda_f * computeCompressionVelocityMassFluxMatrixContribution();

  return (_advected_interp_coeffs.second * _total_adv_mass_face_flux + comp_mass_flux) *
         _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeElemRightHandSideContribution()
{
  Real rhs = 0;
  if (_dim > 1 && _use_nonorthogonal_correction && _c_alpha > 1e-42)
    rhs += _lambda_f * computeCompressionVelocityMassFluxRHSContribution();

  Real tol_mass_flux = _advected_interp_coeffs.first * _total_adv_mass_face_flux;
  Real alpha_holo = this->getHighOrderFaceValue(_var) - this->getLowOrderFaceValue(_var);
  rhs -= tol_mass_flux * alpha_holo * _current_face_area;

  return rhs;
}

Real
LinearFVMultiPhaseFractionAdvection::computeNeighborRightHandSideContribution()
{
  Real rhs = 0;
  if (_dim > 1 && _use_nonorthogonal_correction && _c_alpha > 1e-42)
    rhs += _lambda_f * computeCompressionVelocityMassFluxRHSContribution();

  Real tol_mass_flux = _advected_interp_coeffs.second * _total_adv_mass_face_flux;
  Real alpha_holo = this->getHighOrderFaceValue(_var) - this->getLowOrderFaceValue(_var);
  rhs -= tol_mass_flux * alpha_holo * _current_face_area;

  return rhs;
}

Real
LinearFVMultiPhaseFractionAdvection::computeBoundaryMatrixContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  const auto boundary_value_matrix_contrib = adv_bc->computeBoundaryValueMatrixContribution();

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM) ? 1.0 : -1.0;

  // Adding compression velocity contribution
  Real compression_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    compression_mass_flux = computeCompressionVelocityMassFlux();

  return boundary_value_matrix_contrib * factor *
         (_total_adv_mass_face_flux + compression_mass_flux) * _current_face_area;
}

Real
LinearFVMultiPhaseFractionAdvection::computeBoundaryRHSContribution(
    const LinearFVBoundaryCondition & bc)
{
  const auto * const adv_bc = static_cast<const LinearFVAdvectionDiffusionBC *>(&bc);
  mooseAssert(adv_bc, "This should be a valid BC!");

  // We support internal boundaries too so we have to make sure the normal points always outward
  const auto factor = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM ? 1.0 : -1.0);

  const auto boundary_value_rhs_contrib = adv_bc->computeBoundaryValueRHSContribution();

  // Adding compression velocity contribution
  Real compression_mass_flux = 0.0;
  if (_c_alpha > 1e-42)
    compression_mass_flux = computeCompressionVelocityMassFlux();

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
  _low_order_face = makeFace(
      *_current_face_info, limiterType(_advected_interp_method), _total_adv_mass_face_flux);

  // MULES
  if(_use_mules)
  {
    const auto total_adv_volume_flux =
        _mass_flux_provider.getVolumetricFaceFlux(*face_info) * _current_face_area * _dt;

    const bool donor_is_elem = _low_order_face.elem_is_upwind;

    auto donor =
        donor_is_elem ? _low_order_face.makeElem() : _low_order_face.makeNeighbor();
    auto acceptor =
        donor_is_elem ? _low_order_face.makeNeighbor() : _low_order_face.makeElem();

    auto donor_info = donor_is_elem ? _current_face_info->elemInfo()
                                    : _current_face_info->neighborInfo();
    auto acceptor_info = donor_is_elem ? _current_face_info->neighborInfo()
                                      : _current_face_info->elemInfo();

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

    const auto donnor_capacity =
        MetaPhysicL::raw_value(_var(donor, determineState())) * donor_info->volume();
    const auto acceptor_capacity =
        (1.0 - MetaPhysicL::raw_value(_var(acceptor, determineState()))) * acceptor_info->volume();

    _lambda_f = std::max(std::min(std::min(1.0, donnor_capacity / std::abs(total_adv_volume_flux)),
                                  acceptor_capacity / std::abs(total_adv_volume_flux)),
                        1e-10);
  }
  else
    _lambda_f = 1.0;

}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFluxMatrixContribution()
{
  return computeCompressionVelocityMassFlux();
}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFluxRHSContribution()
{
  // Get the gradients from the adjacent cells
  const auto grad_elem = _var.gradSln(*_current_face_info->elemInfo());
  const auto & grad_neighbor = _var.gradSln(*_current_face_info->neighborInfo());

  // Interpolate the two gradients to the face
  const auto interp_coeffs =
      interpCoeffs(Moose::FV::InterpMethod::Average, *_current_face_info, true);

  const auto correction_vector =
      _current_face_info->normal() -
      1 / (_current_face_info->normal() * _current_face_info->eCN()) * _current_face_info->eCN();

  //  Compute compression velocity
  Real compression_velocity = computeCompressionVelocityMassFlux();

  return compression_velocity *
         (interp_coeffs.first * grad_elem + interp_coeffs.second * grad_neighbor) *
         correction_vector;
}

Real
LinearFVMultiPhaseFractionAdvection::computeCompressionVelocityMassFlux()
{

  // const auto alpha_f = this->getHighOrderFaceValue(_var);
  const auto alpha_f = this->getLowOrderFaceValue(_var);
  if (alpha_f <= 1e-12 || alpha_f >= 1.0 - 1e-12) // pure phase ⇒ no compression
    return 0.0;

  const auto grad = MetaPhysicL::raw_value(_var.gradSln(*_current_face_info->elemInfo()));
  const auto grad_mag = grad.norm();

  if (grad_mag < 1e-14)
    return 0.0;

  const Real Un = std::fabs(_total_adv_mass_face_flux);

  const Real u_c = _c_alpha * Un;

  const auto compression_dir = (grad / grad_mag) * _current_face_info->normal();

  const auto rho = (*_rho)(_low_order_face, determineState());

  const auto compression_mass_flux = rho * u_c * alpha_f * (1.0 - alpha_f) * compression_dir;

  return compression_mass_flux;
}

Real
LinearFVMultiPhaseFractionAdvection::getLowOrderFaceValue(MooseLinearVariableFV<Real> & variable)
{
  return MetaPhysicL::raw_value(variable(_low_order_face, determineState()));
}

Real
LinearFVMultiPhaseFractionAdvection::getHighOrderFaceValue(MooseLinearVariableFV<Real> & variable)
{

  //---------------------------------------------------------------------------
  // 1. Donor / acceptor bookkeeping
  //---------------------------------------------------------------------------
  const bool donor_is_elem = _low_order_face.elem_is_upwind;

  auto donor    = donor_is_elem ? _low_order_face.makeElem()
                                : _low_order_face.makeNeighbor();
  auto acceptor = donor_is_elem ? _low_order_face.makeNeighbor() 
                                : _low_order_face.makeElem();

  const auto * donor_info    = donor_is_elem ? _current_face_info->elemInfo()
                                             : _current_face_info->neighborInfo();
  const auto * acceptor_info = donor_is_elem ? _current_face_info->neighborInfo()
                                             : _current_face_info->elemInfo();

  // Handle boundaries where one side is missing
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
    
  //---------------------------------------------------------------------------
  // 2. Cell-centred values and donor increment \Delta \phi_P (=\nabla \phi_P \cdot dP)
  //---------------------------------------------------------------------------
  const Real phi_P = MetaPhysicL::raw_value(variable(donor, determineState()));
  const Real phi_N = MetaPhysicL::raw_value(variable(acceptor, determineState()));
  const auto  gradP = variable.gradSln(*donor_info);               // \phi_P
  const Point face_c   = _current_face_info->faceCentroid();
  const Point donor_c  = donor_info->centroid();
  const auto  dP       = face_c - donor_c;                         // dP
  const Real delta_P   = gradP * dP;                               // \nabla \phi_P \cdot dP
  
  //---------------------------------------------------------------------------
  // 3. Slope ratio  r  and limiter \psi(r)
  //---------------------------------------------------------------------------

  constexpr Real tiny = 1.0e-14;

  const Real deltaPhi = phi_N - phi_P;
  const Real r        = delta_P / (deltaPhi + (deltaPhi >= 0 ? tiny : -tiny));

  Real psi = 1.0;   // default = second-order (\psi=1)

  if (_limiter_method == "min_mod")
    psi = std::max(0.0, std::min(1.0, r));
  else if (_limiter_method == "vanLeer")
    psi = (r + std::fabs(r)) / (1.0 + std::fabs(r));
  else if (_limiter_method == "vanAlbada")
    psi = (r * r + r) / (r * r + 1.0);
  else if (_limiter_method == "quick")           // Koren QUICK
    psi = std::max(0.0,
                   std::min({ 2.0 / 3.0 * r + 1.0 / 6.0,  // bounded cubic
                              2.0 / 3.0,                  // upper plateaux
                              r }));                      // monotone
  else if (_limiter_method == "venkatakrishnan")
    psi = (r * r + 2.0 * r) / (r * r + r + 2.0);
  else if (_limiter_method == "average" || _limiter_method == "upwind")      // retain first-order
    psi = 0.0;
  // 'sou' (second-order upwind) and anything unrecognised fall back to \psi = 1

  //---------------------------------------------------------------------------
  // 4. High-order face value \phi_f  and storage for later access
  //---------------------------------------------------------------------------
  const Real phi_f = phi_P + psi * delta_P;  // Eq.  \phi_f = \phi_P + \psi(r)·\Delta \phi_P
  return phi_f;
}