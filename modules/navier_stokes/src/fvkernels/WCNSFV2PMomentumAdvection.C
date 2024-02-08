//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMomentumAdvection.h"
#include "NS.h"
#include "FVUtils.h"
#include "INSFVRhieChowInterpolator.h"
#include "SystemBase.h"
#include "NSFVUtils.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMomentumAdvection);

InputParameters
WCNSFV2PMomentumAdvection::validParams()
{
  InputParameters params = INSFVMomentumAdvection::validParams();
  params.addRequiredParam<MooseFunctorName>("fd", "The phase fraction functor.");
  params.addClassDescription("Object for two-phase advecting momentum, e.g. alpha*rho*u");
  return params;
}

WCNSFV2PMomentumAdvection::WCNSFV2PMomentumAdvection(const InputParameters & params)
  : INSFVMomentumAdvection(params),
    _alpha(getFunctor<ADReal>("fd"))
{
}

void
WCNSFV2PMomentumAdvection::initialSetup()
{
  INSFVAdvectionKernel::initialSetup();

  _rho_u = std::make_unique<PiecewiseByBlockLambdaFunctor<ADReal>>(
      "rho_u",
      [this](const auto & r, const auto & t) -> ADReal
      { return _rho(r, t) * _alpha(r, t) * _var(r, t) / epsilon()(r, t); },
      std::set<ExecFlagType>({EXEC_ALWAYS}),
      _mesh,
      this->blockIDs());
}

void
WCNSFV2PMomentumAdvection::computeResidualsAndAData(const FaceInfo & fi)
{
  mooseAssert(!skipForBoundary(fi),
              "We shouldn't get in here if we're supposed to skip for a boundary");

  _face_info = &fi;
  _normal = fi.normal();
  _face_type = fi.faceType(std::make_pair(_var.number(), _var.sys().number()));
  const auto state = determineState();

  using namespace Moose::FV;

  const bool correct_skewness = _advected_interp_method == InterpMethod::SkewCorrectedAverage;
  const bool on_boundary = onBoundary(fi);

  _elem_residual = 0, _neighbor_residual = 0, _ae = 0, _an = 0;

  const auto v_face = velocity();
  const auto vdotn = _normal * v_face;

  if (on_boundary)
  {
    const auto ssf = singleSidedFaceArg();
    const Elem * const sided_elem = ssf.face_side;
    const auto dof_number = sided_elem->dof_number(_sys.number(), _var.number(), 0);
    const auto rho_alpha_face = _rho(ssf, state) * _alpha(ssf, state);
    const auto eps_face = epsilon()(ssf, state);
    const auto u_face = _var(ssf, state);
    const Real d_u_face_d_dof = u_face.derivatives()[dof_number];
    const auto coeff = vdotn * rho_alpha_face / eps_face;

    if (sided_elem == fi.elemPtr())
    {
      _ae = coeff * d_u_face_d_dof;
      _elem_residual = coeff * u_face;
      if (_approximate_as)
        _ae = _cs / fi.elem().n_sides() * rho_alpha_face / eps_face;
    }
    else
    {
      _an = -coeff * d_u_face_d_dof;
      _neighbor_residual = -coeff * u_face;
      if (_approximate_as)
        _an = _cs / fi.neighborPtr()->n_sides() * rho_alpha_face / eps_face;
    }
  }
  else // we are an internal fluid flow face
  {
    const bool elem_is_upwind = MetaPhysicL::raw_value(v_face) * _normal > 0;
    const Moose::FaceArg advected_face_arg{
        &fi, limiterType(_advected_interp_method), elem_is_upwind, correct_skewness, nullptr};
    if (const auto [is_jump, eps_elem_face, eps_neighbor_face] =
            NS::isPorosityJumpFace(epsilon(), fi, state);
        is_jump)
    {
      const auto rho_alpha_face = _rho(advected_face_arg, state) * _alpha(advected_face_arg, state);

      // We set the + and - sides of the superficial velocity equal to the interpolated value
      const auto & var_elem_face = v_face(_index);
      const auto & var_neighbor_face = v_face(_index);

      const auto elem_dof_number = fi.elem().dof_number(_sys.number(), _var.number(), 0);
      const auto neighbor_dof_number = fi.neighbor().dof_number(_sys.number(), _var.number(), 0);

      const auto d_var_elem_face_d_elem_dof = var_elem_face.derivatives()[elem_dof_number];
      const auto d_var_neighbor_face_d_neighbor_dof =
          var_neighbor_face.derivatives()[neighbor_dof_number];

      // We allow a discontintuity in the advective momentum flux at the jump face. Mainly the
      // advective flux is:
      const auto elem_coeff = vdotn * rho_alpha_face / eps_elem_face;
      const auto neighbor_coeff = vdotn * rho_alpha_face / eps_neighbor_face;
      _ae = elem_coeff * d_var_elem_face_d_elem_dof;
      _elem_residual = elem_coeff * var_elem_face;
      _an = -neighbor_coeff * d_var_neighbor_face_d_neighbor_dof;
      _neighbor_residual = -neighbor_coeff * var_neighbor_face;
      if (_approximate_as)
      {
        _ae = _cs / fi.elem().n_sides() * rho_alpha_face / eps_elem_face;
        _an = _cs / fi.neighborPtr()->n_sides() * rho_alpha_face / eps_elem_face;
      }
    }
    else
    {
      const auto [interp_coeffs, advected] =
          interpCoeffsAndAdvected(*_rho_u, advected_face_arg, state);

      const auto elem_arg = elemArg();
      const auto neighbor_arg = neighborArg();

      const auto rho_alpha_elem = _rho(elem_arg, state) * _alpha(elem_arg, state);
      const auto rho_alpha_neighbor = _rho(neighbor_arg, state) * _alpha(elem_arg, state);
      const auto eps_elem = epsilon()(elem_arg, state),
                 eps_neighbor = epsilon()(neighbor_arg, state);
      const auto var_elem = advected.first / rho_alpha_elem * eps_elem,
                 var_neighbor = advected.second / rho_alpha_neighbor * eps_neighbor;

      _ae = vdotn * rho_alpha_elem / eps_elem * interp_coeffs.first;
      // Minus sign because we apply a minus sign to the residual in computeResidual
      _an = -vdotn * rho_alpha_neighbor / eps_neighbor * interp_coeffs.second;

      _elem_residual = _ae * var_elem - _an * var_neighbor;
      _neighbor_residual = -_elem_residual;

      if (_approximate_as)
      {
        _ae = _cs / fi.elem().n_sides() * rho_alpha_elem / eps_elem;
        _an = _cs / fi.neighborPtr()->n_sides() * rho_alpha_neighbor / eps_neighbor;
      }
    }
  }

  _ae *= fi.faceArea() * fi.faceCoord();
  _an *= fi.faceArea() * fi.faceCoord();
  _elem_residual *= fi.faceArea() * fi.faceCoord();
  _neighbor_residual *= fi.faceArea() * fi.faceCoord();
}
