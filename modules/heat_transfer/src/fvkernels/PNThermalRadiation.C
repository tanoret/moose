//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PNThermalRadiation.h"
#include "MathUtils.h"
#include "HeatConductionNames.h"

registerMooseObject("HeatTransferApp", PNThermalRadiation);

InputParameters
PNThermalRadiation::validParams()
{
  InputParameters params = FVDiffusion::validParams();
  params.addClassDescription(
      "Kernel that assembles the equations for PN thermal radiation transport.");

  params.addRequiredParam<unsigned int>("n", "N-th Order of the Equation.");

  params.addRequiredParam<MooseFunctorName>(
      "phi_n_minus_2", "Lagged order of the radiation flux.");
  params.addParam<MooseFunctorName>(
      "phi_n_plus_2", 0.0, "Leading order of the radiation flux.");

  params.addRequiredParam<MooseFunctorName>(
        "sigma_n_minus_1", "Lagged order of the radiation flux.");
  params.addParam<MooseFunctorName>(
        "sigma_n_plus_1", 1.0, "Leading order of the radiation flux.");

  return params;
}

PNThermalRadiation::PNThermalRadiation(const InputParameters & params)
  : FVDiffusion(params),
    _n(getParam<unsigned int>("n")),
    _phi_n_minus_2(getFunctor<ADReal>("phi_n_minus_2")),
    _phi_n_plus_2(getFunctor<ADReal>("phi_n_plus_2")),
    _sigma_n_minus_1(getFunctor<ADReal>("sigma_n_minus_1")),
    _sigma_n_plus_1(getFunctor<ADReal>("sigma_n_plus_1"))
{
}

ADReal
PNThermalRadiation::computeQpResidual()
{
  using namespace Moose::FV;
  const auto state = determineState();
  const auto elem_arg = elemArg();

  // Building the residual for the leading order
  const auto n_order_plus_2 = (_n+1)*(_n+2)/((2*_n+1)*(2*_n+3));
  ADReal sigma_n_plus_1_interpolated;
  interpolate(_coeff_interp_method,
    sigma_n_plus_1_interpolated,
    _sigma_n_plus_1(elemArg(), state),
    _sigma_n_plus_1(neighborArg(), state),
    *_face_info,
    true);
  const auto grad_phi_n_plus_2_normal = _phi_n_plus_2.gradient(elem_arg, state) * _face_info->normal();
  const auto res_leading_order = -n_order_plus_2 / sigma_n_plus_1_interpolated * grad_phi_n_plus_2_normal;

  // Building the residual for the leagged order
  const auto n_order_minus_2 = _n*(_n-1)/((2*_n+1)*(2*_n-1));
  ADReal sigma_n_minus_1_interpolated;
  interpolate(_coeff_interp_method,
    sigma_n_minus_1_interpolated,
    _sigma_n_minus_1(elemArg(), state),
    _sigma_n_minus_1(neighborArg(), state),
    *_face_info,
    true);
  const auto grad_phi_n_minus_2_normal = _phi_n_minus_2.gradient(elem_arg, state) * _face_info->normal();
  const auto res_lagged_order = -n_order_minus_2 / sigma_n_minus_1_interpolated * grad_phi_n_minus_2_normal;

  // Building the residual for the current term
  const auto current_n_order_plus_2 = Utility::pow<2>(_n+1)/((2*_n+1)*(2*_n+3));
  const auto current_n_order_minus_2 = Utility::pow<2>(_n)/((2*_n+1)*(2*_n-1));
  const auto res_order = -(current_n_order_plus_2/sigma_n_plus_1_interpolated + 
                           current_n_order_minus_2/sigma_n_minus_1_interpolated) * gradUDotNormal(state);

  return res_leading_order + res_lagged_order + res_order;
}
