//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVSP3TemperatureBC.h"
#include "Function.h"

registerMooseObject("HeatTransferApp", FVSP3TemperatureBC);

InputParameters
FVSP3TemperatureBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params.addClassDescription("Semi-Transparent Boundary Condition for SP3 Radiation Temperature");

  params.addRequiredParam<MooseFunctorName>("Tb", "The temperature of the boundary");
  params.addParam<MooseFunctorName>("n1", 1.0, "The refraction coefficient of medium 1.");
  params.addParam<MooseFunctorName>("n2", 1.0, "The refraction coefficient of medium 2.");
  params.addRequiredParam<MooseFunctorName>("h", "The convective heat transfer coefficient of medium.");
  params.addRequiredParam<MooseFunctorName>("k", "The conductive heat transfer coefficient of medium.");

  params.addRequiredParam<MooseFunctorName>("alpha", "The hemispheric emissivity of the medium.");
  params.addRequiredParam<Real>("nu1", "The maximum opaque frequency of the medium.");
  params.addParam<unsigned int> ("Nint", "The number of Plank function integral interval");

  return params;
}

FVSP3TemperatureBC::FVSP3TemperatureBC(const InputParameters & parameters)
  : FVFluxBC(parameters),
    _Tb(getFunctor<ADReal>("Tb")),
    _n1(getFunctor<ADReal>("n1")),
    _n2(getFunctor<ADReal>("n2")),
    _h(getFunctor<ADReal>("h")),
    _k(getFunctor<ADReal>("k")),
    _alpha(getFunctor<ADReal>("alpha")),
    _nu1(getParam<Real>("nu1")),
    _Nint(getParam<unsigned int>("Nint"))
{
}

ADReal
FVSP3TemperatureBC::computeQpResidual()
{
  // Allow the functors to pick their side evaluation
  const Moose::FaceArg face{
      _face_info, Moose::FV::LimiterType::CentralDifference, true, false, nullptr, nullptr};
  const auto state = determineState();

  // Build the convective source at the boundary
  const auto T = _var(face, state);
  const auto Tb = _Tb(face, state);
  const auto thermal_conv_source = _h(face, state) * (Tb - T);

  // Build the radiative source at the boundary
  ADReal Plank_integral(0.0);
  const auto d_nu = _nu1 / _Nint;
  const auto n1_pow_2 = Utility::pow<2>(_n1(face, state));

  const auto prefactor_coeff = n1_pow_2 * 2 * HeatConduction::Constants::hp  / (Utility::pow<2>(HeatConduction::Constants::c));
  const auto source_coeff = HeatConduction::Constants::hp / (HeatConduction::Constants::kb);

  for (int itr = 0; itr < (int)_Nint; itr++){
    const auto local_nu = (itr + 0.5) * d_nu;
    const auto local_nu_pow_3 = Utility::pow<3>(local_nu);

    const auto local_prefactor = prefactor_coeff * local_nu_pow_3;
    const auto local_boundary_source = std::exp(source_coeff * local_nu / Tb) - 1.;
    const auto local_cell_source = std::exp(source_coeff * local_nu / T) - 1.;

    Plank_integral += local_prefactor * (1 / local_boundary_source - 1 / local_cell_source) * d_nu;
  }

  const auto thermal_rad_source = _alpha(face, state) * libMesh::pi * Utility::pow<2>(_n2(face, state)/_n1(face, state)) * Plank_integral;
  
  return -1. * (thermal_conv_source + thermal_rad_source) / _k(face, state);
}
