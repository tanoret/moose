//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVSP3ThermalRadiationBC.h"
#include "Function.h"

registerMooseObject("HeatTransferApp", FVSP3ThermalRadiationBC);

InputParameters
FVSP3ThermalRadiationBC::validParams()
{
  InputParameters params = FVFluxBC::validParams();
  params.addClassDescription("Marshak Boundary Condition for SP3 Radiation Transport");

  params.addRequiredRangeCheckedParam<MooseFunctorName>("T", "T>0", "The temperature of the medium.");
  params.addRequiredRangeCheckedParam<MooseFunctorName>("nu", "T>0", "The mean frequency of the thermal radiation band.");
  params.addParam<MooseFunctorName>("refraction_index", 1.0, "The refraction index in the spectral band.");
  params.addRequiredParam<MooseFunctorName>("kappa", "The absorptivity of the medium.");

  params.addRequiredParam<MooseFunctorName>("epsilon", "The optical thickness of the medium.");
  params.addRequiredParam<MooseFunctorName>("psi", "The incoming momentum flux from the conjugated SP3 order.");

  MooseEnum order("first second", "first");
  params.addParam<MooseEnum>("order", order, "The order of the diffusion term.");

  return params;
}

FVSP3ThermalRadiationBC::FVSP3ThermalRadiationBC(const InputParameters & parameters)
  : FVFluxBC(parameters),
    _T(getFunctor<ADReal>("T")),
    _nu(getFunctor<ADReal>("nu")),
    _n1(getFunctor<ADReal>("refraction_index")),
    _absorptivity(getFunctor<ADReal>("kappa")),
    _optical_thickness(getFunctor<ADReal>("epsilon")),
    _psi(getFunctor<ADReal>("psi")),
    _order(getParam<MooseEnum>("order"))
{
  if (_order == "first")
  {
    _alpha_order = _alpha_1;
    _beta_order = _beta_2;
    _eta_order = _eta_1;
  }
  else
  {
    _alpha_order = _alpha_2;
    _beta_order = _beta_1;
    _eta_order = _eta_2;
  }
}

ADReal
FVSP3ThermalRadiationBC::computeQpResidual()
{
  // Negative means incoming flux

  // Allow the functors to pick their side evaluation
  const Moose::FaceArg face{
      _face_info, Moose::FV::LimiterType::CentralDifference, true, false, nullptr, nullptr};
  const auto state = determineState();

  // Build the emission source at the boundary
  const auto n1_pow_2 = Utility::pow<2>(_n1(face, state));
  const auto nu = _nu(face, state);
  const auto nu_pow_3 = Utility::pow<3>(nu);
  const auto T = _T(face, state);

  const auto pre_factor = n1_pow_2 * 2.0 * HeatConduction::Constants::hp * nu_pow_3 / (Utility::pow<2>(HeatConduction::Constants::c));
  const auto inv_thermal_source = std::exp(HeatConduction::Constants::hp*nu/(HeatConduction::Constants::kb * T)) - 1.0;
  const auto thermal_rad_source = -_eta_order * pre_factor / inv_thermal_source;

  // Radiation leaking the boundary
  const auto thermal_rad_sink = _beta_order * _psi(face, state) + _alpha_order * _var(face, state);
  
  // Compute Marshak flux
  const auto protected_epsilon = (_optical_thickness(face, state) > 1e-12 ? _optical_thickness(face, state) : 1e-12);
  const auto flux = (thermal_rad_source + thermal_rad_sink) * _absorptivity(face, state) / protected_epsilon;
  return flux;
}
