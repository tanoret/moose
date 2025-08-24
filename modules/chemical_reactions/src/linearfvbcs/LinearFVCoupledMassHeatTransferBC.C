//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVCoupledMassHeatTransferBC.h"

registerMooseObject("ChemicalReactionsApp", LinearFVCoupledMassHeatTransferBC);

InputParameters
LinearFVCoupledMassHeatTransferBC::validParams()
{
  InputParameters params = LinearFVAdvectionDiffusionBC::validParams();

  params.addRequiredParam<MooseFunctorName>("c_fluid", "Liquid concentration variable [mol/m^3]");
  params.addRequiredParam<MooseFunctorName>("c_solid", "Solid/wall variable");
  params.addRequiredParam<MooseFunctorName>("c_eq",
                                            "Equilibrium liquid concentration at eta=0 [mol/m^3]");

  // Electrochemistry (default: symmetric transfer, no driving)
  params.addParam<Real>("alpha_anode", 0.5, "Anodic transfer coefficient [-]");
  params.addParam<Real>("alpha_cathode", 0.5, "Cathodic transfer coefficient [-]");
  params.addParam<MooseFunctorName>("E_reaction", 0.0, "Reaction potential E_eq [V]");
  params.addParam<MooseFunctorName>("phi", 0.0, "Background electric potential phi [V]");
  params.addParam<Real>("n", 1.0, "Electrons transferred n [-]");
  params.addParam<MooseFunctorName>("T", 300.0, "Temperature [K]");
  params.addParam<Real>("i0", 0.0, "Exchange current density i0 [A/m^2]");

  MooseEnum mass_transfer_treatment("constant correlation resolved", "constant");
  params.addParam<MooseEnum>("mass_transfer_treatment",
                             mass_transfer_treatment,
                             "Mass-transfer model: 'constant', 'correlation', or 'resolved'");

  // Film models
  params.addParam<MooseFunctorName>("km", "Direct mass-transfer coefficient k_m [m/s]");
  params.addParam<MooseFunctorName>("Re", "Reynolds number [-]");
  params.addParam<MooseFunctorName>("Sc", "Schmidt number [-]");
  params.addParam<MooseFunctorName>("D", "Molecular diffusivity [m^2/s]");
  params.addParam<MooseFunctorName>("dh", "Hydraulic diameter [m]");
  params.addParam<MooseFunctorName>("k", "Turbulent kinetic energy k [m^2/s^2]");
  params.addParam<MooseFunctorName>("vel_bulk", "Bulk velocity [m/s]");

  params.addClassDescription(
      "Robin mass-transfer BC with electrochemical correction (plating/corrosion).");
  return params;
}

LinearFVCoupledMassHeatTransferBC::LinearFVCoupledMassHeatTransferBC(
    const InputParameters & parameters)
  : LinearFVAdvectionDiffusionBC(parameters),
    _c_fluid(getFunctor<Real>("c_fluid")),
    _c_solid(getFunctor<Real>("c_solid")),
    _c_eq(getFunctor<Real>("c_eq")),
    _var_is_fluid("wraps_" + _var.name() == _c_fluid.functorName() ||
                  "wraps_" + _var.name() + "_raw_value" == _c_fluid.functorName()),
    _alpha_anode(getParam<Real>("alpha_anode")),
    _alpha_cathode(getParam<Real>("alpha_cathode")),
    _reaction_potential(getFunctor<Real>("E_reaction")),
    _phi(getFunctor<Real>("phi")),
    _n(getParam<Real>("n")),
    _T(getFunctor<Real>("T")),
    _i0(getParam<Real>("i0")),
    _mass_transfer_treatment(getParam<MooseEnum>("mass_transfer_treatment")),
    _km(parameters.isParamValid("km") ? &(getFunctor<Real>("km")) : nullptr),
    _Re(parameters.isParamValid("Re") ? &(getFunctor<Real>("Re")) : nullptr),
    _Sc(parameters.isParamValid("Sc") ? &(getFunctor<Real>("Sc")) : nullptr),
    _D(parameters.isParamValid("D") ? &(getFunctor<Real>("D")) : nullptr),
    _dh(parameters.isParamValid("dh") ? &(getFunctor<Real>("dh")) : nullptr),
    _k(parameters.isParamValid("k") ? &(getFunctor<Real>("k")) : nullptr),
    _u_bulk(parameters.isParamValid("vel_bulk") ? &(getFunctor<Real>("vel_bulk")) : nullptr)
{
  // Parameter checks for k_m models
  if (_mass_transfer_treatment == "constant" && !_km)
    paramError("km", "Provide 'km' for constant mass-transfer treatment.");
  else if (_mass_transfer_treatment == "correlation")
  {
    if (!_Re)
      paramError("Re", "Provide 'Re' for correlation treatment.");
    if (!_Sc)
      paramError("Sc", "Provide 'Sc' for correlation treatment.");
    if (!_D)
      paramError("D", "Provide 'D' for correlation treatment.");
    if (!_dh)
      paramError("dh", "Provide 'dh' for correlation treatment.");
  }
  else if (_mass_transfer_treatment == "resolved")
  {
    if (!_k)
      paramError("k", "Provide 'k' for resolved treatment.");
    if (!_u_bulk)
      paramError("vel_bulk", "Provide 'vel_bulk' for resolved treatment.");
    if (!_Re)
      paramError("Re", "Provide 'Re' for resolved treatment.");
    if (!_Sc)
      paramError("Sc", "Provide 'Sc' for resolved treatment.");
    if (!_D)
      paramError("D", "Provide 'D' for resolved treatment.");
    if (!_dh)
      paramError("dh", "Provide 'dh' for resolved treatment.");
  }
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryValue() const
{
  const auto elem_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                             ? _current_face_info->elemInfo()
                             : _current_face_info->neighborInfo();
  return _var.getElemValue(*elem_info, determineState());
}

Real
LinearFVCoupledMassHeatTransferBC::computeCorrectedEquilibriumConcentration() const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  // Overpotential: eta = phi - E_reaction
  const Real eta = _phi(face, state) - _reaction_potential(face, state);
  const Real Tval = _T(face, state);

  // c_eq*(eta) = c_eq * exp(+alpha_a F eta / RT)  (clamped for stability)
  const Real expo = std::max(std::min(_alpha_anode * F * _n * eta / (R * Tval), 10.0), -10.0);
  return _c_eq(face, state) * std::exp(expo);
}

Real
LinearFVCoupledMassHeatTransferBC::computeMassTransferCoefficient() const
{
  Real km_val = 0.0;
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  if (_mass_transfer_treatment == "constant")
    km_val = (*_km)(face, state);

  else if (_mass_transfer_treatment == "correlation")
  {
    const Real Re = (*_Re)(face, state);
    const Real Sc = (*_Sc)(face, state);

    Real Sh_lam = 1.86 * std::pow(Re * Sc, 1.0 / 3.0);             // developing laminar
    Real Sh_tur = 0.023 * std::pow(Re, 0.83) * std::pow(Sc, 0.33); // fully turbulent (fixed)

    Real Sh = 0.0;
    if (Re <= 3e3)
      Sh = Sh_lam;
    else if (Re >= 1e4)
      Sh = Sh_tur;
    else
    {
      const Real w = (Re - 3e3) / (1e4 - 3e3); // smooth blend
      Sh = (1.0 - w) * Sh_lam + w * Sh_tur;
    }

    km_val = Sh * (*_D)(face, state) / std::max((*_dh)(face, state), 1e-10);
  }

  else // "resolved" via Chilton–Colburn using wall friction from k
  {
    // u_*^2 ~ C_mu^{1/2} k;   f = 2 u_*^2 / u_b^2 = 2 sqrt(C_mu) k / u_b^2
    const Real f =
        2.0 * std::sqrt(C_mu) * (*_k)(face, state) / Utility::pow<2>((*_u_bulk)(face, state));

    const Real Sh = 0.5 * f * (*_Re)(face, state) * std::pow((*_Sc)(face, state), 1.0 / 3.0);
    km_val = Sh * (*_D)(face, state) / std::max((*_dh)(face, state), 1e-10);
  }

  return km_val;
}

Real
LinearFVCoupledMassHeatTransferBC::computeEffectiveMassTransferCoefficient() const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  // k0 from i0: k0 = i0 / (n F c_eq)  (units m/s)   *** FIXED UNITS ***
  const Real c_eq0 = std::max(_c_eq(face, state), 1e-30);
  const Real k0 = (_i0 / std::max(_n * F * c_eq0, 1e-30));

  // Overpotential eta = phi - E ; reaction-limited exponent uses -alpha_c
  const Real eta = _phi(face, state) - _reaction_potential(face, state);
  const Real Tval = _T(face, state);
  const Real expo = std::max(std::min(-_alpha_cathode * F * _n * eta / (R * Tval), 10.0), -10.0);

  const Real k0_eff = k0 * std::exp(expo);
  const Real km = this->computeMassTransferCoefficient();

  // Series resistance
  return km * k0_eff / std::max(km + k0_eff, 1e-12);
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryNormalGradient() const
{
  // Provide the physical molar flux J as a normal gradient equivalent for FV:
  // J = k_eff * (c_fluid - c_eq*)
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  const Real keff = this->computeEffectiveMassTransferCoefficient();
  const Real ceqcor = this->computeCorrectedEquilibriumConcentration();
  const Real cL = _c_fluid(face, state);

  const Real J = keff * (cL - ceqcor); // [mol m^-2 s^-1], positive = plating (liq->solid)

  // Orient with outward normal of *this* variable's side
  const auto elem_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                             ? _current_face_info->elemInfo()
                             : _current_face_info->neighborInfo();
  const auto neighbor_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                                 ? _current_face_info->neighborInfo()
                                 : _current_face_info->elemInfo();

  const auto this_side = _var_is_fluid ? elem_info : neighbor_info;

  const int sgn = (_current_face_info->normal() *
                   (_current_face_info->faceCentroid() - this_side->centroid())) > 0
                      ? 1
                      : -1;

  return sgn * J;
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryValueMatrixContribution() const
{
  // First-order FV face value approximation -> coefficient on c at the face
  return 1.0;
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryValueRHSContribution() const
{
  // No explicit RHS from the value part in this Robin implementation
  return 0.0;
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryGradientMatrixContribution() const
{
  // For the liquid equation, assemble the k_eff coefficient on c_liq
  if (_var_is_fluid)
    return this->computeEffectiveMassTransferCoefficient();
  else
    return 0.0; // solid-side variable is fed via RHS only (acts as accumulator/source)
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryGradientRHSContribution() const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  const Real keff = this->computeEffectiveMassTransferCoefficient();
  const Real ceqcor = this->computeCorrectedEquilibriumConcentration();

  if (_var_is_fluid)
  {
    // Fluid equation: -n·(D∇c) = k_eff (c - c_eq*) -> RHS = k_eff * c_eq*
    return keff * ceqcor;
  }
  else
  {
    // Solid-side: add the same interfacial molar flux J to the solid variable
    // J = k_eff * (c_liq - c_eq*)
    // Note: this makes plating (J>0) increase solid content; corrosion (J<0) decreases it.
    const Real cL = _c_fluid(face, state);
    return keff * (cL - ceqcor);
  }
}
