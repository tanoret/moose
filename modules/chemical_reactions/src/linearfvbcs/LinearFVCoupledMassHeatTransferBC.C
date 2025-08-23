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
  params.addRequiredParam<MooseFunctorName>("c_fluid", "The fluid conentration variable");
  params.addRequiredParam<MooseFunctorName>("c_solid", "The solid/wall concentration variable");
  params.addRequiredParam<MooseFunctorName>("c_eq", "The fluid equilibrium concentration variable");

  params.addParam<Real>("alpha_anode", 0.5, "Anodic charge transfer coefficient");
  params.addParam<Real>("alpha_cathode", 0.5, "Cathodic charge transfer coefficient");
  params.addParam<MooseFunctorName>("E_reaction", 0.0, "Reaction potential");
  params.addParam<MooseFunctorName>("phi", 0.0, "Brackground electric potential");
  params.addParam<Real>("n", 1.0, "Number of transferred electrons in the reaction");
  params.addParam<MooseFunctorName>("T", 300.0, "Temperature");
  params.addParam<Real>("i0", 0.0, "Exchange current density.");

  MooseEnum mass_transfer_treatment("constant correlation resolved", "constant");

  params.addParam<MooseEnum>("mass_transfer_treatment",
                             mass_transfer_treatment,
                             "The method used for computing the mass transfer coefficient "
                             "'constant', 'correlation', 'resolved'");
  params.addParam<MooseFunctorName>("km", "The convective heat transfer coefficient");
  params.addParam<MooseFunctorName>("Re", "The Reynolds number");
  params.addParam<MooseFunctorName>("Sc", "The Schmidt number");
  params.addParam<MooseFunctorName>("D", "The molecular diffusion coefficient");
  params.addParam<MooseFunctorName>("dh", "The hydraulics diameter");
  params.addParam<MooseFunctorName>("k", "The turbulent kinetic energy");
  params.addParam<MooseFunctorName>("vel_bulk", "The bulk velocity");

  params.addClassDescription("Class describing a convective heat transfer between two domains.");
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

  // Check for parameter errors in the mass transport models
  if(_mass_transfer_treatment == "constant" && !_km)
    paramError("km", "The mass transfer coefficient should be specified for the constant mass transfer treatment.");
  else if (_mass_transfer_treatment == "correlation")
  {
    if(!_Re)
      paramError("Re", "The Reynolds number should be specified for the correlation treatment of the mass transfer coefficient.");

    if(!_Sc)
      paramError("Sc", "The Schmidt number should be specified for the correlation treatment of the mass transfer coefficient.");

    if(!_D)
      paramError("D", "The molecular diffusion coefficient should be specified for the correlation treatment of the mass transfer coefficient.");

    if(!_dh)
      paramError("dh", "The hydraulic diameter should be specified for the correlation treatment of the mass transfer coefficient.");
  }
  else if (_mass_transfer_treatment == "resolved")
  {
    if(!_k)
      paramError("k", "The turbulent kinetic energy should be specified for the correlation treatment of the mass transfer coefficient.");

    if(!_u_bulk)
      paramError("vel_bulk", "The bulk velocity should be specified for the correlation treatment of the mass transfer coefficient.");

    if(!_D)
      paramError("D", "The molecular diffusion coefficient should be specified for the correlation treatment of the mass transfer coefficient.");

    if(!_dh)
      paramError("dh", "The hydraulic diameter should be specified for the correlation treatment of the mass transfer coefficient.");
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

  const Real E = _reaction_potential(face, state);
  const Real phi = _phi(face, state);
  const Real T = _T(face, state);

  const Real exponent = std::max(std::min(_alpha_anode * F * _n * (E - phi) / (R * T), 10.0), -10.0);

  return _c_eq(face, state) * std::exp(exponent);
  
}

Real
LinearFVCoupledMassHeatTransferBC::computeMassTransferCoefficient() const
{
  Real mtc;
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  if(_mass_transfer_treatment == "constant")
    mtc = (*_km)(face, state);
  else if(_mass_transfer_treatment == "correlation")
  {
    const auto Re = (*_Re)(face, state);
    const auto Sc = (*_Sc)(face, state);
    Real Sh; // Sherwood number
    if(Re <= 3e3)
      Sh = 1.86 * std::pow(Re * Sc, 1./3.);
    else if (Re >= 1e4)
      Sh = 1.86 * std::pow(Re, 0.83) * std::pow(Sc, 0.33);
    else
    {
      const auto Sc_lam = 1.86 * std::pow(Re * Sc, 1./3.);
      const auto Sc_tur = 1.86 * std::pow(Re, 0.83) * std::pow(Sc, 0.33);
      const auto w_tur = (Re - 1e3) / (1e4 - 1e3);
      Sh = (1. - w_tur) * Sc_lam + w_tur * Sc_tur;
    }
    mtc = Sh * (*_D)(face, state) / std::max((*_dh)(face, state), 1e-10);
  }
  else // _mass_transfer_treatment == "resolved"
  {
    const auto f = 2.0 * std::sqrt(C_mu) * (*_k)(face, state) / Utility::pow<2>((*_u_bulk)(face, state));
    const auto Sh = f/2.0 * (*_Re)(face, state) * std::pow((*_Sc)(face, state), 1./3.);
    mtc = Sh * (*_D)(face, state) / std::max((*_dh)(face, state), 1e-10);
  }
  return mtc;
}

Real
LinearFVCoupledMassHeatTransferBC::computeEffectiveMassTransferCoefficient() const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  const Real k0 = _i0 / (_n * F);

  const Real E = _reaction_potential(face, state);
  const Real phi = _phi(face, state);
  const Real T = _T(face, state);

  const Real exponent = std::max(std::min(-_alpha_cathode * F * _n * (E - phi) / (R * T), 10.0), -10.0);
  const Real k0_eff = k0 * std::exp(exponent);
  const Real km = this->computeMassTransferCoefficient();

  return km*k0_eff / std::max(km + k0_eff, 1e-10);
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryNormalGradient() const
{
  const auto face = singleSidedFaceArg(_current_face_info);
  const auto state = determineState();

  const auto elem_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                             ? _current_face_info->elemInfo()
                             : _current_face_info->neighborInfo();

  const auto neighbor_info = (_current_face_type == FaceInfo::VarFaceNeighbors::ELEM)
                                 ? _current_face_info->neighborInfo()
                                 : _current_face_info->elemInfo();

  const auto fluid_side_elem_info = _var_is_fluid ? elem_info : neighbor_info;

  // All this fuss is just for cases when we have an internal boundary, then the flux will change
  // signs depending on which side of the face we are at.
  const auto multiplier = _current_face_info->normal() * (_current_face_info->faceCentroid() -
                                                          fluid_side_elem_info->centroid()) >
                                  0
                              ? 1
                              : -1;

  return multiplier * this->computeEffectiveMassTransferCoefficient() * (_c_fluid(face, state) - this->computeCorrectedEquilibriumConcentration());
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryValueMatrixContribution() const
{
  // We approximate the face value with the cell value here.
  // TODO: we can extend this to a 2-term expansion at some point when the need arises.
  return 1.0;
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryValueRHSContribution() const
{
  // We approximate the face value with the cell value, we
  // don't need to add anything to the right hand side.
  // TODO: we can extend this to a 2-term expansion at some point when the need arises.
  return 0.0;
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryGradientMatrixContribution() const
{
  // We just put the heat transfer coefficient on the diagonal (multiplication with the
  // surface area is taken care of in the kernel).
  if (_var_is_fluid)
    return this->computeEffectiveMassTransferCoefficient();
  else
    return 0.0;
}

Real
LinearFVCoupledMassHeatTransferBC::computeBoundaryGradientRHSContribution() const
{
  // We check where the functor contributing to the right hand side lives. We do this
  // because this functor lives on the domain where the variable of this kernel doesn't.
  if (_var_is_fluid)
    return this->computeEffectiveMassTransferCoefficient() * this->computeCorrectedEquilibriumConcentration();
  else
  {
    auto face_cl = singleSidedFaceArg(_current_face_info);
    const auto state = determineState();

    if (_c_fluid.hasFaceSide(*_current_face_info, true))
      face_cl.face_side = _current_face_info->elemPtr();
    else
      face_cl.face_side = _current_face_info->neighborPtr();
    return this->computeEffectiveMassTransferCoefficient() * (_c_fluid(face_cl, state) - this->computeCorrectedEquilibriumConcentration());
  }
}
