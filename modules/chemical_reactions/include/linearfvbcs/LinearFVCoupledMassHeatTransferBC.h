//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVAdvectionDiffusionBC.h"

/**
 * Mass-transfer (not heat) Robin boundary condition coupling a liquid-side
 * species to a solid/wall-side variable with optional electrochemical
 * correction via Butler–Volmer linearization.
 *
 * Flux (molar, from liquid to solid) is:
 *   J = k_eff(eta) * ( c_liq - c_eq*(eta) )   [mol m^-2 s^-1]
 * with k_eff = (k_m k0*) / (k_m + k0*),  k0* = k0 exp(-alpha_c F eta / RT),
 *      c_eq*(eta) = c_eq exp(+alpha_a F eta / RT), and  eta = phi - E_reaction.
 *
 * This BC is sign-symmetric: if c_liq < c_eq*(eta) it yields corrosion
 * (J<0, solid-to-liquid), otherwise plating (J>0).
 */
class LinearFVCoupledMassHeatTransferBC : public LinearFVAdvectionDiffusionBC
{
public:
  static InputParameters validParams();
  LinearFVCoupledMassHeatTransferBC(const InputParameters & parameters);

  // FV hooks
  Real computeBoundaryValue() const override;
  Real computeBoundaryNormalGradient() const override;
  Real computeBoundaryValueMatrixContribution() const override;
  Real computeBoundaryValueRHSContribution() const override;
  Real computeBoundaryGradientMatrixContribution() const override;
  Real computeBoundaryGradientRHSContribution() const override;

  bool includesMaterialPropertyMultiplier() const override { return true; }

protected:
  // Electrochemically corrected equilibrium concentration c_eq*(eta)
  Real computeCorrectedEquilibriumConcentration() const;

  // Film mass-transfer coefficient k_m
  Real computeMassTransferCoefficient() const;

  // Effective (series) mass-transfer coefficient k_eff(eta)
  Real computeEffectiveMassTransferCoefficient() const;

  // ---- Inputs / Functors ----
  const Moose::Functor<Real> & _c_fluid; // liquid concentration [mol/m^3]
  const Moose::Functor<Real> & _c_solid; // solid-side variable (accumulator or conc)
  const Moose::Functor<Real> & _c_eq;    // base equilibrium conc (eta=0) [mol/m^3]

  bool _var_is_fluid;

  // Electrochemistry
  const Real & _alpha_anode;                        // [-]
  const Real & _alpha_cathode;                      // [-]
  const Moose::Functor<Real> & _reaction_potential; // E_reaction [V]
  const Moose::Functor<Real> & _phi;                // background potential [V]
  const Real & _n;                                  // electrons transferred [-]
  const Moose::Functor<Real> & _T;                  // temperature [K]
  const Real & _i0;                                 // exchange current density [A/m^2]

  // Mass-transfer modeling
  const MooseEnum _mass_transfer_treatment;
  const Moose::Functor<Real> * _km;     // direct k_m [m/s]
  const Moose::Functor<Real> * _Re;     // Reynolds [-]
  const Moose::Functor<Real> * _Sc;     // Schmidt [-]
  const Moose::Functor<Real> * _D;      // molecular diffusivity [m^2/s]
  const Moose::Functor<Real> * _dh;     // hydraulic diameter [m]
  const Moose::Functor<Real> * _k;      // turbulent kinetic energy [m^2/s^2]
  const Moose::Functor<Real> * _u_bulk; // bulk velocity [m/s]

private:
  static constexpr Real R = 8.314;      // J/(mol K)
  static constexpr Real F = 96485.3321; // C/mol
  static constexpr Real C_mu = 0.09;    // k–ε constant
};
