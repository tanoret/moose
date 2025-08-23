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
 * Class describing a convective heat transfer between two domains.
 * The heat flux is described  by:
 * h * (T_solid - T_liquid),
 * where h is the heat transfer coefficient, while T_solid and T_liquid
 * denote the solid and liquid temperatures, respectively.
 */
class LinearFVCoupledMassHeatTransferBC : public LinearFVAdvectionDiffusionBC
{
public:
  static InputParameters validParams();

  /**
   * Class constructor.
   * @param parameters The InputParameters for the object
   */
  LinearFVCoupledMassHeatTransferBC(const InputParameters & parameters);

  virtual Real computeBoundaryValue() const override;

  virtual Real computeBoundaryNormalGradient() const override;

  virtual Real computeBoundaryValueMatrixContribution() const override;

  virtual Real computeBoundaryValueRHSContribution() const override;

  virtual Real computeBoundaryGradientMatrixContribution() const override;

  virtual Real computeBoundaryGradientRHSContribution() const override;

  virtual bool includesMaterialPropertyMultiplier() const override { return true; }

protected:

  /// Function that computes the electro-chemistry corrected equilibrium concentration
  Real computeCorrectedEquilibriumConcentration() const;

  /// Function that computes the mass transfer coefficient
  Real computeMassTransferCoefficient() const;

  /// Function that computes the electro-chemistry corrected equilibrium concentration
  Real computeEffectiveMassTransferCoefficient() const;

  /// The fluid temperature, we use the functor form to enable situations when
  /// the user wants to supply a solution-independent form for this.
  const Moose::Functor<Real> & _c_fluid;

  /// The solid/wall temperature, we use the functor form to enable situations when
  /// the user wants to supply a solution-independent form for this.
  const Moose::Functor<Real> & _c_solid;

  /// The solid/wall temperature, we use the functor form to enable situations when
  /// the user wants to supply a solution-independent form for this.
  const Moose::Functor<Real> & _c_eq;

  /// Helper boolean to see if the variable we have is the fluid variable
  bool _var_is_fluid;

  /// Electrochemistry related variables
  const Real & _alpha_anode;
  const Real & _alpha_cathode;
  const Moose::Functor<Real> & _reaction_potential;
  const Moose::Functor<Real> & _phi; // Background potential (V)
  const Real & _n; // Electrons involved in the reaction
  const Moose::Functor<Real> & _T; // Temperature (K)
  const Real & _i0; // Exchange current density (A)

  /// Parameters to define the mass transfer coefficient via modeling
  const MooseEnum _mass_transfer_treatment;
  const Moose::Functor<Real> * _km; /// The convective heat transfer coefficient
  const Moose::Functor<Real> * _Re; // Reynolds number
  const Moose::Functor<Real> * _Sc; // Schmidt number
  const Moose::Functor<Real> * _D; // Liquid diffusivity coefficient - m2/s 
  const Moose::Functor<Real> * _dh; // Hydraulics diameter - m
  const Moose::Functor<Real> * _k; // Turbulent kintic energy - m2/s2
  const Moose::Functor<Real> * _u_bulk; // Bulk velocity - m/s

  /// The temperature which will contribute to the right hand side.
  /// When this is the fluid domain, the solid temperature will go to the
  /// right hand side. When it is the solid domain, the fluid temperature
  /// will contribute to the right hand side.
  // const Moose::Functor<Real> * _rhs_temperature;

private:

  static constexpr Real R = 8.314; // Universal Gas Constant - J/(mol.K)
  static constexpr Real F = 96485.3321; // Faraday Constant - C/mol
  static constexpr Real C_mu = 0.09; // Kolmogorov Turbulence Constant
};
