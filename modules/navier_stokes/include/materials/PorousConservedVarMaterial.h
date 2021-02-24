//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

class PorousConservedVarMaterial : public Material
{
public:
  PorousConservedVarMaterial(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual void computeQpProperties() override;

  /// fluid properties
  const SinglePhaseFluidProperties & _fluid;

  const ADVariableValue & _var_rho;
  const ADVariableValue & _var_rho_ud;
  const ADVariableValue & _var_rho_vd;
  const ADVariableValue & _var_rho_wd;
  const ADVariableValue & _var_total_energy_density;
  const MaterialProperty<Real> & _epsilon;
  ADMaterialProperty<Real> & _rho;
  ADMaterialProperty<RealVectorValue> & _mass_flux;
  ADMaterialProperty<Real> & _total_energy_density;
  ADMaterialProperty<RealVectorValue> & _velocity;
  ADMaterialProperty<Real> & _vel_x;
  ADMaterialProperty<Real> & _vel_y;
  ADMaterialProperty<Real> & _vel_z;
  ADMaterialProperty<Real> & _v;
  ADMaterialProperty<Real> & _specific_internal_energy;
  ADMaterialProperty<Real> & _pressure;
  ADMaterialProperty<Real> & _specific_total_enthalpy;
  ADMaterialProperty<Real> & _rho_ht;
};