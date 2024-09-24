//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVInterfaceKernel.h"

class FVButlerBolmerInterface : public FVInterfaceKernel
{
public:
  static InputParameters validParams();
  FVButlerBolmerInterface(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  /// Model constants
  const Real _F;
  const Real _z;
  const Real _R;

  /// Temperature field
  const Moose::Functor<ADReal> & _T;

  /// Electric potential
  const Moose::Functor<ADReal> & _phi;

  /// The main tracked variable
  const Moose::Functor<ADReal> & _c;

  /// The main tracked variable on the solid side
  const Moose::Functor<ADReal> * _c_solid;

  /// Standard electrode potential
  const Real _c_E;

  /// Standard electrode potential temperature derivative
  const Real _c_dE_dT;

  /// Reference electrode temperature
  const Real _T_ref;

  /// Electrode reaction charge
  const unsigned int _c_Z;

  /// Formation energy
  const Real _c_G;

  /// Formulation used for electro-kinetics
  const MooseEnum & _interface_formulation;

  /// Reaction kinetics constant
  const Real * _c_k0;

  /// Partial current for single partial current
  const Real _c_alpha;

  /// Coupled species
  std::vector<MooseFunctorName> _coupled_c;

  /// Coupled species on the solid side
  const std::vector<MooseFunctorName> * _coupled_c_solid;

  /// Coupled species electrode potential
  std::vector<Real> _coupled_E;

  /// Standard electrode potential temperature derivative
  const std::vector<Real> * _coupled_dE_dT;

  /// Coupled species electrode reaction charge
  std::vector<unsigned int> _coupled_Z;

  /// Coupled species formation energy
  const std::vector<Real> * _coupled_G;

  /// Reactions constant
  const std::vector<Real> * _coupled_k0;

  /// Partial current for reaction
  const std::vector<Real> * _coupled_alpha;

  /// Vector of the functors
  std::vector<const Moose::Functor<GenericReal<true>> *> _coupled_c_functors;
  std::vector<const Moose::Functor<GenericReal<true>> *> _coupled_c_solid_functors;

  /// The distance from the wall before evaluating the bulk temperature
  const Real _bulk_distance;

  /// Whether to use the wall cell for the bulk fluid temperature
  const bool _use_wall_cell;

  /// libmesh object to find points in the mesh
  std::unique_ptr<PointLocatorBase> _pl;

  /// Boolean to see if variable1 is the fluid
  const bool _var1_is_fluid;
};
