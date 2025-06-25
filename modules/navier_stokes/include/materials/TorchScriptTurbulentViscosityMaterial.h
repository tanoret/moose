//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef LIBTORCH_ENABLED

#pragma once

#include "Material.h"
#include "TorchScriptUserObject.h"

/**
 * This material declares properties which are evaluated as
 * based on a torch script neural network.
 */
class TorchScriptTurbulentViscosityMaterial : public Material
{
public:
  static InputParameters validParams();

  TorchScriptTurbulentViscosityMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Dimmension of the mesh
  const unsigned int _mesh_dimension;

  /// x-velocity
  const Moose::Functor<ADReal> & _u_var;
  /// y-velocity
  const Moose::Functor<ADReal> * _v_var;
  /// z-velocity
  const Moose::Functor<ADReal> * _w_var;
  /// k for viscosity
  const Moose::Functor<ADReal> & _k;
  /// epsilon for viscosity
  const Moose::Functor<ADReal> & _eps;
  /// Debug output flag
  const bool & _debug;
  /// Use NN or ASRM
  const bool & _use_NN;
  /// Minimum output for nu_t
  const Real & _mu_t_min;
  /// Max output for nu_t
  const Real & _mu_t_max;
  /// The user object that holds the torch module
  const TorchScriptUserObject & _torch_script_userobject;

  /// Old prop functor
  const Moose::Functor<ADReal> * _mu_t_old;

  /// Relaxation Factor
  const Real _rf;

  /// Place holder for the inputs to the neural network
  torch::Tensor _input_tensor;

  /// Vector of all the properties, for now we don't support AD
  GenericMaterialProperty<Real, false> * _properties;


private:
  /**
   * A helper method for evaluating the torch script module and populating the
   * material properties.
   */

  void arsm(const Real& eta1, const Real& eta2, Real& G1, Real& G2, Real& G3);
  void computeQpValues();
};

#endif
