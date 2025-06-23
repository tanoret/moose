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
class TorchScriptTurbulentAnisotropyMaterial : public Material
{
public:
  static InputParameters validParams();

  TorchScriptTurbulentAnisotropyMaterial(const InputParameters & parameters);

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
  /// Density
  const Moose::Functor<ADReal> & _rho;
  /// k for viscosity
  const Moose::Functor<ADReal> & _k;
  /// epsilon for viscosity
  const Moose::Functor<ADReal> & _eps;
  /// Debug output flag
  const bool & _debug;
  /// Use NN or ASRM
  const bool & _use_NN;

  /// The user object that holds the torch module
  const TorchScriptUserObject & _torch_script_userobject;

  /// Place holder for the inputs to the neural network
  torch::Tensor _input_tensor;

  /// Vector of all property prefix
  const MooseFunctorName & _property_prefix;
  
  /// Vector of all the properties, for now we don't support AD
  std::vector<GenericMaterialProperty<Real, false>*>  _properties;

private:
  /**
   * A helper method for evaluating the torch script module and populating the
   * material properties.
   */

  void arsm(const Real& eta1, const Real& eta2, Real& G1, Real& G2, Real& G3);
  void computeQpValues();
};

#endif
