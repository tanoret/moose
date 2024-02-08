//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVFunctionDirichletBC.h"
#include "INSFVFullyDevelopedFlowBC.h"

/**
 * A class for setting the value of the pressure at an outlet of the system.
 * It may not be used with a mean/pinned-pressure approach
 */
template <class T>
class INSFVOutflowPressureBCTempl : public FVDirichletBCBase, public T
{
public:
  static InputParameters validParams();
  INSFVOutflowPressureBCTempl(const InputParameters & params);

  ADReal boundaryValue(const FaceInfo & /* fi */) const override;

protected:
  /// the dimension of the domain
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<ADReal> & _u;
  /// y-velocity
  const Moose::Functor<ADReal> * const _v;
  /// z-velocity
  const Moose::Functor<ADReal> * const _w;

  /// density
  const Moose::Functor<ADReal> & _rho;

  /// dynamic viscosity
  const Moose::Functor<ADReal> & _mu;

  using FVDirichletBCBase::_t;
  using FVDirichletBCBase::_var;
  using FVDirichletBCBase::determineState;
  using FVDirichletBCBase::getFunction;
  using FVDirichletBCBase::getPostprocessorValue;
  using FVDirichletBCBase::isParamValid;
  using FVDirichletBCBase::mooseError;
  using FVDirichletBCBase::paramError;
  using FVDirichletBCBase::singleSidedFaceArg;
};

typedef INSFVOutflowPressureBCTempl<INSFVFullyDevelopedFlowBC> INSFVOutflowPressureBC;
