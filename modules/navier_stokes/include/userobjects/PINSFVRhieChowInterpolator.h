//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVRhieChowInterpolator.h"
#include "CellCenteredMapFunctor.h"
#include <unordered_map>

class PINSFVRhieChowInterpolator : public INSFVRhieChowInterpolator
{
public:
  static InputParameters validParams();
  PINSFVRhieChowInterpolator(const InputParameters & params);
  VectorValue<ADReal>
  getVelocity(Moose::FV::InterpMethod m, const FaceInfo & fi, THREAD_ID tid) const override;

  void meshChanged() override;
  void residualSetup() override;

protected:
  Moose::Functor<ADReal> & _eps;
  std::vector<const Moose::Functor<ADReal> *> _epss;
  const unsigned short _smoothing_layers;
  std::vector<const FaceInfo *> _geometric_fi;

  CellCenteredMapFunctor<ADReal, std::unordered_map<dof_id_type, ADReal>> _smoothed_eps;

  /// Whether to force the Rhie-Chow correction of velocity in regions of porosity gradients. This
  /// data member will probably be removed soon since we now have the ability to smooth the porosity
  /// using successive interpolations and reconstructions such that we should always do the RC
  /// correction
  const bool _force_rc_correction;

private:
  void pinsfvSetup();
  bool _initial_setup_done = false;
};
