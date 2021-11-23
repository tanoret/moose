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
  const unsigned short _rec;
  std::vector<const FaceInfo *> _geometric_fi;

  CellCenteredMapFunctor<ADReal, std::unordered_map<dof_id_type, ADReal>> _reconstructed_eps;

  /// Whether the porosity field is smooth or has discontinuities
  const bool _smooth_porosity;

private:
  void pinsfvSetup();
  bool _initial_setup_done = false;
};
