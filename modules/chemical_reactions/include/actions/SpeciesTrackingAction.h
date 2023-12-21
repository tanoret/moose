//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AddFunctionAction.h"

/**
 * The SpeciesTrackingAction sets up user objects, aux kernels, and aux variables
 * for a thermochemistry calculation using Thermochimica.
 */
class SpeciesTrackingAction : public Action
{
public:
  static InputParameters validParams();
  SpeciesTrackingAction(const InputParameters & params);

  virtual void act() override;

  // protected:
};
