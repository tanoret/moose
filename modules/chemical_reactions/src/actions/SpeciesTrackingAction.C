//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SpeciesTrackingAction.h"
#include "ThermochimicaUtils.h"
#include "FEProblemBase.h"
#include "MooseMesh.h"
#include "MooseUtils.h"
#include "MooseUtils.h"
#include "AddVariableAction.h"
#include "libmesh/string_to_enum.h"

registerMooseAction("ChemicalReactionsApp", SpeciesTrackingAction, "add_variable");
registerMooseAction("ChemicalReactionsApp", SpeciesTrackingAction, "add_aux_variable");
registerMooseAction("ChemicalReactionsApp", SpeciesTrackingAction, "add_ic");
registerMooseAction("ChemicalReactionsApp", SpeciesTrackingAction, "add_user_object");
registerMooseAction("ChemicalReactionsApp", SpeciesTrackingAction, "add_aux_kernel");

InputParameters
SpeciesTrackingAction::validParams()
{
  InputParameters params = Action::validParams();

  ThermochimicaUtils::addClassDescription(params, "Sets up species tracking action.");
  return params;
}

SpeciesTrackingAction::SpeciesTrackingAction(const InputParameters & params) : Action(params) {}

void
SpeciesTrackingAction::act()
{
}
