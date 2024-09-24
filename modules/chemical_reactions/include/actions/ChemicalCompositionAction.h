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
#include "CommonChemicalCompositionAction.h"

/**
 * The ChemicalCompositionAction sets up user objects, aux kernels, and aux variables
 * for a thermochemistry calculation using Thermochimica.
 */
class ChemicalCompositionAction : public Action
{
public:
  static InputParameters validParams();
  ChemicalCompositionAction(const InputParameters & params);

  const std::vector<unsigned int> & elementIDs() const { return _element_ids; }

  const std::vector<std::string> & phases() const { return _phases; }
  const std::vector<std::string> & elementPotentials() const
  {
    return _tokenized_element_potentials;
  }

  const std::vector<std::pair<std::string, std::string>> & speciesPhasePairs() const
  {
    return _tokenized_species;
  }

  const std::vector<std::pair<std::string, std::string>> & vaporPhasePairs() const
  {
    return _tokenized_vapor_species;
  }

  const std::vector<std::pair<std::string, std::string>> & phaseElementPairs() const
  {
    return _tokenized_chemical_potential;
  }

  const std::vector<std::pair<std::string, std::string>> & chemicalPotentialPairs() const
  {
    return _tokenized_phase_elements;
  }

  virtual void act();

protected:
  void readCSV();

  /// Element names
  std::vector<std::string> _elements;

  /// Initial conditions for each element: [element name] => initial condition value
  std::map<std::string, Real> _initial_conditions;

  /// List of phases tracked by Thermochimica
  std::vector<std::string> _phases;

  /// Atomic numbers of the selected elements
  std::vector<unsigned int> _element_ids;

  /// Tokenized versions of the output variables to avoid redoing tokenization
  std::vector<std::pair<std::string, std::string>> _tokenized_species;
  std::vector<std::pair<std::string, std::string>> _tokenized_chemical_potential;
  std::vector<std::string> _tokenized_element_potentials;
  std::vector<std::pair<std::string, std::string>> _tokenized_vapor_species;
  std::vector<std::pair<std::string, std::string>> _tokenized_phase_elements;
};
