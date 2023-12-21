//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "ElementUserObject.h"
#include "NonADFunctorInterface.h"

// Forward declaration
// class SpeciesTrackingAction;

/**
 * User object that performs a Gibbs energy minimization at each node by calling
 * the Thermochimica code.
 */
class MagicBookSpecies : public ElementUserObject, public NonADFunctorInterface
{
public:
  static InputParameters validParams();
  MagicBookSpecies(const InputParameters & parameters);

  virtual void initialSetup() override;
  virtual void initialize() override{};
  virtual void execute() override;
  virtual void finalize() override;
  virtual void threadJoin(const UserObject &) override {}

  void setSpeciesValue(const Real & set_value,
                       const std::vector<std::string> & species_name,
                       std::map<std::string, MooseVariableFV<Real> *> & species_map);
  void setSpeciesValue(const std::vector<Real> & set_values,
                       const std::vector<std::string> & species_name,
                       std::map<std::string, MooseVariableFV<Real> *> & species_map);
  void updateElementConcentrationFromIsotopes();
  void updateIsotopeConecentrationFromElements();

private:
  /// Naming storage
  std::vector<std::string> _isotopes_name;
  std::vector<std::string> _elements_name;

  /// Size Storage
  size_t _number_of_isotopes;
  size_t _number_of_elements;

  /// Map containing all isotopes
  std::map<std::string, MooseVariableFV<Real> *> _isotopes;

  // Map containign all elements
  std::map<std::string, MooseVariableFV<Real> *> _elements;

  // Releationship mappers
  std::map<std::string, size_t> _isotopeToElementIndex;
  std::map<std::string, std::vector<size_t>> _elementToIsotopeIndex;

  // Intial values
  Real _initial_value_isotopes;
  bool _initial_value_provided_isotopes;
  std::vector<Real> _initial_values_isotopes;
  bool _initial_values_provided_isotopes;

  // Enumerator to determine UO function
  MooseEnum _uo_function;
};
