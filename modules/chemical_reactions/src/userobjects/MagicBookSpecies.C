//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MagicBookSpecies.h"
#include "MooseVariableFV.h"
#include "SystemBase.h"

registerMooseObject("ChemicalReactionsApp", MagicBookSpecies);

InputParameters
MagicBookSpecies::validParams()
{
  InputParameters params = ElementUserObject::validParams();
  params += NonADFunctorInterface::validParams();
  params.addClassDescription("UO keeping track of species and elements.");
  params.addParam<std::vector<std::string>>("isotopes_name", "Name of all isotopes trackied.");
  params.addParam<Real>("initial_value_isotopes", "Initial value for all isotopes.");
  params.addParam<std::vector<Real>>("intial_values_isotopes",
                                     "Vector of initial values for all isotopes.");

  MooseEnum uo_function("isotopes_to_elements elements_to_isotopes", "isotopes_to_elements");

  params.addParam<MooseEnum>(
      "uo_function",
      uo_function,
      "Function to call in the user object. "
      "'isotopes_to_elements' computes elements concentrations from isotopes concentration. "
      "'elements_to_isotopes' computes isotopes concentrations from updates in element "
      "concentrations.");
  return params;
}

MagicBookSpecies::MagicBookSpecies(const InputParameters & parameters)
  : ElementUserObject(parameters),
    NonADFunctorInterface(this),
    _isotopes_name(parameters.isParamValid("isotopes_name")
                       ? getParam<std::vector<std::string>>("isotopes_name")
                       : std::vector<std::string>{}),
    _initial_value_isotopes(parameters.isParamValid("initial_value_isotopes")
                                ? getParam<Real>("initial_value_isotopes")
                                : 0.0),
    _initial_values_isotopes(parameters.isParamValid("intial_values_isotopes")
                                 ? getParam<std::vector<Real>>("intial_values_isotopes")
                                 : std::vector<Real>{}),
    _uo_function(getParam<MooseEnum>("uo_function"))
{
  if (_isotopes_name.empty())
    mooseWarning("No isotope names were provided");

  // Define the lambda function to extract unique letters
  auto extractElements = [this]() -> std::vector<std::string>
  {
    std::set<std::string> uniqueLetters; // Using unordered_set to store unique letters

    for (const std::string & isotope : this->_isotopes_name)
    {
      // Extract letters from each isotope (first one or two characters)
      std::string letters =
          isotope.substr(0, std::min(isotope.find_first_of("0123456789"), isotope.size()));
      uniqueLetters.insert(letters); // Insert the letters into the set
    }

    return std::vector<std::string>(uniqueLetters.begin(), uniqueLetters.end());
  };

  // Extract unique letters using the lambda function
  _elements_name = extractElements();

  // Populating mappers
  // Populating isotope to element mapper
  for (const std::string & isotope : _isotopes_name)
  {
    auto it = std::find(_elements_name.begin(),
                        _elements_name.end(),
                        isotope.substr(0, isotope.find_first_of("0123456789")));

    // If the related element is found, store its index for the isotope
    if (it != _elements_name.end())
    {
      size_t elementIndex = std::distance(_elements_name.begin(), it);
      _isotopeToElementIndex[isotope] = elementIndex;
    }
  }

  // Populating element to isotopes mapper
  for (const std::string & element : _elements_name)
    _elementToIsotopeIndex[element] = {};

  for (size_t i = 0; i < _isotopes_name.size(); ++i)
  {
    auto it = std::find(_elements_name.begin(),
                        _elements_name.end(),
                        _isotopes_name[i].substr(0, _isotopes_name[i].find_first_of("0123456789")));

    if (it != _elements_name.end())
    {
      size_t elementIndex = std::distance(_elements_name.begin(), it);
      _elementToIsotopeIndex[_elements_name[elementIndex]].push_back(i);
    }
  }

  // Initializing isotopes functors
  _number_of_isotopes = _isotopes_name.size();
  for (auto isotope : _isotopes_name)
  {
    InputParameters loc_params = MooseVariableFVReal::validParams();
    loc_params.addPrivateParam("_moose_app", &getMooseApp());
    auto variable_name = static_cast<const VariableName &>(isotope);
    UserObject::_fe_problem.addAuxVariable("MooseVariableFVReal", variable_name, loc_params);
  }

  // Initializing elements functors
  _number_of_elements = _elements_name.size();
  for (auto element : _elements_name)
  {
    InputParameters loc_params = MooseVariableFVReal::validParams();
    loc_params.addPrivateParam("_moose_app", &getMooseApp());
    UserObject::_fe_problem.addAuxVariable("MooseVariableFVReal", element, loc_params);
  }

  // Error contorl for initial conditions
  _initial_value_provided_isotopes = parameters.isParamValid("initial_value_isotopes");
  _initial_values_provided_isotopes = parameters.isParamValid("intial_values_isotopes");
  if (_initial_value_provided_isotopes && _initial_values_provided_isotopes)
    mooseError("Only 'initial_value_isotopes' or 'intial_values_isotopes' may be provided for the "
               "initial conditions.");
  if (!(_initial_value_provided_isotopes || _initial_values_provided_isotopes))
  {
    mooseWarning("No 'initial_value_isotopes' or 'intial_values_isotopes' provided, initializing "
                 "all isotopes to zero.");
    _initial_value_provided_isotopes = true;
  }
  if (_initial_values_provided_isotopes)
    if (_isotopes_name.size() != _initial_values_isotopes.size())
      mooseError("Number of specified isotopes ",
                 _isotopes_name.size(),
                 " is different that the number of ",
                 _initial_values_isotopes.size(),
                 " initial values provided");
}

void
MagicBookSpecies::initialSetup()
{
  // Linking isotopes dictoinary to functors
  for (auto isotope : _isotopes_name)
  {
    // Determine if the functor is already in the system
    // We need to do this because `getFunctor` will create a template null functor if it is not
    bool functor_type = UserObject::_subproblem.hasFunctorWithType<ADReal>(isotope, _tid);
    if (!functor_type)
      mooseError("The isotope functor named: '",
                 isotope,
                 "' does not exist in the system ",
                 "or returns other type than 'Real' or 'ADReal'");

    const MooseVariableFieldBase * loc_isotope_functor_ptr =
        &UserObject::_subproblem.getVariable(_tid, isotope);
    _isotopes[isotope] = const_cast<MooseVariableFV<Real> *>(
        dynamic_cast<const MooseVariableFV<Real> *>(loc_isotope_functor_ptr));
  }

  // Set all elements to zero
  if (_initial_value_provided_isotopes)
    this->setSpeciesValue(_initial_value_isotopes, _isotopes_name, _isotopes);

  if (_initial_values_provided_isotopes)
    this->setSpeciesValue(_initial_values_isotopes, _isotopes_name, _isotopes);

  // Setting up element functors
  for (auto element : _elements_name)
  {
    // Determine if the functor is already in the system
    // We need to do this because `getFunctor` will create a template null functor if it is not
    bool functor_type = UserObject::_subproblem.hasFunctorWithType<ADReal>(element, _tid);
    if (!functor_type)
      mooseError("The element functor named: '",
                 element,
                 "' does not exist in the system ",
                 "or returns other type than 'Real' or 'ADReal'");

    const MooseVariableFieldBase * loc_element_functor_ptr =
        &UserObject::_subproblem.getVariable(_tid, element);
    _elements[element] = const_cast<MooseVariableFV<Real> *>(
        dynamic_cast<const MooseVariableFV<Real> *>(loc_element_functor_ptr));
  }

  // Elment concentrations are initially set from isotopic concentrations
  this->updateElementConcentrationFromIsotopes();
}

void
MagicBookSpecies::execute()
{
  if (_uo_function == "isotopes_to_elements")
    this->updateElementConcentrationFromIsotopes();

  if (_uo_function == "elements_to_isotopes")
    this->updateIsotopeConecentrationFromElements();
}

void
MagicBookSpecies::finalize()
{
}

void
MagicBookSpecies::updateElementConcentrationFromIsotopes()
{
  // Compute element concentration from isotope concentration
  this->setSpeciesValue(0.0, _elements_name, _elements);
  for (auto element : _elements_name)
  {
    for (const auto * const elem :
         UserObject::_fe_problem.mesh().getMesh().active_local_element_ptr_range())
    {
      for (auto isotope : _isotopes_name)
      {
        size_t element_index = _isotopeToElementIndex[isotope];
        std::string loc_element_name = _elements_name[element_index];
        if (loc_element_name == element)
        {
          auto elemArg = makeElemArg(elem, false);
          Real element_value = raw_value((*_isotopes[isotope])(elemArg, determineState()));
          _elements[element]->sys().solution().add(
              elem->dof_number(_elements[element]->sys().number(),
                               _elements[element]->number(),
                               /*comp=*/0),
              element_value);
        }
      }
    }
    _elements[element]->sys().solution().close();
    _elements[element]->sys().update();
  }
}

void
MagicBookSpecies::updateIsotopeConecentrationFromElements()
{
  // Compute element concentration from isotope concentration
  auto new_state = Moose::currentState();
  auto old_state = Moose::oldState();
  for (const auto & isotope : _isotopes_name)
  {
    auto element_index = _isotopeToElementIndex[isotope];
    auto element = _elements_name[element_index];

    for (const auto * const elem :
         UserObject::_fe_problem.mesh().getMesh().active_local_element_ptr_range())
    {
      auto elemArg = makeElemArg(elem, false);
      Real new_element_value = raw_value((*_elements[element])(elemArg, new_state));
      Real old_element_value = raw_value((*_elements[element])(elemArg, old_state));
      auto delta_element = new_element_value - old_element_value;
      _isotopes[isotope]->sys().solution().add(elem->dof_number(_isotopes[isotope]->sys().number(),
                                                                _isotopes[isotope]->number(),
                                                                /*comp=*/0),
                                               delta_element);
    }
    _isotopes[isotope]->sys().solution().close();
    _isotopes[isotope]->sys().update();
  }
}

void
MagicBookSpecies::setSpeciesValue(const Real & set_value,
                                  const std::vector<std::string> & species_name,
                                  std::map<std::string, MooseVariableFV<Real> *> & species_map)
{
  for (const auto & species : species_name)
  {
    for (const auto * const elem :
         UserObject::_fe_problem.mesh().getMesh().active_local_element_ptr_range())
    {
      species_map[species]->sys().solution().set(
          elem->dof_number(species_map[species]->sys().number(),
                           species_map[species]->number(),
                           /*comp=*/0),
          set_value);
    }
    // Do global assembly
    species_map[species]->sys().solution().close();
    // Scatter from parallel solution vector into ghosted vector which is what things read from
    species_map[species]->sys().update();
  }
}

void
MagicBookSpecies::setSpeciesValue(const std::vector<Real> & set_values,
                                  const std::vector<std::string> & species_name,
                                  std::map<std::string, MooseVariableFV<Real> *> & species_map)
{
  if (set_values.size() != species_name.size())
    mooseError(
        "Set species value called with different size for the number of species and values.");
  for (auto i : make_range(set_values.size()))
  {
    for (const auto * const elem :
         UserObject::_fe_problem.mesh().getMesh().active_local_element_ptr_range())
    {
      species_map[species_name[i]]->sys().solution().set(
          elem->dof_number(species_map[species_name[i]]->sys().number(),
                           species_map[species_name[i]]->number(),
                           /*comp=*/0),
          set_values[i]);
    }
    // Do global assembly
    species_map[species_name[i]]->sys().solution().close();
    // Scatter from parallel solution vector into ghosted vector which is what things read from
    species_map[species_name[i]]->sys().update();
  }
}
