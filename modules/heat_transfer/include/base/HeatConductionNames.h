//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <string>

namespace HeatConduction
{
static const std::string emissivity = "emissivity";
static const std::string T_ambient = "T_ambient";

namespace Constants
{
static const Real sigma = 5.670374419e-8;
static const Real hp = 6.62607015e-34; // J.s
static const Real c = 2.99792458e8; // m/s
static const Real kb = 1.380649e-23; // J/K
}
}
