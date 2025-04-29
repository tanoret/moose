import numpy as np

# input parameters

number_of_bands = 8
wavelengths = [0.2, 3.0, 3.5, 4.0, 4.5, 5.5, 6.0, 7.0] # Unit: micro-meter
# opacities = [0.4, 0.5, 7.7, 15.45, 27.98, 267.98, 567.32, 7136.06] # Unit: 1/m
opacities = [1.0, 1.0, 1.0, 1.0 ,1.0 ,1.0 ,1.0 ,1.0]

# Check and calculate

if(len(wavelengths) != number_of_bands):
    print("Wrong input for wavelengths")
    exit(-1)

if(len(opacities) != number_of_bands):
    print("Wrong input for opacities")
    exit(-1)

c = 299792458

nu = [0] * number_of_bands
kappa = [0] * number_of_bands
nuhalf = [0] * number_of_bands

for i in range(1, number_of_bands+1):
    nu[i-1] = c / (wavelengths[number_of_bands - i] * 1e-6)
    kappa[i-1] = opacities[number_of_bands -i]
for i in range(2, number_of_bands+1):
    nuhalf[i-2] = 2 * c / ((wavelengths[number_of_bands -i] + wavelengths[number_of_bands +1 -i]) * 1e-6)
nuhalf[number_of_bands-1] = 2 * c / (wavelengths[number_of_bands - 1] * 1e-6)

# Print for .i file
print(f"nu1 = {nu[0]}")
print("")

for i in range (0, number_of_bands):
    print(f"nu{i+1}half = {nuhalf[i]}")
print("")

for i in range (0, number_of_bands):
    print(f"kappa{i+1} = {kappa[i]}")
print("")

for i in range(0, number_of_bands-1):
    print(f"nu{i+1}to{i+2} = {nu[i+1] - nu[i]}")
print(f"nu{number_of_bands}to{number_of_bands+1} = {nu[number_of_bands-1] - nu[number_of_bands-2]}")
print("")

print("[Variables]")
print("  [T]")
print("    type = MooseVariableFVReal")
print("    initial_condition = ${T0}")
print("  []")
print("")

for band in range(1, number_of_bands+1):
    for order in range(1,3):
        print(f"  [psi{order}{band}]")
        print(f"      type = MooseVariableFVReal")
        print("  []")
        print("")
print("[]")
print("")

print("[FVKernels]")

for band in range(1, number_of_bands+1):
    for order in range(1,3):
        print(f"  [diffusion{order}{band}]")
        print("      type = FVSP3ThermalRadiationDiffusion")
        print(f"      variable = psi{order}{band}")
        print("      epsilon = ${epsilon}")
        print(f"      kappa = ${{kappa{band}}}")
        if order == 1:
            print("      order = first")
        else:
            print("      order = second")
        print("  []")
        print("")

for band in range(1, number_of_bands+1):
    for order in range(1,3):
        print(f"  [source{order}{band}]")
        print("      type = FVSP3ThermalRadiationSourceSink")
        print(f"      variable = psi{order}{band}")
        print("      T = 'T'")
        print(f"      nu = ${{nu{band}half}}")
        print(f"      refraction_index = ${{n1}}")
        print(f"      kappa = ${{kappa{band}}}")
        print("  []")
        print("")
        
print("  [energy_source]")
print("      type = FVSP3TemperatureSourceSink")
print("      variable = T")

print("      absorptivities = '", end = "")
for i in range (1, number_of_bands+1):
    print(f"${{kappa{i}}} ", end = "")
print("'")

print("      psi_1 = '", end = "")
for i in range (1, number_of_bands+1):
    print(f"psi1{i} ", end = "")
print("'")

print("      psi_2 = '", end = "")
for i in range (1, number_of_bands+1):
    print(f"psi2{i} ", end = "")
print("'")

print("      band_frequency_width = '", end = "")
for i in range (1, number_of_bands+1):
    print(f"${{nu{i}to{i+1}}} ", end = "")
print("'")

print("  []")
print("")

print("    [energy_time]")
print("        type = FVTimeKernel")
print("        variable = T")
print("    []")
print("")

print("    [energy_diffusion]")
print("        type = FVDiffusion")
print("        variable = T")
print("        coeff = ${k}")
print("  []")

print("[]")
print("")

print("[FVBCs]")

for band in range(1, number_of_bands + 1):
    for order in range(1,3):
        print(f"  [BC{order}{band}]")
        print("      type = FVSP3ThermalRadiationBC")
        print("      boundary = outer")
        print(f"      variable = psi{order}{band}")
        print("      T = 'T'")
        print(f"      nu = ${{nu{band}half}}")
        print(f"      refraction_index = ${{n1}}")
        print(f"      kappa = ${{kappa{band}}}")
        print("      epsilon = ${epsilon}")
        
        if order == 1:
            print(f"      psi = 'psi2{band}'")
            print(f"      order = first")
        else:
            print(f"      psi = 'psi1{band}'")
            print(f"      order = second")
        print("  []")
        print("")


print("    [BC_temperature]")
print("        type = FVSP3TemperatureBC")
print("        boundary =  outer")
print("        variable = T")
print("        Tb = ${Tb}")
print("        n1 = ${n1}")
print("        n2 = ${n2}")
print("        h = ${h}")
print("        k = ${k}")
print("        alpha = ${alpha}")
print("        nu1 = ${nu1}")
print("        Nintegral = 100")
print("    []")
print("[]")
print("")