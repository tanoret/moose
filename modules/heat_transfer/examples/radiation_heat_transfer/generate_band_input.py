import numpy as np
import sys

sys.stdout = open('input.i', 'w', encoding = 'utf-8')

# input parameters

number_of_bands = 8
wavelengths = [0.2, 3.0, 3.5, 4.0, 4.5, 5.5, 6.0, 7.0] # Unit: micro-meter
opacities = [0.4, 0.5, 7.7, 15.45, 27.98, 267.98, 567.32, 7136.06] # Unit: 1/m
# opacities = [1.0, 1.0, 1.0, 1.0 ,1.0 ,1.0 ,1.0 ,1.0]

Nints = [1, 1, 1, 1, 1, 1, 1, 1]
# Nints = [2, 2, 2, 2, 2, 2, 2, 2]
# Nints = [3, 3, 3, 3, 3, 3, 3, 3]
min_wavelength = 0.05 # Unit: 1/m

# Check and calculate

if(len(wavelengths) != number_of_bands):
    print("Wrong input for wavelengths")
    exit(-1)

if(len(opacities) != number_of_bands):
    print("Wrong input for opacities")
    exit(-1)

if(len(Nints) != number_of_bands):
    print("Wrong input for integral numbers")
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

prev = -1
for i in range (2, number_of_bands+1):
    N = Nints[number_of_bands +1 -i]
    for split in range(0, N):
        calc =  c / ((wavelengths[number_of_bands -i] * split / N + wavelengths[number_of_bands +1 -i]* (1-(split / N ))) * 1e-6)
        if(prev != -1):
            if(split == 0):
                print(f"dnu{i-2}split{Nints[number_of_bands +1 -i]} = {calc - prev}")
            else:
                print(f"dnu{i-1}split{split} = {calc - prev}")
        print("")
        print(f"nu{i-1}split{split+1} = {calc}")
        prev = calc

for split in range(0, Nints[0]):
    calc =  c / ((min_wavelength* split / Nints[0] + wavelengths[0] * (1-(split / Nints[0]))) * 1e-6)
    if(split == 0):
                print(f"dnu{number_of_bands-1}split{Nints[number_of_bands - 1]} = {calc - prev}")
    else:
        print(f"dnu{number_of_bands}split{split} = {calc - prev}")
    print("")
    print(f"nu{number_of_bands}split{split+1} = {calc}")
    prev = calc

print(f"dnu{number_of_bands}split{split+1} = {c/(min_wavelength*1e-6) - calc}")
print("")

for i in range (0, number_of_bands):
    print(f"kappa{i+1} = {kappa[i]}")
print("")

# for i in range(0, number_of_bands-1):
#     print(f"nu{i+1}to{i+2}split = {(nu[i+1] - nu[i])/Nints[i]}")
# print(f"nu{number_of_bands}to{number_of_bands+1} = {nu[number_of_bands-1] - nu[number_of_bands-2]}")
# print("")

print("[Variables]")
print("  [T]")
print("    type = MooseVariableFVReal")
print("    initial_condition = ${T0}")
print("  []")
print("")

for band in range(1, number_of_bands+1):
    for order in range(1,3):
        for split in range(1, Nints[number_of_bands - band]+1):
            if(split < 10 ):
                print(f"  [psi{order}{band}{split}]")
                print(f"      type = MooseVariableFVReal")
                print("  []")
                print("")
print("[]")
print("")

print("[FVKernels]")

for band in range(1, number_of_bands+1):
    for order in range(1,3):
        for split in range(1, Nints[number_of_bands - band]+1):
            print(f"  [diffusion{order}{band}{split}]")
            print("      type = FVSP3ThermalRadiationDiffusion")
            print(f"      variable = psi{order}{band}{split}")
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
        for split in range(1, Nints[number_of_bands - band]+1):
            print(f"  [source{order}{band}{split}]")
            print("      type = FVSP3ThermalRadiationSourceSink")
            print(f"      variable = psi{order}{band}{split}")
            print("      T = 'T'")
            print(f"      nu = ${{nu{band}split{split}}}")
            print(f"      refraction_index = ${{n1}}")
            print(f"      kappa = ${{kappa{band}}}")
            print("  []")
            print("")
        
print("  [energy_source]")
print("      type = FVSP3TemperatureSourceSink")
print("      variable = T")

print("      absorptivities = '", end = "")
for i in range (1, number_of_bands+1):
    for split in range(1, Nints[number_of_bands - i]+1):
        print(f"${{kappa{i}}} ", end = "")
print("'")

print("      psi_1 = '", end = "")
for i in range (1, number_of_bands+1):
    for split in range(1, Nints[number_of_bands - i]+1):
        print(f"psi1{i}{split} ", end = "")
print("'")

print("      psi_2 = '", end = "")
for i in range (1, number_of_bands+1):
    for split in range(1, Nints[number_of_bands - i]+1):
        print(f"psi2{i}{split} ", end = "")
print("'")

print("      band_frequency_width = '", end = "")
for i in range (1, number_of_bands+1):
    for split in range(1, Nints[number_of_bands - i]+1):
        print(f"${{dnu{i}split{split}}} ", end = "")
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
        for split in range(1, Nints[number_of_bands - band]+1):
            print(f"  [BC{order}{band}{split}]")
            print("      type = FVSP3ThermalRadiationBC")
            print("      boundary = 'left right'")
            print(f"      variable = psi{order}{band}{split}")
            print("      T = 'T'")
            print(f"      nu = ${{nu{band}split{split}}}")
            print(f"      refraction_index = ${{n1}}")
            print(f"      kappa = ${{kappa{band}}}")
            print("      epsilon = ${epsilon}")
            
            if order == 1:
                print(f"      psi = 'psi2{band}{split}'")
                print(f"      order = first")
            else:
                print(f"      psi = 'psi1{band}{split}'")
                print(f"      order = second")
            print("  []")
            print("")


print("    [BC_temperature]")
print("        type = FVSP3TemperatureBC")
print("        boundary = 'left right'")
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