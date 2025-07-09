from scipy.integrate import quad
import re
import numpy as np

file_path = "Larsen.i"

def const_internal(n1, n2):
    muc = np.sqrt(1 - (n2/n1)**2)
    
    def rho_internal(mu):
        if abs(mu) < abs(muc):
            return 1.0
        
        sin1 = np.sqrt(1.0 - mu**2)
        sin2 = (n1/n2) * sin1
        mu2 = np.sqrt(1.0 - sin2 **2)

        # Rs = ((n1*mu - n2*mu2)/(n1*mu2 + n2*mu)) ** 2
        # Rp = ((n1*mu2 - n2*mu)/(n1*mu2 + n2*mu)) ** 2

        # return 0.5*(Rs+Rp)
        prefix = 0.5 * ((n1*mu-n2*mu2)/(n1*mu + n2*mu2))**2
        parse = 1 + ((n2*mu*mu2-n1*(1-mu*mu))/(n2*mu*mu2+n1*(1-mu*mu)))**2

        return prefix * parse
    
    def P2(mu):
        return (3/2*mu**2 - 1/2)
    
    def P3(mu):
        return (5/2*mu**3 - 3/2*mu)
    
    alpha_integrand = lambda mu : 1.0 - rho_internal(mu)
    alpha_I, _ = quad(alpha_integrand, 0.0,1.0)
    alpha = 2.0 * n1 * alpha_I

    r1_integrand = lambda mu : mu * rho_internal(mu)
    r1, _ = quad(r1_integrand, 0.0, 1.0)

    r2_integrand = lambda mu : mu**2 * rho_internal(mu)
    r2, _ = quad(r2_integrand, 0.0, 1.0)

    r3_integrand = lambda mu : mu**3 * rho_internal(mu)
    r3, _ = quad(r3_integrand, 0.0, 1.0)

    r4_integand = lambda mu : mu * P3(mu) * rho_internal(mu)
    r4, _ = quad(r4_integand, 0.0, 1.0)

    r5_integrand = lambda mu: P3(mu) * rho_internal(mu)
    r5, _ = quad(r5_integrand, 0.0, 1.0)

    r6_integrand = lambda mu: P2(mu) * P3(mu) * rho_internal(mu)
    r6, _ = quad(r6_integrand, 0.0, 1.0)

    r7_integrand = lambda mu: P3(mu) * P3(mu) * rho_internal(mu)
    r7, _ = quad(r7_integrand, 0.0, 1.0)

    rho1 = (1 - 2*r1)*np.pi
    rho3 = -(1/4 + 2*r5)*np.pi

    gamma1 = 5/7*(1 - 3*np.sqrt(6/5))
    gamma2 = 5/7*(1 + 3*np.sqrt(6/5))

    w0 = 1/(gamma2 - gamma1)

    A1 = (1-2*r1)/4
    A2 = (1-8*r3)*5/16
    A3 = (1+3*r2)/6
    A4 = ((1+3*r2)/3 + (3*r4/2))  #*2/3

    B1 = -(1+8*r5)/16
    B2 = (1-8*r6)*5/16
    B3 = 3*r4/6
    B4 = r4 + 3/14*(1+7*r7)
    
    C1 = w0 * (gamma2*A1 - A2)
    C2 = w0 * (A2 - gamma1*A1)
    C3 = w0 * (gamma2*A3 - A4)
    C4 = w0 * (A4 - gamma1*A3)

    D1 = w0 * (gamma2*B1 - B2)
    D2 = w0 * (B2 - gamma1*B1)
    D3 = w0 * (gamma2*B3 - B4)
    D4 = w0 * (B4 - gamma1*B3)

    D = C3*D4 - D3*C4

    alpha1 = (C1*D4 - D1*C4) / D
    alpha2 = (C3*D2 - C2*D3) / D
    beta1 = (C3*D1 - D3*C1) / D
    beta2 = (C2*D4 - D2*C4) / D
    eta1 = (D4*rho1 - C4*rho3) / D
    eta2 = (C3*rho3 - D3*rho1) / D

    return alpha, alpha1, alpha2, beta1, beta2, eta1, eta2

n1 = None
n2 = None

with open(file_path, 'r') as f:
    lines = f.readlines()

for line in lines:
    if line[:1] in (" ", "\t"):
        continue

    stripped = line.lstrip()
    if stripped.startswith("n1 ="):
        n1 = float(stripped.split("=")[1].strip())
    elif stripped.startswith("n2 = "):
        n2 = float(stripped.split("=")[1].strip())

if n1 is None or n2 is None:
    raise ValueError("Cannot find n1 or n2")

alpha, alpha1, alpha2, beta1, beta2, eta1, eta2 = const_internal(n1, n2)

def replace_line(name, value):
    for i, line in enumerate(lines):
        if line[:1] in (" ", "\t"):
            continue

        if line.lstrip().startswith(f"{name} = "):
            lines[i] = f"{name} = {value}\n"

replace_line("alpha", alpha)
replace_line("alpha1", alpha1)
replace_line("alpha2", alpha2)
replace_line("beta1", beta1)
replace_line("beta2", beta2)
replace_line("eta1", eta1)
replace_line("eta2", eta2)

with open(file_path, 'w') as f:
    f.writelines(lines)

print(f"{file_path} has updated with n1 {n1}, n2 {n2}")