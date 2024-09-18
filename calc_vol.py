import numpy as np

# Concentration calculator
def calc_conc(n, l, spacing):
    return n/((l*spacing)**3)*(1/(1E-8))**3*1000/(6.022E23)

def calc_box_length(n, conc, spacing):
    return int(((n / ((conc / 1000) * (6.022E23))) ** (1 / 3) / spacing)*(1E8))

def calc_n(l, conc, spacing):
    return int(((conc / 1000) * (6.022E23) * (l * spacing) ** 3) / (1 / (1E-8)) ** 3)

# N corresponds to total # of associated species
# So for example, n=300 corresponds to 300 CaCO3 or 300 NaCl or 300 ZnCl4
n = 300
# Spacing corresponds to the lattice spacing you are using, for now it's just 1
spacing = 1
# Concentrations are the concentrations you want in units of  !!!! --> mM <-- !!!!
concentrations = [25, 50, 100, 200]
print("Conc (mM) | Box Length (Spacings) -- Spacing = 1")
for concentration in concentrations:
    c = concentration / 1000
    print(f"{c*1000:<10} {calc_box_length(n, c, spacing):>10}")

