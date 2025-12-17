import numpy as np
import matplotlib.pyplot as plt

# Load data
dft_data = np.loadtxt('data_DFT.dat')
lammps_data = np.loadtxt('data_LAMMPS.dat')
lammps_data0 = np.loadtxt('data_LAMMPS0.dat')

# Extract columns
x_dft, y_dft = dft_data[:, 0], dft_data[:, 1]
x_lammps, y_lammps = lammps_data[:, 0], lammps_data[:, 1]
x_lammps0, y_lammps0 = lammps_data0[:, 0], lammps_data0[:, 1]


# Plot
plt.figure(figsize=(4,4))

plt.plot(x_dft, y_dft, 'o-', label='DFT')
plt.plot(x_lammps0, y_lammps0, 'b-', label='FF-LB')
plt.plot(x_lammps, y_lammps, 's-', label='FF-Opt')


plt.xlabel(r'$\mathrm{R_{CM} (\AA)}$')
plt.ylabel('Energy (kcal/mol)')
plt.ylim(-5,5)
plt.xlim(3,10)

plt.legend()
plt.tight_layout()

# Save as PDF
plt.savefig('comparison_plot.pdf')

# (Optional) Show the plot
# plt.show()
