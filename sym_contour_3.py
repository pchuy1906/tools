import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse

plt.rcParams.update({'font.size': 16})

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate and plot 2D distribution.')
    parser.add_argument("--file_input_density", default='density.dat', help='Input file for density profile')
    parser.add_argument("--file_input", default='contour.dat', help='Input file for contour data')
    parser.add_argument("--isym", type=int, default=0, help='Symmetry flag')
    parser.add_argument('--xlim', nargs=2, type=float, required=True, help='X-axis limits (min max)')
    return parser.parse_args()

def read_density_profile(filename):
    F = np.loadtxt(filename)
    X0 = F[:, 0]
    Y0 = F[:, 1]
    indices = np.where(Y0 > 0)[0]
    idmin = indices[0]
    idmax = indices[-1]
    xshift = 0.5 * (X0[idmin] + X0[idmax])
    print(f"Density profile is symmetric with the center z = {xshift:.3f}")
    return xshift, idmin, idmax, (idmin + idmax) // 2

def read_contour_data(filename):
    with open(filename, "r") as f:
        Ny, Nx = map(int, f.readline().split())
        y = np.zeros(Ny)
        x = np.zeros(Nx)
        fyx = np.zeros((Ny, Nx))
        for iy in range(Ny):
            for ix in range(Nx):
                tmp = f.readline().split()
                y[iy] = float(tmp[0])
                x[ix] = float(tmp[1])
                fyx[iy, ix] = float(tmp[2])
    return x, y, fyx, Ny, Nx

def apply_symmetry(fyx, isym, Ny, idmin, idmax, idmiddle):
    if isym == 0:
        return fyx
    for iy in range(Ny):
        for i in range(idmiddle - idmin + 1):
            if isym == -1:
                f1 = 0.5 * (fyx[iy, idmin + i] + fyx[Ny - 1 - iy, idmax - i])
                fyx[iy, idmin + i] = f1
                fyx[Ny - 1 - iy, idmax - i] = f1
            elif isym == 1:
                f1 = 0.5 * (fyx[iy, idmin + i] + fyx[iy, idmax - i])
                fyx[iy, idmin + i] = f1
                fyx[iy, idmax - i] = f1
    return fyx

def compute_2D_distribution(x, y, fyx, xlim, file_input):
    collect_ave = []
    for i in range(len(y)):
        mask = (x >= xlim[0]) & (x <= xlim[1])
        a_tmp = np.sum(fyx[i, mask])
        collect_ave.append(a_tmp)
    res = np.vstack((y, collect_ave)).T
    newfile = "DISTRIBUTION_" + file_input
    np.savetxt(newfile, res, fmt='%15.9f %15.9f')

def main():
    args = parse_arguments()
    xshift, idmin, idmax, idmiddle = read_density_profile(args.file_input_density)
    x, y, fyx, Ny, Nx = read_contour_data(args.file_input)
    x = x - xshift
    fyx = apply_symmetry(fyx, args.isym, Ny, idmin, idmax, idmiddle)

    Y, X = np.meshgrid(x, y)
    maxval = 0.15
    fyx[0,0] = maxval
    contour_min = 0.0
    contour_max = np.amax(fyx)
    v = np.linspace(contour_min, contour_max, 2000)
    v2 = np.linspace(0.0, maxval, 4)

    plt.figure(figsize=(8, 6))
    cf = plt.contourf(Y, X, fyx, levels=v, cmap=cm.jet, vmin=0.0, vmax=maxval)
    cbar = plt.colorbar(cf, shrink=1.00, ticks=v2)
    cbar.set_ticklabels([f"{tick:.2f}" for tick in v2])

    plt.xlim(args.xlim[0], args.xlim[1])
    xtick_positions = np.arange(args.xlim[0], args.xlim[1], 4) + 1
    plt.xticks(xtick_positions)
    plt.ylim(y[0], y[-1])
    plt.xlabel(r'$z$ $(\AA)$')
    plt.ylabel(r'$q$')

    # Output file naming
    file_base = args.file_input.split('.')[0]
    pre_out = f"sym_{file_base}" if args.isym != 0 else file_base
    plt.savefig(pre_out + '.png', bbox_inches='tight')

    compute_2D_distribution(x, y, fyx, args.xlim, args.file_input)

if __name__ == "__main__":
    main()
