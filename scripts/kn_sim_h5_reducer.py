from pathlib import Path
import sys
import numpy as np
import h5py as h5
from astropy import units
from astropy import constants

def run(filenames, output_filename, Nθ_new, Nλ_new):

    dL = (10 * units.pc).cgs

    for filename in filenames:
        with h5.File(filenames[0], "r") as f:
            fλ = f['fla_cgs_per_angstrom'][...]
            t = f['t_days'][...]
            λe = f['lambda_cm'][...]
            θe = f['theta_rad'][...]
            md = f['md_Msolar'][...]
            mw = f['mw_Msolar'][...]
            vd = f['vd_c'][...]
            vw = f['vw_c'][...]
            topo = f['topo'][...]
            wind = f['wind'][...]
            
        Nt, Nλ, Nθ = fλ.shape
       
        idx_λ = (np.arange(Nλ_new+1)*Nλ) // Nλ_new
        idx_θ = (np.arange(Nθ_new+1)*Nθ) // Nθ_new

        λe_new = np.empty((Nλ_new, 2), dtype=np.float32)
        θe_new = np.empty((Nθ_new, 2), dtype=np.float32)
        fλ_new = np.zeros((Nt, Nλ_new, Nθ_new), dtype=np.float32)

        λe_new[:, 0] = λe[idx_λ[:-1], 0]
        λe_new[:, 1] = λe[idx_λ[1:]-1, 1]
        
        θe_new[:, 0] = θe[idx_θ[:-1], 0]
        θe_new[:, 1] = θe[idx_θ[1:]-1, 1]

        dλ = λe[:, 1] - λe[:, 0]
        dλ_new = λe_new[:, 1] - λe_new[:, 0]

        dΩ = 2*np.pi * (np.cos(θe[:, 0]) - np.cos(θe[:, 1]))
        dΩ_new = 2*np.pi * (np.cos(θe_new[:, 0]) - np.cos(θe_new[:, 1]))


        for i in range(Nλ_new):
            for j in range(Nθ_new):
                fλ_new[:, i, j] = Nθ * (
                        fλ[:, idx_λ[i]:idx_λ[i+1], idx_θ[j]:idx_θ[j+1]]
                        * dλ[idx_λ[i]:idx_λ[i+1]][None, :, None] 
                        * dΩ[idx_θ[j]:idx_θ[j+1]][None, None, :]
                        ).sum(axis=(1, 2)) / (
                            dλ_new[i] * dΩ_new[j])

        Lλ_new = 4*np.pi * dL**2 * fλ_new

        output_name = output_dir / (filename.stem + ".reduced.h5")

        print("Saving", output_name)

        with h5.File(output_name, "w") as f:
            f.create_dataset('t_days', data=t)
            f.create_dataset('lambda_cm', data=λe_new)
            f.create_dataset('theta_rad', data=θe_new)
            f.create_dataset('Lla_erg_per_s_angstrom', data=Lλ_new)

            f.create_dataset('md_Msolar', data=md)
            f.create_dataset('mw_Msolar', data=mw)
            f.create_dataset('vd_c', data=vd)
            f.create_dataset('vw_c', data=vw)
            f.create_dataset('topo', data=topo)
            f.create_dataset('wind', data=wind)


if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("usage: python3 kn_sim_h5_reducer [input.h5 ...] output_dir")
        sys.exit()

    filenames = [Path(x) for x in sys.argv[1:-1]]
    output_dir = Path(sys.argv[-1])

    run(filenames, output_dir, 10, 128)

