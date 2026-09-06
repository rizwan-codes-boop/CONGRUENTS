from congruents import Context

with Context(threads=4) as model:
    print("C version:", model.version, "OpenMP:", model.openmp_enabled)
    print("Ionisation loss [GeV/s]:",
          model.ionisation_loss([0.001, 1., 100.], density_cm3=1.))

