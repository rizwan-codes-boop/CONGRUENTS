/* Independent executable calls the original function, not the new API. */
#include <stdio.h>
#include "../CR_spectra/ionisation.h"
int main(void) {
    double energies[] = {0.001, 0.01, 1., 100., 1.e5};
    double densities[] = {0., 1.e-3, 1., 1.e3};
    for (unsigned d=0; d<4; ++d)
        for (unsigned e=0; e<5; ++e)
            printf("%.17g\n", dEdtm1_ion__GeVsm1(energies[e], densities[d]));
    return 0;
}

