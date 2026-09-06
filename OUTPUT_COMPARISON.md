# Output comparison with CONGRUENTS-c

This comparison uses the completed 11-galaxy production outputs in each
folder. Spectral component changes below are ratios of the L2 norm across the
tabulated energy grid. Radio changes refer to the sum of all five luminosity
columns at 1.49 GHz. Positive values mean the corrected result is larger.

| Galaxy | IC | Synchrotron | Bremsstrahlung | Total radio | Gamma luminosity |
|---|---:|---:|---:|---:|---:|
| M31 | +109.45% | -15.16% | +0.32% | +6.97% | +44.20% |
| Arp 220 | -9.12% | +1.83% | +1.20% | +2.63% | +0.30% |
| NGC 253 | +24.34% | -3.82% | +0.92% | +4.50% | +8.30% |
| NGC 2403 | +40.41% | -4.01% | +0.29% | +2.62% | +9.77% |
| NGC 2146 | +2.40% | +1.40% | +2.00% | +5.96% | +1.62% |
| M82 | +24.55% | -0.94% | +1.51% | +4.63% | +4.69% |
| Arp 220 corrected | -9.46% | +0.62% | +0.92% | -3.44% | +0.27% |
| M33 | +76.91% | -8.13% | +0.19% | +3.10% | +23.92% |
| LMC | +33.94% | -5.25% | +0.32% | +4.31% | +14.87% |
| SMC | +41.88% | -4.59% | +0.17% | +3.63% | +15.13% |
| LMC corrected | +1125.75% | +776.29% | +835.99% | +155.06% | +595.85% |

The extreme corrected-LMC changes arise mainly because the legacy source
reduced primary-electron injection by a factor of ten in that special catalogue
entry. That suppression is not stated in Section 5.2.2 of the paper and has
been removed here.

Across every galaxy, removing the erroneous extra `1/f_e` factor increases
the free-free optical-depth norm and the free-free-emission norm by exactly
10 per cent. The 1 GHz Gaunt-factor branch correction adds a much smaller,
frequency-dependent change around 1-10 GHz.

The IC changes include both the doubled non-CMB radiation-field normalization
and correction of the CMB/FIR 3-D kernel-selection bug. The final IC response
is therefore galaxy-dependent: the stronger target field changes both photon
production and the steady-state electron population.
