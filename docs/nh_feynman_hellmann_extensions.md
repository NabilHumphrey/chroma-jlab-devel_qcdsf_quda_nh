# Feynman-Hellmann Extensions for Chroma

This document describes the new inline measurements and modifications added to
support Feynman-Hellmann (FH) calculations for hadronic structure functions,
including single-hadron and multi-hadron (e.g., deuteron) observables.

## Overview

The extensions add the following new inline measurements:

| Measurement | Description |
|-------------|-------------|
| `BARYON_SPECTRUM_FH` | Single-baryon correlators with FH support |
| `MESON_SPECTRUM_FH` | Single-meson correlators with FH support |
| `MULTI_HADRON_FH` | Multi-hadron (dibaryon, dimeson) correlators with spin projection |

These measurements are designed to work with perturbed and unperturbed
propagators to extract matrix elements via the Feynman-Hellmann theorem.

---

## New Files

### `lib/meas/inline/hadron/inline_baryon_spectrum_fh.{cc,h}`

Single-baryon spectrum measurement supporting:
- Proton and neutron correlators
- Flexible flavor assignment (perturbed/unperturbed for up/down quarks)
- Multiple contraction combinations in a single measurement

### `lib/meas/inline/hadron/inline_meson_spectrum_fh.{cc,h}`

Single-meson spectrum measurement supporting:
- Pseudoscalar mesons (pion, kaon)
- Vector mesons (rho) with all three polarizations (x, y, z)
- Flexible flavor assignment for quark/antiquark

### `lib/meas/inline/hadron/inline_multi_hadron_fh.{cc,h}`

Multi-hadron correlator measurement using the product ansatz:

```
C_AB(t) = (1/2) * sum_p C_A(p,t) * C_B(-p,t)
```

Supports:
- **DIBARYON**: Proton-proton, proton-neutron (deuteron), neutron-neutron
- **DIMESON**: Rho-rho (di-rho)

All states are computed as spin-1 in the **Cartesian polarization basis** (x, y, z).

### `lib/meas/inline/hadron/inline_deuteron_correlator.{cc,h}`

Earlier implementation of deuteron correlator (may be deprecated in favor of
the more general `MULTI_HADRON_FH`).

---

## Modified Files

### `lib/meas/inline/hadron/inline_hadron_aggregate.cc`

Registers the new inline measurements with the Chroma factory system.

### `lib/CMakeLists.txt`

Adds the new source files to the build system.

### `lib/meas/hadron/simple_meson_seqsrc_w.{cc,h}`

Modifications to support momentum-projected meson sequential sources.

### `lib/meas/sources/mom_source_const.cc`

Modifications to momentum source construction.

---

## MULTI_HADRON_FH Details

### Cartesian Polarization Basis

The measurement outputs spin-1 correlators in the Cartesian basis rather than
the spherical (m = +1, 0, -1) basis. This simplifies the physics interpretation:

**For dibaryons (1/2 x 1/2 -> 1):**
```
C_x(t) = C_y(t) = 0.5 * [G_up1 * G_up2 + G_down1 * G_down2]
C_z(t) = 0.5 * [G_up1 * G_down2 + G_down1 * G_up2]
```

Note: C_x = C_y identically (bit-for-bit) for the product ansatz, but both are
output for uniformity with other states.

**For di-rho (1 x 1 -> 1, antisymmetric):**

Uses Levi-Civita (cross product) coupling:
```
C_x(t) = 0.5 * sum_p [C_y^rho1(p) * C_z^rho2(-p) + C_z^rho1(p) * C_y^rho2(-p)]
C_y(t) = 0.5 * sum_p [C_z^rho1(p) * C_x^rho2(-p) + C_x^rho1(p) * C_z^rho2(-p)]
C_z(t) = 0.5 * sum_p [C_x^rho1(p) * C_y^rho2(-p) + C_y^rho1(p) * C_x^rho2(-p)]
```

### Sum Relations

The total unpolarized correlator is:
```
C_total = C_x + C_y + C_z
```

### Tensor Polarization

The tensor polarization observable b1 can be extracted from:
```
b1 ~ C_z - 0.5 * (C_x + C_y)
```

### Momentum Structure

- Multi-hadron states are computed at **total momentum zero**
- The product ansatz sums over all relative momenta up to `mom2_max`
- Single-hadron correlators are momentum-projected at all |p|^2 <= mom2_max

---

## XML Input Format

### MULTI_HADRON_FH Example

```xml
<elem>
    <Name>MULTI_HADRON_FH</Name>
    <Frequency>1</Frequency>
    <Param>
        <mom2_max>4</mom2_max>
        <avg_equiv_mom>true</avg_equiv_mom>
        <t0>0</t0>

        <states>
            <elem>
                <name>deuteron</name>
                <type>DIBARYON</type>
                <hadron1>PROTON</hadron1>
                <hadron2>NEUTRON</hadron2>
            </elem>
            <elem>
                <name>pp_dibaryon</name>
                <type>DIBARYON</type>
                <hadron1>PROTON</hadron1>
                <hadron2>PROTON</hadron2>
            </elem>
            <elem>
                <name>dirho</name>
                <type>DIMESON</type>
                <hadron1>RHO</hadron1>
                <hadron2>RHO</hadron2>
            </elem>
        </states>

        <contractions>
            <elem>
                <up_is_perturbed>false</up_is_perturbed>
                <down_is_perturbed>false</down_is_perturbed>
                <label>unperturbed</label>
            </elem>
            <elem>
                <up_is_perturbed>true</up_is_perturbed>
                <down_is_perturbed>false</down_is_perturbed>
                <label>up_perturbed</label>
            </elem>
        </contractions>
    </Param>
    <NamedObject>
        <gauge_id>default_gauge_field</gauge_id>
        <perturbed_prop_id>fh_prop</perturbed_prop_id>
        <unperturbed_prop_id>std_prop</unperturbed_prop_id>
        <output_file>multi_hadron.dat</output_file>
    </NamedObject>
</elem>
```

### State Types

| Type | Supported hadron1/hadron2 |
|------|---------------------------|
| `DIBARYON` | `PROTON`, `NEUTRON` |
| `DIMESON` | `RHO`, `RHO_X`, `RHO_Y`, `RHO_Z` |

### Notes

- The `total_spin` parameter has been **removed**. All states are spin-1.
- Spin-0 (singlet) and spin-2 states are no longer supported.

---

## Output Format

### XML Output

```xml
<state>
    <name>deuteron</name>
    <type>DIBARYON</type>
    <hadron1>PROTON</hadron1>
    <hadron2>NEUTRON</hadron2>
    <multi_hadron>
        <elem><t_eff>0</t_eff><re>...</re><im>...</im></elem>
        ...
    </multi_hadron>
    <pol_x>
        <elem><t_eff>0</t_eff><re>...</re><im>...</im></elem>
        ...
    </pol_x>
    <pol_y>...</pol_y>
    <pol_z>...</pol_z>
</state>
```

### Plain Text Output

```
# Columns: state flavor t_eff re im
deuteron unperturbed 0 1.234e-05 5.678e-10
deuteron unperturbed 1 1.123e-05 4.567e-10
...
deuteron_x unperturbed 0 4.112e-06 1.892e-10
deuteron_y unperturbed 0 4.112e-06 1.892e-10
deuteron_z unperturbed 0 4.114e-06 1.893e-10
...
```

---

## Building

The new files are integrated into the CMake build system. No special build
flags are required. Simply rebuild Chroma after pulling the changes:

```bash
cd build
cmake ..
make -j8
```

---

## Physics Background

### Feynman-Hellmann Theorem

For a Hamiltonian H(lambda) with perturbation lambda*O:

```
d E_n / d lambda |_{lambda=0} = <n|O|n>
```

In lattice QCD, this is implemented by computing correlators with both
standard propagators and propagators computed with a modified Dirac operator
that includes the perturbation.

### Product Ansatz

The multi-hadron correlators use the product ansatz, which approximates the
two-hadron correlator as a product of single-hadron correlators:

```
C_AB(p_tot=0, t) = sum_p C_A(p, t) * C_B(-p, t)
```

This is exact in the non-interacting limit and provides the leading
contribution in interacting theories when the hadrons are well-separated.

### Tensor Structure of Deuteron

The deuteron is a spin-1 state. In the Cartesian basis:
- C_x and C_y probe transverse polarizations
- C_z probes longitudinal polarization

The tensor polarization b1 structure function is sensitive to:
```
b1 ~ C_z - (C_x + C_y)/2
```

which vanishes for unpolarized or purely vector-polarized states.

---

## Authors and Version History

- Initial implementation: 2024-2025
- Cartesian basis conversion: February 2025
