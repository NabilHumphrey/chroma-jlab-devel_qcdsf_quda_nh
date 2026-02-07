#!/usr/bin/env python3
"""
Compare hadron correlators from two Chroma HADRON_SPECTRUM XML output files.

Usage:
    python compare_correlators.py ref_hadspec.xml quda_hadspec.xml

This implements Level 3 validation from the QUDA SLRC integration guide:
compare pion (and other hadron) correlators computed from the BiCGStab
reference propagator and the QUDA SLRC propagator.

Expected agreement:
  - Double precision (no mixed precision): ~1e-12 or better
  - Mixed precision (single sloppy):       ~1e-6 to ~1e-8
  - With multigrid:                         ~1e-8 to ~1e-10

Exit codes:
  0 = PASS (all correlators agree within tolerance)
  1 = FAIL (significant disagreement found)
  2 = ERROR (file parsing or other error)
"""

import sys
import xml.etree.ElementTree as ET
import numpy as np


def extract_correlators(xml_file):
    """
    Extract correlator data from a Chroma HADRON_SPECTRUM XML output file.

    Returns a dict: { (hadron_name, momentum) : np.array of complex values }
    """
    try:
        tree = ET.parse(xml_file)
    except ET.ParseError as e:
        print(f"ERROR: Failed to parse {xml_file}: {e}")
        sys.exit(2)

    root = tree.getroot()
    correlators = {}

    # Navigate the Chroma HADRON_SPECTRUM output structure
    # The structure varies slightly but typically contains:
    #   <hadron_spectrum>
    #     <... meson/baryon sections ...>
    #       <momenta ...>
    #         <correlator>
    #           <re> ... </re>
    #           <im> ... </im>

    # Try multiple possible XML structures

    # Structure 1: Direct Wilson_hadron_measurements path
    for had_meas in root.iter('Wilson_hadron_measurements'):
        # Mesons
        for meson_elem in had_meas.iter('meson'):
            name = "meson"
            gamma_tag = meson_elem.find('.//gamma_value')
            if gamma_tag is not None:
                name = f"meson_g{gamma_tag.text.strip()}"

            mom_tag = meson_elem.find('.//sink_mom')
            mom_str = "0 0 0"
            if mom_tag is not None:
                mom_str = mom_tag.text.strip()

            re_tag = meson_elem.find('.//re')
            im_tag = meson_elem.find('.//im')

            if re_tag is not None and im_tag is not None:
                re_vals = [float(x) for x in re_tag.text.split()]
                im_vals = [float(x) for x in im_tag.text.split()]
                corr = np.array(re_vals) + 1j * np.array(im_vals)
                correlators[(name, mom_str)] = corr

        # Baryons
        for baryon_elem in had_meas.iter('baryon'):
            name = "baryon"
            name_tag = baryon_elem.find('.//baryon_name')
            if name_tag is not None:
                name = name_tag.text.strip()

            mom_tag = baryon_elem.find('.//sink_mom')
            mom_str = "0 0 0"
            if mom_tag is not None:
                mom_str = mom_tag.text.strip()

            re_tag = baryon_elem.find('.//re')
            im_tag = baryon_elem.find('.//im')

            if re_tag is not None and im_tag is not None:
                re_vals = [float(x) for x in re_tag.text.split()]
                im_vals = [float(x) for x in im_tag.text.split()]
                corr = np.array(re_vals) + 1j * np.array(im_vals)
                correlators[(name, mom_str)] = corr

    # Structure 2: Try iterating over all elements with re/im children
    # (fallback for different Chroma output versions)
    if not correlators:
        idx = 0
        for elem in root.iter():
            re_tag = elem.find('re')
            im_tag = elem.find('im')
            if re_tag is not None and im_tag is not None:
                try:
                    re_vals = [float(x) for x in re_tag.text.split()]
                    im_vals = [float(x) for x in im_tag.text.split()]
                    if len(re_vals) > 1:  # Skip scalar values
                        corr = np.array(re_vals) + 1j * np.array(im_vals)
                        name = elem.tag if elem.tag else f"correlator_{idx}"

                        # Try to get momentum info from parent or attributes
                        mom_str = "0 0 0"
                        mom_elem = elem.find('.//sink_mom')
                        if mom_elem is not None:
                            mom_str = mom_elem.text.strip()

                        correlators[(name, mom_str)] = corr
                        idx += 1
                except (ValueError, AttributeError):
                    continue

    return correlators


def compare_correlators(ref_corrs, test_corrs):
    """
    Compare two sets of correlators and report results.

    Returns (passed, max_relative_diff).
    """
    if not ref_corrs:
        print("ERROR: No correlators found in reference file.")
        return False, float('inf')

    if not test_corrs:
        print("ERROR: No correlators found in test file.")
        return False, float('inf')

    print(f"Reference file: {len(ref_corrs)} correlator(s)")
    print(f"Test file:      {len(test_corrs)} correlator(s)")
    print()

    # Find common correlator keys
    common_keys = set(ref_corrs.keys()) & set(test_corrs.keys())

    if not common_keys:
        # Try matching by position if keys don't match
        print("WARNING: No matching correlator keys found. Comparing by position...")
        ref_list = sorted(ref_corrs.items())
        test_list = sorted(test_corrs.items())
        n_compare = min(len(ref_list), len(test_list))
        common_pairs = [(ref_list[i], test_list[i]) for i in range(n_compare)]
    else:
        common_pairs = [((k, ref_corrs[k]), (k, test_corrs[k])) for k in sorted(common_keys)]

    overall_max_diff = 0.0
    overall_pass = True
    n_compared = 0

    for (ref_key, ref_corr), (test_key, test_corr) in common_pairs:
        if len(ref_corr) != len(test_corr):
            print(f"WARNING: Length mismatch for {ref_key}: ref={len(ref_corr)}, test={len(test_corr)}")
            continue

        n_t = len(ref_corr)
        name = ref_key if isinstance(ref_key, str) else f"{ref_key[0]} (p={ref_key[1]})"

        print(f"--- {name} (T={n_t}) ---")
        print(f"{'t':>4} {'|ref|':>16} {'|test|':>16} {'rel_diff':>14} {'status':>8}")
        print("-" * 64)

        max_diff_this = 0.0
        for t in range(n_t):
            r = ref_corr[t]
            q = test_corr[t]

            abs_r = abs(r)
            abs_q = abs(q)

            if abs_r > 0:
                rel_diff = abs(q - r) / abs_r
            else:
                rel_diff = abs(q - r)

            max_diff_this = max(max_diff_this, rel_diff)

            # Flag problematic timeslices
            if rel_diff > 1e-6:
                status = "FAIL"
            elif rel_diff > 1e-10:
                status = "WARN"
            else:
                status = "OK"

            # Only print every 4th timeslice for brevity, plus any flagged ones
            if t % 4 == 0 or status != "OK":
                print(f"{t:4d} {abs_r:16.10e} {abs_q:16.10e} {rel_diff:14.6e} {status:>8}")

        print(f"\n  Max relative difference: {max_diff_this:.3e}")

        if max_diff_this < 1e-10:
            print("  -> PASS (double precision agreement)")
        elif max_diff_this < 1e-6:
            print("  -> PASS (mixed precision agreement)")
        else:
            print("  -> FAIL (significant disagreement)")
            overall_pass = False

        overall_max_diff = max(overall_max_diff, max_diff_this)
        n_compared += 1
        print()

    return overall_pass, overall_max_diff, n_compared


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)

    ref_file = sys.argv[1]
    test_file = sys.argv[2]

    print("=" * 64)
    print("QUDA SLRC Correlator Comparison (Level 3 Validation)")
    print("=" * 64)
    print(f"Reference: {ref_file}")
    print(f"Test:      {test_file}")
    print()

    ref_corrs = extract_correlators(ref_file)
    test_corrs = extract_correlators(test_file)

    passed, max_diff, n_compared = compare_correlators(ref_corrs, test_corrs)

    print("=" * 64)
    print("SUMMARY")
    print("=" * 64)
    print(f"Correlators compared: {n_compared}")
    print(f"Max relative diff:    {max_diff:.3e}")
    print()

    if n_compared == 0:
        print("ERROR: No correlators could be compared.")
        print("Check that the XML files contain HADRON_SPECTRUM output.")
        sys.exit(2)
    elif passed:
        if max_diff < 1e-10:
            print("RESULT: PASS - Double precision agreement")
            print("  The QUDA SLRC solver matches the reference solver")
            print("  to ~1e-10 or better at all timeslices.")
        else:
            print("RESULT: PASS - Mixed precision agreement")
            print("  The QUDA SLRC solver matches the reference solver")
            print(f"  to ~{max_diff:.0e} (consistent with mixed precision rounding).")
        sys.exit(0)
    else:
        print("RESULT: FAIL - Significant disagreement detected")
        print("  Investigate the SLRC gauge/clover field decomposition.")
        print("  Common causes:")
        print("    - Thin/fat links swapped")
        print("    - Boundary conditions not applied to one gauge field")
        print("    - Clover coefficient applied twice (Chroma + QUDA)")
        print("  See docs/quda_slrc_integration_guide.md Section 5.7")
        sys.exit(1)


if __name__ == "__main__":
    main()
