"""Run the standard reproducibility suite with per-command logs and no silent skips.

All preparations are analytic. Numerical checks consume the resulting formulas
and never choose their coefficients, integration constants or physical branches.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
# Regenerate disposable inputs needed by the checks, in dependency order.
PREPARE = [
    'compton_real_equivalence_allfour.py',
    'compton_real_allfour_endpoint.py',
    'compton_matching_simplify.py',
    'compton_real_four_finite_collect.py',
    'compton_real_four_finite_shorten.py',
    'compton_real_four_matched_collect.py',
    'compton_total_energy_high.py',
    'compton_total_energy_high_constants.py',
    'compton_total_energy_low.py',
]
CHECKS = [
    'check_compton_repository.py',
    'check_compton_born_analytic.py',
    'compton_dilog_identities.py',
    'compton_real_equivalence_u.py',
    'compton_ck_simplification.py',
    'compton_real_combined_kernel.py',
    'compton_total_regular.py',
    'check_compton_generated_sources.py',
    'check_compton_amplitudes.py',
    'check_compton_scalar_analytic.py',
    'check_compton_vertex_zero.py',
    'check_compton_vertex_analytic.py',
    'check_compton_virtual_trace_analytic.py',
    'check_compton_virtual_angular_analytic.py',
    'check_compton_infrared_analytic.py',
    'check_compton_real_soft_matching.py',
    'check_compton_real_integrated_ir.py',
    'compton_box_analytic_check.py',
    'compton_box_analytic_angle_check.py',
    'compton_box_analytic_gpl_check.py',
    'compton_box_analytic_bose_check.py',
    'compton_box_analytic_bose_crosscheck.py',
    'compton_box_analytic_matching_check.py',
    'compton_real_analytic_check.py',
    'compton_real_analytic_angular_check.py',
    'compton_real_analytic_gpl_check.py',
    'compton_real_analytic_rational_check.py',
    'compton_real_analytic_rational_fieldcheck.py',
    'compton_real_allfour_endpoint_check.py',
    'compton_real_qpair_endpoint_check.py',
    'compton_real_four_finite_check.py',
    'compton_real_four_matched_check.py',
    'compton_ck_simplification_check.py',
    'compton_ck_classical_check.py',
    'compton_ck_energy_check.py',
    'compton_soft_uv_pole_certificate.py',
    'compton_soft_uv_uw_origin_certificate.py',
    'compton_soft_uv_uw_endpoint_certificate.py',
    'compton_uv_primitive_certificate.py',
    'compton_total_energy_low_check.py',
    'compton_total_energy_high_check.py',
    'check_compton_dilog_identities.py',
    'check_compton_total_regular.py',
    'check_compton_total_analytic.py',
    'compton_total_reference_identity.py',
    'check_compton_nlo.py',
    'check_compton_virtual.py',
]

# Some checks consume another check's independently evaluated result.
# Always regenerate that input, including when selecting checks with --only.
CHECK_DEPENDENCIES = {
    'compton_real_four_matched_check.py': ('compton_real_four_finite_check.py',),
}


def resolve_checks(requested):
    result, active, done = [], set(), set()

    def visit(name):
        if name in done:
            return
        if name not in CHECKS:
            raise ValueError(f'unknown check: {name}')
        if name in active:
            raise ValueError(f'cyclic check dependency: {name}')
        active.add(name)
        for dependency in CHECK_DEPENDENCIES.get(name, ()):
            visit(dependency)
        active.remove(name)
        done.add(name)
        result.append(name)

    for name in requested:
        visit(name)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--list', action='store_true', help='list commands without running them')
    parser.add_argument('--only', nargs='+', help='run named checks with their preparation and prerequisite checks')
    args = parser.parse_args()
    checks = CHECKS if not args.only else args.only
    if any(name not in CHECKS for name in checks):
        parser.error('--only must name a check from --list')
    checks = resolve_checks(checks)
    commands = [('prepare', name) for name in PREPARE] + [('check', name) for name in checks]
    if args.list:
        for stage, name in commands:
            print(stage, name)
        print(f'{len(PREPARE)} preparation steps, {len(checks)} checks')
        return
    out = ROOT / 'build' / 'checks'
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for stage, name in commands:
        print(f'{stage.upper()} {name}', flush=True)
        start = time.monotonic()
        log = out / (Path(name).stem + '.log')
        with log.open('w') as handle:
            result = subprocess.run([sys.executable, str(ROOT/'scripts'/name)],
                                    cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)
        record = dict(stage=stage, script=name, returncode=result.returncode,
                      seconds=round(time.monotonic()-start, 3), log=str(log.relative_to(ROOT)))
        records.append(record)
        (out/'summary.json').write_text(json.dumps(records, indent=2)+'\n')
        print(f"{'PASS' if result.returncode == 0 else 'FAIL'} {name} ({record['seconds']} s)", flush=True)
        if result.returncode:
            print(log.read_text()[-6000:], flush=True)
            # A failed preparation invalidates all dependent checks.
            if stage == 'prepare':
                raise SystemExit(1)
    failures = [record for record in records if record['returncode']]
    print(f'{len(checks)-len(failures)}/{len(checks)} checks passed; logs: {out}', flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
