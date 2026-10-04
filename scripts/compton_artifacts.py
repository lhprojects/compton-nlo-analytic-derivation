"""Locations of analytic inputs, final reports and disposable intermediates."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYTIC_FILES = frozenset({
    'ck_energy_collected.json',
    'compton_box_analytic_angle.json',
    'compton_box_analytic_bose_gpl.json',
    'compton_box_analytic_bose_kernels.json',
    'compton_box_analytic_trace.json',
    'compton_real_analytic_UW_endpoints.json',
    'compton_real_analytic_UW_finite.json',
    'compton_real_analytic_angles.json',
    'compton_real_analytic_logs_gpl.json',
    'compton_real_analytic_logs_inner.json',
    'compton_real_analytic_rational_final_kernels.json',
    'compton_real_analytic_rational_gpl.json',
    'compton_real_analytic_single_gpl.json',
    'compton_real_analytic_single_sector.json',
    'compton_real_analytic_trace.json',
    'compton_real_bose_subtraction.json',
    'compton_real_soft_matching_integrated.json',
    'compton_real_uv_analytic.json',
    'compton_virtual_trace_analytic.json',
    'real_equivalence_combined.json',
    'real_four_matched_coefficients.json',
    'real_four_matched_short.json',
    'total_energy_collected.json',
    'total_energy_compact.json',
    'total_simplification_virtual_uv.json',
    'uw_inverse_basis.json',
})
REPORT_FILES = frozenset({
    'amplitude_checks.json',
    'ck_classical_full_check.json',
    'ck_energy_closed_check.json',
    'ck_energy_masters_check.json',
    'ck_hermite_check.json',
    'ck_rational_weight3_check.json',
    'ck_simplification_check.json',
    'comparison_checks.json',
    'compton_box_analytic_reference_check.json',
    'compton_real_analytic_numerical_audit.json',
    'compton_real_analytic_numerical_wide_audit.json',
    'compton_soft_uv_numeric_virtual.json',
    'experimental_curve_analytic.json',
    'inclusive_checks.json',
    'independent_closed_reference_check.json',
    'independent_closed_s10.json',
    'independent_closed_s100_55.json',
    'independent_closed_s100_75.json',
    'independent_closed_s10_70.json',
    'independent_closed_s1p0002_70.json',
    'independent_closed_s1p0002_90.json',
    'independent_closed_s1p02_70.json',
    'independent_closed_s1p02_90.json',
    'independent_closed_s1p6_50.json',
    'independent_closed_s1p6_70.json',
    'independent_closed_s200001_70.json',
    'independent_closed_s200001_90.json',
    'independent_closed_s2001_70.json',
    'independent_closed_s2001_90.json',
    'independent_total_table_check.json',
    'real_allfour_endpoint_check.json',
    'real_four_finite_check.json',
    'real_four_matched_check.json',
    'real_integrated_ir_check.json',
    'real_qpair_endpoint_check.json',
    'real_soft_matching_check.json',
    'root_classical_check.json',
    'soft_pole_certificate.json',
    'total_compact_table_values.json',
    'total_energy_closed_check.json',
    'total_energy_high_check.json',
    'total_energy_low_check.json',
    'total_reference_identity.json',
    'uv_primitive_certificate_finite.json',
    'uv_primitive_certificate_pole.json',
    'uw_endpoint_certificate.json',
    'uw_energy_full_check.json',
    'uw_origin_certificate.json',
    'uw_simplified_kernel_check.json',
    'vertex_zero_check.json',
    'virtual_checks.json',
})


def artifact_path(path):
    """Resolve legacy root/script paths while preserving explicit external paths."""
    path = Path(path)
    resolved = path.resolve()
    if not resolved.is_relative_to(ROOT) or path.suffix not in {'.json', '.txt'}:
        return path
    if path.name == 'requirements.txt':
        return path
    if path.name in ANALYTIC_FILES:
        return ROOT / 'data' / 'analytic' / path.name
    if path.name in REPORT_FILES or path.name.startswith('independent_closed_s'):
        return ROOT / 'reports' / path.name
    return ROOT / 'build' / 'intermediate' / path.name


def read_artifact(path, *args, **kwargs):
    target = artifact_path(path)
    if not target.exists() and target.parent == ROOT / 'build' / 'intermediate':
        raise FileNotFoundError(
            f'{target.name} is a disposable intermediate. Regenerate it with '
            'the corresponding derivation script before this step. '
            f'Expected location: {target}')
    return target.read_text(*args, **kwargs)


def read_artifact_bytes(path):
    return artifact_path(path).read_bytes()


def write_artifact(path, content, *args, **kwargs):
    target = artifact_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    return target.write_text(content, *args, **kwargs)
