"""Compare independently assembled coefficients with the reference, exactly.

The reference is an output-side comparison only. No independent coefficient
or integration constant is set by this script.
"""
from compton_artifacts import read_artifact, read_artifact_bytes, write_artifact
import hashlib
import json
from pathlib import Path
import sympy as S
from compton_box_analytic_equivalence_basis import reference_x_words, x, Y

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT

def main():
    own_path = FOLDER / 'total_energy_collected.json'
    source = json.loads(read_artifact(own_path))
    own = {tuple(row['word']): S.sympify(row['coefficient'], locals={'X': x, 'Y': Y})
           for row in source['terms']}
    reference = {tuple(map(str, word)): coefficient
                 for word, coefficient in reference_x_words().items()}
    rows = []
    for word in sorted(own.keys() | reference.keys(), key=lambda w: (len(w), w)):
        difference = S.factor(S.cancel(own.get(word, 0) - reference.get(word, 0)))
        rows.append({'word': list(word), 'difference': str(difference)})
    failures = [row for row in rows if row['difference'] != '0']
    report = {
        'reference': 'Lee, Schwartz and Zhang, arXiv:2102.06718, Eq. (6)',
        'normalization': 'sigma_1 = alpha_T^3/m^2 C1(X), X=s_physical/m^2-1',
        'method': 'Exact shuffle expansion followed by rational coefficient subtraction in X and Y.',
        'own_sha256': hashlib.sha256(read_artifact_bytes(own_path)).hexdigest(),
        'reference_transcription_sha256': hashlib.sha256(read_artifact_bytes(ROOT / 'scripts/compton_reference_high_precision.py')).hexdigest(),
        'own_word_count': len(own), 'reference_word_count': len(reference),
        'differences': rows, 'nonzero_differences': len(failures),
        'status': 'passed' if not failures else 'failed',
        'scope': 'Exact coefficient identity. Physical branches and zero-endpoint normalization are established in analytic_equivalence.tex.',
        'transcription_audit': 'All 21 groups in Eq. (6) and the weights in Eqs. (7)-(8) checked against arXiv:2102.06718v2 on 2026-10-04.'
    }
    write_artifact(FOLDER / 'total_reference_identity.json' , json.dumps(report, indent=2) + '\n')
    print(f'{len(rows)} coefficient differences, {len(failures)} nonzero', flush=True)
    assert not failures, failures

if __name__ == '__main__':
    main()
