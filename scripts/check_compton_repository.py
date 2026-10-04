"""Check the standalone source graph and separation of derivation from comparison."""
import ast
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'scripts'
TEX = ROOT / 'tex'


def check():
    entries = [p for p in TEX.rglob('*.tex') if r'\documentclass' in p.read_text()]
    assert entries == [TEX / 'compton_calculation_details.tex'], entries
    seen = set()

    def expand(path, ancestors=()):
        assert path.is_relative_to(ROOT) and path.is_file(), path
        assert path not in ancestors, ('cyclic input', path)
        seen.add(path)
        source = re.sub(r'(?<!\\)%[^\n]*', '', path.read_text())
        def include(match):
            name = match[1]
            child = (TEX / name).with_suffix('.tex').resolve()
            return expand(child, ancestors + (path,))
        return re.sub(r'\\input\{([^{}]+)\}', include, source)

    source = expand(entries[0])
    chapters = re.findall(r'\\chapter\{([^{}]+)\}', source)
    expected = ['Physical setup and principal results',
                'From the Feynman rules to the NLO integrands',
                'Integration of the amplitudes', 'Simplified analytic cross section',
                'Asymptotic expansions', 'Comparison with the published analytic total',
                'Numerical verification']
    assert chapters == expected, chapters

    # The overview must reproduce the later derivation's final formulas.
    overview = (TEX / 'overview.tex').read_text()
    final = (TEX / 'total_simplification.tex').read_text()
    def displayed_equation(text, label):
        marker = r'\label{' + label + '}'
        matches = [equation for equation in
                   re.findall(r'\\begin\{equation\}(.*?)\\end\{equation\}', text, re.S)
                   if marker in equation]
        assert len(matches) == 1, ('missing or duplicate formula', label)
        return re.sub(r'\s+', '', matches[0].replace(marker, ''))
    for name in ('total', 'born'):
        assert displayed_equation(overview, 'eq:comp-overview-' + name) == displayed_equation(
            final, 'eq:comp-simplified-' + name), ('overview formula out of sync', name)
    shared = r'\input{total_energy_product_body}'
    assert overview.count(shared) == 1, 'overview must use the shared final formula'
    assert (TEX / 'total_energy_compact_formula.tex').read_text().count(shared) == 1, \
        'later chapter must use the shared final formula'
    labels = [v for v in re.findall(r'\\label\{([^{}]+)\}', source) if not v.startswith('\\')]
    assert len(labels) == len(set(labels)), 'duplicate equation or section label'
    labels += re.findall(r'\\def\\comptonproductlabel\{([^{}]+)\}', source)
    references = re.findall(r'\\(?:eqref|ref|pageref)\{([^{}]+)\}', source)
    missing = set(references) - set(labels)
    assert not missing, ('unresolved references', missing)
    assert r'\mathsf r' not in source
    assert not re.search(r'(?<!mathcal )C_K', source), 'inconsistent soft coefficient'
    for name in ['asymptotic_soft', 'asymptotic_real', 'high_energy_constants']:
        notation = (TEX/(name+'.tex')).read_text()
        if name == 'asymptotic_real':
            # The explicit Dalitz alias and its regulator use the earlier x,y,
            # not the incident-energy X governed by this notation check.
            notation = notation.replace('X_1=x=', 'X_1=')
            notation = notation.replace(r'(xy)^\eta', r'(X_1X_2)^\eta')
        assert not re.search(r'(?<![A-Za-z\\])x(?:\b|[wztU])', notation), name
    trees = {p.stem: ast.parse(p.read_text(), filename=str(p))
             for p in SCRIPTS.glob('*.py')}
    graph = {}
    for name, tree in trees.items():
        dependencies = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                dependencies.update(n.name.split('.')[0] for n in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                dependencies.add(node.module.split('.')[0])
        graph[name] = dependencies & trees.keys()
    # Include local imports inside functions as well as top-level imports.
    forbidden = {'compton_reference_high_precision', 'compton_nlo',
                 'compton_box_analytic_equivalence_basis',
                 'compton_box_analytic_reference_check'}
    roots = ['compton_total_independent_closed', 'compton_total_energy_closed',
             'compton_box_analytic_equivalence_inverse_basis',
             'compton_real_four_matched_collect', 'compton_total_energy_collect',
             'compton_total_energy_low', 'compton_total_energy_high',
             'compton_total_energy_high_constants',
             'compton_real_equivalence_u', 'compton_ck_simplification',
             'compton_real_combined_kernel', 'compton_total_regular']
    for root in roots:
        assert root in graph, root
        visited, pending = set(), [root]
        while pending:
            name = pending.pop()
            if name in visited:
                continue
            assert name not in forbidden, (root, 'imports reference', name)
            visited.add(name)
            pending.extend(graph[name] - visited)
        print(f'{root}: {len(visited)} local modules, no reference dependency')
    print(f'PASS: sole entry, seven chapters, {len(seen)} included TeX files, '
          f'{len(labels)} unique labels and {len(trees)} Python sources')


if __name__ == '__main__':
    check()
