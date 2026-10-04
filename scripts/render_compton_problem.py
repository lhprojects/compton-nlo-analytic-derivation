"""Render the README reaction diagram from TikZ to PNG using pdfLaTeX and Poppler."""
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    for command in ('pdflatex', 'pdftoppm'):
        if shutil.which(command) is None:
            raise SystemExit(f'{command} is required to render the README diagram')
    target = ROOT / 'build' / 'readme-diagram'
    target.mkdir(parents=True, exist_ok=True)
    source = target / 'diagram.tex'
    source.write_text(
        r'\documentclass[tikz,border=8pt]{standalone}'+'\n'
        r'\usetikzlibrary{arrows.meta}'+'\n'
        r'\begin{document}'+'\n'
        r'\input{'+(ROOT/'tex/figures/compton_problem.tex').as_posix()+'}\n'
        r'\end{document}'+'\n')
    subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                    source.name], cwd=target, check=True, stdout=subprocess.PIPE,
                   stderr=subprocess.STDOUT, text=True)
    output = ROOT / 'tex' / 'figures' / 'compton_problem'
    subprocess.run(['pdftoppm', '-png', '-singlefile', '-r', '180',
                    str(target/'diagram.pdf'), str(output)], check=True)
    print(output.with_suffix('.png'))


if __name__ == '__main__':
    main()
