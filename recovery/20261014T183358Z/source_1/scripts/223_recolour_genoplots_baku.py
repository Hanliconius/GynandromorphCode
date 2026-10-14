#!/usr/bin/env python3
"""Make isolated palette-only variants of the existing collapsed genoplot script."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import sys

PALETTES = ('baku_solar', 'baku_polaire')
OLD = ('#2C9BCB', '#F28E2B', '#B61F72')


def transform(source, destination, palette):
    tree = ast.parse(source)
    replacements = dict(zip(OLD, [palette[0], palette[len(palette)//2], palette[-1]]))
    counts = {'output': 0, 'colormap': 0}

    class Change(ast.NodeTransformer):
        def visit_Assign(self, node):
            node = self.generic_visit(node)
            if any(isinstance(t, ast.Name) and t.id == 'OUTPUT_DIR' for t in node.targets):
                node.value = ast.Name(id='__baku_destination', ctx=ast.Load())
                counts['output'] += 1
            return node

        def visit_FunctionDef(self, node):
            if node.name == 'make_colormap':
                node.body = [ast.Return(value=ast.Call(
                    func=ast.Name(id='__baku_colormap', ctx=ast.Load()), args=[], keywords=[]))]
                counts['colormap'] += 1
                return node
            return self.generic_visit(node)

        def visit_Constant(self, node):
            if isinstance(node.value, str) and node.value.upper() in replacements:
                return ast.copy_location(ast.Constant(replacements[node.value.upper()]), node)
            return node

    tree = Change().visit(tree)
    ast.fix_missing_locations(tree)
    if counts != {'output': 1, 'colormap': 1}:
        raise RuntimeError(f'Expected one OUTPUT_DIR and make_colormap; found {counts}. Original script unchanged.')
    return tree


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--palette', choices=PALETTES, required=True)
    p.add_argument('--source', type=Path, default=Path('scripts/212_plot_collapsed_Krista_loci.py'))
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.colors import LinearSegmentedColormap
    project = Path.cwd()
    source = args.source.resolve()
    recipe_path = Path(__file__).resolve().parents[1] / 'palettes' / (args.palette + '.json')
    recipe = json.loads(recipe_path.read_text())
    colours = recipe['colours']
    if len(colours) != 200:
        raise RuntimeError('Expected the complete 200-colour palette')
    destination = project / 'analysis/figures' / ('genoplots_' + source.stem + '_' + args.palette)
    code = source.read_text()
    tree = transform(code, destination, colours)
    compiled = compile(tree, str(source), 'exec')

    def make_colormap():
        cmap = LinearSegmentedColormap.from_list(args.palette, colours, N=256)
        cmap.set_bad('#F2F2F2')
        return cmap

    cmap = make_colormap()
    if not np.isfinite(cmap(np.linspace(0, 1, 201))).all():
        raise RuntimeError('Invalid colormap')
    print(f'PALETTE\t{args.palette}\tblue_to={colours[0]}\tmidpoint={colours[100]}\tend={colours[-1]}')
    print(f'OUTPUT_DIR\t{destination}')
    if destination.exists():
        raise RuntimeError(f'Output directory already exists; preserved: {destination}')
    if args.check:
        print('PREFLIGHT_PASS\tsource syntax, palette hooks, dependencies; data execution deferred to job')
        return
    destination.mkdir(parents=True)
    provenance = {'source': str(source), 'source_sha256': hashlib.sha256(code.encode()).hexdigest(),
                  'palette': recipe, 'output_directory': str(destination),
                  'changes': ['OUTPUT_DIR', 'make_colormap', 'three old literal genotype colours'],
                  'status': 'RUNNING'}
    marker = destination / 'palette_variant.audit.json'
    marker.write_text(json.dumps(provenance, indent=2) + '\n')
    namespace = {'__name__': 'baku_genoplot_source', '__file__': str(source),
                 '__baku_destination': destination, '__baku_colormap': make_colormap}
    saved_argv = sys.argv
    sys.argv = [str(source)]
    try:
        exec(compiled, namespace)
        # Guard known output constants before invoking main.
        for name, value in namespace.items():
            if isinstance(value, Path) and (name.endswith('_PDF') or name.endswith('_PNG') or name.endswith('_TSV') or name.endswith('_JSON')):
                if not value.resolve().is_relative_to(destination.resolve()) and 'analysis/figures' in str(value):
                    raise RuntimeError(f'Output {name} escaped new directory: {value}')
        if not callable(namespace.get('main')):
            raise RuntimeError('Source has no callable main()')
        namespace['main']()
        pdfs = sorted(destination.rglob('*.pdf'))
        if not pdfs:
            raise RuntimeError('Source produced no PDFs in isolated destination')
        provenance.update(status='COMPLETE', pdfs=[str(p) for p in pdfs])
        print('COMPLETE')
        for path in pdfs:
            print(path)
    except BaseException:
        provenance['status'] = 'FAILED'
        raise
    finally:
        sys.argv = saved_argv
        marker.write_text(json.dumps(provenance, indent=2) + '\n')


if __name__ == '__main__':
    main()
