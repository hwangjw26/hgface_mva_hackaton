"""Prepare portable Exomiser/LIRICAL commands; never launches genomic analysis."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex


def prepare(config, out):
    required = ['java', 'exomiser_jar', 'lirical_jar', 'exomiser_data', 'remm',
                'lirical_data', 'vcf', 'sample_id', 'hpo_ids']
    if not all(k in config for k in required):
        raise ValueError('Missing native configuration fields')
    hpo = config['hpo_ids']
    if not hpo or len(hpo) != len(set(hpo)) or not all(re.fullmatch(r'HP:[0-9]{7}', h) for h in hpo):
        raise ValueError('HPO terms must be unique, nonempty and valid identifiers')
    if any('\n' in str(config[k]) or '\r' in str(config[k]) for k in required if k != 'hpo_ids'):
        raise ValueError('Configuration paths and sample identifier must be single-line')
    template = Path(__file__).resolve().parents[1] / 'config'
    out.mkdir(parents=True, exist_ok=False)
    os.chmod(out, 0o700)
    def save(name, value):
        p = out / name
        p.write_text(json.dumps(value, indent=2) + '\n')
        return str(p.resolve())
    sample = save('sample.json', {
        'id': config['sample_id'], 'subject': {'id': config['sample_id']},
        'phenotypicFeatures': [{'type': {'id': h}} for h in hpo],
        'metaData': {'created': '2026-10-04T00:00:00Z', 'createdBy': 'reproducible-workflow', 'phenopacketSchemaVersion': '2.0'}})
    properties = out / 'application.properties'
    properties.write_text('exomiser.data-directory=' + str(Path(config['exomiser_data']).resolve()) + '\n'
                          'exomiser.hg38.data-version=2602\nexomiser.phenotype.data-version=2602\n'
                          'exomiser.hg38.remm-path=' + str(Path(config['remm']).resolve()) + '\n'
                          'spring.cache.type=none\nlogging.level.com.zaxxer.hikari=ERROR\n')
    commands = []
    java, vcf = str(Path(config['java']).resolve()), str(Path(config['vcf']).resolve())
    for preset in ['exome', 'genome']:
        commands.append([java, '-Xmx24g', '-XX:ActiveProcessorCount=8',
                         '-Dspring.config.location=file:' + str(properties.resolve()),
                         '-jar', str(Path(config['exomiser_jar']).resolve()), 'analyse', '--assembly', 'GRCh38',
                         '--vcf', vcf, '--sample', sample,
                         '--analysis', str(template / f'preset-{preset}-analysis.yml'),
                         '--output', str(template / 'output-options.json'),
                         '--output-directory', str(out.resolve() / preset), '--output-filename', f'singleton_{preset}'])
    for i in range(len(hpo) + 1):
        case = 'baseline' if i == 0 else f'loo_{i:02d}'
        inp = save(case + '.yaml', {'sampleId': config['sample_id'], 'hpoIds': hpo if i == 0 else hpo[:i-1]+hpo[i:],
                                  'negatedHpoIds': [], 'sex': 'UNKNOWN', 'vcf': vcf})
        commands.append([java, '-XX:-UsePerfData', '-XX:ActiveProcessorCount=8', '-Xmx16g',
                         '-jar', str(Path(config['lirical_jar']).resolve()), 'yaml', '-y', inp, '--assembly', 'hg38',
                         '-ed38', str(Path(config['exomiser_data']).resolve() / '2602_hg38'),
                         '-d', str(Path(config['lirical_data']).resolve()), '--transcript-db', 'ENSEMBL',
                         '--parallelism', '8', '--pathogenicity-threshold', '0.8', '--validation-policy', 'STRICT',
                         '-f', 'json', '-f', 'tsv', '-f', 'html', '--display-all-variants',
                         '-o', str(out.resolve() / case), '-x', case])
    save('commands.json', commands)
    lines = ['#!/usr/bin/env bash', 'set -euo pipefail', 'umask 077',
             ': "${TMPDIR:?Set TMPDIR to an existing local scratch directory}"', 'test -d "$TMPDIR"']
    for i, cmd in enumerate(commands):
        output_flag = '--output-directory' if '--output-directory' in cmd else '-o'
        paths = [Path(cmd[cmd.index(output_flag) + 1]), out.resolve() / f'native_{i:02d}.log']
        for path in paths:
            lines.append('if test -e ' + shlex.quote(str(path)) + '; then echo "Refusing existing native output" >&2; exit 1; fi')
    lines.append('set -o noclobber')
    for i, cmd in enumerate(commands):
        # Insert TMPDIR as a separately quoted shell argument, never as user code.
        line = shlex.quote(cmd[0]) + ' "-Djava.io.tmpdir=$TMPDIR" ' + shlex.join(cmd[1:])
        lines.append(line + ' >' + shlex.quote(str(out.resolve() / f'native_{i:02d}.log')) + ' 2>&1')
    (out / 'run_native.sh').write_text('\n'.join(lines) + '\n')
    return commands


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    os.umask(0o077)
    commands = prepare(json.loads(args.config.read_text()), args.output)
    print(f'Prepared {len(commands)} commands; no native analysis was executed.')
