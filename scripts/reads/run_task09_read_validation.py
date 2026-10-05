"""Thin wrapper: run the unchanged Task06 targeted-read orchestrator on the Task09 frozen targets.

Sets the expected pair/allele cardinality from the frozen target handoff. If data/sensitive/probe_excluded_alleles.tsv
exists, the listed alleles' probes are removed after probe design (attempt02 smoke: low-complexity probes exceeded the
frozen smoke gate). Those alleles keep probe_ok=False, so Task06 rules label them READ_NONCALL, never contradicted.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_targeted_fastq_validation as task06  # noqa: E402  byte-identical Task06 snapshot


def key(allele):
    c = str(allele['chrom']).removeprefix('chr')
    return f"{'MT' if c == 'M' else c}:{int(allele['pos1'])}:{allele['ref']}:{allele['alt']}"


def main():
    handoff = Path(sys.argv[sys.argv.index('--candidate-handoff') + 1])
    with handoff.open() as s:
        table = list(csv.DictReader(s, delimiter='\t'))
    task06.EXPECTED_PAIR_COUNT = len(table)
    task06.EXPECTED_ALLELE_COUNT = len({r[f'allele_{i}_normalized_allele_key'] for r in table for i in (1, 2)})
    exclusion = handoff.parent / 'probe_excluded_alleles.tsv'
    if exclusion.exists():
        with exclusion.open() as s:
            excluded = {r['allele_key'] for r in csv.DictReader(s, delimiter='\t')}
        original = task06.build_probes

        def build_probes(alleles, reference, tools, sensitive, runner, probe_contract):
            probes = original(alleles, reference, tools, sensitive, runner, probe_contract)
            drop = {str(a['allele_id']) for a in alleles if key(a) in excluded}
            assert len(drop) == len(excluded), 'exclusion list does not match frozen targets'
            kept = [p for p in probes if str(p['allele_id']) not in drop]
            for a in alleles:
                if str(a['allele_id']) in drop:
                    a.update(left_probe_count=0, right_probe_count=0, unique_probe_count=0, probe_ok=False, mapping_caveat=True)
            task06.write_tsv(sensitive / 'probes' / 'probe_manifest.tsv',
                             ['probe_id', 'allele_id', 'reference_chrom', 'side', 'start0', 'end0', 'orientation', 'sequence',
                              'exact_hit_count', 'maximum_mapq', 'is_unique'], kept)
            fasta = sensitive / 'probes' / 'allele_neutral_probes.fa'
            with fasta.open('w', encoding='ascii') as handle:
                for index, sequence in enumerate(sorted({str(p['sequence']) for p in kept}), start=1):
                    handle.write(f'>T{index:05d}\n{sequence}\n')
            fasta.chmod(0o600)
            return kept
        task06.build_probes = build_probes
    return task06.main()


if __name__ == '__main__':
    raise SystemExit(main())
