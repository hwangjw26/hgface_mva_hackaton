"""Extract the unchanged native two-contributing-heterozygote AR hypotheses."""
import subprocess
from canonical.identity import allele_key, pair_key, canonical_pair
from canonical.primary import primary_genotypes
from io_utils import index, require, rows, write_tsv


def normalize(variants, stem, bcftools, fasta):
    source, dest = stem.with_suffix(".input.vcf"), stem.with_suffix(".normalized.vcf")
    chrom = lambda value: "MT" if value.removeprefix("chr") == "M" else value.removeprefix("chr")
    with source.open("x") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        for contig in sorted({chrom(row["CONTIG"]) for row in variants}):
            handle.write(f"##contig=<ID={contig}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for n, row in enumerate(variants):
            require(bool(row["REF"]) and bool(row["ALT"]) and "," not in row["ALT"], "Invalid native allele")
            handle.write(f'{chrom(row["CONTIG"])}\t{row["START"]}\tROW_{n}\t{row["REF"]}\t{row["ALT"]}\t.\t.\t.\n')
    with stem.with_suffix(".stderr").open("x") as log:
        subprocess.run([bcftools, "norm", "--check-ref", "e", "-f", fasta, "-m", "-any",
                        "-Ov", "-o", str(dest), str(source)], stdout=log, stderr=log, check=True)
    keys = {}
    for line in dest.read_text().splitlines():
        if line.startswith("#"):
            continue
        c = line.split("\t")
        n = int(c[2].removeprefix("ROW_"))
        require(n not in keys, "Normalization duplicated one native row")
        keys[n] = allele_key(c[0], c[1], c[3], c[4])
    require(set(keys) == set(range(len(variants))), "Native normalization row coverage mismatch")
    return [keys[n] for n in range(len(variants))]


def extract(normalized_variants, gene_rows):
    modes = index(gene_rows, ["preset", "ENTREZ_GENE_ID", "MOI"])
    selected = primary_genotypes(normalized_variants)
    selected.sort(key=lambda x: (["exome", "genome"].index(x[0]),
                                int(modes[(x[0], x[1].removeprefix("NCBIGene:"), "AR")]["#RANK"]),
                                x[1], x[2]))
    ar_rank = {}
    for preset in ("exome", "genome"):
        genes = sorted((g for g in gene_rows if g["preset"] == preset and g["MOI"] == "AR"),
                       key=lambda g: int(g["#RANK"]))
        require(len({g["ENTREZ_GENE_ID"] for g in genes}) == len(genes), "Duplicate native AR gene")
        ar_rank.update({(preset, g["ENTREZ_GENE_ID"]): n for n, g in enumerate(genes, 1)})
    pairs, presets = {}, []
    for preset, gene, keys in selected:
        key = pair_key(*keys)
        g = modes[(preset, gene.removeprefix("NCBIGene:"), "AR")]
        if (gene, key) not in pairs:
            a, b = canonical_pair(*keys)
            pairs[(gene, key)] = {"review_id": f"ER{len(pairs)+1:04d}", "gene_id": gene,
                                  "gene_symbol": g["GENE_SYMBOL"], "pair_key": key,
                                  "allele_1_key": a, "allele_2_key": b,
                                  "exome_native_ar_rank": "", "genome_native_ar_rank": ""}
        pair = pairs[(gene, key)]
        rank = ar_rank[(preset, g["ENTREZ_GENE_ID"])]
        pair[preset + "_native_ar_rank"] = str(rank)
        presets.append({"preset": preset, "review_id": pair["review_id"], "gene_id": gene,
                        "pair_key": key, "native_gene_rank": g["#RANK"], "native_ar_rank": rank,
                        "combined_score": g["EXOMISER_GENE_COMBINED_SCORE"],
                        "genotype_score": g["EXOMISER_GENE_VARIANT_SCORE"],
                        "phenotype_score": g["EXOMISER_GENE_PHENO_SCORE"]})
    require(bool(pairs), "No native two-heterozygote AR pairs")
    return list(pairs.values()), presets


def run(cfg, out):
    variants, genes = [], []
    for preset in ("exome", "genome"):
        native = rows(cfg[preset + "_variants"])
        keys = normalize(native, out / (preset + "_alleles"), cfg["bcftools"], cfg["reference_fasta"])
        variants.extend({"preset": preset, "allele_key": key, **row} for key, row in zip(keys, native))
        genes.extend({"preset": preset, **row} for row in rows(cfg[preset + "_genes"]))
    pairs, presets = extract(variants, genes)
    write_tsv(out / "native_pairs.tsv", pairs)
    write_tsv(out / "native_presets.tsv", presets)
    return pairs, presets
