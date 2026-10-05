package org.monarchinitiative.lirical.cli.cmd;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.monarchinitiative.lirical.core.model.*;
import org.monarchinitiative.svart.CoordinateSystem;
import org.monarchinitiative.svart.Strand;
import picocli.CommandLine;
import java.nio.file.*;
import java.util.*;

/** Observe the entire canonical stream; never subset it for disease scoring. */
public class FullMetadata extends YamlCommand {
    static final ObjectMapper JSON = new ObjectMapper();
    static List<String> strings(JsonNode n) {
        List<String> out = new ArrayList<>(); n.forEach(x -> out.add(x.asText())); return out;
    }
    static String key(LiricalVariant v) {
        var g = v.variant();
        return g.contigName().replaceFirst("^chr", "") + ":"
            + g.startOnStrandWithCoordinateSystem(Strand.POSITIVE, CoordinateSystem.oneBased())
            + ":" + g.ref() + ":" + g.alt();
    }
    public static void main(String[] args) throws Exception {
        JsonNode input = JSON.readTree(Path.of(args[0]).toFile());
        FullMetadata command = new FullMetadata();
        List<String> cli = strings(input.get("cli_args"));
        if (cli.contains("--sdwndv")) throw new IllegalStateException("Inverted flag must be omitted");
        new CommandLine(command).parseArgs(cli.toArray(String[]::new));
        command.dataSection.parallelism = 4;
        var lirical = command.bootstrapLirical(GenomeBuild.HG38);
        var options = command.prepareAnalysisOptions(lirical, GenomeBuild.HG38, TranscriptDatabase.ENSEMBL);
        if (!options.includeDiseasesWithNoDeleteriousVariants() || options.targetDiseases() != null)
            throw new IllegalStateException("Expected unrestricted native configuration");
        Map<String, Set<String>> wanted = new TreeMap<>();
        for (var c : input.get("contexts")) wanted.computeIfAbsent(c.get("key").asText(), k -> new TreeSet<>()).add(c.get("gene").asText());
        Set<String> scoped = new HashSet<>(strings(input.get("task19_scope_genes")));
        String sample = input.get("sample").asText();
        Map<String, long[]> geneCounts = new TreeMap<>();
        List<Map<String,Object>> observed = new ArrayList<>(), contexts = new ArrayList<>();
        long[] counts = new long[2];
        try (var parser = lirical.variantParserFactory().forPath(Path.of(input.get("vcf").asText()), GenomeBuild.HG38, TranscriptDatabase.ENSEMBL).orElseThrow()) {
            if (!parser.sampleNames().contains(sample)) throw new IllegalStateException("Sample identity mismatch");
            parser.variantStream().sequential().forEach(v -> {
                counts[0]++;
                var identifiers = new TreeMap<String, org.monarchinitiative.phenol.annotations.formats.GeneIdentifier>();
                v.annotations().forEach(a -> identifiers.put(a.getGeneId().id().getValue(), a.getGeneId()));
                for (String g : identifiers.keySet()) {
                    long[] c = geneCounts.computeIfAbsent(g, unused -> new long[3]);
                    c[0]++; c[v.passedFilters() ? 1 : 2]++;
                }
                if (identifiers.keySet().stream().anyMatch(scoped::contains)) counts[1]++;
                String k = key(v);
                if (!wanted.containsKey(k)) return;
                Map<String,Object> row = new TreeMap<>();
                row.put("key", k); row.put("genes", new ArrayList<>(identifiers.keySet()));
                row.put("passed", v.passedFilters()); row.put("score", v.pathogenicityScore().orElse(null));
                row.put("clinvar", v.clinvarClnSig().name());
                row.put("pathogenic_clinvar_count", v.pathogenicClinVarAlleleCount(sample));
                row.put("ref_count", v.alleleCount(sample).map(a -> (int)a.ref()).orElse(null));
                row.put("alt_count", v.alleleCount(sample).map(a -> (int)a.alt()).orElse(null));
                observed.add(row);
                for (String g : wanted.get(k)) {
                    Map<String,Object> context = new TreeMap<>();
                    context.put("key", k); context.put("gene", g);
                    context.put("mapped", identifiers.containsKey(g));
                    context.put("score_eligible_allele_count", null); context.put("plp_allele_count", null);
                    if (identifiers.containsKey(g)) {
                        var single = Gene2Genotype.of(identifiers.get(g), v.passedFilters() ? List.of(v) : List.of(), v.passedFilters() ? 0 : 1);
                        context.put("score_eligible_allele_count", single.deleteriousAlleleCount(sample, options.variantDeleteriousnessThreshold()));
                        context.put("plp_allele_count", single.pathogenicClinVarCount(sample));
                    }
                    contexts.add(context);
                }
            });
        }
        observed.sort(Comparator.comparing(r -> (String)r.get("key")));
        contexts.sort(Comparator.comparing(r -> r.get("gene") + "|" + r.get("key")));
        Map<String,Object> report = new TreeMap<>();
        report.put("parsed_alleles", counts[0]); report.put("task19_scope_alleles", counts[1]);
        report.put("all_gene_annotation_counts", geneCounts); report.put("allele_observations", observed);
        report.put("contexts", contexts); report.put("threshold", options.variantDeleteriousnessThreshold());
        report.put("gene_retention_filter_used", false);
        JSON.writerWithDefaultPrettyPrinter().writeValue(Path.of(args[1]).toFile(), report);
    }
}
