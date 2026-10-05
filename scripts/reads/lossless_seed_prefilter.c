#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define KMER_LENGTH 15
#define TABLE_SIZE 131072u
#define TABLE_MASK (TABLE_SIZE - 1u)
#define KMER_MASK ((1u << (2 * KMER_LENGTH)) - 1u)
#define MAX_PROBES 2048
#define PROBE_LENGTH 51
#define MAX_SEED_NODES (MAX_PROBES * (PROBE_LENGTH - KMER_LENGTH + 1))
#define MAX_LINE 1048576
#define MAX_PROBE 4096

static uint32_t keys[TABLE_SIZE];
static unsigned char occupied[TABLE_SIZE];
static int32_t buckets[TABLE_SIZE];
static char probes[MAX_PROBES][PROBE_LENGTH + 1];
static int probe_count = 0;
static struct {
    uint32_t key;
    int32_t next;
    uint16_t probe_index;
    uint8_t probe_offset;
} seed_nodes[MAX_SEED_NODES];
static int seed_node_count = 0;
static char r1_lines[4][MAX_LINE];
static char r2_lines[4][MAX_LINE];

static uint32_t hash_key(uint32_t value) {
    value ^= value >> 16;
    value *= 0x7feb352du;
    value ^= value >> 15;
    value *= 0x846ca68bu;
    value ^= value >> 16;
    return value;
}

static int base_code(char base) {
    switch (toupper((unsigned char)base)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

static int insert_key(uint32_t key) {
    uint32_t index = hash_key(key) & TABLE_MASK;
    while (occupied[index]) {
        if (keys[index] == key) return 0;
        index = (index + 1u) & TABLE_MASK;
    }
    occupied[index] = 1;
    keys[index] = key;
    return 1;
}

static int add_seed_node(uint32_t key, int probe_index, int probe_offset) {
    if (seed_node_count >= MAX_SEED_NODES) return 0;
    uint32_t bucket = hash_key(key) & TABLE_MASK;
    seed_nodes[seed_node_count].key = key;
    seed_nodes[seed_node_count].next = buckets[bucket];
    seed_nodes[seed_node_count].probe_index = (uint16_t)probe_index;
    seed_nodes[seed_node_count].probe_offset = (uint8_t)probe_offset;
    buckets[bucket] = seed_node_count;
    seed_node_count++;
    return 1;
}

static size_t sequence_length(const char *line) {
    size_t length = strlen(line);
    while (length > 0 && (line[length - 1] == '\n' || line[length - 1] == '\r')) length--;
    return length;
}

static uint64_t add_probe_kmers(const char *sequence) {
    size_t length = sequence_length(sequence);
    if (length != PROBE_LENGTH || probe_count >= MAX_PROBES) return 0;
    memcpy(probes[probe_count], sequence, length);
    probes[probe_count][length] = '\0';
    uint32_t code = 0;
    int valid = 0;
    uint64_t added = 0;
    for (size_t i = 0; i < length; i++) {
        int value = base_code(sequence[i]);
        if (value < 0) {
            valid = 0;
            code = 0;
            continue;
        }
        code = ((code << 2) | (uint32_t)value) & KMER_MASK;
        valid++;
        if (valid >= KMER_LENGTH) {
            int offset = (int)i - KMER_LENGTH + 1;
            if (!add_seed_node(code, probe_count, offset)) return 0;
            added += (uint64_t)insert_key(code);
        }
    }
    probe_count++;
    return added;
}

static uint64_t load_probe_kmers(const char *path) {
    for (uint32_t i = 0; i < TABLE_SIZE; i++) buckets[i] = -1;
    FILE *handle = fopen(path, "r");
    if (handle == NULL) return 0;
    char line[MAX_PROBE];
    char sequence[MAX_PROBE] = "";
    size_t used = 0;
    uint64_t added = 0;
    int saw_record = 0;
    while (fgets(line, sizeof(line), handle) != NULL) {
        if (line[0] == '>') {
            if (used > 0) {
                sequence[used] = '\0';
                added += add_probe_kmers(sequence);
            }
            used = 0;
            saw_record = 1;
            continue;
        }
        size_t length = sequence_length(line);
        if (used + length >= sizeof(sequence)) {
            fclose(handle);
            return 0;
        }
        memcpy(sequence + used, line, length);
        used += length;
    }
    if (used > 0) {
        sequence[used] = '\0';
        added += add_probe_kmers(sequence);
    }
    fclose(handle);
    return saw_record ? added : 0;
}

static int alignment_matches(const char *read, size_t read_length, const char *probe, int relative_start) {
    int read_start = relative_start > 0 ? relative_start : 0;
    int probe_start = relative_start < 0 ? -relative_start : 0;
    int read_available = (int)read_length - read_start;
    int probe_available = PROBE_LENGTH - probe_start;
    int overlap = read_available < probe_available ? read_available : probe_available;
    if (overlap < 31) return 0;
    int maximum_errors = (overlap * 4) / 100;
    int errors = 0;
    for (int offset = 0; offset < overlap; offset++) {
        if (toupper((unsigned char)read[read_start + offset]) != toupper((unsigned char)probe[probe_start + offset])) {
            errors++;
            if (errors > maximum_errors) return 0;
        }
    }
    return 1;
}

static int sequence_has_verified_match(const char *sequence) {
    size_t length = sequence_length(sequence);
    uint32_t code = 0;
    int valid = 0;
    for (size_t i = 0; i < length; i++) {
        int value = base_code(sequence[i]);
        if (value < 0) {
            valid = 0;
            code = 0;
            continue;
        }
        code = ((code << 2) | (uint32_t)value) & KMER_MASK;
        valid++;
        if (valid >= KMER_LENGTH) {
            int read_offset = (int)i - KMER_LENGTH + 1;
            uint32_t bucket = hash_key(code) & TABLE_MASK;
            for (int32_t node = buckets[bucket]; node >= 0; node = seed_nodes[node].next) {
                if (seed_nodes[node].key != code) continue;
                int relative_start = read_offset - (int)seed_nodes[node].probe_offset;
                if (alignment_matches(sequence, length, probes[seed_nodes[node].probe_index], relative_start)) return 1;
            }
        }
    }
    return 0;
}

static size_t mate_core_length(const char *header) {
    size_t length = 0;
    const char *start = header[0] == '@' ? header + 1 : header;
    while (start[length] != '\0' && start[length] != ' ' && start[length] != '\t' && start[length] != '\r' && start[length] != '\n') length++;
    if (length >= 2 && start[length - 2] == '/' && (start[length - 1] == '1' || start[length - 1] == '2')) length -= 2;
    return length;
}

static int mate_names_match(const char *header1, const char *header2) {
    const char *start1 = header1[0] == '@' ? header1 + 1 : header1;
    const char *start2 = header2[0] == '@' ? header2 + 1 : header2;
    size_t length1 = mate_core_length(header1);
    size_t length2 = mate_core_length(header2);
    return length1 == length2 && strncmp(start1, start2, length1) == 0;
}

static int valid_record(char lines[4][MAX_LINE]) {
    return lines[0][0] == '@' && lines[2][0] == '+' && sequence_length(lines[1]) == sequence_length(lines[3]);
}

int main(int argc, char **argv) {
    if (argc != 7) {
        fprintf(stderr, "usage: lossless_seed_prefilter PROBES R1 R2 OUT1 OUT2 STATS\n");
        return 2;
    }
    uint64_t seed_count = load_probe_kmers(argv[1]);
    if (seed_count == 0) {
        fprintf(stderr, "probe seed construction failed\n");
        return 3;
    }
    gzFile input1 = gzopen(argv[2], "rb");
    gzFile input2 = gzopen(argv[3], "rb");
    gzFile output1 = gzopen(argv[4], "wb1");
    gzFile output2 = gzopen(argv[5], "wb1");
    if (input1 == NULL || input2 == NULL || output1 == NULL || output2 == NULL) {
        fprintf(stderr, "gzip open failed\n");
        return 4;
    }
    gzbuffer(input1, 1 << 20);
    gzbuffer(input2, 1 << 20);
    gzbuffer(output1, 1 << 20);
    gzbuffer(output2, 1 << 20);
    uint64_t input_pairs = 0;
    uint64_t retained_pairs = 0;
    int failed = 0;
    while (1) {
        char *first1 = gzgets(input1, r1_lines[0], MAX_LINE);
        char *first2 = gzgets(input2, r2_lines[0], MAX_LINE);
        if (first1 == NULL || first2 == NULL) {
            if (first1 == NULL && first2 == NULL && gzeof(input1) && gzeof(input2)) break;
            failed = 1;
            break;
        }
        for (int line = 1; line < 4; line++) {
            if (gzgets(input1, r1_lines[line], MAX_LINE) == NULL || gzgets(input2, r2_lines[line], MAX_LINE) == NULL) {
                failed = 1;
                break;
            }
        }
        if (failed || !valid_record(r1_lines) || !valid_record(r2_lines) || !mate_names_match(r1_lines[0], r2_lines[0])) {
            failed = 1;
            break;
        }
        input_pairs++;
        if (sequence_has_verified_match(r1_lines[1]) || sequence_has_verified_match(r2_lines[1])) {
            for (int line = 0; line < 4; line++) {
                if (gzputs(output1, r1_lines[line]) < 0 || gzputs(output2, r2_lines[line]) < 0) {
                    failed = 1;
                    break;
                }
            }
            if (failed) break;
            retained_pairs++;
        }
    }
    int close_status = gzclose(input1) | gzclose(input2) | gzclose(output1) | gzclose(output2);
    if (failed || close_status != Z_OK) {
        fprintf(stderr, "FASTQ validation or gzip I/O failed\n");
        return 5;
    }
    FILE *stats = fopen(argv[6], "w");
    if (stats == NULL) return 6;
    fprintf(stats, "metric\tvalue\n");
    fprintf(stats, "seed_length\t%d\n", KMER_LENGTH);
    fprintf(stats, "distinct_seeds\t%llu\n", (unsigned long long)seed_count);
    fprintf(stats, "input_pairs\t%llu\n", (unsigned long long)input_pairs);
    fprintf(stats, "retained_pairs\t%llu\n", (unsigned long long)retained_pairs);
    fprintf(stats, "retained_fraction\t%.12g\n", input_pairs ? (double)retained_pairs / (double)input_pairs : 0.0);
    fprintf(stats, "invalid_records\t0\n");
    fclose(stats);
    return 0;
}
