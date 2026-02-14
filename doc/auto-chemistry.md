# Chemistry Auto-Detection Logic Description

## 1. Overall Detection Workflow

1. **Sampling**:
   For each FASTQ file, the program extracts the first `max_read` reads (default: 10,000) for detection.

2. **Per-read Detection**:
   For each sampled read, the following matching logic is attempted sequentially:

   * **Linker/Offset Matching** (for V3 or flv-V2):
     A sliding window search (allowing 0–4 bp offset) is used to identify specific linker sequences.
     If a match is found, the chemistry is immediately determined.

   * **Whitelist Matching** (for V1/V2, etc.):
     Barcodes are extracted according to predefined patterns and compared against a whitelist
     (allowing up to one base mismatch).

3. **Counting & Thresholds**:

   * The program counts the frequency of detected chemistry types across all reads and selects the most frequent one as the final result.

   * **Threshold checks**:

     * If the proportion is **< 10%**:
       The program raises an error and terminates:
       `Exception: Auto chemistry detection failed!`
     * If the proportion is **10% – 50%**:
       The program proceeds but prints a low-confidence warning:
       `Valid chemistry read counts percent < 0.5`

4. **Consistency Check**:
   For multiple FASTQ file inputs, the program verifies that all files yield a consistent chemistry detection result.
   If inconsistent, the program exits with an error.

---

## 2. Distinguishing GEXSCOPE-V1 and flv_rna

These two chemistries share highly similar barcode and initial linker structures:

* **GEXSCOPE-V1**: `C8L16C8L16C8` + **`L1(C)`** + `U12`
* **flv_rna**: `C8L16C8L16C8` + `U9`

Since `flv_rna` lacks the one-base linker (fixed as C) present in V1, the program uses a **match-first, then refine** strategy.

### 1. Structural Matching

The program first calls:

```
is_chemistry(seq, "GEXSCOPE-V1")
```

Because the upstream structure is identical, `flv_rna` sequences can also pass the V1 barcode whitelist validation.

### 2. Position-Based Disambiguation (Position 57 Check)

After a successful structural match, the program determines the final chemistry by examining the **57th base** (index `seq[56]`):

| Detection Result | Decision Logic                                |
| ---------------- | --------------------------------------------- |
| **GEXSCOPE-V1**  | Matches V1 structure **and** `seq[56] == 'C'` |
| **flv_rna**      | Matches V1 structure **and** `seq[56] != 'C'` |

**Reference code:**
[https://github.com/singleron-RD/CeleScope/blob/master/celescope/tools/parse_chemistry.py#L380](https://github.com/singleron-RD/CeleScope/blob/master/celescope/tools/parse_chemistry.py#L380)
