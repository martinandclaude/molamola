# Filters

This page covers three mechanisms: **SV-mode noise flags** (events stay in the report, just greyed), **rearrangement tiers** (how BND junctions are paired and which pairs become candidates — display emphasis only) and **karyotype-mode QC** (bins / sites are dropped before plotting). Jump to [Rearrangement tiers](#rearrangement-tiers-sv-mode) or [Karyotype-mode QC](#karyotype-mode-qc).

## SV-mode noise flags

Each event in the VCF is checked against a few independent noise-flag rules. A flag does not drop the event from the output — flagged events still appear in the HTML report, but render grey and faint in the figures so the eye goes to the unflagged signal first.

The directional sense of the **raise** / **lower** columns is consistent: raising a threshold makes the filter looser (fewer events flagged), lowering makes it stricter (more events flagged).

## Defaults

| flag | default | raise → | lower → |
|---|---|---|---|
| `--cov-ratio` | `auto` (= `max(2.0, p99)`) | passing a higher fixed number (e.g. `--cov-ratio 4`) flags fewer events; only extreme coverage spikes get caught. | passing a lower fixed number (e.g. `--cov-ratio 2`) flags more events; even mild coverage anomalies are excluded. The default `auto` adapts per sample. |
| `--cov-vaf-max` | `0.35` | flags more high-coverage events as noise (some real mosaic / subgermline events may get greyed). | flags fewer; only unambiguously low-VAF artefacts. |
| `--mark-acrocentric` | on for hg38, off for t2t | binary — on means chr13/14/15/21/22 p-arm-only BNDs render grey. Off when investigating real acrocentric biology, or for T2T (where those p-arms are real, fully-resolved sequence). |  |
| `--min-svlen` | `50` | shows only larger SVs everywhere (density rings and stdout summaries); cleaner picture, but small events (Alu insertions, micro-deletions) drop out. | includes very small events (10–50 bp); polymorphism noise dominates. Set `0` to disable. |
| `--focus-window` | `1000` | accepts more candidate matches around a focus position; useful for IMPRECISE breakpoints. | tighter exact-position match. |

## What each filter does

### Coverage anomaly (`--cov-ratio` + `--cov-vaf-max`)

A long-read SV with `max(COVERAGE) > N × genome-median-coverage` AND `VAF < threshold` is the canonical repeat-collapse / mismapping signature. Read pile-up is high (multiple repeat copies aligning to the same locus) but the variant fraction stays low (most reads agree with reference). Default flag rule: VAF < 0.35 *and* the per-event `max(COVERAGE) / median_coverage` is at or above the threshold described next.

`--cov-ratio` defaults to `auto`, which sets the threshold per sample to `max(2.0, p99 of the in-sample max-coverage / median-coverage distribution)`. This adapts to each sample's coverage profile (high-noise samples get a stricter cutoff; clean samples fall to the 2.0× floor). The chosen value is printed at run start, e.g. `--cov-ratio auto: 2.72x (p99 of 25,167 PASS events; floor 2.0x)`.

Pass a fixed number to override (int or float — e.g. `--cov-ratio 3` or `--cov-ratio 2.5` — reproduces a literature-style constant for samples where you'd rather pin the threshold).

### Acrocentric short arms (`--mark-acrocentric`)

On hg38, the p-arms of chr13/14/15/21/22 are largely unresolved (rDNA arrays, satellite DNA). BNDs with both endpoints in those p-arms are almost always mapping artefacts on hg38. The flag is off by default for T2T-CHM13v2.0, where those regions are properly resolved sequence.

### Minimum size (`--min-svlen`)

Non-BND SVs (INS / DEL / DUP / INV) shorter than `--min-svlen` (default 50 bp) are dropped from every downstream consumer: density rings and stdout summaries. BND records have no SVLEN and are not affected. Set `--min-svlen 0` to keep every event.

### Focus window (`--focus-window`)

When `--focus CHR:POS` is given, BNDs whose endpoints fall within `±--focus-window bp` of the requested coordinate are kept and the rest are filtered out for that figure. Rearrangements are still classified on the whole call set first, so a focused BND is drawn as part of its event — paired with its reciprocal junction even when that one lies outside the window.

## Rearrangement tiers (SV mode)

Separately from the noise flags, molamola groups BND junctions into rearrangement events, so that a balanced translocation or inversion is not just one arc among a hundred. The circos draws one arc per event, emphasised by tier and with candidates numbered and named (see [Outputs](OUTPUTS.md#circos-plot-sv-mode)). The run prints one line per sample, e.g. `Rearrangements: 1 candidate (1 translocation), 11 paired in repeats, 83 single junctions`, which the report's metadata block repeats, and `--focus` names each matched BND by the event it belongs to.

A BND record is one **junction**: two breakpoints joined, each keeping the sequence on one side of it. A balanced translocation leaves **two** junctions, one per derivative chromosome. These sit within `100 kb` of each other at both breakpoints and keep opposite sides. Artefacts usually leave one. Records describing the same junction twice (cuteSV writes each mate separately) are merged first.

| event | what it is |
|---|---|
| translocation | reciprocal junction pair between two chromosomes, named `t(A;B)(band;band)` |
| inversion | reciprocal pair on one chromosome with inversion sides, or an `INV` record ≥ 1 Mb, named `inv(A)(band band)` |
| insertion | pair sharing one breakpoint, with the other ends bounding a donor segment, named `ins(A;B)(...)` |
| single | an unpaired junction, named `8q21.3::21q22.12 (single junction)`. Never `t()`: one junction does not show the event is balanced. Unbalanced events are singles by nature; the karyotype report shows their copy-number side. |

Intrachromosomal events shorter than 1 Mb are left out; small inversions are common polymorphisms.

Each paired event is then a **candidate** unless one of these applies, in which case it is **paired in repeats**:

- a junction is noise-flagged (see above);
- a breakpoint is within 2 Mb of an `acen` / `gvar` / `stalk` band (whole-arm and Robertsonian events land here, as does most normal-genome background);
- every breakpoint falls inside the bundled exclusion mask (no padding: the mask is stored as 500 bp runs and covers 40 % of hg38 at 0 bp, 84 % at 1 kb);
- the pair has the **mobile-element shape**: within 100 bp at one end (the target-site duplication of a non-reference L1 / Alu / SVA) and bounding a 100 bp – 15 kb shared segment at the other (the reference copy it aligned to).

One exception: a pair demoted **only** because every breakpoint is masked is promoted back to candidate when every breakpoint also falls inside a gene (the panel says so). The mask covers 40 % of hg38 and would otherwise demote ~7–24 % of real recurrent fusions; in the ten normal call sets no mask-demoted pair had a gene at every breakpoint, so the promotion costs no measurable specificity. A synthetic inv(16) whose breakpoints happen to fall in masked runs comes out as a candidate, `CBFB::MYH11`, because of it.

Tiers set display emphasis only; nothing is dropped. The trade-off is measured, not assumed. On ten normal call sets (Sniffles2 2.7.5 / 2.8.0 and cuteSV 2.1.3, hg38 and T2T), false translocation candidates fall from 13 to 1 with the mask rule. But with breakpoints uniform over MANE gene bodies, the same rule would demote 7–24 % of typical recurrent leukaemia fusions (median 11 %), all P2RY8::CRLF2 and all DUX4 events. Most false candidates that remain in normals are large Sniffles2 `INV` records, 0–6 per genome. Several recur across unrelated people, so they are either common polymorphic inversions or recurrent caller artefacts. cuteSV caps SV length at ~100 kb by default, so a cuteSV VCF may carry no large inversion at all.

## Karyotype-mode QC

Unlike the SV noise flags, these **drop** data before plotting — a masked bin or a failed BAF site is gone from the figure, not greyed. The intent is a clean CN / BAF surface for cytogenetic review rather than an exhaustive dump.

### Coverage bins: exclusion mask + centromere pad

Before CN is computed, every mosdepth bin is tested against the bundled exclusion mask (`--mask`, or `--no-mask` to disable) — the union of low-mappability regions and a polymorphic-TR catalog. A bin is dropped when **more than `--mask-overlap` (default 0.5) of it** falls inside the mask. Overlap is measured in base pairs rather than as a yes/no touch: the bundled masks are built from 500 bp runs, so an any-overlap test would scale the exclusion with your mosdepth bin size instead of with the mask — on a 1 kb run it excluded 57.8 % of bins for a mask covering 39.8 % of hg38. Half-masked 1 kb bins measure like clean sequence (SD 0.230 log2 vs 0.225), while fully-masked bins sit at 1.072, so majority-masked is the empirically right cut. `--mask-overlap 0` restores the pre-v0.5.0 behaviour. Masked bins are dropped from the depth scatter, the rolling-median smooth, **and** the autosomal-median normalisation anchor (so a few residual high-copy repeats can't skew the CN 2.0 baseline). The mask is additionally widened by `--centromere-pad-kb` (default 1000) either side of every `acen` band, because centromeric depth stays unreliable on ONT even after the TR mask.

Non-uniform mosdepth inputs are refused outright: molamola expects `mosdepth --by <int>` (a single fixed bin size). A `regions.bed.gz` from `--by some.bed` (variable-width target bins) exits 1 with a clear error rather than producing a meaningless aggregation.

### GC correction (`--gc` / `--no-gc`)

Depth is corrected by a per-1 % GC-bucket median ratio (bundled 10 kb GC table) before normalisation. This is a mild correction on well-prepared ONT data; `--no-gc` disables it. It is not a filter (no bins are dropped) but it shapes the CN values, so it is listed here for completeness.

### BAF sites (only when `--vcf` is given)

The optional BAF panel needs a small-variant VCF. The whole input is refused upfront if its header carries `##INFO=<ID=SVTYPE,...>` — that catches an SV / CNV / BND VCF (Sniffles2 / cuteSV / SVIM / pbsv / NanoVar / Spectre / hificnv) passed by mistake; the error points at the small-variant call set instead.

Each record then has to pass, or it is silently skipped:

| rule | default | why |
|---|---|---|
| `FILTER == PASS` | — | only confident calls. |
| biallelic het GT | — | BAF is only meaningful at heterozygous sites (`0/1`, `1/0`, `0\|1`, `1\|0`). |
| SNV-only | — | `REF` and every `ALT` must be a single base. Indel-het AF is noisier on ONT (alignment ambiguity around the breakpoint inflates the estimate); the panel reads cleaner with SNVs only. Also drops symbolic `<...>` ALTs. |
| `FORMAT/DP >= --min-baf-dp` | `10` | below ~10 reads the alt-fraction estimate is too noisy to place on the panel. |
| `FORMAT/GQ >= --min-baf-gq` | `20` | de-facto floor for a confident het call across Clair3 / DeepVariant / GATK. Skipped when the VCF doesn't emit GQ. |
| canonical chrom | — | chr1–22, X, Y only. |

Raising `--min-baf-dp` or `--min-baf-gq` thins the panel to higher-confidence sites; lowering them (e.g. `--min-baf-dp 5 --min-baf-gq 0`) keeps more, at the cost of a noisier cloud. `--min-baf-gq 0` effectively disables the GQ filter.

### Adaptive-sampling detection (informational, not a filter)

If more than 5 % of autosomal non-masked bins sit above 5× the autosomal median, the depth distribution is treated as bimodal (adaptive-sampling-like) and an **"AS suspected"** chip is added to the run-metadata, plus a one-line stderr warning. Nothing is dropped — but CN is anchored on the autosomal median, which on an AS sample sits at the off-target background, so the CN scale is biased. Stratum-aware AS normalisation is a future feature; for now the flag is a read-with-care signal.
