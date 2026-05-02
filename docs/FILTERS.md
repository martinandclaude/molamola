# Filters

Each event in the VCF is checked against a few independent noise-flag rules. A flag does not drop the event from the output — flagged events still appear in the HTML report, but render greyed/dashed in the figures so the eye goes to the unflagged signal first.

The directional sense of the **raise** / **lower** columns is consistent: raising a threshold makes the filter looser (fewer events flagged), lowering makes it stricter (more events flagged).

## Defaults

| flag | default | raise → | lower → |
|---|---|---|---|
| `--cov-ratio` | `auto` (= `max(2.0, p99)`) | passing a higher fixed number (e.g. `--cov-ratio 4`) flags fewer events; only extreme coverage spikes get caught. | passing a lower fixed number (e.g. `--cov-ratio 2`) flags more events; even mild coverage anomalies are excluded. The default `auto` adapts per sample. |
| `--cov-vaf-max` | `0.35` | flags more high-coverage events as noise (some real mosaic / subgermline events may get greyed). | flags fewer; only unambiguously low-VAF artefacts. |
| `--mark-acrocentric` | on for hg38, off for t2t | binary — on means chr13/14/15/21/22 p-arm-only BNDs render grey. Off when investigating real acrocentric biology, or for T2T (where those p-arms are real, fully-resolved sequence). |  |
| `--min-svlen` | `50` | shows only larger SVs everywhere (density tracks and stdout summaries); cleaner picture, but small events (Alu insertions, micro-deletions) drop out. | includes very small events (10–50 bp); polymorphism noise dominates. Set `0` to disable. |
| `--focus-window` | `1000` | accepts more candidate matches around a focus position; useful for IMPRECISE breakpoints. | tighter exact-position match. |

## What each filter does

### Coverage anomaly (`--cov-ratio` + `--cov-vaf-max`)

A long-read SV with `max(COVERAGE) > N × genome-median-coverage` AND `VAF < threshold` is the canonical repeat-collapse / mismapping signature. Read pile-up is high (multiple repeat copies aligning to the same locus) but the variant fraction stays low (most reads agree with reference). Default flag rule: VAF < 0.35 *and* the per-event `max(COVERAGE) / median_coverage` is at or above the threshold described next.

`--cov-ratio` defaults to `auto`, which sets the threshold per sample to `max(2.0, p99 of the in-sample max-coverage / median-coverage distribution)`. This adapts to each sample's coverage profile (high-noise samples get a stricter cutoff; clean samples fall to the 2.0× floor). The chosen value is printed at run start, e.g. `--cov-ratio auto: 2.72x (p99 of 25,167 PASS events; floor 2.0x)`.

Pass a fixed number to override (int or float — e.g. `--cov-ratio 3` or `--cov-ratio 2.5` — reproduces a literature-style constant for samples where you'd rather pin the threshold).

### Acrocentric short arms (`--mark-acrocentric`)

On hg38, the p-arms of chr13/14/15/21/22 are largely unresolved (rDNA arrays, satellite DNA). BNDs with both endpoints in those p-arms are almost always mapping artefacts on hg38. The flag is off by default for T2T-CHM13v2.0, where those regions are properly resolved sequence.

### Minimum size (`--min-svlen`)

Non-BND SVs (INS / DEL / DUP / INV) shorter than `--min-svlen` (default 50 bp) are dropped from every downstream consumer: density tracks and stdout summaries. BND records have no SVLEN and are not affected. Set `--min-svlen 0` to keep every event.

### Focus window (`--focus-window`)

When `--focus CHR:POS` is given, BNDs whose endpoints fall within `±--focus-window bp` of the requested coordinate are kept and the rest are filtered out for that figure.
