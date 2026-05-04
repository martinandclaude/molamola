<p align="center"><img src="header_fish.png" alt="molamola" width="128"></p>

# molamola

A Python plotting tool for Oxford Nanopore variation data.
**One VCF in, one self-contained HTML report out.**

molamola inspects the VCF header and picks the right plot type
automatically — no flags or subcommands to remember:

- **SV / cytogenetics report** for long-read SV VCFs
  (Sniffles2, cuteSV, SVIM, pbsv, NanoVar).
- **Per-gene phased-haplotype panels** for compound-het workup
  from phased + VEP-annotated small-variant VCFs
  (WhatsHap, HiPhase).

## Install

```sh
pip install molamola
```

Or via conda from the bioconda channel:

```sh
conda install -c bioconda molamola
```

## Quick start

```sh
molamola --vcf sample.vcf
```

The plot type is auto-detected. Output is a single self-contained
HTML report — figures embedded as base64, no external assets, opens
offline.

## Preparing a phased VCF for compound-het mode

Compound-het mode needs both phasing (`PS` FORMAT field) and VEP
annotation (`CSQ` INFO field). A raw phased small-variant VCF —
e.g. straight Clair3 output — has the first but not the second,
and molamola will refuse it. Annotate with
[Ensembl VEP](https://github.com/Ensembl/ensembl-vep) first.

**1. Download the matching VEP cache once** (one-time, ~20 GB),
on a machine with internet:

```sh
wget https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz
```

If your VEP isn't 105, swap `105` for your release number (it
appears twice in the URL); the cache release must match the VEP
release exactly. Transfer to wherever you run VEP if that's a
different machine.

**2. Unpack into a stable cache directory:**

```sh
mkdir -p VEP_cache && cd VEP_cache
tar -xzf ../homo_sapiens_vep_105_GRCh38.tar.gz
# creates VEP_cache/homo_sapiens/105_GRCh38/
```

**3. Run VEP fully offline**, with `--canonical --symbol --pick`
so the CSQ shape matches what molamola consumes:

```sh
vep --input_file sample.phased.vcf \
    --output_file sample.phased.vep.vcf \
    --vcf --offline \
    --cache --dir_cache /path/to/VEP_cache \
    --assembly GRCh38 \
    --fasta /path/to/hg38.fa \
    --canonical --symbol --pick \
    --force_overwrite
```

Then `molamola --vcf sample.phased.vep.vcf` picks it up as
compound-het mode.

**Notes on VEP.** VEP is third-party software (Ensembl); molamola
does not bundle or wrap it. The cache release and VEP binary
release must match exactly — a mismatch leads to silent
mis-annotation rather than a clean error. Compound-het mode reads
VEP's `Consequence`, `SYMBOL`, and `CANONICAL` fields as-is; any
quirks of a particular VEP build are inherited. `--pick` reduces
multi-transcript CSQ entries to one per variant.

## Documentation

- [CLI reference](CLI.md) — every flag, with defaults and meanings.
- [Examples](EXAMPLES.md) — worked commands for each plot type.
- [Filters](FILTERS.md) — noise heuristics, focus windows, ClinVar gating.
- [Output spec](OUTPUTS.md) — what's in the HTML; how figures are encoded.

## Source

- Repo: [github.com/martinandclaude/molamola](https://github.com/martinandclaude/molamola)
- PyPI: [pypi.org/project/molamola](https://pypi.org/project/molamola/)
- Changelog: [CHANGELOG.md](https://github.com/martinandclaude/molamola/blob/main/CHANGELOG.md)
- Issues: [github.com/martinandclaude/molamola/issues](https://github.com/martinandclaude/molamola/issues)
