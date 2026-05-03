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

## Quick start

```sh
molamola --vcf sample.vcf
```

The plot type is auto-detected. Output is a single self-contained
HTML report — figures embedded as base64, no external assets, opens
offline.

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
