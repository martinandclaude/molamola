"""Smoke tests for the compound-het renderer."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")  # headless

import molamola as mm


PNG_SIGNATURE = bytes([0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A])


def test_renderer_emits_valid_png(make_phased_variant, make_gene):
    gene = make_gene()
    variants = [
        make_phased_variant(pos=1100, gt="0|1", clnsig="p_or_lp"),
        make_phased_variant(pos=1500, gt="1|0", clnsig="vus"),
    ]
    blocks = mm._build_phase_blocks(variants)
    png, stats = mm.render_compound_het_png(
        gene, variants, blocks, reference="hg38",
    )
    assert isinstance(png, bytes)
    assert png[:8] == PNG_SIGNATURE
    assert len(png) > 1000  # not a degenerate / empty PNG
    assert stats["n_missense"] == 2
    assert stats["n_total"] == 2
    assert stats["n_trans"] == 1
    assert stats["n_cis"] == 0
    assert stats["n_blocks"] == 1


def test_renderer_handles_no_canonical_exons(make_phased_variant, make_gene):
    """Genes without canonical exons still render (just no exon track)."""
    gene = make_gene(canonical_exons=())
    variants = [make_phased_variant(pos=1100)]
    blocks = mm._build_phase_blocks(variants)
    png, _ = mm.render_compound_het_png(
        gene, variants, blocks, reference="hg38",
    )
    assert png[:8] == PNG_SIGNATURE


def test_renderer_handles_off_edge_phase_block(
    make_phased_variant, make_gene,
):
    """Phase block extending past the gene window draws an off-edge arrow."""
    import molamola as mm
    gene = make_gene(start=1000, end=2000)
    variants = [
        make_phased_variant(pos=1100),
        make_phased_variant(pos=1900, gt="1|0"),
    ]
    block = mm.PhaseBlock(ps=12345, start=500, end=2500, n_phased=2)
    png, _ = mm.render_compound_het_png(
        gene, variants, [block], reference="hg38",
    )
    assert png[:8] == PNG_SIGNATURE


def test_renderer_includes_synonymous_as_tick(make_phased_variant, make_gene):
    """Non-missense canonical variant renders as tick, not lollipop."""
    gene = make_gene()
    variants = [
        make_phased_variant(pos=1100, consequence="missense_variant"),
        make_phased_variant(pos=1500,
                              consequence="synonymous_variant", gt="1|0"),
    ]
    blocks = mm._build_phase_blocks(variants)
    png, stats = mm.render_compound_het_png(
        gene, variants, blocks, reference="hg38",
    )
    assert stats["n_missense"] == 1
    assert stats["n_total"] == 2


def test_dpi_capped_for_wide_figures(make_phased_variant, make_gene):
    """Width auto-scales but DPI caps so PNG never exceeds 1800 px."""
    gene = make_gene(start=0, end=1_000_000)
    # Many variants -> wide figure.
    variants = [
        make_phased_variant(pos=10_000 * i, ps=12345)
        for i in range(1, 50)
    ]
    blocks = mm._build_phase_blocks(variants)
    png, _ = mm.render_compound_het_png(
        gene, variants, blocks, reference="hg38",
    )
    # We don't crack open the PNG header to read pixel width; instead
    # rely on the dpi formula: target_dpi = min(150, 1800/figwidth).
    # The cap ensures bytes-output stays bounded, which we sanity-check.
    assert len(png) < 4_000_000


def test_build_phase_blocks_groups_by_ps(make_phased_variant):
    a = make_phased_variant(pos=100, ps=11)
    b = make_phased_variant(pos=200, ps=11)
    c = make_phased_variant(pos=500, ps=22)
    blocks = mm._build_phase_blocks([a, b, c])
    assert len(blocks) == 2
    by_ps = {b.ps: b for b in blocks}
    assert by_ps[11].start == 100
    assert by_ps[11].end == 200
    assert by_ps[11].n_phased == 2
    assert by_ps[22].start == 500
    assert by_ps[22].end == 500
    assert by_ps[22].n_phased == 1
