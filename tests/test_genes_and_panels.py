"""Gene tables, breakpoint gene labels, fusion naming, derivative panels.

The gene tables are NCBI RefSeq (RS_2025_08), native on both builds.
Fusion names follow from the kept sides and the gene strands alone; the
derivative panels are built from the same junction sides. A known
answer - t(8;21) makes RUNX1::RUNX1T1 on der(8), inv(16) makes
CBFB::MYH11 - pins both.
"""

from __future__ import annotations

import gzip
import importlib.util

import numpy as np
import pytest

import molamola as mm

P8, P21 = 92_050_000, 34_900_000
INV16 = (15_800_000, 67_070_000)


@pytest.fixture(scope="module")
def cyto(bundled_cytoband):
    return mm.load_cytobands(bundled_cytoband)


@pytest.fixture(scope="module")
def genes():
    return mm.load_gene_table(mm.find_gene_file("hg38"))


def J(ca, pa, sa, cb, pb, sb, ident="j", source="BND"):
    return mm.Junction(ca, pa, sa, cb, pb, sb, support=20, vaf=0.4,
                       sv_ids=(ident,), source=source)


T821 = (J("chr8", P8, "L", "chr21", P21, "R", "der8"),
        J("chr8", P8 + 10, "R", "chr21", P21 - 10, "L", "der21"))
INV = (J("chr16", INV16[0], "L", "chr16", INV16[1], "L", "i1"),
       J("chr16", INV16[0], "R", "chr16", INV16[1], "R", "i2"))


# --- bundled tables --------------------------------------------------------

@pytest.mark.parametrize("build", ["hg38", "t2t"])
def test_bundled_tables_carry_genes_and_ig_tr_loci(build):
    t = mm.load_gene_table(mm.find_gene_file(build))
    names = {g.name for c in t.values() for g in c[2]}
    assert {"RUNX1", "RUNX1T1", "CBFB", "MYH11", "KMT2A", "CRLF2"} <= names
    assert {"IGH", "IGK", "IGL", "TRA/TRD", "TRB", "TRG"} <= names
    assert len(names) > 19_000


def test_t2t_coordinates_are_native_and_match_nasvar():
    """NASVAR's own T2T leukaemia config places RUNX1 at
    33170406-33432142 (1-based). The native RefSeq table must agree -
    a liftover would not, exactly."""
    t = mm.load_gene_table(mm.find_gene_file("t2t"))
    runx1 = [g for g in t["chr21"][2] if g.name == "RUNX1"]
    assert [(g.start + 1, g.end) for g in runx1] == [(33_170_406, 33_432_142)]


# --- derivation script -----------------------------------------------------

def _load_script(repo_root):
    spec = importlib.util.spec_from_file_location(
        "derive_gene_tables", repo_root / "scripts" / "derive_gene_tables.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_locus_of_does_not_sweep_in_ordinary_genes(repo_root):
    d = _load_script(repo_root)
    assert d.locus_of("IGHV3-21") == "IGH"
    assert d.locus_of("IGHM") == "IGH"
    assert d.locus_of("TRDV2") == "TRA/TRD"
    assert d.locus_of("TRBC1") == "TRB"
    assert d.locus_of("IGHMBP2") is None
    assert d.locus_of("TRAF1") is None


def test_derivation_keeps_primary_protein_coding_and_drops_orphons(
        repo_root, tmp_path):
    d = _load_script(repo_root)
    cols = ["feature", "class", "assembly", "assembly_unit", "seq_type",
            "chromosome", "genomic_accession", "start", "end", "strand",
            "product_accession", "non-redundant_refseq", "related_accession",
            "name", "symbol", "GeneID", "locus_tag",
            "feature_interval_length", "product_length", "attributes"]

    def row(cls, chrom, start, end, sym, seq_type="chromosome", feat="gene"):
        r = dict.fromkeys(cols, "")
        r.update(feature=feat, **{"class": cls}, seq_type=seq_type,
                 chromosome=chrom, start=str(start), end=str(end),
                 strand="+", symbol=sym)
        return "\t".join(r[c] for c in cols)

    lines = ["# " + "\t".join(cols),
             row("protein_coding", "21", 101, 200, "GENEA"),
             row("protein_coding", "21", 101, 200, "ALTCOPY",
                 seq_type="alternate scaffold"),
             row("protein_coding", "21", 101, 200, "NOTGENE", feat="mRNA"),
             row("lncRNA", "21", 101, 200, "LNC"),
             row("V_segment", "14", 1_000_001, 1_010_000, "IGHV1-1"),
             row("C_region", "14", 1_200_001, 1_210_000, "IGHM"),
             row("V_segment", "14", 20_000_001, 20_010_000, "IGHV1-99"),
             row("V_segment", "15", 1_000_001, 1_010_000, "IGHV1-OR15")]
    ft = tmp_path / "ft.txt.gz"
    with gzip.open(ft, "wt") as fh:
        fh.write("\n".join(lines) + "\n")
    rows = d.derive(ft)
    assert ("chr21", 100, 200, "GENEA", 0, "+") in rows
    assert not any(r[3] in {"ALTCOPY", "NOTGENE", "LNC"} for r in rows)
    igh = [r for r in rows if r[3] == "IGH"]
    assert igh == [("chr14", 1_000_000, 1_210_000, "IGH", 0, ".")]


# --- gene lookup -----------------------------------------------------------

def test_breakpoint_labels(genes):
    assert mm.breakpoint_gene_label(genes, "chr21", P21) == "RUNX1"
    assert mm.breakpoint_gene_label(genes, "chr14", 106_000_000) == "IGH"
    # 157 kb downstream of MYC, which has no protein-coding neighbour closer
    assert mm.breakpoint_gene_label(genes, "chr8", 127_900_000) == \
        "near MYC (157 kb)"
    assert mm.breakpoint_gene_label(genes, "chrZZ", 1) == ""


def test_an_ig_locus_wins_over_genes_inside_it(genes):
    hits = mm.genes_at(genes, "chr14", 106_000_000)
    assert [g.name for g in hits] == ["IGH"]


# --- fusion naming ---------------------------------------------------------

def test_t821_fuses_runx1_to_runx1t1_on_der8(genes):
    assert mm.junction_fusion(genes, T821[0]) == "RUNX1::RUNX1T1"
    assert mm.junction_fusion(genes, T821[1]) == "RUNX1T1::RUNX1"


def test_inv16_fuses_cbfb_to_myh11(genes):
    assert mm.junction_fusion(genes, INV[0]) == "CBFB::MYH11"
    assert mm.junction_fusion(genes, INV[1]) == "MYH11::CBFB"


def test_head_to_head_junction_is_not_a_fusion(genes):
    """RUNX1T1 and RUNX1 are both on the minus strand; keeping the left of
    both puts them head to head - no fusion transcript."""
    j = J("chr8", P8, "L", "chr21", P21, "L")
    assert mm.junction_fusion(genes, j) == ""


def test_ig_locus_is_named_first(genes):
    j = J("chr14", 106_000_000, "L", "chrX", 1_200_000, "R")
    assert mm.junction_fusion(genes, j) == "IGH::CRLF2"


def test_breakpoint_outside_genes_makes_no_fusion(genes):
    assert mm.junction_fusion(genes, J("chr8", 127_900_000, "L",
                                       "chr21", P21, "R")) == ""


# --- promotion -------------------------------------------------------------

def test_masked_pair_with_genes_at_every_breakpoint_is_promoted(cyto, genes):
    mask = mm.load_mask_intervals(mm.find_mask_file("hg38"))
    # both INV16 positions fall in masked 500 bp runs (measured)
    assert mm.in_mask(mask, "chr16", INV16[0])
    assert mm.in_mask(mask, "chr16", INV16[1])
    without = mm.classify_rearrangements(list(INV), cyto, mask)
    with_genes = mm.classify_rearrangements(list(INV), cyto, mask, genes)
    assert without[0].tier == "repeat"
    assert with_genes[0].tier == "candidate"
    assert with_genes[0].notes


def test_masked_pair_outside_genes_stays_demoted(cyto, genes):
    # 8q24.21 gene desert beside MYC, on chr21 inside RUNX1
    js = [J("chr8", 127_900_000, "L", "chr21", P21, "R"),
          J("chr8", 127_900_010, "R", "chr21", P21 - 10, "L")]
    mask = {"chr8": (np.array([127_899_000]), np.array([127_901_000])),
            "chr21": (np.array([P21 - 1_000]), np.array([P21 + 1_000]))}
    ev = mm.classify_rearrangements(js, cyto, mask, genes)
    assert ev[0].tier == "repeat"


# --- derivatives -----------------------------------------------------------

def test_translocation_derivatives(cyto):
    ev = mm.Rearrangement("translocation", T821, "candidate")
    cols = mm.panel_columns(ev, cyto)
    assert [c[0] for c in cols] == ["8", "der(8)", "21", "der(21)"]
    der8, der21 = cols[1][1], cols[3][1]
    L8, L21 = mm._chrom_len(cyto, "chr8"), mm._chrom_len(cyto, "chr21")
    assert der8 == [("chr8", 0, P8, False), ("chr21", P21, L21, False)]
    assert der21 == [("chr21", 0, P21 - 10, False),
                     ("chr8", P8 + 10, L8, False)]


def test_inversion_derivative_reverses_the_middle(cyto):
    ev = mm.Rearrangement("inversion", INV, "candidate")
    cols = mm.panel_columns(ev, cyto)
    assert [c[0] for c in cols] == ["16", "inv(16)"]
    assert cols[1][1][1] == ("chr16", INV16[0], INV16[1], True)
    assert [b for b, _ in cols[1][2]] == [0, 1]


def test_insertion_derivatives(cyto):
    js = (J("chr8", P8, "L", "chr21", P21, "R"),
          J("chr8", P8 + 20, "R", "chr21", P21 + 3_000_000, "L"))
    cols = mm.panel_columns(mm.Rearrangement("insertion", js, "candidate"),
                            cyto)
    assert [c[0] for c in cols] == ["8", "der(8)", "21", "der(21)"]
    assert cols[1][1][1] == ("chr21", P21, P21 + 3_000_000, False)
    assert len(cols[3][1]) == 2  # donor chromosome minus the segment


def test_orient_flips_a_derivative_read_from_its_qter(cyto):
    L = mm._chrom_len(cyto, "chr8")
    pieces = [("chr21", 40_000_000, 46_000_000, True),
              ("chr8", 0, L, True)]
    oriented, name = mm._orient(pieces, cyto)
    assert name == "der(8)"
    assert oriented[0] == ("chr8", 0, L, False)


def test_two_centromeres_make_a_dicentric(cyto):
    pieces = [("chr8", 0, 100_000_000, False),
              ("chr21", 0, 30_000_000, True)]
    assert mm._orient(pieces, cyto)[1] == "dic(8;21)"


def test_panel_renders(cyto, genes):
    png = mm.render_rearrangement_panel(
        mm.Rearrangement("translocation", T821, "candidate"), 1, cyto,
        genes, "hg38")
    assert png.startswith(b"\x89PNG")


def test_event_title_carries_the_fusion(cyto, genes):
    ev = mm.Rearrangement("translocation", T821, "candidate")
    assert mm.event_title(ev, cyto, genes) == \
        "t(8;21)(q21.3;q22.12)  RUNX1::RUNX1T1"


# --- report ----------------------------------------------------------------

def _t821_vcf(tmp_path):
    header = (
        "##fileformat=VCFv4.2\n##source=Sniffles2_2.8.0\n"
        "##contig=<ID=chr1,length=248956422>\n"
        "##contig=<ID=chr8,length=145138636>\n"
        "##contig=<ID=chr21,length=46709983>\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"t\">\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
    info = "SUPPORT=20;COVERAGE=30,30,30,30,30;VAF=0.400"
    recs = [
        f"chr8\t{P8}\tB1\tG\tG[chr21:{P21}[\t60\tPASS\tSVTYPE=BND;{info}\tGT\t0/1",
        f"chr8\t{P8 + 10}\tB2\tG\t]chr21:{P21 - 10}]G\t60\tPASS\tSVTYPE=BND;{info}\tGT\t0/1",
    ]
    p = tmp_path / "t821.vcf"
    p.write_text(header + "\n".join(recs) + "\n")
    return p


def test_report_embeds_one_panel_per_candidate(tmp_path):
    out = tmp_path / "o"
    assert mm.main(["--vcf", str(_t821_vcf(tmp_path)), "--out", str(out),
                    "--png"]) == 0
    html = (out / "t821.report.html").read_text()
    assert 'id="fig-rearrangements"' in html
    assert "rearrangement 1: t(8;21)(q21.3;q22.12)  RUNX1::RUNX1T1" in html
    assert (out / "t821.report.rearrangement_1.png").read_bytes()[:4] == \
        b"\x89PNG"


def test_report_without_candidates_has_no_panel_section(tiny_vcf, tmp_path):
    mm.main(["--vcf", str(tiny_vcf), "--out", str(tmp_path)])
    html = (tmp_path / "tiny.report.html").read_text()
    assert 'id="fig-rearrangements"' not in html


def test_panels_html_counts_what_the_cap_left_out():
    html = mm._panels_html([(1, "t(1;2)", b"x")], n_candidates=3)
    assert "2 further candidates not drawn" in html
    assert mm._panels_html([], n_candidates=0) == ""
