# qPCR Primer & Probe Design for Pathogen Detection

Designs TaqMan primer/probe sets from a FASTA target, scoring candidates on
nearest-neighbour thermodynamics and secondary structure.

Shipped target: SARS-CoV-2 nucleocapsid (N) gene, `NC_045512.2:28274-29533`, 1,260 bp.

## Quick start

No packages to install. Python 3.8+.

```bash
python qpcr_design.py
```

Output on the bundled target: **442** primers pass filtering, **830** complete
assay sets form, **8** remain after redundancy removal.

| | Sequence (5'→3') | Tm | GC | Pos |
|---|---|---|---|---|
| Forward | `GCATTACGTTTGGTGGACCCTCAG` | 59.9 °C | 54.2 % | 41 |
| Probe | `TCGCGCCCCACTGCGTTCTCCATTCTG` | 69.1 °C | 63.0 % | 82 |
| Reverse | `TATTGGGTAAACCTTGGGGCCGA` | 60.1 °C | 52.2 % | 120 |

102 bp amplicon · ΔTm 0.14 °C between primers · probe melts 9.0 °C above them.

## Files

| File | Purpose |
|---|---|
| `qpcr_design.py` | The designer. FASTA in, ranked CSV out. |
| `qpcr_report.py` | Renders the same run as a self-contained HTML page. |
| `pathogen.fasta` | Target sequence. |
| `example_output_assays.csv` | Sample output from the bundled target. |
| `qPCR_Assay_Design_Report.docx` | Formatted write-up of the results. |
| `.vscode/launch.json` | Run configurations for VS Code (F5). |

Both scripts resolve `pathogen.fasta` against their own folder, so the working
directory does not matter.

## Options

```
--fasta PATH          target FASTA, one or more records   (default: pathogen.fasta)
--out PATH            output CSV                          (default: qpcr_assays.csv)
--report-fasta PATH   also write oligos as FASTA, for BLAST
--top N               assays to report                     (default: 10)
--min-amplicon N      shortest product, bp                 (default: 70)
--max-amplicon N      longest product, bp                  (default: 150)
--na X                monovalent cation, mol/L             (default: 0.05)
--primer-conc X       oligo concentration, mol/L           (default: 250e-9)
```

Multi-record FASTA files are supported; each record is designed separately, so no
chimeric primer can span a record junction.

## Visual report

```bash
python qpcr_report.py
```

Writes `qpcr_report.html` — an assay map drawn to scale along the target, a Tm
plot, and the results table. Open it in any browser. No plotting libraries used.

## How it works

1. **Enumerate** every 18–25 nt window on the plus strand.
2. **Filter** on GC (40–60 %), Tm (58–62 °C), homopolymer runs, 3' GC clamp,
   hairpin stem, self-dimer and 3' self-dimer. Windows containing any non-ACGT
   base are skipped.
3. **Pair** forward and reverse primers, bucketed by start position so only
   candidates within the amplicon window are compared — linear in sequence
   length rather than quadratic in primer count.
4. **Reject** pairs with cross-dimer or extendable 3' complementarity.
5. **Place a probe** in the amplicon interior, on either strand, that melts
   6–12 °C above both primers, does not start with G, and prefers C over G.
6. **Score** on primer Tm, ΔTm, GC, probe Tm offset, amplicon length and
   cross-dimer, then drop assays whose forward primers sit within 20 bp of a
   higher-scoring one.

### Thermodynamics

Melting temperatures use the SantaLucia (1998) unified nearest-neighbour
parameters with the Owczarzy (2004) monovalent-salt correction, at configurable
Na⁺ and oligo concentration.

The GC-content approximation common in teaching scripts —
`64.9 + 41 × (GC − 16.4) / N` — is sequence-order blind and ignores salt and
strand concentration. On this target it misplaces Tm by a mean of 2.3 °C and up
to 6.1 °C, which is enough to let primers through a 58–62 °C filter whose true
Tm is 66 °C.

## Limitation

**Specificity is not checked.** Candidates are selected on thermodynamics and
structure alone. Before any wet-lab use, screen every oligo against the host
genome, near-neighbour organisms and the target's strain diversity with BLASTn
or Primer-BLAST, and confirm the amplicon region is conserved across circulating
variants. `--report-fasta` writes the oligos in FASTA form for that step.

## References

- SantaLucia J. (1998) *PNAS* 95:1460–1465 — unified nearest-neighbour parameters.
- Owczarzy R. *et al.* (2004) *Biochemistry* 43:3537–3554 — salt correction.
