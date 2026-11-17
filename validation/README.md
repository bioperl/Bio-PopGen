# Validation of Bio::PopGen statistics against PopGenome

This directory checks the numbers Bio::PopGen produces against an
independent implementation, the R package PopGenome 2.7.5 (Pfeifer et
al. 2014, Mol Biol Evol 31:1929), on the small datasets used by the test
suite.

## Setup

    cd validation
    pixi install          # R + PopGenome from conda-forge, see pixi.toml

The Bio::PopGen side needs a perl with BioPerl and Bio::MolEvol::CodonModel
on its include path (Bio::PopGen::Statistics loads the latter).

## Run

    perl prepare_data.pl      # writes data/ and bioperl_values.tsv
    pixi run reference        # writes popgenome_values.tsv
    perl compare.pl           # joins the two and reports differences

or `pixi run all` after `prepare_data.pl`. `compare.pl` exits non-zero
if any value differs by more than 1e-3 (relative). `data/` and the two
TSV files are generated and not committed.

## Datasets

| dataset    | source                               | statistics                                                    |
|------------|--------------------------------------|---------------------------------------------------------------|
| popstats   | t/data/popstats.prettybase (haploid) | S, pi, theta, Tajima D, Fu and Li D, F, D*, F*                |
| diploid    | the fixture in t/Diploid.t           | same; Bio::PopGen works on the diploid Population, PopGenome on its six haplotypes |
| mk_CG2865  | t/data/CG2865.fasaln                 | McDonald-Kreitman 2x2 counts, sim ingroup vs mel outgroup     |
| mk_CG11099 | t/data/CG11099.fasaln                | same                                                          |

The first two are written as FASTA alignments with one column per marker,
padded with 30 monomorphic columns because PopGenome cannot read an
alignment of only a few columns. Monomorphic columns change no compared
statistic.

For the MK alignments, `prepare_data.pl` keeps mel and the two sim
sequences with the fewest N, then drops every codon in which a kept
sequence has a character other than A, C, G or T, and every codon with a
three-allele position. Reasons:

- Bio::PopGen skips a codon when fewer than two ingroup alleles are
  called. PopGenome counts any site with at least one call. On the raw
  files, which have long runs of N in most sim sequences, the two count
  different sites.
- PopGenome analyses biallelic sites only and lists others under
  `polyallelic.sites`. Bio::PopGen counts a codon position with three
  alleles (two in the ingroup, a third in the outgroup) as one fixed and
  one polymorphic change. Before this filter CG11099 differed by exactly
  that one codon (codon 132, GTA/GTC vs GTT) in fixed_S and poly_S.

## Result (2026-10-30)

All compared values agree to at least six digits, with one documented
exception:

- Fu and Li F*: Bio::PopGen uses the numerator pi - eta_s (n-1)/n and the
  Simonsen et al. 1995 variance, and reproduces the values of
  `fu_and_li_F_star_counts` exactly. PopGenome's `init_coef` has the same
  variance terms, but its `fl_f2` multiplies `thetaFL`, which already
  carries the (n-1)/n factor (see `ap_pop`), by (n-1)/n again.
  `compare.pl` checks that PopGenome's F* equals the Bio::PopGen value
  scaled by (pi - eta k^2)/(pi - eta k), k = (n-1)/n, and reports it as
  "ok after correction". This is reported as a PopGenome issue, not a
  Bio::PopGen one.

Statistics with no PopGenome equivalent (singleton count, external
mutation count, codon count) are listed as "n/a".

## Definitions that had to match

- Bio::PopGen `pi` and `theta` without a site count are totals over the
  region, as are PopGenome `Pi` (from `diversity.stats`) and
  `theta_Watterson`.
- Fu and Li D and F need an outgroup. PopGenome reports D* and F* under
  the same column names when no outgroup is set, which is how the script
  obtains them.
- For an alignment read with a GFF, PopGenome classifies sites as
  synonymous or not while reading. `set.synnonsyn` is only for SNP data
  and fails on alignment input.
