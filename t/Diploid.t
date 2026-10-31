# -*-Perl-*- Test Harness script for Bioperl
#
# Fu and Li statistics and haploid_population with diploid genotypes.
# A Bio::PopGen::Population of diploid individuals must give the same
# answers as the array of haploid individuals it expands to.

use strict;

BEGIN {
    use lib '.';
    use Bio::Root::Test;

    test_begin(-tests => 44);

    use_ok('Bio::PopGen::Individual');
    use_ok('Bio::PopGen::Genotype');
    use_ok('Bio::PopGen::Population');
    use_ok('Bio::PopGen::Statistics');
    use_ok('Bio::PopGen::IO');
    use_ok('Bio::PopGen::Simulation::Coalescent');
    use_ok('Bio::PopGen::Simulation::GeneticDrift');
}

sub mk {
    my ($id,%g) = @_;
    my $i = Bio::PopGen::Individual->new(-unique_id => $id);
    for my $m ( sort keys %g ) {
        $i->add_Genotype(Bio::PopGen::Genotype->new(-marker_name => $m,
                                                    -alleles     => $g{$m}));
    }
    return $i;
}

# 3 diploid individuals = 6 chromosomes, 5 sites. Per site allele counts:
#  S1 A5 G1 (G derived singleton -> external)
#  S2 C3 T3 (T derived, count 3   -> internal)
#  S3 T5 A1 (A derived singleton -> external)
#  S4 A4 C2 (C derived, count 2   -> internal)
#  S5 monomorphic
# So n = 6, S = 4, external = 2, singletons = 2.
my @ingroup = (
    mk('i1', S1 => ['A','G'], S2 => ['C','C'], S3 => ['T','T'], S4 => ['A','A'], S5 => ['G','G']),
    mk('i2', S1 => ['A','A'], S2 => ['C','T'], S3 => ['T','T'], S4 => ['A','C'], S5 => ['G','G']),
    mk('i3', S1 => ['A','A'], S2 => ['T','T'], S3 => ['T','A'], S4 => ['A','C'], S5 => ['G','G']),
);
my @outgroup = ( mk('o1', S1 => ['A'], S2 => ['C'], S3 => ['T'], S4 => ['A'], S5 => ['G']) );

my $ingroup  = Bio::PopGen::Population->new(-individuals => \@ingroup);
my $outgroup = Bio::PopGen::Population->new(-individuals => \@outgroup);
my $stats    = Bio::PopGen::Statistics->new();

# haploid_population: one new individual per chromosome
my $hap = $ingroup->haploid_population;
is($hap->get_number_individuals, 6, 'haploid_population has 2 individuals per diploid');
my @hapinds = sort { $a->unique_id cmp $b->unique_id } $hap->get_Individuals;
is_deeply([ map { $_->unique_id } @hapinds ], [qw(i1.1 i1.2 i2.1 i2.2 i3.1 i3.2)],
          'haploid individuals are named <id>.<chromosome>');
is(scalar( grep { my $h = $_; grep { $h == $_ } @ingroup } @hapinds ), 0,
   'haploid individuals are new objects, not the diploid originals');
is_deeply([ map { join('', $_->get_Genotypes(-marker => 'S1')->get_Alleles) } @hapinds ],
          [qw(A G A A A A)], 'each haploid individual carries one allele per site');
is(join('', $hapinds[1]->get_Genotypes(-marker => 'S2')->get_Alleles), 'C',
   'chromosome 2 of i1 carries the second allele of every site');

# derived mutations: population input must match arrayref input
my @from_array = $stats->derived_mutations(\@ingroup, \@outgroup);
my @from_pop   = $stats->derived_mutations($ingroup, $outgroup);
is_deeply(\@from_array, [2, 2], 'derived_mutations(arrayref): 2 external, 2 internal');
is_deeply(\@from_pop,   [2, 2], 'derived_mutations(population) agrees');

is($stats->segregating_sites_count($ingroup), 4, 'segregating sites');
is($stats->singleton_count($ingroup), 2, 'singletons');
is(sprintf("%.4f", $stats->pi($ingroup)), '1.8000', 'pi');
is(sprintf("%.4f", $stats->pi(\@hapinds)), '1.8000', 'pi on haploid individuals');

# Fu and Li statistics. Expected values come from the *_counts forms
# with n = 6 sequences, S = 4, external = 2, singletons = 2, pi = 1.8.
my $D  = $stats->fu_and_li_D_counts(6, 4, 2);
my $Ds = $stats->fu_and_li_D_star_counts(6, 4, 2);
my $F  = $stats->fu_and_li_F_counts(6, 1.8, 4, 2);
my $Fs = $stats->fu_and_li_F_star_counts(6, 1.8, 4, 2);
my $T  = $stats->tajima_D_counts(6, 4, 1.8);
ok(defined $D, 'fu_and_li_D_counts defined for n=6');

for my $case ( ['population', $ingroup, $outgroup],
               ['haploid arrayref', \@hapinds, \@outgroup] ) {
    my ($label, $in, $out) = @$case;
    is(sprintf("%.5f", $stats->fu_and_li_D($in, $out)),  sprintf("%.5f", $D),  "fu_and_li_D ($label)");
    is(sprintf("%.5f", $stats->fu_and_li_D($in, 2)),     sprintf("%.5f", $D),  "fu_and_li_D with external count ($label)");
    is(sprintf("%.5f", $stats->fu_and_li_D_star($in)),   sprintf("%.5f", $Ds), "fu_and_li_D_star ($label)");
    is(sprintf("%.5f", $stats->fu_and_li_F($in, $out)),  sprintf("%.5f", $F),  "fu_and_li_F ($label)");
    is(sprintf("%.5f", $stats->fu_and_li_F_star($in)),   sprintf("%.5f", $Fs), "fu_and_li_F_star ($label)");
    is(sprintf("%.5f", $stats->tajima_D($in)),           sprintf("%.5f", $T),  "tajima_D ($label)");
    is(sprintf("%.5f", $stats->theta($in)),              sprintf("%.5f", 4 / (1 + 1/2 + 1/3 + 1/4 + 1/5)), "theta ($label)");
}

# add_Individual keeps only IndividualI objects
{
    my $pop = Bio::PopGen::Population->new;
    my $w = 0;
    local $SIG{__WARN__} = sub { $w++ };
    $pop->add_Individual($ingroup[0], 'not an individual', undef);
    is($pop->get_number_individuals, 1, 'add_Individual drops non-individuals');
    is($w, 1, 'add_Individual warned once');
}

# composite_LD accepts an arrayref
{
    my %ld = $stats->composite_LD(\@ingroup);
    ok(exists $ld{S1}, 'composite_LD(arrayref) returned results');
}

# Coalescent maxcount is honored
{
    my $sim = Bio::PopGen::Simulation::Coalescent->new(-sample_size => 4, -maxcount => 2);
    is($sim->maxcount, 2, 'maxcount getter');
    my $n = 0;
    while ( $n < 10 && $sim->next_tree ) { $n++ }
    is($n, 2, 'next_tree stops after maxcount trees');
}

# GeneticDrift can return a Population
{
    my $sim = Bio::PopGen::Simulation::GeneticDrift->new(-popsize => 50,
                                                          -alleles => { A => 0.5, B => 0.5 });
    my $pop = $sim->next_generation('population');
    isa_ok($pop, 'Bio::PopGen::PopulationI');
    my %af = $pop->get_Marker('Marker1')->get_Allele_Frequencies;
    is(sprintf("%.3f", ($af{A} || 0) + ($af{B} || 0)), '1.000', 'frequencies sum to 1');
}

# format guessing
is(Bio::PopGen::IO->_guess_format('x.dat'), 'csv', 'guess .dat');
is(Bio::PopGen::IO->_guess_format('x.prettybase'), 'prettybase', 'guess .prettybase');
is(Bio::PopGen::IO->_guess_format('x.phase'), 'phase', 'guess .phase');
is(Bio::PopGen::IO->_guess_format('x.hapmap'), 'hapmap', 'guess .hapmap');
