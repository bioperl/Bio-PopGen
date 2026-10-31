#!/usr/bin/env perl
# Write the small test datasets as FASTA alignments for PopGenome and
# record the Bio::PopGen values for the same data in bioperl_values.tsv.
# Run from the validation/ directory with a perl that has BioPerl and
# Bio::MolEvol::CodonModel on its include path:
#   perl -I../lib prepare_data.pl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Bio::PopGen::IO;
use Bio::PopGen::Individual;
use Bio::PopGen::Genotype;
use Bio::PopGen::Population;
use Bio::PopGen::Statistics;
use Bio::PopGen::Utilities;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;

my $datadir = "$FindBin::Bin/data";
mkdir $datadir unless -d $datadir;
my $stats   = Bio::PopGen::Statistics->new;
my @rows;   # dataset, statistic, value

sub mk {
    my ($id, %g) = @_;
    my $i = Bio::PopGen::Individual->new(-unique_id => $id);
    $i->add_Genotype(Bio::PopGen::Genotype->new(-marker_name => $_, -alleles => $g{$_}))
        for sort keys %g;
    return $i;
}

# write one FASTA per dataset: ingroup haplotypes then outgroup haplotypes.
# Each marker becomes one alignment column, markers in sorted order. The
# alignment is padded with monomorphic columns because PopGenome cannot
# read an alignment of only a few columns. Monomorphic columns change no
# statistic here except the site count.
my $PAD = 30;
sub write_fasta {
    my ($name, $in, $out) = @_;
    mkdir "$datadir/$name" unless -d "$datadir/$name";
    my @markers = sort { $a cmp $b } $in->[0]->get_marker_names;
    open my $fh, '>', "$datadir/$name/$name.fa" or die $!;
    for my $ind ( @$in, @$out ) {
        my $seq = join '', map {
            my ($g) = $ind->get_Genotypes(-marker => $_);
            my @a = $g ? $g->get_Alleles(1) : ('N');
            @a == 1 ? $a[0] : die "$name: ".$ind->unique_id." is not haploid at $_";
        } @markers;
        print $fh ">", $ind->unique_id, "\n$seq", "A" x $PAD, "\n";
    }
    close $fh;
    return scalar(@markers) + $PAD;
}

sub record_stats {
    my ($name, $in, $out, $nsites) = @_;
    my $pop  = Bio::PopGen::Population->new(-individuals => $in);
    my $opop = Bio::PopGen::Population->new(-individuals => $out);
    my ($ext, $int) = $stats->derived_mutations($pop, $opop);
    my %v = (
        n_sequences       => $stats->_sample_size($pop),
        n_sites           => $nsites,
        segregating_sites => $stats->segregating_sites_count($pop),
        singletons        => $stats->singleton_count($pop),
        external_mutations=> $ext || 0,
        pi                => $stats->pi($pop),
        theta             => $stats->theta($pop),
        tajima_D          => $stats->tajima_D($pop),
        fu_li_D           => $stats->fu_and_li_D($pop, $opop),
        fu_li_F           => $stats->fu_and_li_F($pop, $opop),
        fu_li_D_star      => $stats->fu_and_li_D_star($pop),
        fu_li_F_star      => $stats->fu_and_li_F_star($pop),
    );
    push @rows, [ $name, $_, $v{$_} ] for sort keys %v;
}

# ---- dataset 1: t/data/popstats.prettybase (haploid, 5 ingroup + 1 outgroup)
{
    my $io = Bio::PopGen::IO->new(-format => 'prettybase', -no_header => 1,
                                  -file => "$FindBin::Bin/../t/data/popstats.prettybase");
    my (@in, @out);
    while ( my $ind = $io->next_individual ) {
        push @{ $ind->unique_id =~ /out/ ? \@out : \@in }, $ind;
    }
    @in = sort { $a->unique_id cmp $b->unique_id } @in;
    my $n = write_fasta('popstats', \@in, \@out);
    record_stats('popstats', \@in, \@out, $n);
}

# ---- dataset 2: the diploid fixture from t/Diploid.t, as its haploid expansion
{
    my @dip = (
        mk('i1', S1 => ['A','G'], S2 => ['C','C'], S3 => ['T','T'], S4 => ['A','A'], S5 => ['G','G']),
        mk('i2', S1 => ['A','A'], S2 => ['C','T'], S3 => ['T','T'], S4 => ['A','C'], S5 => ['G','G']),
        mk('i3', S1 => ['A','A'], S2 => ['T','T'], S3 => ['T','A'], S4 => ['A','C'], S5 => ['G','G']),
    );
    my @out = ( mk('o1', S1 => ['A'], S2 => ['C'], S3 => ['T'], S4 => ['A'], S5 => ['G']) );
    my $dpop = Bio::PopGen::Population->new(-individuals => \@dip);
    my @hap  = sort { $a->unique_id cmp $b->unique_id } $dpop->haploid_population->get_Individuals;
    my $n = write_fasta('diploid', \@hap, \@out);
    # statistics are computed on the diploid Population; the FASTA holds its haplotypes
    record_stats('diploid', \@dip, \@out, $n);
}

# ---- datasets 3 and 4: McDonald-Kreitman on the two test alignments
# (sim ingroup, mel outgroup). The alignments have long runs of N in most
# sim sequences. Bio::PopGen skips a codon when fewer than two ingroup
# alleles are called and PopGenome counts any site with at least one call,
# so the two cannot be compared on the raw files. Keep mel and the two sim
# sequences with the fewest N, then drop every codon in which any kept
# sequence has a character other than A, C, G or T, and every codon with a
# position that has three alleles.
for my $gene ( qw(CG2865 CG11099) ) {
    my $name = "mk_$gene";
    my $raw = Bio::AlignIO->new(-format => 'fasta',
                                -file => "$FindBin::Bin/../t/data/$gene.fasaln")->next_aln;
    my ($mel) = grep { $_->display_id =~ /^mel/ } $raw->each_seq;
    my @sim = sort { ($a->seq =~ tr/Nn//) <=> ($b->seq =~ tr/Nn//) }
              grep { $_->display_id =~ /^sim/ } $raw->each_seq;
    my @keep = ($mel, @sim[0,1]);
    my @codons;
    CODON: for ( my $i = 0; $i + 3 <= $raw->length; $i += 3 ) {
        for my $s ( @keep ) {
            next CODON if substr($s->seq, $i, 3) =~ /[^ACGTacgt]/;
        }
        # PopGenome analyses biallelic sites only; Bio::PopGen counts a
        # three-allele codon position as one fixed and one polymorphic change.
        # Skip such codons so both tools see the same sites.
        for my $p ( 0..2 ) {
            my %b = map { uc(substr($_->seq, $i + $p, 1)) => 1 } @keep;
            next CODON if keys %b > 2;
        }
        push @codons, $i;
    }
    my $aln = Bio::SimpleAlign->new;
    for my $s ( @keep ) {
        (my $id = $s->display_id) =~ s/[^\w.+-]/_/g;
        $aln->add_seq(Bio::LocatableSeq->new(-id => $id, -alphabet => 'dna',
            -seq => join('', map { substr($s->seq, $_, 3) } @codons)));
    }
    my $pop = Bio::PopGen::Utilities->aln_to_population(-alignment => $aln, -site_model => 'codon',
                                                        -include_monomorphic => 1);
    my (@in, @out);
    for my $ind ( $pop->get_Individuals ) {
        push @in,  $ind if $ind->unique_id =~ /^sim/;
        push @out, $ind if $ind->unique_id =~ /^mel/;
    }
    my @c = $stats->mcdonald_kreitman(-ingroup => \@in, -outgroup => \@out);
    # PopGenome wants the alignment in its own directory and the GFF (same
    # file name) in another; the whole alignment is one CDS in frame 0
    mkdir "$datadir/$name" unless -d "$datadir/$name";
    mkdir "$datadir/${name}_gff" unless -d "$datadir/${name}_gff";
    open my $fh, '>', "$datadir/$name/$gene.fa" or die $!;
    print $fh ">", $_->display_id, "\n", $_->seq, "\n" for $aln->each_seq;
    close $fh;
    open my $gff, '>', "$datadir/${name}_gff/$gene.gff" or die $!;
    print $gff join("\t", "$gene.fa", 'test', 'gene', 1, $aln->length, '.', '+', '.', 'ID=g1'), "\n";
    print $gff join("\t", "$gene.fa", 'test', 'CDS',  1, $aln->length, '.', '+', 0,   'Parent=g1'), "\n";
    close $gff;
    push @rows, [ $name, 'n_codons', scalar @codons ],
                [ $name, 'poly_N',  $c[0] ], [ $name, 'fixed_N', $c[1] ],
                [ $name, 'poly_S',  $c[2] ], [ $name, 'fixed_S', $c[3] ];
}

open my $out, '>', "$FindBin::Bin/bioperl_values.tsv" or die $!;
print $out join("\t", qw(dataset statistic bioperl)), "\n";
print $out join("\t", @$_), "\n" for @rows;
close $out;
print "wrote $FindBin::Bin/bioperl_values.tsv (", scalar @rows, " values)\n";
