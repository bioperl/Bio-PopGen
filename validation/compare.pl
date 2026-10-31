#!/usr/bin/env perl
# Join bioperl_values.tsv and popgenome_values.tsv and report differences.
use strict; use warnings; use FindBin;
my $tol = 1e-3;
my %v;
for my $src (qw(bioperl popgenome)) {
    open my $fh, '<', "$FindBin::Bin/${src}_values.tsv" or die "run the $src step first: $!";
    my $hdr = <$fh>;
    while (<$fh>) { chomp; my ($ds, $st, $val) = split /\t/; $v{$ds}{$st}{$src} = $val }
}
# PopGenome fl_f2() multiplies thetaFL, which already carries the (n-1)/n factor,
# by (n-1)/n again. Its F* therefore equals Bio::PopGen F* scaled by
# (pi - eta k^2) / (pi - eta k) with k = (n-1)/n and eta the singleton count.
sub fstar_explained {
    my ($set, $b, $p) = @_;
    my ($pi, $eta, $n) = map { $set->{$_}{bioperl} } qw(pi singletons n_sequences);
    return 0 unless defined $pi && defined $eta && defined $n && $n > 1;
    my $k = ($n - 1) / $n;
    my $pred = $b * ($pi - $eta * $k**2) / ($pi - $eta * $k);
    return abs($pred - $p) <= $tol * (1 + abs($p));
}
my ($n, $bad) = (0, 0);
printf "%-10s %-20s %14s %14s  %s\n", qw(dataset statistic bioperl popgenome status);
for my $ds ( sort keys %v ) {
    for my $st ( sort keys %{ $v{$ds} } ) {
        my ($b, $p) = @{ $v{$ds}{$st} }{qw(bioperl popgenome)};
        my $status;
        if ( $st =~ /^(singletons|external_mutations|n_codons)$/ ) { $status = 'n/a (no PopGenome equivalent)' }
        elsif ( !defined $b || !defined $p || $b eq '' || $p eq 'NA' || $p eq '' ) { $status = 'missing' }
        elsif ( $st =~ /^(singletons|external_mutations|n_codons)$/ ) { $status = 'n/a (no PopGenome equivalent)' }
        elsif ( abs($b - $p) <= $tol * (1 + abs($p)) ) { $status = 'ok' }
        elsif ( $st eq 'fu_li_F_star' && fstar_explained($v{$ds}, $b, $p) ) {
            $status = 'ok after correction (PopGenome applies (n-1)/n to the singleton term twice)';
        }
        else { $status = sprintf 'DIFF %.4g', $b - $p; $bad++ }
        $n++;
        printf "%-10s %-20s %14s %14s  %s\n", $ds, $st, map({ defined $_ ? $_ : 'NA' } $b, $p), $status;
    }
}
print "\n$n comparisons, $bad differences beyond tolerance $tol\n";
exit($bad ? 1 : 0);
