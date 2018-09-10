use strict;

use Bio::Root::Test;

use Bio::PopGen::Genotype;
use Bio::Tree::AlleleNode;

my $allele_node = Bio::Tree::AlleleNode->new();
$allele_node->add_Genotype(
  Bio::PopGen::Genotype->new(
    -marker_name => 'm1',
    -alleles     => [0]
  )
);
$allele_node->add_Genotype(
  Bio::PopGen::Genotype->new(
    -marker_name => 'm3',
    -alleles     => [ 1, 1 ]
  )
);
$allele_node->add_Genotype(
  Bio::PopGen::Genotype->new(
    -marker_name => 'm4',
    -alleles     => [ 0, 4 ]
  )
);
ok($allele_node);
my @mkrs = $allele_node->get_marker_names;

is( @mkrs, 3 );
my ($m3) = $allele_node->get_Genotypes( -marker => 'm3' );
is( $m3->get_Alleles, 2 );
my ($a1) = $allele_node->get_Genotypes( -marker => 'm1' )->get_Alleles;
is( $a1, 0 );

my ( $a2, $a3 ) = $allele_node->get_Genotypes( -marker => 'm4' )->get_Alleles;
is( $a2, 0 );
is( $a3, 4 );

done_testing();
