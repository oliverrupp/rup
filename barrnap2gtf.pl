use strict;

my $c = 1;

while(<>) {
	chomp;
	next if(/^#/);
	my @x = split(/\t/, $_);

	my ($t) = $x[8] =~ /Name=([^;]*)/;
	my $id = "rRNA_".$c."_$t";
	#$x[2] = "mRNA";
	#$x[8] = "ID=".$id;
	#print join("\t", @x), "\n";

	$x[2] = "exon";
	$x[8] = "transcript_id \"".$id."\"; gene_id \"".$id."\"";
	print join("\t", @x), "\n";

	$c++;
}
