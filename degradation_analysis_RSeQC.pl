use strict;

sub gtf_to_bed12 {
    my($in, $out) = @_;
    
    my @f;
    my $pid;

    open(IN, $in) or die;
    open(OUT, ">$out") or die;
    
    while(<IN>) {
	chomp;
	my @d = split(/\t/, $_);
	
	if($d[2] eq 'exon') {
	    my ($p) = $d[8] =~ /Parent=([^;]*)/;
	    ($p) = $d[8] =~ /transcript_id "([^"]*)"/ unless($p);
	    unless($pid && $p eq $pid) {
		print OUT dump_line(\@f) if(@f);
		@f = ();
		$pid = $p;
	    }
	    push @f, [@d];
	}
    }
    
    print OUT dump_line(\@f) if(@f);

    close IN;
    close OUT;
}

sub dump_line {
    my ($f) = @_;
    my @f = @$f;
    
    return unless(@f);
    
    @f = sort { $a->[3] <=> $b->[3] } @f;

    my $chr = $f[0][0];
    my ($min) = sort { $a <=> $b } map { $_->[3] } @f;
    my ($max) = sort { $b <=> $a } map { $_->[4] } @f;

    my ($p) = $f[0][8] =~ /Parent=([^;]*)/;
    ($p) = $f[0][8] =~ /transcript_id "([^"]*)"/ unless($p);
    my ($d) = $f[0][6];

    return join("\t", $chr, $min-1, $max, $p, 0, $d, $min-1, $max, 0, scalar(@f), join(",", map { $_->[4]-$_->[3]+1 } @f), join(",", map { $_->[3] - $min } @f)), "\n";    
}

my $workingdir = ".";

for(my $i = 0; $i < @ARGV; $i++) {
    my $arg = $ARGV[$i];
    if($arg =~ /-d/) {
	$workingdir = $ARGV[$i+1];
    }
}

my $geneBody_coverage = `which geneBody_coverage.py`;
my $tin = `which tin.py`;

my $bedfile = "$workingdir/reference/annotation.bed12";
my $gtffile = "$workingdir/reference/annotation.gtf";

chomp $geneBody_coverage;
chomp $tin;

die "geneBody_coverage.py not found!" if(! -x $geneBody_coverage);
die "tin.py not found!" if(! -x $tin);

if(! -e $bedfile) {
    die "annotation.gtf not found!" if(! -e $gtffile);
    
    print "creating BED12 file ...\n\n";
    gtf_to_bed12($gtffile, $bedfile);
}

my @bams = <$workingdir."/results/sorted_bam/*.bam">;
die "no BAM file found!\n\nPlease run rup.R first!\n\n" unless(@bams);

mkdir $workingdir."/results/rseqc/" if(! -d $workingdir."/results/rseqc/");

chdir($workingdir."/results/rseqc/");

print "WD: $workingdir\n";
print "RSeQC: $geneBody_coverage\n";
print "RSeQC: $tin\n\n";

my $bamfiles = join(",", <../sorted_bam/*.bam>);

$bedfile = "../../reference/annotation.bed12";

print "BAM: $bamfiles\n\n";

system($geneBody_coverage, "-i", $bamfiles, "-r", $bedfile, "-o", "rseqc") if(-x $geneBody_coverage);

system($tin, "-i", $bamfiles, "-r", $bedfile) if(-x $tin);
 
