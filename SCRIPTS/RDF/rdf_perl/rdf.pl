#!/usr/bin/perl

# An example call: ./rdf.pl OW OW HISTORY > output_file &

##############################    Options    ###################################
#
# To skip frames, set $read_every to the interval you want to skip. e.g., 
# $read_every = 1 reads every frame, $read_every = 2 reads every other frame
#
my $read_every = 1;
#my $read_every = 1;
#my $stop_time = 10.0; #only bin RDF up to this simulation time
my $stop_time = undef; #if undef, read the whole HISTORY file
#my $start_time =  10.0; #begin benning when simulation time >= $start_time, or
my $start_time = undef; #if undef, read whole history file
#
# Verbose = 1 prints each frame number
my $verbose = 0;
#
my $rmax = 11.000;
my $delta = 0.005;
#
#
################################################################################

#
# 5 November 2011 - Greg Medders - gmedders@ucsd.edu
#
# Added option to read every set number of frames.
#
# 11 August 2001 - Jeffrey Sung jcsung@ucsd.edu
# This script no longer requires HistoryFile.pm
#
# 4 August 2011 - Greg Medders - gmedders@ucsd.edu
# This script computes the radial distribution function for two given atom
# types using dlpoly2 or dlpoly4 HISTORY files. If you want the RDF over
# multiple simulations (multiple history files), you simply have to list the
# HISTORY_FILEs in the command line. 
#
# rdf.pl depends on HistoryFile.pm. Include it in the directory in which you 
# call rdf.pl
#
# By default, the RDF is collected over a pair separation from 0.05 - 9.7A.
# This script can take a long time to run, especially if you have a huge box
# or if you have tons many, many frames in a history file or many histories.
# For water, I generally average over 100 frames with 256 molecules/box.
# Depending on the size/orderedness of your system, you may need more or 
# fewer frames.
#
# Note, this is a modification of a script by Volodymyr

use strict;
use POSIX;
use Math::Trig;

use FindBin;
use lib $FindBin::Bin;

#use HistoryFile;

die "usage: rdf.pl TYPE1 TYPE2 HISTORY_FILE ... HISTORY_FILE\n"
    unless defined(@ARGV) and scalar @ARGV > 2;

my $type1 = shift @ARGV;
my $type2 = shift @ARGV;

my $oldfh = select(STDERR); $| = 1; select($oldfh);

my $invvolm = 0.0; # sum of 1/volm
my $nframes = 0;
my $pairs_per_frame = -1;

my %bins; # bins are of width $delta : [0, delta), [delta, 2*delta), ...

print STDERR "<> computing RDF for atoms '$type1' and '$type2'\n";
print STDERR "<> rmax = $rmax, bin size = $delta\n";

my $filename="HISTORY";
my $natom=0;
my $lineno=0;
my $frame=undef;
my $this_ppf;

foreach (@ARGV) {
    print STDERR "* processing '$_' ...";
#    my $traj = &HistoryFile::load($_);

    #BEGIN READ HERE

    $filename=$_;
    $frame=undef;
    
    open FH, $filename
	or die " ** Error ** : could not open '$filename' for reading: $!\n";

    $natom = 0;
    my $every_frame = 0;
    my $okay = 0;
    $lineno = 0;

    while (<FH>){
	++$lineno;
	if (/timestep/){
	    ++$every_frame;
	    if($every_frame % $read_every == 0) {
		$okay = 1;
	    }else{
		$okay = 0;
	    }
	    if($okay){
		if (defined($frame)){
#		    die " ** Error ** : '$filename': incomplete frame ending at line $lineno,\n", $natom," != ", scalar @{$frame->{ATOMS}}
#		    unless $natom  == scalar @{$frame->{ATOMS}};
		    if($verbose){
			my $frame_number = $every_frame - $read_every;
			print STDERR "\nProcessing Frame #$frame_number";
		    }
		    if((not defined($start_time)) 
			|| ($frame->{TIME} >= $start_time)){
			$this_ppf = &process_frame($frame);
			if ($pairs_per_frame < 0) {
			    $pairs_per_frame = $this_ppf;
			} else {
			    die unless $pairs_per_frame == $this_ppf; 
# this constraint won't work with GCMC (fluctuating number of pairs per frame).
# Can lift this problem by makin $normfac vary for each frame
    			}
			$nframes++;
		    }
		}

		my @toks = split;
		if(defined $stop_time){
		    if($toks[5] > $stop_time){
			last;
		    }
		}

		$frame = {};
		die " ** Error ** : '$filename': expected six tokens in line $lineno\n"
	    	unless scalar @toks >= 6;#dlpoly4 has 7 tokens per line

		$frame->{STEP}   = $toks[1];
		$natom           = $toks[2];
		$frame->{LEVCFG} = $toks[3];
		$frame->{IMCON}  = $toks[4];
		$frame->{TIME}   = $toks[5];


		die " ** FIXME **\n" unless $frame->{LEVCFG} == 0;

		my @cell = ();
		$frame->{CELL} = \@cell;

		my @box;
		$box[0] = <FH>;
		$box[1] = <FH>;
		$box[2] = <FH>;	

		$lineno += 3;


		for (my $n = 0; $n < 3; ++$n) {
	    	    $box[$n] =~ s/^\s+//;
    		    chomp($box[$n]);
		    my @dim = split(/\s+/,$box[$n]);
		    push @{$frame->{CELL}}, $dim[0];
		    push @{$frame->{CELL}}, $dim[1];
		    push @{$frame->{CELL}}, $dim[2];
		}

		if($frame->{IMCON}==3){
		    #invert matrix
		    $frame->{RCELL}->[0] = $frame->{CELL}->[4]*$frame->{CELL}->[8]
		                          -$frame->{CELL}->[5]*$frame->{CELL}->[7];
		    $frame->{RCELL}->[1] = $frame->{CELL}->[2]*$frame->{CELL}->[7]
		                          -$frame->{CELL}->[1]*$frame->{CELL}->[8];
		    $frame->{RCELL}->[2] = $frame->{CELL}->[1]*$frame->{CELL}->[5]
		                          -$frame->{CELL}->[2]*$frame->{CELL}->[4];
		    $frame->{RCELL}->[3] = $frame->{CELL}->[5]*$frame->{CELL}->[6]
		                          -$frame->{CELL}->[3]*$frame->{CELL}->[8];
		    $frame->{RCELL}->[4] = $frame->{CELL}->[0]*$frame->{CELL}->[8]
		                          -$frame->{CELL}->[2]*$frame->{CELL}->[6];
		    $frame->{RCELL}->[5] = $frame->{CELL}->[2]*$frame->{CELL}->[3]
		                          -$frame->{CELL}->[0]*$frame->{CELL}->[5];
		    $frame->{RCELL}->[6] = $frame->{CELL}->[3]*$frame->{CELL}->[7]
		                          -$frame->{CELL}->[4]*$frame->{CELL}->[6];
		    $frame->{RCELL}->[7] = $frame->{CELL}->[1]*$frame->{CELL}->[6]
		                          -$frame->{CELL}->[0]*$frame->{CELL}->[7];
		    $frame->{RCELL}->[8] = $frame->{CELL}->[0]*$frame->{CELL}->[4]
		                          -$frame->{CELL}->[1]*$frame->{CELL}->[3];

		    my $determinant = $frame->{CELL}->[0]*$frame->{RCELL}->[0]
		                    + $frame->{CELL}->[3]*$frame->{RCELL}->[1]
				    + $frame->{CELL}->[6]*$frame->{RCELL}->[2];
		    my $detinv=0;
		    if (abs($determinant)>0){
			$detinv=1/$determinant;
		    }
		    else{
			die " ** Error ** : Singular matrix!!!!  (Determinant 0)\n";
		    }
		    for (my $count = 0; $count < 9; $count++){
			$frame->{RCELL}->[$count]*=$detinv;
		    }
		    #inversion complete
		}
    		my @atoms = ();
		$frame->{ATOMS} = \@atoms;
	    }
	} elsif (($okay) and (/^(\w+)\s+(\d+)\s+/)){
	    my $atom_name = $1;
	    my $atom_index = $2;

	    my $crd = <FH>;
    	    if ($crd) {
		++$lineno;

		$crd =~ s/^\s+//;
 		chomp($crd);
		my @toks = split(/\s+/,$crd);
		my @pos = ();

		push @pos, $toks[0];
		push @pos, $toks[1];
		push @pos, $toks[2];

	    	my $atom = {};
		$atom->{NAME} = $atom_name;
    		$atom->{INDEX} = $atom_index;
		$atom->{POS} = \@pos;

	    	die unless defined($frame) and defined($frame->{ATOMS});
		if($atom->{NAME} =~ m/$type1/i || $atom->{NAME} =~ m/$type2/i){
		    push @{$frame->{ATOMS}}, $atom;
		}
    	    }
	} elsif ($okay and $lineno < 3) {#this is okay: HISTORY_PIMD has two lines
	    #before "timestep".
	} elsif($okay) {
	    print STDERR $_;
	    die "** Error ** : unexpected input in '$filename', line $lineno\n";
	}
    }

    if (($okay) and defined($frame) and ((not defined($start_time)) || ($frame->{TIME} >= $start_time))) {
#	die " ** Error ** : '$filename': incomplete frame ending at line $lineno\n"
#	unless $natom == scalar @{$frame->{ATOMS}};
	if($verbose){
	    my $frame_number = $every_frame;
	    print STDERR "\nProcessing Frame #$frame_number";
	}
	$this_ppf = &process_frame($frame);
	if ($pairs_per_frame < 0) {
	    $pairs_per_frame = $this_ppf;
	} else {
	    die unless $pairs_per_frame == $this_ppf; # this constraint won't
	    # work with GCMC (fluctuating number of pairs per frame). Can
	    # lift this problem by makin $normfac vary for each frame
	}
	$nframes++;	
    }

    print STDERR ' ', $nframes, ' frames ...';
    print STDERR ' ', ($pairs_per_frame*($nframes)), " samples\n";
}

exit(0) unless $nframes > 0;

my $normfac = 1.0/$invvolm/$pairs_per_frame;

my @sortedbins = sort {$a <=> $b} keys(%bins);
foreach (@sortedbins) {
    my $r = $delta*($_ + 0.5);
    last if ($r + $delta > $rmax);

    my $r_lo = $r - $delta/2.0;
    my $r_hi = $r + $delta/2.0;

    my $dvol = (4.0/3.0)*&pi()*($r_hi**3 - $r_lo**3);
    my $rdf = $normfac*$bins{$_}/$dvol;

    printf "%16.12f %16.12f\n", $r, $rdf;
}

################################################################################

   sub process_frame {
    my $frame = shift;

    my $cell = $frame->{CELL};
    my $volm = 0;
# determine the cell volume
    if($cell->[1] == 0.0 && $cell->[2] == 0.0 &&
    	$cell->[3] == 0.0 && $cell->[5] == 0.0 &&
	$cell->[6] == 0.0 && $cell->[7] == 0.0){
	$volm = $cell->[0]*$cell->[4]*$cell->[8];
    } else { #if not orthorhombic, then use general case: [a] . [b] x [c]
	my $cross_x =  $cell->[4]*$cell->[8] - $cell->[7]*$cell->[5];
	my $cross_y = -$cell->[1]*$cell->[8] + $cell->[7]*$cell->[2];
	my $cross_z =  $cell->[1]*$cell->[5] - $cell->[4]*$cell->[2];
	$volm = $cell->[0]*$cross_x + $cell->[3]*$cross_y
                + $cell->[6]*$cross_z;	
    }
    $invvolm += 1.0/$volm;

    my $npairs = 0;
    my $natom = scalar @{$frame->{ATOMS}};

    for (my $i = 0; $i < $natom; ++$i) {
        for (my $j = $i + 1; $j < $natom; ++$j) {
            my $atomi = $frame->{ATOMS}->[$i];
            my $atomj = $frame->{ATOMS}->[$j];

            if (($atomi->{NAME} eq $type1 and $atomj->{NAME} eq $type2)
             or ($atomi->{NAME} eq $type2 and $atomj->{NAME} eq $type1)) {

		my $rij=-1;
		if($frame->{IMCON} == 3){
			$rij = &distance3($cell, $atomi->{POS}, $atomj->{POS}, $frame->{RCELL});
		}
		else{
			$rij = &distance1($cell, $atomi->{POS}, $atomj->{POS});
		}

                if ($rij < $rmax) {
                    my $binid = &floor($rij/$delta);
                    if (exists($bins{$binid})) {
                        $bins{$binid} += 1;
                    } else {
                        $bins{$binid} = 0;
                    }
                }

                $npairs++;
            }
        }
    }

    return $npairs;
}

################################################################################

sub distance1 {
    my ($cell, $r1, $r2) = @_;

    my ($dx, $dy, $dz) =
       ($r1->[0] - $r2->[0], $r1->[1] - $r2->[1], $r1->[2] - $r2->[2]);

    $dx -= $cell->[0]*&floor($dx/$cell->[0] + 0.5);
    $dy -= $cell->[4]*&floor($dy/$cell->[4] + 0.5);
    $dz -= $cell->[8]*&floor($dz/$cell->[8] + 0.5);

    return sqrt($dx*$dx + $dy*$dy + $dz*$dz);
}

################################################################################

sub distance3 {
# Calculate the distance with imcon = 3

    my ($cell, $r1, $r2, $rcell) = @_;

    my ($dx, $dy, $dz) =
       ($r1->[0] - $r2->[0], $r1->[1] - $r2->[1], $r1->[2] - $r2->[2]);

    my ($ssx,$ssy,$ssz)=($rcell->[0]*$dx+$rcell->[3]*$dy+$rcell->[6]*$dz,
                         $rcell->[1]*$dx+$rcell->[4]*$dy+$rcell->[7]*$dz,
                         $rcell->[2]*$dx+$rcell->[5]*$dy+$rcell->[8]*$dz);

    my ($xss,$yss,$zss)=($ssx-&floor($ssx+0.5),$ssy-&floor($ssy+0.5),$ssz-&floor($ssz+0.5));

    $dx = $cell->[0]*$xss+$cell->[3]*$yss+$cell->[6]*$zss;
    $dy = $cell->[1]*$xss+$cell->[4]*$yss+$cell->[7]*$zss;
    $dz = $cell->[2]*$xss+$cell->[5]*$yss+$cell->[8]*$zss;

    return sqrt($dx*$dx + $dy*$dy + $dz*$dz);

}

################################################################################
