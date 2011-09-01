#!/usr/bin/perl -w
use strict;
use warnings;

use Text::CSV_XS;

my $csv_parser = Text::CSV_XS->new;

my $first_line = 1;
my $genes = {};
while (<>){
    if ($first_line){
        $first_line = 0;
        next;
    }
    chomp;
    $csv_parser->parse($_);
	# counting loop
    my ($marker, $epName, $epd, $distribute, $tt, $gt_count) = $csv_parser->fields;
    $genes->{$marker}->{gt_count} = 0;
    if($epd){
        $genes->{$marker}->{epd_count}++;
    }
    if($distribute eq 'yes'){
        $genes->{$marker}->{distribute_count}++;
    }
    if($tt eq 'yes'){
        $genes->{$marker}->{tt_count}++;
    }
    if($gt_count){
        $genes->{$marker}->{gt_count} = $gt_count;
    }
	$genes->{$marker}->{ep}->{$epName}++;
}

	# fetching and printing loop
print join("\t", "MARKER_SYMBOL","EP_PLATE_NAME","EPD_WELL_NAME","EPD_DISTRIBUTE","TARGETED_TRAP","MGI_GT_COUNT"."\n");
while ( my ( $marker, $marker_data ) = each %{$genes} ) {
        my $epd_count = $marker_data->{epd_count} || 0;
        my $distribute_count = $marker_data->{distribute_count} || 0;
        my $tt_count = $marker_data->{tt_count} || 0;
        my $gt_count = $marker_data->{gt_count} || 0;
		my $eps = join q{,}, sort keys %{ $marker_data->{ep} };
        print join( "\t", $marker,$eps,$epd_count,$distribute_count,$tt_count,$gt_count) . "\n"; 
}

