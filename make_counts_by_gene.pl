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
    my ($marker, $epName, $epd, $distribute, $tt, $gt_count) = $csv_parser->fields;
    $genes->{$marker}->{$epName}->{gt_count} = 0;
    if($epd){
        $genes->{$marker}->{$epName}->{epd_count}++;
    }
    if($distribute eq 'yes'){
        $genes->{$marker}->{$epName}->{distribute_count}++;
    }
    if($tt eq 'yes'){
        $genes->{$marker}->{$epName}->{tt_count}++;
    }
    if($gt_count){
        $genes->{$marker}->{$epName}->{gt_count} = $gt_count;
    }
}

print join("\t", "MARKER_SYMBOL","EP_PLATE_NAME","EPD_WELL_NAME","EPD_DISTRIBUTE","TARGETED_TRAP","MGI_GT_COUNT"."\n");

	while ( my ( $marker, $marker_data ) = each %{$genes} ) {
	    while ( my ( $epName, $ep_data ) = each %{$marker_data} ) {
		my $epd_count = $ep_data->{epd_count} || 0;
		my $distribute_count = $ep_data->{distribute_count} || 0;
		my $tt_count = $ep_data->{tt_count} || 0;
		my $gt_count = $ep_data->{gt_count} || 0;
		print join("\t",$marker,$epName,$epd_count,$distribute_count,$tt_count,$gt_count)."\n"; 
	    }
	}


