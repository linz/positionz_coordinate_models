#!/usr/bin/perl
use strict;
use POSIX;

die unless @ARGV;
foreach my $f (@ARGV)
{
    my $mtime=(stat($f))[9];
    my $timestr=strftime("%Y-%m-%dT%H:%M:%S",gmtime($mtime));
    open( my $xml, "<$f");
    my $data=join('',<$xml>);
    close($xml);
    $data =~ s/start_date/version_date="$timestr" start_date/;
    open( $xml, ">$f");
    print $xml $data;
    close($xml);
    utime $mtime, $mtime, $f;
}
