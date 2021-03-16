#This Perl code identifies the first base before the 5' of each read on the opposite strand (indirect mode)
#!/usr/bin/perl
use strict;
use warnings;

my $usage = "usage perl tmp.pl <filename>";
my $filename = $ARGV[0];
die $usage unless $filename;

open(FH,"<$filename");
my @lines = <FH>;

foreach my $line(@lines){
  chomp $line;
  my @elements = split(/\t/,$line);

  if($elements[5]=~s/\+/-/){
    $elements[2]=$elements[1];
    $elements[1]--;
  }
  elsif($elements[5]=~s/\-/+/){
    $elements[1]=$elements[2];
    $elements[2]++;
  }

  my $out = join("\t", @elements);
  print $out."\n";




} 
