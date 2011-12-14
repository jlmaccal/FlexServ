#! /usr/bin/perl
use strict;
#need 2 files: first is CA atom pdb (trimmed)
#second is the all_atom pdb, possibly not trimmed

#Assign the beta factors of the CA pdb to the other
#Trim to same size (for comparison with different structures

my ($CA,$pdb) = @ARGV;

open (CA,"<$CA");
  
#Read CA atom structure wt beta factors
my (@crystal,$res,$allxal,$beta,@beta);
my $count = 0;
while (<CA>) {
  if (/^ATOM/){ 
   $res = substr($_,17,3);
   $beta = substr($_,60,6);
   push(@beta,$beta);
   push(@crystal,$res);
   $allxal .= $res;
   $count++;
  }
} 
close (CA);

#Read ALl atom structure
my (@template,$res,$alltem);
open (PDB,"<$pdb");
while (<PDB>) {
  if ((/^ATOM/)&&(/ CA /)){
   $res = substr($_,17,3);
   push(@template,$res);
   $alltem .= $res;
  }
}
close (PDB);

#Align
my $alltem =~ s/$allxal/999/g;
my @parts = split("999",$alltem);
my @before = split("",$parts[0]);
my @after = split("",$parts[1]);
my $res_before = int(($#before+1)/3);
my $res_after = int(($#after+1)/3);

open (IN,"<$pdb");
my ($bef,$aft);
my $new_count = 0;
my $ini_res = 0;
while (<IN>) {
  if (/^ATOM/){
   $res = substr($_,22,4); #if get 21,5 then if there's chain information present we do badly
   if (!$ini_res) {
	$ini_res = $res;
   }
   if ($res-$ini_res+1 != $new_count){
   	$new_count++;
   }
   next if ($new_count <= $res_before); 
   last if ($new_count > $count-$res_after); 
   $bef = substr($_,0,20);
   $res = sprintf("%6d",$new_count);
   $aft = substr($_,26,34);
   $beta = $beta[$new_count-1];
   print "$bef$res$aft$beta\n";
  }
}
close (IN);


