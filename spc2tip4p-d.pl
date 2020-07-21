#!/usr/bin/perl
# converts PDB file with 3-center water (SPC) to 
# Ivan’s tip4pd.idp water
# requires atoms in order O H H

if ($#{ARGV} < 0) {
  die("at least input file name is required!\n");
}
open(IN, $ARGV[0]) ||
  die("cannot open file $ARGV[0] -- $!\n");

# numbers to start counting
$natom = 1;
$nres = 1;

if ($#{ARGV} >= 2) {
  $natom = $ARGV[1];
  $nres = $ARGV[2];
}

# deg to rad conversion factor
$cf = 3.141592654 / 180;

# oxygen - hydrogen distance
$doh = 0.9572; 
# oxygen - extra point distance (consistent w/tip4pd.idp)
$doe = 0.1546;
# H-O-H angle  (consistent w/tip4pd.idp)
$ahoh = 52.26*2; 
$ahoh = $ahoh * $cf;
# H-O-E anle
$ahoe = 52.26; # (consistent w/tip4pd.idp)
$ahoe = $ahoe * $cf;
$cahoe = cos($ahoe);
$sahoe = sin($ahoe);

$ac = 0;	# atom count
while (<IN>) {
  if ($_ =~ (/^ATOM/)) {
    @line = split();
    # get x, y, z
    $crd[$ac][0] = $line[5];
    $crd[$ac][1] = $line[6];
    $crd[$ac][2] = $line[7];
    $ac++;
    if ($ac == 3) { # we have the whole water
      # coordinate of oxygen doesn't change
      $o[0] = $crd[0][0]; $o[1] = $crd[0][1]; $o[2] = $crd[0][2];

      # get both (original) o->h unit vectors 
      $oh1[0] = $crd[1][0] - $crd[0][0]; $oh2[0] = $crd[2][0] - $crd[0][0];
      $oh1[1] = $crd[1][1] - $crd[0][1]; $oh2[1] = $crd[2][1] - $crd[0][1];
      $oh1[2] = $crd[1][2] - $crd[0][2]; $oh2[2] = $crd[2][2] - $crd[0][2];
      $oh1l = sqrt($oh1[0]**2 + $oh1[1]**2 + $oh1[2]**2);
      $oh2l = sqrt($oh2[0]**2 + $oh2[1]**2 + $oh2[2]**2);
      $oh1[0] = $oh1[0] / $oh1l; $oh1[1] = $oh1[1] / $oh1l; $oh1[2] = $oh1[2] / $oh1l;
      $oh2[0] = $oh2[0] / $oh2l; $oh2[1] = $oh2[1] / $oh2l; $oh2[2] = $oh2[2] / $oh2l;
#print("oh1 coords: $oh1[0] $oh1[1] $oh1[2]\n");
#print("oh2 coords: $oh2[0] $oh2[1] $oh2[2]\n");

      # get bisector unit vector o->e
      $bis[0] = $oh1[0] + $oh2[0]; $bis[1] = $oh1[1] + $oh2[1]; $bis[2] = $oh1[2] + $oh2[2];
      $bisl = sqrt($bis[0]**2 + $bis[1]**2 + $bis[2]**2);
      $bis[0] = $bis[0] / $bisl; $bis[1] = $bis[1] / $bisl; $bis[2] = $bis[2] / $bisl;
      # set extra point coords
      $ep[0] = $o[0] + $bis[0] * $doe; $ep[1] = $o[1] + $bis[1] * $doe; $ep[2] = $o[2] + $bis[2] * $doe;
#$bisl = sqrt($bis[0]**2 + $bis[1]**2 + $bis[2]**2);
#print("bis: $bis[0] $bis[1] $bis[2] $bisl\n");

      # set new o->h unit vectors - rotate o->e around normal vector by +- HOE
      # get water molecule plane unit normal vector
      $wn[0] = $oh1[1] * $oh2[2] - $oh1[2] * $oh2[1];
      $wn[1] = (-1) * $oh1[0] * $oh2[2] + $oh1[2] * $oh2[0];
      $wn[2] = $oh1[0] * $oh2[1] - $oh1[1] * $oh2[0];
      $wnl = sqrt($wn[0]**2 + $wn[1]**2 + $wn[2]**2);
      $wn[0] = $wn[0] / $wnl; $wn[1] = $wn[1] / $wnl; $wn[2] = $wn[2] / $wnl;
#$wnl = sqrt($wn[0]**2 + $wn[1]**2 + $wn[2]**2);
#print("wnl: $wn[0] $wn[1] $wn[2] $wnl\n");
      
      # get third lattice vector perpendicular to BIS and WN - WP
      $wp[0] = $bis[1] * $wn[2] - $bis[2] * $wn[1];
      $wp[1] = $bis[2] * $wn[0] - $bis[0] * $wn[2];
      $wp[2] = $bis[0] * $wn[1] - $bis[1] * $wn[0];
      $wpl = sqrt($wp[0]**2 + $wp[1]**2 + $wp[2]**2);
      $wp[0] = $wp[0] / $wpl; $wp[1] = $wp[1] / $wpl; $wp[2] = $wp[2] / $wpl;
#$wpl = sqrt($wp[0]**2 + $wp[1]**2 + $wp[2]**2);
#print("wp: $wp[0] $wp[1] $wp[2] $wpl\n");
      
      # perform transformation of the coordinate system
      # x, y, z -> wp, bis, wn
      # build transformation matrix
      $tm[0][0] = $wp[0]; $tm[0][1] = $wp[1]; $tm[0][2] = $wp[2];
      $tm[1][0] = $bis[0]; $tm[1][1] = $bis[1]; $tm[1][2] = $bis[2];
      $tm[2][0] = $wn[0]; $tm[2][1] = $wn[1]; $tm[2][2] = $wn[2];

#print("trans matrix\n");
#for ($c = 0; $c <= 2; $c++) {
#  print("$tm[$c][0], $tm[$c][1], $tm[$c][2]\n");
#}
      
      # build inverse transformation matrix
      for ($c = 0; $c <= 2; $c++) {
        for ($d = 0; $d <= 2; $d++) {
	  $itm[$c][$d] = $tm[$d][$c];
	}
      }

      # transform BIS to new coordinates
      for ($c = 0; $c <= 2 ; $c++) {
        $bist[$c] = $tm[$c][0] * $bis[0] + $tm[$c][1] * $bis[1] + $tm[$c][2] * $bis[2];
#print("$bist[$c]: $tm[$c][0] $tm[$c][1] $tm[$c][2]\n");
#$bist[$c] = $tm[$c][0] * $oh1[0] + $tm[$c][1] * $oh1[1] + $tm[$c][2] * $oh1[2];
      }
#print("bist: $bist[0], $bist[1], $bist[2]\n");
      $bistl = sqrt($bist[0]**2 + $bist[1]**2 + $bist[2]**2);
      $bist[0] = $bist[0] / $bistl; $bist[1] = $bist[1] / $bistl; $bist[2] = $bist[2] / $bistl;
#$bistl = sqrt($bist[0]**2 + $bist[1]**2 + $bist[2]**2);
#print("bist $bist[0], $bist[1], $bist[2] $bistl\n");
  
      # rotate bist by +- ahoe around wn to get new o->h vectors
      $oh1t[0] = $bist[0] * $cahoe - $bist[1] * $sahoe; 
      $oh2t[0] = $bist[0] * $cahoe + $bist[1] * $sahoe;
      $oh1t[1] = $bist[1] * $cahoe + $bist[0] * $sahoe;
      $oh2t[1] = $bist[1] * $cahoe - $bist[0] * $sahoe;
      $oh1t[2] = $bist[2];
      $oh2t[2] = $bist[2];

#print("oh1t coords: $oh1t[0] $oh1t[1] $oh1t[2]\n");
#print("oh2t coords: $oh2t[0] $oh2t[1] $oh2t[2]\n");

      # rotate bist by +- aeol around wp to get o->lp vectors
      $ol1t[0] = $bist[0];
      $ol2t[0] = $bist[0];
      $ol1t[1] = $bist[1] * $caeol - $bist[2] * $saeol;
      $ol2t[1] = $bist[1] * $caeol + $bist[2] * $saeol;
      $ol1t[2] = $bist[2] * $caeol + $bist[1] * $saeol;
      $ol2t[2] = $bist[2] * $caeol - $bist[1] * $saeol;

      # transform back to old coordinate system and get new unit o->h vectors
      for ($c = 0; $c <= 2 ; $c++) {
        $oh1[$c] = $itm[$c][0] * $oh1t[0] + $itm[$c][1] * $oh1t[1] + $itm[$c][2] * $oh1t[2];
        $oh2[$c] = $itm[$c][0] * $oh2t[0] + $itm[$c][1] * $oh2t[1] + $itm[$c][2] * $oh2t[2];
        $ol1[$c] = $itm[$c][0] * $ol1t[0] + $itm[$c][1] * $ol1t[1] + $itm[$c][2] * $ol1t[2];
        $ol2[$c] = $itm[$c][0] * $ol2t[0] + $itm[$c][1] * $ol2t[1] + $itm[$c][2] * $ol2t[2];
#$bis[$c] = $itm[$c][0] * $bist[0] + $itm[$c][1] * $bist[1] + $itm[$c][2] * $bist[2];
      }
#print("BIS coordinates after 2 transformations: $bis[0], $bis[1], $bis[2]\n");

      $oh1l = sqrt($oh1[0]**2 + $oh1[1]**2 + $oh1[2]**2);
      $oh2l = sqrt($oh2[0]**2 + $oh2[1]**2 + $oh2[2]**2);
      $oh1[0] = $oh1[0] / $oh1l; $oh1[1] = $oh1[1] / $oh1l; $oh1[2] = $oh1[2] / $oh1l;
      $oh2[0] = $oh2[0] / $oh2l; $oh2[1] = $oh2[1] / $oh2l; $oh2[2] = $oh2[2] / $oh2l;

#print("oh1 coords: $oh1[0] $oh1[1] $oh1[2]\n");
#print("oh2 coords: $oh2[0] $oh2[1] $oh2[2]\n");

      # set new hydrogen positions
      $h1[0] = $o[0] + $oh1[0] * $doh; $h1[1] = $o[1] + $oh1[1] * $doh; $h1[2] = $o[2] + $oh1[2] * $doh;
      $h2[0] = $o[0] + $oh2[0] * $doh; $h2[1] = $o[1] + $oh2[1] * $doh; $h2[2] = $o[2] + $oh2[2] * $doh;

      # set lone pairs...
      #$l1[0] = $o[0] + $ol1[0] * $dol; $l1[1] = $o[1] + $ol1[1] * $dol; $l1[2] = $o[2] + $ol1[2] * $dol;
      #$l2[0] = $o[0] + $ol2[0] * $dol; $l2[1] = $o[1] + $ol2[1] * $dol; $l2[2] = $o[2] + $ol2[2] * $dol;
      
      # write the water
      printf("%-6s%5d%4s  %3s  %4d    %8.3f%8.3f%8.3f\n", ATOM, $natom + 0, OW, SOL, $nres, $o[0], $o[1], $o[2]);
      printf("%-6s%5d%4s  %3s  %4d    %8.3f%8.3f%8.3f\n", ATOM, $natom + 1, HW2, SOL, $nres, $h1[0], $h1[1], $h1[2]);
      printf("%-6s%5d%4s  %3s  %4d    %8.3f%8.3f%8.3f\n", ATOM, $natom + 2, HW3, SOL, $nres, $h2[0], $h2[1], $h2[2]);
      printf("%-6s%5d%4s  %3s  %4d    %8.3f%8.3f%8.3f\n", ATOM, $natom + 3, MW4, SOL, $nres, $ep[0], $ep[1], $ep[2]);
      print("TER\n");
      $nres++;
      $natom += 4;
      $ac = 0;
    }
  } else {
    # skip the line
  }
}
print("END\n");
