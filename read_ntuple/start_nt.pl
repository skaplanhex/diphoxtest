#! /usr/local/bin/perl -w
# Ce programme genere:
# 1) les trois fichiers de parametres
# 2) les directoires pour stocker les fichiers .bs
# 3) un fichier executable qui permet de tout lancer a la fois
##############################################################
# cette routine prend deux listes @liste1 et @liste2 et en fait un
# tableau associatif avec @list1 comme cles et @list2 comme valeurs
# @list1 et @list2 doivent avoir le meme nombre d'elements
sub list_to_hash {
  my ($rl_list1,$rl_list2) = @_;
  my ($nb_element1,$nb_element2,$i,%hash_temp,$rh_hash_temp);
  $nb_element1 = @$rl_list1;
  $nb_element2 = @$rl_list2;
  if ($nb_element1 == $nb_element2) {
    for ($i=0;$i<$nb_element1;$i++) {
      $hash_temp{"$rl_list1->[$i]"} = $rl_list2->[$i];
    }
  }
  else {
    print "Erreur dans list_to_hash\n";
    print "les deux listes n'ont pas le memes nombres d'elements\n";
    print "nombre d'elements de la liste 1: $nb_element1\n";
    print "nombre d'elements de la liste 2: $nb_element2\n";
 }
  $rh_hash_temp = {%hash_temp};
  return($rh_hash_temp);
}
######################## partie pour ntuple ############################ 
##############################################################
use Getopt::Long;
GetOptions("histo=s" => \$fichier_histo);
# valeur par defaut du fichier de parametres
unless ($fichier_histo) {
  $fichier_histo = "param_ntuple.indat";
}
# initialisation
$nom_fichier_fortran_init = "init_paw.f";
$nom_fichier_fortran_end = "end_paw.f";
$nom_fichier_histogrammation = "further_param.f";
$nom_fichier_remplissage = "remplissage.f";
# $nb_event = 10;
$ivar = 0;
$ihisto_equi = 0;
$ihisto_nonequi = 0;
$iscatter_equi = 0;
%hash_var_histo = ();
@list_var_histo_clef = qw/ name_ntuple name_histo identique cut_grand_pt cut_petit_pt y_grand_max y_grand_min y_petit_max y_petit_min mass_min mass_max degree_or_rad/;
@list_titre_equi = ();
@list_titre_nonequi = ();
@list_variable_equi = ();
@list_variable_nonequi = ();
@list_order_equi = ();
@list_order_nonequi = ();
@list_nb_bin_equi = ();
@list_nb_bin_nonequi = ();
@list_ref_bin_value = ();
@list_xmin_equi = ();
@list_xmax_equi = ();
@list_ref_histo_equi = ();
%hash_variable = ();
%hash_variable_gamma_gamma = ( "mass_gamma_gamma" => ["xmgg"],
		   "pt_gamma" => ["pt3","pt4"],
		   "pt_pair" => ["qt"],
		   "phi_gamma_gamma" => ["fi34"],
		   "y_gamma" => ["y3","y4"],	
		   "y_boost" => ["yboost"], 	
		   "y_star" => ["ystar"],  	
		   "y_gamma_gamma" => ["ygg"],	
		   "cos_theta_star" => ["cos_thetas"],  
		   "p_out" => ["p_out3","p_out4"],		
		   "pt_balance" => ["z3_trig","z4_trig"],	
		   "x_fragmentation_1" => ["x3"],
		   "x_fragmentation_2" => ["x4"]);
%hash_variable_pion_gamma = ( "mass_pion_gamma" => ["xmgg"],
		   "pt_pion" => ["pt3"],
		   "pt_gamma" => ["pt4"],
		   "pt_pair" => ["qt"],
		   "phi_pion_gamma" => ["fi34"],
		   "y_pion" => ["y3"],	
		   "y_gamma" => ["y4"],	
		   "y_boost" => ["yboost"], 	
		   "y_star" => ["ystar"],  	
		   "y_pion_gamma" => ["ygg"],	
		   "cos_theta_star" => ["cos_thetas"],  
		   "p_out_pion" => ["p_out3"],		
		   "p_out_gamma" => ["p_out4"],		
		   "pt_balance_pion" => ["z3_trig"],	
		   "pt_balance_gamma" => ["z4_trig"],	
		   "x_fragmentation_pion" => ["x3"],
		   "x_fragmentation_gamma" => ["x4"]);
$name_histo = "histo.res";
# on lit le fichier de donnees et on remplit le hashage %hash_var_histo
open (EREAD,$fichier_histo) || die "cannot open $fichier_histo";
open (EWRITE,">$name_histo") || die "cannot open $name_histo";
$deja_fait = 1;
while(<EREAD>) {
  chomp;
# il faut absolument enlever les blancs en debut de ligne sinon l'ordre 
# unless ne joue pas son role
  s/^\s+//;
  unless (/^#/ || /^histo_equi/ || /^histo_nonequi/) {  
# on enleve les blancs qui pourraient se trouver dans ou autour des 
# variables que l'on lit
    s/\s+//g; 
    $hash_var_histo{$list_var_histo_clef[$ivar]} = "$_";
    print EWRITE "$list_var_histo_clef[$ivar] $_\n";
    $ivar++;
  }
  if ($hash_var_histo{"identique"} && $deja_fait) {
    $deja_fait = 0;
# on test si les particules produites sont identiques ou non
    if ($hash_var_histo{"identique"} eq "identical") {
      $rh_variable = \%hash_variable_gamma_gamma;
      %hash_variable = %hash_variable_gamma_gamma;
      $identique = ".TRUE.";
    }
    elsif ($hash_var_histo{"identique"} eq "nonidentical") {
      $rh_variable = \%hash_variable_pion_gamma;
      %hash_variable = %hash_variable_pion_gamma;
      $identique = ".FALSE.";
    }
  }
  if (/^histo_equi/) { 
# on enleve les blancs qui sont devant et derriere la chaine et l'on fait en 
# sorte qu'il n'y ait qu'un blanc entre les champs 
    s/\s$//;
    s/\s+/ /g;
    my ($temp_dummy,$temp_variable,$temp_order,$temp_cut,$temp_titre,$temp_nb_bin,
    $temp_xmin,$temp_xmax);
    print EWRITE "$_\n";
    ($temp_dummy,$temp_variable,$temp_order,$temp_cut,$temp_titre,
    $temp_nb_bin,$temp_xmin,$temp_xmax) = split(/ /,$_);
    unless (exists($hash_variable{$temp_variable})) {
      die "The variable $temp_variable does not exist\nPlease check the equidistant histogram variables\n";
    }
    push (@list_variable_equi,$temp_variable);
    push (@list_order_equi,$temp_order);
    push (@list_cut_equi,$temp_cut);
    push (@list_titre_equi,$temp_titre);
    push (@list_nb_bin_equi,$temp_nb_bin);
    push (@list_xmin_equi,$temp_xmin);
    push (@list_xmax_equi,$temp_xmax);
    $ihisto_equi++;
  }
  elsif (/^histo_nonequi/) { 
# on enleve les blancs qui sont devant et derriere la chaine et l'on fait en 
# sorte qu'il n'y ait qu'un blanc entre les champs 
    s/\s$//;
    s/\s+/ /g;
    my ($temp_dummy,$temp_variable,$temp_order,$temp_cut,$temp_titre,
    $temp_nb_bin,@temp_bin_value,$rl_temp);
    print EWRITE "$_\n";
    ($temp_dummy,$temp_variable,$temp_order,$temp_cut,$temp_titre,
    $temp_nb_bin,@temp_bin_value) = split(/ /,$_);
    unless (exists($hash_variable{$temp_variable})) {
      die "The variable $temp_variable does not exist\nPlease check the non equidistant histogram variables\n";
    }
    push (@list_variable_nonequi,$temp_variable);
    push (@list_order_nonequi,$temp_order);
    push (@list_cut_nonequi,$temp_cut);
    push (@list_titre_nonequi,$temp_titre);
    push (@list_nb_bin_nonequi,$temp_nb_bin);
    $rl_temp = \@temp_bin_value;
    push (@list_ref_bin_value,$rl_temp);
    $ihisto_nonequi++;
  }
  elsif (/^scatter_equi/) { 
# on enleve les blancs qui sont devant et derriere la chaine et l'on fait en 
# sorte qu'il n'y ait qu'un blanc entre les champs 
    s/\s$//;
    s/\s+/ /g;
    my ($temp_dummy,$temp_variablex,$temp_variabley,$temp_order,$temp_cut,$temp_titre,
    $temp_nb_binx,$temp_xmin,$temp_xmax,$temp_nb_biny,$temp_ymin,
    $temp_ymax);
    print EWRITE "$_\n";
    ($temp_dummy,$temp_variablex,$temp_variabley,$temp_order,$temp_cut,$temp_titre,
    $temp_nb_binx,$temp_xmin,$temp_xmax,$temp_nb_biny,$temp_ymin,
    $temp_ymax) = split(/ /,$_);
    unless (exists($hash_variable{$temp_variablex}) ||
     exists($hash_variable{$temp_variabley})) {
      die "One of the variable $temp_variablex or $temp_variabley does not exist\nPlease check the scatter-plot variables\n";
    }
    push (@list_scatter_variablex_equi,$temp_variablex);
    push (@list_scatter_variabley_equi,$temp_variabley);
    push (@list_scatter_order_equi,$temp_order);
    push (@list_scatter_cut_equi,$temp_cut);
    push (@list_scatter_titre_equi,$temp_titre);
    push (@list_scatter_nb_binx_equi,$temp_nb_binx);
    push (@list_scatter_xmin_equi,$temp_xmin);
    push (@list_scatter_xmax_equi,$temp_xmax);
    push (@list_scatter_nb_biny_equi,$temp_nb_biny);
    push (@list_scatter_ymin_equi,$temp_ymin);
    push (@list_scatter_ymax_equi,$temp_ymax);
    $iscatter_equi++;
  }
}
close (EWRITE);
close (EREAD);
#
#@temp = split(/\./,$hash_var_histo{"name_ntuple"});
$temp[0] = $hash_var_histo{"name_ntuple"};
$temp[0] =~ s/\.rz//;
print "coucou1 ".$hash_var_histo{"name_ntuple"}." $temp[0]\n";
$name_res = $temp[0].".res";
if ($temp[0] =~ m/ggd/) {
  @name_dot = ("dir = .true.","oneb = .false.","twob = .false.");
}
elsif ($temp[0] =~ m/ggo/) {
  @name_dot = ("dir = .false.","oneb = .true.","twob = .false.");
}
elsif ($temp[0] =~ m/ggt/) {
  @name_dot = ("dir = .false.","oneb = .false.","twob = .true.");
}
# pour chaque histo avec bins equidistants, on donne:
# un entier
# un identificateur
# la variable qui remplit l'histogramme
# un flag pour savoir si on remplit l'histogramme avec les resultats Leading Log ou Next-to-Leading Log
# un tableau (en reference) pour les coupures specifiques a cet histogramme
# un titre
# le nombre de bin
# le maximum et minimum  de la variable
# ici on cree une reference a tous les objets du package Histo_equi
$starting = 20;
for ($ih=0;$ih < $ihisto_equi;$ih++) {
    my $r_histo = Histo_equi::new(eval($ih+1),eval($starting+$ih),
		    $list_variable_equi[$ih],$list_order_equi[$ih],
		    $list_cut_equi[$ih],
		    $list_titre_equi[$ih],$list_nb_bin_equi[$ih],
		    $list_xmin_equi[$ih],$list_xmax_equi[$ih]);
    push (@list_ref_histo_equi,$r_histo);
}
# pour chaque histo avec bins non equidistants, on donne:
# un entier
# un identificateur
# la variable qui remplit l'histogramme
# un flag pour savoir si on remplit l'histogramme avec les resultats Leading Log ou Next-to-Leading Log
# un tableau (en reference) pour les coupures specifiques a cet histogramme
# un titre
# le nombre de bin
# un tableau donnant les valeurs inferieurs de chaque bin + la valeur superieure pour le dernier bin
# ici on cree une reference a tous les objets du package Histo_nonequi
$starting = $starting+$ihisto_equi;
for ($ih=0;$ih < $ihisto_nonequi;$ih++) {
    my $r_histo = Histo_nonequi::new(eval($ih+1),eval($starting+$ih),
		    $list_variable_nonequi[$ih],$list_order_nonequi[$ih],
		    $list_cut_nonequi[$ih],
		    $list_titre_nonequi[$ih],$list_nb_bin_nonequi[$ih],
		    $list_ref_bin_value[$ih]);
    push (@list_ref_histo_nonequi,$r_histo);
}
# pour chaque scatterplot avec bins equidistants, on donne:
# un entier
# un identificateur
# les variables qui remplissent le scatterplot
# un flag pour savoir si on remplit le scatterplot avec les resultats Leading Log ou Next-to-Leading Log
# un tableau (en reference) pour les coupures specifiques a ce scatterplot
# un titre
# le nombre de bin pour la variable x
# le maximum et minimum  de la variable x
# le nombre de bin pour la variable y
# le maximum et minimum  de la variable y
# ici on cree une reference a tous les objets du package Scatter_equi
$starting = $starting+$ihisto_nonequi;
for ($ih=0;$ih < $iscatter_equi;$ih++) {
    my $r_scatter = Scatter_equi::new(eval($ih+1),eval($starting+$ih),
		    $list_scatter_variablex_equi[$ih],
		    $list_scatter_variabley_equi[$ih],
		    $list_scatter_order_equi[$ih],
		    $list_scatter_cut_equi[$ih],
		    $list_scatter_titre_equi[$ih],
		    $list_scatter_nb_binx_equi[$ih],
		    $list_scatter_xmin_equi[$ih],
		    $list_scatter_xmax_equi[$ih],
		    $list_scatter_nb_biny_equi[$ih],
		    $list_scatter_ymin_equi[$ih],
		    $list_scatter_ymax_equi[$ih]);
    push (@list_ref_scatter_equi,$r_scatter);
}
#
# creation du fichier d'initialisation de paw
open (INITCREATE,">$nom_fichier_fortran_init")|| 
die "cannot create  $nom_fichier_fortran_init";
$valeur_de_ipawc = 3000000;
$initcreate = *INITCREATE;
print INITCREATE "\tsubroutine histo_init(dir,oneb,twob)\n";
print INITCREATE "\tlogical dir,oneb,twob\n";
print INITCREATE "\tparameter(iwpawc = $valeur_de_ipawc)\n";
print INITCREATE "\tcommon/pawc/hmemor(iwpawc)\n";
unless ($ihisto_equi == 0) {
  print INITCREATE "\tparameter (inbhisto_equi=".$ihisto_equi.")\n";
}
unless ($ihisto_nonequi == 0) {
  print INITCREATE "\tparameter (inbhisto_nonequi=".$ihisto_nonequi.")\n";
}
unless ($iscatter_equi == 0) {
  print INITCREATE "\tparameter (inbscatter_equi=".$iscatter_equi.")\n";
}
unless ($ihisto_equi == 0) {
  print INITCREATE "\tcommon/init_close_equi/ibin(inbhisto_equi),wxmin(inbhisto_equi),\n";
  print INITCREATE "     #\twxmax(inbhisto_equi)\n";
}
unless ($ihisto_nonequi == 0) {
  Histo_nonequi::print_common($initcreate,\@list_ref_histo_nonequi);
}
unless ($iscatter_equi == 0) {
  print INITCREATE "\tcommon/init_close_scatter_equi/ibinx(inbscatter_equi),\n";
  print INITCREATE "     #\tibiny(inbscatter_equi),\n";
  print INITCREATE "     #\twsxmin(inbscatter_equi),wsxmax(inbscatter_equi),\n";
  print INITCREATE "     #\twsymin(inbscatter_equi),wsymax(inbscatter_equi)\n";
}
for (my $i=0;$i<=2;$i++) {
  print INITCREATE "\t$name_dot[$i]\n";
}
print INITCREATE "\tiwpawt = 150000*7*inbevent/1000000\n";
print INITCREATE "\tif(iwpawt.gt.iwpawc) then\n";
print INITCREATE "\t  write(6,*)' you must have iwpawc at least ',iwpawt\n"; 
print INITCREATE "\tstop\n";  
print INITCREATE "\tendif\n";
print INITCREATE "\tcall hlimit(iwpawc)\n";
print INITCREATE "\tcall hropen(1,'histogramme',\n";
print INITCREATE "     #\t\'$hash_var_histo{\"name_histo\"}\',\n";
print INITCREATE "     #\t'n',1024,istat)\n";
print INITCREATE "\tcall hropen(2,'ntuple',\n";
print INITCREATE "     #\t\'$hash_var_histo{\"name_ntuple\"}\',\n";
print INITCREATE "     #\t' ',8192,istat)\n";
print INITCREATE "\tcall hrin(0,9999,0) ! lecture\n";
foreach $r_histo (@list_ref_histo_equi) {
    $r_histo->print_value($initcreate);
}
foreach $r_histo (@list_ref_histo_equi) {
    $r_histo->print_hbook($initcreate);
}
foreach $r_histo (@list_ref_histo_nonequi) {
    $r_histo->print_value($initcreate);
}
foreach $r_histo (@list_ref_histo_nonequi) {
    $r_histo->print_hbook($initcreate);
}
foreach $r_scatter (@list_ref_scatter_equi) {
    $r_scatter->print_value($initcreate);
}
foreach $r_scatter (@list_ref_scatter_equi) {
    $r_scatter->print_hbook($initcreate);
}
print_return($initcreate);
close (INITCREATE);
#
# creation du fichier de fermeture paw
open (ENDCREATE,">$nom_fichier_fortran_end")|| 
die "cannot create  $nom_fichier_fortran_end";
$endcreate = *ENDCREATE;
print ENDCREATE "\tsubroutine histo_close\n";
unless ($ihisto_equi == 0) {
  print ENDCREATE "\tparameter (inbhisto_equi=".$ihisto_equi.")\n";
}
unless ($ihisto_nonequi == 0) {
  print ENDCREATE "\tparameter (inbhisto_nonequi=".$ihisto_nonequi.")\n";
}
unless ($ihisto_equi == 0) {
  print ENDCREATE "\tcommon/init_close_equi/ibin(inbhisto_equi),wxmin(inbhisto_equi),\n";
  print ENDCREATE "     #\twxmax(inbhisto_equi)\n";
  print ENDCREATE "\tdimension wbin_size(inbhisto_equi),wnorma(inbhisto_equi)\n";
}
unless ($ihisto_nonequi == 0) {
  Histo_nonequi::print_common($endcreate,\@list_ref_histo_nonequi);
  Histo_nonequi::print_dimension($endcreate,\@list_ref_histo_nonequi);
}
unless ($iscatter_equi == 0) {
  print ENDCREATE "\tcommon/init_close_scatter_equi/ibinx(inbscatter_equi),\n";
  print ENDCREATE "     #\tibiny(inbscatter_equi),\n";
  print ENDCREATE "     #\twsxmin(inbscatter_equi),wsxmax(inbscatter_equi),\n";
  print ENDCREATE "     #\twsymin(inbscatter_equi),wsymax(inbscatter_equi)\n";
  print ENDCREATE "\tdimension wsbin_size(inbscatter_equi),wsnorma(inbscatter_equi)\n";
}
print ENDCREATE "\topen(unit=30,file='$name_res',\n";
print ENDCREATE "     #\tstatus='unknown')\n";
print ENDCREATE "\tread(30,100) inbevent,wfirst_int\n";
print ENDCREATE "100\tformat(1x,i12,e12.5)\n";
unless ($ihisto_equi == 0) {
  print ENDCREATE "\tdo i=1,inbhisto_equi\n";
  print ENDCREATE "\t  wbin_size(i) = (wxmax(i)-wxmin(i))/float(ibin(i))\n";
  print ENDCREATE "\t  wnorma(i) = wfirst_int/(float(inbevent)*wbin_size(i))\n"; 
  print ENDCREATE "\tenddo\n";
}  
unless ($iscatter_equi == 0) {
  print ENDCREATE "\tdo i=1,inbscatter_equi\n";
  print ENDCREATE "\t  wsbin_size(i) = (wsxmax(i)-wsxmin(i))/float(ibinx(i))\n";
  print ENDCREATE "     #\t*(wsymax(i)-wsymin(i))/float(ibiny(i))\n";
  print ENDCREATE "\t  wsnorma(i) = wfirst_int/(float(inbevent)*wsbin_size(i))\n";   
  print ENDCREATE "\tenddo\n";  
}
foreach $r_histo (@list_ref_histo_equi) {
    $r_histo->print_opera($endcreate);
}
foreach $r_histo (@list_ref_histo_nonequi) {
    $r_histo->print_opera($endcreate);
}
foreach $r_scatter (@list_ref_scatter_equi) {
    $r_scatter->print_opera($endcreate);
}
foreach $r_histo (@list_ref_histo_equi) {
    $r_histo->print_out($endcreate);
}
foreach $r_histo (@list_ref_histo_nonequi) {
    $r_histo->print_out($endcreate);
}
foreach $r_scatter (@list_ref_scatter_equi) {
    $r_scatter->print_out($endcreate);
}
print ENDCREATE "\tcall hrend('ntuple')\n";  
print ENDCREATE "\tcall hrend('histogramme')\n";  
print_return($endcreate);
close (ENDCREATE);
#
# creation du fichier des parametres pour l'histogrammation
open (HISTO,">$nom_fichier_histogrammation")|| 
die "cannot create  $nom_fichier_histogrammation";
$histo = *HISTO;
# print HISTO "\tsubroutine further_param(pt1min,pt2min,ymax,ymin,etmin,\n";
# print HISTO "     #\tr,xmasmax,xmasmin)\n";
print HISTO "\tsubroutine further_param(pt1min,pt2min,y1max,y1min,\n";
print HISTO "     #\ty2max,y2min,xmasmax,xmasmin,ident)\n";
print HISTO "\timplicit real*4 (a-h,l-z)\n";
print HISTO "\tlogical ident\n";
print HISTO "\tident = $identique\n";
print HISTO "\tpt1min = $hash_var_histo{\"cut_grand_pt\"}\n";
print HISTO "\tpt2min = $hash_var_histo{\"cut_petit_pt\"}\n";
# print HISTO "\tymax = $hash_var_histo{\"ymax\"}\n";
# print HISTO "\tymin = $hash_var_histo{\"ymin\"}\n";
print HISTO "\ty1max = $hash_var_histo{\"y_grand_max\"}\n";
print HISTO "\ty1min = $hash_var_histo{\"y_grand_min\"}\n";
print HISTO "\ty2max = $hash_var_histo{\"y_petit_max\"}\n";
print HISTO "\ty2min = $hash_var_histo{\"y_petit_min\"}\n";
# print HISTO "\tetmin = $hash_var_histo{\"etmax\"}\n";
# print HISTO "\tr = $hash_var_histo{\"r_isol\"}\n";
print HISTO "\txmasmin = $hash_var_histo{\"mass_min\"}\n";
print HISTO "\txmasmax = $hash_var_histo{\"mass_max\"}\n";
print_return($histo);
close (HISTO);
#
# creation du fichier de remplissage des histogrammes
open (REMPLISSAGE,">$nom_fichier_remplissage")|| 
die "cannot create  $nom_fichier_remplissage";
$remplissage = *REMPLISSAGE;
print REMPLISSAGE "\tsubroutine remplissage(iprov,xmgg,pt3,pt4,qt,fi34,\n";
print REMPLISSAGE "     #\ty3,y4,yboost,ystar,ygg,p_out3,p_out4,z3_trig,z4_trig,\n";
print REMPLISSAGE "     #\tcos_thetas,wx3,wx4,weight)\n";
print REMPLISSAGE "\timplicit real*4 (a-h,l-z)\n";
$order = "lo";
print REMPLISSAGE "\t  if (iprov.eq.11) then\n";
foreach $r_histo (@list_ref_histo_equi) {
    $r_histo->print_remplissage($remplissage,$order,$rh_variable);
}
foreach $r_histo (@list_ref_histo_nonequi) {
    $r_histo->print_remplissage($remplissage,$order,$rh_variable);
}
foreach $r_scatter (@list_ref_scatter_equi) {
    $r_scatter->print_remplissage($remplissage,$order,$rh_variable);
}
print REMPLISSAGE "\t  endif\n";
$order = "nlo";
foreach $r_histo (@list_ref_histo_equi) {
    $r_histo->print_remplissage($remplissage,$order,$rh_variable);
}
foreach $r_histo (@list_ref_histo_nonequi) {
    $r_histo->print_remplissage($remplissage,$order,$rh_variable);
}
foreach $r_scatter (@list_ref_scatter_equi) {
    $r_scatter->print_remplissage($remplissage,$order,$rh_variable);
}
print_return($remplissage);
close (REMPLISSAGE);
#
sub print_return{
  my ($df) = @_;
  print $df "\treturn\n";
  print $df "\tend\n";
}
#
 system("gmake");
#
########################################################################
package Histo_equi;
sub new {
  my ($i,$id,$variable,$order,$cut,$title,$nb_bin,$xmin,$xmax) = @_;
  my $r_hash_histo = {
    "id" => $id,
    "title" => $title,
    "variable" => $variable,
    "order" => $order,
    "cut" => $cut,
    "nb_bin" => "ibin($i)",
    "value_nb_bin" => $nb_bin,
    "min_bin" => "wxmin($i)",
    "value_min_bin" => $xmin,
    "max_bin" => "wxmax($i)",
    "value_max_bin" => $xmax,
    "normalisation" => "wnorma($i)"
  };
  bless $r_hash_histo,'Histo_equi';
  return $r_hash_histo;
}
#
sub print_value {
  my ($r_histo,$df) = @_;
  print $df "\t".$r_histo->{"nb_bin"}." = ".$r_histo->{"value_nb_bin"}."\n";
  print $df "\t".$r_histo->{"min_bin"}." = ".$r_histo->{"value_min_bin"}."\n";
  print $df "\t".$r_histo->{"max_bin"}." = ".$r_histo->{"value_max_bin"}."\n";
  print $df "c\n";
}
#
sub print_hbook {
  my ($r_histo,$df) = @_;
  print $df "\tcall hbook1(".$r_histo->{"id"}.",'".
  $r_histo->{"title"}."',".
  $r_histo->{"nb_bin"}.",".
  $r_histo->{"min_bin"}.",".
  $r_histo->{"max_bin"}.",".
  "0.)\n";
  print $df "\tcall hbook1(9".$r_histo->{"id"}.",'dummy',".
  $r_histo->{"nb_bin"}.",".
  $r_histo->{"min_bin"}.",".
  $r_histo->{"max_bin"}.",".
  "0.)\n";
  print $df "\tcall hbarx(".$r_histo->{"id"}.")\n";
  print $df "\tcall hbarx(9".$r_histo->{"id"}.")\n";
  print $df "c\n";
}
#
sub print_opera {
  my ($r_histo,$df) = @_;
  print $df "\tcall hopera(".$r_histo->{"id"}.",'+',9".
  $r_histo->{"id"}.",".
  $r_histo->{"id"}.",".
  $r_histo->{"normalisation"}.",1.)\n";
}
#
sub print_remplissage {
  my ($r_histo,$df,$order,$rh_variable) = @_;
  my ($rl_temp,$first_parenthese,$last_parenthese,$machin,$truc,$variable);
  if ($r_histo->{"order"} eq $order) {
    eval('$rl_temp = '.$r_histo->{"cut"});
    $first_parenthese = "";
    $last_parenthese = "";
    if (scalar(@$rl_temp) > 3){
      $first_parenthese = "(";
      $last_parenthese = ")";
    }
    if ($rl_temp->[0]) {
      for ($i=1;$i<=scalar(@$rl_temp);$i=$i+3) {
        $machin = ($i==1)? "\t  if  $first_parenthese":"     #\t  .and.";
        print $df "$machin";
        for ($j=1;$j<=scalar(@{$rh_variable->{$rl_temp->[$i]}});$j++) {
	  $truc = ($j == 1)? "":"     #\t  .and.";
          print $df "$truc($rl_temp->[$i-1].le.".
	  @{$rh_variable->{"$rl_temp->[$i]"}}[$j-1]
	  .".and.".
	  @{$rh_variable->{"$rl_temp->[$i]"}}[$j-1]
	  .".le.$rl_temp->[$i+1])\n";
	}
      }
      print $df "     #\t  $last_parenthese then\n";
      foreach $variable (@{$rh_variable->{$r_histo->{"variable"}}}) {
        print $df "\t    call hfill (".$r_histo->{"id"}.",$variable,0.,weight)\n";
      }
      print $df "\t  endif\n";
    }
    else {
      foreach $variable (@{$rh_variable->{$r_histo->{"variable"}}}) {
      print $df "\t  call hfill (".$r_histo->{"id"}.",$variable,0.,weight)\n";
      }
    }
  }
}
#
sub print_out {
  my ($r_histo,$df) = @_;
  print $df "\tcall hrout(".$r_histo->{"id"}.",ICYCLE,' ')\n";
}
#
########################################################################
package Histo_nonequi;
sub new {
  my ($i,$id,$variable,$order,$cut,$title,$nb_bin,$rl_bin_value) = @_;
  my $r_hash_histo = {
    "id" => $id,
    "title" => $title,
    "variable" => $variable,
    "order" => $order,
    "cut" => $cut,
    "nb_bin" => "jbin($i)",
    "value_nb_bin" => $nb_bin,
    "bin" => "wbin".$i,
    "value_bin" => $rl_bin_value,
    "normalisation" => "wnorma".$i
  };
  bless $r_hash_histo,'Histo_nonequi';
  return $r_hash_histo;
}
#
sub print_common {
  my ($df,$r_list) = @_;
  my ($chaine_ini,$chaine_fin,$nb_term);
  $nb_term = @$r_list;
  $nb_term--;
  my ($i) = 0;
  foreach my $r_histo (@$r_list) {
    if ($nb_term ==  $i) {
      $chaine_fin = "";
    } 
    else {
      $chaine_fin = ",";
    } 
    if ($i ==  0) {
      $chaine_ini = "\tcommon/init_close_nonequi/jbin(inbhisto_nonequi),";
    } 
    else {
      $chaine_ini = "     #\t";
    } 
    print $df $chaine_ini.$r_histo->{"bin"}."(".eval($r_histo->{"value_nb_bin"}+1).")".$chaine_fin."\n";
    $i++;
  }
  print $df "c\n";
}
#
sub print_dimension {
  my ($df,$r_list) = @_;
  my ($chaine_ini,$chaine_fin,$nb_term);
  $nb_term = @$r_list;
  $nb_term--;
  my ($i) = 0;
  foreach my $r_histo (@$r_list) {
    if ($nb_term ==  $i) {
      $chaine_fin = "";
    } 
    else {
      $chaine_fin = ",";
    } 
    if ($i ==  0) {
      $chaine_ini = "\tdimension ";
    } 
    else {
      $chaine_ini = "     #\t";
    } 
    print $df $chaine_ini.$r_histo->{"normalisation"}."(".eval($r_histo->{"value_nb_bin"}+1).")".$chaine_fin."\n";
    $i++;
  }
  print $df "c\n";
}
#
sub print_value {
  my ($r_histo,$df) = @_;
  print $df "\t".$r_histo->{"nb_bin"}." = ".$r_histo->{"value_nb_bin"}."\n";
  my ($i) = 0;
  foreach my $value (@{$r_histo->{"value_bin"}}) {
    $i++;
    print $df "\t".$r_histo->{"bin"}."($i) = ".$value."\n";
  }
  print $df "c\n";
}
#
sub print_hbook {
  my ($r_histo,$df) = @_;
  print $df "\tcall hbookb(".$r_histo->{"id"}.",'".
  $r_histo->{"title"}."',".
  $r_histo->{"nb_bin"}.",".
  $r_histo->{"bin"}.",".
  "0.)\n";
  print $df "\tcall hbookb(9".$r_histo->{"id"}.",'dummy',".
  $r_histo->{"nb_bin"}.",".
  $r_histo->{"bin"}.",".
  "0.)\n";
  print $df "\tcall hbarx(".$r_histo->{"id"}.")\n";
  print $df "\tcall hbarx(9".$r_histo->{"id"}.")\n";
  print $df "c\n";
}
#
sub print_opera {
  my ($r_histo,$df) = @_;
  print $df "\tdo j=1,".$r_histo->{"nb_bin"}."\n";
  print $df "\t  ".$r_histo->{"normalisation"}."(j) = wfirst_int/(float(inbevent)*\n";
  print $df "     #\t ( ".$r_histo->{"bin"}."(j+1)-".$r_histo->{"bin"}."(j)))\n";
  print $df "\tenddo\n";
  print $df "\tcall hpak(9".$r_histo->{"id"}.",".$r_histo->{"normalisation"}.")\n";
  print $df "\tcall hopera(".$r_histo->{"id"}.",'*',9".
  $r_histo->{"id"}.",".
  $r_histo->{"id"}.",1.,1.)\n";
}
#
sub print_remplissage {
  my ($r_histo,$df,$order,$rh_variable) = @_;
  my ($rl_temp,$first_parenthese,$last_parenthese,$machin,$truc,$variable);
  if ($r_histo->{"order"} eq $order) {
    eval('$rl_temp = '.$r_histo->{"cut"});
    $first_parenthese = "";
    $last_parenthese = "";
    if (scalar(@$rl_temp) > 3){
      $first_parenthese = "(";
      $last_parenthese = ")";
    }
    if ($rl_temp->[0]) {
      for ($i=1;$i<=scalar(@$rl_temp);$i=$i+3) {
        $machin = ($i==1)? "\t  if  $first_parenthese":"     #\t  .and.";
        print $df "$machin";
        for ($j=1;$j<=scalar(@{$rh_variable->{$rl_temp->[$i]}});$j++) {
	  $truc = ($j == 1)? "":"     #\t  .and.";
          print $df "$truc($rl_temp->[$i-1].le.".
	  @{$rh_variable->{"$rl_temp->[$i]"}}[$j-1]
	  .".and.".
	  @{$rh_variable->{"$rl_temp->[$i]"}}[$j-1]
	  .".le.$rl_temp->[$i+1])\n";
	}
      }
      print $df "     #\t  $last_parenthese then\n";
      foreach $variable (@{$rh_variable->{$r_histo->{"variable"}}}) {
        print $df "\t    call hfill (".$r_histo->{"id"}.",$variable,0.,weight)\n";
      }
      print $df "\t  endif\n";
    }
    else {
      foreach $variable (@{$rh_variable->{$r_histo->{"variable"}}}) {
      print $df "\t  call hfill (".$r_histo->{"id"}.",$variable,0.,weight)\n";
      }
    }
  }
}
#
sub print_out {
  my ($r_histo,$df) = @_;
  print $df "\tcall hrout(".$r_histo->{"id"}.",ICYCLE,' ')\n";
}
#
########################################################################
package Scatter_equi;
sub new {
  my ($i,$id,$variablex,$variabley,$order,$cut,$title,$nb_binx,$xmin,$xmax,
  $nb_biny,$ymin,$ymax) = @_;
  my $r_hash_scatter = {
    "id" => $id,
    "title" => $title,
    "variablex" => $variablex,
    "variabley" => $variabley,
    "order" => $order,
    "cut" => $cut,
    "nb_binx" => "ibinx($i)",
    "value_nb_binx" => $nb_binx,
    "min_binx" => "wsxmin($i)",
    "value_min_binx" => $xmin,
    "max_binx" => "wsxmax($i)",
    "value_max_binx" => $xmax,
    "nb_biny" => "ibiny($i)",
    "value_nb_biny" => $nb_biny,
    "min_biny" => "wsymin($i)",
    "value_min_biny" => $ymin,
    "max_biny" => "wsymax($i)",
    "value_max_biny" => $ymax,
    "normalisation" => "wsnorma($i)"
  };
  bless $r_hash_scatter,'Scatter_equi';
  return $r_hash_scatter;
}
#
sub print_value {
  my ($r_scatter,$df) = @_;
  print $df "\t".$r_scatter->{"nb_binx"}." = ".$r_scatter->{"value_nb_binx"}."\n";
  print $df "\t".$r_scatter->{"nb_biny"}." = ".$r_scatter->{"value_nb_biny"}."\n";
  print $df "\t".$r_scatter->{"min_binx"}." = ".$r_scatter->{"value_min_binx"}."\n";
  print $df "\t".$r_scatter->{"max_binx"}." = ".$r_scatter->{"value_max_binx"}."\n";
  print $df "\t".$r_scatter->{"min_biny"}." = ".$r_scatter->{"value_min_biny"}."\n";
  print $df "\t".$r_scatter->{"max_biny"}." = ".$r_scatter->{"value_max_biny"}."\n";
  print $df "c\n";
}
#
sub print_hbook {
  my ($r_scatter,$df) = @_;
  print $df "\tcall hbook2(".$r_scatter->{"id"}.",'".
  $r_scatter->{"title"}."',".
  $r_scatter->{"nb_binx"}.",\n";
  print $df "     #\t".$r_scatter->{"min_binx"}.",".
  $r_scatter->{"max_binx"}.",".
  $r_scatter->{"nb_biny"}.",".
  $r_scatter->{"min_biny"}.",".
  $r_scatter->{"max_biny"}.",".
  "0.)\n";
  print $df "\tcall hbook2(9".$r_scatter->{"id"}.",'dummy',".
  $r_scatter->{"nb_binx"}.",\n";
  print $df "     #\t".$r_scatter->{"min_binx"}.",".
  $r_scatter->{"max_binx"}.",".
  $r_scatter->{"nb_biny"}.",".
  $r_scatter->{"min_biny"}.",".
  $r_scatter->{"max_biny"}.",".
  "0.)\n";
  print $df "c\n";
}
#
sub print_opera {
  my ($r_scatter,$df) = @_;
  print $df "\tcall hopera(".$r_scatter->{"id"}.",'+',9".
  $r_scatter->{"id"}.",".
  $r_scatter->{"id"}.",".
  $r_scatter->{"normalisation"}.",1.)\n";
}
#
sub print_remplissage {
  my ($r_scatter,$df,$order,$rh_variable) = @_;
  my ($rl_temp,$first_parenthese,$last_parenthese,$machin,$truc,$variable);
  if ($r_scatter->{"order"} eq $order) {
    eval('$rl_temp = '.$r_scatter->{"cut"});
    $first_parenthese = "";
    $last_parenthese = "";
    if (scalar(@$rl_temp) > 3){
      $first_parenthese = "(";
      $last_parenthese = ")";
    }
    if ($rl_temp->[0]) {
      for ($i=1;$i<=scalar(@$rl_temp);$i=$i+3) {
        $machin = ($i==1)? "\t  if  $first_parenthese":"     #\t  .and.";
        print $df "$machin";
        for ($j=1;$j<=scalar(@{$rh_variable->{$rl_temp->[$i]}});$j++) {
	  $truc = ($j == 1)? "":"     #\t  .and.";
          print $df "$truc($rl_temp->[$i-1].le.".
	  @{$rh_variable->{"$rl_temp->[$i]"}}[$j-1]
	  .".and.".
	  @{$rh_variable->{"$rl_temp->[$i]"}}[$j-1]
	  .".le.$rl_temp->[$i+1])\n";
	}
      }
      print $df "     #\t  $last_parenthese then\n";
      print $df "\t  call hfill (".$r_scatter->{"id"}.",".
      @{$rh_variable->{$r_scatter->{"variablex"}}}[0].",".
      @{$rh_variable->{$r_scatter->{"variabley"}}}[0].",weight)\n";
      print $df "\t  endif\n";
    }
    else {
      print $df "\t  call hfill (".$r_scatter->{"id"}.",".
      @{$rh_variable->{$r_scatter->{"variablex"}}}[0].",".
      @{$rh_variable->{$r_scatter->{"variabley"}}}[0].",weight)\n";
    }
  }
}
#
sub print_out {
  my ($r_scatter,$df) = @_;
  print $df "\tcall hrout(".$r_scatter->{"id"}.",ICYCLE,' ')\n";
}
#
