#!/usr/bin/perl
# Modish, version 0.159
# Author: Gian Luca Brunetti, Politecnico di Milano - gianluca.brunetti@polimi.it.
# The subroutine createconstrdbfile has been modified by ESRU, University of Strathclyde, Glasgow (2018) to adapt it to the new ESP-r construction database format.
# All rights reserved, 2015-18.
# This is free software.  You can redistribute it and/or modify it under the terms of the
# GNU General Public License, version 3, as published by the Free Software Foundation.

use v5.14;
use Math::Trig;
use List::Util qw[ min max reduce shuffle any];
use List::MoreUtils qw(uniq);
use List::AllUtils qw(sum);
use Statistics::Basic qw(:all);
use Data::Dump qw(dump);
use Data::Dumper;
$Data::Dumper::Terse  = 1;
$Data::Dumper::Purity  = 1;
use Vector::Object3D::Polygon;
use Math::Polygon::Tree;
use Storable qw(lock_store lock_nstore lock_retrieve);
use feature 'say';
no strict;
no warnings;

# Modish is a program for modifying the shading factors in the ISH (shading and insolation) files of the ESP-r building performance simulation suite in order to make it take into account the solar reflections from obstructions.
# The program, more precisely,  brings into account the reflective effect of solar obstructions on solar gains in the ESP-r building models on the basis of irradiance ratios. Those ratios are obtained combining the direct radiation on a surface, calculated by the means of ESP-r and by the means of a raytracer (Radiance), and the total radiation on the same surface calculated by the means of the raytracer. Using proportions, the values of the total radiation to be input to ESP-r, and from it, the modifications to the shading coefficients needed to obtain that, are calculated.
#
# How the program works
# The effect of solar reflections is taken into account at each hour on the basis of the ratios between the irradiances measured at the models' surfaces in two transitional, fictious model derived from the primary model.
# The procedure in question has to variants, at the user's choice.

# In the first variant of the procedure, the irradiances calculated by the means of Radiance can be derived from the following two alternative sets of models (the choice between them has to be done in the configuration file "modish_defaults.pl"):
#
# 1)
# a) a model in which all the surfaces are reflective, excepted the obstructions, which are black;
# b) a model in which everything is reflective.
#
# 2)
# a) a model in which everything is black;
# b) a model in which all the surfaces are black, excepted the obstructions, which are reflective.
#
# The value given by 1 minus the irradiance ratios gives the diffuse shading factors that are put in the ISH file of the ESP-r model in place of the original values. Solution 1) is the most correct one, but there are cases in which the results of the two solutions converge.
# The value given by 1 minus the irradiance ratios gives the diffuse shading factors that are put in the ISH file of the ESP-r model in place of the original values. Solution 1) is the most correct one, but there are cases in which the results of the two solutions converge.

# In the other variant of the procedure, the irradiances calculated by the means of Radiance can be derived from the following two alternative sets of models:
# 3)
# a) a model in which all the surfaces are reflective, and the obstruction are absent;
# b) a model in which everything is reflective, and the obstructions are present.
# In the case of this variant, the "new" shading factors are not calculated as a correction of the original ESP-r's ISH shading factors, but with the only means of Radiance, again, on the basis of irradiance ratios.
#
# Variant 3) can be combined with variant 1) OR 2).

# The original ISH's ".shda" files are not substituted. Two new files are added in the "zone" folder of the ESP-r model: the ".mod.shda" file is usable by ESP-r. It features the newly calculated shading factors; the ".report.shda" file lists the original shading factors and, at the bottom, the irradiance ratios from which the new shading factors in the ".mod.shda" file have been derived. Note that when the radiation on a surface is increased, instead of decreased, as an effect of reflections on obstructions, the shading factor are negative.
#
# To launch Modish the following command has to be issued:
#
# perl ./modish PATH_TO_THE_ESP-r_CONFIGURATION_FILE.cfg zone_number surface_1_number surface_2_number surface_n_number
#
# For example:
#
# perl ././Modish.pm /home/x/model/cfg/model.cfg 1 7 9 (which means: calculate for zone 1, surfaces 7 and 9.)
#
# In calculating the irradiance ratios, the program defaults to the following settings: diffuse reflections: 1 ; direct reflections: 7; surface grid: 2 x 2; direction vectors for each surface: 1 ; distance from the surface for calculating the irradiances: 0.01 (metres); ratio of the of the original shading factor to the "new" shading factor under which the new shading factor is used to substitute the original one in the ".shda" file. If this value is 0, it is inactive, there is no threshold.
# These defaults are a compromise between quality and speed. They can be overridden by preparing a "modish_defaults.pl" file and placing it in the same directory from which modish is called. In that directory,
# the files "fix.sh" and "perlfix.pl" must also be present, and "fix.sh" should be chmodded 755.
#
# The content of a configuration file for the present defaults, for example, would be constituted by the following line (note that after it a line of comment follows):
#
# @defaults = ( [  2, 2 ], 1, 1, 7, 0.01, 1 );### i.e ( [ resolution_x, resolution_y ], $dirvectorsnum, $bounceambnum, $bouncemaxnum, $distgrid, $threshold )
#
# The value "$dirvectorsnum" controls the numbers of direction vectors that are used for computing the irradiances at each point of each grid over the surfaces. The values that currently can be chosen are 1, 5 and 17. When the points are more than 1, they are evenly distributed on a hemisphere following a precalculated geodetic pattern.
# Modish works with no limitation of number of surfaces, it can be used with complex geometries and takes into account parent and child surfaces.
#
# For the program to work correctly, the materials, construction and optics databases
# must be local to the ESP-r model.
# Note that in the present version of the program, the materials used for obstructions should be not shared
# by objects which are not reflective obstructions. If it is, that material should be duplicated with a different name in the ESP-r material database and be given a different name, so that one variant of the material can be assigned to the obstructions, and another to the zones surfaces.
#
# The value @calcprocedures in a configuration file controls the algorithmic variant that is used in the calculations.
# If @calcprocedures is set to "diluted", the algorithmic strategy 1 listed above will be used.
# @calcprocedures = ( "diluted" );
# Otherwise, the algorithmic strategy 2 will be used.
# If @calcprocedures is set to "radical", strategy 3 will be used. Together with the option "radical", the options "diluted" or "complete" can be specified.
# If calcprocedures is set to "keepdirshdf", the direct shading factors calculated by ESP-r's ISH will be kept.
# Experience suggests that the most appropriate settings are:
# @calcprocedures = ( "diluted", "radical" );
# and
# @calcprocedures = ( "diluted" );
#
# In a configuration file, values of specular ratio and rougnness can be specified
# ovverriding the values specified in the Radiance material database.
# To obtain that, values of the kind "construction:specularratio:roughnessvalue"
# should be specified in the configuration file, with, for example:
# < @specularratios = ( "polishedmetal:0.03:0.05" );
# like shown in the example configuration case listed at the bottom of this writing.
# This will enable the program to take into account the specular reflective components.
#
# Another manner for specifying a complete (ideal) specular reflection for a material
# in a configuration file, of the kind that can be obtained by settings
# @specularratios = ( "constructiontype:1:0" );
# is that of giving to the material the "mirror" property, in the following manner:
# @specularratios = ( "constructiontype:mirror" );
#
# Considerations on the efficiency of the program.
# The speed of the program largely depends on the number of times that the Radiance raytracer is called, which in turn depends on the resolution of the grid on the external surface which is being considered.
#
# One drawback of the procedure in question from the viewpoint of execution speed may seem to be that the number of calls to the raytracer is double the number of the grid points defined on the considered external surface(s) for taking into account the solar reflections from obstructions. But another implication of this strategy is that it makes possible to decouple the gridding resolution on the considered external surface(s) regarding the effect of direct and diffuse reflections from obstruction from those on: (a) the considered external surface(s), for what is aimed to calculating direct radiation; (b) the internal surfaces, as regards the insolation. This makes possible to adopt a low gridding resolution for the considered external surface(s) relative to the diffuse and specular solar reflections from obstructions while adopting a higher resolution for (a) and (b). Which entails that the calculations regarding the direct radiation, which are likely to be the most important quantitatively for determining the solar gains in the thermal zones, and which are much quicker to calculate than the ones performed by the raytracer (which are necessary for determining the amount of solar radiation reflected from obstructions) can be carried out with a higher resolution than those involved in the calculations of the raytracer, so as to avoid to slow down the calculations themselves by a considerable amount. The amount of computations spared in the described manner may be significant, because the gridding entailed in the calculations not requiring the raytracer is commonly in the order of tens (for example, 20 x 20), whilst a gridding suitable for the use of a raytracer in this kind of operation is commonly in the order of units (for example, 2 x 2).
#
# The alternative to this strategy would be that of calculating all the solar radiation explicitly by defining one only gridding density for each surface; one only for all the radiation components entailed: the direct one, the diffuse one, and the one (diffuse and specular) reflected from obstruction. But this would require a gridding resolution of compromise between the components. For this reason, the calculation efficiency of the Modish procedure is likely to be most of the times not lower, but rather higher, than the alternative one entirely relying on calls to a raytracer.
#
# Modish should work with Linux and the Mac.
#
#
######### EXAMPLE OF CONFIGURATION FILE TO BE NAMED "modish_defaults.pl", ###################
######### TO BE PLACED IN THE DIRECTORY FROM WHICH MODISH IS CALLED, ########################
######### TOGETHER WITH "fix.pl" AND "perlfix.pl". ##########################################
######### THE FIRST "#" DIGIT OF EACH LINE MUST BE UNCHECKED FOR IT TO WORK. ################



# Example of configuration file for modish, version 1.59.
#
#@defaults = ( [ 2, 2 ], 5, 1, 7, 0.01, 0.01 );
#
### The line above means: ( [ resolution_x, resolution_y ], $dirvectorsnum, $bounceambnum, $bouncemaxnum, $distgrid, $threshold )
### resolution_x, resolution_y are the gridding values
### The value "$dirvectorsnum" controls the numbers of direction vectors that are used for
### computing the irradiances at each point of each grid over the surfaces. The values that currently
### can be chosen are 1, 5 and 9. When the points are more than 1, they are evenly distributed
### on a hemisphere following a geodetic pattern.
### $bounceambnum are the number of the bounces of the diffuse light which are taken into account.
# $bouncemaxnum are the number of the bounces of direct light which are taken into account.
# $distgrid is the distance of the grid in meters outside the surfaces which are taken into account.
# $threshold is the threshold under which the changes of shading value are not taken into account.
# A value of "1" means that if the new shading value is increased instead of decreased in place of the old one,
# the change is not executed.
#
#
#@calcprocedures = ( "diluted", "gendaylit", "radical", "coexistent", "getweather", "getsimple", "espdirres" );
## Quick instructions:
## The best groups of settings of @calcprocedures for calculating the shading factors are likely to be
## overall the following ones in most cases:
## 1)
## @calcprocedures = ( "diluted", "gensky", "directlydirect" ); ## FOR USING A CIE SKY AND ANALOGUE MODELS BASED ON DIFFERENCES OF REFLECTIVITY. THIS IS THE DEFAULT SETTING, IF NOTHING IS SPECIFIED.
## 2)
## @calcprocedures = ( "diluted", "gendaylit", "getweather", "getsimple" ); ## FOR USING A PEREZ SKY BUILT WITH MONTHLY AVERAGE WEATHER DATA AND ANALOGUE MODELS BASED ON DIFFERENCES OF REFLECTIVITY.
## 3)
## @calcprocedures = ( "diluted", "radical", "gendaylit", "getweather", "getsimple", "keepdirshdf" ); ## FOR USING A PEREZ SKY BUILT WITH MONTHLY AVERAGE WEATHER DATA AND ANALOGUE MODELS BASED ON DIFFERENCES OF GEOMETRY.
## 4)
## @calcprocedures = ( "diluted", "gensky", "radical", "coexistent", "directlydirect" );
## OR:
## @calcprocedures = ( "diluted", "gendaylit", "radical", "coexistent", "getweather", "getsimple" ); ## FOR USING
## STRATEGY 1) OR 2) FOR DIRECT SHADING FACTORS AND STRATEGY 3) FOR DIFFUSE ONES.
##5)
## @calcprocedures = ( "diluted", "gensky", "radical", "directlydirect", "espdirres" ); # THE "espdirres"
## OPTION MAY BE USED WHEN LOTS OF DIRECT REFLECTIONS COME INTO PLAY
## Explanations follow.
## "diluted" means that the two models from which the shading ratios are derived
## are going to be the following:
## 1)
## a) a model in which all the surfaces are reflective,
## excepted the obstructions, which are black;
## b) a model in which everything is reflective.
## 2)
## if "complete" is specified, the two models from which the shading ratios
## are derived are going to be the following:
## a) a model in which everything is black, and
## b) a model in which all the surfaces are black, excepted the obstructions,
## which are reflective. The settings "diluted" and "complete" are alternatives.
## 3)
## If the strategy "radical" is specified, the following sets of models are going to be used:
## a) a model in which all the surfaces are reflective, and the obstruction are absent;
## b) a model in which everything is reflective, and the obstructions are present.
## The options "diluted" or "complete" can be specified together with
## the option "radical".
## The option "coexistent" specified together with the option "radical" allows
## to utilize strategy 3 ("radical") for the diffuse calculations and strategy 1) or 2)
## for the direct calculations.
## If calcprocedures is set to "keepdirshdf", the direct shading factors calculated by ESP-r's ISH will be kept.
## If the setting "altcalcdiff" is specified, the diffuse irradiances are calculated
## from skies with no sun (be careful). If "directlydirect" is specified, the direct irradiances
## are calculates with skies with suns and imposing 0 diffuse bounces. This is the
## fastest method for direct calculations.
## If either "altcalcdiff" or "directlydirect" are specified, the total irradiances are
## calculated with skies with sun and complete diffuse bounces. But if both
## "altcalcdiff" or "directlydirect" are specified, the total irradiances are calculated
## as sum of direct and diffuse irradiances.
## If "gensky" is specified, the irradiances are calculated using the gensky program
## of Radiance, entailing the use of the CIE standard skies, for both the diffuse and direct
## calculations, and the result is sensible to the setting of sky condition for each month (below:
## clear, cloudy, or overcast).
## If "gendaylit" is specified, the irradiances are calculated using the gendaylit program
## of Radiance, entailing the use of the Perez sky model for the diffuse
## calculations and the direct ones. If the "getweather" setting is not specified,
## the direct calculations are performed by the means of gensky. If "getweather" is specified,
## the both the direct and the diffuse calculations are used with gendaylit by the means
## of averages of the weather data about direct normal and horizontal diffuse irradiances.
## For the setting "gendaylit" to work, it has to be specified together with the "altcalcdiff" setting.
## The setting "sparecycles" is specified together with the "radical" setting, it makes
## possible to spare a certain amount of computations (calls to Radiance).
## It ensures that the direct calculations are calculated directly, without
## requesting total irradiance calculations. With gendaylit it is more correct,
## and in any case it is quicker.
## The setting "getweather" used with "gendaylit" ("it can't be used without it)
## makes possible that the average radiation values of the weather data are utilized
## when calling gendaylit.
## The option "getsimple" used with "getweather" (it can't be used without it)
## determines the fact that the proportion of direct to diffuse radiation
## is determined directly from the shading data and overriding the other methods
## (determined by "altcalcdiff" and "directlydirect", or their absence) for defining that ratio.
## Note that the materials used in the obstructions should be not shared
## by objects which are not obstructions. If necessary, to obtain that,
## some materials may have to be suitably duplicated and renamed.
## The option "espdirres" makes possible to adopt the esp-r setting for the resolution
## of the direct calculations with Radiance with the "directlydirect" option,
## or "getweather", "getsimple".
## By specifying in @calcprocedure items of the kind "light/infrared-ratio:materialname:ratio"
## (for example: "light/infrared-ratio:gypsum:1.2" ) it is possible to model
## obstruction material which are selective in reflection - i.e. having different
## reflectivities in the range of light and solar infrared.
#
#@specularratios = (  );
#
##@specularratios = ( "reflector:mirror" );
##@specularratios = ( "reflector:0.03:0.05" );
## Here values of the kind "construction:specularratio:roughnessvalue"
## may be specified. For example, "reflector:0.03:0.05".
## The textual element ("reflector") is the name
## of a construction. The first number is the specular ratio
## for that construction. The second number is the roughness value.
## Specifying those values here makes possible
## to override the values specified in a Radiance database.
## (for example, the "0"s that may be in the database
## by defaul as regards specular ratios and roughness values).
## As an alternative, a material can be declared to be of the "mirror" type.
## This is done by specifying a value "construction:mirror".
## For example: reflector:mirror (see Radiance documentation
## about the properties of the "mirror" material type).
#
#%skycondition = ( 1=> "clear", 2=> "clear", 3=> "clear", 4=> "clear", 5=> "clear", 6=> "clear", 7=> "clear", 8=> "clear", 9=> "clear", 10=> "clear", 11=> "clear", 12=> "clear" );
## PREVAILING CONDITION OF THE SKY FOR EACH MONTH, EXPRESSED WITH ITS NUMBER, IN THE CASE IN WHICH
## CIE SKIES (WITH GENSKY) ARE UTILIZED.
## THE OPTIONS ARE: "clear", "cloudy" and "overcast".
## IF NO VALUE IS SPECIFIED, THE DEFAULT IS "clear".


######### END OF EXAMPLE CONFIGURATION LINES ###############################################


######## BEGINNING OF THE CONTENT OF THE FILE "perlfix.pl", ################################
######### TO BE PLACED IN THE DIRECTORY FROM WHICH MODISH IS CALLED ########################

##!/usr/bin/perl
#
#if ( -e "./fixl.pl" )
#{
#
#	open ( FIXL, "./fixl.pl" ) or die( $! );
#	my @files = <FIXL>;
#	close FIXL;
#	$" = " ";
#	print "FILES: @files\n";
#	my $to = $files[0];
#	chomp( $to );
#	print "TO: $to\n";
#
#	my $from = $files[1];
#	chomp( $from );
#	print "FROM: $from\n";
#
#	`cp -f $from $to`;
###	print "WHAT I AM DOING IS: cp -f $from $to \n";
#}

######## END OF THE CONTENT OF THE FILE "perlfix.pl" #######################################


######## BEGINNING OF THE CONTENT OF THE FILE "fix.sh" #####################################
######### TO BE PLACED IN THE DIRECTORY FROM WHICH MODISH IS CALLED ########################

#perl ./perlfix.pl

######## END OF THE CONTENT OF THE FILE "fix.sh" ###########################################


############################################################################################
######### BEGINNING OF MODISH ##############################################################


my $max_processes = $main::max_processes;
if ( not ( defined( $max_processes ) ) ) { $max_processes = 1; }

if ( ( "$^O" eq "MSWin32" ) or ( "$^O" eq "MSWin64" ) )
{
  say "\nSorry, this function presently works only on Linux and OSX." and die;
}

my ( @zoneshds, @winsdata );
my ( %surfslist, %shdfileslist, %obsinfo );

my %days_inmonths = ( Jan => 16, Feb => 15, Mar => 16, Apr => 16, May => 16, Jun => 16, Jul => 16, Aug => 16, Sep => 16, Oct => 16, Nov => 16, Dec => 16 );

my %monthsnum = ( Jan => 1, Feb => 2, Mar => 3, Apr => 4, May => 5, Jun => 6, Jul => 7, Aug => 8, Sep => 9, Oct => 10, Nov => 11, Dec => 12 );

sub getmonthname
{
  my $monthnum = shift;

  my %monthsnames = ( 1 => "Jan" , 2 => "Feb", 3 => "Mar", 4 => "Apr", 5 => "May", 6 => "Jun", 7 => "Jul", 8 => "Aug", 9 => "Sep", 10 => "Oct", 11 => "Nov", 12 => "Dec" );
  my $monthname = $monthsnames{ "$monthnum" };
  return $monthname;
}

sub getmonthnum
{
  my $monthname = shift;
  my %monthsnums = ( Jan => 1, Feb => 2, Mar => 3, Apr => 4, May => 5, Jun => 6, Jul => 7, Aug => 8, Sep => 9, Oct => 10, Nov => 11, Dec => 12 );
  my $monthnum = $monthsnums{ "$monthname" };
  return $monthnum;
}

sub getconffilenames
{  # THIS GETS THE CONSTRUCTION AND MATERIALS FILES FROM THE CFG FILE. IT IS CALLED BY sub createfictitious
  my ( $conffile, $path, $askedzonenum ) = @_;
  open ( CONFFILE, "$conffile") or die;
  my @lines = <CONFFILE>;
  close CONFFILE;

  my ( $constrdbfile, $matdbfile );
  my @zonedata;

  my ( $zonepath, $netpath, $ctlpath, $aimpath, $radpath, $imgpath,
    $docpath, $dbspath, $hvacpath, $bsmpath, $matfile, $constfile, $long );
  my %paths;

  my $semaphore = "no";
  my $countline = 0;
  foreach my $line (@lines)
  {
    my ($geofile, $constrfile, $shdfile, $zonenum_cfg );
    my @row = split(/\s+|,/, $line);

    if ( ( "#" ~~ @row ) and ( "Site" ~~ @row ) and ( "exposure" ~~ @row ) and ( "ground" ~~ @row ) and ( "reflectivity" ~~ @row ) )
    {
      my @minis;
      foreach my $mini ( @row )
      {
        unless ( $mini eq "" )
        {
          push( @minis, $mini );
        }
      }
      $groundrefl = $minis[1]; #say "\$groundrefl $groundrefl";
    }

    if ( ( "#" ~~ @row ) and ( "Latitude" ~~ @row ) and ( "Longitude" ~~ @row ) )
    {
      my @minis;
      foreach my $mini ( @row )
      {
        unless ( $mini eq "" )
        {
          push( @minis, $mini );
        }
      }
      $lat = $minis[0]; #say "\$lat $lat";
      $longdiff = $minis[1]; #say "\$longdiff $longdiff";
    }

    if ( $row[0] eq "*zonpth" )
    {
      $zonepath = $row[1];
      $zonepath =~ s/^\.//;
      $zonepath =~ s/^\.//;
      $zonepath =~ s/\/$//;
      if ( not ( $zonepath =~ /^\// ) )
      {
        $zonepath = "/" . $zonepath;
      }
      $zonepath = $path . $zonepath;
    }
    elsif ( $row[0] eq "*netpth" )
    {
      $netpath = $row[1];
      $netpath =~ s/^\.//;
      $netpath =~ s/^\.//;
      $netpath =~ s/\/$//;
      if ( not ( $netpath =~ /^\// ) )
      {
        $netpath = "/" . $netpath;
      }
      $netpath = $path . $netpath;
    }
    elsif ( $row[0] eq "*ctlpth" )
    {
      $ctlpath = $row[1];
      $ctlpath =~ s/^\.//;
      $ctlpath =~ s/^\.//;
      $ctlpath =~ s/\/$//;
      if ( not ( $ctlpath =~ /^\// ) )
      {
        $ctlpath = "/" . $ctlpath;
      }
      $ctlpath = $path . $ctlpath;
    }
    elsif ( $row[0] eq "*aimpth" )
    {
      $aimpath = $row[1];
      $aimpath =~ s/^\.//;
      $aimpath =~ s/^\.//;
      $aimpath =~ s/\/$//;
      if ( not ( $aimpath =~ /^\// ) )
      {
        $aimpath = "/" . $aimpath;
      }
      $aimpath = $path . $aimpath;
    }
    elsif ( $row[0] eq "*radpth" )
    {
      $radpath = $row[1];
      $radpath =~ s/^\.//;
      $radpath =~ s/^\.//;
      $radpath =~ s/\/$//;
      if ( not ( $radpath =~ /^\// ) )
      {
        $radpath = "/" . $radpath;
      }
      $radpath = $path . $radpath;
    }
    elsif ( $row[0] eq "*imgpth" )
    {
      $imgpath = $row[1];
      $imgpath =~ s/^\.//;
      $imgpath =~ s/^\.//;
      $imgpath =~ s/\/$//;
      if ( not ( $imgpath =~ /^\// ) )
      {
        $imgpath = "/" . $imgpath;
      }
      $imgpath = $path . $imgpath;
    }
    elsif ( $row[0] eq "*docpth" )
    {
      $docpath = $row[1];
      $docpath =~ s/^\.//;
      $docpath =~ s/^\.//;
      $docpath =~ s/\/$//;
      if ( not ( $docpath =~ /^\// ) )
      {
        $docpath = "/" . $docpath;
      }
      $docpath = $path . $docpath;
    }
    elsif ( $row[0] eq "*dbspth" )
    {
      $dbspath = $row[1];
      $dbspath =~ s/^\.//;
      $dbspath =~ s/^\.//;
      $dbspath =~ s/\/$//;
      if ( not ( $dbspath =~ /^\// ) )
      {
        $dbspath = "/" . $dbspath;
      }
      $dbspath = $path . $dbspath;
    }
    elsif ( $row[0] eq "*hvacpth" )
    {
      $hvacpath = $row[1];
      $hvacpath =~ s/^\.//;
      $hvacpath =~ s/^\.//;
      $hvacpath =~ s/\/$//;
      if ( not ( $hvacpath =~ /^\// ) )
      {
        $hvacpath = "/" . $hvacpath;
      }
      $hvacpath = $path . $hvacpath;
    }
    elsif ( $row[0] eq "*bsmpth" )
    {
      $bsmpath = $row[1];
      $bsmpath =~ s/^\.//;
      $bsmpath =~ s/^\.//;
      $bsmpath =~ s/\/$//;
      if ( not ( $bsmpath =~ /^\// ) )
      {
        $bsmpath = "/" . $bsmpath;
      }
      $bsmpath = $path . $bsmpath;
    }
    elsif ( $row[0] eq "*mat" )
    {
      $matdbfile = $row[1];
      $matdbfile =~ s/^\.//;
      $matdbfile =~ s/^\.//;
      if ( not ( $matdbfile =~ /^\// ) )
      {
        $matdbfile = "/" . $matdbfile;
      }
      $matdbfile = $path . $matdbfile; #say "\$matdbfile $matdbfile";
    }
    elsif ( $row[0] eq "*mlc" )
    {
      $constrdbfile = $row[1];
      $constrdbfile =~ s/^\.//;
      $constrdbfile =~ s/^\.//;
      if ( not ( $constrdbfile =~ /^\// ) )
      {
        $constrdbfile = "/" . $constrdbfile;
      }
      $constrdbfile = $path . $constrdbfile; #say "\$constrdbfile $constrdbfile";
    }
    elsif ( $row[0] eq "*clm" )
    {
      $clmfile = $row[1];
      $clmfile =~ s/^\.//;
      $clmfile =~ s/^\.//;
      if ( not ( $clmfile =~ /^\// ) )
      {
        $clmfile = "/" . $clmfile;
      }
      $clmfile = $path . $clmfile; #say "\$clmfile $clmfile";
    }
    #say "\$constrdbfile: " .  $constrdbfile ; ###.
    #say "\$matdbfile: " .  $matdbfile ; ###.
    #say "\$path $path";

    if ($row[0] eq "*zon")
    {
      $countzone++;
      my $zonenum = $row[1];
      if ( $zonenum eq $askedzonenum )
      {
        $semaphore = "yes";
        push ( @zonedata, $zonenum );
      }
    }

    if ( $semaphore eq "yes" )
    {
      if ($row[0] eq "*geo")
      {
        $geofile = $row[1];
        $geofile =~ s/^\.//;
        $geofile =~ s/^\.//;
        if ( not ( $geofile =~ /^\// ) )
        {
          $geofile = "/" . $geofile;
        }
        $geofile = $path . $geofile; #say "\$geofile $geofile";
        push ( @zonedata, $geofile );
      }

      if ( $row[0] eq "*con" )
      {
        $constrfile = $row[1];
        $constrfile =~ s/^\.//;
        $constrfile =~ s/^\.//;
        if ( not ( $constrfile =~ /^\// ) )
        {
          $constrfile = "/" . $constrfile;
        }
        $constrfile = $path . $constrfile; #say "\$constrfile $constrfile";
        push ( @zonedata, $constrfile );
      }
      #say "\$constrfile: " .  $constrfile ;###.


      if ( $row[0] eq "*isi")
      {
        $shdfile = $row[1];
        $shdfile =~ s/^\.//;
        $shdfile =~ s/^\.//;
        if ( not ( $shdfile =~ /^\// ) )
        {
          $shdfile = "/" . $shdfile;
        }
        $shdfile = $path . $shdfile; #say "\$shdfile $shdfile";
        push ( @zonedata, $shdfile );
        $semaphore = "no";
      }
    }
  }

  my $clmfilea = $clmfile . ".a"; #say "\$clmfilea: $clmfilea";

  if ( "getweather" ~~ @calcprocedures )
  {
    if ( not ( -e $clmfilea ) )
    {
`prj -file $conffile -mode script<<YYY
b
b
a
f

-
-
-
YYY
`;
    }
  }



#  if ( $longitude eq "" )
#  {
#    open( CLMFILE, "$clmfilea");
#    my @clmlines = <CLMFILE>;
#    close CLMFILE;
#    my @elts = split( ",", $clmlines[0] );
#    my $long = $elts[6];
#    if ( $long eq "" )
#    {
#      $long = $longitude;
#    }
#  }

  my $standardmeridian = 0; ################################ WORKING, BUT UNELEGANT FOR REPORTS
  #if ( $long eq "" )
  #{
    $long = $standardmeridian + $longdiff; ################################ WORKING, BUT UNELEGANT FOR REPORTS
  #}

  my $clmavgs = $clmfilea . "_avgs";

  $paths{zonepath} = $zonepath;
  $paths{netpath} = $netpath;
  $paths{ctlpath} = $ctlpath;
  $paths{aimpath} = $aimpath;
  $paths{radpath} = $radpath;
  $paths{imgpath} = $imgpath;
  $paths{docpath} = $docpath;
  $paths{dbspath} = $dbspath;
  $paths{hvacpath} = $hvacpath;
  $paths{bsmpath} = $bsmpath;
  $paths{matdbfile} = $matdbfile;
  $paths{constrdbfile} = $constrdbfile;
  $paths{conffile} = $conffile;
  $paths{lat} = $lat;
  $paths{longdiff} = $longdiff;
  $paths{long} = $long; #say "LONGITUDE: $long" ;
  $paths{groundrefl} = $groundrefl;
  $paths{standardmeridian} = $standardmeridian;
  $paths{clmfile} = $clmfile;
  $paths{clmfilea} = $clmfilea;
  $paths{clmavgs} = $clmavgs;

  my $cfgpath = $conffile;
  $cfgpath =~ s/\.cfg$// ;

  while ( not ( $cfgpath =~ /\/$/ ) )
  {
    $cfgpath =~ s/(\w+)$// ;
  }
  $cfgpath =~ s/\/$// ;

  $paths{cfgpath} = $cfgpath;

  return ( $constrdbfile, $matdbfile, \@zonedata, \@lines, \%paths );
}

sub createfictitiousfiles
{
  # THIS CREATES THE FILES FOR THE MODELS FEATURING FICTITIOUS QUALITIES AIMED TO THE MAIN Modish PROGRAM,
  # MODIFIES THE MATERIALS DB AS REQUESTED
  # _AND_ PREPARES THE CONFIGURATION FILES FOR THE FICTITIOUS MODELS
  my ($conffile, $path, $zonenum, $calcprocedures_ref ) = @_;
  my $conffile_f1 = $conffile;
  my ($flaggeo, $flagconstrdb, $flagmatdb, $flagconstr);
  $conffile_f1 =~ s/\.cfg/\_f1\.cfg/;
  my $conffile_f2 = $conffile;
  $conffile_f2 =~ s/\.cfg/\_f2\.cfg/;
  my $conffile_f3 = $conffile;
  $conffile_f3 =~ s/\.cfg/\_f3\.cfg/;
  my @calcprocedures = @{ $calcprocedures_ref };

  print REPORT "cp -R -f $conffile $conffile_f1\n";
  print REPORT "cp -R -f $conffile $conffile_f2\n";
  if ( "coexistent" ~~ @calcprocedures)
  {
    print REPORT "cp -R -f $conffile $conffile_f3\n";
  }
  #if ( not ( -e $conffile_f1 )
  #{
	`cp -R -f $conffile $conffile_f1\n`;
  #}
  `cp -R -f $conffile $conffile_f2\n`;
  if ( "coexistent" ~~ @calcprocedures)
  {
    `cp -R -f $conffile $conffile_f3\n`;
  }

  my ($constrdbfile, $matdbfile, $zonedataref, $conflinesref, $paths_ref ) = getconffilenames( $conffile, $path, $zonenum );
  my @zonedata = @$zonedataref;
  my $geofile = $zonedata[1];
  #say "MATDBFILE2: " .  $matdbfile ; ###.

  my @conflines = @$conflinesref;
  my %paths = %{ $paths_ref }; ################################say "PATHS RECEIVED: " . dump( \%paths ); ###.
  my (@originals, @fictitia1, @fictitia2, @fictitia3 );

  push ( @originals, $constrdbfile);

  my $constrdbfile_f = $constrdbfile;
  $constrdbfile_f = $constrdbfile . "_f" ;
  push ( @fictitia1, $constrdbfile_f);
  push ( @fictitia2, $constrdbfile_f);

  push ( @originals, $matdbfile);

  @fictitia3 = @fictitia2;

  my $matdbfile_f1 = $matdbfile;
  $matdbfile_f1 = $matdbfile . "_f1";
  push ( @fictitia1, $matdbfile_f1 );

  my $matdbfile_f2 = $matdbfile;
  $matdbfile_f2 = $matdbfile . "_f2";
  push ( @fictitia2, $matdbfile_f2 );

  push ( @fictitia3, $matdbfile_f1 );


  my ( @tempbox_original, @tempbox_fictitia1, @tempbox_fictitia2, @tempbox_fictitia3 );
  my $geofile = $zonedata[1];
  push ( @tempbox_originals, $geofile );

  my $geofile_f = $geofile;
  $geofile_f =~ s/\.geo/_f\.geo/;

  push ( @tempbox_fictitia1, $geofile_f); ######## IF ( "radical" ~~ @calcprocedures ) this is risky )
  push ( @tempbox_fictitia2, $geofile_f);
  push ( @tempbox_fictitia3, $geofile_f);

  my $constrfile = $zonedata[2];
  push ( @tempbox_originals, $constrfile );

  my $constrfile_f1 = $constrfile;
  $constrfile_f1 =~ s/\.con/_f1\.con/;
  push ( @tempbox_fictitia1, $constrfile_f1 );

  my $constrfile_f2 = $constrfile;
  $constrfile_f2 =~ s/\.con/_f2\.con/;
  push ( @tempbox_fictitia2, $constrfile_f2 );

  push ( @tempbox_fictitia3, $constrfile_f1 );


  print REPORT "cp -R -f $constrfile $constrfile_f1\n";
  print REPORT "cp -R -f $constrfile $constrfile_f2\n";

  `cp -R -f $constrfile $constrfile_f1\n`; $flagconstr = "y";
  `cp -R -f $constrfile $constrfile_f2\n`; $flagconstr = "y";


  my $shdfile = $zonedata[3];
  push ( @tempbox_originals, $shdfile);
  push ( @tempbox_fictitia1, $shdfile);
  push ( @tempbox_fictitia2, $shdfile);
  push ( @tempbox_fictitia3, $shdfile);

  my $zonenum_cfg = $zonedata[0];
  push ( @tempbox_originals, $zonenum_cfg);
  push ( @tempbox_fictitia1, $zonenum_cfg);
  push ( @tempbox_fictitia2, $zonenum_cfg);
  push ( @tempbox_fictitia3, $zonenum_cfg);

  push ( @originals, [ @tempbox_originals ] );
  push ( @fictitia1, [ @tempbox_fictitia1 ] );
  push ( @fictitia2, [ @tempbox_fictitia2 ] );
  push ( @fictitia3, [ @tempbox_fictitia3 ] );

  my ( @correctlines, $addline );
  open ( CONFFILE_F1, ">$conffile_f1");
  my @conflines2 = @conflines;

  foreach my $line (@conflines)
  {
    my $counter = 0;
    foreach my $elt ( @fictitia1 )
    {
      if ( not ( ref($elt) ) )
      {
        my $original = $originals[$counter];
        $elt =~ s/$path//;
        $original =~ s/$path//;
        if ( $elt )
        {
          $line =~ s/$original/$elt/;
        }
      }
      else
      {
        my @elts = @$elt;
        my @originalelts = @{$originals[$counter]};
        my $count = 0;
        foreach my $el ( @elts )
        {
          my $original = $originalelts[$count];
          $el =~ s/$path//;
          $original =~ s/$path//;
          if ( $el )
          {
            $line =~ s/$original/$el/;
          }
          $count++;
        }
      }

      if ( ( $counter == 0 ) and ( not ( $line =~ /^\*/ ) ) )
      {
        open( CORRECTCONF, $conffile ) or die;
        @correctlines = <CORRECTCONF>;
        close CORRECTCONF;
      }

      if ( @correctlines )
      {
        $addline = /^(.)$correctlines[ $counter ]/;
        if ( $addline )
        {
          $line = $addline . $line;
        }
      }
      $counter++;
    }
    print CONFFILE_F1 $line;
  }
  close CONFFILE_F1;


  open ( CONFFILE_F2, ">$conffile_f2");
  foreach my $line (@conflines2)
  {
    my $counter = 0;
    foreach my $elt ( @fictitia2 )
    {
      if ( not ( ref($elt) ) )
      {
        my $original = $originals[$counter];
        $elt =~ s/$path//;
        $original =~ s/$path//;
        if ( $elt )
        {
          $line =~ s/$original/$elt/;
        }
      }
      else
      {
        my @elts = @$elt;
        my @originalelts = @{$originals[$counter]};
        my $count = 0;
        foreach my $el ( @elts )
        {
          my $original = $originalelts[$count];
          $el =~ s/$path//;
          $original =~ s/$path//;
          if ( $el )
          {
            $line =~ s/$original/$el/;
          }
          $count++;
        }
      }
      if ( @correctlines )
      {
        $addline = /^(.)$correctlines[ $counter ]/;
        if ( $addline )
        {
          $line = $addline . $line;
        }
      }
      $counter++;
    }
    print CONFFILE_F2 $line;
  }
  close CONFFILE_F2;


  open ( CONFFILE_F1, "$conffile_f1");
  my @conflines1 = <CONFFILE_F1>;
  close CONFFILE_F1;

  my $take;
  foreach my $line ( @conflines1 )
  {
    if ( $line =~ /^\*mat/ )
    {
      chomp $line;
      my @elts = split( " +", $line );
      $take = $elts[1];
    }
  }

  open ( CONFFILE_F2, "$conffile_f2" );
  my @conflines2 = <CONFFILE_F2>;
  close CONFFILE_F2;

  open ( CONFFILE_F3, ">$conffile_f3" );

  foreach my $line (@conflines2)
  {
    if ( $line =~ /^\*mat/ )
    {
      chomp $line;
      my @elts = split( " +", $line );
      $line = "$elts[0]  $take\n";
    }
    print CONFFILE_F3 $line;
  }
  close CONFFILE_F3;


  if ( "bouncemethod" ~~ @calcprocedures )
  {
    $conffile_f1 = $conffile_f2;
    $matdbfile_f1 = $matdbfile_f2;
    @fictitia1 = @fictitia2;
  }

  return ( $conffile, $conffile_f1, $conffile_f2, $conffile_f3, $constrdbfile, $constrdbfile_f,
  $matdbfile, $matdbfile_f1, $matdbfile_f2, $flagconstrdb, $flagmatdb,
  $flaggeo, $flagconstr, [ @originals ], [ @fictitia1], [ @fictitia2 ], [ @fictitia3 ], \%paths );
}

sub definepath
{
  # THIS DEFINES THE PATH STARTING FROM THE PATH OF THE CONFIGURATION FILE. IT IS CALLED FROM sub modish
  my $launchfile = shift;
  my $path = $launchfile;

  if ( $path =~ /\/cfg.+/ )
  {
    $path =~ s/\/cfg.+// ;
  }
  else
  {
    $path =~ s/\.cfg// ;
    $path =~ s/\.cfg// ;

    while ( not ( $path =~ /\/$/ ) )
    {
      $path =~ s/(\w+)$// ;
    }
    $path =~ s/\/$// ;

  }
  return ( $path );
}

sub readgeofile
{  # THIS READS A GEO FILE TO GET THE DATA OF THE REQUESTED SURFACES
  my ( $geofile, $transpsurfs_ref, $zonenum, $calcprocedures_ref ) = @_;
  my @transpsurfs = @{ $transpsurfs_ref };
  my @calcprocedures = @{ $calcprocedures_ref };
  open ( GEOFILE, "$geofile") or die;
  my @lines = <GEOFILE>;
  close GEOFILE;
  my ( @geofilestruct, @transpelts, @obs );
  my %datalist;
  my $countsurf = 0;

  my @orresolutions;

  foreach my $surfnum ( @transpsurfs )
  {
    foreach my $line (@lines)
    {
      my @elts = split(/\s+|,/, $line);
      if ( $elts[0] eq "\*surf")
      {
        my $surfname = $elts[1];
        my $surfnum = $elts[12];
        my $parent = $elts[3];
        my $constr = $elts[6];
        # THIS POPULATES THE VARIABLE %surfslist (HASH - DICTIONARY - ASSOCIATIVE ARRAY) LINKING ZONES, SURFACES NAMES AND SURFACE NUMBER:

        $surfslist{ $zonenum }{ $surfnum }{surfname} = $surfname;
        $surfslist{ $zonenum }{ $surfname }{surfnum} = $surfnum;
        $datalist{ $zonenum }{ $surfnum }{ surfname } = $surfname;
        $datalist{ $zonenum }{ $surfnum }{ parent } = $parent;
        $datalist{ $zonenum }{ $surfnum }{ constr } = $constr;
        $datalist{ $zonenum }{ $surfnum }{ surfnum } = $surfnum;
        $datalist{ $zonenum }{ $surfnum }{ geofile } = $geofile;
        unless ( $parent eq "-" )
        {
          my $parentnum = $surfslist{ $zonenum }{ $parent }{surfnum};
          push ( @{ $datalist{ $zonenum }{ children }{ $parentnum } }, $surfnum );
          @{ $datalist{ $zonenum }{ children }{ $parentnum } } = uniq( @{ $datalist{ $zonenum }{ children }{ $parentnum } } );
        }
      }

      if ( $elts[0] eq "\*vertex")
      {
        my $x =  $elts[1];
        my $y =  $elts[2];
        my $z =  $elts[3];
        my $vertnum =  $elts[5];
        $datalist{ $zonenum }{ vertex }{ $vertnum } = [ $x, $y, $z ];
      }

      if ( $elts[0] eq "\*edges")
      {
        my $surfnum = $elts[ $#surfnum ];
        my $border = scalar( @elts - 3 );
        my @vertnums = @elts[ 1..$border ];
        $datalist{ $zonenum }{ $surfnum }{ edges }{ $surfnum } = [ @vertnums ];
      }

      if ( ($elts[0] eq "\*surf") and ( $surfnum == $elts[12] ) )
      {
        my $surfname = $elts[1];
        my $parent = $elts[3];
        my $constr = $elts[6];
        my $surfnum = $elts[12];
        push (@transpelts, [ $surfname, $parent, $constr, $surfnum, $geofile, $zonenum ] );
      }

      if ( ( $elts[0] eq "\*block_start" ) and ( "espdirres" ~~ @calcprocedures ) )
      {
        if ( $elts[1] eq "" )
        {
          push ( @orresolutions, $elts[2], $elts[3] );
        }
        else
        {
          push ( @orresolutions, $elts[1], $elts[2] );
        }
      }

      if ($elts[0] eq "\*obs")
      {
        my $obsconstr = $elts[10];
        my $obsname = $elts[9];
        my $obsnum = $elts[13];
        push (@obs, [ $obsname, $obsconstr, $obsnum ] );
      }

      if ( $countsurf == 0 )
      {
        push ( @geofilestruct, [ @elts ] );
      }
    }
    $countsurf++;
  }

  foreach ( @geofilestruct )
  {
    my $obsmaterial = $_->[9];
    my $obsnumber = $_->[12];
    if ( ( $_->[0] eq "#" ) and ( $_->[1] eq "*obs" ) )
    {
      $semaphore = 1;
    }
    if ( ( $semaphore == 1) and ( $_->[0] eq "*obs" ) )
    {
      push ( @obsconstr, $obsmaterial  );
    }
    $obsinfo{ $obsnumber } = $obsnumber;
    $obsinfo{ $obsmaterial } = $obsmaterial;
  }
  my @obsconstrset = uniq( @obsconstr );

  return ( \@transpelts, \@geofilestruct, \%surfslist, \@obs, \@obsconstrset, \%datalist, \@obsmaterials, \@orresolutions );
}

sub readverts
{
  # THIS READS THE VERTEX NUMBERS OF THE REQUESTED SURFACES IN A GEO FILE
  my @transpelts = @{$_[0]};
  my $geofile = $_[1];
  my @geodata = @{$_[2]};
  my %datalist = %{$_[3]};
  my @winselts;
  foreach my $transpelt (@transpelts)
  {
    my $surfname = $transpelt->[0];
    my $parent = $transpelt->[1];
    my $constr = $transpelt->[2];
    my $surfnum = $transpelt->[3];
    my $geofile = $transpelt->[4];
    my $zonenum = $transpelt->[5];
    my @winelts;
    foreach my $datum (@geodata)
    {
      my @data = @$datum;
      if ( ($data[0] eq "*edges") and ( $data[$#data] == $surfnum ) )
      {
        push ( @winelts, [ [ @data[ 2..( $#data - 2 ) ] ], $surfnum ] );
        my @surfverts = @data[ 2..( $#data - 2 ) ];
        $datalist{ $zonenum }{ $surfnum }{vertnums} = [ @surfverts ];
      }
    }
    push ( @winselts, [ @winelts ] );
  }
  return ( \@winselts, \%datalist );
}

sub readcoords
{
  # THIS READS THE COORDINATES OF THE REQUESTED VERTEX NUMBERS
  my ( $winseltsref, $geofile, $geodataref, $datalistref, $transpeltsref ) = @_;
  my @winselts = @$winseltsref;
  my @geodata = @$geodataref;
  my %datalist = %$datalistref;
  my @transpelts = @$transpeltsref;
  my @allcoords;
  my $count = 1;
  foreach my $winseltref (@winselts)
  {
    my @transpelt = @{ $transpelts[ $count -1 ] };
    my $zonenum = $transpelt[5];

    my @winselt = @$winseltref;
    my @vertnums = @{ $winselt[0][0] };
    my $surfnum = $winselt[0][1];
    my @coords;
    foreach my $num (@vertnums)
    {
      foreach my $datum (@geodata)
      {
        my @data = @$datum;
        if ( ($data[0] eq "*vertex") and ( $data[5] == $num ) )
        {
          push ( @coords, [ [ @data[ 1..3 ] ], $num ] );
          $datalist{ $zonenum }{ $num }{vertcoords} = [ @data[ 1..3 ] ];
        }
      }
    }
    push ( @allcoords, [ @coords ] );
    $count++;
  }
  return (\@allcoords, \%datalist );
}

sub getcorners
{
  # THIS PACKS THE X, Y, AND Z COORDINATES OF THE VERTICES OF THE REQUESTED SURFACES INTO SUBARRAYS
  my ( $winscoordsref, $winseltsref ) = @_;
  my @winscoords = @$winscoordsref;
  my @winselts = @$winseltsref;
  my @packsurfsdata;
  my $countsurf = 0;
  foreach $surfcoordsref ( @winscoords )
  {
    my @surfcoords = @$surfcoordsref;
    my ( @xdata, @ydata, @zdata );
    my @packsurfdata;
    my $surfnum = $winselts[$countsurf][0][1];
    foreach my $coordsetref (@surfcoords)
    {
      my @coordset = @$coordsetref;
      push (@xdata, $coordset[0][0]);
      push (@ydata, $coordset[0][1]);
      push (@zdata, $coordset[0][2]);
    }
    push (@packsurfdata, [ @xdata ], [ @ydata ], [ @zdata ], $surfnum  );
    push ( @packsurfsdata, [ @packsurfdata ] );
    $countsurf++;
  }
  return ( @packsurfsdata );
}

sub findextremes
{
  # THIS FINDS THE MAXIMA AND THE MINIMA FOR EACH COORDINATE FOR THE REQUESTED SURFACE
  my @xyzcoords = @_;
  my @surfsdata;
  foreach my $coordsdataref ( @xyzcoords )
  {
    my @coordsdata = @$coordsdataref;
    my $count = 0;
    my @surfdata;
    foreach $coordstyperef (@coordsdata)
    {
      if ($count < 3)
      {
        my @coordstype = @$coordstyperef;
        my $extreme1 = max(@coordstype);
        my $extreme2 = min(@coordstype);
        my $countpos = 0;
        my (@extreme1positions, @extreme2positions);
        foreach my $elt ( @coordstype )
        {
          if ( $elt ~~ $extreme1 )
          {
            push ( @extreme1positions, $countpos );
          }
          if ( $elt ~~ $extreme2 )
          {
            push ( @extreme2positions, $countpos );
          }
          $countpos++;
        }
        push ( @surfdata, [ [ $extreme1, [ @extreme1positions ] ], [ $extreme2, [ @extreme2positions ] ] ] );
        $count++;
      }
      else
      {
        if ( $surfdata[0][0][1] ~~ $surfdata[1][1][1] )
        {
          my $swap = $surfdata[1][1];
          $surfdata[1][1] = $surfdata[1][0];
          $surfdata[1][0] = $swap;
        }

        my $surfnum = $coordstyperef;
        push ( @surfdata, $surfnum );
      }
    }
    push (@surfsdata, [ @surfdata ] );
  }
  return ( @surfsdata );
}

sub makecoordsgrid
{
  # THIS FORMS A GRID OVER EACH REQUESTED SURFACE
  my ($extremesref, $resolutionsref, $dirsvectorsrefsref) = @_;
  my @extremesdata = @$extremesref;
  my @resolutions = @$resolutionsref;
  my @dirsvectorsrefs = @$dirsvectorsrefsref;
  my @wholegrid;
  my $countsurf = 0;
  foreach my $surfcase ( @extremesdata )
  {
    my $dirsvectorsref = $dirsvectorsrefs[$countsurf];
    my @surfdata = @$surfcase;
    my $surf = pop @surfdata;
    my @coordspoints;
    my $count = 0;
    foreach ( @surfdata )
    {
      my $extreme1 = $_->[0][0];
      my $extreme2 = $_->[1][0];
      my @extreme1positions = @{$_->[0][1]};
      my @extreme2positions = @{$_->[1][1]};
      my $resolution = $resolutions[$counter];
      my $diffextremes = ( $extreme1 - $extreme2 );
      my $variation = ( $diffextremes / ( $resolution + 1) );
      my @coordpoints;
      my $othercount = 1;
      while ( $othercount < ( $resolution +1 ) )
      {
        my $compoundvariation = ( $variation * $othercount );
        my $coordvalue = ( $extreme2 + $compoundvariation );
        push ( @coordpoints, $coordvalue );
        $othercount++;
      }
      push ( @coordspoints, [ @coordpoints ] );
      $count++;
    }
    push ( @coordspoints, $surf, $dirsvectorsref );
    push ( @wholegrid, [ @coordspoints ] );
    $countsurf++;
  }
  return(@wholegrid);
}

sub makegrid
{ # THIS CONVERTS THE GRID DATA IN VERTEX FORM
  my @gridcoords = @_;
  my @gridsvertices;
  foreach my $surfdataref ( @gridcoords )
  {
    my @xyzcoords;
    my @surfdata = @$surfdataref;
    my @xdata = @{$surfdata[0]};
    my @ydata = @{$surfdata[1]};
    my @zdata = @{$surfdata[2]};
    my $surf = $surfdata[3];
    my $dirvectorsref = $surfdata[4];
    my $counter = 0;
    my @gridvertices;
    my ( @groups, @xyzdata );
    foreach my $xdatum (@xdata)
    {
      my $ydatum = $ydata[$counter];
      push ( @xyzdata, [ $xdatum, $ydatum ] );
      $counter++;
    }
    foreach my $elt (@xyzdata)
    {
      foreach my $zdatum ( @zdata )
      {
        my @group = @$elt;
        push ( @group, $zdatum );
        push ( @groups, [ @group ] );
      }
    }
    push ( @gridvertices, [ @groups ], $surf, $dirvectorsref );
    push ( @gridsvertices, [ @gridvertices ] );
  }
  return ( @gridsvertices );
}

sub adjust_dirvector
{  # THIS SCALES THE DIRECTION VECTORS TO EASE THE MANIPULATION OF THE GRIDS IN DETACHING THEM FROM THE SURFACES.
  my ( $vectorref, $distgrid ) = @_;
  my @vector = @$vectorref;
  my $denominator = ( 1 / $distgrid );
  my @adjusted_vector;
  foreach my $elt ( @vector )
  {
    my $adjusted_component = ( $elt / $denominator );
    $adjusted_component = sprintf ( "%.3f", $adjusted_component );
    push ( @adjusted_vector, $adjusted_component );
  }
  return ( @adjusted_vector );
}

sub adjustgrid
{  # THIS ADJUSTS THE GRIDS OF POINTS OVER THE REQUESTED SURFACES BY DETACHING THEM OF ABOUT 1 CM TOWARDS THE OUTSIDE.
  my ( $griddataref, $distgrid )  = @_;
  my @griddata = @$griddataref;
  my @adjustedsurfs;
  foreach my $elt ( @griddata )
  {
    my @surfdata = @$elt;
    my @vertexdatarefs = @{$surfdata[0]};
    my $surfnum = $surfdata[1];
    my @dirvector = @{$surfdata[2]};
    my @adjusted_dirvector = adjust_dirvector( \@dirvector, $distgrid );
    my @adjustedsurfs;
    foreach my $vertexref ( @vertexdatarefs )
    {
      my @vertexcoords = @$vertexref;
      my @adjustedvertex;
      $countcomp = 0;
      foreach my $el ( @vertexcoords )
      {
        my $component = $adjusted_dirvector[$countcomp];
        my $newel = ( $el + $component );
        push ( @adjustedvertex, $newel );
        $countcomp++;
      }
      push ( @adjustedsurfs, [ @adjustedvertex ] );
    }
    push ( @adjusteddata, [ [ @adjustedsurfs ], $surfnum, [ @dirvector ] ] );
  }
  return ( @adjusteddata );
}

sub treatshdfile
{ # THIS PREPARES THE SHDA FILES IN MEMORY FOR USE.
  my @lines = @_;
  my @newlines;
  my $count = 0;
  foreach my $line ( @lines )
  {
    my $lineafter = $lines[ $count + 1 ];
    my $linebefore = $lines[ $count - 1 ];
    my $linecontrol = $line;
    if ( ( $lineafter =~ /# direct - surface/ ) or ( $lineafter =~ /# diffuse - surface/ ) )
    {
      $linecontrol = "";
    }
    elsif ( ( $line =~ /# direct - surface/ ) or ( $line =~ /# diffuse - surface/ ) )
    {
      chomp $linebefore;
      $line = "$linebefore" . " " . "$line" ;
    }

    unless ( $linecontrol eq "" )
    {
      push ( @newlines, $line );
    }
    $count++;
  }
  return ( @newlines );
}


sub readshdfile
{ # THIS READS THE RELEVANT CONTENT OF THE SHDA FILE.
  my ( $shdfile, $calcprocedures_ref, $conffile_f2, $message ) = @_;
  my @calcprocedures = @{ $calcprocedures_ref };

  if ( ( "keepdirshdf" ~~ @calcprocedures ) and ( $message eq "go" ) )
  {
`prj -file $conffile_f2 -mode script<<YYY
b
m
c
f
*
b
a
-
-
-
-
-
-
YYY
`;

say REPORT "prj -file $conffile_f2 -mode script<<YYY
b
m
c
f
*
b
a
-
-
-
-
-
-
YYY
";
  }

  my $shdafile = $shdfile . "a";
  if ( not ( -e $shdafile ) )#
  {
    say "\nExiting. A file \".shda\" must be present in the model folders for the operation to be performed.
    Now it isn't. To obtain that, a shading and insolation calculation must have been performed." and die;
  }
  my $tempfile = $shdafile;
  $tempfile =~ s/\.shda/\.temp\.shda/ ;

  open ( SHDAFILE, "$shdafile") or die;
  my @shdalines = <SHDAFILE>;
  close SHDAFILE,

  my (@filearray, @rawlines, @months);
  foreach my $line ( @shdalines )
  {
    push ( @rawlines, $line );
  }
  my @treatedlines = treatshdfile ( @rawlines );

  foreach my $line ( @treatedlines )
  { # THIS READS THE ".shda" FILES.
    my @elts = split(/\s+|,/, $line);
    if ( $line =~ /\* month:/ )
    {
      my $month = $elts[2];
      push ( @months, $month );
    }
    push ( @filearray, [ @elts ] );
  }

  open ( TEMP , ">$tempfile" ) or die;
  foreach my $line ( @treatedlines )
  {
    print TEMP $line;
  }
  close TEMP;
  return ( \@treatedlines, \@filearray, \@months );
}

sub tellsurfnames
{ # THIS RETURNS THE NAMES OF THE SURFACES.
  my ( $transpsurfsref, $geodataref ) = @_;
  my @transpsurfs = @$transpsurfsref;
  my @geodata = @$geodataref;
  my ( @containernums, @containernames, @nums, @names );
  my $count = 0;
  foreach my $surf ( @transpsurfs )
  {
    foreach my $rowref ( @geodata )
    {
      my @row = @$rowref;
      if ( ( $surf eq $row[12] ) and ( $row[0] eq "*surf" ) )
      {
        push ( @nums, $surf );
        push ( @names, $row[1] );
      }
      $count++;
    }
  }
  return ( \@nums, \@names );
}

sub getsurfshd
{ # THIS RETURN ALL THE DATA NEEDED FROM THE ".shda" FILE.
  my ( $shdfilearrayref, $monthsref, $surfnumsref, $surfnamesref ) = @_;
  my @shdfilearray = @$shdfilearrayref;
  my @months = @$monthsref;
  my @surfnums = @$surfnumsref;
  my @surfnames = @$surfnamesref;

  my @yearbag;
  foreach my $month ( @months )
  {
    my $semaphore = 0;
    my @monthbag;
    foreach my $rowref ( @shdfilearray )
    {
      my @row = @$rowref;
      if ( ( $row[0] eq "*") and ( $row[1] eq "month:" ) and ( $row[2] eq "$month" ) )
      {
        $semaphore = 1;
      }
      elsif ( ( $row[0] eq "*") and ( $row[1] eq "month:" ) and ( not ( $row[2] eq "$month" ) ) )
      {
        $semaphore = 0;
      }
      foreach my $surfname ( @surfnames )
      {
        if ( ( $row[25] eq "diffuse") and ( $row[27] eq "surface") and ( $row[28] eq "$surfname" ) and ( $semaphore == 1 ) )
        {
          push ( @monthbag, [ [ @row[0..23] ], $surfname ] );
        }
      }
    }
    push ( @yearbag, [ [ @monthbag ], $month ] );
  }
  return ( @yearbag );
}

sub checklight
{ # THIS LOOKS INTO THE "shda" DATA AND SEES WHEN DAYLIGHTING IS PRESENT.
  my ( $shdfilearrayref, $monthsref ) = @_;
  my @shdfilearray = @$shdfilearrayref;
  my @months = @$monthsref;

  my @yearbag;
  foreach my $month ( @months )
  {
    my @monthbag;
    my $countrow = 0;
    my $semaphore = 0;
    foreach my $rowref ( @shdfilearray )
    {
      my @row = @$rowref;
      if ( ( $row[0] eq "*") and ( $row[1] eq "month:" ) and ( $row[2] eq "$month" ) )
      {
        $semaphore = 1;
      }
      elsif ( ( $row[0] eq "*") and ( $row[1] eq "month:" ) and ( not ( $row[2] eq "$month" ) ) )
      {
        $semaphore = 0;
      }

      if ( ( $row[0] eq "surfaces") and ( $row[1] eq "insolated") and ( $semaphore == 1 ) )
      {
        my @bag;

        foreach my $el ( @{ $shdfilearray[ $countrow + 1 ] } )
        {
          if ( $el < 0 )
          {
            $el = 1;
            push ( @bag, $el );
          }
          else
          {
            $el = 0;
            push ( @bag, $el );
          }
        }
        push ( @monthbag, [ @bag ] );
      }
      $countrow++;
    }
    push ( @yearbag, [ [ @monthbag ], $month ] );
  }
  return ( @yearbag );
}

sub tellradfilenames
{ # THIS RETURNS THE NAMES OF THE RADIANCE FILES.
  my ( $path, $conffile_f1, $conffile_f2, $conffile_f3, $conffile_f4, $paths_ref ) = @_;
  my %paths =%{ $paths_ref };
  my @confs = ( $conffile_f1, $conffile_f2, $conffile_f3, $conffile_f4 );
  my @result;
  foreach my $conf ( @confs )
  {
    my $confstripped = $conf;

    if ( $confstripped =~ /$path\/cfg\// )
    {
      $confstripped =~ s/$path\/cfg\///;
    }
    else
    {
      $confstripped =~ s/$path\///;
    }

    $confstripped =~ s/.cfg//;
    my $radoctfile = "$confstripped" . "_Extern.oct";
    my $rcffile = "$confstripped" . ".rcf" ;
    push ( @result, [ $conf, $radoctfile, $rcffile ] );
  }
  return( @result );
}

sub tellradnames
{
  my ( $conffile, $path, $radpath ) = @_;
  my $confroot = $conffile;

  if ( $confroot =~ /$path\/cfg\// )
  {
    $confroot =~ s/$path\/cfg\/// ;
  }
  else
  {
    $confroot =~ s/$path\/// ;
  }

  $confroot =~ s/\.cfg$// ;
  my $fileroot = "$path/$confroot";
  my $rcffile = "$radpath/$confroot.rcf" ;
  my $radoctfile = "$radpath/$confroot" . "_Extern.oct";
  my $riffile = $radoctfile;
  $riffile =~ s/\.oct$/\.rif/ ;
  my $skyfile = $radoctfile;
  $skyfile =~ s/\.oct$/\.sky/ ;
  my $radmatfile = $radoctfile;
  $radmatfile =~ s/\.oct$/\.mat/ ;
  my $radmatcopy = $radmatfile . ".copy";
  my $diffskyfile = $skyfile;
  $diffskyfile =~ s/\.sky$/_diff\.sky/ ;
  return ( $fileroot, $rcffile, $radoctfile, $riffile, $skyfile, $radmatfile, $radmatcopy, $diffskyfile );
}

sub adjustlaunch
{
  my ( $skyfile, $diffskyfile, $path, $radpath ) = @_;

  $skyfile_short = $skyfile;
  $skyfile_short =~ s/$radpath\///;
  $diffskyfile_short = $diffskyfile;
  $diffskyfile_short =~ s/$radpath\///;

  open( SKYFILE, "$skyfile" ) or die "Can't open $skyfile $!";
  my @lines = <SKYFILE>;
  close SKYFILE;
  open( DIFFSKYFILE, ">$diffskyfile" ) or die "$!";
  foreach my $line ( @lines )
  {
    $line =~ s/^3 (.+)$/3 0 0 0/ ;
    $line =~ s/^4 (.+)$/4 0 0 0 0/ ;
    print DIFFSKYFILE $line;
  }
  close DIFFSKYFILE;

  my $oldskyfile = $skyfile . ".old";
  `mv -f $skyfile $oldskyfile`;
  print REPORT "mv -f $skyfile $oldskyfile\n";
  `mv -f $diffskyfile $skyfile`;
  print REPORT "mv -f $diffskyfile $skyfile\n";
}

sub setrad
{
  # THIS CREATES THE RADIANCE SCENES.
  my ( $conffile, $radoctfile, $rcffile, $path, $radpath, $monthnum, $day, $hour, $countfirst, $exportconstrref,
  $exportreflref, $skycondition_ref, $countrad, $specularratios_ref, $calcprocedures_ref, $debug, $paths_ref ) = @_;
  my %paths = %{ $paths_ref }; #say "\%paths_in_setrad: " . dump ( %paths );###.
  #say "\$radpath_in_setrad: " . dump ( $radpath ); ###.
  my %skycondition = %$skycondition_ref; #print REPORT "\%skycondition: " . dump ( %skycondition );
  my @calcprocedures = @$calcprocedures_ref;
  #say REPORT "IN SETRAD"; say "IN SETRAD";
  if ( $debug == 1 )
  {
    $debugstr = ">>out.txt";
  }
  else
  {
    $debugstr = "";
  }

  my $skycond = $skycondition{$monthnum}; #print REPORT "\$skycond: " . dump ( $skycond );

  my $radoctroot = $radoctfile;
  $radoctroot =~ s/$radoctfile/\.oct/ ;

  my $shortrcffile = $rcffile; #say  "\$rcffile: $rcffile";
  $shortrcffile =~ s/$radpath\/// ; #say REPORT "\$shortrcffile: $shortrcffile"; #say  "\$shortrcffile: $shortrcffile";

  my $skyfile = $rcffile;
  $skyfile =~ s/rif$/sky/ ;

  my $riffile = $rcffile;
  $riffile =~ s/\.rcf$/_Extern\.rif/ ;

  my $shortriffile = $riffile;
  $shortriffile =~ s/$radpath\/// ; #say REPORT "\$shortriffile: $shortriffile";

  my $add;
  if ( $skycond eq "cloudy" ) { $add = "\nf"; }
  if ( $skycond eq "overcast" ) { $add = "\nf\nf"; } #print REPORT "\$add: " . dump ( $add );

  my $moment;

  if ( ( ( $monthnum == 12 ) or ( $monthnum == 1 ) or ( $monthnum == 11 ) or ( $monthnum == 2 ) ) and ( $hour < 11 ) )
  { $moment = "a"; }
  elsif ( ( ( $monthnum == 12 ) or ( $monthnum == 1 ) or ( $monthnum == 11 ) or ( $monthnum == 2 ) ) and ( ( $hour == 11 ) or ( $hour == 12 ) or ( $hour == 13 ) ) )
  { $moment = "b"; }
  elsif  ( ( ( $monthnum == 12 ) or ( $monthnum == 1 ) or ( $monthnum == 11 ) or ( $monthnum == 2 ) ) and ( $hour > 13 ) )
  { $moment = "c"; }
  elsif ( ( ( $monthnum == 3 ) or ( $monthnum == 4 ) or ( $monthnum == 9 ) or ( $monthnum == 10 ) ) and ( $hour < 11 ) )
  { $moment = "d"; }
  elsif ( ( ( $monthnum == 3 ) or ( $monthnum == 4 ) or ( $monthnum == 9 ) or ( $monthnum == 10 ) ) and ( ( $hour == 11 ) or ( $hour == 12 ) or ( $hour == 13 ) ) )
  { $moment = "e"; }
  elsif  ( ( ( $monthnum == 3 ) or ( $monthnum == 4 ) or ( $monthnum == 9 ) or ( $monthnum == 10 ) ) and ( $hour > 13 ) )
  { $moment = "f"; }
  elsif  ( ( ( $monthnum == 5 ) or ( $monthnum == 6 ) or ( $monthnum == 7 ) or ( $monthnum == 8 ) ) and ( $hour < 11 ) )
  { $moment = "g"; }
  elsif  ( ( ( $monthnum == 5 ) or ( $monthnum == 6 ) or ( $monthnum == 7 ) or ( $monthnum == 8 ) ) and ( ( $hour == 11 ) or ( $hour == 12 ) or ( $hour == 13 ) ) )
  { $moment = "h"; }
  elsif  ( ( ( $monthnum == 5 ) or ( $monthnum == 6 ) or ( $monthnum == 7 ) or ( $monthnum == 8 ) ) and ( $hour > 13 ) )
  { $moment = "i"; }

  if ( not ( -e "$paths{cfgpath}/fix.sh" ) ) { `cp ./fix.sh $paths{cfgpath}/fix.sh`; }
  if ( not ( -e "$paths{cfgpath}/perlfix.pl" ) ) { `cp ./perlfix.pl $paths{cfgpath}/perlfix.pl`; }
  if ( not ( -e "$paths{radpath}/fix.sh" ) ) { `cp ./fix.sh $paths{radpath}/fix.sh`; }
  if ( not ( -e "$paths{radpath}/perlfix.pl" ) ) { `cp ./perlfix.pl $paths{radpath}/perlfix.pl`; }

if ( not ( -e "$paths{cfgpath}" ) )
{
  `mkdir $paths{cfgpath}`;
  print REPORT "mkdir $paths{cfgpath}";
}

  say REPORT "cd $paths{cfgpath}
e2r -file $conffile -mode text $debugstr <<YYY
c
$rcffile
a
a
$moment
1
n
d

d
$day $monthnum $hour$add
g
-
h
c
a
d
a
f
c
h
y
>
$riffile
-
-
YYY
.Done this.
";

`cd $paths{cfgpath}
e2r -file $conffile -mode text $debugstr <<YYY
c
$rcffile
a
a
$moment
1
n
d

d
$day $monthnum $hour
g
-
h
c
a
d
a
f
c
h
y
>
$riffile
-
-
YYY
`;

  if ( $countrad == 0 )
  {
    adjust_radmatfile( $exportconstrref, $exportreflref, $conffile, $path, \@specularratios,
    \%obslayers, \@selectives, \%paths );
  }
}

sub setroot
{ # THIS SETS THE MODELS' ROOT NAME.
  my ( $conffile, $path, $debug, $paths_ref ) = @_;
  my %paths = %{ $paths_ref };
  my $rootname = $conffile;
  if ( $rootname =~ /$path\/cfg\// )
  {
    $rootname =~ s/$path\/cfg\///;
  }
  else
  {
    $rootname =~ s/$path\///;
  }

  $rootname =~ s/\.cfg//;
  if ( $debug == 1 )
  {
    $debugstr = ">>out.txt";
  }
  else
  {
    $debugstr = "";
  }

  if ( not ( -e "$paths{cfgpath}/fix.sh" ) ) { `cp ./fix.sh $paths{cfgpath}/fix.sh`; }
  if ( not ( -e "$paths{cfgpath}/perlfix.pl" ) ) { `cp ./perlfix.pl $paths{cfgpath}/perlfix.pl`; }
  if ( not ( -e "$paths{radpath}/fix.sh" ) ) { `cp ./fix.sh $paths{radpath}/fix.sh`; }
  if ( not ( -e "$paths{radpath}/perlfix.pl" ) ) { `cp ./perlfix.pl $paths{radpath}/perlfix.pl`; }

print REPORT "cd $paths{radpath}
prj -file $conffile -mode text $debugstr <<YYY
b

s

$rootname

m
c
b
#
y
-
-
-
-

YYY
";
`cd $paths{radpath}
prj -file $conffile -mode text $debugstr <<YYY
b

s

$rootname

m
c
b
#
y
-
-
-
-

YYY
`;
}

sub populatelight
{ # THIS POPULATES THE DATA STRUCTURE DEDICATED TO SIGNAL THE DAYLIT HOURS.
  my @daylighthoursarr = @_;
  my %daylighthours;
  my $count = 0;
  foreach my $monthref ( @daylighthoursarr )
  {
    my @monthdata = @$monthref;
    my $month = $monthdata[1];
    $month =~ s/`//g;
    my @lithours = @{$monthdata[0][0]};
    $daylighthours{$month} = [ @lithours ] ;
    $count++;
  }
  return ( %daylighthours );
}


sub deg2rad
{
	my $degrees = shift;
	return ( ( $degrees / 180 ) * 3.14159265358979 );
}

sub rad2deg
{
	my $radians = shift;
	return ( ( $radians / 3.14159265358979 ) * 180 ) ;
}

sub rotate2d
{   # SELF-EXPLAINING.
    my ( $x, $y, $angle ) = @_;
    $angle = deg2rad( $angle );
    my $x_new = cos($angle)*$x - sin($angle)*$y;
    my $y_new = sin($angle)*$x + cos($angle)*$y;
  return ( $x_new, $y_new);
}

sub getdirvectors
{ # THIS GETS THE NEEDED DIRECTION VECTORS AT EACH GRID POINT DEPENDING FROM THE LAUNCH SETTINGS.
   ( $basevectorsref, $dirvectorref, $pointcoordsref ) = @_;
   my @basevectors = @$basevectorsref;;
   my @dirvector = @$dirvectorref;
   my @topcoords = @{$basevectors[0]};
   my @newdirvectors;
   my $xbase = $topcoords[0];
   my $ybase = $topcoords[1];
   my $zbase = $topcoords[2];
   my $xnew = $dirvector[0];
   my $ynew = $dirvector[1];
   my $znew = $dirvector[2];
   my $anglebasexz = acos($xbase);
   my $anglebaseyz = acos($zbase);
   my $anglenewxz = acos($xnew);
   my $anglenewyz = acos($znew);
   my $anglediffxz = ( $anglenewxz - $anglebasexz );
   my $anglediffyz = ( $anglenewyz - $anglebaseyz );
   foreach my $eltsref ( @basevectors )
   {
     my @elts = @$eltsref;
     my ( $x, $y, $z ) = @elts ;
     my ( $x_ok, $tempy ) = rotate2d( $x, $y, $anglediffxz );
     my ( $y_ok, $z_ok ) = rotate2d( $tempy, $z, $anglediffyz );
     $x_ok = sprintf ( "%.3f", $x_ok );
     $y_ok = sprintf ( "%.3f", $y_ok );
     $z_ok = sprintf ( "%.3f", $z_ok );
     push ( @newdirvectors, [ $x_ok, $y_ok, $z_ok ] );
   }
   return ( @newdirvectors );
}

sub pursue
{ # THIS CALCULATES THE IRRADIANCES BY THE MEANS OF RADIANCE.
  # RADIANCE EXAMPLE: echo 1 dat-0.01 2 0 -1 0 | rtrace  -I -ab 2 -lr 7 -h /home/luca/boxform/rad/boxform_f1_Extern.oct | rcalc -e '$1=179*(.265*$1+.670*$2+.065*$3)'
  $" = " ";
  my $dat = shift;
  my %d = %$dat;
  my $zonenum = $d{zonenum};
  my $geofile = $d{geofile};
  my %paths = %{ $d{paths} }; #say "PATHS IN PURSUE " . dump ( \%paths ); ###.
  my $constrfile = $d{constrfile};
  my $shdfile = $d{shdfile};
  my @gridpoints = @{ $d{gridpoints} };
  my @shdsurfdata = @{ $d{shdsurfdata} };
  my @daylighthoursarr = @{ $d{daylighthoursarr} };
  my %daylighthours =  %{ $d{daylighthours} };
  my @shdfilearray = @{ $d{shdfilearray} };
  my $exportconstrref = $d{exportconstrref};
  my $exportreflref = $d{exportreflref};
  my $conffile = $d{conffile};
  my $path = $d{path};
  my $radpath = $paths{radpath};
  my @basevectors = @{ $d{basevectors} };
  my @resolutions = @{ $d{resolutions} };
  my $dirvectorsnum = $d{dirvectorsnum};
  my @calcprocedures = @{ $d{calcprocedures} };
  my @specularratios = @{ $d{specularratios} };
  my $bounceambnum = $d{bounceambnum};
  my $bouncemaxnum = $d{bouncemaxnum};
  my ( $fict1ref, $fict2ref, $fict3ref ) = @{ $d{radfilesrefs} };
  my @surfsnums = @{ $d{transpsurfs} };
  my @selectives = @{ $d{selectives} };
  my $lat = $paths{lat};
  my $longdiff = $paths{longdiff};
  my $long = $paths{long};
  my $groundrefl = $paths{groundrefl};
  my $standardmeridian = $paths{standardmeridian};
  my $clmavgs = $paths{clmavgs};
  my @orgridpoints = @{ $d{orgridpoints} };

  my ( $conffile_f1, $radoctfile_f1, $rcffile_f1 ) = @$fict1ref;
  my ( $conffile_f2, $radoctfile_f2, $rcffile_f2 ) = @$fict2ref;
  my ( $conffile_f3, $radoctfile_f3, $rcffile_f3 ) = @$fict3ref;
  my ( $conffile_f4, $radoctfile_f4, $rcffile_f4 ) = @$fict4ref;
  my @conffiles = ( $conffile_f1, $conffile_f2, $conffile_f3, $conffile_f4 );
  my ( @radoctfiles, @rcffiles );

  if ( scalar( @selectives ) == 0 )
  {
    @radoctfiles = ( $radoctfile_f1, $radoctfile_f2  );
    @rcffiles = ( $rcffile_f1, $rcffile_f2 );
  }
  elsif ( scalar( @selectives ) > 0 )
  {
    @radoctfiles = ( $radoctfile_f1, $radoctfile_f2, $radoctfile_f4  );
    @rcffiles = ( $rcffile_f1, $rcffile_f2, $rcffile_f4 );
  }

  if ( "coexistent" ~~ @calcprocedures )
  {
    @radoctfiles = ( $radoctfile_f1, $radoctfile_f2, $radoctfile_f3  );
    @rcffiles = ( $rcffile_f1, $rcffile_f2, $rcffile_f3 );
  }


  my $resolnumber = ( $resolutions[0] * $resolutions[1] );

  my $fileroot = $radoctfile_f1;
  $fileroot =~ s/_f1_Extern.oct//;
  my ( %irrs, %irrs_dir, %irrs_amb);
  my $setoldfiles = "on";
  my ( %surftests, %surftestsdiff, %surftestsdir, %irrs );

  my ( $totaldirect, $totalrad, $directratio, $diffuseratio );
  my $countfirst = 0;
  my $countrad = 0;

  my ( $t_ref, $clmlines_ref, %t, @clmlines );
  if ( ( not( -e $clmavgs ) ) and ( "getweather" ~~ @calcprocedures ) )
  {
    ( $t_ref, $clmlines_ref ) = getsolar( \%paths );
    %t = %{ $t_ref };
    @clmlines = @{ $clmlines_ref };
  }

  if ( ( scalar( @clmlines ) == 0 ) and ( "getweather" ~~ @calcprocedures ) )
  {
    open( CLMAVGS, "$clmavgs" ) or die;
    @clmlines = <CLMAVGS>;
    close CLMAVGS;
  }

  my %avgs;
  foreach my $clmline ( @clmlines )
  {
    chomp $clmline;
    my @es = split( ",", $clmline );
    $avgs{dir}{$es[0]}{$es[2]} = $es[3];
    $avgs{diff}{$es[0]}{$es[2]} = $es[4];
    $avgs{alt}{$es[0]}{$es[2]} = $es[5];
    $avgs{azi}{$es[0]}{$es[2]} = $es[6];
  } #say REPORT "\%avgs: " . dump( \%avgs );

  foreach my $radoctfile ( @radoctfiles )
  { #say "\$radoctfile $radoctfile";
    my $conffile = $conffiles[$countrad];

    my ( $fileroot, $rcffile, $radoctfile, $riffile, $skyfile, $radmatfile, $radmatcopy, $diffskyfile ) =
    tellradnames( $conffile, $path, $radpath );
    #say "\$fileroot, $fileroot, \$rcffile, $rcffile, \$radoctfile, $radoctfile, \$riffile, $riffile, \$skyfile, $skyfile, \$radmatfile, $radmatfile, \$radmatcopy, $radmatcopy, \$diffskyfile, $diffskyfile\n\n";

    my $conffile_a = $conffile;
    $conffile_a =~ s/\.cfg$/_a\.cfg/ ;
    print REPORT "cp -f $conffile $conffile_a\n";
    `cp -f $conffile $conffile_a`;

    my ( $fileroot_a, $rcffile_a, $radoctfile_a, $riffile_a, $skyfile_a, $radmatfile_a, $radmatcopy_a, $diffskyfile_a ) =
    tellradnames( $conffile_a, $path, $radpath );
    #say "\$fileroot_a, $fileroot_a, \$rcffile_a, $rcffile_a, \$radoctfile_a, $radoctfile_a, \$riffile_a, $riffile_a, \$skyfile_a, $skyfile_a, \$radmatfile_a, $radmatfile_a, \$radmatcopy_a, $radmatcopy_a, \$diffskyfile_a, $diffskyfile_a\n\n";

    my $countmonth = 0;
    foreach my $monthsdataref ( @daylighthoursarr )
    {
      my @monthsdata = @$monthsdataref;
      my @hourslight = @{$monthsdata[0][0]};
      my $month = $monthsdata[1];
      $month =~ s/`//g;
      my $monthnum = getmonthnum( $month );
      my $day = $days_inmonths{ $month };

      my $countsurf = 0;
      foreach my $surfref ( @gridpoints )
      {
        my @surfdata = @$surfref;
        my @pointrefs = @{ $surfdata[0] };
        my $surfnum = $surfdata[1];
        my @dirvector = @{ $surfdata[2] };
        my ( $dirvectorx, $dirvectory, $dirvectorz ) = @dirvector;
        my $countlithour = 0;
        my $counthour = 0;
        foreach my $hourlight ( @hourslight )
        {
          if ( $hourlight != 1 )
          {
            $countlithour++;
            my $hour = ( $counthour + 1) ;

            my $dir = $avgs{dir}{$monthnum}{$hour}; say REPORT "HERE \$dir: $dir";
            my $diff = $avgs{diff}{$monthnum}{$hour}; say REPORT "HERE \$diff: $diff";
            my $alt = $avgs{alt}{$monthnum}{$hour}; say REPORT "HERE \$alt: $alt";
            my $azi = $avgs{azi}{$monthnum}{$hour}; say REPORT "HERE \$azi: $azi";

            my $countpoint = 0;

            if ( ( $countmonth == 0 ) and ( $countsurf == 0 ) and ( $countlithour == 1 ) )
            {
              if ( $countrad >= 1 )
              {
                my ( $conffileother1, $filerootother1, $rcffileother1, $radoctfileother1, $riffileother1, $skyfileother1, $radmatfileother1, $radmatcopyother1, $diffskyfileother1 );
                #my ( $conffileother2, $filerootother2, $rcffileother2, $radoctfileother2, $riffileother2, $skyfileother2, $radmatfileother2, $radmatcopyother2, $diffskyfileother2 );
                if ( ( $countrad == 2 ) and ( "coexistent" ~~ @calcprocedures ) )
                {
                  $conffileother1 = $conffiles[0];
                  $conffileother1 = $conffiles[1];

                  ( $filerootother1, $rcffileother1, $radoctfileother1, $riffileother1, $skyfileother1, $radmatfileother1, $radmatcopyother1, $diffskyfileother1 ) = tellradnames( $conffileother1, $path, $radpath );
                  #( $filerootother2, $rcffileother2, $radoctfileother2, $riffileother2, $skyfileother2, $radmatfileother2, $radmatcopyother2, $diffskyfileother2 ) = tellradnames( $conffileother2, $path, $radpath );
                  `cp -f $radmatfileother1 $radmatfile`;
                  say REPORT "$radmatfileother1 $radmatfile";
                }
              }

              setroot( $conffile, $path, $debug, \%paths );

              setrad( $conffile, $radoctfile, $rcffile, $path, $radpath, $monthnum,
              $day, $hour, $countfirst, $exportconstrref, $exportreflref, \%skycondition, $countrad,
              \@specularratios, \@calcprocedures, $debug, \%paths );

              if ( $countrad >= 1 )
              {
                my ( $conffileother1, $filerootother1, $rcffileother1, $radoctfileother1, $riffileother1, $skyfileother1, $radmatfileother1, $radmatcopyother1, $diffskyfileother1 );
                #my ( $conffileother2, $filerootother2, $rcffileother2, $radoctfileother2, $riffileother2, $skyfileother2, $radmatfileother2, $radmatcopyother2, $diffskyfileother2 );
                unless ( ( $countrad == 2 ) and ( "coexistent" ~~ @calcprocedures ) )
                {
                  adjust_radmatfile( $exportconstrref, $exportreflref, $conffile, $path, \@specularratios,
                  \%obslayers, \@selectives, \%paths );
                  say REPORT "cp -f $radmatfile $radmatcopy";
                  `cp -f $radmatfile $radmatcopy`;

                  open ( FIXLIST, ">$path/rad/fixl.pl" ) or die( $! );
                  print FIXLIST "$radmatfile\n";
                  print FIXLIST "$radmatcopy\n";
                  close FIXLIST;
                }
                else
                {
                  $conffileother1 = $conffiles[0];
                  $conffileother1 = $conffiles[1];

                  ( $filerootother1, $rcffileother1, $radoctfileother1, $riffileother1, $skyfileother1, $radmatfileother1, $radmatcopyother1, $diffskyfileother1 ) = tellradnames( $conffileother1, $path, $radpath );
                  #( $filerootother2, $rcffileother2, $radoctfileother2, $riffileother2, $skyfileother2, $radmatfileother2, $radmatcopyother2, $diffskyfileother2 ) = tellradnames( $conffileother2, $path, $radpath );
                  `cp -f $radmatfileother1 $radmatfile`;
                  say REPORT "$radmatfileother1 $radmatfile";
                }
              }
              $setoldfiles = "off";
            }


            ## HERE THE CALLS FOR THE DIFFUSE IRRADIANCES FOLLOW.
            unless ( ( "gendaylit" ~~ @calcprocedures ) or ( "gensky" ~~ @calcprocedures ) )
            {
              setrad( $conffile, $radoctfile, $rcffile, $path, $radpath, $monthnum, $day, $hour,
              $countfirst, $exportconstrref, $exportreflref, \%skycondition, $countrad, \@specularratios,
              \@calcprocedures, $debug, \%paths );
            }

            my $skycond = $skycondition{$monthnum};
            my ( @returns, $altreturns );

            #say "\$month:$monthnum, \$hour:$hour, \$dir:$dir, \$diff:$diff, \$alt:$alt, \$azi:$azi, \$lat:$lat, \$long:$long, \$standardmeridian:$standardmeridian";
            if ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) )
            {
              if ( $alt <= 0 ) { $alt = 0.0001; say "IMPOSED \$alt = 0.0001;"; say REPORT "IMPOSED \$alt = 0.0001;"; } # IMPORTANT: THIS SETS THE ALTITUDE > 0 OF A TINY AMOUNT IF IT IS < 0 DUE TO THE FACT
              # THAT THE MAJORITY OF THAT HOUR THE SUN WAS BELOW THE HORIZON, WHILE THE NET GAINED AMOUNT OF RADIATION WAS STILL > 0.
              if ( $dir == 0 ){ $dir = 0.0001; say "IMPOSED \$dir = 0.0001;"; say REPORT "IMPOSED \$dir = 0.0001;"; }
              if ( $diff == 0 ){ $diff = 0.0001; say "IMPOSED \$diff = 0.0001;"; say REPORT "IMPOSED \$diff = 0.0001;"; }
            }
            # IMPORTANT: THE TWO LINES ABOVE SET THE DIFFUSE AND DIRECT IRRADIANCE > 0 OF A TINY AMOUNT IF THERE ARE 0
            # TO AVOID ERRORS IN THE rtrace CALLS WHEN THE ALTITUDE IS > 0.

            my $altreturn;


            if ( ( "getweather" ~~ @calcprocedures )
              or ( ( "gensky" ~~ @calcprocedures ) and ( "altcalcdiff" ~~ @calcprocedures ) )
              or ( ( "gendaylit" ~~ @calcprocedures ) and ( "altcalcdiff" ~~ @calcprocedures )  ) ) # IF CONDITION, PRODUCE LIGHT WITH NO SUN TO CHECK DIFFUSE RADIATION
            {
              if ( "gendaylit" ~~ @calcprocedures )
              {
                if ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) )
                {
                  @returns = `gendaylit -ang $alt $azi +s -g $groundrefl -W 0 $diff -a $lat -o $long -m $standardmeridian`;
                  say REPORT "gendaylit -ang $alt $azi +s -g $groundrefl -W 0 $diff -a $lat -o $long -m $standardmeridian";
                }
                else
                {
                  if ( not( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) ) )
                  {
                    @returns = `gendaylit $monthnum $day $hour -s -g $groundrefl -a $lat -o $long -m $standardmeridian`;
                    say REPORT "gendaylit $monthnum $day $hour -s -g $groundrefl -a $lat -o $long -m $standardmeridian";
                  }
                  elsif ( ( "getweather" ~~ @calcprocedures ) and ( not( "getsimple" ~~ @calcprocedures ) ) )
                  {
                    @returns = `gendaylit -ang $alt $azi -s -g $groundrefl -W $dir $diff -a $lat -o $long -m $standardmeridian`;
                    say REPORT "gendaylit -ang $alt $azi -s -g $groundrefl -W $dir $diff -a $lat -o $long -m $standardmeridian";
                  }
                }
              }
              elsif ( "gensky" ~~ @calcprocedures )
              {
                if ( $skycond eq "clear" )
                {
                  @returns = `gensky $monthnum $day $hour -s -g $groundrefl -a $lat -o $long -m $standardmeridian`;
                  say REPORT "gensky $monthnum $day $hour -s -g $groundrefl -a $lat -o $long -m $standardmeridian";
                }
                elsif ( $skycond eq "cloudy" )
                {
                  @returns = `gensky $monthnum $day $hour -i -g $groundrefl -a $lat -o $long -m $standardmeridian`;
                  say REPORT "gensky $monthnum $day $hour -i -g $groundrefl -a $lat -o $long -m $standardmeridian";
                }
                elsif ( $skycond eq "overcast" )
                {
                  @returns = `gensky $monthnum $day $hour -c -g $groundrefl -a $lat -o $long -m $standardmeridian`;
                  say REPORT "gensky $monthnum $day $hour -c -g $groundrefl -a $lat -o $long -m $standardmeridian";
                }
                #say "gensky $monthnum $day $hour -s -g $groundrefl -a $lat -o $long -m $standardmeridian";
              }

              my $newline;
              my @alts = split( " +", $returns[2] ); #say REPORT "ALTS: " . dump ( @alts );
              $altreturn = $alts[5]; #say REPORT "\$altreturn: " . dump ( $altreturn );
              if ( $altreturn <= 0 )
              {
                $newline = "$alts[0] $alts[1] $alts[2] $alts[3] $alts[4] 0.0001 $alts[6]";
              }


              if ( ( "gendaylit" ~~ @calcprocedures ) or ( "gensky" ~~ @calcprocedures ) )
              {
                say REPORT "rm $skyfile";
                open( SKYFILE, ">$skyfile" ) or die;

                #say "ALT: $alt\n";
                my $counter = 0;
                foreach my $line ( @returns )
                {
                  if ( $counter == 2 )
                  {
                    if ( $newline ne "" )
                    {
                      print SKYFILE $newline;
                      #print REPORT "SUBSTITUTED:" . $newline;
                    }
                    else
                    {
                      print SKYFILE $line;
                      #print REPORT "CONFIRMED:" . $line;
                    }
                  }
                  else
                  {
                    print SKYFILE $line;
                    #print REPORT "CONFIRMED:" . $line;
                  }
                  $counter++;
                }

                if ( not ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) ) )
                {

print SKYFILE "

void light solar
0
0
3 0.0 0.0 0.0

solar source sun
0
0
4 0.0 0.0 0.0 0.0
";
close SKYFILE;


print REPORT"

void light solar
0
0
3 0.0 0.0 0.0

solar source sun
0
0
4 0.0 0.0 0.0 0.0
";
                }
              }

              say REPORT "IN CALCULATIONS FOR DIFFUSE RADIATION. \$hour: $hour, \$surfnum: $surfnum, \$month: $month";

              open( SKYFILE, "$skyfile" ) or die;
              my @lines = <SKYFILE>;
              foreach my $line ( @lines )
              {
                print REPORT "$line";
              }
              close SKYFILE;

              my @shuffled;
              if ( ( $countrad == 0 ) and ( "radical" ~~ @calcprocedures ) and ( "sparecycles" ~~ @calcprocedures ) )
              {
                @shuffleds = shuffle( @pointrefs );
                @shuffleds = @shuffleds[0..1];
              }

              #say "\$countrad $countrad";
              unless ( ( "coexistent" ~~ @calcprocedures ) and ( $countrad == 2 ) )
              {
                my $countpoint = 0;
                #my $pm = Parallel::ForkManager->new( $max_processes ); #Sets up the possibility of opening child processes
                unless ( ( $altreturn < 0 ) and ( ( "gensky" ~~ @calcprocedures ) or ( "gendaylit" ~~ @calcprocedures ) )
                    and ( not ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) ) ) )
                {
                  foreach my $pointref ( @pointrefs )
                  {
                    my $shuffled = $shuffleds[$countpoint];
                    if ( $shuffled ne "" )
                    {
                      $pointref = $shuffled;
                    }
                    #$pm->start and next; # Begins the child process
                    my @pointcoords = @$pointref;
                    my ( $xcoord, $ycoord, $zcoord ) = @pointcoords;

                    #my $raddir = "$path/rad/"; #OVERRIDDEN BELOW
                    my $raddir = $paths{radpath};
                    #my $cfgpath = "$path/cfg/"; #OVERRIDDEN BELOW
                    my $cfgpath = $paths{cfgpath};

                    my @dirvgroup = getdirvectors ( \@basevectors, \@dirvector );


                    my $countdirvec = 0;
                    foreach my $dirvector ( @dirvgroup )
                    {
                      my ( $valstring, $valstring1, $valstring2, $irr, $irr1, $irr2 );
                      my ( $dirvx, $dirvy, $dirvz ) = @{ $dirvector };

                      $valstring = `cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab $bounceambnum -lr $bouncemaxnum -h $radoctfile`;
                      say REPORT "5TO SHELL: cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab $bounceambnum -lr $bouncemaxnum -h $radoctfile";
                      if ( $valstring =~ /fatal/ ){ die ; };
                      my ( $x, $y, $z ) = ( $valstring =~ m/(.+)\t(.+)\t(.+)\t/ );
                      $irr = ( 179 * ( ( .265 * $x ) + ( .670 * $y ) + ( .065 * $z ) ) );
                      push ( @{ $surftestsdiff{$radoctfile}{$monthnum}{$surfnum}{$hour} }, $irr );
                      $countdirvec++;
                    }
                    $countpoint++;
                    #$pm->finish;
                  }
                }
                else
                {
                  push ( @{ $surftestsdiff{$radoctfile}{$monthnum}{$surfnum}{$hour} }, 0 );
                }
              }
            }


            ## HERE FOLLOW THE OPERATIONS FOR THE TOTAL AND DIRECT IRRADIANCES.

            ## HERE FOLLOWS THE CREATION OF THE SKIES FOR OBTAINING THE TOTAL AND DIRECT IRRADIANCES.
            unless ( "keepdirshdf" ~~ @calcprocedures )
            {
              if ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) ) # IF CONDITION, CHECK DIRECT RADIATION
              {
                @returns = `gendaylit -ang $alt $azi +s -g $groundrefl -W $dir 0 -a $lat -o $long -m $standardmeridian`;
                say REPORT "gendaylit -ang $alt $azi +s -g $groundrefl -W $dir 0 -a $lat -o $long -m $standardmeridian";
              }
              elsif ( ( "gendaylit" ~~ @calcprocedures ) and ( "getweather" ~~ @calcprocedures ) and ( not ( ( "getsimple" ~~ @calcprocedures ) ) ) )
              {
                @returns = `gendaylit -ang $alt $azi +s -g $groundrefl -W $dir $diff -a $lat -o $long -m $standardmeridian`;
                say REPORT "gendaylit -ang $alt $azi +s -g $groundrefl -W $dir $diff -a $lat -o $long -m $standardmeridian";
              }
              elsif ( ( ( not( "getweather" ~~ @calcprocedures ) ) and ( not( "getsimple" ~~ @calcprocedures ) ) )
                or ( "gensky" ~~ @calcprocedures ) )
              {
                if ( $skycond eq "clear" )
                {
                  @returns = `gensky $monthnum $day $hour +s -g $groundrefl -a $lat -o $long -m $standardmeridian`;
                  say REPORT "gensky $monthnum $day $hour +s -g $groundrefl -a $lat -o $long -m $standardmeridian";
                }
                elsif ( $skycond eq "cloudy" )
                {
                  @returns = `gensky $monthnum $day $hour +i -g $groundrefl -a $lat -o $long -m $standardmeridian`;
                  say REPORT "gensky $monthnum $day $hour +i -g $groundrefl -a $lat -o $long -m $standardmeridian";
                }
                elsif ( $skycond eq "overcast" )
                {
                  @returns = `gensky $monthnum $day $hour -c -g $groundrefl -a $lat -o $long -m $standardmeridian`;
                  say REPORT "gensky $monthnum $day $hour -c -g $groundrefl -a $lat -o $long -m $standardmeridian";
                }
              }

              say REPORT "IN CALCULATIONS FOR TOTAL AND/OR DIRECT RADIATION. \$hour: $hour, \$surfnum: $surfnum, \$month: $month";

              my $newline;
              my @alts = split( " +", $returns[2] ); #say REPORT "ALTS: " . dump ( @alts );
              $altreturn = $alts[5]; #say REPORT "\$altreturn: " . dump ( $altreturn );
              if ( $altreturn <= 0 )
              {
                $newline = "$alts[0] $alts[1] $alts[2] $alts[3] $alts[4] 0.0001 $alts[6]";
              }

              if ( ( ( "gendaylit" ~~ @calcprocedures ) ) or ( "gensky" ~~ @calcprocedures ) )
              {
                say REPORT "rm -f $skyfile";

                open( SKYFILE, ">$skyfile" ) or die;

                my $counter = 0;
                my $counter = 0;
                foreach my $line ( @returns )
                {
                  if ( $counter == 2 )
                  {
                    if ( $newline ne "" )
                    {
                      print SKYFILE $newline;
                      #print REPORT "SUBSTITUTED:" . $newline;
                    }
                    else
                    {
                      print SKYFILE $line;
                      #print REPORT "CONFIRMED:" . $line;
                    }
                  }
                  else
                  {
                    print SKYFILE $line;
                    #print REPORT "CONFIRMED:" . $line;
                  }
                  $counter++;
                }

                close SKYFILE;
              }
            }


            ## HERE FOLLOW THE OPERATIONS FOR THE TOTAL IRRADIANCES.
            #say "\$countrad $countrad";
            unless ( ( ( "coexistent" ~~ @calcprocedures ) and ( $countrad == 0 ) )
              or ( "keepdirshdf" ~~ @calcprocedures ) )
            {
              unless ( ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) )
                or ( ( "directlydirect" ~~ @calcprocedures ) and ( "altcalcdiff" ~~ @calcprocedures ) ) )
              {
                my $countpoint = 0;
                #my $pm4 = Parallel::ForkManager->new( $max_processes ); #Sets up the possibility of opening child processes
    	           #DATA_LOOP:
                unless ( ( $altreturn < 0 ) and ( ( "gensky" ~~ @calcprocedures ) or ( "gendaylit" ~~ @calcprocedures ) ) )
                {
                  foreach my $pointref ( @pointrefs )
                  {
                    #my $pid4 = $pm4->start and next DATA_LOOP; # Begins the child process
                    my @pointcoords = @$pointref;
                    my ( $xcoord, $ycoord, $zcoord ) = @pointcoords;
                    #my $raddir = "$path/rad/"; # GOING TO BE OVERRIDDEN BELOW
                    my $raddir = $radpath;
                    #my $cfgpath = "$path/cfg/";  # GOING TO BE OVERRIDDEN BELOW
                    my $cfgpath = $paths{cfgpath};
                    my @dirvgroup = getdirvectors ( \@basevectors, \@dirvector );

                    my $countdirvec = 0;
                    foreach my $dirvector ( @dirvgroup )
                    { #say REPORT "IN RTRACE 5";
                      my ( $valstring, $valstring1, $valstring2, $irr, $irr1, $irr2 );
                      my ( $dirvx, $dirvy, $dirvz ) = @{ $dirvector };
                      $valstring = `cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab $bounceambnum -lr $bouncemaxnum -h $radoctfile`;
                      say REPORT "5TO SHELL: cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab $bounceambnum -lr $bouncemaxnum -h $radoctfile";
      		            my ( $x, $y, $z ) = ( $valstring =~ m/(.+)\t(.+)\t(.+)\t/ );
                      $irr = ( 179 * ( ( .265 * $x ) + ( .670 * $y ) + ( .065 * $z ) ) );
                      push ( @{ $surftests{$radoctfile}{$monthnum}{$surfnum}{$hour} }, $irr );

                      $countdirvec++;
                    }
                    $countpoint++;
      	            #$pm4->finish;
                  }
                  #say "Obtained total @{ $surftests{$radoctfile}{$monthnum}{$surfnum}{$hour} }";
                }
                else
                {
                  push ( @{ $surftests{$radoctfile}{$monthnum}{$surfnum}{$hour} }, 0 );
                }
              }
            }



            ## HERE FOLLOW THE OPERATIONS FOR THE DIRECT IRRADIANCES.
            unless ( ( "coexistent" ~~ @calcprocedures ) and ( $countrad == 0 ) )
            {
              if ( "espdirres" ~~ @calcprocedures )
              {
                @gridpoints = @orgridpoints;
              }
              unless ( "keepdirshdf" ~~ @calcprocedures ) # IF CONDITION, CHECK DIRECT RADIATION
              {
                if ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) )
                {
                  my $countpoint = 0;
                  #my $pm4 = Parallel::ForkManager->new( $max_processes ); #Sets up the possibility of opening child processes
      	          #DATA_LOOP:
                  unless ( $alt <= 0 )
                  {
                    foreach my $pointref ( @pointrefs )
                    {
                      #my $pid4 = $pm4->start and next DATA_LOOP; # Begins the child process
                      my @pointcoords = @$pointref;
                      my ( $xcoord, $ycoord, $zcoord ) = @pointcoords;
                      #my $raddir = "$path/rad/"; # GOING TO BE OVERRIDDEN BELOW
                      my $raddir = $radpath;
                      #my $cfgpath = "$path/cfg/";  # GOING TO BE OVERRIDDEN BELOW
                      my $cfgpath = $paths{cfgpath};
                      my @dirvgroup = getdirvectors ( \@basevectors, \@dirvector );

                      my $countdirvec = 0;
                      foreach my $dirvector ( @dirvgroup )
                      { #say REPORT "IN RTRACE 4";
                        my ( $valstring, $valstring1, $valstring2, $irr, $irr1, $irr2 );
                        my ( $dirvx, $dirvy, $dirvz ) = @{ $dirvector };
                        $valstring = `cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab $bounceambnum -lr $bouncemaxnum -h $radoctfile`;
                        say REPORT "5TO SHELL: cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab $bounceambnum -lr $bouncemaxnum -h $radoctfile";
                        if ( $valstring =~ /fatal/ ){ die; };
                        my ( $x, $y, $z ) = ( $valstring =~ m/(.+)\t(.+)\t(.+)\t/ );
                        $irr = ( 179 * ( ( .265 * $x ) + ( .670 * $y ) + ( .065 * $z ) ) );
                        push ( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} }, $irr );

                        $countdirvec++;
                      }
                      $countpoint++;
        	            #$pm4->finish;
                    }
                    #say "Obtained direct @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} }";
                  }
                  else
                  {
                    push ( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} }, 0 );
                  }
                }
                elsif ( ( "directlydirect" ~~ @calcprocedures ) or ( ( not( "directlydirect" ~~ @calcprocedures ) ) and ( not( "altcalcdiff" ~~ @calcprocedures ) ) ) )
                { #say "DIRECT CALCDIFF";
                  my $countpoint = 0;
                  #my $pm = Parallel::ForkManager->new( $max_processes ); #Sets up the possibility of opening child processes
                  unless ( ( $altreturn < 0 ) and ( ( "gensky" ~~ @calcprocedures ) or ( "gendaylit" ~~ @calcprocedures ) ) )
                  {
                    foreach my $pointref ( @pointrefs )
                    {
                      #$pm->start and next; # Begins the child process
                      my @pointcoords = @$pointref;
                      my ( $xcoord, $ycoord, $zcoord ) = @pointcoords;

                      #my $raddir = "$path/rad/"; #OVERRIDDEN BELOW
                      my $raddir = $paths{radpath};
                      #my $cfgpath = "$path/cfg/"; #OVERRIDDEN BELOW
                      my $cfgpath = $paths{cfgpath};

                      my @dirvgroup = getdirvectors ( \@basevectors, \@dirvector );

                      my $countdirvec = 0;
                      foreach my $dirvector ( @dirvgroup )
                      { #say REPORT "IN RTRACE 7a";
                        my ( $valstring, $valstring1, $valstring2, $irr, $irr1, $irr2 );
                        my ( $dirvx, $dirvy, $dirvz ) = @{ $dirvector };
                        $valstring = `cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab 0 -av 0 0 0 -lr $bouncemaxnum -h $radoctfile`;
                        say REPORT "5TO SHELL: cd $raddir \n echo $xcoord $ycoord $zcoord $dirvx $dirvy $dirvz | rtrace  -I -ab 0 -av 0 0 0 -lr $bouncemaxnum -h $radoctfile";
                        my ( $x, $y, $z ) = ( $valstring =~ m/(.+)\t(.+)\t(.+)\t/ );
                        $irr = ( 179 * ( ( .265 * $x ) + ( .670 * $y ) + ( .065 * $z ) ) );
                        push ( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} }, $irr );
                        $countdirvec++;
                      }
                      $countpoint++;
                      #$pm->finish;
                    }
                    #$pm->wait_all_children;
                  }
                  else
                  {
                    push ( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} }, 0 );
                  }
                  #say "Obtained direct @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} }";
                }
              }
            }


            my ( $meanvaluesurf, $meanvaluesurf_diff, $meanvaluesurf_dir, $meantoptot, $meantopdir, $meantopdiff );


            if ( ( ( "getweather" ~~ @calcprocedures ) and ( "getsimple" ~~ @calcprocedures ) ) or ( ( "altcalcdiff" ~~ @calcprocedures ) and ( "directlydirect" ~~ @calcprocedures ) ) )
            {
              if ( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} } )
              {
                $meanvaluesurf_dir = mean( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} } );
              }

              if ( @{ $surftestsdiff{$radoctfile}{$monthnum}{$surfnum}{$hour} } )
              {
                $meanvaluesurf_diff = mean( @{ $surftestsdiff{$radoctfile}{$monthnum}{$surfnum}{$hour} } );
              }

              $meanvaluesurf = ( $meanvaluesurf_diff + $meanvaluesurf_dir );
            }
            elsif ( ( "altcalcdiff" ~~ @calcprocedures ) and ( not ( "directlydirect" ~~ @calcprocedures ) ) )
            {
              if ( @{ $surftests{$radoctfile}{$monthnum}{$surfnum}{$hour} } )
              {
                $meanvaluesurf = mean( @{ $surftests{$radoctfile}{$monthnum}{$surfnum}{$hour} } );
              }

              if ( @{ $surftestsdiff{$radoctfile}{$monthnum}{$surfnum}{$hour} } )
              {
                $meanvaluesurf_diff = mean( @{ $surftestsdiff{$radoctfile}{$monthnum}{$surfnum}{$hour} } );
              }

              $meanvaluesurf_dir = ( $meanvaluesurf - $meanvaluesurf_diff );
            }
            elsif ( ( ( "directlydirect" ~~ @calcprocedures ) and ( not ( "altcalcdiff" ~~ @calcprocedures ) ) )
              or ( ( not( "directlydirect" ~~ @calcprocedures ) ) and ( not( "altcalcdiff" ~~ @calcprocedures ) ) ) )
            {
              if ( @{ $surftests{$radoctfile}{$monthnum}{$surfnum}{$hour} } )
              {
                $meanvaluesurf = mean( @{ $surftests{$radoctfile}{$monthnum}{$surfnum}{$hour} } );
              }

              if ( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} } )
              {
                $meanvaluesurf_dir = mean( @{ $surftestsdir{$radoctfile}{$monthnum}{$surfnum}{$hour} } );
              }

              $meanvaluesurf_diff = ( $meanvaluesurf - $meanvaluesurf_dir );
            }


            if ( $meanvaluesurf_diff < 0 ) { $meanvaluesurf_diff = 0; };
            if ( $meanvaluesurf_dir < 0 ) { $meanvaluesurf_dir = 0; };
            if ( $meantopdiff < 0 ) { $meantopdiff = 0; };
            if ( $meantopdir < 0 ) { $meantopdir = 0; };

            if ( ( $meanvaluesurf_diff + $meanvaluesurf_dir ) > $meanvaluesurf )
            {
              $meanvaluesurf = ( $meanvaluesurf_diff + $meanvaluesurf_dir );
            }

            say "Surface $surfnum, zone $zonenum, month $monthnum, day $day, hour $hour, octree $radoctfile.";
            say "Diffuse irradiance: $meanvaluesurf_diff. ";
            say "Direct irradiance: $meanvaluesurf_dir. ";
            say "Total irradiance: $meanvaluesurf.\n";

            say REPORT "Surface $surfnum, zone $zonenum, month $monthnum, day $day, hour $hour, octree $radoctfile.";
            say REPORT "Diffuse irradiance: $meanvaluesurf_diff. ";
            say REPORT "Direct irradiance: $meanvaluesurf_dir. ";
            say REPORT "Total irradiance: $meanvaluesurf.\n";

            if ( $meanvaluesurf_diff and $surfnum and $hour )
            {
              $irrs{ $zonenum }{ $countrad + 1 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr } = $meanvaluesurf_diff;
            }

            if ( $meanvaluesurf_dir and $surfnum and $hour )
            {
              $irrs{ $zonenum }{ $countrad + 1 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr } = $meanvaluesurf_dir;
            }
          }
          $counthour++;
        }
        $countsurf++;
      }
      $countmonth++;
    }
    $countrad++;
  }
  $" = ",";
  return ( \%irrs );
}

sub cleanblanks
{ # SELF-EXPLAINING.
  my @arr = @_;
  my @box;
  my $count = 0;
  foreach my $el ( @arr )
  {
    unless ( $el eq "" )
    {
      push ( @box, $el );
    }
    $count++;
  }
  return( @box );
}

sub createconstrdbfile
{ # THIS CREATES THE CONSTRUCTION DB FILE OF THE FICTITIOUS ESP-r MODELS.
  my ( $constrdbfile, $constrdbfile_f, $obsconstrsetref, $calcprocedures_ref ) = @_;
  my @obsconstrset = @$obsconstrsetref;
  my @calcprocedures = @{ $calcprocedures_ref };
  @obsconstrset = uniq( @obsconstrset );
  my ( @bigcopy, @updatedlines );
  open ( DBFILE, "$constrdbfile" ) or die;
  my @lines = <DBFILE>;
  close DBFILE;

#  my $topline; this is not used
  my $countline = 0;

if ( "oldconstrdb" ~~ @calcprocedures )
{
  # --- OLD CONSTR DATABASE ---
     foreach my $line ( @lines )
     { # THIS PUSHES IN @UPDATEDLINES THE CONSTR DATABASE EXCEPTED THE FIRST LINES (HEADER LINES)
       # This actually pushes all lines, including header lines, which explains the "if ( $countline > 2 )" below.
       my @row = split( /\s+|,/ , $line);
       @row = cleanblanks( @row );
       my ( $oldnumber, $newnumber );
       my $atline;
       if ( $line =~ /\# no of composites/ )
       {
         $atline == $countline;
         $oldnumber = $row[0];
         $newnumber = ( $oldnumber + scalar( @obsconstrset ) );
         $line =~ s/$oldnumber/$newnumber/;
         push ( @updatedlines, $line );
       }
       else
       {
         push ( @updatedlines, $line );
       }
       $countline++;
     } # --- END OLD CONSTR DATABASE ---
  }
  elsif ( not ( "oldconstrdb" ~~ @calcprocedures ) )
  {
    # --- NEW CONSTR DATABASE ---
    # Push database contents into @updatedlines.
    # Add a category called "Modish_fict" while doing this.
    foreach my $line ( @lines )
    {
      push ( @updatedlines, $line );
      if ( $line =~ /^\*date,/ )
      {
        push ( @updatedlines, "*Category,Modish_fict,Modish fictitious constructions,fictitious versions of existing constructions used for shading factor modifier script Modish\n" );
      }
    } # --- END NEW CONSTR DATABASE ---
  }

  if ( "oldconstrdb" ~~ @calcprocedures )
  {
    # --- OLD CONSTR DATABASE ---
    my $coun = 0;
    foreach my $el ( @obsconstrset )
    { #FOREARCH MATERIAL USED IN THE OBSTRUCTIONS, PUSHES THE CONSTRUCTION SOLUTIONS IN WHICH IT IS USED IN @COPY, AND PUSHES EACH [ @COPY ] IN @BIGCOPY
      my @copy;
      my $semaphore = 0;
      $countel = 0;
      $countline = 0;
      foreach my $line ( @updatedlines )
      {
        my @row = split( /\s+|,/ , $line);
        @row = cleanblanks( @row );
       # if ( $countline > 2 ) #WHY IS THIS?
       # {
          if ( $el eq $row[1] )
          {
            $semaphore = 1;
          }

          if ( ( $semaphore == 1 ) and ( $countel == 0) )
          {
            push ( @copy, "# layers  description   optics name   symmetry tag\n" );
            push ( @copy, $line );
            $countel++;
          }
          elsif ( ( $semaphore == 1 ) and ( $countel > 0) )
          {
            push ( @copy, $line );
            if (  ( $row[0] eq "#" ) and ( $row[1] eq "layers" ) and ( $row[2] eq "description" )  )
            {
              pop(@copy);
              $semaphore = 0;
            }
            $countel++;
          }

        #}
        $countline++;
      }
      #say REPORT "\@copy " . dump(@copy);
      #if ( $coun == $#obsconstrset )
      #{
      #  @bigcopy = [ @copy ];
      #}
      push ( @bigcopy, [ @copy ] );
      $coun++;
    } # --- END OLD CONSTR DATABASE ---
  }
  elsif ( not ( "oldconstrdb" ~~ @calcprocedures ) )
  {
    # --- NEW CONSTR DATABASE ---
    # If there are obstruction constructions, loop through each @updatedlines.
    # If an obstruction construction is found, push this into @copy, and push each [ @copy ] into @bigcopy.
    my $semaphore = 0;
    my $nummatches = 0;
    my $numobsconstr = scalar @obsconstrset;
    my @copy;
    if ( $numobsconstr > 0 )
    {
      foreach my $line ( @updatedlines )
      {
        my @row = split( /,/ , $line);
        @row = cleanblanks( @row );
        if ( ( $semaphore == 0 ) and ( $row[0] == "*item" ) and ( any { $_ eq $row[1] } @obsconstrset ) )
        {
          $semaphore = 1;
        }
        if ( $semaphore == 1 )
        {
          push ( @copy, $line );
          if ( $row[0] =~ /\*end_item/ )
          {
            $semaphore = 0;
            push ( @bigcopy, [ @copy ] );
            undef ( @copy );
            $nummatches++;
            if ( $nummatches == $numobsconstr ) { last }
          }
        }
      }
    }# --- END NEW CONSTR DATABASE ---
  }


  if ( "oldconstrdb" ~~ @calcprocedures )
  {
  # --- OLD CONSTR DATABASE ---
    my $cn = 0;
     my ( @materials, @newmaterials, @newcopy, @newbigcopy );
     my ( $newmatinssurf, $newmatextsurf );
     my %exportconstr;
     foreach my $copyref ( @bigcopy )
     {
       my @constrlines = @$copyref;
       my $firstrow = $constrlines[1];
       my @row = split ( /\s+|,/ , $firstrow );
       @row = cleanblanks( @row );
       my $constrname = $row[1];
       my $newconstrname = $constrname;
       $newconstrname =~ s/\w\b// ;
       $newconstrname =~ s/\w\b// ;
       $newconstrname = "f_" . "$newconstrname";

       my $intlayer = $constrlines[3];
       my @row = split ( /\s+|,/ , $intlayer );
       @row = cleanblanks( @row );
       my $matintlayer = $row[2];
       my $newmatintlayer = $matintlayer;
       $newmatintlayer =~ s/\w\b// ;
       $newmatintlayer =~ s/\w\b// ;
       $newmatintlayer = "f_" . "$newmatintlayer";
       my $extlayer = $constrlines[$#constrlines];
       my @row = split ( /\s+|,/ , $extlayer );
       @row = cleanblanks( @row );
       my $matextlayer = $row[2];
       my $newmatextlayer = $matextlayer;
       $newmatextlayer =~ s/\w\b// ;
       $newmatextlayer =~ s/\w\b// ;
       $newmatextlayer = "f_" . "$newmatextlayer";
       push ( @materials, $matintlayer, $matextlayer );
       push ( @newmaterials, $newmatintlayer, $newmatextlayer );

       $constrlines[1] =~ s/$constrname/$newconstrname/g;
       $constrlines[3] =~ s/$matintlayer/$newmatintlayer/g;
       $constrlines[$#constrlines] =~ s/$matextlayer/$newmatextlayer/g;
       foreach my $line ( @constrlines )
       {
         push ( @newcopy, $line );
       }
       @newbigcopy = [ @newcopy ] ;
       #push ( @newbigcopy, [ @newcopy ] );
       $exportconstr{ $newconstrname }{ extlayer } = $newmatextlayer;
       $exportconstr{ $newconstrname }{ intlayer } = $newmatintlayer;
       $cn++;
     } # --- END OLD CONSTR DATABASE ---
  }
  elsif ( not ( "oldconstrdb" ~~ @calcprocedures ) )
  {
    # --- NEW CONSTR DATABASE ---
    # In each [ @copy ], modify the:
    # construction name and documentation,
    # category,
    # internal material name, and
    # external material name.
    # Store the old and new material names, and form a hash relating new materials to the new constructions.
    my ( @materials, @newmaterials, @newcopy, @newbigcopy );
    my %exportconstr;
    foreach my $copyref ( @bigcopy )
    {
      my @constrlines = @$copyref;
      my $onlyonelayer = 0;
      if ( $#constrlines == 5 ) { $onlyonelayer = 1 }

      my $firstrow = $constrlines[0];
      my @row = split ( /,/ , $firstrow );
      @row = cleanblanks( @row );
      my $constrname = $row[1];
      my $newconstrname = $constrname;
      if ( length($constrname) > 30 )
      {
        $newconstrname =~ s/\w\b// ;
        $newconstrname =~ s/\w\b// ;
      }
      $newconstrname = "f_" . "$newconstrname";

      my $intlayer = $constrlines[4];
      my @row = split ( / : |,/ , $intlayer );
      @row = cleanblanks( @row );
      my $matintlayer = $row[3];
      my $newmatintlayer = $matintlayer;
      if ( length($matintlayer) > 30 )
      {
        $newmatintlayer =~ s/\w\b// ;
        $newmatintlayer =~ s/\w\b// ;
      }
      $newmatintlayer = "f_" . "$newmatintlayer";
      push ( @materials, $matintlayer );
      push ( @newmaterials, $newmatintlayer );

      my ( $matextlayer, $newmatextlayer );
      unless ( $onlyonelayer == 1 )
      {
        my $extlayer = $constrlines[$#constrlines-1];
        my @row = split ( / : |,/ , $extlayer );
        @row = cleanblanks( @row );
        $matextlayer = $row[3];
        $newmatextlayer = $matextlayer;
        if ( length($matextlayer) > 30 )
        {
          $newmatextlayer =~ s/\w\b// ;
          $newmatextlayer =~ s/\w\b// ;
        }
        $newmatextlayer = "f_" . "$newmatextlayer";
        push ( @materials, $matextlayer );
        push ( @newmaterials, $newmatextlayer );
      }

      $constrlines[0] =~ s/$constrname/$newconstrname/;
      $constrlines[1] = "*itemdoc,fictitious version of construction " . $constrname . " created by Modish script\n";
      $constrlines[2] = "*incat,Modish_fict\n";
      $constrlines[4] =~ s/$matintlayer/$newmatintlayer/;
      unless ( $onlyonelayer == 1 ) { $constrlines[$#constrlines-1] =~ s/$matextlayer/$newmatextlayer/ }
      foreach ( @constrlines )
      {
        push ( @newcopy, $_ );
      }
      @newbigcopy = [ @newcopy ] ;
      $exportconstr{ $newconstrname }{ extlayer } = $newmatextlayer;
      $exportconstr{ $newconstrname }{ intlayer } = $newmatintlayer;
    }

    @materials = uniq( @materials );
    @newmaterials = uniq( @newmaterials );

    my %newmatnums;
    my $countmat = 1;
    foreach ( @newmaterials )
    {
      $newmatnums{$_} = $countmat;
      $countmat++;
    }

    my %matnums;
    $countmat = 1;
    foreach ( @materials )
    {
      $matnums{$_} = $countmat ;
      $countmat++;
    } # --- END NEW CONSTR DATABASE ---
  }

  if ( "oldconstrdb" ~~ @calcprocedures )
  {
     # --- OLD CONSTR DATABASE ---
     my ( @lastbigcopy );
     $countmat = 1;
     foreach my $copyref ( @newbigcopy )
     {
       my ( @lastcopy );
       my @constrlines = @$copyref;
       my $intlayer = $constrlines[3];
       my @row = split ( /\s+|,/ , $intlayer );
       @row = cleanblanks( @row );
       my $matintlayernum = $row[0];
       my $matintlayer = $row[2];

       my $extlayer = $constrlines[$#constrlines];
       my @row = split ( /\s+|,/ , $extlayer );
       @row = cleanblanks( @row );
       my $matextlayernum = $row[0];
       my $matextlayer = $row[2];

       my $newmatnumint = $newmatnums{$matintlayer};
       my $newmatnumext = $newmatnums{$matextlayer};
       $constrlines[3] =~ s/$matintlayernum/$newmatnumint/g;
       $constrlines[$#constrlines] =~ s/$matextlayernum/$newmatnumext/g;
       foreach my $line ( @constrlines )
       {
         push ( @lastcopy, $line );
       }
       push ( @lastbigcopy, [ @lastcopy ] );
     } # --- END OLD CONSTR DATABASE ---
  }
  elsif ( not ( "oldconstrdb" ~~ @calcprocedures ) )
  {
  # --- NEW CONSTR DATABASE ---
    my ( @lastbigcopy );
    $countmat = 1;
    foreach my $copyref ( @newbigcopy )
    {
      my @lastcopy;
      my @constrlines = @$copyref;
      my $onlyonelayer = 0;
      if ( $#constrlines == 5 ) { $onlyonelayer = 1 }

      my $intlayer = $constrlines[4];
      my @row = split ( / : |,/ , $intlayer );
      @row = cleanblanks( @row );
      my $matintlayernum = $row[1];
      my $matintlayer = $row[3];

      unless ( $onlyonelayer == 1 )
      {
        my $extlayer = $constrlines[$#constrlines-1];
        my @row = split ( / : |,/ , $extlayer );
        @row = cleanblanks( @row );
        my $matextlayernum = $row[1];
        my $matextlayer = $row[3];
      }

      my $newmatnumint = $newmatnums{$matintlayer};
      my $newmatnumext = $newmatnums{$matextlayer};
      $constrlines[4] =~ s/$matintlayernum/$newmatnumint/;
      unless ( $onlyonelayer == 1 ) { $constrlines[$#constrlines-1] =~ s/$matextlayernum/$newmatnumext/ }
      foreach my $line ( @constrlines )
      {
        push ( @lastcopy, $line );
      }
      push ( @lastbigcopy, [ @lastcopy ] );
    }
  } # --- END NEW CONSTR DATABASE ---


  foreach ( @lastbigcopy )
  {
    splice ( @updatedlines, $#updatedlines, 0, @$_ );
  }

  open ( CONSTRDBFILE_F, ">$constrdbfile_f" ) or die;
  foreach ( @updatedlines )
  {
    print CONSTRDBFILE_F $_;
  }
  close CONSTRDBFILE_F;

  return ( \@materials, \@newmaterials, \%matnums, \%newmatnums, \%exportconstr );
}

sub compareirrs
{ # THIS COMPARES THE IRRADIANCES TO OBTAIN THE IRRADIANCE RATIOS.
  my ( $zonefilelistsref, $irrsref, $computype, $calcprocedures_ref, $selectives_ref ) = @_;
  my %zonefilelists = %$zonefilelistsref;
  my %irrs = %$irrsref;
  my @calcprocedures = @$calcprocedures_ref;
  my @selectives = @{ $selectives_ref };

  my %irrvars;
  foreach my $zonenum ( sort {$a <=> $b} ( keys %irrs ) )
  {
    my $shdfile = $zonefilelists{ $zonenum }{ shdfile };
    foreach my $monthnum ( sort {$a <=> $b} ( keys %{ $irrs{ $zonenum }{ 1 } } ) )
    {
      foreach my $surfnum ( sort {$a <=> $b} ( keys %{ $irrs{ $zonenum }{ 1 }{ $monthnum } } ) )
      {
        foreach my $hour ( sort {$a <=> $b} ( keys %{ $irrs{ $zonenum }{ 1 }{ $monthnum }{ $surfnum } } ) )
        {

          my $diffsurfirr = $irrs{ $zonenum }{ 1 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr };
          my $whitediffsurfirr;

          if ( scalar( @selectives == 0 ) )
          {
            $whitediffsurfirr = $irrs{ $zonenum }{ 2 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr };
          }
          elsif ( scalar( @selectives > 0 ) )
          {
            $whitediffsurfirr = ( ( $irrs{ $zonenum }{ 2 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr }
              + $irrs{ $zonenum }{ 4 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr } ) / 2 );
          }

          #my $surfirr_amb = $irrs{ $zonenum }{ 3 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr };
          my $diffirrratio;

            if ( $diffsurfirr == 0 ) { $diffsurfirr = 0.0001 }
            if ( $whitediffsurfirr == 0 ) { $whitediffsurfirr = 0.0001 }

            if ( ( $whitediffsurfirr == 0.0001 ) and ( $diffsurfirr == 0.0001 ) )
            {
              $diffirrratio = "1.0000";
            }
            elsif ( $whitediffsurfirr == 0.0001 )
            {
              $diffirrratio = "0.0000";
            }
            elsif ( $diffsurfirr == 0.0001 )
            {
              $diffirrratio = "1.0000";
            }
            else
            {
              $diffirrratio = ( $whitediffsurfirr / $diffsurfirr );
            }


            #say REPORT"\$zonenum $zonenum \$monthnum $monthnum \$surfnum $surfnum \$hour $hour";
            #say REPORT "\$diffsurfirr: $diffsurfirr"; say REPORT "\$whitediffsurfirr: $whitediffsurfirr"; say REPORT "\$diffirrratio: $diffirrratio\n";



          $irrvars{ $zonenum }{ $monthnum }{ $surfnum }{ $hour }{ diffirrvar } = $diffirrratio;


          my $dirsurfirr;

          unless ( "coexistent" ~~ @calcprocedures )
          {
            $dirsurfirr = $irrs{ $zonenum }{ 1 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr };
          }
          else
          {
            $dirsurfirr = $irrs{ $zonenum }{ 3 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr };
          }

          my $dirwhitesurfirr;
          if ( scalar( @selectives == 0 ) )
          {
            $dirwhitesurfirr = $irrs{ $zonenum }{ 2 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr };
          }
          elsif ( scalar( @selectives > 0 ) )
          {
            $dirwhitesurfirr = ( ( $irrs{ $zonenum }{ 2 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr }
              + $irrs{ $zonenum }{ 4 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr } ) / 2 );
          }

          my $dirirrratio;

          if ( $dirsurfirr == 0 ) { $dirsurfirr = 0.0001 }
          if ( $dirwhitesurfirr == 0 ) { $dirwhitesurfirr = 0.0001 }

          if ( ( $dirwhitesurfirr == 0.0001 ) and ( $dirsurfirr == 0.0001 ) )
          {
            $dirirrratio = "1.0000";
          }
          elsif ( $dirwhitesurfirr == 0.0001 )
          {
            $dirirrratio = "0.0000";
          }
          elsif ( $dirsurfirr == 0.0001 )
          {
            $dirirrratio = "1.0000";
          }
          else
          {
            $dirirrratio = ( $dirwhitesurfirr / $dirsurfirr );
          }
          #say REPORT"\$zonenum $zonenum \$monthnum $monthnum \$surfnum $surfnum \$hour $hour";
          #say REPORT "\$dirsurfirr: $dirsurfirr, \$whitedirsurfirr: $whitedirsurfirr, \$dirirrratio: $dirirrratio\n";

          $irrvars{ $zonenum }{ $monthnum }{ $surfnum }{ $hour }{ dirirrvar } = $dirirrratio;

          $countsurf++;
        }
        $counthour++;
      }
      $countmonth++;
    }
    $countzone++;
  }
  return ( \%irrvars );
}

sub fillhours
{ # THIS COMPLETES THE FILLING OF THE DATA STRUCTURES GIVING INFORMATION ABOUT THE DAYLIT HOURS.
  my ( $newhourvalsref, $monthname, $daylighthoursref ) = @_;
  my @hadhours = @$newhourvalsref;
  my %lithours = %$daylighthoursref;
  my @monthhours = @{ $lithours{ $monthname } };
  my @values;

  my $sunhoursnum = 0;
  foreach my $lightcond ( @monthhours )
  {
    unless ( $lightcond == 1 )
    {
      $sunhoursnum++;
    }
  }

  if ( $sunhoursnum == scalar( @hadhours ) )
  {
    my $counthr = 1;
    my $countlit = 0;
    foreach my $lightcond ( @monthhours )
    {
      if ( $lightcond == 1 )
      {
        push ( @values, "1.0000" );
      }
      else
      {
        push ( @values, $hadhours[ $countlit ] );
        $countlit++;
      }
      $counthr++;
    }
    return ( @values );
  }
}


sub modifyshda
{ # THIS MODIFIES THE ".shda" FILE ON THE BASIS OF THE IRRADIANCE RATIOS.
  my ( $comparedirrsref, $surfslistref, $zonefilelistsref, $shdfileslistref, $daylighthoursref, $irrvarsref, $threshold, $tempmod, $tempreport, $tempmoddir, $tempreportdir, $elm, $radtype, $calcprocedures_ref, $irrs_ref, $conffile_f2, $shdfile ) = @_; ##### CONDITION! "diffuse" AND "direct".
  my %surfslist = %$surfslistref;
  my %zonefilelists = %$zonefilelistsref;
  my %shdfileslist = %$shdfileslistref;
  my %daylighthours = %$daylighthoursref;
  my %irrvars = %$irrvarsref;
  my @calcprocedures = @$calcprocedures_ref;
  my %irrs = %{ $irrs_ref };
  my ( @printcontainer, @monthnames, @pushmodline, @pushreportline, @mainbag, @mainoriginal );

  foreach my $zonenum ( sort {$a <=> $b} ( keys %irrvars ) )
  {
    my $inlinesref = $shdfileslist{ $zonenum };
    my @inlines = @$inlinesref;

    my ( @surfsdo, @vehicle );
    my ( @insertlines, @inserts, $bringline ) ;
    if ( ( "keepdirshdf" ~~ @calcprocedures ) and ( $radtype eq "direct" ) )
    {
      foreach my $monthnum ( sort {$a <=> $b} ( keys %{ $irrvars{ $zonenum } } ) )
      {
        foreach my $surfnum ( sort {$a <=> $b} ( keys %{ $irrvars{ $zonenum }{ $monthnum } } ) )
        {
          my $surfname = $surfslist{$zonenum}{$surfnum}{surfname};
          push( @surfsdo, $surfname );
        }
        last;
      }

      my ( $insertlines_ref ) = readshdfile( $shdfile, \@calcprocedures, $conffile_f2, "go" );
      @insertlines = @{ $insertlines_ref }; #say "INSERTLINES: " . dump( @insertlines );

      foreach my $lin ( @insertlines )
      {
        my $signal = "no";
        foreach my $surf ( @surfsdo )
        {
          if ( $lin =~ /$surf/ )
          {
            $signal = "yes";
          }
        }

        if ( ( $lin =~ /direct - surface / ) and ( $signal eq "yes" ) )
        {
          my @elts = split( "#", $lin );
          $bringline = $elts[0];
          push ( @vehicle, $bringline );
        }
      }
    }

    my $semaphore = 0;
    my ( $readmonthname, $readmonthnum );
    my $c = 0;
    foreach my $line ( @inlines )
    {
      my $line2;
      my @row = split( /\s+|,/ , $line);
      my ( $readsurfname, $readsurfnum );

      if ( ( $row[0] eq "*" ) and ( $row[1] eq "month:" ) )
      {
        $semaphore = 1;
        $readmonthname = $row[2];
        $readmonthname =~ s/`//g;
        $readmonthnum = getmonthnum( $readmonthname );
      }

      if ( ( ( $row[0] eq "24" ) and ( $row[1] eq "hour" ) and ( $row[1] eq "surface" ) ) or ( $row[0] eq "*end" ) )
      {
        $semaphore = 0;
      }

      my ( @newhourvals, @newhourvals2, @were );
      foreach my $monthnum ( sort {$a <=> $b} ( keys %{ $irrvars{ $zonenum } } ) )
      {
        my $monthname = getmonthname( $monthnum );
        push ( @monthnames, $monthname );
        push ( @monthnames2, $monthname );

        foreach my $surfnum ( sort {$a <=> $b} ( keys %{ $irrvars{ $zonenum }{ $monthnum } } ) )
        {
          my $surfname = $surfslist{$zonenum}{$surfnum}{surfname};
          push( @surfsdo, $surfname );
          foreach my $hour ( sort {$a <=> $b} ( keys %{ $irrvars{ $zonenum }{ $monthnum }{ $surfnum } } ) )
          {
            my $irrvariation;
            #say REPORT"\$radtype $radtype \$zonenum $zonenum \$monthnum $monthnum \$surfnum $surfnum \$hour $hour";
            if ( $radtype eq "diffuse" )
            {
              $irrvariation = $irrvars{ $zonenum }{ $monthnum }{ $surfnum }{ $hour }{ diffirrvar };
              #say REPORT "DIFF VARIATION: $irrvariation";
            }
            elsif ( $radtype eq "direct" )
            {
              $irrvariation = $irrvars{ $zonenum }{ $monthnum }{ $surfnum }{ $hour }{ dirirrvar };
              #say REPORT "DIR VARIATION: $irrvariation";
            }

            if ( ( $zonenum ) and ( $monthnum ) and ( $hour ) and ( $surfnum ) and ( $irrvariation ) and ( $surfname ) and ( $monthname ) and ( $surfnum eq $elm ) ) # I.E. IF ALL THE NEEDED DATA EXIST
            {
              if ( $semaphore == 1 )
              {
                if ( $row[27] eq "surface" )
                {
                  $readsurfname = $row[28];
                  $readsurfnum = $surfslist{$zonenum}{$readsurfname}{surfnum} ;
                  my @filledhourvals;
                  if ( ( $row[25] eq $radtype ) and ( $readsurfname eq $surfname ) )
                  {
                    my @hourvals = ( @row[ 0..23 ] );
                    my $counthour = 1;
                    foreach my $el ( @hourvals )
                    { # $el IS THE SHADING FACTOR IN THE ORIGINAL SHDA FILE.
                      # %irrvariation IS THE IRRADIANCE DIFFERENCE BETWEEN THE "WHITE" MODEL AND THE "BLACK" MODEL.
                      # $ambase IS THE AMBIENT RADIATION WITHOUT SHADINGS

                      my ( $calcamount, $improvedguess, $newshadingvalue );


                      if ( $radtype eq "diffuse" )
                      {
                        unless ( ( ( "coexistent" ~~ @calcprocedures ) and ( "radical" ~~ @calcprocedures ) )
                          or ( "radical" ~~ @calcprocedures ) )
                        {
                          $calcamount = ( 1 - $el ); # THIS IS THE RATIO OF NON-SHADED IRRADIATION AS CALCULATED BY THE ESP-r's ISH MODULE
                          $improvedguess = ( $calcamount * $irrvariation ); # THIS IS THE RATIO ABOVE CORRECTED BY MULTIPLYING IT BY THE IRRADIANCE RATIO TO TAKE REFLECTIONS INTO ACCOUNT.
                          $newshadingvalue = ( 1 - $improvedguess ); # AS THE NAME SAYS, THIS IS THE NEW SHADING FACTOR.
                          if ( $newshadingvalue > $el ) { $newshadingvalue = $el }; # IF THE SHADING VALUE IS INCREASING, KEEP THE OLD ONE.
                          #say REPORT "IN MOD CD DIFF \$radtype: $radtype, \$calcamount: $calcamount, \$improvedguess: $improvedguess, \$newshadingvalue: $newshadingvalue.";
                          #say REPORT "\$radtype: $radtype, \$calcamount: $calcamount, \$improvedguess: $improvedguess, \$newshadingvalue: $newshadingvalue.";
                        }
                      }
                      elsif ( $radtype eq "direct" )
                      {
                        unless ( ( ( not( "coexistent" ~~ @calcprocedures ) ) and ( "radical" ~~ @calcprocedures ) ) )
                        {
                          $calcamount = ( 1 - $el ); # THIS IS THE RATIO OF NON-SHADED IRRADIATION AS CALCULATED BY THE ESP-r's ISH MODULE
                          $improvedguess = ( $calcamount * $irrvariation ); # THIS IS THE RATIO ABOVE CORRECTED BY MULTIPLYING IT BY THE IRRADIANCE RATIO TO TAKE REFLECTIONS INTO ACCOUNT.
                          $newshadingvalue = ( 1 - $improvedguess ); # AS THE NAME SAYS, THIS IS THE NEW SHADING FACTOR.
                          if ( $newshadingvalue > $el ) { $newshadingvalue = $el }; # IF THE SHADING VALUE IS INCREASING, KEEP THE OLD ONE.
                          #say REPORT "IN MOD CD DIR \$radtype: $radtype, \$calcamount: $calcamount, \$improvedguess: $improvedguess, \$newshadingvalue: $newshadingvalue.";
                          #say REPORT "\$radtype: $radtype, \$calcamount: $calcamount, \$improvedguess: $improvedguess, \$newshadingvalue: $newshadingvalue.";
                        }
                      }


                      if ( "radical" ~~ @calcprocedures )
                      {
                        if ( $radtype eq "diffuse" )
                        {
                          my $noshads = $irrs{ $zonenum }{ 1 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr };
                          my $withshads = $irrs{ $zonenum }{ 2 }{ $monthnum }{ $surfnum }{ $hour }{ meandiffirr };
                          unless ( $noshads == 0 )
                          {

                            if ( ( $withshads == 0.0001 ) and ( $noshads == 0.0001 ) )
                            {
                              $newshadingvalue = "1.0000";
                            }
                            elsif ( $withshads == 0.0001 )
                            {
                              $newshadingvalue = "0.0000";
                            }
                            elsif ( $noshads == 0.0001 )
                            {
                              $newshadingvalue = "1.0000";
                            }
                            else
                            {
                              $newshadingvalue = ( 1 - ( $withshads / $noshads ) );
                            }
                          }
                          else
                          {
                            $newshadingvalue = "1.0000";
                          }
                          #say REPORT "IN MOD RADICAL DIFF \$radtype: $radtype, \$noshads: $noshads, \$withshads: $withshads, \$newshadingvalue2: $newshadingvalue2";
                        }
                        elsif ( ( $radtype eq "direct" ) and ( not( "coexistent" ~~ @calcprocedures ) ) )
                        {
                          if ( "keepdirshdf" ~~ @calcprocedures )
                          {
                            $newshadingvalue = $el;
                          }
                          else
                          {
                            my $nodirshads = $irrs{ $zonenum }{ 1 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr };
                            my $withdirshads = $irrs{ $zonenum }{ 2 }{ $monthnum }{ $surfnum }{ $hour }{ meandirirr };
                            unless ( $nodirshads == 0 )
                            {

                              if ( ( $withdirshads == 0.0001 ) and ( $nodirshads == 0.0001 ) )
                              {
                                $newshadingvalue = "1.0000";
                              }
                              elsif ( $withdirshads == 0.0001 )
                              {
                                $newshadingvalue = "0.0000";
                              }
                              elsif ( $nodirshads == 0.0001 )
                              {
                                $newshadingvalue = "1.0000";
                              }
                              else
                              {
                                $newshadingvalue = ( 1 - ( $withdirshads / $nodirshads ) );
                              }
                            }
                            else
                            {
                              $newshadingvalue = "1.0000";
                            }
                            #say REPORT "IN MOD RADICAL DIR \$radtype: $radtype, \$nodirshads: $nodirshads, \$withdirshads: $withdirshads, \$newshadingvalue2: $newshadingvalue2."
                          }
                        }
                      }


                      if ( ( $counthour == $hour ) and ( $readmonthname eq $monthname ) and ( $readsurfnum == $surfnum ) ) # I.E.: IF THIS LINE IS THE RIGHT ONE...
                      {

                        if ( $newshadingvalue == 1 )
                        {
                          $newshadingvalue = "1.0000";
                        }

                        if ( $newshadingvalue == 0 )
                        {
                          $newshadingvalue = "0.0000";
                        }

                        if ( $newshadingvalue > 1 )
                        {
                          $newshadingvalue = "1.0000";
                        }

                        if ( ( $newshadingvalue > 0 ) and ( $newshadingvalue < 1 ) ) # IF THE VARIATION OF IRRADIANCE FROM MODEL A AND MODEL B IS NEGATIVE...
                        {  # ...INCREASE THE SHADING FACTOR ACCORDINGLY.
                          $newshadingvalue = sprintf ( "%.4f", $newshadingvalue ); # FORMAT THE NUMBER SO THAT IT HAS FOUR DECIMALS
                        }

                        if ( ( $newshadingvalue > -10 ) and ( $newshadingvalue < 0 ) )
                        {
                          $newshadingvalue = sprintf ( "%.3f", $newshadingvalue ); # IF THE NUMBER IS COMPRISED BETWEEN -10 AND = 0 FORMAT IT SO THAT IT HAS 3 DECIMALS
                        }

                        if ( ( $newshadingvalue > -100 ) and ( $newshadingvalue <= -10 ) )
                        {
                          $newshadingvalue = sprintf ( "%.2f", $newshadingvalue ); # IT THE NUMBER IS COMPRISED BETWEEN -100 AND = -10 FORMAT IT SO THAT IT HAS 2 DECIMALS
                        }

                        my $irrvariation = sprintf ( "%.4f", $irrvariation );
                        print REPORT "OBTAINED \$radtype,$radtype, \$zonenum $zonenum, \$monthnum $monthnum, \$surfnum $surfnum, \$hour $hour, OLD SHADING VALUE $el, NEW SHADING VALUE $newshadingvalue, \$irrvariation $irrvariation\n";

                        push ( @newhourvals, $newshadingvalue );
                        push ( @newhourvals2, $irrvariation );
                        push ( @were, $el );
                      }
                      $counthour++;
                    }
                    #say REPORT "IN MODIFISHDA \@newhourvals (NEW SHADING VALUES) " . dump(@newhourvals) . ", in monthnum: $monthnum, monthname: $monthname, surfnum: $surfnum, hour: $hour";
                    #say REPORT "IN MODIFISHDA \@newhourvals2 (IRRADIANCE RATIOS)" . dump(@newhourvals2) . ", in monthnum: $monthnum, monthname: $monthname, surfnum: $surfnum, hour: $hour";
                    #say REPORT "IN MODIFISHDA OLD SHADING VALUES WERE : " . dump(@were) . ", in monthnum: $monthnum, monthname: $monthname, surfnum: $surfnum, hour: $hour";

                    my @filledhourvals = fillhours( \@newhourvals, $monthname, \%daylighthours );

                    my @filledhourvals2 = fillhours( \@newhourvals2, $monthname, \%daylighthours );

                    if ( ( scalar ( @filledhourvals ) == 24 ) and ( $monthname eq $monthnames[0] ) )
                    {
                      shift @monthnames;
                      my @firstarr = @filledhourvals[ 0..11 ];
                      my @secondarr = @filledhourvals[ 12..$#filledhourvals ];
                      my $joinedfirst = join ( ' ' , @firstarr );
                      my $joinedsecond = join ( ' ' , @secondarr );
                      if ( $radtype eq "diffuse" )
                      {
                        my $newline = "$joinedfirst " . "$joinedsecond" . " # diffuse - surface " . "$readsurfname $monthname\n";
                        print TEMPMOD $newline;
                      }
                      elsif ( $radtype eq "direct" )#
                      {
                        my $newline;
                        if ( "keepdirshdf" ~~ @calcprocedures )
                        {
                          my $bringline = shift( @vehicle );
                          $newline = "$bringline" . "# direct - surface " . "$readsurfname $monthname\n";
                        }
                        else
                        {
                          $newline = "$joinedfirst " . "$joinedsecond" . " # direct - surface " . "$readsurfname $monthname\n";
                        }
                        print TEMPMODDIR $newline;
                      }
                    }

                    if ( ( scalar ( @filledhourvals2 ) == 24 ) and ( $monthname eq $monthnames2[0] ) )
                    {
                      shift @monthnames2;
                      my @firstarr2 = @filledhourvals2[ 0..11 ];
                      my @secondarr2 = @filledhourvals2[ 12..$#filledhourvals2 ];
                      my $joinedfirst2 = join ( ' ' , @firstarr2 );
                      my $joinedsecond2 = join ( ' ' , @secondarr2 );
                      if ( $radtype eq "diffuse" )
                      {
                        my $newline2 = "$joinedfirst2 " . "$joinedsecond2" . " # diffuse for surface " . "$readsurfname in $monthname\n";
                        print TEMPREPORT $newline2;
                      }
                      elsif ( $radtype eq "direct" )#
                      {
                        my $newline2 = "$joinedfirst2 " . "$joinedsecond2" . " # direct for surface " . "$readsurfname in $monthname\n";
                        print TEMPREPORTDIR $newline2;#
                      }
                    }
                  }
                }
              }
            }
            $countref++;
          }
        }
      }
      $c++;
    }
  }
}


sub getbasevectors
{ # THIS GETS THE PRE-COMPUTED EVENLY DISTRIBUTED N POINTS ON THE SURFACE OF A HEMISPHERE.
  my ( $dirvectorsnum ) = @_;
  my @basevectors;
  if ( $dirvectorsnum == 1 )
  {
    @basevectors = ( #[ 0, 0, 0 ] , # origin, base point of direction vector
        [ 0, 0, 1 ], # direction vector of high, central, vertical point
         ); # lowest vertices
  }
  if ( $dirvectorsnum == 5 )
  {
    @basevectors = ( #[ 0, 0, 0 ] , # origin, base point of direction vector
        [ 0, 0, 1 ], # direction vector of high, central, vertical point
        [ 0.7071, -0.7071, 0.4472 ] , [ 0.7071, 0.7071, 0.4472 ], [ -0.7071, 0.7071, 0.4472 ], [ -0.7071, -0.7071, 0.4472 ] ); # lowest vertices
  }
  elsif ( $dirvectorsnum == 17 )
  {
    @basevectors = ( #[ 0, 0, 0 ] , # origin, base point of direction vector
    [ 0, 0, 1 ], # direction vector of high, central, vertical point
    [ 0.1624, 0.4999, 0.8506 ], [ 0.1624, -0.4999, 0.8506 ], [ -0.2628, 0.8090, 0.5257 ], [ -0.2628, -0.8090, 0.5257 ],
    [ 0.2763, 0.8506, 0.4472 ], [ 0.2763, -0.8506, 0.4472 ], [ -0.4253, 0.3090, 0.8506 ], [ -0.4253, -0.3090, 0.8506,  ],
    [ 0.5257, 0.0, 0.8506 ], [ 0.5877, 0.8090, 0.0 ], [ 0.6881, 0.4999, 0.5257 ], [ 0.6881, -0.4999, 0.5257 ],
    [ -0.7236, 0.5257, 0.4472 ], [ -0.7236, -0.5257, 0.4472 ], [ -0.8506, 0.0, 0.5257 ], [ 0.8944, 0.0, 0.4472 ]
    ); # lowest vertices
  }
  return ( @basevectors );
}

sub createfictgeofile
{  # THIS MANAGES THE MODIFICATION OF THE FICTITIOUS GEO FILES FOR THE ZONE BY ADJUSTING THE OBSTRUCTION CONSTRUCTIONS TO FICTITIOUS EQUIVALENTS
  my ( $geofile, $obsconstrsetref, $geofile_f, $paths_ref, $calcprocedures_ref, $conffile_f1 ) = @_;
  my @obsconstrset = @$obsconstrsetref;
  my %paths = %{ $paths_ref }; #say "\%paths: " . dump( \%paths );
  my @calcprocedures = @{ $calcprocedures_ref }; #say "\@calcprocedures: " . dump( \@calcprocedures );
  open ( GEOFILE, "$geofile" ) or die;
  my @lines = <GEOFILE>;
  close GEOFILE;

  open ( GEOFILE_F, ">$geofile_f" ) or die;

  foreach my $line ( @lines )
  {
    if ( $line =~ /^\*obs/ )
    {
      foreach my $obsconstr ( @obsconstrset )
      {
        my $newobsconstr = $obsconstr;
        if ( length($obsconstr) > 30 )
        {
          $newobsconstr =~ s/\w\b// ;
          $newobsconstr =~ s/\w\b// ;
        }
        $newobsconstr = "f_" . $newobsconstr;
        $line =~ s/$obsconstr/$newobsconstr/;
      }
      print GEOFILE_F $line;
    }
    else
    {
      print GEOFILE_F $line;
    }
  }
  close ( GEOFILE_F);
  ###### PART ADDED IN VERSION 1.09
  my ( $shortgeofile_f, $shortgeofile_f1, $geofile_f1 );
  if ( "radical" ~~ @calcprocedures )
  {
    $geofile_f1 = $geofile_f;
    $geofile_f1 =~ s/\.geo$// ;
    $geofile_f1 = $geofile_f1 . "1";
    $geofile_f1 = $geofile_f1 . ".geo";
    #say "GEOFILE_F1: $geofile_f1";
    #say "GEOFILE_F: $geofile_f";

    open( GEOFILE_F, "$geofile_f" ) or die;
    my @lines_f = <GEOFILE_F>;
    close GEOFILE_F;

    open( GEOFILE_F1, ">$geofile_f1" );
    foreach my $line_f ( @lines_f )
    {
      if ( $line_f =~ /^\*obs/ )
      {
        $line =~ s/(\s+)$// ;
        my @elts = split( ",", $line_f );
        $elts[4] = 0.01;
        $elts[5] = 0.01;
        $elts[6] = 0.01;
        $line_f = "$elts[0],$elts[1],$elts[2],$elts[3],$elts[4],$elts[5],$elts[6],$elts[7],$elts[8],$elts[9],$elts[10]";
        $line_f =~ s/^,//;
        $line_f =~ s/,$//;
      }
      print GEOFILE_F1 $line_f;
    }
    close GEOFILE_F1;
    my $zonepath = $paths{zonepath};
    $shortgeofile_f1 = $geofile_f1;
    $shortgeofile_f1 =~ s/^$zonepath// ;
    $shortgeofile_f1 =~ s/^\/// ;
    $shortgeofile_f = $geofile_f;
    $shortgeofile_f =~ s/^$zonepath// ;
    $shortgeofile_f =~ s/^\/// ;

    #say "\$shortgeofile_f1: " . dump( $shortgeofile_f1 );
    #say "\$shortgeofile_f: " . dump( $shortgeofile_f );
    #say "\$geofile_f1: " . dump( $geofile_f1 );
    #say "\$geofile_f: " . dump( $geofile_f );
    open ( CONFFILE_F1, "$conffile_f1" ) or die;
    my @lines_old = <CONFFILE_F1>;
    close CONFFILE_F1;

    my $conffile_f1_old = $conffile_f1 . ".old"; #say "\$conffile_f1_old: $conffile_f1_old" ;
    `mv -f $conffile_f1 $conffile_f1_old`;
    say REPORT "mv -f $conffile_f1 $conffile_f1_old";

    open( CONFFILE_F1, ">$conffile_f1" ) or die;
    foreach my $line ( @lines_old )
    {
      if ( $line =~ /^\*geo/ )
      {
        $line =~ s/$shortgeofile_f/$shortgeofile_f1/ ; #say "changed $shortgeofile_f in $shortgeofile_f1 .";
      }
      print CONFFILE_F1 $line;
    }

    if ( "radical1" ~~ @calcprocedures )
    {
`prj -file $conffile_f1 -mode script<<YYY
b
m
c
f
*
b
a
-
-
-
-
-
-
YYY
`;

say REPORT "prj -file $conffile_f1 -mode script<<YYY
b
m
c
f
*
b
a
-
-
-
-
-
-
YYY
";
    }
  }
  ###### END PART ADDED IN VERSION 1.09
}


sub creatematdbfiles
{ # THIS MANAGES THE CREATION OF THE TWO FICTITIOUS MATERIALS DATABASES:
  # ONE FOR THE THE "UNREFLECTIVE" MODEL AND THE OTHER FOR THE "REFLECTIVE" ONE.
  my ( $matdbfile,  $matdbfile_f1, $matdbfile_f2, $calcprocedures_ref, $constrdbfile_f, $obsdata_ref ) = @_;

  my @calcprocedures = @{ $calcprocedures_ref };
  my @obsdata = @{ $obsdata_ref }; #say "OBSDATA: " .dump( @obsdata );

  my ( @box );
  my ( %exportrefl, %obslayers );

  #say "MATDBFILE: " .  $matdbfile ;
  open ( MATDBFILE, "$matdbfile" ) or die;
  my @matlines = <MATDBFILE>;
  close MATDBFILE;

  open( CONSTRDBFILE_F, "$constrdbfile_f" ) or die;
  my @constrlines = <CONSTRDBFILE_F>;
  close CONSTRDBFILE_F;

  my @obs = uniq( map{ $_->[1] } @obsdata );

  my $count = 0;
  foreach my $ob ( @obs )
  {
    my $semaphore = "off";
    foreach my $constrline ( @constrlines )
    {
      my @row = split( ",", $constrline );
      if ( $row[0] eq "*item" )
      {
        if ( $row[1] eq $ob )
        {
          $semaphore = "on";
        }
      }
      if ( $row[0] eq "*end_item" )
      {
        $semaphore = "off";
        $count++;
      }
      if ( $semaphore eq "on" )
      {
        if ( $row[0] eq "*layer" )
        {
          my @els = split( "\s+|,|:", $row[3] );
          $els[0] =~ s/(\s+)// ;
          push( @{ $obslayers{$ob} }, $els[0] );
        }
      }
    }
  }

  my @obsmats;
  foreach my $obkey ( keys %obslayers )
  {
    my @ob = @{ $obslayers{$obkey} };
    push( @obsmats, $ob[0], $ob[-1] );
  }
  @obsmats = uniq( @obsmats );

  my ( @bag, @row, @firstloop, @secondloop );
  my $semaphore = "off";
  foreach my $matline ( @matlines )
  {
    chomp $matline;
    $matline =~ s/\s+// ;
    my @row = split( ",", $matline );
    if ( $row[0] eq "*item" )
    {
      if ( $row[1] ~~ @obsmats )
      {
        $semaphore = "on";
      }
      else
      {
        $semaphore = "off";
      }
    }

    my @e = split( ",", $matline );
    if ( ( $e[0] =~ /^\d/ ) and ( $e[-1] =~ /\D$/ ) )
    {
      if ( $semaphore eq "off" )
      {
        if ( not( "diluted" ~~ @calcprocedures ) )
        {
          my $lin = "$e[0],$e[1],$e[2],$e[3],$e[4],0.999,0.999,$e[7],$e[8],$e[9]";
          push( @firstloop, $lin );
          my $linn = "$e[0],$e[1],$e[2],$e[3],$e[4],0.999,0.999,$e[7],$e[8],$e[9]";
          push( @secondloop, $linn );
        }
        elsif ( "diluted" ~~ @calcprocedures )
        {
          my $lin = "$e[0],$e[1],$e[2],$e[3],$e[4],$e[5],$e[6],$e[7],$e[8],$e[9]";
          push( @firstloop, $lin );
          my $linn = "$e[0],$e[1],$e[2],$e[3],$e[4],$e[5],$e[6],$e[7],$e[8],$e[9]";
          push( @secondloop, $linn );
        }
      }
      elsif ( $semaphore eq "on" )
      {
        $exportrefl{ $row[1] }{ absout } =  $e[5];
        $exportrefl{ $row[1] }{ absin } = $e[6];
        if ( not( "diluted" ~~ @calcprocedures ) )
        {
          my $lin = "$e[0],$e[1],$e[2],$e[3],$e[4],0.999,0.999,$e[7],$e[8],$e[9]";
          push( @firstloop, $lin );
          my $linn = "$e[0],$e[1],$e[2],$e[3],$e[4],$e[5],$e[6],$e[7],$e[8],$e[9]";
          push( @secondloop, $linn );
        }
        elsif ( "diluted" ~~ @calcprocedures )
        {
          my $lin = "$e[0],$e[1],$e[2],$e[3],$e[4],0.999,0.999,$e[7],$e[8],$e[9]";
          push( @firstloop, $lin );
          my $linn = "$e[0],$e[1],$e[2],$e[3],$e[4],$e[5],$e[6],$e[7],$e[8],$e[9]";
          push( @secondloop, $linn );
        }
      }
    }
    else
    {
      push( @firstloop, $matline );
      push( @secondloop, $matline );
    }
  }

  open( my $MATDBFILE_F1, ">$matdbfile_f1" ) or die;
  foreach my $line ( @firstloop )
  {
    say $MATDBFILE_F1 $line ;
  }
  close $MATDBFILE_F1;

  open( my $MATDBFILE_F2, ">$matdbfile_f2" ) or die;
  foreach my $line ( @secondloop )
  {
    say $MATDBFILE_F2 $line ;
  }
  close $MATDBFILE_F2;

  return ( \%exportrefl, \%obslayers );
}


sub adjust_radmatfile
{ # THIS CHECKS IF THE RADIANCE MATERIALS FILE HAS BEEN PROPERLY MODIFIED.
# IF NOT, THIS DOES THE MODIFICATION. THIS IS USED WHEN THE SPECULAR FRACTIONS
# AND ROUGHNESSES IN THE MATERIALS DATABASE ARE ALL SET TO 0.
  my ( $exportconstrref, $exportreflref, $conffile, $path, $specularratios_ref,
  $obslayers_ref, $selectives_ref, $paths_ref ) = @_;
  my %exportconstr = %$exportconstrref;
  my %exportrefl = %$exportreflref;
  my %paths = %{ $paths_ref };

  my @specularratios = @$specularratios_ref;
  my %obslayers = %{ $obslayers_ref };
  my @selectives = @{ $selectives_ref };

  my %hs;
  foreach $el ( @specularratios )
  {
    my @row = split( ":", $el );
    {
      $row[0] = "_" . $row[0] ;
      $hs{$row[0]}{spec} = $row[1];
      $hs{$row[0]}{roughn} = $row[2];
    }
  }

  my $radmat_f2 = $conffile;
  if ( $radmat_f2 =~ /$path\/cfg\// )
  {
    $radmat_f2 =~ s/$path\/cfg\///;
  }
  else
  {
    $radmat_f2 =~ s/$path\///;
  }

  my $radpath = $paths{radpath}; #say "\$radpath " . $radpath ; ###.
  $radmat_f2 =~ s/.cfg//;
  $radmat_f2 = $radmat_f2 . "_Extern.mat";
  #$radmat_f2 = "$path/rad/$radmat_f2"; #OVERRIDDEN below
  $radmat_f2 = "$radpath/$radmat_f2"; #say "\$radmat_f2 " . $radmat_f2 ; ###.

  my $radmattemp = $radmat_f2 . ".temp"; #say "\$radmattemp " . $radmattemp ; ###.
  `mv -f $radmat_f2 $radmattemp`;
  open( RADMATTEMP, "$radmattemp" ) or die;
  my @lines = <RADMATTEMP>;
  close RADMATTEMP;
  open( RADMAT_F2, ">$radmat_f2" ) or die;
  my $count = 0;
  my @constrs = keys %exportconstr;
  foreach ( @lines )
  {
    my ( $spec, $roughn );
    my $lin = $lines[ $count + 4 ];
    my @arr = split( /\s+/, $lin );
    if ( ( $_ =~ /^#/ ) and ( $_ =~ /ternal MLC Colours.../ ) )
    {
      my $description = $lines[ $count + 1 ] ;

      foreach my $const ( keys %hs )
      {
        if ( $description =~ /$const/ )
        {
          $spec = $hs{$const}{spec};
          $roughn = $hs{$const}{roughn};
          $lines[ $count + 4 ] = "5  $arr[1] $arr[2] $arr[3] $spec $roughn \n";
          last;
        }
      }
    }
    print RADMAT_F2 $lines[ $count ];
    $count++;
  }
  close RADMAT_F2;

  if ( scalar( @selectives ) > 0 )
  {
    my $radmat_f3 = $conffile;
    if (  $radmat_f3 =~ /$path\/cfg\// )
    {
      $radmat_f3 =~ s/$path\/cfg\///;
    }
    else
    {
      $radmat_f3 =~ s/$path\///;
    }

    $radmat_f3 =~ s/.cfg//;
    $radmat_f3 = $radmat_f3 . "_Extern.mat";
    $radmat_f3 = "$path/rad/$radmat_f3";
    my $radmattemp3 = $radmat_f3 . ".temp";
    `mv -f $radmat_f3 $radmattemp3`;
    open( RADMATTEMP3, "$radmattemp3" ) or die;
    my @lines = <RADMATTEMP3>;
    close RADMATTEMP3;
    open( RADMAT_F3, ">$radmat_f3" ) or die;
    my $count = 0;
    my @constrs = keys %exportconstr;
    foreach ( @lines )
    {
      my ( $spec, $roughn );
      my $lin = $lines[ $count + 4 ];
      my @arr = split( /\s+/, $lin );
      if ( ( $_ =~ /^#/ ) and ( $_ =~ /ternal MLC Colours.../ ) )
      {
        my $description = $lines[ $count + 1 ] ;

        foreach my $const ( keys %hs )
        {
          if ( $description =~ /$const/ )
          {
            $spec = $hs{$const}{spec};
            $roughn = $hs{$const}{roughn};
            $lines[ $count + 4 ] = "5  $arr[1] $arr[2] $arr[3] $spec $roughn \n";
            last;
          }
        }
      }
      print RADMAT_F3 $lines[ $count ];
      $count++;
    }
    close RADMAT_F3;
  }
}


sub calcdirvectors
{ # THIS CALCULATES THE NEEDED DIRECTION VECTORS AT EACH GRID POINT.
  my @winscoords = @_;
  my ( @groupbag );
  foreach my $surf ( @winscoords )
  {
    my ( @surfbag );
    foreach my $v ( @$surf )
    {
      my @fields = @$v;
      my $coordsref = $fields[0];
      my @coords = @$coordsref;
      my $vertex = Vector::Object3D::Point->new( x => $coords[0], y => $coords[1], z => $coords[2], );
      push ( @surfbag, $vertex);
    }

    my $polygon = Vector::Object3D::Polygon->new(vertices => [ @surfbag ]);
    my $normal_vector = $polygon->get_normal_vector;
    my ($x_, $y_, $z_) = $normal_vector->array;
    my ( $x, $y, $z ) = ( -$x_, -$y_, -$z_ );
    my $max = max( abs($x), abs($y), abs($z) );
    my @dirvector;
    unless ( $max == 0 )
    {
      $x = ( $x / $max );
      $y = ( $y / $max );
      $z = ( $z / $max );
      @dirvector =  ( $x, $y, $z );
    }
    else
    {
      @dirvector =  ( $x, $y, $z );
    }
    push ( @groupbag, [ @dirvector ] )
  }
  return ( @groupbag);
}

sub prunepoints
{ # IT PRUNES AWAY THE GRID POINTS FALLING OUTSIDE THE SURFACE.
  my ( $gridpoints_transitionalref, $xyzcoordsref ) = @_;
  my @gridpoints = @$gridpoints_transitionalref;
  my @vertstaken = @$xyzcoordsref;

  my @verts;
  foreach ( @vertstaken )
  {
    my @fields = @$_;
    my @xs = @{ $fields[0] };
    my @ys = @{ $fields[1] };
    my @zs = @{ $fields[2] };
    my $i = 0;
    my @bag;
    foreach ( @xs )
    {
      push ( @bag, [ $xs[ $i ], $ys[ $i ], $zs[ $i ] ] );
      $i++;
    }
    push ( @verts, [ @bag ] );
  }

  my ( @coords, @prunegridpoints );

  foreach my $v ( @gridpoints )
  {
    my @fields = @$v;
    my $coordsref = $fields[0];
    push( @coords, [ @$coordsref ] );
  }

  my ( @boxpointxy, @boxpointxz, @boxpointyz );
  foreach my $gridpoint ( @coords )
  {
    my ( @point_xys, @point_xzs, @point_yzs );
    foreach ( @$gridpoint )
    {
      push ( @point_xys, [ $_->[0], $_->[1] ] );
      push ( @point_xzs, [ $_->[0], $_->[2] ] );
      push ( @point_yzs, [ $_->[1], $_->[2] ] );
    }
    push ( @boxpointxy, [ @point_xys ] );
    push ( @boxpointxz, [ @point_xzs ] );
    push ( @boxpointyz, [ @point_yzs ] );
  }

  my ( @boxvertxy, @boxvertxz, @boxvertyz );
  foreach my $surf ( @verts )
  {
    my ( @vert_xys, @vert_xzs, @vert_yzs );
    foreach my $vert ( @$surf )
    {
      push ( @vert_xys, [ $vert->[0], $vert->[1] ] );
      push ( @vert_xzs, [ $vert->[0], $vert->[2] ] );
      push ( @vert_yzs, [ $vert->[1], $vert->[2] ] );
    }
    push ( @boxvertxy, [ @vert_xys ] );
    push ( @boxvertxz, [ @vert_xzs ] );
    push ( @boxvertyz, [ @vert_yzs ] );
  }

  my $count = 0;
  my ( $vert_xys, $vert_xzs, $vert_yzs, $polyxy, $polyxz, $polyyz );
  my ( @verts_xys, @verts_xzs, @verts_yzs);
  foreach my $case ( @boxvertxy )
  {
    $vert_xys = $boxvertxy[$count];
    $vert_xzs = $boxvertxz[$count];
    $vert_yzs = $boxvertyz[$count];
    $polyxy = Math::Polygon::Tree->new( $vert_xys );
    $polyxz = Math::Polygon::Tree->new( $vert_xzs );
    $polyyz = Math::Polygon::Tree->new( $vert_yzs );
    $count++;
  }

  my $count = 0;
  my @newbox;
  foreach my $caseref ( @gridpoints )
  {
    my @case = @$caseref;
    my $surfnum = $case[ 1 ];
    my $dirvector = $case[ 2 ];
    my @bag;
    foreach my $vert ( @{ $case[ 0 ] } )

    {
      my $xyref = $boxpointxy[ $count ][ 0 ];
      my $xzref = $boxpointxz[ $count ][ 0 ];
      my $yzref = $boxpointyz[ $count ][ 0 ];
      unless ( ( ( $polyxy->contains( $xyref ) ) == 0 ) and ( ( $polyxz->contains( $xzref ) ) == 0 ) and ( ( $polyyz->contains( $yzref ) ) == 0 ) )
      {
        push( @bag, $vert );
      }
    }
    push ( @newbox, [ [ @bag ], $surfnum, $dirvector ] );
    $count++;
  }
  return ( @newbox );
}


sub solveselective
{
  my ( $matdbfile_f2, $selectives_ref, $conffile, $conffile_f2, $path ) = @_;
  my @selectives = @{ $selectives_ref };

  my $matdbfile_f3 = $matdbfile_f2;
  $matdbfile_f3 =~ s/_f2/_f3/ ; #say "\$matdbfile_f3 $matdbfile_f3";
  my $matdbfile_f4 = $matdbfile_f2;
  $matdbfile_f4 =~ s/_f2/_f4/ ; #say "\$matdbfile_f4 $matdbfile_f4";

  my $shortmatdbfile_f2 = $matdbfile_f2;
  my $shortmatdbfile_f3 = $matdbfile_f3;
  my $shortmatdbfile_f4 = $matdbfile_f4;
  $shortmatdbfile_f2 =~ s/$path\/dbs\/// ;
  $shortmatdbfile_f3 =~ s/$path\/dbs\/// ;
  $shortmatdbfile_f4 =~ s/$path\/dbs\/// ;

  open( my $MATDBFILE_F2, "$matdbfile_f2" ) or die;
  my @matlines = <$MATDBFILE_F2>;
  close $MATDBFILES_F2;

  my @mats = map{ $_->[0] } @selectives; #say "\@mats @mats";

  my ( @row, @thirdloop, @fourthloop );
  my $semaphore = "off";
  foreach my $matline ( @matlines )
  {
    chomp $matline;
    $matline =~ s/\s+// ;
    my @e = split( ",", $matline );
    if ( $e[0] eq "*item" )
    {
      if ( $e[1] ~~ @mats )
      {
        $semaphore = "on";
      }
      else
      {
        $semaphore = "off";
      }
    }

    if ( ( $e[0] =~ /^\d/ ) and ( $e[-1] =~ /\D$/ ) )
    {
      if ( $semaphore eq "on" )
      {
        my $ratio = $selective->[1];
        my $changeratio3 = ( 1 - ( ( $ratio - 1 ) / 2 ) ); #say "\$changeratio3 $changeratio3";
        my $changeratio4 = ( 1 + ( ( $ratio - 1 ) / 2 ) ); #say "\$changeratio4 $changeratio4";
        my $e5_3 = $e[5] * $changeratio3;
        my $e5_4 = $e[5] * $changeratio4;
        my $e6_3 = $e[6]  * $changeratio3;
        my $e6_4 = $e[6] * $changeratio4;

        my $lin = "$e[0],$e[1],$e[2],$e[3],$e[4],$e5_3,$e6_3,$e[7],$e[8],$e[9]";
        push( @thirdloop, $lin );
        my $linn = "$e[0],$e[1],$e[2],$e[3],$e[4],$e5_4,$e6_4,$e[7],$e[8],$e[9]";
        push( @fourthloop, $linn );
      }
      elsif ( $semaphore eq "off" )
      {
        push( @thirdloop, $matline );
        push( @fourthloop, $matline );
      }
    }
    else
    {
      push( @thirdloop, $matline );
      push( @fourthloop, $matline );
    }

    open( my $MATDBFILE_F3, ">$matdbfile_f3" ) or die;
    foreach my $line ( @thirdloop )
    {
      say $MATDBFILE_F3 $line ;
    }
    close $MATDBFILE_F3;

    open( my $MATDBFILE_F4, ">$matdbfile_f4" ) or die;
    foreach my $line ( @fourthloop )
    {
      say $MATDBFILE_F4 $line ;
    }
    close $MATDBFILE_F4;
  }

  my $conffile_f4 = $conffile;
  $conffile_f4 =~ s/\.cfg/\_f4\.cfg/;
  print REPORT "cp -R -f $conffile_f2 $conffile_f4\n";
  `cp -R -f $conffile_f2 $conffile_f4\n`;

  open( my $CONFFILE_F2, "$conffile_f2" ) or die;
  my @lines2 =<$CONFFILE_F2>;
  close $CONFFILE_F2;

  open( my $CONFFILE_F2, ">$conffile_f2" ) or die;
  foreach my $line2 ( @lines2 )
  {
    $line2 =~ s/$shortmatdbfile_f2/$shortmatdbfile_f3/ ;
    print $CONFFILE_F2 $line2;
  }
  close $CONFFILE_F2;

  open( my $CONFFILE_F4, "$conffile_f4" ) or die;
  my @lines4 =<$CONFFILE_F4>;
  close $CONFFILE_F4;

  open( my $CONFFILE_F4, ">$conffile_f4" ) or die;
  foreach my $line4 ( @lines4 )
  {
    $line4 =~ s/$shortmatdbfile_f2/$shortmatdbfile_f4/ ;
    print $CONFFILE_F4 $line4;
  }
  close $CONFFILE_F4;

  return( $conffile_f4 );
}


sub getsolar
{
  my ( $paths_ref ) = @_;
  my %paths = %{ $paths_ref };
  my $clma = $paths{clmfilea};
  my $lat = $paths{lat};
  my $longdiff = $paths{longdiff};
  my $clmavgs = $paths{clmavgs};
  my $lstm = $paths{standardmeridian};
  my $long = $lstm + $longdiff;

  my %daymonths = ( 1 => 16, 2 => 15, 3 => 16, 4 => 16, 5 => 16, 6 => 16, 7 => 16, 8 => 16, 9 => 16, 10 => 16, 11 => 16, 12 => 16 );
  my $daynumber;

  open ( WFILE, "$clma" ) or die;
  my @wlines = <WFILE>;
  close WFILE;

  my $countline = 0;
  my $count = -1;
  my $c = -1;
  my ( %ts, $month, $day, $hour, $dir, $diff );
  my $sem = 1;
  foreach my $line ( @wlines )
  {
    chomp $line;
    if ( $countline >= 12 )
    {
      my @elts;
      my $sem = 0;
      $line =~ s/  / / ;
      $line =~ s/  / / ;
      $line =~ s/ /,/ ;
      $line =~ s/^\,// ; #say $line;
      my @elts = split( "  | |,", $line ); #say "ELTS: " . dump ( @elts );

      if ( $elts[0] eq "*" )
      {
        $sem = 0;
        $month = $elts[4]; #say "MONTH: " . dump ( $month );
        $day = $elts[2]; #say "\$day: " . dump ( $day );
        $count = 0;
        if ( $countline == 12 )
        {
          $c = 0;
        }
      }
      else
      {
        $sem = 1;
      }

      if ( ( $sem == 1 ) and ( $month ne "" ) )
      {
        $diff = $elts[0];
        $dir = $elts[2];
        $hour = $count;
        $ts{or}{dir}{$month}{$day}{$hour} = $dir;
        $ts{or}{diff}{$month}{$day}{$hour} = $diff; #say "\$diff: $diff";
        $count++;
      }
      $c++;
    }
    $countline++;
  } #say "TS: " . dump ( \%ts );

  if ( $c <= 8776 )
  {
    %daynums = ( 1 => 16, 2 => 46, 3 => 75, 4 => 106, 5 => 136, 6 => 167, 7 => 197, 8 => 228, 9 => 259, 10 => 289, 11 => 320, 12 => 351 );
  }
  elsif ( $c >= 8790 )
  {
    %daynums = ( 1 => 16, 2 => 46, 3 => 76, 4 => 107, 5 => 137, 6 => 168, 7 => 198, 8 => 229, 9 => 260, 10 => 290, 11 => 321, 12 => 352 );
  }

  my @dds = ( "dir", "diff" );
  my ( %t, %ti );
  foreach my $dd ( @dds )
  { #say "\$dd: " . dump ( $dd );
    foreach my $m ( sort { $a <=> $b } ( keys %{ $ts{or}{$dd} } ) )
    {
      my @bag;
      foreach my $d ( sort { $a <=> $b } ( keys %{ $ts{or}{$dd}{$m} } ) )
      { #say "\$d: " . dump ( $d );
        foreach my $h ( sort { $a <=> $b } ( keys %{ $ts{or}{$dd}{$m}{$d} } ) )
        {
          push ( @{ $t{vals}{$dd}{$m}{$h} }, $ts{or}{$dd}{$m}{$d}{$h} );
        }
      }
    }
  } #say "T: " . dump ( \%t );


  foreach my $dd ( @dds )
  {
    foreach my $m ( sort { $a <=> $b} ( keys %{ $t{vals}{$dd} } ) )
    {
      foreach my $h ( sort { $a <=> $b} ( keys %{ $t{vals}{$dd}{$m} } ) )
      {
        #say "h: " . dump ( $h );
        my @ddvals = @{ $t{vals}{$dd}{$m}{$h} };
        my $ddval = mean( @ddvals );
        $t{avg}{$dd}{$m}{$h} = $ddval;
      }
    }
  } #say "\T: " . dump ( $t{avg} );


  $lat, $longdiff, $long, $localtime;
  open( NEWCLM, ">$clmavgs" ) or die;
  foreach my $m ( sort { $a <=> $b} ( keys %{ $t{avg}{dir} } ) )
  {
    foreach my $h ( sort { $a <=> $b} ( keys %{ $t{avg}{dir}{$m} } ) )
    { #say "\T: " . dump ( $t{avg}{dir}{$m}{$h} );

      my $decl = 23.45 * sin( deg2rad( 280.1 + 0.9863 * $daynums{$m} ) );
      my $declrad = deg2rad($decl);

      my $timeq = ( 9.87 * sin( deg2rad( ( 1.978 * $daynums{$m}  )- 160.22 ) ) ) - ( 7.53 * cos( deg2rad( ( 0.989 * $daynums{$m} )- 80.11 ) ) )
        - ( 1.5 * sin( deg2rad( ( 0.989 * $daynums{$m}  ) - 80.11 ) ) );  #say "\$timeq: $timeq";

      my $tcf = ( 4 * ( $lstm - $long ) ) + $timeq; #say "\$tcf: $tcf";

      #my $solartime = ( $h + ( $longdiff / 15 ) + $timeq ); #say "\$solartime1: $solartime";

      my $solartime = ( $h + ( $tcf / 60 ) ); #say "\$solartime2: $solartime";

      my $hourangle = ( 15 * ( 12 - $solartime ) );
      my $houranglerad = deg2rad($hourangle);

      my $latrad = deg2rad($lat);

      my $altrad = asin( ( cos( $latrad ) * cos( $declrad ) * cos( $houranglerad ) ) + ( sin( $latrad ) * sin( $declrad ) ) );
      my $alt = rad2deg($altrad);
      $alt = sprintf ( "%.3f", $alt );
      $t{avg}{alt}{$m}{$h} = $alt;

      my $azirad = ( asin( cos( $declrad ) * ( sin( $houranglerad ) / cos( $altrad ) ) ) );
      my $azi = rad2deg($azirad);
      $azi = sprintf ( "%.3f", $azi );
      $t{avg}{azi}{$m}{$h} = $azi;

      say NEWCLM "$m,$daymonths{$m},$h,$t{avg}{dir}{$m}{$h},$t{avg}{diff}{$m}{$h},$alt,$azi";
      #say NEWCLM "\$m:$m,\$daymonths{\$m}:$daymonths{$m},\$h:$h,\$t{avg}{dir}{\$m}{\$h}:$t{avg}{dir}{$m}{$h},\$t{avg}{diff}{\$m}{\$h}:$t{avg}{diff}{$m}{$h},\$alt:$alt,\$azi:$azi";
    }
  }
  close NEWCLM;
  return( \%t );
}


sub modish
{ # MAIN PROGRAM
  my @things = @_;
  my $modishdefpath;
  my %paths;

  my $launchfile = shift( @things ); #say "\$launchfile: " . $launchfile; ###.

  my $path = $launchfile;
  $path =~ s/\.cfg$// ;
  while ( not ( $path =~ /\/$/ ) )
  {
    $path =~ s/(\w+)$// ;
  }
  $path =~ s/\/$// ;

  if ( $path =~ /\/cfg$/ )
  {
    $path =~ s/\/cfg$// ;
  }
  #say "\$path: " . $path; ###.

  my ( @restpars, @settings, @received );
  my @received = @things;

  if ( not ( @ARGV) )
  {
    foreach ( @received )
    {
      if ( not ( ref( $_ ) ) )
      {
        push( @restpars, $_ );
      }
      else
      {
        @settings = @{ $_ };
      }
    }
  }
  else
  {
    foreach ( @received )
    {
      if ( not ( ref( $_ ) ) )
      {
        push( @restpars, $_ );
      }
      else
      {
        push( @settings );
      }
    }
  }

  if ( scalar( @restpars ) == 0 ) { say "NO ZONE HAVE BEEN SPECIFIED. EXITING." and die; }

  my ( $zonenum, $dirvectorsnum, $bounceambnum, $bouncemaxnum, $distgrid, $threshold );
  my ( @transpdata, @surfaces, @dirvectorsrefs, @transpsurfs, @resolutions, @treatedlines );

  say "Setting things up...\n";

  #my $radpath = $path . "/rad";


  ##################################################

  my $zonenum = $restpars[0];
  my @transpsurfs = @restpars[ 1..$#restpars ];

  if ( scalar( @{ $settings } ) == 0 )
  {
    if ( -e "$modishdefpath" )
    {
      require "$modishdefpath";
      @resolutions = @{ $defaults[0] };
      $dirvectorsnum = $defaults[1];
      $bounceambnum = $defaults[2];
      $bouncemaxnum = $defaults[3];
      $distgrid = $defaults[4];
      $threshold = $defaults[5];
    }
    elsif ( -e "./modish_defaults.pl" )
    {
      require "./modish_defaults.pl";
      @resolutions = @{ $defaults[0] };
      $dirvectorsnum = $defaults[1];
      $bounceambnum = $defaults[2];
      $bouncemaxnum = $defaults[3];
      $distgrid = $defaults[4];
      $threshold = $defaults[5];
    }
  }
  else
  {
    @resolutions = @{ $settings[0] };
    $dirvectorsnum = $settings[1];
    $bounceambnum = $settings[2];
    $bouncemaxnum = $settings[3];
    $distgrid = $settings[4];
    $threshold = $settings[5];
  }
  my $writefile = "$path/writefile.txt";
  open ( REPORT, ">$writefile" ) or die "Can't open $writefile !";

  if ( not ( @resolutions ) ) { @resolutions = ( 2, 2 ); };
  if ( not defined( $dirvectorsnum ) ) { $dirvectorsnum = 1; };
  if ( not defined( $bounceambnum ) ) { $bounceambnum = 1; };
  if ( not defined( $bouncemaxnum ) ) { $bouncemaxnum = 7; };
  if ( not defined( $distgrid ) ) { $distgrid = 0.01; };
  if ( not defined( $threshold ) ) { $threshold = 0.99; };
  if ( not ( @calcprocedures ) ) { @calcprocedures = ( "diluted", "gensky" ); }
  if ( not ( @specularratios ) ) { @specularratios = ( ) }
  if ( !keys %skycondition )
  {
    %skycondition = ( 1=> "clear", 2=> "clear", 3=> "clear", 4=> "clear", 5=> "clear", 6=> "clear", 7=> "clear", 8=> "clear", 9=> "clear", 10=> "clear", 11=> "clear", 12=> "clear" );
  }

  #if ( not defined( $max_processes ) ) { $max_processes = 1; };
  push ( @calcprocedures, "besides", "extra" ); # THESE SETTINGS WERE ONCE SPECIFIABLE IN THE CONFIGURATION FILE.


# Debug output from ESP-r (out.txt in /cfg and /rad), $debug = 1 to enable.
  my $debug = 0;
  if ( $debug == 1 )
  {
    say REPORT "ESP-r debug output activated.";
    if ( -e "$path/cfg/out.txt" )
    {
      say REPORT "rm $path/cfg/out.txt";
      `rm $path/cfg/out.txt`;
    }
    say REPORT "touch $path/cfg/out.txt";
    `touch $path/cfg/out.txt`;
    if ( -e "$path/rad/out.txt" )
    {
      say REPORT "rm $path/rad/out.txt";
      `rm $path/rad/out.txt`;
    }
    say REPORT "touch $path/rad/out.txt";
    `touch $path/rad/out.txt`;
  }

  my ( $conffile, $conffile_f1, $conffile_f2, $conffile_f3, $constrdbfile, $constrdbfile_f,
  $matdbfile, $matdbfile_f1, $matdbfile_f2, $flagconstrdb, $flagmatdb, $flaggeo, $flagconstr, $originalsref,
  $fictitia1ref, $fictitia2ref, $fictitia3ref, $paths_ref ) = createfictitiousfiles( $launchfile, $path, $zonenum, \@calcprocedures );
  my %paths = %{ $paths_ref }; #say "PATH RECEIVED2: " . dump( \%paths ); ###.
  my $radpath = $paths{radpath};

  my @basevectors = getbasevectors( $dirvectorsnum );

  my @originals = @$originalsref;
  my @fictitia1 = @$fictitia1ref;
  my @fictitia2 = @$fictitia2ref;
  my @fictitia3 = @$fictitia3ref;

  my ( @daylighthours);
  my %actiondata;
  my @zoneoriginals = @originals;
  shift(@zoneoriginals); shift(@zoneoriginals);
  my @zonefictitia1 = @fictitia1; # "BLACK" MODEL
  shift(@zonefictitia1); shift(@zonefictitia1);
  my @zonefictitia2 = @fictitia2; # "REFLECTIVE" MODEL
  shift(@zonefictitia2); shift(@zonefictitia2);
  my @zonefictitia3 = @fictitia3; # "REFLECTIVE" MODEL
  shift(@zonefictitia3); shift(@zonefictitia3);

  my ( %zonefilelists, %fict1filelists, %fict2filelists, %fict3filelists, %fict4filelists );
  my @daylighthoursarr;
  my %daylighthours;
  my ( $exportreflref__, $exportconstrref__ );

  my $tempmod = "$launchfile.mod.temp";
  my ( $tempmoddir, $tempreportdir );
  if ( $^O eq "linux" ) # THESE LINES DEPEND FROM THE OPERATING SYSTEM.
  {
    if ( $tempmod =~ /$path\/cfg\// )
    {
      $tempmod =~ s/$path\/cfg\///;
    }
    else
    {
      $tempmod =~ s/$path\///;
    }

    $tempmod = "$path/tmp/$tempmod";
    unless ( -e "$path/tmp" )
    {
      `mkdir $path/tmp`;
      say REPORT "mkdir $path/tmp";
    }
  }
  elsif ( $^O eq "darwin" )
  { ; }

  say REPORT "\$tempmod $tempmod";
  open ( TEMPMOD, ">$tempmod" ) or die "$!";

  my $tempreport = "$launchfile.report.temp";
  if ( $^O eq "linux" ) # THESE LINES DEPEND FROM THE OPERATING SYSTEM.
  {
    if ( $tempreport =~ /$path\/cfg\// )
    {
      $tempreport =~ s/$path\/cfg\///;
    }
    else
    {
      $tempreport =~ s/$path\///;
    }

    $tempreport = "$path/tmp/$tempreport";
  }
  elsif ( $^O eq "darwin" )
  { ; }

  open ( TEMPREPORT, ">$tempreport" ) or die "$!";

  $tempmoddir = $tempmod . ".dir";
  open ( TEMPMODDIR, ">$tempmoddir" ) or die "$!";

  $tempreportdir = $tempreport . ".dir";
  open ( TEMPREPORTDIR, ">$tempreportdir" ) or die "$!";

  my @treatedlines;

  `cp -f ./fix.sh $radpath/fix.sh`;
  say REPORT "cp -f ./fix.sh $radpath/fix.sh\n";
  `chmod 755 $radpath/fix.sh`;
  `chmod 755 ./fix.sh`;
  say REPORT "chmod 755 $radpath/fix.sh\n";
  `cp -f ./perlfix.pl $radpath/perlfix.pl`;
  say REPORT "cp -f ./perlfix.pl $radpath/perlfix.pl\n";
  #say "PATHS BEFORE BEGINNING: " . dump ( \%paths );###. #############################################################
  my $countzone = 1;
  foreach my $elt (@zoneoriginals)
  {
    my @zonefiles = @$elt;
    my @fict1files = @{ $zonefictitia1[ $countzone - 1 ] };
    my @fict2files = @{ $zonefictitia2[ $countzone - 1 ] };
    my @fict3files = @{ $zonefictitia3[ $countzone - 1 ] };
    my $geofile = $zonefiles[0];
    my $constrfile = $zonefiles[1];
    my $shdfile = $zonefiles[2];
    my $zonenum_cfg = $zonefiles[3];
    my $geofile_f = $fict1files[0];
    my $constrfile_fict = $fict1files[1];
    $zonefilelists{ $zonenum }{ geofile } = $geofile;
    $zonefilelists{ $zonenum }{ geofile_f } = $geofile_f;
    $zonefilelists{ $zonenum }{ constrfile } = $constrfile;
    $zonefilelists{ $zonenum }{ constrfile_f } = $constrfile_f;
    $zonefilelists{ $zonenum }{ shdfile } = $shdfile;

    my ( $transpeltsref, $geofilestructref, $surfslistref, $obsref, $obsconstrsetref, $datalistref,
      $obsmaterialsref, $orresolutions_ref ) =
      readgeofile( $geofile, \@transpsurfs, $zonenum, \@calcprocedures );

    my @transpelts = @$transpeltsref;
    my @geodata = @$geofilestructref;
    my %surfslist = %$surfslistref;
    my @obsdata = @$obsref;
    my @obsconstrset = @$obsconstrsetref;
    my %datalist = %$datalistref;
    my @obsmaterials = @{ $obsmaterialsref };
    my @orresolutions = @{ $orresolutions_ref };

    createfictgeofile( $geofile, \@obsconstrset, $geofile_f, \%paths, \@calcprocedures, $conffile_f1 );

    setroot( $conffile_f1, $path, $debug, \%paths );
    setroot( $conffile_f2, $path, $debug, \%paths );
    if ( "coexistent" ~~ @calcprocedures )
    {
      setroot( $conffile_f3, $path, $debug, \%paths );
    }
    if ( scalar( @selectives) > 0 )
    {
      setroot( $conffile_f4, $path, $debug, \%paths );
    }

    my ( $materialsref, $newmaterialsref, $matnumsref, $newmatnumsref, $exportconstrref ) =
    createconstrdbfile( $constrdbfile, $constrdbfile_f, \@obsconstrset );

    my ( $exportreflref, $obslayers_ref, $selectives_ref ) = creatematdbfiles( $matdbfile,
      $matdbfile_f1, $matdbfile_f2, \@calcprocedures, $constrdbfile_f, \@obsdata );
    my %obslayers = %{ $obslayers_ref };

    my @selectives;
    foreach my $item ( @calcprocedures )
    {
      my @els = split( ":", $item );
      if ( $els[0] eq "light/infrared-ratio" )
      {
        my $mat = $els[1];
        my $ratio = $els[2];
        push( @selectives, [ $mat, $ratio ] );
      }
    }
    @selectives = uniq( @selectives );

    my $conffile_f4;
    if ( scalar( @selectives ) > 0 )
    {
      $conffile_f4 = solveselective( $matdbfile_f2, \@selectives, $conffile, $conffile_f2, $path );
    }
    my ( $surfnumsref, $surfnamesref ) = tellsurfnames( \@transpsurfs, \@geodata );
    my @surfnums = @$surfnumsref;
    my @surfnames = @$surfnamesref;
    my ( $winseltsref, $datalistref ) = readverts( \@transpelts, $geofile, \@geodata, \%datalist );
    my @winselts = @$winseltsref;
    my %datalist = %$datalistref;
    my ( $winscoordsref, $datalistref ) = readcoords( \@winselts, $geofile, \@geodata, \%datalist, \@transpelts );
    my @winscoords = @$winscoordsref;
    my %datalist = %$datalistref;
    my @dirvectorsrefs = calcdirvectors( @winscoords );
    my @xyzcoords = getcorners( \@winscoords, \@winselts );
    my @extremes = findextremes( @xyzcoords );
    my @gridcoords = makecoordsgrid( \@extremes, \@resolutions, \@dirvectorsrefs );
    my @gridpoints_transitional = makegrid( @gridcoords );
    my @gridpoints_newtransitional = prunepoints( \@gridpoints_transitional, \@xyzcoords );
    my @gridpoints = adjustgrid( \@gridpoints_newtransitional, $distgrid );

    my ( @orgridcoords, @orgridpoints_transitional, @orgridpoints_newtransitional, @orgridpoints );
    if ( "espdirres" ~~ @calcprocedures )
    {
      @orgridcoords = makecoordsgrid( \@extremes, \@orresolutions, \@dirvectorsrefs );
      @orgridpoints_transitional = makegrid( @orgridcoords );
      @orgridpoints_newtransitional = prunepoints( \@orgridpoints_transitional, \@xyzcoords );
      @orgridpoints = adjustgrid( \@orgridpoints_newtransitional, $distgrid );
    }

    my ( $treatedlinesref, $filearrayref, $monthsref ) = readshdfile( $shdfile, \@calcprocedures, $conffile_f2 );
    @treatedlines = @$treatedlinesref;
    my @shdfilearray = @$filearrayref;
    my @months = @$monthsref;
    my @shdsurfdata = getsurfshd( \@shdfilearray, \@months, \@surfnums, \@surfnames );
    @daylighthoursarr = checklight( \@shdfilearray, \@months );
    %daylighthours = populatelight( @daylighthoursarr );
    $shdfileslist{ $zonenum } = \@treatedlines;
    $countzone++;
    my @radfilesrefs = tellradfilenames( $path, $conffile_f1, $conffile_f2, $conffile_f3, $conffile_f4, \%paths );
    my $hashirrsref = pursue( { zonenum => $zonenum, geofile => $geofile, constrfile => $constrfile,
      shdfile => $shdfile, gridpoints => \@gridpoints, shdsurfdata => \@shdsurfdata,
      daylighthoursarr => \@daylighthoursarr, daylighthours=> \%daylighthours,
      shdfilearray => \@shdfilearray, exportconstrref => $exportconstrref,
      exportreflref => $exportreflref, conffile => $conffile,  path => $path,
      radpath => $radpath, basevectors => \@basevectors, resolutions => \@resolutions,
      dirvectorsnum => $dirvectorsnum, calcprocedures => \@calcprocedures,
      specularratios => \@specularratios, bounceambnum => $bounceambnum,
      bouncemaxnum => $bouncemaxnum, radfilesrefs => \@radfilesrefs,
      conffile_f1 => $conffile_f1, conffile_f2 => $conffile_f2, conffile_f3 => $conffile_f3, conffile_f4 => $conffile_f4,
      transpsurfs=> \@transpsurfs, selectives => \@selectives, paths => \%paths, orgridpoints => \@orgridpoints } );

    my $irrvarsref = compareirrs( \%zonefilelists, $hashirrsref, $computype, \@calcprocedures, \@selectives );

    foreach my $elm ( @transpsurfs )
    {
      my @transpsurfs;
      push ( @transpsurfs, $elm );
      say "Closing calculations for surface " . dump( @transpsurfs );
      modifyshda( \@comparedirrs, \%surfslist, \%zonefilelists, \%shdfileslist, \%daylighthours, $irrvarsref, $threshold, $tempmod, $tempreport,  $tempmoddir, $tempreportdir, $elm, "diffuse", \@calcprocedures, $hashirrsref, $conffile_f2, $shdfile );

      modifyshda( \@comparedirrs, \%surfslist, \%zonefilelists, \%shdfileslist, \%daylighthours, $irrvarsref, $threshold, $tempmod, $tempreport,  $tempmoddir, $tempreportdir, $elm, "direct", \@calcprocedures, $hashirrsref, $conffile_f2, $shdfile );
    }
  }

  close TEMPMOD;
  open ( TEMPMOD, "$tempmod" ) or die;
  my @tempmodlines = <TEMPMOD>;
  close TEMPMOD;
  @tempmodlines = uniq( @tempmodlines );

  close TEMPREPORT;
  open ( TEMPREPORT, "$tempreport" ) or die;
  my @tempreportlines = <TEMPREPORT>;
  close TEMPREPORT;
  @tempreportlines = uniq( @tempreportlines );

  close TEMPMODDIR;

  open ( TEMPMODDIR, "$tempmoddir" ) or die;
  my @tempmoddirlines = <TEMPMODDIR>;
  close TEMPMODDIR;
  @tempmoddirlines = uniq( @tempmoddirlines );

  close TEMPREPORTDIR;
  open ( TEMPREPORTDIR, "$tempreportdir" ) or die;
  my @tempreportdirlines = <TEMPREPORTDIR>;
  close TEMPREPORTDIR;
  @tempreportdirlines = uniq( @tempreportdirlines );

  setroot( $launchfile, $path, $debug, \%paths );

  my $shdfile = $zonefilelists{ $zonenum }{ shdfile };
  my $shdafile = "$shdfile" . "a";
  my $shdafilemod = $shdafile;
  $shdafilemod =~ s/.shda/.mod.shda/;
  open ( SHDAMOD, ">$shdafilemod" ) or die;

  my $shdafilereport = $shdafile;
  $shdafilereport =~ s/.shda/.report.shda/;
  print REPORT "cp -R -f $shdafile $shdafilereport";
  `cp -R -f $shdafile $shdafilereport`;
  open ( SHDAREPORT, ">>$shdafilereport" ) or die;
  print SHDAREPORT "# FOLLOWING, THE VERIFIED VARIATIONS (AS RATIOS) OF IRRADIANCES DUE TO REFLECTIONS BY OBSTRUCTIONS.\n";

  my $counter = 0;;
  foreach my $lin ( @tempreportlines )
  {
    my $lindir = $tempreportdirlines[ $counter ];
    print SHDAREPORT $lindir;
    print SHDAREPORT $lin;
    $counter++;
  }
  close SHDAREPORT;

  my @container;
  my ($currentmonth, $signal );
  foreach my $lin ( @treatedlines )
  {
    my @arr = split(/\s+|,/, $lin);
    if ( $lin =~ /^\* month:/ )
    {
      $currentmonth = $arr[2];
      $currentmonth =~ s/`//g;
    }

    if ( ( $lin =~ /24 hour external surface shading/ ) or ( $lin =~ /\* month:/ ) or ( $lin =~ /\* end/ ) or ( $lin =~ /\* end/ ) or ( $lin =~ /Shading and insolation data in db/ ) )
    {
      $signal = "on";
    }
    elsif ( $lin =~ /24 hour internal surface insolation/ )
    {
      $signal = "off";
    }

    my $count = 0;
    my $i = 0;
    foreach my $el ( @tempmodlines )
    {
      my $eldir = $tempmoddirlines[ $i ];
      my @modarrdir = split( /\s+|,/, $eldir );
      my @modarr = split( /\s+|,/, $el );
      if ( $signal eq "on" )
      {
        push( @arr, $currentmonth );
        if ( ( $arr[25] eq $modarrdir[25] ) and ( $arr[28] eq $modarrdir[28] ) and ( $arr[29] eq $modarrdir[29] ) )
        {
          $eldir =~ s/ $currentmonth?//;
          push( @container, $eldir );
          $count++;
          last;
        }
        elsif ( ( $arr[25] eq $modarr[25] ) and ( $arr[28] eq $modarr[28] ) and ( $arr[29] eq $modarr[29] ) )
        {
          $el =~ s/ $currentmonth?//;
          push( @container, $el );
          $count++;
          last;
        }
      }
      $i++;
    }

    if ( $count == 0 )
    {
      push( @container, $lin );
    }

  }


  my $signalins;

  foreach my $lin ( @container )
  {
    if ( ( $lin =~ / # diffuse - / ) or ( $lin =~ / # direct - / ) )
    {
      my  @arr = split(/\s+|,/, $lin);
      my @firstarr = @arr[ 0..11 ];
      my @secondarr = @arr[ 12..$#arr ];
      my $joinedfirst = join ( ' ' , @firstarr );
      my $joinedsecond = join ( ' ' , @secondarr );
      $lin = "$joinedfirst\n" . "$joinedsecond\n";
    }

    if ( "noins" ~~ @calcprocedures )
    {
      if ( ( $lin =~ /24 hour external surface shading/ ) or
        ( $lin =~ /\* month:/ ) or ( $lin =~ /\* end/ ) or
        ( $lin =~ /Shading and insolation data in db/ ) )
      {
        $signalins = "off";
        #say REPORT "signalins off";
      }
      elsif ( $lin =~ /24 hour internal surface insolation/ )
      {
        $signalins = "on";
        #say REPORT "signalins on";
      }
    }

    unless ( $signalins eq "on" )
    {
      print SHDAMOD $lin;
      #say REPORT "I AM GOING TO PRINT THIS IN $shdafilemod: " . "$lin, because \$signalins is $signalins." ;
    }

  }
  close SHDAMOD;
}

if ( @ARGV )
{
  modish( @ARGV );
}

1;
