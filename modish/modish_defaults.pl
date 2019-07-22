# Configuration file for modish, version 1.59.

@defaults = ( [ 2, 2 ], 5, 1, 7, 0.01, 0.01 );

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


@calcprocedures = ( "diluted", "gensky" );
# Quick instructions:
# The best groups of settings of @calcprocedures for calculating the shading factors are likely to be 
# overall the following ones in most cases:
# 1) 
# @calcprocedures = ( "diluted", "gensky", "directlydirect" ); # FOR USING A CIE SKY AND ANALOGUE MODELS BASED ON DIFFERENCES OF REFLECTIVITY. THIS IS THE DEFAULT SETTING, IF NOTHING IS SPECIFIED.
# 2)
# @calcprocedures = ( "diluted", "gendaylit", "getweather", "getsimple" ); # FOR USING A PEREZ SKY BUILT WITH MONTHLY AVERAGE WEATHER DATA AND ANALOGUE MODELS BASED ON DIFFERENCES OF REFLECTIVITY.
# 3) 
# @calcprocedures = ( "diluted", "radical", "gendaylit", "getweather", "getsimple", "keepdirshdf" ); # FOR USING A PEREZ SKY BUILT WITH MONTHLY AVERAGE WEATHER DATA AND ANALOGUE MODELS BASED ON DIFFERENCES OF GEOMETRY.
# 4)
# @calcprocedures = ( "diluted", "gensky", "radical", "coexistent", "directlydirect" ); 
# OR:
# @calcprocedures = ( "diluted", "gendaylit", "radical", "coexistent", "getweather", "getsimple" ); # FOR USING
# STRATEGY 1) OR 2) FOR DIRECT SHADING FACTORS AND STRATEGY 3) FOR DIFFUSE ONES.
# 5)
# @calcprocedures = ( "diluted", "gensky", "radical", "directlydirect", "espdirres" ); # THE "espdirres"
# OPTION MAY BE USED WHEN LOTS OF DIRECT REFLECTIONS COME INTO PLAY
# Explanations follow.
# "diluted" means that the two models from which the shading ratios are derived
# are going to be the following: 
# 1)
# a) a model in which all the surfaces are reflective,
# excepted the obstructions, which are black;
# b) a model in which everything is reflective.
# 2) 
# if "complete" is specified, the two models from which the shading ratios 
# are derived are going to be the following:
# a) a model in which everything is black, and
# b) a model in which all the surfaces are black, excepted the obstructions,
# which are reflective. The settings "diluted" and "complete" are alternatives.
# 3)
# If the strategy "radical" is specified, the following sets of models are going to be used:
# a) a model in which all the surfaces are reflective, and the obstruction are absent;
# b) a model in which everything is reflective, and the obstructions are present.
# The options "diluted" or "complete" can be specified together with
# the option "radical".
# The option "coexistent" specified together with the option "radical" allows
# to utilize strategy 3 ("radical") for the diffuse calculations and strategy 1) or 2)
# for the direct calculations.
# If calcprocedures is set to "keepdirshdf", the direct shading factors calculated by ESP-r's ISH will be kept.
# If the setting "altcalcdiff" is specified, the diffuse irradiances are calculated
# from skies with no sun (be careful). If "directlydirect" is specified, the direct irradiances
# are calculates with skies with suns and imposing 0 diffuse bounces. This is the
# fastest method for direct calculations.
# If either "altcalcdiff" or "directlydirect" are specified, the total irradiances are 
# calculated with skies with sun and complete diffuse bounces. But if both
# "altcalcdiff" or "directlydirect" are specified, the total irradiances are calculated
# as sum of direct and diffuse irradiances.
# If "gensky" is specified, the irradiances are calculated using the gensky program 
# of Radiance, entailing the use of the CIE standard skies, for both the diffuse and direct
# calculations, and the result is sensible to the setting of sky condition for each month (below: 
# clear, cloudy, or overcast).
# If "gendaylit" is specified, the irradiances are calculated using the gendaylit program 
# of Radiance, entailing the use of the Perez sky model for the diffuse
# calculations and the direct ones. If the "getweather" setting is not specified, 
# the direct calculations are performed by the means of gensky. If "getweather" is specified,
# the both the direct and the diffuse calculations are used with gendaylit by the means
# of averages of the weather data about direct normal and horizontal diffuse irradiances.
# For the setting "gendaylit" to work, it has to be specified together with the "altcalcdiff" setting.
# The setting "sparecycles" is specified together with the "radical" setting, it makes
# possible to spare a certain amount of computations (calls to Radiance).
# It ensures that the direct calculations are calculated directly, without
# requesting total irradiance calculations. With gendaylit it is more correct,
# and in any case it is quicker.
# The setting "getweather" used with "gendaylit" ("it can't be used without it) 
# makes possible that the average radiation values of the weather data are utilized 
# when calling gendaylit. 
# The option "getsimple" used with "getweather" (it can't be used without it) 
# determines the fact that the proportion of direct to diffuse radiation 
# is determined directly from the shading data and overriding the other methods
# (determined by "altcalcdiff" and "directlydirect", or their absence) for defining that ratio.
# Note that the materials used in the obstructions should be not shared
# by objects which are not obstructions. If necessary, to obtain that,
# some materials may have to be suitably duplicated and renamed.
# The option "espdirres" makes possible to adopt the esp-r setting for the resolution
# of the direct calculations with Radiance with the "directlydirect" option,
# or "getweather", "getsimple".  
# By specifying in @calcprocedure items of the kind "light/infrared-ratio:materialname:ratio"
# (for example: "light/infrared-ratio:gypsum:1.2" ) it is possible to model
# obstruction material which are selective in reflection - i.e. having different
# reflectivities in the range of light and solar infrared.

@specularratios = (  );

#@specularratios = ( "reflector:mirror" );
#@specularratios = ( "reflector:0.03:0.05" );
# Here values of the kind "construction:specularratio:roughnessvalue"
# may be specified. For example, "reflector:0.03:0.05".
# The textual element ("reflector") is the name
# of a construction. The first number is the specular ratio
# for that construction. The second number is the roughness value.
# Specifying those values here makes possible
# to override the values specified in a Radiance database.
# (for example, the "0"s that may be in the database
# by defaul as regards specular ratios and roughness values).
# As an alternative, a material can be declared to be of the "mirror" type.
# This is done by specifying a value "construction:mirror".
# For example: reflector:mirror (see Radiance documentation
# about the properties of the "mirror" material type).

%skycondition = ( 1=> "clear", 2=> "clear", 3=> "clear", 4=> "clear", 5=> "clear", 6=> "clear", 7=> "clear", 8=> "clear", 9=> "clear", 10=> "clear", 11=> "clear", 12=> "clear" );
# PREVAILING CONDITION OF THE SKY FOR EACH MONTH, EXPRESSED WITH ITS NUMBER, IN THE CASE IN WHICH
# CIE SKIES (WITH GENSKY) ARE UTILIZED.
# THE OPTIONS ARE: "clear", "cloudy" and "overcast".
# IF NO VALUE IS SPECIFIED, THE DEFAULT IS "clear".

