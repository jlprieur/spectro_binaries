#!/bin/csh
############################################################
# plot_histo_log.com
# To display the ascii mass_histo list on a graphic device 
# in log10/log10 corrdinates 
# 
# Syntax:
# plot_histo.com filename device
#      0 = Histogram-like
#      1 = Small dot
#      2 = Open triangle
#      3 = Filled triangle
#      4 = Cross +
#      5 = Cross X
#      6 = Open square
#      7 = Filled square
#      8 = Open circle
#      9 = Filled circle
#
# JLP
# Version 12/03/2007
############################################################
#
#if ($#argv != 3) then
#  echo "Syntax is: plot_histo.com filename graphic_device symbol"
#  goto end
#endif
############################################################ 
# Output of values:
echo "y" >! jlp_lu5.tmp 
# Input of a file:
echo "1" >> jlp_lu5.tmp 
# Profile format:
echo "1" >> jlp_lu5.tmp 
# File name:
echo "histo_mass.dat" >> jlp_lu5.tmp 
# X,Y columns:
echo "1,2,0" >> jlp_lu5.tmp 
# Curve type (L0=solid line)
# echo "L0" >> jlp_lu5.tmp 
# Curve type (43= crosses)
echo "93" >> jlp_lu5.tmp 
# Plot X,Y after transformations:
echo "5" >> jlp_lu5.tmp 
# Log10
echo "2" >> jlp_lu5.tmp 
# Constant to be substracted
echo "0." >> jlp_lu5.tmp 
# Log10
echo "3" >> jlp_lu5.tmp 
# Constant to be substracted
echo "0." >> jlp_lu5.tmp 
# Title:
echo "Histogram" >> jlp_lu5.tmp 
# X Label:
echo "log10 mass" >> jlp_lu5.tmp 
# Y label:
echo "log10 Nstars" >> jlp_lu5.tmp 
# Output device (Upper case => TeX)
echo "&square" >> jlp_lu5.tmp 
# Same scale in X and Y:
echo "n" >> jlp_lu5.tmp 
# Scale: automatic (default), -Y=Yaxis inversion ?
echo "-1.,0.5,4.6,5.2" >> jlp_lu5.tmp 
#echo " " >> jlp_lu5.tmp
# Postcript file? 
# echo "n" >> jlp_lu5.tmp 
# Change window parameters?
echo "n" >> jlp_lu5.tmp 
# Exit:
echo "10" >> jlp_lu5.tmp 
$EXEC/plot1.exe < jlp_lu5.tmp
end:
