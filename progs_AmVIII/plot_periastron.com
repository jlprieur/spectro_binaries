#!/bin/csh
############################################################
# plot_periastron.com
# To display an ascii list on a graphic device 
# X,Y
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
# History:
# runs latex_to_ascii ../submit/tab_orbitals.tex vsini_asini_e.dat 11,8,5
#
############################################################ 
# Output of values:
echo "y" >! jlp_lu5.tmp 
# Input of a file:
echo "1" >> jlp_lu5.tmp 
# Profile format:
echo "1" >> jlp_lu5.tmp 
# File name:
echo "periastron.dat" >> jlp_lu5.tmp 
# X,Y columns:
echo "4,1,0" >> jlp_lu5.tmp 
# Curve type (L0=solid line)
echo "82" >> jlp_lu5.tmp 
# Curve type (43= crosses)
#echo "43" >> jlp_lu5.tmp 
# Simple X,Y plot
echo "4" >> jlp_lu5.tmp 
# Title:
echo "Histogram\0" >> jlp_lu5.tmp 
# X Label:
echo "a sin i (1 - e)" >> jlp_lu5.tmp 
# Y label:
echo " v sin i" >> jlp_lu5.tmp 
# Output device (Upper case => TeX)
echo "&xterm" >> jlp_lu5.tmp 
# Same scale in X and Y:
echo "n" >> jlp_lu5.tmp 
# Scale: automatic (default), -Y=Yaxis inversion ?
# echo " " >> jlp_lu5.tmp 
echo "0.,20.,0.,40." >> jlp_lu5.tmp 
# Postcript file? 
# echo "n" >> jlp_lu5.tmp 
# Change window parameters?
echo "n" >> jlp_lu5.tmp 
# Exit:
echo "10" >> jlp_lu5.tmp 
$EXEC/plot1.exe < jlp_lu5.tmp
end:
