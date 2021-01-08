# Build SF profile using the following command
./Build-Structure_factor.iso.py -f COLVAR -o sf.dat

# plot sf.dat in GNUPLOT
#p 'sf.dat' u ($1**2):2 w lp
