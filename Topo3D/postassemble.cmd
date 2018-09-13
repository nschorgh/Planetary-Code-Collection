#!/bin/csh

# merges and sorts output files from parallel runs of programs shadows and fieldofviews

#/bin/rm -f horizons.dat fieldofviews.dat
#touch horizons.dat fieldofviews.dat
/bin/rm -f horizons.dat 
touch horizons.dat

foreach i ( horizon.* )
  cat $i >> horizons.dat
end
#foreach i ( fieldofview.* )
#  cat $i >> fieldofviews.dat
#end

sort -k1n -k2n horizons.dat >! tmptmp
/bin/mv -f tmptmp horizons.dat

#sort -k1n -k2n fieldofviews.dat >! tmptmp
#/bin/mv -f tmptmp fieldofviews.dat




