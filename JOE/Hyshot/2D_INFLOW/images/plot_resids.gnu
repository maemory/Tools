# gnuplot script to generate residual plots
set term postscript eps color enhanced "Helvetica"
set output "temp.eps"

# plotting info
set grid
set title "2D Hyshot: Shot 810 Cold"
set logscale y
set xlabel 'Iteration #'
set ylabel 'Normalized Residual Value'
set key right box

# scale factor information (user defined for now)
s3=3.1195e-02
s4=6.9337e+01
s5=4.0655e+01
s7=8.4829e+04
s8=6.1869e+02
s9=4.5666e+06

# marker type
marker=1

plot '<grep "RESID:" ../output.SST.log*' u 2:($3/s3)   ti "{/Symbol r}"   w p 0 marker, \
     '<grep "RESID:" ../output.SST.log*' u 2:($4/s4)   ti "{/Symbol r}-u" w p 1 marker, \
     '<grep "RESID:" ../output.SST.log*' u 2:($5/s5)   ti "{/Symbol r}-v" w p 2 marker, \
     '<grep "RESID:" ../output.SST.log*' u 2:($7/s7)   ti "{/Symbol r}-E" w p 3 marker, \
     '<grep "RESID:" ../output.SST.log*' u 2:($8/s8)   ti "{/Symbol k}"   w p 4 marker, \
     '<grep "RESID:" ../output.SST.log*' u 2:($9/s9)   ti "{/Symbol w}"   w p 5 marker
