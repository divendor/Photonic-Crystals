reset 
set terminal pngcairo size 1280,1280
#set terminal postscript enhanced
set output "TM_TE_M_GAMMA.png"

set lmargin 5
set rmargin 2
set tmargin 5
set bmargin 5
set size square

set border lw 2

set ytics font ", 20"

set grid linecolor rgb "black" lt -1 lw 0.5
set format x""
set label "{/Symbol M}" at -0.001,-0.04 font ", 35"
set label "{/Symbol G}" at 0.434,-0.04 font ", 35"

set xrange[0:0.44428828]
set yrange[0:0.8]

#set style rect back fs border lc rgb "light-cyan" 
#set object 1 rect from 0,0.46240515 to  0.31415927,0.49132657 fc rgb "light-cyan" behind
#set object 2 rect from 0,0.26198202 to  0.31415927,0.44386894 fc rgb "pink" behind

m = "./TM_M_GAMMA.dat"
n = "./TE_M_GAMMA.dat"
plot for [col=2:6] m using 1:col with lines linecolor rgb "blue" lw 3 notitle, for [col=2:6] n using 1:col with lines linecolor rgb "red" lw 3 notitle



# interesting colours: 
#	red: dark-salmon, pink 
#	blue: light-blue, light-cyan
