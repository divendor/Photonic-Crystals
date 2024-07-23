reset 
set terminal pngcairo size 1280,1280
#set terminal postscript enhanced
set output "TM_TE.png"

set lmargin 6
set rmargin 3
set tmargin 3
set bmargin 3
set size ratio 0.8

set border lw 3.5

set ytics font ", 20"

set grid linecolor rgb "black" lt -1 lw 0.5 

set format x""

set label 1 "TM MODES" at 0.33,0.17 font ", 35" tc rgb "blue"
set label 2 "TE MODES" at 0.3399,0.13 font ", 35" tc rgb "red"
set label "{/Symbol G}" at -0.001,-0.04 font ", 35"
set label "{/Symbol C}" at 0.3012,-0.04 font ", 35"
set label "{/Symbol M}" at 0.6099,-0.04 font ", 35"
set label "{/Symbol G}" at 1.049,-0.04 font ", 35"



set xrange[0:1.0726068]
set yrange[0:0.8]

#set style rect back fs border lc rgb "light-cyan" 
#set object 1 rect from 0,0.46240515 to  0.31415927,0.49132657 fc rgb "light-cyan" behind
#set object 1 rect from 0,0.46240515 to  0.31415927,0.50741357 fc rgb "light-cyan" behind
#set object 2 rect from 0,0.26198202 to  0.31415927,0.44386894 fc rgb "pink" behind

set arrow from 0.31415927,0 to 0.31415927,0.8 nohead lc rgb "black" lw 3.5 front
set arrow from 0.62831855,0 to 0.62831855,0.8 nohead lc rgb "black" lw 3.5 front


m = "./TM.dat"
n = "./TE.dat"
plot for [col=2:5] m using 1:col every ::0::180 with lines linecolor rgb "blue" lw 3 notitle, for [col=2:5] n using 1:col every ::0::180 with lines linecolor rgb "red" lw 3 notitle
#plot for [col=2:5] m using 1:col every ::0::180 with lp linecolor rgb "blue" lw 1 notitle, for [col=2:5] n using 1:col every ::0::180 with lines linecolor rgb "red" lw 3 notitle


# interesting colours: 
#	red: dark-salmon, pink 
#	blue: light-blue, light-cyan
