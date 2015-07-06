reset
set terminal x11
set style data lines

set style line 15 pt 7 lc rgb "#AF5000" ps 0.9 lt 1 lw 5 #lu brown
set style line 14 pt 7 lc rgb "#0030FF" ps 0.9 lt 1 lw 5 #lu blue

set style line 1 lt 1 lc rgb "blue" lw 5
set style line 2 lt 1 lc rgb "red" lw 5 pt 7 ps 2
set style line 3 lt 1 lc rgb "green" lw 5
set style line 4 lt 1 lc rgb "yellow" lw 5
set style line 5 lt 1 lc rgb "orange" lw 5
set style line 6 lt 1 lc rgb "purple" lw 5
set style line 7 lt 1 lc rgb "magenta" lw 5
set style line 8 lt 1 lc rgb "cyan" lw 5
set style line 9 lt 1 lc rgb "violet" lw 5
set style line 10 lt 1 lc rgb "black" lw 3 pt 7 ps 2
set style line 11 lt 1 lc rgb "white" lw 3
set style line 12 lt 1 lc rgb "black" lw 15
set style line 13 lt 1 lc rgb "black" lw 1

set style line 16 lt 1 lc rgb "black" lw 3 pt 7 ps 1
set style line 18 lt 1 lc rgb "#FFDDAA" lw 3 pt 7 ps 1
set style line 17 lt 1 lc rgb "#BBCCFF" lw 3 pt 7 ps 1

set xrange [-7:7]
set yrange [-7:7]

pTscale = 0.1

set terminal postscript enhanced
set size square

# 
# "VirtLeftppB1128State.dat" using (($4==0 && $8==1) ? ($2):1/0):3 ls 11 title "DGLAPsafe dipoles",\

set out "AARealState.eps"

#eps_2 is 0.108

plot "AAB55State.dat" using (($4==0) ? (($5==1) ? ($2/5-5.5/2):1/0):1/0):($3/5) with points ls 17 title "Partons from left Nucleus",\
"AAB55State.dat" using (($4==0) ? (($5==0) ? ($2/5-5.5/2):1/0):1/0):($3/5) with points ls 18 title "Partons from right Nucleus",\
"AAB55State.dat" using (($4==0 && $9**2 < 1) ? ($2/5-5.5/2):1/0):($3/5) with points ls 16 title "Partons in {/Symbol h} [-1,1]"

set out "VirtRandCascade.eps"

set xrange [-1.5:1.5]
set yrange [-1.5:1.5]

plot "VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 title "Dipoles, Left Cascade",\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==2) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==3) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==4) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==5) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==6) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==7) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($1==8) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0 && $8==1) ? (($2/5-1.28/2)):1/0):($3/5) ls 11 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==3) ? ($2/5-1.28/2):1/0):($3/5) ls 12 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($5==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 10 title "Partons, Left Cascade",\
"VirtRandLeftppB1128State.dat" using (($4==0) ? (($5==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 2 notitle,\
"VirtRandLeftppB1128State.dat" using (($4==2 && $11==1) ? ($2/5-1.28/2):1/0):($3/5) with points pt 6 lt 1 ps 4 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 title "Dipoles, Right Cascade",\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==2) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==3) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==4) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==5) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==6) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==7) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($1==8) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0 && $8==1) ? (($2/5-1.28/2)):1/0):($3/5) ls 11 notitle,\
"VirtRandRightppB1128State.dat" using (($4==3) ? ($2/5-1.28/2):1/0):($3/5) ls 12 notitle,\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($5==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 2 title "Partons, Right Cascade",\
"VirtRandRightppB1128State.dat" using (($4==0) ? (($5==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 10 notitle,\
"VirtRandRightppB1128State.dat" using (($4==2 && $11==1) ? ($2/5-1.28/2):1/0):($3/5) with points pt 6 lt 1 ps 4 title "Valence Partons"

set out "VirtWeightedCascade.eps"

plot "VirtLeftppB1128State.dat" using (($4==0) ? (($1==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 title "Dipoles, Left Cascade",\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==2) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==3) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==4) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==5) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==6) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==7) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($1==8) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"VirtLeftppB1128State.dat" using (($4==0 && $8==1) ? (($2/5-1.28/2)):1/0):($3/5) ls 11 notitle,\
"VirtLeftppB1128State.dat" using (($4==3) ? ($2/5-1.28/2):1/0):($3/5) ls 12 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($5==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 2 notitle,\
"VirtLeftppB1128State.dat" using (($4==0) ? (($5==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 10 title "Partons, Left Cascade",\
"VirtLeftppB1128State.dat" using (($4==2 && $11==1) ? ($2/5-1.28/2):1/0):($3/5) with points pt 6 lt 1 ps 4 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 title "Dipoles, Right Cascade",\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==2) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==3) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==4) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==5) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==6) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==7) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($1==8) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 15 notitle,\
"VirtRightppB1128State.dat" using (($4==0 && $8==1) ? (($2/5-1.28/2)):1/0):($3/5) ls 11 notitle,\
"VirtRightppB1128State.dat" using (($4==3) ? ($2/5-1.28/2):1/0):($3/5) ls 12 notitle,\
"VirtRightppB1128State.dat" using (($4==0) ? (($5==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 2 title "Partons, Right Cascade",\
"VirtRightppB1128State.dat" using (($4==0) ? (($5==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 10 notitle,\
"VirtRightppB1128State.dat" using (($4==2 && $11==1) ? ($2/5-1.28/2):1/0):($3/5) with points pt 6 lt 1 ps 4 title "Valence Partons"

set out "RealCascade.eps"

plot "RealppB1128State.dat" using (($4==0) ? (($1==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 title "Dipoles, Real Cascade",\
"RealppB1128State.dat" using (($4==0) ? (($1==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($1==2) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($1==3) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($1==4) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($1==5) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($1==6) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($1==7) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($1==8) ? ($2/5-1.28/2):1/0):1/0):($3/5) ls 14 notitle,\
"RealppB1128State.dat" using (($4==0 && $8==1) ? (($2/5-1.28/2)):1/0):($3/5) ls 11 notitle,\
"RealppB1128State.dat" using (($4==3) ? ($2/5-1.28/2):1/0):($3/5) ls 12 notitle,\
"RealppB1128State.dat" using (($4==0) ? (($5==1) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 10 title "Partons, Real Cascade",\
"RealppB1128State.dat" using (($4==0) ? (($5==0) ? ($2/5-1.28/2):1/0):1/0):($3/5) with points ls 2 notitle,\
"VirtRightppB1128State.dat" using (($4==2 && $11==1) ? ($2/5-1.28/2):1/0):($3/5) with points pt 6 lt 1 ps 4 title "Valence Partons"

set terminal x11
set xrange [*:*]
set yrange [*:*]

plot "state.dat" using (($4==0) ? (($1==0) ? ($2/5):1/0):1/0):($3/5) ls 14 title "Dipoles, Real Cascade",\
"state.dat" using (($4==0) ? (($1==1) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0) ? (($1==2) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0) ? (($1==3) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0) ? (($1==4) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0) ? (($1==5) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0) ? (($1==6) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0) ? (($1==7) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0) ? (($1==8) ? ($2/5):1/0):1/0):($3/5) ls 14 notitle,\
"state.dat" using (($4==0 && $8==1) ? (($2/5)):1/0):($3/5) ls 11 notitle,\
"state.dat" using (($4==3) ? ($2/5):1/0):($3/5) ls 12 notitle,\
"state.dat" using (($4==0) ? (($5==1) ? ($2/5):1/0):1/0):($3/5) with points ls 10 title "Partons, Real Cascade",\
"state.dat" using (($4==0) ? (($5==0) ? ($2/5):1/0):1/0):($3/5) with points ls 2 notitle,\
"state.dat" using (($4==0 && $7==1) ? (($2/5)):1/0):($3/5) ls 13 notitle,\
"state.dat" using ((($4==1) || ($4==2)) ? (($4==1) ? ($2/5):(($2/5)+$5*pTscale)):1/0):(($4==1) ? ($3/5):(($3/5)+$6*pTscale)) ls 10 notitle,\
"state.dat" using (($4==2 && $11==1) ? ($2/5):1/0):($3/5) with points pt 6 lt 1 ps 4 title "Valence Partons"

# show interacting dipoles
#"state.dat" using (($4==0 && $7==1) ? (($2/5)):1/0):($3/5) ls 13 notitle,\

# shows the PT
# "state.dat" using ((($4==1) || ($4==2)) ? (($4==1) ? ($2/5):(($2/5-1.28/2)+$5*pTscale)):1/0):(($4==1) ? ($3/5):(($3/5)+$6*pTscale)) ls 10 notitle

#marking the real gluons
#"state.dat" using (($4==0) ? (($6==1) ? ($2/5):1/0):1/0):($3/5) with points pt 6 lt 1 ps 3 notitle

#marks non-participants
#"state.dat" using (($4==3) ? ($2/5):1/0):($3/5) ls 12 notitle,\