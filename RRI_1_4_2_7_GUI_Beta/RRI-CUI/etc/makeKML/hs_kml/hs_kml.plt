reset

#set terminal gif medium size 336,204 crop
set terminal gif medium size 672, 408 crop

set lmargin 0
set bmargin 0
set rmargin 0
set tmargin 0
set notics
set nokey
unset colorbox

set pm3d map

set xrange [0:]
set yrange [:] reverse
set zrange [0:] reverse

set cbrange[0:5]
set zrange[0:]

set notics

#set palette defined (0.0 "white", 1.0 "green", 5.0 "blue")
set palette defined (0.0 "gray", 1.5 "blue", 3 "green")

set output "./hs/hs_number.gif"
splot "./out/hs_number.out" matrix t "tstep / 96"

