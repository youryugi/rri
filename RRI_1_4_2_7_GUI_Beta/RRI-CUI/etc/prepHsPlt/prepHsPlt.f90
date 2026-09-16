program prepHsPlt
implicit none

character*256 infile
parameter(infile = "hs_header.plt")

character*256 outfile
parameter(outfile = "hs.plt" )

integer maxstep
parameter(maxstep = 96)

integer ios, i
character*6 cnum, cnum2
character*256 ctemp, line

open(10, file = infile, status = "old")
open(20, file = outfile)

do
 read(10, '(a)', iostat = ios) ctemp
 if(ios.ne.0) exit
 write(20, '(a)') trim(ctemp)
enddo

call int2char(maxstep, cnum2)
do i = 1, maxstep
 call int2char(i, cnum)

 line = 'set output "./hs/hs_' // trim(cnum) // '.gif"'
 write(20, '(a)') trim(line)

 line = 'splot "./out/hs_' // trim(cnum) // '.out" matrix t "' & 
        // trim(cnum) // ' / ' // trim(cnum2) // '"'
 write(20, '(a)') trim(line)

 write(20, '("")')

enddo

line = 'set output "temp.gif"'
write(20, '(a)') trim(line)

line = 'set terminal windows'
write(20, '(a)') trim(line)

end program prepHsPlt


! numbers to characters
subroutine int2char( num, cwrk )
implicit none

integer num, j
character*6 cwrk
cwrk = ' '
write( cwrk, '( I6 )' ) num
do j = 1, 6
if( cwrk( j:j ) .eq. ' ') cwrk(j:j) = '0'
enddo
end subroutine int2char