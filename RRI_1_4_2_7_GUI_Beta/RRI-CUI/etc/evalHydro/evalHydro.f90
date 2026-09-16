implicit none

character*256 filefile
parameter( filefile = "evalHydro.txt" )

character*256 obsfile, calcfile
real, allocatable :: obs(:), calc(:)
integer num, itemp, ios, i, j
real eff, rmse

open(1, file = filefile, status = "old")
read(1, '(a)') obsfile
read(1, '(a)') calcfile
close(1)

open(2, file = obsfile, status = "old")
open(3, file = calcfile, status = "old")

num = 0
do
 read(2, *, iostat = ios) itemp
 if(ios .ne. 0) exit
 num = num + 1
enddo
!write(*,*) num
rewind(2)

allocate( obs(num), calc(num) )

do i = 1, num
 read(2,*) itemp, obs(i)
 read(3,*) itemp, calc(i)
enddo

call eff_calc(obs, calc, num, eff)
call rmse_calc(obs, calc, num, rmse)

write(*,*) "nash: ", eff
write(*,*) "rmse: ", rmse

pause

end


subroutine eff_calc(obs, calc, num, eff)
implicit none

integer num
real obs(num), calc(num)

real eff, f0, f, ave
integer i

f = 0
do i = 1, num
 f = f + (obs(i) - calc(i)) ** 2.
enddo

ave = 0
do i = 1, num
 ave = ave + obs(i)
enddo
ave = ave / real(num)

f0 = 0
do i = 1, num
 f0 = f0 + (obs(i) - ave) ** 2.
enddo

eff = 1. - f / f0

return
end



subroutine rmse_calc(obs, calc, num, rmse)
implicit none

integer num
real obs(num), calc(num), rmse

real eff
integer i

rmse = 0.
do i = 1, num
 rmse = rmse + (obs(i) - calc(i)) ** 2.
enddo
rmse = sqrt(rmse / real(num))

return
end
