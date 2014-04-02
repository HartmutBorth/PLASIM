program f90check
integer :: i = 258
real    :: r = 1.0
open(21,file='F90_INTEGER',form='unformatted')
write (21) i
close(21)
open(21,file='F90_REAL',form='unformatted')
write (21) r
close(21)
!call csub(i,r)
stop
end program f90check
