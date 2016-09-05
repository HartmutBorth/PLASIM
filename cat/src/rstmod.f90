! **********
! * RSTMOD *
! **********
module rstmod
use catmod, only: nurstini,cat_rstini, nurstfin, nudiag

integer, parameter :: nrstdim  = 200   ! Max number of records
integer, parameter :: nrstdiag = 5     ! 1D restart arrays up
                                       ! to size nrstdiag will be
                                       ! written to cat_diag
integer            :: nexcheck =   1   ! Extended checks
integer            :: nrstnum  =   0   ! Actual number of records
integer            :: nlastrec =   0   ! Last read record
                                          
character (len=16) :: yrstnam(nrstdim) ! Array of record names

end module rstmod


! ************************
! * SUBROUTINE CHECK_RES *
! ************************
subroutine check_rst
use rstmod
implicit none

integer :: iostat,j
character (len=16) :: yn ! variable name


write(nudiag, &
 '(/," *************************************************")')
write(nudiag,'("  Checking for restart file <cat_rstini>")')
write(nudiag, &
 '(" *************************************************",/)')

do
   read (nurstini,IOSTAT=iostat) yn
   if (iostat /= 0) exit
   nrstnum = nrstnum + 1
   yrstnam(nrstnum) = yn
   read (nurstini,IOSTAT=iostat)
   if (iostat /= 0) exit
   if (nrstnum >= nrstdim) then
      write(nudiag,*) 'Too many variables in restart file'
      write(nudiag,*) 'Increase NRESDIM in module rstmod'
      write(nudiag,*) '*** Error Stop ***'
      stop
   endif
enddo

write(nudiag,'(a,i4,3a/)') 'Found ',nrstnum, &
      ' variables in file <',trim(cat_rstini),'>'
do j = 1 , nrstnum
   write(nudiag,'(i4," : ",8x,1x,a)') j,yrstnam(j)
enddo
nlastrec = nrstnum

return
end subroutine check_rst

! ********************************
! * SUBROUTINE GET_RESTART_ARRAY *
! ********************************
subroutine get_restart_array(yn,pa,k1,k2)
use rstmod
implicit none

character (len=*) :: yn
integer :: k1,k2,j,m
real(8) :: pa(k1,k2)

do j = 1 , nrstnum
   if (trim(yn) == trim(yrstnam(j))) then
      call fileseek(yn,j)
      read (nurstini) pa(:,:)
      do m = 1 , min(k1,nrstdiag)
         if (k1 .eq. 1 .and. k2 .eq. 1) then
            write(nudiag,'(" ", a," = ",e15.9)')  &
            trim(yn),pa(m,:)
         elseif (k1 .ne. 1 .and. k2 .eq. 1) then
            write(nudiag,'(" ",a, "(",i2,") = ", e15.9)')  &
            trim(yn),m, pa(m,:)
         else
            write(nudiag,'(" ",a, "(",i2,",1) = ", e15.9)')  &
            trim(yn),m, pa(m,:)
         endif
      enddo
      nlastrec = nlastrec + 1
      return
   endif
enddo
if (nexcheck == 1) then
   write(nudiag,*) '*** Error in get_restart_array ***'
   write(nudiag,*) 'Requested array {',yn,'} was not found'
   stop
endif

return
end subroutine get_restart_array


! *********************************
! * SUBROUTINE GET_RESTART_CARRAY *
! *********************************
subroutine get_restart_carray(yn,ca,k1,k2)
use rstmod
implicit none


character (len=*) :: yn
integer    :: k1,k2,j,m
complex(8) :: ca(k1,k2)

do j = 1 , nrstnum
   if (trim(yn) == trim(yrstnam(j))) then
      call fileseek(yn,j)
      read (nurstini) ca(:,:)
      do m = 1 , min(k1,nrstdiag)
         if (k1 .eq. 1 .and. k2 .eq. 1) then
            write(nudiag,'(" ",a," = ", e15.9," / ",e15.9)') &
               trim(yn), real(ca(m,1)), aimag(ca(m,1))
         elseif (k1 .ne. 1 .and. k2 .eq. 1) then
            write(nudiag,'(" ",a, "(",i2,") = ", e15.9," / ",e15.9)') &
               trim(yn),m, real(ca(m,1)), aimag(ca(m,1))
         else
            write(nudiag,'(" ",a, "(",i2,",1) = ", e15.9," / ",e15.9)') &
               trim(yn),m, real(ca(m,1)), aimag(ca(m,1))
         endif
      enddo
      nlastrec = nlastrec + 1
      return
   endif
enddo
if (nexcheck == 1) then
   write(nudiag,*) '*** Error in get_restart_array ***'
   write(nudiag,*) 'Requested array {',yn,'} was not found'
   stop
endif

return
end subroutine get_restart_carray


! *********************************
! * SUBROUTINE GET_RESTART_IARRAY *
! *********************************
subroutine get_restart_iarray(yn,ia,k1,k2)
use rstmod
implicit none


character (len=*) :: yn
integer :: k1,k2,j,m
integer :: ia(k1,k2)

do j = 1 , nrstnum

   if (trim(yn) == trim(yrstnam(j))) then
      call fileseek(yn,j)
      read (nurstini) ia(:,:)
      do m = 1 , min(k1,nrstdiag)
         if (k1 .eq. 1 .and. k2 .eq. 1) then
            write(nudiag,'(" ",a," = ", i16)') trim(yn), ia(:,:)
         elseif (k1 .ne. 1 .and. k2 .eq. 1) then
            write(nudiag,'(" ",a, "(",i2,") = ", i16)')  &
               trim(yn),m, ia(m,:)
         else
            write(nudiag,'(" ",a, "(",i2,",1) = ", i16)')  &
               trim(yn),m, ia(m,:)
         endif
      enddo
      nlastrec = nlastrec + 1
      return
   endif
enddo
if (nexcheck == 1) then
   write(nudiag,*) '*** Error in get_restart_array ***'
   write(nudiag,*) 'Requested array {',yn,'} was not found'
   stop
endif

return
end subroutine get_restart_iarray


! ********************************
! * SUBROUTINE PUT_RESTART_ARRAY *
! ********************************
subroutine put_restart_array(yn,pa,k1,k2)
use rstmod
implicit none

character (len=*)  :: yn
character (len=16) :: yy
integer :: k1,k2
real(8) :: pa(k1,k2)

yy = yn
write(nurstfin) yy
write(nurstfin) pa(1:k1,1:k2)

return
end subroutine put_restart_array


! *********************************
! * SUBROUTINE PUT_RESTART_CARRAY *
! *********************************
subroutine put_restart_carray(yn,ca,k1,k2)
use rstmod
implicit none

character (len=*)  :: yn
character (len=16) :: yy
integer    :: k1,k2
complex(8) :: ca(k1,k2)

yy = yn
write(nurstfin) yy
write(nurstfin) ca(1:k1,1:k2)

return
end subroutine put_restart_carray


! *********************************
! * SUBROUTINE PUT_RESTART_IARRAY *
! *********************************
subroutine put_restart_iarray(yn,ia,k1,k2)
use rstmod
implicit none

character (len=*)  :: yn
character (len=16) :: yy
integer :: k1,k2
integer :: ia(k1,k2)

yy = yn
write(nurstfin) yy
write(nurstfin) ia(1:k1,1:k2)

return
end subroutine put_restart_iarray


! ***********************
! * SUBROUTINE FILESEEK *
! ***********************

subroutine fileseek(yn,k)
use rstmod
implicit none

character (len=*)  :: yn
character (len=16) :: yy
integer :: iostat, k


if (k <= nlastrec) then
   rewind nurstini
   nlastrec = 0
endif

do
   read (nurstini,iostat=iostat) yy
   if (iostat /= 0) exit
   if (trim(yn) == trim(yy)) return ! success
   read (nurstini,iostat=iostat)    ! skip data
   if (iostat /= 0) exit
   nlastrec = nlastrec + 1
enddo
write(nudiag,*) 'Variable <',trim(yn),'> not in restart file'

return
end


! *****************************
! * SUBROUTINE CHECK_EQUALITY *
! *****************************
subroutine check_equality(yn,pa,pb,k1,k2)
use rstmod
implicit none

character (len=*) :: yn
integer :: k1,k2,j1,j2
real(8) :: pa(k1,k2)
real(8) :: pb(k1,k2)

do j2 = 1 , k2
   do j1 = 1 , k1
      if (pa(j1,j2) /= pb(j1,j2)) then
         write(nudiag,*) 'No Equality on ',yn,'(',j1,',',j2,')',pa(j1,j2),pb(j1,j2)
         return
      endif
   enddo
enddo
 write(nudiag,*) 'Array {',yn,'} is OK'

return
end


! **********************
! * SUBROUTINE VARSEEK *
! **********************
subroutine varseek(yn,knum)
use rstmod
implicit none

character (len=*)  :: yn
character (len=16) :: ytmp
integer :: k, knum

knum = 0
do k = 1,nrstdim
   ytmp = yrstnam(k)
   if (trim(yn) == trim(ytmp)) then 
      knum = k
   endif
enddo

return
end
