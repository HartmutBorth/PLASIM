program pusca
implicit none

integer :: nstep
integer :: nlev
integer :: nloni,nlono
integer :: nlati,nlato
integer :: nrspi,nrspo
integer :: ntrui,ntruo

integer :: narg,io

character (len=1024) :: yarg,yfni,yfno
character (len=16)   :: yid

narg = iargc()

if (narg /= 3) then
   print *,"Usage: pusca.x <newlats> <inputfile> <outputfile>"
   stop
endif

call getarg(1,yarg)
call getarg(2,yfni)
call getarg(3,yfno)

read(yarg,*) nlato

nlono = 2 * nlato
ntruo = (nlono - 1) / 3
nrspo = (ntruo + 1) * (ntruo + 2)

open(10,file=yfni,form='unformatted',status='old',iostat=io)
if (io /= 0) then
   print *,"Error opening file <",trim(yfni),">"
   stop
endif

open(20,file=yfno,form='unformatted')

call get_restart_integer('nstep',nstep)
call get_restart_integer('nlat' ,nlati)
call get_restart_integer('nlon' ,nloni)
call get_restart_integer('nlev' ,nlev )
call get_restart_integer('nrsp' ,nrspi)

ntrui = (nloni - 1) / 3

write (*,'(/,56("*"))')
write (*,'("*        * ",a20," * ",a20," *")') trim(yfni),trim(yfno)
write (*,'(56("*"))')
write (*,'("* NLAT   * ",i20," * ",i20," *")') nlati,nlato
write (*,'("* NLON   * ",i20," * ",i20," *")') nloni,nlono
write (*,'("* NTRU   * ",i20," * ",i20," *")') ntrui,ntruo
write (*,'("* NRSP   * ",i20," * ",i20," *")') nrspi,nrspo
write (*,'("* NLEV   * ",i20," * ",i20," *")') nlev,nlev
write (*,'("* NSTEP  * ",i20," * ",i20," *")') nstep,nstep
write (*,'(56("*"))')

call put_restart_integer('nstep',nstep)
call put_restart_integer('nlat' ,nlato)
call put_restart_integer('nlon' ,nlono)
call put_restart_integer('nlev' ,nlev )
call put_restart_integer('nrsp' ,nrspo)

call convert_restart_array('sz' ,ntrui,ntruo,nlev)
call convert_restart_array('sd' ,ntrui,ntruo,nlev)
call convert_restart_array('st' ,ntrui,ntruo,nlev)
call convert_restart_array('st2',ntrui,ntruo,nlev)
call convert_restart_array('st3',ntrui,ntruo,nlev)
call convert_restart_array('sr1',ntrui,ntruo,nlev)
call convert_restart_array('sr2',ntrui,ntruo,nlev)
call convert_restart_array('sp' ,ntrui,ntruo,   1)   
call convert_restart_array('sp2',ntrui,ntruo,   1)   
call convert_restart_array('sp3',ntrui,ntruo,   1)   
call convert_restart_array('so' ,ntrui,ntruo,   1)   
call convert_restart_array('szm',ntrui,ntruo,nlev)
call convert_restart_array('sdm',ntrui,ntruo,nlev)
call convert_restart_array('stm',ntrui,ntruo,nlev)
call convert_restart_array('spm',ntrui,ntruo,   1)   

close (10)
close (20)
stop
end


subroutine convert_restart_array(yn,ntrui,ntruo,nlev)
character (len=*)  :: yn
character (len=16) :: yi

real, allocatable :: fi(:,:)
real, allocatable :: fo(:,:)
real, allocatable :: fx(:,:)
real, allocatable :: fy(:,:)

ntrux = max(ntrui,ntruo) + 1
nrspi = (ntrui + 1) * (ntrui + 2)
nrspo = (ntruo + 1) * (ntruo + 2)

allocate(fi(nrspi,nlev))
allocate(fo(nrspo,nlev))
allocate(fx(ntrux,ntrux))
allocate(fy(ntrux,ntrux))

read (10,iostat=io) yi
if (io /= 0) then
   print *,"I/O Error reading ",trim(yn)
   stop
endif
if (trim(yi) /= trim(yn)) then
   print *,"Looking for: <",trim(yn),">  found: <",trim(yi),">"
   stop
endif
read (10) fi(1:nrspi,1:nlev)

do jlev = 1 , NLEV
   fx(:,:) = 0.0
   fy(:,:) = 0.0
   k = 0
   do m = 1 , ntrui+1
      do n = m , ntrui+1
         fx(n,m) = fi(k+1,jlev)
         fy(n,m) = fi(k+2,jlev)
         k = k + 2
      enddo
   enddo
   k = 0
   do m = 1 , ntruo+1
      do n = m , ntruo+1
         fo(k+1,jlev) = fx(n,m)
         fo(k+2,jlev) = fy(n,m)
         k = k + 2
      enddo
   enddo
enddo

write (20) yi
write (20) fo

print *,"Converted array: ",trim(yi)

deallocate(fi)
deallocate(fo)
deallocate(fx)
deallocate(fy)

return
end


! ==============================
! SUBROUTINE GET_RESTART_INTEGER
! ==============================

subroutine get_restart_integer(yn,kv)

character (len=* ) :: yn
character (len=16) :: yi
integer :: kv

read (10) yi
read (10) kv
if (trim(yi) /= trim(yn)) then
   print *,'*** Error in get_restart_integer ***'
   print *,'Requested integer {',yn,'} was not found'
   stop
endif
return
end subroutine get_restart_integer


! ==============================
! SUBROUTINE PUT_RESTART_INTEGER
! ==============================

subroutine put_restart_integer(yn,kv)

character (len=* ) :: yn
character (len=16) :: yi
integer :: kv

yi = yn

write (20) yi
write (20) kv
return
end subroutine put_restart_integer



