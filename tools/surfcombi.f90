module splitmod
implicit none
integer :: nlon = 0
integer :: nlat,ntru
integer :: icemon = 0
integer :: sstmon = 0
integer :: nrec   = 0
character(len=1024) :: ofile,idir,yres
integer :: head(8) = 0
real, parameter :: GA = 9.81 ! g
real, parameter :: TMELT = 273.16 
real, allocatable :: f(:,:)
character(len=1), allocatable :: dls(:,:)
character(len=1), allocatable :: dglac(:,:)
character(len=1), allocatable :: dice(:,:)
real,  allocatable :: doro(:,:)

end module splitmod

program surfcombi
use splitmod
implicit none

interface
   integer function iargc()
   end function iargc
end interface

integer :: nargc,iostat

nargc = iargc()

if (nargc /= 3) then
   print '("Usage surfcombi <t-res>  <input directory> <surface file>")'
   stop
endif

call getarg(1,yres)
call getarg(2,idir)
call getarg(3,ofile)
read(yres,*) ntru

open(10,file=ofile)

print '(48("*"))'
print '("* Surface file generator for Planet Simulator *")'
print '(48("*"))'

call submask ! Read combined Land/Sea - Glacier mask

if (nlon == 0) stop
head(5) = nlon
head(6) = nlat

allocate(f(nlon,nlat))
allocate(dls(nlon,nlat))
allocate(dice(nlon,nlat))
allocate(doro(nlon,nlat))

! 1.0 Process binary masks

call subdls   ! Land/Sea mask [172]
call subdglac ! Glacier mask  [232]

! 2.0 Process constant arrays

call subdoro  ! Orography [129]
call subcommon(173 ,'z0')
call subcommon(1730,'z0noveg')
call subcommon(174 ,'albedo')
call subcommon(209 ,'soil')

! 3.0 Process arrays with annual cycle (1-12) or (0-13) months

call subyear(169,'sst')
call subyear(210,'ice')

stop
end

subroutine subdls
use splitmod
implicit none
integer :: jlon,jlat
character(len=80) yform

head(1) = 172
f(:,:) = 0
where (dglac(:,:) == '*' .or. dglac(:,:) == '=') f(:,:) = 1.0
write(10,'(8i10)') head
write(10,'(8e12.6)') f(:,:)

nrec = nrec + 1
write(*,'("Record",i3,"   Code",i5,3x,a)') &
          nrec,head(1),'dls'
return
end

subroutine subdglac
use splitmod
implicit none
integer :: jlon,jlat
character(len=80) yform

head(1) = 232
f(:,:) = 0
where (dglac(:,:) == '=') f(:,:) = 1.0
write(10,'(8i10)') head
write(10,'(8e12.6)') f(:,:)

nrec = nrec + 1
write(*,'("Record",i3,"   Code",i5,3x,a)') &
          nrec,head(1),'dglac'
return
end

subroutine subdoro
use splitmod
implicit none
integer :: jlon,jlat,ilat
character (len=1024) :: ifile
character (len=80) :: yline
character (len=80) :: yform = ' '
integer :: itopo(nlon,nlat)

head(1) = 129

write (ifile,'(a,"/t",i3.3,"_topo.txt")') trim(idir),ntru

! write(*,'("Open file <",a,">")') trim(ifile)
open(11,file=ifile)
do
read(11,'(a)') yline
if (trim(yline) == '[BEGIN ARRAY]') exit
if (yline(1:4) == 'NLAT') then
   read(yline(10:),*) ilat
   if (ilat /= nlat) then
      write(*,'("*** Resolution mismatch ***")')
      write(*,'("<mask> NLAT    = ",I4)') nlat
      write(*,'("<topo> NLAT    = ",I4)') ilat
      stop 
   endif
endif
if (yline(1:9) == 'FORMAT  =') then
   read(yline(11:),'(a)') yform
endif
enddo
read (11,yform) itopo(:,:)
close(11)

f(:,:) = GA * itopo(:,:)

write(10,'(8i10)') head
write(10,'(8e12.6)') f(:,:)

nrec = nrec + 1
write(*,'("Record",i3,"   Code",i5,3x,a)') &
          nrec,head(1),'dls'
return
end

subroutine subcommon(kcode,ktext)
use splitmod
implicit none
integer :: kcode
integer :: ilat
character (len=*) :: ktext
character (len=1024) :: ifile
character (len=80) :: yline
character (len=80) :: yform = ' '
character (len=80) :: yarray = ' '
logical :: lc2k

lc2k = .false.
head(1) = kcode

write (ifile,'(a,"/t",i3.3,"_",a,".txt")') trim(idir),ntru,trim(ktext)

! write(*,'("Open file <",a,">")') trim(ifile)
open(11,file=ifile)
do
read(11,'(a)') yline
if (trim(yline) == '[BEGIN ARRAY]') exit
if (yline(1:4) == 'NLAT') then
   read(yline(10:),*) ilat
   if (ilat /= nlat) then
      write(*,'("*** Resolution mismatch ***")')
      write(*,'("<mask> NLAT    = ",I4)') nlat
      write(*,'("<",a,"> NLAT    = ",I4)') trim(ktext),ilat
      stop 
   endif
endif
if (yline(1:9) == 'FORMAT  =') then
   read(yline(11:),'(a)') yform
endif

! Check units
  
if (yline(1:4) == 'UNIT') then
   if (trim(yline(11:)) == 'Celsius') then
      lc2k = .true.
   endif
endif
   
! Read array name
  
if (yline(1:5) == 'ARRAY') then
   yarray = trim(yline(11:))
endif
   
enddo
read (11,yform) f(:,:)
close(11)

if (lc2k) f(:,:) = f(:,:) + TMELT

write(10,'(8i10)') head
write(10,'(8e12.6)') f(:,:)

nrec = nrec + 1
write(*,'("Record",i3,"   Code",i5,3x,a)') &
          nrec,kcode,trim(yarray)
return
end

subroutine subyear(kcode,ktext)
use splitmod
implicit none
integer :: kcode
integer :: ilat,icode,imonth,ios
character (len=*) :: ktext
character (len=1024) :: ifile
character (len=80) :: yline
character (len=80) :: yform
character (len=1)  :: yice(nlon,nlat)
character (len=80) :: yarray
logical :: lc2k ! Celsius -> Kelvin flag
logical :: leof ! EOF flag

lc2k = .false.
leof = .false.
yform = ' '
yarray = ' '

head(1) = kcode

write (ifile,'(a,"/t",i3.3,"_",a,".txt")') trim(idir),ntru,trim(ktext)

! write(*,'("Open file <",a,">")') trim(ifile)
open(11,file=ifile,iostat=ios)
if (ios /= 0) then
   write(*,'("*** Error IOSTAT=",i5," opening <",a,"> ***")') ios,trim(ifile)
   stop 
endif


! Read first line

read(11,'(a)',iostat=ios) yline
if (ios /= 0) then
   write(*,'("*** Could not read from file <",a,">")') trim(ifile)
   stop 
endif

! loop until end-of-file

do while (.not. leof)

   ! Check header

   if (yline(1:22) /= '[PLASIM SURFACE ARRAY]') then
      write(*,'("*** No [PLASIM SURFACE ARRAY] found ***")')
      stop 
   endif

   do
   read(11,'(a)') yline
   if (trim(yline) == '[BEGIN ARRAY]') exit ! Data section
   
   ! Check code
   
   if (yline(1:4) == 'CODE') then
      read(yline(10:),*) icode
      if (icode /= kcode) then
         write(*,'("*** Code mismatch ***")')
         write(*,'("expected CODE    = ",I4)') kcode
         write(*,'("found    CODE    = ",I4)') icode
         stop 
      endif
   endif
   
   ! Check resolution
   
   if (yline(1:4) == 'NLAT') then
      read(yline(10:),*) ilat
      if (ilat /= nlat) then
         write(*,'("*** Resolution mismatch ***")')
         write(*,'("expected NLAT    = ",I4)') nlat
         write(*,'("found    NLAT    = ",I4)') ilat
         stop 
      endif
   endif
   
   ! Check units
   
   if (yline(1:4) == 'UNIT') then
      if (trim(yline(11:)) == 'Celsius') then
         lc2k = .true.
      endif
   endif
   
   ! Read array name
   
   if (yline(1:5) == 'ARRAY') then
      yarray = trim(yline(11:))
   endif
   
   ! Set header
   
   if (yline(1:5) == 'MONTH') then
      read(yline(10:),*) imonth
      if (imonth < 0 .or. imonth > 13) then
         write(*,'("*** Illegal value for month [0-13] ***")')
         write(*,'("MONTH   = ",I10)') imonth
         stop 
      endif
      head(3) = 990000 + imonth * 100
      head(4) = -1
   endif
   
   ! Read format
   
   if (yline(1:9) == 'FORMAT  =') then
      read(yline(11:),'(a)') yform
   endif
   enddo

   if (kcode == 210) then ! sea ice mask
      read (11,yform) yice(:,:)
   !  write(*,yform) yice
      f(:,:) = 0.0
      where (yice(:,:) == '.') f(:,:) = 1.0
   else
      read (11,yform) f(:,:)
      if (lc2k) f(:,:) = f(:,:) + TMELT
   endif
   
   write(10,'(8i10)') head
   write(10,'(8e12.6)') f(:,:)
   nrec = nrec + 1
   write(*,'("Record",i3,"   Code",i5,"   Month",i3,3x,a)') &
             nrec,kcode,imonth,trim(yarray)

   read(11,'(a)') yline

   ! Check end of array marker

   if (yline(1:11) /= '[END ARRAY]') then
      write(*,'("*** No [END ARRAY] found ***")')
      stop 
   endif

   read(11,'(a)',iostat=ios) yline
   leof = ios < 0

enddo ! while (.not. leof)

close(11)

return
end

subroutine submask
use splitmod
implicit none
integer :: jlon,jlat
character(len=80) :: yform
character(len=80) :: yline

character (len=1024) :: ifile
write (ifile,'(a,"/t",i3.3,"_mask.txt")') trim(idir),ntru

! write(*,'("Open file <",a,">")') trim(ifile)
open(11,file=ifile)
do
read(11,'(a)') yline
if (trim(yline) == '[BEGIN ARRAY]') exit
if (yline(1:4) == 'NLON') then
   read(yline(10:),*) nlon
   write(*,'("NLON    = ",I4)') nlon
endif
if (yline(1:4) == 'NLAT') then
   read(yline(10:),*) nlat
   write(*,'("NLAT    = ",I4)') nlat
endif
enddo
allocate(dglac(nlon,nlat))
write (yform,'("(",i,"a1)")') nlon
read (11,yform) dglac(:,:)
write (*,yform) dglac(:,:)
close(11)
return
end

