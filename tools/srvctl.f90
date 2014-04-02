 program srvctl
!     SERVICE to GrADS format converter
!     INPUT (on commandline):
!     arg#1: name of inputfile with data in SERVICE-format
!     The program assumes that the records of the file in SERVICE format
!     are in the order produced by the PUMABURNER or AFTERBURNER program,
!     i. e. the sequence and the number of codes for each timestep must not
!     vary in the inputfile, different levels of one code at one timestep
!     must be ordered and in sequence. Usually this is not a problem.
!
!     The order is (because GrADS assumes this): level - code - time
!
!     Author: Simon Blessing (blessing(at)dkrz.de), 27.11.2000
!     Last maintenance: Oct. 27, 2004, got rid of complaint by sunf90
!                       Jan. 1st, 2005, added recunit detection
!
implicit none
type levlst
   integer :: level
   type (levlst), pointer :: next
end type levlst
type codelst
   integer :: code
   integer :: nlev
   type (codelst), pointer :: next
end type codelst
type (codelst), pointer :: cl_first, cl_help, codelist
type (levlst), pointer :: ll_first, levlist, ll_help
integer :: nlon, nlat, nlev,j, deltat,fdate,fhour,flev, nlevmax
integer :: dd,mm,mi,yyyy
integer :: ddl,mml,mil,yyyyl,lhour,ldate
integer :: tdef=1
integer :: nvar=1
logical :: l_vsect=.false.
real, allocatable :: si(:) !sinus of gaussian latitude
real, parameter :: PI     = 3.14159265359  ! Pi
character*2 :: ut
character*16 :: vform = '(a1,i3,i6,a,3i6)'
integer,parameter :: maxdpm(12) =(/ 31,29,31,30,31,30,31,31,30,31,30,31/)
character*3,parameter :: month(12) = (/'jan','feb','mar','apr','may','jun',&
                                       'jul','aug','sep','oct','nov','dec'/)
integer irec, i, narg
integer ih(8)             !header of the SERVICE file
REAL, ALLOCATABLE :: x(:) !data array
Character*255 filein, filedes !filenames
!     external iargc
!==========================================================================
      narg=iargc()
      if (narg < 1) call usage
      call getarg(1,filein)    
      if ( filein == '' .or. filein(1:1) == '-' ) call usage
      i = len_trim(filein)
      if (i > 4 .and. filein(i-3:i) == '.srv') i = i - 4
      filedes = filein(1:i) // '.ctl'
      if (narg > 1) call getarg(2,filedes)
      print*,'reading binary service file <',trim(filein),'>'
      print*,'writing grads  control file <',trim(filedes),'>'
!
!     Open the input file
!
      open(1,File=filein,FORM='UNFORMATTED',status='OLD')
!
!     Read the first header to know the length of the data array
!
      irec=0
      read (1,End=999) ih

      
      nlon = ih(5)
      nlat = ih(6)
      nlev = ih(2)               !very first guess
      if (ih(2)==-1) l_vsect=.true.
      flev    = ih(2)                !remember first level in file
      fdate   = ih(3)               !first date in dataset
      fhour   = ih(4)               !first hour in dataset
      deltat  = 1                  !first guess
      nlevmax = 1                 !first guess
      if ( NLAT .lt. 32 .or. NLAT .gt. 1000) then
         print*,ih
         print*,'Unusual number of Latitudes in header of inputfile:',  &
     &        ' "" ',filein(1:len_trim(filein)),'"" : ',NLAT
         if (nlat.gt.1000) then
            print*,'Exiting.'
            stop
         end if
      end if
!
      ALLOCATE( x(ih(5)*ih(6)) )
!
!     Open the descriptor file (formatted, sequential)
!
      Open (3,File=filedes,FORM='FORMATTED',ACCESS='SEQUENTIAL')
!
!     Calculate Gaussian latitudes
!
      if (l_vsect) then
         print*,'Presuming ',nlon,' gaussian latitudes.'
         allocate (si(nlon))     
         call inigau(nlon,si)
      else
         print*,'Presuming ',nlat,' gaussian latitudes.'
         allocate (si(nlat))     
         call inigau(nlat,si)
      end if
!
!     Write to descriptor file what we know so far
!
      write(3,'(2a)')'DSET ^',trim(filein)
!      write(3,'(2a)')'TITLE ',name
      if (kind(x)==8) then
         write(3,'(a)')'UNDEF 9e+99'
      else
         write(3,'(a)')'UNDEF 9e+09'
      end if
!      write(3,'(a)')'OPTIONS BIG_ENDIAN'
      if (l_vsect) then
         write(3,'(a,i6,a,f8.4,a,f8.4)')'XDEF ',1,                      &
     &        ' LEVELS ',1.
            write(3,'(a)')'OPTIONS YREV'
            write(3,'(a,i6,a/(8f10.4))')'YDEF ',nlon,                   &
     &           ' LEVELS ',-180./PI*asin(si(:))
      else
         write(3,'(a,i6,a,f8.4,a,f8.4)')'XDEF ',nlon,                   &
     &        ' LINEAR ',0.,' ',(360./float(nlon))
         write(3,'(a)')'OPTIONS YREV'
         if (nlat.gt.1) then
            write(3,'(a,i6,a/(8f10.4))')'YDEF ',nlat,                   &
     &           ' LEVELS ',-180./PI*asin(si(:))
         else
            write(3,'(a,i6,a,f8.4,a,f8.4)')'YDEF ',nlat,                &
     &           ' LINEAR ',0.,' ',180.
         end if
      end if
      write(3,'(A)') 'OPTIONS SEQUENTIAL'
      write(3,'(A)') 'XYHEADER 40'
!
!     Start list about codes and levels
!
      nvar=1
      allocate(cl_first)
      cl_first%code=ih(1)
      cl_first%nlev=min(ih(2),1) !special treatment of level zero
      nullify(cl_first%next)    !last in list so far; Lahey and VAST need this
      
      allocate(ll_first)
      ll_first%level=ih(2)
      nullify(ll_first%next)    !last in list so far
!
!     Do the real work (read & write data)
!     
      irec = 1                  !start record count
      do
!
!     Read the input data
!
         read (1,END=999) x
!
!     Do for the next record
!
         read (1,End=999) ih
         irec = irec+1          !count records
!
!     count time
!
         if ( ih(1)==cl_first%code.and.ih(2)==flev ) then
                                !first code first level again?
            tdef=tdef+1         !so time must have advanced! 
         end if
!
!     get first non-zero level
!
         if (ll_first%level==0 .and. ih(2).ne.0) then
            ll_first%level=ih(2)
         end if
!
!     Work on codelist
!      
         codelist=>cl_first
         do while (ih(1).ne.codelist%code)
            if (.not. associated(codelist%next)) then !add new code to list
               cl_help=>codelist
               allocate(codelist)
               cl_help%next=>codelist
               nvar=nvar+1      !count codes
               codelist%code=ih(1)
               codelist%nlev=min(ih(2),1) !special treatment of level zero
               nullify(codelist%next) !this is last in list
            else
               codelist=>codelist%next !try next code in list
            end if
         end do
!     count levels for that code
         if ( ih(2).ne.ll_first%level .and. ih(2).ne.0 .and. tdef==1 )  &
     &        then
                                !level information should be complete in
                                !first time level so we are only looking there
            codelist%nlev=codelist%nlev+1 !code repeats at first~
                          ! timelevel => must be new level and can't be zero
            levlist=>ll_first
!     maintain list of levels
            do while (levlist%level.ne.ih(2)) 
               if (.not. associated(levlist%next)) then
                  ll_help=>levlist
                  allocate(levlist)
                  ll_help%next=>levlist
                  nullify(levlist%next)
                  levlist%level=ih(2)
                  nlevmax=nlevmax+1
               else
                  levlist=>levlist%next !try next level in list
               end if
            end do
         end if
      end do
 999  Continue
      if (irec==0) then
         print*,"*** No records processed. Input file was empty! ***"
         stop
      end if
!
!     write to descriptorfile levels and codes
!
      if (l_vsect) then
         write(3,'(a,i4,a,i8,a)')'ZDEF ',nlat,' LINEAR ', &
     &        abs(ll_first%level),' 1'
      else
         if (nlevmax==1) then
            write(3,'(a,i8,a)')'ZDEF 1 LINEAR ',ll_first%level,' 1'
         else
            write(3,'(a,i6,a)')'ZDEF ',nlevmax,' LEVELS '
            levlist=>ll_first
            do while (associated(levlist))
               write(3,'(i10)')levlist%level
               levlist=>levlist%next
            end do
         endif
      end if
!     make a date
!      fdate and ih(3) are of the form yyyymmdd
!      extract first time:
      if (fhour.gt.99) then
         mi=mod(fhour,100)
         fhour=max(min(fhour/100,23),0)
      else
         fhour=max(min(fhour,23),0)
         mi=0
      end if
!      extract first date
      mm=min(12,max(mod(fdate,10000)/100,1))
      dd=min(maxdpm(mm),max(mod(fdate,100),1))
      yyyy=max(1,min((fdate-100*mm-dd)/10000,9999)) !There is no year zero
!      extract last time
      if (ih(4).gt.99) then
         mil=mod(ih(4),100)
         lhour=max(min(ih(4)/100,23),0)
      else
         lhour=max(min(ih(4),23),0)
         mil=0
      end if
!      extract last date
      ddl=min(31,max(mod(ih(3),100),1))
      mml=min(12,max(mod(ih(3),10000)/100,1))
      yyyyl=max(1,min((ih(3)-100*mml-ddl)/10000,9999)) ! There is no year zero
!      convert everything to minutes
      fdate=mi +60*(fhour+24*(dd +30*(mm +12*(yyyy )))) !  a modeler's month
      ldate=mil+60*(lhour+24*(ddl+30*(mml+12*(yyyyl)))) !  has 30 days :-)
!      calculate time increment
      deltat=max(1,(ldate-fdate)/(max(tdef-1,1))) ! just let's be robust
!      find appropriate unit for time increment
      ut='mn'                   ! fallback
      if (mod(deltat,60).eq.0) then
         deltat=deltat/60
         ut='hr'
         if (mod(deltat,24).eq.0) then
            deltat=deltat/24
            ut='dy'
            if (mod(deltat,30).eq.0) then
               deltat=deltat/30
               ut='mo'
               if (mod(deltat,12).eq.0) then
                  deltat=deltat/12
                  ut='yr'
                  if (deltat.eq.0) then
                     deltat=6
                     ut='hr'
                     print*,'Could not extract proper time increment.'
                     print*,'Assuming 6-hourly data...'
                  end if
               end if
            end if
         end if
      end if
!    Was this complicated?
            
!    We need a date string of the form hh:mmZddmmmyyyy
      write(3,'(a4,i6,a8,i2.2,a1,i2.2,a1,i2.2,a3,i4.4,i9,a2)')          &
     &     'TDEF',tdef,' LINEAR ',fhour,':',mi,'Z',dd,month(mm),        &
     &     yyyy,deltat,ut
!     write codelist
      write(3,'(a,i6)')'VARS ',nvar
      codelist=>cl_first
      do while(associated(codelist))
         write(vform(6:6),'(i1)')                                       &
     &        (min(9,1+int(log10(max(1.,real(codelist%code))))))!varying format
         if (l_vsect) codelist%nlev=nlat
         write(3,vform)'c',codelist%code,                               &
     &        codelist%nlev,' 99 ',codelist%code,ih(7),ih(8)
         if (l_vsect) then
            print*,"level -1 detected, assuming vertical section with"
            print*,nlon,"latitudes and ",nlat," linear levels."
         else
            if (codelist%nlev.ne.0 .and. codelist%nlev.ne.nlevmax) then
               print*,'***CAUTION!*** code ',codelist%code,                &
     &              ' is not defined on all levels. '
               print*,'GrADS will assume that they are the first few',     &
     &              ' from the descriptor file: ',filedes
            end if
         end if
         codelist=>codelist%next
      end do
      write(3,'(a)')'ENDVARS'
      Write (*,*) 'Processed ', irec, ' record(s)'
      End

subroutine usage
implicit none
print*,'Usage: srvctl <service file>'
print*,'Example: srvctl z500.srv will create z500.ctl'
stop
end subroutine usage


subroutine inigau(kl,pga)
implicit none
integer :: kl             ! Number of Gaussian latitudes
integer :: j0,j1          ! Loop index
real          :: pga(kl)  ! Gaussian abscissas
real (kind=8) :: pz0(kl)  ! Gaussian abscissas
real (kind=8), parameter :: pi = 3.14159265358979_8
real (kind=8) :: z0,z1,z2,z3
real (kind=8) :: ql,qld

! 1. Compute Gaussian abscissas

z0 = pi / (2*kl+1)
z1 = 1.0_8 / (kl*kl*8)

do j0 = 1 , kl/2
   z2 = z0 * (2*j0 - 0.5_8)
   z2 = cos(z2 + z1 / tan(z2))
   do j1 = 1 , 50
      z3 = ql(kl,z2) / qld(kl,z2)
      z2 = z2 - z3
      if (abs(z3) < 1.0e-16) exit
   enddo ! j1
!  print *,j1,' iter for lat= ',j0,' eps= ',z3
   pz0(j0) = z2
   pz0(kl-j0+1) = -z2
enddo ! j0
pga(:) = pz0(:) ! Double precision to native
return
end

real (kind=8) function ql(k,p)
implicit none
integer :: k
real (kind=8) :: p
real (kind=8) :: z0,z1,z2,z3,z4
integer :: j
z0 = acos(p)
z1 = 1.0
z2 = 0.0
do j = k , 0 , -2
   z3 = z1 * cos(z0 * j)
   z2 = z2 + z3
   z4 = (k-j+1) * (k+j) * 0.5
   z1 = z1 * z4 / (z4 + (j-1))
enddo ! j
if (mod(k,2) == 0) z2 = z2 - 0.5 * z3

z0 = sqrt(2.0_8)
do j = 1 ,k
   z0 = z0 * sqrt(1.0_8 - 0.25_8 / (j*j))
enddo ! j
ql = z0 * z2
return
end


real (kind=8) function qld(k,p)
implicit none
integer :: k
real (kind=8) :: p
real (kind=8) :: z0,z1,z2,z3,z4
real (kind=8) :: ql
integer :: j

if (k == 0) then
   qld = 0.0
   return
endif

z0 = 1.0 / (p*p - 1.0)
z1 = k
z2 = 2.0 * z1
z3 = sqrt((z2 + 1.0) / (z2 - 1.0))
z4 = p * ql(k,p) - z3 * ql(k-1,p)
qld = z0 * z1 * z4

return
end
