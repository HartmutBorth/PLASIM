      module restartmod
      integer, parameter :: nresdim  = 200   ! Max number of records
      integer, parameter :: nreaunit =  33   ! FORTRAN unit for reading
      integer, parameter :: nwriunit =  34   ! FORTRAN unit for writing
      integer, parameter :: ndebunit =  35   ! FORTRAN unit for debugging
      integer            :: nexcheck =   1   ! Extended checks
      integer            :: nresnum  =   0   ! Actual number of records
      integer            :: nlastrec =   0   ! Last read record
      integer            :: nud      =   6   ! Standard output
      integer            :: nresdbg  =   0   ! Debug mode
      character (len=16) :: yresnam(nresdim) ! Array of record names
      end module restartmod

!     ======================
!     SUBROUTINE RESTART_INI
!     ======================

      subroutine restart_ini(lrestart,yrfile)
      use restartmod

      logical :: lrestart
      character (len=*)  :: yrfile
      character (len=16) :: yn ! variable name

      inquire(file=yrfile,exist=lrestart)
      if (lrestart) then
         open(nreaunit,file=yrfile,form='unformatted')
         do
            read (nreaunit,IOSTAT=iostat) yn
            if (iostat /= 0) exit
            nresnum = nresnum + 1
            yresnam(nresnum) = yn
            read (nreaunit,IOSTAT=iostat)
            if (iostat /= 0) exit
            if (nresnum >= nresdim) then
               write(nud,*) 'Too many variables in restart file'
               write(nud,*) 'Increase NRESDIM in module restartmod'
               write(nud,*) '*** Error Stop ***'
               stop
            endif
         enddo
   
         write(nud,'(a,i4,3a/)') 'Found ',nresnum, &
               ' variables in file <',trim(yrfile),'>'
         do j = 1 , nresnum
            write(nud,'(i4," : ",8x,1x,a)') j,yresnam(j)
         enddo
         nlastrec = nresnum
      endif ! (lrestart)

!     file must be left open for further access

      return
      end subroutine restart_ini     


!     ==========================
!     SUBROUTINE RESTART_PREPARE
!     ==========================

      subroutine restart_prepare(ywfile)
      use restartmod

      character (len=*) :: ywfile

      open(nwriunit,file=ywfile,form='unformatted')
      if (nresdbg > 0) then
         open(ndebunit,file='status.txt',form='formatted')
      endif

      return
      end subroutine restart_prepare


!     =======================
!     SUBROUTINE RESTART_STOP
!     =======================

      subroutine restart_stop
      use restartmod

      close (nreaunit)
      close (nwriunit)

      return
      end subroutine restart_stop


!     ==============================
!     SUBROUTINE GET_RESTART_INTEGER
!     ==============================

      subroutine get_restart_integer(yn,kv)
      use restartmod

      character (len=*) :: yn
      integer :: kv

      do j = 1 , nresnum
         if (trim(yn) == trim(yresnam(j))) then
            call fileseek(yn,j)
            read (nreaunit) kv
            nlastrec = nlastrec + 1
            return
         endif
      enddo
      if (nexcheck == 1) then
         write(nud,*) '*** Error in get_restart_integer ***'
         write(nud,*) 'Requested integer {',yn,'} was not found'
         stop
      endif
      return
      end subroutine get_restart_integer


!     ============================
!     SUBROUTINE GET_RESTART_ARRAY
!     ============================

      subroutine get_restart_array(yn,pa,k1,k2,k3)
      use restartmod

      character (len=*) :: yn
      real :: pa(k2,k3)

      do j = 1 , nresnum
         if (trim(yn) == trim(yresnam(j))) then
            call fileseek(yn,j)
            read (nreaunit) pa(1:k1,:)
            nlastrec = nlastrec + 1
            return
         endif
      enddo
      if (nexcheck == 1) then
         write(nud,*) '*** Error in get_restart_array ***'
         write(nud,*) 'Requested array {',yn,'} was not found'
         stop
      endif
      return
      end subroutine get_restart_array


!     ==============================
!     SUBROUTINE PUT_RESTART_INTEGER
!     ==============================

      subroutine put_restart_integer(yn,kv)
      use restartmod

      character (len=*)  :: yn
      character (len=16) :: yy
      integer :: kv

      yy = yn
      write(nwriunit) yy
      write(nwriunit) kv
      if (nresdbg > 0) then
         write(ndebunit,'(a,a,i10)') yy," = ",kv
      endif
      return
      end subroutine put_restart_integer


!     ============================
!     SUBROUTINE PUT_RESTART_ARRAY
!     ============================

      subroutine put_restart_array(yn,pa,k1,k2,k3)
      use restartmod

      character (len=*)  :: yn
      character (len=16) :: yy
      integer :: k1,k2,k3
      real :: pa(k2,k3)

      yy = yn
      write(nwriunit) yy
      write(nwriunit) pa(1:k1,1:k3)
      if (nresdbg > 0) then
         do jl = 1 , k3
         do jm = 1 , k1 , 6
            if (jm < k2 - 5) then
               write(ndebunit,'(a,"(",i4,",",i2,") = ",6e11.4)') &
               yy,jm,jl,pa(jm:jm+5,jl)
            endif
         enddo
         enddo
      endif
      return
      end subroutine put_restart_array


!     ===================
!     SUBROUTINE FILESEEK
!     ===================

      subroutine fileseek(yn,k)
      use restartmod

      character (len=*)  :: yn
      character (len=16) :: yy

!     write(nud,*) 'Pos:',nlastrec,'   Want:',k
      if (k <= nlastrec) then
!        write(nud,*) 'Rewinding'
         rewind nreaunit
         nlastrec = 0
      endif

      do
         read (nreaunit,iostat=iostat) yy
         if (iostat /= 0) exit
         if (trim(yn) == trim(yy)) return ! success
         read (nreaunit,iostat=iostat)    ! skip data
         if (iostat /= 0) exit
         nlastrec = nlastrec + 1
      enddo
      write(nud,*) 'Variable <',trim(yn),'> not in restart file'
      return
      end


!     =========================
!     SUBROUTINE CHECK_EQUALITY
!     =========================

      subroutine check_equality(yn,pa,pb,k1,k2)
      character (len=*) :: yn
      real :: pa(k1,k2)
      real :: pb(k1,k2)

      do j2 = 1 , k2
      do j1 = 1 , k1
         if (pa(j1,j2) /= pb(j1,j2)) then
            write(nud,*) 'No Equality on ',yn,'(',j1,',',j2,')',pa(j1,j2),pb(j1,j2)
            return
         endif
      enddo
      enddo
      write(nud,*) 'Array {',yn,'} is OK'
      return
      end


!     ==================
!     SUBROUTINE VARSEEK
!     ==================

      subroutine varseek(yn,knum)
      use restartmod

      character (len=*)  :: yn
      character (len=16) :: ytmp
      integer :: k, knum

      knum = 0
      do k = 1,nresdim
         ytmp = yresnam(k)
         if (trim(yn) == trim(ytmp)) then 
            knum = k
         endif
      enddo
      return
      end
