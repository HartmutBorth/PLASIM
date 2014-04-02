!     ===================
!     SUBROUTINE ATM_DUMP
!     ===================

      subroutine atm_dump
      use pumamod
      character (len=256) :: yfn

      call step2cal(nstep,ntspd,iyea,imon,iday,ihou,imin)
      write(yfn,'("atm_dump_",I4.4,4I2.2)') iyea,imon,iday,ihou,imin
      open(72,file=yfn)

  100 format(a16,i16)
      write (72,100) 'nice',nice
      write (72,100) 'mocd',mocd
      write (72,100) 'n_start_year',n_start_year
      write (72,100) 'n_start_month',n_start_month
      write (72,100) 'ntspd',ntspd
      write (72,100) 'noutput',noutput
      write (72,100) 'nout',nout
      write (72,100) 'n_start_step',n_start_step
      write (72,100) 'n_days_per_month',n_days_per_month
      write (72,100) 'n_days_per_year',n_days_per_year
      write (72,100) 'nprint',nprint
      write (72,100) 'nprhor',nprhor
      write (72,100) 'nstep',nstep
      write (72,100) 'nrestart',nrestart
      write (72,100) 'naccuout',naccuout
      write (72,100) 'mpstep',mpstep
      write (72,100) 'nwpd',nwpd
      write (72,100) 'kick',kick
      write (72,100) 'nafter',nafter
      write (72,100) 'ncoeff',ncoeff
      write (72,100) 'ndiag',ndiag
      write (72,100) 'ndebug',ndebug
      write (72,100) 'nexp',nexp
      write (72,100) 'nexper',nexper
      write (72,100) 'ngui',ngui
      write (72,100) 'npacksp',npacksp
      write (72,100) 'npackgp',npackgp
      write (72,100) 'ndiaggp',ndiaggp
      write (72,100) 'ndiagsp',ndiagsp
      write (72,100) 'ndiagcf',ndiagcf
      write (72,100) 'ndiaggp2d',ndiaggp2d
      write (72,100) 'ndiaggp3d',ndiaggp3d
      write (72,100) 'ndiagsp2d',ndiagsp2d
      write (72,100) 'ndiagsp3d',ndiagsp3d
      write (72,100) 'nhdiff',nhdiff
      write (72,100) 'ntime',ntime
      write (72,100) 'nperpetual',nperpetual
      write (72,100) 'n_sea_points',n_sea_points
      write (72,100) 'mars',mars

  200 format(a16,e20.10)
      write (72,200) 'dawn',dawn
      write (72,200) 'delt',delt
      write (72,200) 'delt2',delt2
      write (72,200) 'dtep',dtep
      write (72,200) 'dtns',dtns
      write (72,200) 'dtrop',dtrop
      write (72,200) 'dttrp',dttrp
      write (72,200) 'tgr',tgr
      write (72,200) 'psurf',psurf
      write (72,200) 'time0',time0
      write (72,200) 'co2',co2
      write (72,200) 'umax',umax
      write (72,200) 't2mean',t2mean
      write (72,200) 'precip',precip
      write (72,200) 'evap',evap
      write (72,200) 'olr',olr
      write (72,200) 'akap',akap
      write (72,200) 'alr',alr
      write (72,200) 'ga',ga
      write (72,200) 'gascon',gascon
      write (72,200) 'plarad',plarad
      write (72,200) 'pnu',pnu
      write (72,200) 'siderial_day',siderial_day
      write (72,200) 'solar_day',solar_day
      write (72,200) 'ww',ww
      write (72,200) 'ra1',ra1
      write (72,200) 'ra2',ra2
      write (72,200) 'ra4',ra4
      write (72,200) 'acpd',acpd
      write (72,200) 'adv',adv
      write (72,200) 'cv',cv
      write (72,200) 'ct',ct
      write (72,200) 'pnu21',pnu21
      write (72,200) 'rdbrv',rdbrv

      close(72)

      end


      subroutine atm_dump_hor(ya,pa)
      use pumamod

      character (len=*) :: ya
      real :: pa(NHOR)

      write(72,'(/A)') trim(ya)
      do j = 1 , NHOR , 4
      write(72,'(i6,2x,4e18.10)') j,pa(j:j+3)
      enddo
      return
      end


      subroutine atm_dump_hor14(ya,pa)
      use pumamod

      character (len=*) :: ya
      real :: pa(NHOR,14)


      do mon = 1 , 14
         write(72,'(/I2,2X,A)') mon-1,trim(ya)
         do j = 1 , NHOR , 4
         write(72,'(i6,2x,4e18.10)') j,pa(j:j+3,mon)
         enddo
      enddo
      return
      end


