!     =====================
!     SUBROUTINE OCEAN_DUMP
!     =====================

      subroutine oce_dump
      use oceanmod
      character (len=256) :: yfn

      call step2cal(nstep,ntspd,iyea,imon,iday,ihou,imin)
      write(yfn,'("ocean_dump_",I4.4,4I2.2)') iyea,imon,iday,ihou,imin
      open(72,file=yfn)

  100 format(a16,i16)
      write (72,100) 'nfluko',nfluko
      write (72,100) 'ndiag',ndiag
      write (72,100) 'noutput',noutput
      write (72,100) 'nout',nout
      write (72,100) 'nocean',nocean
      write (72,100) 'newsurf',newsurf
      write (72,100) 'ntspd',ntspd
      write (72,100) 'nperpetual_ocean',nperpetual_ocean
      write (72,100) 'nprint',nprint
      write (72,100) 'nprhor',nprhor
      write (72,100) 'nstep',nstep
      write (72,100) 'nrestart',nrestart
      write (72,100) 'naccuout',naccuout

  200 format(a16,e20.10)
      write (72,200) 'taunc',taunc
      write (72,200) 'dlayer(1)',dlayer(1)
      write (72,200) 'vdiffk',vdiffk
      write (72,200) 'dtmix',dtmix
      write (72,200) 'solar_day',solar_day

      call oce_dump_hor('gw',gw)
      call oce_dump_hor('ysst',ysst)
      call oce_dump_hor('ymld',ymld)
      call oce_dump_hor('yls',yls)
!     call oce_dump_hor('yicec',yicec)
!     call oce_dump_hor('yiced',yiced)
!     call oce_dump_hor('yheat',yheat)
!     call oce_dump_hor('yfresh',yfresh)
!     call oce_dump_hor('ytaux',ytaux)
!     call oce_dump_hor('ytauy',ytauy)
!     call oce_dump_hor('yust3',yust3)
      call oce_dump_hor('yiflux',yiflux)
      call oce_dump_hor('ydsst',ydsst)
!     call oce_dump_hor('yclsst2',yclsst2)
!     call oce_dump_hor('yfsst2',yfsst2)
      call oce_dump_hor('yheata',yheata)
      call oce_dump_hor('yfssta',yfssta)
      call oce_dump_hor('yifluxa',yifluxa)
      call oce_dump_hor('ydssta',ydssta)

      call oce_dump_hor14('yclsst',yclsst)
      call oce_dump_hor14('yfsst',yfsst)

      close(72)

      end


      subroutine oce_dump_hor(ya,pa)
      use oceanmod

      character (len=*) :: ya
      real :: pa(NHOR)

      write(72,'(/A)') trim(ya)
      do j = 1 , NHOR , 4
      write(72,'(i6,2x,4e18.10)') j,pa(j:j+3)
      enddo
      return
      end


      subroutine oce_dump_hor14(ya,pa)
      use oceanmod

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


