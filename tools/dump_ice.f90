!     ===================
!     SUBROUTINE ICE_DUMP
!     ===================

      subroutine ice_dump
      use icemod
      character (len=256) :: yfn

      call step2cal(nstep,ntspd,iyea,imon,iday,ihou,imin)
      write(yfn,'("ice_dump_",I4.4,4I2.2)') iyea,imon,iday,ihou,imin
      open(72,file=yfn)

  100 format(a16,i16)
      write (72,100) 'nice',nice
      write (72,100) 'newsurf',newsurf
      write (72,100) 'nsnow',nsnow
      write (72,100) 'ntskin',ntskin
      write (72,100) 'ntspd',ntspd
      write (72,100) 'noutput',noutput
      write (72,100) 'nout',nout
      write (72,100) 'nfluko',nfluko
      write (72,100) 'ncpl_ice_ocean',ncpl_ice_ocean
      write (72,100) 'nperpetual_ice',nperpetual_ice
      write (72,100) 'nprint',nprint
      write (72,100) 'nprhor',nprhor
      write (72,100) 'nstep',nstep
      write (72,100) 'nrestart',nrestart
      write (72,100) 'naccuout',naccuout

  200 format(a16,e20.10)
      write (72,200) 'taunc',taunc
      write (72,200) 'xmind',xmind
      write (72,200) 'xdt',xdt
      write (72,200) 'cicemin',cicemin
      write (72,200) 'solar_day',solar_day

!     call ice_dump_hor('xshfl',xshfl)
!     call ice_dump_hor('xshdt',xshdt)
!     call ice_dump_hor('xlhfl',xlhfl)
!     call ice_dump_hor('xlhdt',xlhdt)
!     call ice_dump_hor('xswfl',xswfl)
!     call ice_dump_hor('xlwfl',xlwfl)
      call ice_dump_hor('xts',xts)
      call ice_dump_hor('xsst',xsst)
      call ice_dump_hor('xmld',xmld)
      call ice_dump_hor('xls',xls)
      call ice_dump_hor('xiced',xiced)
      call ice_dump_hor('xicec',xicec)
      call ice_dump_hor('xsnow',xsnow)
!     call ice_dump_hor('xsmelt',xsmelt)
!     call ice_dump_hor('xsmflx',xsmflx)
!     call ice_dump_hor('ximelt',ximelt)
!     call ice_dump_hor('xsndch',xsndch)
      call ice_dump_hor('xqmelt',xqmelt)
!     call ice_dump_hor('xheat',xheat)
!     call ice_dump_hor('xcflux',xcflux)
!     call ice_dump_hor('xprs',xprs)
!     call ice_dump_hor('xfresh',xfresh)
      call ice_dump_hor('xroff',xroff)
!     call ice_dump_hor('xtaux',xtaux)
!     call ice_dump_hor('xtauy',xtauy)
!     call ice_dump_hor('xust3',xust3)
      call ice_dump_hor('xoflux',xoflux)
!     call ice_dump_hor('xtsflux',xtsflux)
!     call ice_dump_hor('xhnew',xhnew)
      call ice_dump_hor('cheat',cheat)
      call ice_dump_hor('cfresh',cfresh)
      call ice_dump_hor('ctaux',ctaux)
      call ice_dump_hor('ctauy',ctauy)
      call ice_dump_hor('cust3',cust3)
      call ice_dump_hor('xflxicea',xflxicea)
      call ice_dump_hor('xheata',xheata)
      call ice_dump_hor('xofluxa',xofluxa)
      call ice_dump_hor('xqmelta',xqmelta)
      call ice_dump_hor('xcfluxa',xcfluxa)
      call ice_dump_hor('xsmelta',xsmelta)
      call ice_dump_hor('ximelta',ximelta)
      call ice_dump_hor('xtsfluxa',xtsfluxa)
      call ice_dump_hor('xhnewa',xhnewa)

      call ice_dump_hor14('xclsst',xclsst)
      call ice_dump_hor14('xclicec',xclicec)
      call ice_dump_hor14('xcliced',xcliced)
      call ice_dump_hor14('xflxice',xflxice)

      close(72)

      end


      subroutine ice_dump_hor(ya,pa)
      use icemod

      character (len=*) :: ya
      real :: pa(NHOR)

      write(72,'(/A)') trim(ya)
      do j = 1 , NHOR , 4
      write(72,'(i6,2x,4e18.10)') j,pa(j:j+3)
      enddo
      return
      end


      subroutine ice_dump_hor14(ya,pa)
      use icemod

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


