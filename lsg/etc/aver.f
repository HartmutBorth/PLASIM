!     ==================================================================
      subroutine aver(knt,ksum)
      use lsgvar
      implicit none
!
!     ------------------------------------------------------------------
!
!**** *aver* averages data sets.
!     When mod(knt,ksum) equals zero, the averaging is performed.
!
!     by U. Mikolajewicz 1/89.
!     optimzed by E. Kirk 2/2008
!
!*    Input.
!     ------
!     *knt*     actual time step (via parameterlist).
!     *ksum*    number of data sets to be averaged (via parameterlist).
!     Various data sets for 1 time step via module lsgvar
!
!*    Output.
!     -------
!     averaged data sets via module lsgvar
!
!*    Interface.
!     ----------
!     *call* *aver(knt,ksum)*
!
!     ------------------------------------------------------------------
      integer, intent(in) :: knt
      integer, intent(in) :: ksum
      integer, save       :: ncount = 0
      real (kind=8) :: fakt

      real (kind=8),save :: tmea(ienjen,ken)
      real (kind=8),save :: smea(ienjen,ken)
      real (kind=8),save :: umea(ienjen,ken)
      real (kind=8),save :: vmea(ienjen,ken)
      real (kind=8),save :: wmea(ienjen,ken)
      real (kind=8),save :: cmea(ienjen,ken)
      real (kind=8),save :: zmea(ienjen),simea(ienjen),fhmea(ienjen)
      real (kind=8),save :: fwmea(ienjen),ubmea(ienjen),vbmea(ienjen)
      real (kind=8),save :: hfmea(ienjen),psimea(ienjen),flwmea(ienjen)
      real (kind=8),save :: txmea(ienjen),tymea(ienjen),tbmea(ienjen)


!
!*    1.        Copy or sum up
!     ---------------
!
!     ncount counts the number of data sets summed up.
!
      if (ncount == 0) then
         tmea(:,:) = ta(:,:)
         smea(:,:) = sa(:,:)
         umea(:,:) = utota(:,:)
         vmea(:,:) = vtota(:,:)
         wmea(:,:) = wa(:,:)
         cmea(:,:) = convaa(:,:)
         flwmea(:) = flwa(:)
         psimea(:) = psia(:)
         simea(:)  = sicea(:)
         ubmea(:)  = uba(:)
         vbmea(:)  = vba(:)
         fhmea(:)  = flha(:)
         fwmea(:)  = fwa(:)
         txmea(:)  = tauxa(:)
         tymea(:)  = tauya(:)
         tbmea(:)  = tba(:)
         hfmea(:)  = hfa(:)
         zmea(:)   = zetaa(:)
      else
         tmea(:,:) = tmea(:,:) + ta(:,:)
         smea(:,:) = smea(:,:) + sa(:,:)
         umea(:,:) = umea(:,:) + utota(:,:)
         vmea(:,:) = vmea(:,:) + vtota(:,:)
         wmea(:,:) = wmea(:,:) + wa(:,:)
         cmea(:,:) = cmea(:,:) + convaa(:,:)
         flwmea(:) = flwmea(:) + flwa(:)
         psimea(:) = psimea(:) + psia(:)
         simea(:)  = simea(:)  + sicea(:)
         ubmea(:)  = ubmea(:)  + uba(:)
         vbmea(:)  = vbmea(:)  + vba(:)
         fhmea(:)  = fhmea(:)  + flha(:)
         fwmea(:)  = fwmea(:)  + fwa(:)
         txmea(:)  = txmea(:)  + tauxa(:)
         tymea(:)  = tymea(:)  + tauya(:)
         tbmea(:)  = tbmea(:)  + tba(:)
         hfmea(:)  = hfmea(:)  + hfa(:)
         zmea(:)   = zmea(:)   + zetaa(:)
      end if
      ncount=ncount+1
!
!*    2.        Final processing.
!     -----------------
      if (mod(knt,ksum)/=0) return

      if (ncount<ksum-2) then
         print *," AVER: Warning !! average from begin of restart"
         print *,"       nstep=",knt," ncount=",ncount
      end if
!
!
!     Divide through *ncount*.
!
      fakt=1.0/float(ncount)
!
      tmea(:,:) = tmea(:,:) * fakt
      smea(:,:) = smea(:,:) * fakt
      umea(:,:) = umea(:,:) * fakt
      vmea(:,:) = vmea(:,:) * fakt
      wmea(:,:) = wmea(:,:) * fakt
      cmea(:,:) = cmea(:,:) * fakt
      flwmea(:) = flwmea(:) * fakt
      psimea(:) = psimea(:) * fakt
      simea(:)  = simea(:)  * fakt
      ubmea(:)  = ubmea(:)  * fakt
      vbmea(:)  = vbmea(:)  * fakt
      fhmea(:)  = fhmea(:)  * fakt
      fwmea(:)  = fwmea(:)  * fakt
      txmea(:)  = txmea(:)  * fakt
      tymea(:)  = tymea(:)  * fakt
      tbmea(:)  = tbmea(:)  * fakt
      hfmea(:)  = hfmea(:)  * fakt
      zmea(:)   = zmea(:)   * fakt

      call outmea(tmea,smea,umea,vmea,zmea,simea,fhmea,flwmea,ubmea     &
     &           ,vbmea,wmea,cmea,txmea,tymea,tbmea,hfmea,psimea,fwmea  &
     &           ,ncount,knt)
      ncount=0
      return
      end subroutine aver
