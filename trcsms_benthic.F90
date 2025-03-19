MODULE trcsms_ben
#include "bamhbi.h90"
   !!======================================================================
   !!                         ***  MODULE trcsms_ben  ***
   !! TOP :   Main module of the Benthic tracers
   !!======================================================================
   !! History :      !  2017  (T. Lovato) Create Benthic infrastructure
   !!----------------------------------------------------------------------
   !! trc_sms_ben       : Benthic system  main routine
   !! trc_sms_ben_alloc : allocate arrays specific to Benthic sms
   !!----------------------------------------------------------------------
   USE oce_trc         ! Ocean variables
   USE trc             ! TOP variables
   USE par_ben         ! TOP Benthic system
   USE bamhbi
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_sms_ben       ! called by trcsms.F90 module

   ! Defined HERE the arrays specific to Benthic system and ALLOCATE them in trc_sms_ben_alloc

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2017)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

#include "do_loop_substitute.h90"
#include "domzgr_substitute.h90"

CONTAINS

   SUBROUTINE trc_sms_ben( kt , kbb, kmm, krhs )
      !!----------------------------------------------------------------------
      !!                     ***  trc_sms_my_trc  ***
      !!
      !! ** Purpose :   main routine of MY_TRC model
      !!
      !! ** Method  : -
      !!----------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kbb, Kmm, Krhs
      INTEGER ::  ji,jj,ikt,l
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('trc_sms_ben')
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) ' trc_sms_ben  : Benthic system'
      IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~~~~'

      ! call to BGC benthic model, the fluxes are saved in bt(:,:,:,1:6,Krhs)
      call bamhbi_benthic ( Kbb, Kmm , Krhs )

      ! Euler time integration in the benthic model, save result in Kmm
      !! rdttrc is valid for euler, 2*rdttrc for leapfrog ; it's a 1D array (length=lvl, but all the same in our case)
      !! L.V. in nemo 4.2, replaced with rn_Dt (Euler) ; rdt_trc would be for Leapfrog and runge-kutta
      bt(:,:,:,:,Kmm) = bt(:,:,:,:,Kbb) + rn_Dt * bt(:,:,:,:,Krhs)
      where(bt(:,:,:,NCrSED,Kmm)<0.0_wp) bt(:,:,:,NCrSED,Kmm)=0.0_wp

      ! Fluxes to the pelagic variables
      DO_2D(1,1,1,1)
            ikt = mbkt(ji,jj)
            if (tmask(ji,jj,ikt).eq.1) then   ! we have at least 1 layer of water
               tr(ji,jj,ikt,CDI,Krhs) = tr(ji,jj,ikt,CDI,Krhs) + fluxbotCDI(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,NDI,Krhs) = tr(ji,jj,ikt,NDI,Krhs) + fluxbotNDI(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,NOS,Krhs) = tr(ji,jj,ikt,NOS,Krhs) + fluxbotNOS(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,NHS,Krhs) = tr(ji,jj,ikt,NHS,Krhs) + fluxbotNHS(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,SIO,Krhs) = tr(ji,jj,ikt,SIO,Krhs) + fluxbotSIO(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,DOX,Krhs) = tr(ji,jj,ikt,DOX,Krhs) + fluxbotDOX(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,DIC,Krhs) = tr(ji,jj,ikt,DIC,Krhs) + fluxbotDIC(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,ODU,Krhs) = tr(ji,jj,ikt,ODU,Krhs) + fluxbotODU(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,POC,Krhs) = tr(ji,jj,ikt,POC,Krhs) + fluxbotPOC(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,PON,Krhs) = tr(ji,jj,ikt,PON,Krhs) + fluxbotPON(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,SID,Krhs) = tr(ji,jj,ikt,SID,Krhs) + fluxbotSID(ji,jj)/e3t(ji,jj,ikt,Kmm)
               tr(ji,jj,ikt,PHO,Krhs) = tr(ji,jj,ikt,PHO,Krhs) + fluxbotPHO(ji,jj)/e3t(ji,jj,ikt,Kmm)
            end if
      END_2D

      call iom_put('botfluxCDI',fluxbotCDI)
      call iom_put('botfluxNDI',fluxbotNDI)
      call iom_put('botfluxNOS',fluxbotNOS)
      call iom_put('botfluxNHS',fluxbotNHS)
      call iom_put('botfluxSIO',fluxbotSIO)
      call iom_put('botfluxDOX',fluxbotDOX)
      call iom_put('botfluxDIC',fluxbotDIC)
      call iom_put('botfluxODU',fluxbotODU)
      call iom_put('botfluxPOC',fluxbotPOC)
      call iom_put('botfluxPON',fluxbotPON)
      call iom_put('botfluxSID',fluxbotSID)
      call iom_put('botfluxPHO',fluxbotPHO)


      call iom_put('burialC',burialC)
      call iom_put('burialS',burialS)
      call iom_put('denitinsed',denitinsed) 
      call iom_put('pdenit2D',pdenit2D)
      call iom_put('pnit2D',pnit2D)
      call iom_put('panox2D',panox2D)
      call iom_put('fluxbotpoc_resusp',fluxbotpoc_resusp)
      call iom_put('fluxbotpon_resusp',fluxbotpon_resusp)
      call iom_put('fluxbotsid_resusp',fluxbotsid_resusp)
      call iom_put('fluxbotpoc_depo',fluxbotpoc_depo)
      call iom_put('fluxbotpon_depo',fluxbotpon_depo)
      call iom_put('fluxbotsid_depo',fluxbotsid_depo)
      
      !
   END SUBROUTINE trc_sms_ben

   !!======================================================================
END MODULE trcsms_ben



