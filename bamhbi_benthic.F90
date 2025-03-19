###########################################################################################################################################################################################################################################################################
!
!                                                 B A M H B I    B E N T H I C    M O D E L
!
!
! Original code by A. Capet for publication: https://www.sciencedirect.com/science/article/pii/S146350031630004X
! NEMO adjustment: L. Vandenbulcke
! NECCTON refinements: M. Choblet
##########################################################################################################################################################################################################################################################################

subroutine benthic_init (Kbb, Kmm)
  implicit none
  integer :: Kbb, Kmm
  integer :: istat
  integer :: ji,jj, ikt
  integer :: benthiccounter
 
  !benthic vars
  fCSED=1
  sCSED=2
  sSSED=3
  fSSED=4
  NCrSED=5
  benthiccounter=5

!first total number check
IF (jp_ben /= benthiccounter) THEN
    WRITE(ctmp1,*) 'BAMHBI BENTHIC TRACER CHECK: Wrong number of tracers - Namelist jp_ben = ', jp_ben, 'Model configuration (bamhbi.h90) requires = ', benthiccounter
    CALL ctl_stop(ctmp1)
ENDIF

  DO_2D(0,0,0,0)
        if ( tmask(ji,jj,1) .eq. 0 ) then
           bt(ji,jj,:,fCSED,Kbb)  = 0.
           bt(ji,jj,:,sCSED,Kbb)  = 0.
           bt(ji,jj,:,sSSED,Kbb)  = 0.
           bt(ji,jj,:,sCSED,Kbb)  = 0.
           bt(ji,jj,:,NCrSED,Kbb) = 0.
  END_2D

  tauxb=0. ; tauyb=0.
  !Namelist top option to compute bottom stress unlike in NEMO with log-profile
  if (cur_botstr.eq.0) then
    DO_2D(0,0,0,0)
        ikt=mbkt(ji,jj)
        if (tmask(ji,jj,ikt)==1) then
           c_db(ji,jj)=((1+log(2.083E-5/e3t(ji,jj,ikt,Kmm)))/0.396)**(-2)
        end if
    END_2D
  end if
  
end subroutine benthic_init


subroutine bamhbi_benthic(Kbb,Kmm,Krhs)
    implicit none
    integer, intent(in) :: Kbb, Kmm,Krhs
    real(wp)    :: Cmin(jpi,jpj),Nmin(jpi,jpj),Smin(jpi,jpj)
    real(wp)    :: CFLUXPERHOUR,tfCdeg,pdenit,panox,pnit,seuilpanox,seuilpdenit,BDOX,BNOS,BNHS
    real(wp)    :: fluxbotDOXMAX, fluxbotNOSMAX, fluxbotNHSMAX
    real(wp)    :: resusp_fS,resusp_sS, resusp_fS_max
    real(wp)    :: resusp_fC,resusp_sC, resusp_fC_max
    real(wp)    :: fRESUSP_ero_f,fresusp_ero_s,fresusp_dep ! fresusp= t_bottom / t_crit
    real(wp)    :: tf
    real(wp)    :: dt_ben
    integer     :: ji,jj,jk,ikt
    
    dt_ben=dt
    bt(:,:,:,:,Krhs) = 0._wp

    ! bottom stress from waves and current
    call bottomstress(Kmm)

    do jk=1,jp_benlvl
      DO_2D(0,0,0,0)
        ikt=mbkt(ji,jj)  
        if (tmask(ji,jj,ikt)==1) then    ! we have at least 1 layer of water

          fRESUSP_ero_f   = botstress(ji,jj) / criticalstress_ero_f
          fresusp_ero_s   = botstress(ji,jj) / criticalstress_ero_S
          fresusp_dep     = botstress(ji,jj) / criticalstress_dep

          ! 1.! compute remineralised quantities (sediments that become dissolved), burial, and stock trends
          tf = Q10Factor(ts(ji,jj,ikt,1,Kmm),Q10CDEG)

          bt(ji,jj,jk,fCSED,Krhs)=-fCdegrate*tf*bt(ji,jj,jk,fcsed,Kbb)
          bt(ji,jj,jk,sCSED,Krhs)=-(sCdegrate*tf*bt(ji,jj,jk,scsed,Kbb))-bt(ji,jj,jk,scsed,Kbb)*sCburialrate
          bt(ji,jj,jk,fSSED,Krhs)=-(fSdisrate*tf*bt(ji,jj,jk,fssed,Kbb))
          bt(ji,jj,jk,sSSED,Krhs)=-(sSdisrate*tf*bt(ji,jj,jk,sssed,Kbb))-bt(ji,jj,jk,sssed,Kbb)*sSburialrate

          Smin(ji,jj)= fSdisrate*tf*bt(ji,jj,jk,fssed,Kbb)+ sSdisrate*tf*bt(ji,jj,jk,sssed,Kbb)
          Cmin(ji,jj)= fCdegrate*tf*bt(ji,jj,jk,fcsed,Kbb)+ sCdegrate*tf*bt(ji,jj,jk,scsed,Kbb)
          Nmin(ji,jj)= Cmin(ji,jj) *bt(ji,jj,jk,ncrsed,Kbb)


          ! 2.! resusp and sed flux following fRESUP
          if (fRESUSP_dep.le.1.0) THEN  ! the bottom stress is weak --> we have deposition
            fluxbotcdi(ji,jj) =(1-fRESUSP_dep)*tr(ji,jj,ikt,CDI,Kbb)*WDIA(ji,jj,ikt)  ! fluxbotcdi is negative and downward (as WDIA is downward and negative)
            fluxbotpoc(ji,jj) =(1-fRESUSP_dep)*tr(ji,jj,ikt,POC,Kbb)*WPOM(ji,jj,ikt)
            fluxbottotC(ji,jj)=fluxbotpoc(ji,jj)+fluxbotcdi(ji,jj)

            fluxbotndi(ji,jj) =(1-fRESUSP_dep)*tr(ji,jj,ikt,NDI,Kbb)*WDIA(ji,jj,ikt)
            fluxbotpon(ji,jj) =(1-fRESUSP_dep)*tr(ji,jj,ikt,PON,Kbb)*WPOM(ji,jj,ikt)
            fluxbottotN(ji,jj)=fluxbotpon(ji,jj)+fluxbotndi(ji,jj)

            fluxbotsid(ji,jj) =(1-fRESUSP_dep)*tr(ji,jj,ikt,SID,Kbb)*(-vsinkingrate_Silicious_Detritus)
            fluxbottotS(ji,jj)=fluxbotsid(ji,jj)+SiNrDiatoms*fluxbotndi(ji,jj)

            resusp_fS=0.0 ; resusp_sS=0.0 ; resusp_fC=0.0 ; resusp_sC=0.0
            
          else                        ! too much friction -> erosion
            fluxbotcdi(ji,jj) = 0.0
            fluxbotpoc(ji,jj) = 0.0
            fluxbottotC(ji,jj)= 0.0

            fluxbotndi(ji,jj) = 0.0
            fluxbotpon(ji,jj) = 0.0
            fluxbottotN(ji,jj)= 0.0

            fluxbotsid(ji,jj) = 0.0
            fluxbottotS(ji,jj)= 0.0

            resusp_fC_max= alphaRESUSP*Me_fC*(fRESUSP_ero_f-1.0)
            resusp_fC    = min( resusp_fC_max , 0.99*bt(ji,jj,jk,fcsed,Kbb)/dt_ben )
            resusp_fS_max= alphaRESUSP*Me_fs*(fRESUSP_ero_f-1.0)
            resusp_fS    = min( resusp_fS_max , 0.99*bt(ji,jj,jk,fssed,Kbb)/dt_ben )

            if (resusp_fS.le.1e-10) resusp_fS=0.
            if (resusp_fC.le.1e-10) resusp_fC=0.

            resusp_sC=min( alphaRESUSP*Me_sC*(fRESUSP_ero_s-1.0)*(1.0-resusp_fC/resusp_fC_max), 0.99*bt(ji,jj,jk,scsed,Kbb)/dt_ben)
            resusp_sS=min( alphaRESUSP*Me_sS*(fRESUSP_ero_s-1.0)*(1.0-resusp_fS/resusp_fS_max), 0.99*bt(ji,jj,jk,sssed,Kbb)/dt_ben)

            if (resusp_sS.le.1e-10) resusp_sS=0.
            if (resusp_sC.le.1e-10) resusp_sC=0.
          endif

          ! 3.! update of stock trend (the solid part; mineral part is done in 1.)
          bt(ji,jj,jk,fCSED,Krhs) = bt(ji,jj,jk,fCSED,Krhs) - (   pfCSED  *(fluxbottotC(ji,jj)) + resusp_fC )
          bt(ji,jj,jk,sCSED,Krhs) = bt(ji,jj,jk,sCSED,Krhs) - ( (1-pfCSED)*(fluxbottotC(ji,jj)) + resusp_sC )
          bt(ji,jj,jk,fSSED,Krhs) = bt(ji,jj,jk,fSSED,Krhs) - (   pfSSED  *fluxbottotS(ji,jj) + resusp_fS   )
          bt(ji,jj,jk,sSSED,Krhs) = bt(ji,jj,jk,sSSED,Krhs) - ( (1-pfSSED)*fluxbottotS(ji,jj) + resusp_sS   )

          ! 4.! in case of erosion (resuspension), update bottom flux and diagnostics
          if (fRESUSP_dep.GT.1.0) THEN
            fluxbotPOC(ji,jj) = resusp_fC+resusp_sC                ! here fluxbotPOC > 0 and goes upward
            fluxbotSID(ji,jj) = resusp_fS+resusp_sS
            fluxbotPON(ji,jj) = (resusp_fC+resusp_sC)*bt(ji,jj,jk,ncrsed,Kbb)
            fluxbottotC(ji,jj)= fluxbotPOC(ji,jj)
            fluxbottotS(ji,jj)= fluxbotSID(ji,jj)
            fluxbottotN(ji,jj)= fluxbotPON(ji,jj)
            ! pour la variable aggregat je raisonne comme suit : un flux de resuspension de POC revient a une wPOC vers le haut
            ! pour que le rapport POC sur AGG, soit la masse des aggregats, soit conservé, il faut que WAGG = WPOC (cf cas eps<1 dans calculatesink)
            ! on calcule donc la pseudo vitesse WPOC et on l'impose pour AGG


          endif

          ! 6.! computation of remineralized fluxes
          ! the computation done here for sediment model is based on the work of Marie Suleau
          ! finding function approaching at the best the sediment model of Wijssman
          CFLUXPERHOUR=max(Cmin(ji,jj)*CfluxUnitConv(1),0.05*CfluxUnitConv(2))
          CFLUXPERHOUR=min(CFLUXPERHOUR,2.*CfluxUnitConv(3))
          
          BDOX=max(tr(ji,jj,ikt,DOX,Kbb),10.0)
          BNOS=min(max(tr(ji,jj,ikt,NOS,Kbb),1.0),30.0)
          BNHS=min(max(tr(ji,jj,ikt,NHS,Kbb),0.001),5.0)
          
          call Pnitfun(CFLUXPERHOUR,CFLUXPERHOUR*bt(ji,jj,jk,ncrsed,Kbb),BDOX,BNOS,BNHS,pnit)

          if (tr(ji,jj,ikt,DOX,Kbb).le.2.0)   pnit=0.0

          call Panoxfun(CFLUXPERHOUR,CFLUXPERHOUR*bt(ji,jj,jk,ncrsed,Kbb),BDOX,BNOS,BNHS,panox)

          call Pdenitfun(CFLUXPERHOUR,CFLUXPERHOUR*bt(ji,jj,jk,ncrsed,Kbb),BDOX,BNOS,BNHS,pdenit)

          if (tr(ji,jj,ikt,NOS,Kbb).le.0.01) pdenit=0.0

          !***compute mineralisation and release***
          fluxbotODU(ji,jj)=0.
          fluxbotNHSmax=-tr(ji,jj,ikt,NHS,Kbb) * e3t(ji,jj,ikt,Kmm) / DT_BEN
          fluxbotNOSmax=-tr(ji,jj,ikt,NOS,Kbb) * e3t(ji,jj,ikt,Kmm) / DT_BEN
          fluxbotDOXmax=-tr(ji,jj,ikt,DOX,Kbb) * e3t(ji,jj,ikt,Kmm) / DT_BEN
          fluxbotNHS(ji,jj)=Nmin(ji,jj)*(1.0-pnit)

          if (fluxbotNHS(ji,jj).lt.fluxbotNHSmax) then
            pnit=min(pnit,1.0)
            !autre option ; pnit=min(pnit,1-(fluxbotNHSmax/nmin(ji,jj))) ! mais ça amene d'office zéro ..
            fluxbotNHS(ji,jj)=Nmin(ji,jj)*(1.0-pnit)
          endif

          fluxbotNOS(ji,jj)=Nmin(ji,jj)*pnit-Cmin(ji,jj)*pdenit*0.8
          if (fluxbotNOS(ji,jj).lt.fluxbotNOSmax) Then
            pdenit=min(pdenit,bt(ji,jj,jk,ncrsed,Kbb)*pnit/0.8)
            fluxbotNOS(ji,jj)=Nmin(ji,jj)*pnit-Cmin(ji,jj)*pdenit*0.8
          endif

          fluxbotDOX(ji,jj)=-Cmin(ji,jj)*(1.0-pdenit-panox*psoliddepo)*OCrdegrad- Nmin(ji,jj)*pnit*ONrnitrif

          if (fluxbotDOX(ji,jj).lt.fluxbotDOXmax) Then
            panox=1
            fluxbotDOX(ji,jj)=-Cmin(ji,jj)*(1.0-pdenit-panox*psoliddepo)*OCrdegrad- Nmin(ji,jj)*pnit*ONrnitrif
            fluxbotODU(ji,jj)=fluxbotDOXmax-fluxbotDOX(ji,jj)
            fluxbotDOX(ji,jj)=fluxbotDOXmax
          endif

          Fluxbotpho(ji,jj)=Nmin(ji,jj)*PNRedfield
          fluxbotSIO(ji,jj)=Smin(ji,jj)
          fluxbotDIC(ji,jj)=Cmin(ji,jj)
        
          burialC(ji,jj)    = bt(ji,jj,jk,scsed,Kbb)*sCburialrate
          burialS(ji,jj)    = bt(ji,jj,jk,sssed,Kbb)*sSburialrate
          denitinsed(ji,jj) = Cmin(ji,jj)*pdenit*0.8
          pdenit2D(ji,jj) = pdenit
          pnit2D(ji,jj)   = pnit
          panox2D(ji,jj)  = panox
          
          ! 5.! update of NCrSED
          if ((bt(ji,jj,jk,fcsed,Kbb)+bt(ji,jj,jk,scsed,Kbb)).lt.0.05) then
            bt(ji,jj,jk,NCrSED,Krhs)=0.
          else
            bt(ji,jj,jk,NCrSED,Krhs)= (1/(bt(ji,jj,jk,fcsed,Kbb)+bt(ji,jj,jk,scsed,Kbb)) ) * ( (-fluxbottotN(ji,jj)-Nmin(ji,jj)) -bt(ji,jj,jk,ncrsed,Kbb)*(-fluxbottotC(ji,jj)-Cmin(ji,jj)) )
          end if
        end if ! tmask==1
    END_2D
  end do     ! jk benthic levels
  
end subroutine bamhbi_benthic

  SUBROUTINE trc_wbstress(kt,kmm)
    !!----------------------------------------------------------------------
    !!                     ***  trc_ini_wbstress  ***
    !!
    !! ** Purpose :  Compute wave bottom stress effect to be used in bottom stress
    !!
    !!----------------------------------------------------------------------
    use oce,     ONLY : rhop
    use phycst,  ONLY : rpi
    !
    INTEGER, INTENT( in ) ::   kt,kmm              ! ocean time-step index
    !
    INTEGER :: ji, jj
    !
    REAL(wp) :: ks, T, H, d , x, b, y, Uw, A, Rw , n, bottomRho
    REAL(wp) :: fw, fw_rough, fw_smooth

    SELECT CASE (waves_botstr)
    CASE (0)
       !wavestress = 0._wp

    CASE (1)
       CALL fld_read(kt=kt, kn_fsbc=1, sd=sf_wbstress)
       wavestress(:,:)=sf_wbstress(1)%fnow(:,:,1) ! stress from file
       
    CASE (2,3)
       CALL fld_read(kt=kt, kn_fsbc=1, sd=sf_wbstress)
       ! COMPUTE bottom STRESS due to WAVES
       ! using Hunt's method (1979), see https://eprints.hrwallingford.com/588/1/TR155.pdf 
       ! input: wave period and significant height, data from namben_wbstress
       ks = 2.5 * 250E-06
       DO_2D(1,1,1,1)
             IF ( tmask(ji,jj,mbkt(ji,jj)) .eq. 1 ) THEN
             d=gdept(ji,jj,mbkt(ji,jj),Kmm)
             IF ( d .lt. 200.0_wp ) THEN

                T = sf_wbstress(1)%fnow(ji,jj,1) ! peak period
                H = sf_wbstress(2)%fnow(ji,jj,1) ! significant height

                if (T.eq.0.0 .or. H.eq.0) then
                   wavestress(ji,jj)=0.0

                else
                  x=((2*rpi/T)**2)*d/9.81
                  b=1/(1+x*(0.66667+x*(0.35550+x*(0.16084+x*(0.06320+x*(0.02174+x*(0.00654+x*(0.00171+x*(0.00039+x*0.00011)))))))))
                  y=sqrt(x**2+b*x)
                  Uw=H/2*sqrt(2*9.81*y/(d*sinh(2*y)))

                  A=Uw*T/(2*rpi)
                  fw_rough=0.237*(A/ks)**(-0.52)
                  Rw=Uw*A/0.0000012  ! kinematic viscosity nu=0.0000012 [m2/s]
                  if (Rw<=5e5) then  ! laminar
                     B=2      ; N=0.5
                  else               ! turbulent
                     B=0.0521 ; N=0.187
                  end if
                  fw_smooth=B*Rw**(-N)
                  fw=max(fw_rough,fw_smooth)
                  bottomRho=rhop(ji,jj,mbkt(ji,jj))

                  wavestress(ji,jj)=0.5*bottomRho*fw*(Uw**2)
                end if
               
             END IF
             END IF
       END_2D
    END SELECT

    IF (waves_botstr .eq. 3) THEN 
       wave_direction = sf_wbstress(3)%fnow(ji,jj,1)  
    END IF 

  END SUBROUTINE trc_wbstress


subroutine BOTTOMSTRESS(Kmm)
!-------------------------------------------------
! total bottom stress from current and total stress. Bottom stress from currents can be computed as in NEMO (namelist top: cur_botstr = 1) or with a log-profile formula. 
!--------------------------------------------------
  IMPLICIT NONE 
  integer, intent(in) :: Kmm
  integer  :: ji,jj
  real(wp) :: u2,v2,sqrtu2C
  real(wp) :: epsilon=1.E-30
  real(wp) :: currentstress
  real(wp) :: tm
  real(wp) :: vec_sum
  real(wp) :: phi !the vector between currents and waves
  real(wp) :: phi_c !angle currents
  real(wp) :: phi_w !angle waves

  ! M.C 05/2024 :additions for computing currentstress as in nemo (namben_curstress)
  ! code follows src/OCE/DIA/diawri.F90
  real(wp):: zztmp, zztmp2
  REAL(wp), DIMENSION(A2D(     0))     ::   z2d  
  
  ! BOTTOM STRESS FROM WAVES is computed in trc_wbstress.F90 
  if (cur_botstr.eq.1) then 
        zztmp = rho0 * 0.25_wp
        z2d(:,:) = 0._wp
  end if

  ! BOTTOM STRESS FROM CURRENT, AND TOTAL STRESS
  DO_2D(0,0,0,0)
      if (tmask(ji,jj,mbkt(ji,jj)).eq.1) then    ! we have at least 1 layer of water
        u2=0.0 ; v2=0.0
         if (cur_botstr.eq.0) then 
             if (umask(ji,jj,mbku(ji,jj)).eq.1)        u2=uu(ji,jj,mbku(ji,jj),Kmm)
             if (ji.gt.1) then                                                                       ! this will never be the case in R/deSolve
                 if (umask(ji-1,jj,mbku(ji-1,jj)).eq.1) u2=(u2+uu(ji-1,jj,mbku(ji-1,jj),Kmm))/2
             end if
             if (vmask(ji,jj,mbkv(ji,jj)).eq.1)        v2=vv(ji,jj,mbkv(ji,jj),Kmm)
             if (jj.gt.1) then 
                 if (vmask(ji,jj-1,mbkv(ji,jj-1)).eq.1) v2=(v2+vv(ji,jj-1,mbkv(ji,jj-1),Kmm))/2
             end if
             sqrtu2=sqrt(u2**2 + v2**2 + epsilon)
             tauxb(ji,jj)=c_db(ji,jj)*sqrtu2*u2
             tauyb(ji,jj)=c_db(ji,jj)*sqrtu2*v2
             currentstress=rhop(ji,jj,mbkt(ji,jj)) * sqrt((tauxb(ji,jj))**2 + (tauyb(ji,jj))**2)
         elseif (cur_botstr.eq.1) then 
             !code follows src/OCE/DIA/diawri.F90
              zztmp2 = (  ( rCdU_bot(ji+1,jj)+rCdU_bot(ji  ,jj) ) * uu(ji  ,jj,mbku(ji  ,jj),Kmm)  )**2   &
              &   + (  ( rCdU_bot(ji  ,jj)+rCdU_bot(ji-1,jj) ) * uu(ji-1,jj,mbku(ji-1,jj),Kmm)  )**2   &
              &   + (  ( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj  ) ) * vv(ji,jj  ,mbkv(ji,jj  ),Kmm)  )**2   &
              &   + (  ( rCdU_bot(ji,jj  )+rCdU_bot(ji,jj-1) ) * vv(ji,jj-1,mbkv(ji,jj-1),Kmm)  )**2 
              z2d(ji,jj) = zztmp * SQRT( zztmp2 ) * tmask(ji,jj,1)
              currentstress=z2d(ji,jj)
         end if

        ! TOTAL BOTTOM STRESS
        if (waves_botstr .eq. 3) then 
            ! sum of vectors
            phi_c = ATAN2(v2,u2)
            phi_w = wave_direction(ji,jj)        ! angle provided
            ! convert phi_w to the same coordinate system as phi_c 
            ! 1. Substract 90 degrees 2. 0-360 to (-180,180) conversion 3. invert angle (from clockwise to counterclockwise)
            ! 4. Convert to radians
            phi_w = - (MODULO((phi_w-90+180),360.)+180)/180*rpi
            ! nonlinear enhancement,  following R. Soulsby, dynamics of marine sands or also Gayer 2006
            ! (https://link.springer.com/article/10.1007/s10236-006-0070-5)
            tm = currentstress*(1 + 1.2 * (wavestress(ji,jj)/(currentstress+wavestress(ji,jj)))**(3.2))
            ! difference between to angles
            phi=ABS((MODULO(phi_c-phi_w+rpi,2*rpi)-rpi))
            vec_sum = SQRT((tm + wavestress(ji,jj)*COS(phi))**2 + (wavestress(ji,jj)*SIN(phi))**2)
            botstress(ji,jj)=min(1.0, vec_sum) !saveguard against too big wavestresses
        else
           botstress(ji,jj)=min(1.0, currentstress+wavestress(ji,jj)) !saveguard against too big wavestresses
        end if

        curstress(ji,jj)=currentstress
        !resuspension counter for fast and slow decaying stock (1==resuspension, 0==deposition) 
        if (botstress(ji,jj).gt. criticalstress_ero_s) then 
            resuspcounter_s(ji,jj) = 1
        else
            resuspcounter_s(ji,jj) = 0
        endif

        if (botstress(ji,jj).gt. criticalstress_ero_f) then 
            resuspcounter_f(ji,jj) = 1
        else
            resuspcounter_f(ji,jj) = 0
        endif
      end if
  END_2D
  
end subroutine bottomstress




subroutine Panoxfun(CFLUXPERHOUR,nmin,BDOX,bN3n,bN4n,Panox)
IMPLICIT NONE
real(wp),intent(IN) :: CFLUXPERHOUR,nmin,BDOX,bN3n,bN4n
real(wp)            :: Panox


panox = 1.056&
     -0.132*log(CFLUXPERHOUR)&
     +0.005*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)&
     +0.057*log(CFLUXPERHOUR)*log(BDOX)&
     -0.017*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BDOX)&
     -0.008*log(BDOX)*log(BDOX)

Panox=max(Panox,0.d0)
Panox=min(Panox,1.d0)
end subroutine Panoxfun

subroutine Pnitfun(CFLUXPERHOUR,nmin,BDOX,BNOS,BNHS,Pnit)
IMPLICIT NONE
real(wp),intent(IN)  :: CFLUXPERHOUR,nmin,BDOX,BNOS,BNHS
real(wp),intent(OUT) :: Pnit

pnit = -6.280&
     -0.286*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)&
     +0.127*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)&
     +0.006*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)&
     +0.098*log(CFLUXPERHOUR)*log(BDOX)&
     +0.202*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BDOX)&
     -0.021*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BDOX)&
     -0.031*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BDOX)*log(BDOX)&
     -0.066*log(CFLUXPERHOUR)*log(BNHS)&
     +0.009*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BNOS)&
     +0.823*log(BDOX)&
     -0.028*log(BDOX)*log(BDOX)&
     -0.012*log(BDOX)*log(BNHS)&
     +0.119*log(BNHS)&
     +0.015*log(BNHS)*log(BNHS)
if (Nmin.eq.0) then
 pnit=0
else
 pnit=exp(pnit)/(Nmin)
end if

end subroutine Pnitfun


subroutine Pdenitfun(CFLUXPERHOUR,nmin,BDOX,BNOS,bN4n,Pdenit)
IMPLICIT NONE
real(wp),intent(IN) :: CFLUXPERHOUR,nmin,BDOX,BNOS,bN4n
real(wp)            :: Pdenit

pdenit = -5.475&
     -0.786*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)&
     +0.662*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)&
     +0.042*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)&
     +0.064*log(CFLUXPERHOUR)*log(BDOX)&
     +0.794*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BDOX)&
     -0.082*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BDOX)&
     -0.122*log(CFLUXPERHOUR)*log(CFLUXPERHOUR)*log(BDOX)*log(BDOX)&
     +0.077*log(BDOX)*log(BDOX)&
     -0.155*log(BDOX)*log(BNOS)&
     +0.875*log(BNOS)&
     +0.046*log(BNOS)*log(BNOS)
pdenit=exp(pdenit)/(CFLUXPERHOUR)

end subroutine Pdenitfun
