C $Header: /u/gcmpack/MITgcm/model/src/external_forcing.F,v 1.71 2014/08/15 19:19:23 jmc Exp $
C $Name: checkpoint65g $

#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

C--  File external_forcing.F:
C--   Contents
C--   o EXTERNAL_FORCING_U
C--   o EXTERNAL_FORCING_V
C--   o EXTERNAL_FORCING_T
C--   o EXTERNAL_FORCING_S

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_U
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_U(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_U
C     | o Contains problem specific forcing for zonal velocity.
C     *==========================================================*
C     | Adds terms to gU for forcing by external sources
C     | e.g. wind stress, bottom friction etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "EOS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid
	  

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
	  _RL alpha1, U, V, Ux, Uy, timeFactor
	  
C Variables for the isopycnal calculation  
      _RL rhoK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKm1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL slopeX, slopeY, dudz, slopeX1,slopeY1
      _RL drhoY, drhoX, drhoZ
CEOP

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_U(
     U                       gU(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_U(
     U                       gU(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_U(
     U                       gU(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

C     Ocean: Add momentum surface forcing (e.g., wind-stress) in surface level
      IF ( kLev .EQ. kSurface ) THEN
c      DO j=1,sNy
C-jmc: Without CD-scheme, this is OK ; but with CD-scheme, needs to cover [0:sNy+1]
       DO j=0,sNy+1
        DO i=1,sNx+1
          gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingU(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacW(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=0,sNy+1
        DO i=1,sNx+1
         IF ( kSurfW(i,j,bi,bj).EQ.kLev ) THEN
          gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingU(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacW(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF
	  
C    Add terms for advection by flow
C    The strain to be applied is
C	  alpha1=0.4*f0* (1.d0-exp(-(myTime/(2*3.14159/f0))**2.d0))
      alpha1=strainAlpha*f0
C    U, V, du/dx, dV/dx;
C      Lx=200.0e3
C Nx*DXF(1,1,bi,bj)
C      Ly=200.0e3
C Ny*DYF(1,1,bi,bj)
     
C    The contribution to gT is -(U dT/dx+V dT/dy)

      DO j=0,sNy+1
       DO i=1,sNx+1
C    Compute relevant background flow derivatives
      U=-Lx*alpha1/(4.d0*PI)*sin(PI*XC(i,j,bi,bj)/Lx)**2.d0
     &       *(sin(PI*YC(i,j,bi,bj)/Ly)
     &        -3.d0*sin(3.d0*PI*YC(i,j,bi,bj)/Ly))
      V=-Ly*alpha1/(PI)*cos(PI*YC(i,j,bi,bj)/Ly)
     &       *sin(2.d0*PI*XC(i,j,bi,bj)/Lx)
     &       *sin(PI*YC(i,j,bi,bj)/Ly)**2.d0
      Ux=(-0.25d0)*alpha1*sin(2.d0*PI*XG(i,j,bi,bj)/Lx)*(
     &     sin(PI*YC(i,j,bi,bj)/Ly)
     &     -3.d0*sin(3.d0*PI*YC(i,j,bi,bj)/Ly))
      Uy=(-Lx)*alpha1/(4.d0*Ly)*sin(PI*XG(i,j,bi,bj)/Lx)**2.d0
     &    *(cos(PI*YC(i,j,bi,bj)/Ly)
     &      -9.d0*cos(3.d0*PI*YC(i,j,bi,bj)/Ly))
     &       *sin(PI*XG(i,j,bi,bj)/Lx)
C    Add the contribution     
          gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
     &      -( uVel(i,j,kLev,bi,bj) * Ux 
     &         + vVel(i,j,kLev,bi,bj) * Uy
     &         + U*0.5*(uVel(i+1,j,kLev,bi,bj)
     &      -uVel(i-1,j,kLev,bi,bj))*recip_DXC(i,j,bi,bj)
     &         + V*0.5*(uVel(i,j+1,kLev,bi,bj)
     &      -uVel(i,j-1,kLev,bi,bj))*recip_DYG(i,j,bi,bj)
     &      )*_maskW(i,j,kLev,bi,bj)
       ENDDO
      ENDDO
      
C Compute slope for isopycnal background flow
      IF (isopycnalBFlow) THEN
C Only do for non-boundaries, else w=0
      IF ( kLev .GE. 2) THEN
      IF ( kLev .LT. Nr ) THEN
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev,
     I              theta(1-OLx,1-OLy,kLev,bi,bj),
     I              salt(1-OLx,1-OLy,kLev,bi,bj),
     O              rhoK,
     I              kLev, bi, bj, myThid )
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev-1,
     I              theta(1-OLx,1-OLy,kLev-1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev-1,bi,bj),
     O              rhoKm1,
     I              kLev, bi, bj, myThid ) 
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev+1,
     I              theta(1-OLx,1-OLy,kLev+1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev+1,bi,bj),
     O              rhoKp1,
     I              kLev+1, bi, bj, myThid )
      
C Derivatives defined on lower cell edges (U,V pts)
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx+1,sNx+OLx
        sigmaX(i,j)=_maskW(i,j,kLev,bi,bj)
     &        *_recip_dxC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i-1,j))
       ENDDO
      ENDDO
      
      DO j=1-OLy+1,sNy+OLy
       DO i=1-OLx,sNx+OLx
        sigmaY(i,j)=_maskS(i,j,kLev,bi,bj)
     &        *_recip_dyC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i,j-1))
       ENDDO
      ENDDO
C Derivative defined on cell center      
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         sigmaR(i,j)= maskC(i,j,kLev+1,bi,bj)
     &                *maskC(i,j,kLev,bi,bj)
     &                *maskC(i,j,kLev-1,bi,bj)
     &                *0.5d0*(recip_drC(kLev)*rkSign
     &                *(rhoK(i,j)-rhoKm1(i,j))    
     &                +recip_drC(kLev+1)*rkSign
     &                *(rhoKp1(i,j)-rhoK(i,j)))
        ENDDO
       ENDDO
      
	  ENDIF
	  ENDIF
	  
	  ENDIF
      
C BACKGROUND flow from file: Ubar, Vbar
      IF (rC(kLev) .GT. zBFlow) THEN
      DO j=1-OLy+1,sNy+OLy-1
       DO i=1-OLx+2,sNx+OLx-1
C    Compute relevant background flow derivatives
          U=0.5d0*(Ubar(i,j,bi,bj)+Ubar(i-1,j,bi,bj))
     &  *maskC(i,j,kLev,bi,bj)*maskC(i-1,j,kLev,bi,bj)
          V=0.5d0*(Vbar(i,j,bi,bj)+Vbar(i-1,j,bi,bj))
     &  *maskC(i,j,kLev,bi,bj)*maskC(i-1,j,kLev,bi,bj)
          Uy=(Ubar(i-1,j+1,bi,bj)+Ubar(i,j+1,bi,bj)
     &      -Ubar(i-1,j-1,bi,bj)-Ubar(i,j-1,bi,bj))
     &      *0.25d0*recip_DYG(i,j,bi,bj)
     &  *maskC(i-1,j+1,kLev,bi,bj)*maskC(i,j+1,kLev,bi,bj)
     &  *maskC(i-1,j-1,kLev,bi,bj)*maskC(i,j-1,kLev,bi,bj)
          Ux=(Ubar(i,j,bi,bj)-Ubar(i-1,j,bi,bj))
     &      *recip_DXC(i,j,bi,bj)
     &  *maskC(i,j,kLev,bi,bj)*maskC(i-1,j,kLev,bi,bj)
C Add time dependent strain if requested.
      IF ( tauBFlow .GT. 0.d0 ) THEN
        timeFactor=1.d0 - exp(-1.d0*(myTime/tauBFlow)**2.d0)
     	U=U*timeFactor
     	V=V*timeFactor
     	Ux=Ux*timeFactor
     	Uy=Uy*timeFactor
      ENDIF
      
C    Add the contribution     
          gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
     &      -( uVel(i,j,kLev,bi,bj) * Ux 
     &         + vVel(i,j,kLev,bi,bj) * Uy
     &         + U*0.5*(uVel(i+1,j,kLev,bi,bj)
     &      -uVel(i-1,j,kLev,bi,bj))*recip_DXC(i,j,bi,bj)
C     &  *_maskW(i-1,j,kLev,bi,bj)*_maskW(i+1,j,kLev,bi,bj)
     &         + V*0.5*(uVel(i,j+1,kLev,bi,bj)
     &      -uVel(i,j-1,kLev,bi,bj))*recip_DYG(i,j,bi,bj)
C     &  *_maskW(i,j-1,kLev,bi,bj)*_maskW(i,j+1,kLev,bi,bj)
     &      )*_maskW(i,j,kLev,bi,bj)
     
C calculate isopycnal slopes, and add (small slope approx)
C     ... - (Sx*U+Sy*V) d/dz
      IF (isopycnalBFlow) THEN
C slope on U points
        IF (kLev .GE. 2) THEN
        IF (kLev .LT. Nr) THEN
C         slopeX=-1.d0*(sigmaX(i,j))
C     &           /(0.5d0*(sigmaR(i,j)+sigmaR(i-1,j)))
C     &           *_maskW(i,j,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C     &           *maskC(i-1,j,kLev,bi,bj)
C         slopeY=-1.d0*(0.25d0*(sigmaY(i,j)+sigmaY(i,j+1)
C     &           +sigmaY(i-1,j)+sigmaY(i-1,j+1)))
C     &           /(0.5d0*(sigmaR(i,j)+sigmaR(i-1,j)))
C     &           *_maskS(i,j,kLev,bi,bj)
C     &           *_maskS(i,j+1,kLev,bi,bj)
C     &           *_maskS(i-1,j,kLev,bi,bj)
C     &           *_maskS(i-1,j+1,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C     &           *maskC(i-1,j,kLev,bi,bj)
     
C This gives the slope on tracer point (i,j)     
       drhoY=rhoNil*(
     &  sBeta*(salt(i,j+1,kLev,bi,bj)-salt(i,j-1,kLev,bi,bj))
     &  -tAlpha*(theta(i,j+1,kLev,bi,bj)-theta(i,j-1,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i,j+1,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i,j-1,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoX=rhoNil*(
     &  sBeta*(salt(i+1,j,kLev,bi,bj)-salt(i-1,j,kLev,bi,bj))
     &  -tAlpha*(theta(i+1,j,kLev,bi,bj)-theta(i-1,j,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i+1,j,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i-1,j,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoZ=rhoNil*(
     &  sBeta*(salt(i,j,kLev-1,bi,bj)-salt(i,j,kLev+1,bi,bj))
     &  -tAlpha*(theta(i,j,kLev-1,bi,bj)-theta(i,j,kLev+1,bi,bj))
     &  +tCab/2.d0*((theta(i,j,kLev-1,bi,bj)-tRef(kLev-1))**2.d0 
     &  -(theta(i,j,kLev+1,bi,bj)-tRef(kLev+1))**2.d0 ))
       slopeY=drhoY
     &       /(dyC(i,j,bi,bj)+dyC(i,j-1,bi,bj))
     &       /drhoZ
     &       *(drC(kLev)+drC(kLev+1))
     &       *maskC(i,j-1,kLev,bi,bj)
     &       *maskC(i,j+1,kLev,bi,bj)
       slopeX=drhoX
     &       *_recip_dxC(i,j,bi,bj)
     &       /drhoZ
     &       *drC(kLev)
     &       *maskC(i-1,j,kLev,bi,bj)
     &       *maskC(i+1,j,kLev,bi,bj)
C This gives the slope on tracer point (i-1,j) 
       drhoY=rhoNil*(
     &  sBeta*(salt(i-1,j+1,kLev,bi,bj)-salt(i-1,j-1,kLev,bi,bj))
     &  -tAlpha*(theta(i-1,j+1,kLev,bi,bj)-theta(i-1,j-1,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i-1,j+1,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i-1,j-1,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoX=rhoNil*(
     &  sBeta*(salt(i,j,kLev,bi,bj)-salt(i-2,j,kLev,bi,bj))
     &  -tAlpha*(theta(i,j,kLev,bi,bj)-theta(i-2,j,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i,j,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i-2,j,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoZ=rhoNil*(
     &  sBeta*(salt(i-1,j,kLev-1,bi,bj)-salt(i-1,j,kLev+1,bi,bj))
     &  -tAlpha*(theta(i-1,j,kLev-1,bi,bj)-theta(i-1,j,kLev+1,bi,bj))
     &  +tCab/2.d0*((theta(i-1,j,kLev-1,bi,bj)-tRef(kLev-1))**2.d0 
     &  -(theta(i-1,j,kLev+1,bi,bj)-tRef(kLev+1))**2.d0 ))
       slopeY=drhoY
     &       /(dyC(i-1,j,bi,bj)+dyC(i-1,j-1,bi,bj))
     &       /drhoZ
     &       *(drC(kLev)+drC(kLev+1))
     &       *maskC(i-1,j-1,kLev,bi,bj)
     &       *maskC(i-1,j+1,kLev,bi,bj)
       slopeX=drhoX
     &       *_recip_dxC(i-1,j,bi,bj)
     &       /drhoZ
     &       *drC(kLev)
     &       *maskC(i-2,j,kLev,bi,bj)
     &       *maskC(i,j,kLev,bi,bj)   
     
C calculate du/dz    
         dudz = 0.5d0* (recip_drC(kLev)*rkSign
     &    *( uVel(i,j,kLev,bi,bj)-uVel(i,j,kLev-1,bi,bj) )
     &    *maskW(i,j,kLev-1,bi,bj)*maskW(i,j,kLev,bi,bj)
     &    + recip_drC(kLev+1)*rkSign
     &    *( uVel(i,j,kLev+1,bi,bj)-uVel(i,j,kLev,bi,bj) )
     &    *maskW(i,j,kLev+1,bi,bj)*maskW(i,j,kLev,bi,bj)
     &     )
C Check addition is sensible drho/dz<0
       IF (sigmaR(i,j) .LT. 0.d0 ) THEN
       IF (sigmaR(i-1,j) .LT. 0.d0) THEN         
         gU(i,j,kLev,bi,bj) = gU(i,j,kLev,bi,bj)
     &      - (slopeX*U+slopeY*V)*dudz
       ENDIF
       ENDIF
         
        ENDIF
        ENDIF
C End if not a boundary point         
       ENDIF
C ENDIF (isopycnalBFlow)
     
       ENDDO
      ENDDO
      
      ENDIF
	   
#ifdef ALLOW_EDDYPSI
         CALL TAUEDDY_TENDENCY_APPLY_U(
     U                 gU(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
#endif

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gU(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, -1,
     I                 myTime, 0, myThid )

      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_U(
     U                   gU(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_U(
     U                 gU(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_V
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_V(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_V
C     | o Contains problem specific forcing for merid velocity.
C     *==========================================================*
C     | Adds terms to gV for forcing by external sources
C     | e.g. wind stress, bottom friction etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "EOS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
	  _RL alpha1, U, V, Vx, Vy, timeFactor
	  
C Variables for the isopycnal calculation  
      _RL rhoK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKm1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL slopeX, slopeY, dvdz, slopeX1, slopeY1
      _RL drhoY, drhoX, drhoZ
CEOP

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_V(
     U                       gV(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_V(
     U                       gV(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_V(
     U                       gV(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

C     Ocean: Add momentum surface forcing (e.g., wind-stress) in surface level
      IF ( kLev .EQ. kSurface ) THEN
       DO j=1,sNy+1
c       DO i=1,sNx
C-jmc: Without CD-scheme, this is OK ; but with CD-scheme, needs to cover [0:sNx+1]
        DO i=0,sNx+1
          gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingV(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacS(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=1,sNy+1
        DO i=0,sNx+1
         IF ( kSurfS(i,j,bi,bj).EQ.kLev ) THEN
          gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &      +foFacMom*surfaceForcingV(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacS(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF
      
C    Add terms for advection by flow
C    The strain to be applied is
C	  alpha1=0.4*f0* (1.d0-exp(-(myTime/(2*3.14159/f0))**2.d0))
      alpha1=strainAlpha*f0
C    U, V, du/dx, dV/dx;
C      Lx=200.0e3
C Nx*DXF(1,1,bi,bj)
C      Ly=200.0e3
C Ny*DYF(1,1,bi,bj)
     
C    The contribution to gT is -(U dT/dx+V dT/dy)

      DO j=1,sNy+1
       DO i=0,sNx+1
C    Compute relevant background flow derivatives
      U=-Lx*alpha1/(4.d0*PI)*sin(PI*XC(i,j,bi,bj)/Lx)**2.d0
     &       *(sin(PI*YC(i,j,bi,bj)/Ly)
     &        -3.d0*sin(3.d0*PI*YC(i,j,bi,bj)/Ly))
      V=-Ly*alpha1/(PI)*cos(PI*YC(i,j,bi,bj)/Ly)
     &       *sin(2.d0*PI*XC(i,j,bi,bj)/Lx)
     &       *sin(PI*YC(i,j,bi,bj)/Ly)**2.d0
      Vx=(-2.d0)*Ly*alpha1/Lx
     &     *cos(2.d0*PI*XC(i,j,bi,bj)/Lx)
     &     *cos(PI*YG(i,j,bi,bj)/Ly)
     &     *sin(PI*YG(i,j,bi,bj)/Ly)**2.d0
      Vy=(-0.5d0)*alpha1
     &     *(1.d0+3.d0*cos(2.d0*PI*YG(i,j,bi,bj)/Ly))
     &     *sin(2.d0*PI*XC(i,j,bi,bj)/Lx)
     &     *sin(PI*YG(i,j,bi,bj)/Ly)
C    Add the contribution     
          gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &      -( uVel(i,j,kLev,bi,bj) * Vx 
     &         + vVel(i,j,kLev,bi,bj) * Vy
     &         + U*0.5*(vVel(i+1,j,kLev,bi,bj)
     &      -vVel(i-1,j,kLev,bi,bj))*recip_DXC(i,j,bi,bj)
     &         + V*0.5*(vVel(i,j+1,kLev,bi,bj)
     &      -vVel(i,j-1,kLev,bi,bj))*recip_DYG(i,j,bi,bj)
     &      )*_maskW(i,j,kLev,bi,bj)
       ENDDO
      ENDDO
      
C Compute slope for isopycnal background flow
      IF (isopycnalBFlow) THEN
C Only do for non-boundaries, else w=0
      IF ( kLev .GE. 2) THEN
      IF ( kLev .LT. Nr ) THEN
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev,
     I              theta(1-OLx,1-OLy,kLev,bi,bj),
     I              salt(1-OLx,1-OLy,kLev,bi,bj),
     O              rhoK,
     I              kLev, bi, bj, myThid )
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev-1,
     I              theta(1-OLx,1-OLy,kLev-1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev-1,bi,bj),
     O              rhoKm1,
     I              kLev, bi, bj, myThid ) 
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev+1,
     I              theta(1-OLx,1-OLy,kLev+1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev+1,bi,bj),
     O              rhoKp1,
     I              kLev+1, bi, bj, myThid )
      
C Derivatives defined on lower cell edges (U,V pts)
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx+1,sNx+OLx
        sigmaX(i,j)=_maskW(i,j,kLev,bi,bj)
     &        *_recip_dxC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i-1,j))
       ENDDO
      ENDDO
      
      DO j=1-OLy+1,sNy+OLy
       DO i=1-OLx,sNx+OLx
        sigmaY(i,j)=_maskS(i,j,kLev,bi,bj)
     &        *_recip_dyC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i,j-1))
       ENDDO
      ENDDO
C Derivative defined on cell center      
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         sigmaR(i,j)= maskC(i,j,kLev+1,bi,bj)
     &                *maskC(i,j,kLev,bi,bj)
     &                *maskC(i,j,kLev-1,bi,bj)
     &                *0.5d0*(recip_drC(kLev)*rkSign
     &                *(rhoK(i,j)-rhoKm1(i,j))    
     &                +recip_drC(kLev+1)*rkSign
     &                *(rhoKp1(i,j)-rhoK(i,j)))
        ENDDO
       ENDDO
      
	  ENDIF
	  ENDIF
	  
	  ENDIF
      
C BACKGROUND flow from file: Ubar, Vbar
	  IF (rC(kLev) .GT. zBFlow) THEN
      DO j=1-OLy+2,sNy+OLy-1
       DO i=1-OLx+1,sNx+OLx-1
          U=0.5d0*(Ubar(i,j,bi,bj)+Ubar(i,j-1,bi,bj))
          V=0.5d0*(Vbar(i,j,bi,bj)+Vbar(i,j-1,bi,bj))
          Vx=(Vbar(i+1,j-1,bi,bj)+Vbar(i+1,j,bi,bj)
     &      -Vbar(i-1,j-1,bi,bj)-Vbar(i-1,j,bi,bj))
     &      *0.25d0*recip_DXG(i,j,bi,bj)
     &  *maskC(i+1,j-1,kLev,bi,bj)*maskC(i+1,j,kLev,bi,bj)
     &  *maskC(i-1,j-1,kLev,bi,bj)*maskC(i-1,j,kLev,bi,bj)
          Vy=(Vbar(i,j,bi,bj)-Vbar(i,j-1,bi,bj))
     &      *recip_DYC(i,j,bi,bj)
     &    *maskC(i,j,kLev,bi,bj)*maskC(i,j-1,kLev,bi,bj)
C Add time dependent strain if requested.
      IF ( tauBFlow .GT. 0.d0 ) THEN
        timeFactor=1.d0 - exp(-1.d0*(myTime/tauBFlow)**2.d0)
     	U=U*timeFactor
     	V=V*timeFactor
     	Vx=Vx*timeFactor
     	Vy=Vy*timeFactor
      ENDIF
C    Add the contribution     
          gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &      -( uVel(i,j,kLev,bi,bj) * Vx 
     &         + vVel(i,j,kLev,bi,bj) * Vy
     &         + U*0.5*(vVel(i+1,j,kLev,bi,bj)
     &      -vVel(i-1,j,kLev,bi,bj))*recip_DXC(i,j,bi,bj)
C     &    *_maskS(i-1,j,kLev,bi,bj)*_maskS(i+1,j,kLev,bi,bj)
     &         + V*0.5*(vVel(i,j+1,kLev,bi,bj)
     &      -vVel(i,j-1,kLev,bi,bj))*recip_DYG(i,j,bi,bj)
C     &    *_maskS(i,j-1,kLev,bi,bj)*_maskS(i,j+1,kLev,bi,bj)
     &      )*_maskS(i,j,kLev,bi,bj)
          
C calculate isopycnal slopes, and add (small slope approx)
C     ... - (Sx*U+Sy*V) d/dz
      IF (isopycnalBFlow) THEN
C slope on U points
        IF (kLev .GE. 2) THEN
        IF (kLev .LT. Nr) THEN
C         slopeY=-1.d0*(sigmaY(i,j))
C     &           /(0.5d0*(sigmaR(i,j-1)+sigmaR(i,j)))
C     &           *_maskS(i,j,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C     &           *maskC(i,j-1,kLev,bi,bj)
C         slopeX=-1.d0*(0.25d0*(sigmaX(i,j)+sigmaX(i+1,j)
C     &           +sigmaX(i+1,j-1)+sigmaX(i,j-1)))
C     &           /(0.5d0*(sigmaR(i,j)+sigmaR(i,j-1)))
C     &           *_maskW(i,j,kLev,bi,bj)
C     &           *_maskW(i+1,j,kLev,bi,bj)
C     &           *_maskW(i+1,j-1,kLev,bi,bj)
C     &           *_maskW(i,j-1,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C     &           *maskC(i,j-1,kLev,bi,bj)

C This gives the slope on tracer point (i,j)     
       drhoY=rhoNil*(
     &  sBeta*(salt(i,j+1,kLev,bi,bj)-salt(i,j-1,kLev,bi,bj))
     &  -tAlpha*(theta(i,j+1,kLev,bi,bj)-theta(i,j-1,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i,j+1,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i,j-1,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoX=rhoNil*(
     &  sBeta*(salt(i+1,j,kLev,bi,bj)-salt(i-1,j,kLev,bi,bj))
     &  -tAlpha*(theta(i+1,j,kLev,bi,bj)-theta(i-1,j,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i+1,j,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i-1,j,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoZ=rhoNil*(
     &  sBeta*(salt(i,j,kLev-1,bi,bj)-salt(i,j,kLev+1,bi,bj))
     &  -tAlpha*(theta(i,j,kLev-1,bi,bj)-theta(i,j,kLev+1,bi,bj))
     &  +tCab/2.d0*((theta(i,j,kLev-1,bi,bj)-tRef(kLev-1))**2.d0 
     &  -(theta(i,j,kLev+1,bi,bj)-tRef(kLev+1))**2.d0 ))
       slopeY=drhoY
     &       /(dyC(i,j,bi,bj)+dyC(i,j-1,bi,bj))
     &       /drhoZ
     &       *(drC(kLev)+drC(kLev+1))
     &       *maskC(i,j-1,kLev,bi,bj)
     &       *maskC(i,j+1,kLev,bi,bj)
       slopeX=drhoX
     &       *_recip_dxC(i,j,bi,bj)
     &       /drhoZ
     &       *drC(kLev)
     &       *maskC(i-1,j,kLev,bi,bj)
     &       *maskC(i+1,j,kLev,bi,bj)
C This gives the slope on tracer point (i,j-1)     
       drhoY=rhoNil*(
     &  sBeta*(salt(i,j,kLev,bi,bj)-salt(i,j-2,kLev,bi,bj))
     &  -tAlpha*(theta(i,j,kLev,bi,bj)-theta(i,j-2,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i,j,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i,j-2,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoX=rhoNil*(
     &  sBeta*(salt(i+1,j-1,kLev,bi,bj)-salt(i-1,j-1,kLev,bi,bj))
     &  -tAlpha*(theta(i+1,j-1,kLev,bi,bj)-theta(i-1,j-1,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i+1,j-1,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i-1,j-1,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoZ=rhoNil*(
     &  sBeta*(salt(i,j-1,kLev-1,bi,bj)-salt(i,j-1,kLev+1,bi,bj))
     &  -tAlpha*(theta(i,j-1,kLev-1,bi,bj)-theta(i,j-1,kLev+1,bi,bj))
     &  +tCab/2.d0*((theta(i,j-1,kLev-1,bi,bj)-tRef(kLev-1))**2.d0 
     &  -(theta(i,j-1,kLev+1,bi,bj)-tRef(kLev+1))**2.d0 ))
       slopeY1=drhoY
     &       /(dyC(i,j-1,bi,bj)+dyC(i,j-2,bi,bj))
     &       /drhoZ
     &       *(drC(kLev)+drC(kLev+1))
     &       *maskC(i,j-2,kLev,bi,bj)
     &       *maskC(i,j,kLev,bi,bj)
       slopeX1=drhoX
     &       *_recip_dxC(i,j-1,bi,bj)
     &       /drhoZ
     &       *drC(kLev)
     &       *maskC(i-1,j-1,kLev,bi,bj)
     &       *maskC(i+1,j-1,kLev,bi,bj)
C Interp to U points
       slopeX=0.5d0*(slopeX+slopeX1)
       slopeY=0.5d0*(slopeY+slopeY1)
     
C calculate du/dz on U points
         dvdz = 0.5d0* (recip_drC(kLev)*rkSign
     &    *( vVel(i,j,kLev,bi,bj)-vVel(i,j,kLev-1,bi,bj) )
     &    *maskS(i,j,kLev-1,bi,bj)*maskS(i,j,kLev,bi,bj)
     &    + recip_drC(kLev+1)*rkSign
     &    *( vVel(i,j,kLev+1,bi,bj)-vVel(i,j,kLev,bi,bj) )
     &    *maskS(i,j,kLev+1,bi,bj)*maskS(i,j,kLev,bi,bj)
     &     )
C Check addition is sensible drho/dz<0
       IF ( (sigmaR(i,j) .LT. 0.d0) 
     &     .AND. (sigmaR(i-1,j) .LT. 0.d0) ) THEN      
         gV(i,j,kLev,bi,bj) = gV(i,j,kLev,bi,bj)
     &      - (slopeX*U+slopeY*V)*dvdz
       ENDIF
         
        ENDIF
        ENDIF
C End if not a boundary point         
       ENDIF
C ENDIF (isopycnalBFlow)
     
       ENDDO
      ENDDO
      ENDIF

#ifdef ALLOW_EDDYPSI
         CALL TAUEDDY_TENDENCY_APPLY_V(
     U                 gV(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
#endif

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gV(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, -2,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_V(
     U                   gV(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_V(
     U                 gV(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_T
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_T(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_T
C     | o Contains problem specific forcing for temperature.
C     *==========================================================*
C     | Adds terms to gT for forcing by external sources
C     | e.g. heat flux, climatalogical relaxation, etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "SURFACE.h"
#include "EOS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
      INTEGER km, kc, kp
      _RL tmpVar(1:sNx,1:sNy)
      _RL tmpFac, delPI
      _RL recip_Cp
	  _RL alpha1, U, V, timeFactor
	  
C Variables for the isopycnal calculation  
      _RL rhoK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKm1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL slopeX, slopeY, dTdz, drhoY, drhoZ, drhoX
CEOP
#ifdef SHORTWAVE_HEATING
      _RL minusone
      PARAMETER (minusOne=-1.)
      _RL swfracb(2)
      INTEGER kp1
#endif

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingZCoords .AND. useShelfIce ) THEN
       kSurface = -1
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF
      recip_Cp = 1. _d 0 / HeatCapacity_Cp

C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

#ifdef ALLOW_ADDFLUID
      IF ( selectAddFluid.NE.0 .AND. temp_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( temp_addMass - theta(i,j,kLev,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ELSE
         DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( temp_addMass - tRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ADDFLUID */

#ifdef ALLOW_FRICTION_HEATING
      IF ( addFrictionHeating ) THEN
        IF ( fluidIsAir ) THEN
C         conversion from in-situ Temp to Pot.Temp
          tmpFac = (atm_Po/rC(kLev))**atm_kappa
C         conversion from W/m^2/r_unit to K/s
          tmpFac = (tmpFac/atm_Cp) * mass2rUnit
        ELSE
C         conversion from W/m^2/r_unit to K/s
          tmpFac = recip_Cp * mass2rUnit
        ENDIF
        DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &         + frictionHeating(i,j,kLev,bi,bj)
     &          *tmpFac*recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDDO
        ENDDO
      ENDIF
#endif /* ALLOW_FRICTION_HEATING */

      IF ( fluidIsAir .AND. atm_Rq.NE.zeroRL .AND. Nr.NE.1 ) THEN
C--   Compressible fluid: account for difference between moist and dry air
C     specific volume in Enthalpy equation (+ V.dP term), since only the
C     dry air part is accounted for in the (dry) Pot.Temp formulation.
C     Used centered averaging from interface to center (consistent with
C     conversion term in KE eq) and same discretisation ( [T*Q]_bar_k )
C     as for Theta_v in CALC_PHI_HYD

C     conversion from in-situ Temp to Pot.Temp
        tmpFac = (atm_Po/rC(kLev))**atm_kappa
C     conversion from W/kg to K/s
        tmpFac = tmpFac/atm_Cp
        km = kLev-1
        kc = kLev
        kp = kLev+1
        IF ( kLev.EQ.1 ) THEN
          DO j=1,sNy
           DO i=1,sNx
            tmpVar(i,j) = 0.
           ENDDO
          ENDDO
        ELSE
          delPI = atm_Cp*( (rC(km)/atm_Po)**atm_kappa
     &                   - (rC(kc)/atm_Po)**atm_kappa )
          DO j=1,sNy
           DO i=1,sNx
            tmpVar(i,j) = wVel(i,j,kc,bi,bj)*delPI*atm_Rq
     &                  *( theta(i,j,km,bi,bj)*salt(i,j,km,bi,bj)
     &                   + theta(i,j,kc,bi,bj)*salt(i,j,kc,bi,bj)
     &                   )*maskC(i,j,km,bi,bj)*0.25 _d 0
           ENDDO
          ENDDO
        ENDIF
        IF ( kLev.LT.Nr ) THEN
          delPI = atm_Cp*( (rC(kc)/atm_Po)**atm_kappa
     &                   - (rC(kp)/atm_Po)**atm_kappa )
          DO j=1,sNy
           DO i=1,sNx
            tmpVar(i,j) = tmpVar(i,j)
     &                  + wVel(i,j,kp,bi,bj)*delPI*atm_Rq
     &                  *( theta(i,j,kc,bi,bj)*salt(i,j,kc,bi,bj)
     &                   + theta(i,j,kp,bi,bj)*salt(i,j,kp,bi,bj)
     &                   )*maskC(i,j,kp,bi,bj)*0.25 _d 0
           ENDDO
          ENDDO
        ENDIF
        DO j=1,sNy
          DO i=1,sNx
            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &         + tmpVar(i,j)*tmpFac
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDDO
        ENDDO
#ifdef ALLOW_DIAGNOSTICS
        IF ( useDiagnostics ) THEN
C     conversion to W/m^2
          tmpFac = rUnit2mass
          CALL DIAGNOSTICS_SCALE_FILL( tmpVar, tmpFac, 1,
     &                     'MoistCor', kc, 1, 3, bi,bj,myThid )
        ENDIF
#endif /* ALLOW_DIAGNOSTICS */
      ENDIF

C     Ocean: Add temperature surface forcing (e.g., heat-flux) in surface level
      IF ( kLev .EQ. kSurface ) THEN
       DO j=1,sNy
        DO i=1,sNx
          gT(i,j,kLev,bi,bj)=gT(i,j,kLev,bi,bj)
     &      +surfaceForcingT(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=1,sNy
        DO i=1,sNx
         IF ( kSurfC(i,j,bi,bj).EQ.kLev ) THEN
          gT(i,j,kLev,bi,bj)=gT(i,j,kLev,bi,bj)
     &      +surfaceForcingT(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

      IF (linFSConserveTr) THEN
       DO j=1,sNy
        DO i=1,sNx
          IF (kLev .EQ. kSurfC(i,j,bi,bj)) THEN
            gT(i,j,kLev,bi,bj)=gT(i,j,kLev,bi,bj)
     &        +TsurfCor*recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDIF
        ENDDO
       ENDDO
      ENDIF
	  
      
C    Add terms for advection by flow
C    The strain to be applied is
C	  alpha1=0.4*f0* (1.d0-exp(-(myTime/(2*3.14159/f0))**2.d0))
      alpha1=strainAlpha*f0
C    U, V, du/dx, dV/dx;
C      Lx=200.0e3
C Nx*DXF(1,1,bi,bj)
C      Ly=200.0e3
C Ny*DYF(1,1,bi,bj)
     
C    The contribution to gT is -(U dT/dx+V dT/dy)

      DO j=1,sNy
       DO i=1,sNx
C    Compute relevant background flow derivatives
      U=-Lx*alpha1/(4.d0*PI)*sin(PI*XC(i,j,bi,bj)/Lx)**2.d0
     &       *(sin(PI*YC(i,j,bi,bj)/Ly)
     &        -3.d0*sin(3.d0*PI*YC(i,j,bi,bj)/Ly))
      V=-Ly*alpha1/(PI)*cos(PI*YC(i,j,bi,bj)/Ly)
     &       *sin(2.d0*PI*XC(i,j,bi,bj)/Lx)
     &       *sin(PI*YC(i,j,bi,bj)/Ly)**2.d0
C    Add the contribution
C        IF ( (i .LT. 3) .AND. (bi .EQ. 1) ) THEN
C      	  U=0
C        ELSEIF ( (j .LT. 3) .AND. (bj .EQ. 1) ) THEN
C          V=0
C        ELSEIF ( (i .GT. sNx-2) .AND. (bi .EQ. nSx*nPx) ) THEN
C		  U=0
C	    ELSEIF ( (j .GT. sNy-2) .AND. (bj .EQ. nSy*nPy) ) THEN
C	      V=0
C	    ENDIF
	  IF ( (i .GT. 2) .OR. (bi .NE. 1) ) THEN
	   IF ( (j .GT. 2) .OR. (bj .NE. 1) ) THEN
	    IF ( (i .LT. sNx-1) .OR. (bi .NE. nSx*nPx) ) THEN
	    IF ( (j .LT. sNy-1) .OR. (bj .NE. nSy*nPy) ) THEN     
          gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &      -( U*0.5*(theta(i+1,j,kLev,bi,bj)
     &      -theta(i-1,j,kLev,bi,bj))*recip_DXF(i,j,bi,bj)
     &         + V*0.5*(theta(i,j+1,kLev,bi,bj)
     &      -theta(i,j-1,kLev,bi,bj))*recip_DYF(i,j,bi,bj)
     &      )*maskC(i,j,kLev,bi,bj)
         ENDIF
        ENDIF
       ENDIF
      ENDIF
       ENDDO
      ENDDO
      
C Compute slope for isopycnal background flow
      IF (isopycnalBFlow) THEN
C Only do for non-boundaries, else w=0
      IF ( kLev .GE. 2) THEN
      IF ( kLev .LT. Nr ) THEN
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev,
     I              theta(1-OLx,1-OLy,kLev,bi,bj),
     I              salt(1-OLx,1-OLy,kLev,bi,bj),
     O              rhoK,
     I              kLev, bi, bj, myThid )
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev-1,
     I              theta(1-OLx,1-OLy,kLev-1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev-1,bi,bj),
     O              rhoKm1,
     I              kLev, bi, bj, myThid ) 
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev+1,
     I              theta(1-OLx,1-OLy,kLev+1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev+1,bi,bj),
     O              rhoKp1,
     I              kLev+1, bi, bj, myThid )
      
C Derivatives defined on lower cell edges (U,V pts)
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx+1,sNx+OLx
        sigmaX(i,j)=_maskW(i,j,kLev,bi,bj)
     &        *_recip_dxC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i-1,j))
       ENDDO
      ENDDO
      
      DO j=1-OLy+1,sNy+OLy
       DO i=1-OLx,sNx+OLx
        sigmaY(i,j)=_maskS(i,j,kLev,bi,bj)
     &        *_recip_dyC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i,j-1))
       ENDDO
      ENDDO
C Derivative defined on cell center      
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         sigmaR(i,j)= maskC(i,j,kLev+1,bi,bj)
     &                *maskC(i,j,kLev,bi,bj)
     &                *maskC(i,j,kLev-1,bi,bj)
     &                *0.5d0*(recip_drC(kLev)*rkSign
     &                *(rhoK(i,j)-rhoKm1(i,j))    
     &                +recip_drC(kLev+1)*rkSign
     &                *(rhoKp1(i,j)-rhoK(i,j)))
        ENDDO
       ENDDO
      
	  ENDIF
	  ENDIF
	  
	  ENDIF

C     BACKGROUND FLOW ADVECTION FROM FILE: Ubar, Vbar
      IF (rC(kLev) .GT. zBFlow) THEN
C      DO j=1-OLy+1,sNy+OLy-1
C       DO i=1-OLx+1,sNx+OLx-1
      DO j=0,sNy+1
      DO i=0,sNx+1
C  Skip edges to prevent issues (flow should be zero here anyway)
C	  IF ( (i .GT. 2) .OR. (bi .NE. 1) ) THENC
C	   IF ( (j .GT. 2) .OR. (bj .NE. 1) ) THEN
C	    IF ( (i .LT. sNx-1) .OR. (bi .NE. nSx*nPx) ) THEN
C	    IF ( (j .LT. sNy-1) .OR. (bj .NE. nSy*nPy) ) THEN  
	      U=Ubar(i,j,bi,bj)
	      V=Vbar(i,j,bi,bj)  
C Add time dependent strain if requested
     	IF ( tauBFlow .GT. 0.d0 ) THEN
        	timeFactor=1.d0 - exp(-1.d0*(myTime/tauBFlow)**2.d0)
     		U=U*timeFactor
     		V=V*timeFactor
     	ENDIF 
          gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &      -( U*0.5*(theta(i+1,j,kLev,bi,bj)
     &      -theta(i-1,j,kLev,bi,bj))*recip_DXF(i,j,bi,bj)
     &      *maskC(i-1,j,kLev,bi,bj)*maskC(i+1,j,kLev,bi,bj)
     &      + V*0.5*(theta(i,j+1,kLev,bi,bj)
     &      -theta(i,j-1,kLev,bi,bj))*recip_DYF(i,j,bi,bj)
     &      *maskC(i,j-1,kLev,bi,bj)*maskC(i,j+1,kLev,bi,bj)
     &      )*maskC(i,j,kLev,bi,bj)
C         ENDIF
C        ENDIF
C       ENDIF
C      ENDIF

C calculate isopycnal slopes, and add (small slope approx)
C     ... - (Sx*U+Sy*V) d/dz
      IF (isopycnalBFlow) THEN
C slope on centre (T,S) points
        IF (kLev .GE. 2) THEN
        IF (kLev .LT. Nr) THEN
C         slopeY=-1.d0*(0.5d0*(sigmaY(i,j)+sigmaY(i,j+1)))
C     &           /(sigmaR(i,j))
C     &           *_maskS(i,j,kLev,bi,bj)
C     &           *_maskS(i,j+1,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C         slopeX=-1.d0*(0.5d0*(sigmaX(i,j)+sigmaX(i+1,j)))
C     &           /(sigmaR(i,j))
C     &           *_maskW(i,j,kLev,bi,bj)
C     &           *_maskW(i+1,j,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C      	  slopeX=0.0001
C      	  slopeY=0.0001
      	  
       drhoY=rhoNil*(
     &  sBeta*(salt(i,j+1,kLev,bi,bj)-salt(i,j-1,kLev,bi,bj))
     &  -tAlpha*(theta(i,j+1,kLev,bi,bj)-theta(i,j-1,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i,j+1,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i,j-1,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoX=rhoNil*(
     &  sBeta*(salt(i+1,j,kLev,bi,bj)-salt(i-1,j,kLev,bi,bj))
     &  -tAlpha*(theta(i+1,j,kLev,bi,bj)-theta(i-1,j,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i+1,j,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i-1,j,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoZ=rhoNil*(
     &  sBeta*(salt(i,j,kLev-1,bi,bj)-salt(i,j,kLev+1,bi,bj))
     &  -tAlpha*(theta(i,j,kLev-1,bi,bj)-theta(i,j,kLev+1,bi,bj))
     &  +tCab/2.d0*((theta(i,j,kLev-1,bi,bj)-tRef(kLev-1))**2.d0 
     &  -(theta(i,j,kLev+1,bi,bj)-tRef(kLev+1))**2.d0 ))
       slopeY=drhoY
     &       /(dyC(i,j,bi,bj)+dyC(i,j-1,bi,bj))
     &       /drhoZ
     &       *(drC(kLev)+drC(kLev+1))
     &       *maskC(i,j-1,kLev,bi,bj)
     &       *maskC(i,j+1,kLev,bi,bj)
       slopeX=drhoX
     &       *_recip_dxC(i,j,bi,bj)
     &       /drhoZ
     &       *drC(kLev)
     &       *maskC(i-1,j,kLev,bi,bj)
     &       *maskC(i+1,j,kLev,bi,bj)
     
C      rhoLoc(i,j)=rhoNil*(
C     &     sBeta*(sFld(i,j)-refSalt)
C     &   -tAlpha*(tFld(i,j)-refTemp) 
C     &   +tCab/2.d0*(tFld(i,j)-refTemp)**2.d0 )
C     &        + dRho
     
C calculate du/dz    
         dTdz = 0.5d0* (recip_drC(kLev)*rkSign
     &    *( theta(i,j,kLev,bi,bj)-theta(i,j,kLev-1,bi,bj) )
     &    *maskC(i,j,kLev-1,bi,bj)*maskC(i,j,kLev,bi,bj)
     &    + recip_drC(kLev+1)*rkSign
     &    *( theta(i,j,kLev+1,bi,bj)-theta(i,j,kLev,bi,bj) )
     &    *maskC(i,j,kLev+1,bi,bj)*maskC(i,j,kLev,bi,bj)
     &     )
C Check addition is sensible drho/dz<0
       IF ( sigmaR(i,j) .LT. 0.d0 ) THEN            
         gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &      - (slopeX*U+slopeY*V)*dTdz
     &      *maskC(i,j,kLev,bi,bj)
       ENDIF
         
        ENDIF
        ENDIF
C End if not a boundary point         
       ENDIF
C ENDIF (isopycnalBFlow)

       ENDDO
      ENDDO
      ENDIF

#ifdef SHORTWAVE_HEATING
C Penetrating SW radiation
c     IF ( usePenetratingSW ) THEN
       swfracb(1)=abs(rF(kLev))
       swfracb(2)=abs(rF(kLev+1))
       CALL SWFRAC(
     I             2, minusOne,
     U             swfracb,
     I             myTime, 1, myThid )
       kp1 = kLev+1
       IF (kLev.EQ.Nr) THEN
        kp1 = kLev
        swfracb(2)=0. _d 0
       ENDIF
       DO j=1,sNy
        DO i=1,sNx
         gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj)
     &   -Qsw(i,j,bi,bj)*(swfracb(1)*maskC(i,j,kLev,bi,bj)
     &                   -swfracb(2)*maskC(i,j,kp1, bi,bj))
     &    *recip_Cp*mass2rUnit
     &    *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
c     ENDIF
#endif

#ifdef ALLOW_FRAZIL
      IF ( useFRAZIL )
     &     CALL FRAZIL_TENDENCY_APPLY_T(
     U                 gT(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
#endif /* ALLOW_FRAZIL */

#ifdef ALLOW_SHELFICE
      IF ( useShelfIce )
     &     CALL SHELFICE_FORCING_T(
     U                   gT(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
#endif /* ALLOW_SHELFICE */

#ifdef ALLOW_ICEFRONT
      IF ( useICEFRONT )
     &     CALL ICEFRONT_TENDENCY_APPLY_T(
     U                   gT(1-OLx,1-OLy,kLev,bi,bj),
     I                   kLev, bi, bj, myTime, 0, myThid )
#endif /* ALLOW_ICEFRONT */

#ifdef ALLOW_SALT_PLUME
      IF ( useSALT_PLUME )
     &     CALL SALT_PLUME_TENDENCY_APPLY_T(
     U                     gT(1-OLx,1-OLy,kLev,bi,bj),
     I                     iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                     myTime, 0, myThid )
#endif /* ALLOW_SALT_PLUME */

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gT(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, 1,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_T(
     U                   gT(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_BBL
      IF ( useBBL ) CALL BBL_TENDENCY_APPLY_T(
     U                       gT(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_BBL */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_T(
     U                 gT(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: EXTERNAL_FORCING_S
C     !INTERFACE:
      SUBROUTINE EXTERNAL_FORCING_S(
     I           iMin,iMax, jMin,jMax, bi,bj, kLev,
     I           myTime, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R EXTERNAL_FORCING_S
C     | o Contains problem specific forcing for merid velocity.
C     *==========================================================*
C     | Adds terms to gS for forcing by external sources
C     | e.g. fresh-water flux, climatalogical relaxation, etc ...
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"
#include "SURFACE.h"
#include "EOS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     kLev      :: Current vertical level index
C     myTime    :: Current time in simulation
C     myThid    :: Thread Id number
      INTEGER iMin, iMax, jMin, jMax, kLev, bi, bj
      _RL myTime
      INTEGER myThid

#ifdef USE_OLD_EXTERNAL_FORCING
C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j       :: Loop counters
C     kSurface  :: index of surface level
      INTEGER i, j
      INTEGER kSurface
      _RL alpha1, U, V, timeFactor
      
C Variables for the isopycnal calculation  
      _RL rhoK(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKm1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoKp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaX(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL slopeX, slopeY, dSdz, drhoY, drhoZ, drhoX
CEOP

      IF ( fluidIsAir ) THEN
       kSurface = 0
      ELSEIF ( usingZCoords .AND. useShelfIce ) THEN
       kSurface = -1
      ELSEIF ( usingPCoords ) THEN
       kSurface = Nr
      ELSE
       kSurface = 1
      ENDIF

C--   Forcing term
#ifdef ALLOW_AIM
      IF ( useAIM ) CALL AIM_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_AIM */

#ifdef ALLOW_ATM_PHYS
      IF ( useAtm_Phys ) CALL ATM_PHYS_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_ATM_PHYS */

#ifdef ALLOW_FIZHI
      IF ( useFIZHI ) CALL FIZHI_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_FIZHI */

#ifdef ALLOW_ADDFLUID
      IF ( selectAddFluid.NE.0 .AND. salt_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx
            gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( salt_addMass - salt(i,j,kLev,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ELSE
         DO j=1,sNy
          DO i=1,sNx
            gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj)
     &        + addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( salt_addMass - sRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ADDFLUID */

C     Ocean: Add salinity surface forcing (e.g., fresh-water) in surface level
      IF ( kLev .EQ. kSurface ) THEN
       DO j=1,sNy
        DO i=1,sNx
          gS(i,j,kLev,bi,bj)=gS(i,j,kLev,bi,bj)
     &      +surfaceForcingS(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
        ENDDO
       ENDDO
      ELSEIF ( kSurface.EQ.-1 ) THEN
       DO j=1,sNy
        DO i=1,sNx
         IF ( kSurfC(i,j,bi,bj).EQ.kLev ) THEN
          gS(i,j,kLev,bi,bj)=gS(i,j,kLev,bi,bj)
     &      +surfaceForcingS(i,j,bi,bj)
     &      *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
         ENDIF
        ENDDO
       ENDDO
      ENDIF

      IF (linFSConserveTr) THEN
       DO j=1,sNy
        DO i=1,sNx
          IF (kLev .EQ. kSurfC(i,j,bi,bj)) THEN
            gS(i,j,kLev,bi,bj)=gS(i,j,kLev,bi,bj)
     &        +SsurfCor*recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
          ENDIF
        ENDDO
       ENDDO
      ENDIF
      
            
C    Add terms for advection by flow
C    The strain to be applied is
C	  alpha1=0.4*f0* (1.d0-exp(-(myTime/(2*3.14159/f0))**2.d0))
      alpha1=strainAlpha*f0
C    U, V, du/dx, dV/dx;
C      Lx=200.0e3
C Nx*DXF(1,1,bi,bj)
C      Ly=200.0e3
C Ny*DYF(1,1,bi,bj)
     
C    The contribution to gS is -(U dS/dx+V dS/dy)
C    At wall, the normal derivative is zero, so omit.

      DO j=1,sNy
       DO i=1,sNx
C    Compute relevant background flow derivatives
C    Make velocity identically zero near domain edges.
C      XC_rel=(XC(i,j,bi,bj)-XC(2,2,1,1))
C     &       /(XC(sNx-1,sNy-1,nSx*nPx,nSy*nPy)-XC(2,2,1,1))
      U=-Lx*alpha1/(4.d0*PI)*sin(PI*XC(i,j,bi,bj)/Lx)**2.d0
     &       *(sin(PI*YC(i,j,bi,bj)/Ly)
     &        -3.d0*sin(3.d0*PI*YC(i,j,bi,bj)/Ly))
      V=-Ly*alpha1/(PI)*cos(PI*YC(i,j,bi,bj)/Ly)
     &       *sin(2.d0*PI*XC(i,j,bi,bj)/Lx)
     &       *sin(PI*YC(i,j,bi,bj)/Ly)**2.d0
C    Add the contribution 
	  IF ( (i .GT. 2) .OR. (bi .NE. 1) ) THEN
	   IF ( (j .GT. 2) .OR. (bj .NE. 1) ) THEN
	    IF ( (i .LT. sNx-1) .OR. (bi .NE. nSx*nPx) ) THEN
	     IF ( (j .LT. sNy-1) .OR. (bj .NE. nSy*nPy) ) THEN
          gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj)
     &      -( U*0.5*(salt(i+1,j,kLev,bi,bj)
     &      -salt(i-1,j,kLev,bi,bj))*recip_DXF(i,j,bi,bj)
     &         + V*0.5*(salt(i,j+1,kLev,bi,bj)
     &      -salt(i,j-1,kLev,bi,bj))*recip_DYF(i,j,bi,bj)
     &      )*maskC(i,j,kLev,bi,bj)
         ENDIF
        ENDIF
       ENDIF
      ENDIF
       ENDDO
      ENDDO
      
C Compute slope for isopycnal background flow
      IF (isopycnalBFlow) THEN
C Only do for non-boundaries, else w=0
      IF ( kLev .GE. 2) THEN
      IF ( kLev .LT. Nr ) THEN
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev,
     I              theta(1-OLx,1-OLy,kLev,bi,bj),
     I              salt(1-OLx,1-OLy,kLev,bi,bj),
     O              rhoK,
     I              kLev, bi, bj, myThid )
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev-1,
     I              theta(1-OLx,1-OLy,kLev-1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev-1,bi,bj),
     O              rhoKm1,
     I              kLev, bi, bj, myThid ) 
      
          CALL FIND_RHO_2D(
     I              iMin, iMax, jMin, jMax, kLev+1,
     I              theta(1-OLx,1-OLy,kLev+1,bi,bj),
     I              salt(1-OLx,1-OLy,kLev+1,bi,bj),
     O              rhoKp1,
     I              kLev+1, bi, bj, myThid )
      
C Derivatives defined on lower cell edges (U,V pts)
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx+1,sNx+OLx
        sigmaX(i,j)=_maskW(i,j,kLev,bi,bj)
     &        *_recip_dxC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i-1,j))
       ENDDO
      ENDDO
      
      DO j=1-OLy+1,sNy+OLy
       DO i=1-OLx,sNx+OLx
        sigmaY(i,j)=_maskS(i,j,kLev,bi,bj)
     &        *_recip_dyC(i,j,bi,bj)
     &        *(rhoK(i,j)-rhoK(i,j-1))
       ENDDO
      ENDDO
      
C Derivative defined on cell center      
       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         sigmaR(i,j)= maskC(i,j,kLev+1,bi,bj)
     &                *maskC(i,j,kLev,bi,bj)
     &                *maskC(i,j,kLev-1,bi,bj)
     &                *0.5d0*(recip_drC(kLev)*rkSign
     &                *(rhoK(i,j)-rhoKm1(i,j))    
     &                +recip_drC(kLev+1)*rkSign
     &                *(rhoKp1(i,j)-rhoK(i,j)))
        ENDDO
       ENDDO
      
	  ENDIF
	  ENDIF
	  
	  ENDIF
      
C     BACKGROUND FLOW ADVECTION FROM FILE: Ubar, Vbar
      IF (rC(kLev) .GT. zBFlow) THEN
C      DO j=1-OLy+1,sNy+OLy-1
C       DO i=1-OLx+1,sNx+OLx-1
      DO j=0,sNy+1
      DO i=0,sNx+1
C  Skip edges to prevent issues (flow should be zero here anyway)
C	  IF ( (i .GT. 2) .OR. (bi .NE. 1) ) THEN
C	   IF ( (j .GT. 2) .OR. (bj .NE. 1) ) THEN
C	    IF ( (i .LT. sNx-1) .OR. (bi .NE. nSx*nPx) ) THEN
C	    IF ( (j .LT. sNy-1) .OR. (bj .NE. nSy*nPy) ) THEN  
	      U=Ubar(i,j,bi,bj)
	      V=Vbar(i,j,bi,bj)
C Add time factor if requested
     	IF ( tauBFlow .GT. 0.d0 ) THEN
        	timeFactor=1.d0 - exp(-1.d0*(myTime/tauBFlow)**2.d0)
     		U=U*timeFactor
     		V=V*timeFactor
     	ENDIF    
          gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj)
     &      -( U*0.5*(salt(i+1,j,kLev,bi,bj)
     &      -salt(i-1,j,kLev,bi,bj))*recip_DXF(i,j,bi,bj)
     &      *maskC(i-1,j,kLev,bi,bj)*maskC(i+1,j,kLev,bi,bj)
     &      + V*0.5*(salt(i,j+1,kLev,bi,bj)
     &      -salt(i,j-1,kLev,bi,bj))*recip_DYF(i,j,bi,bj)
     &      *maskC(i,j-1,kLev,bi,bj)*maskC(i,j+1,kLev,bi,bj)
     &      )*maskC(i,j,kLev,bi,bj)
C         ENDIF
C        ENDIF
C       ENDIF
C      ENDIF

C calculate isopycnal slopes, and add (small slope approx)
C     ... - (Sx*U+Sy*V) d/dz
      IF (isopycnalBFlow) THEN
C slope on centre (T,S) points
        IF (kLev .GE. 2) THEN
        IF (kLev .LT. Nr) THEN
C         slopeY=-1.d0*(0.5d0*(sigmaY(i,j)+sigmaY(i,j+1)))
C     &           /(sigmaR(i,j))
C     &           *_maskS(i,j,kLev,bi,bj)
C     &           *_maskS(i,j+1,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C         slopeX=-1.d0*(0.5d0*(sigmaX(i,j)+sigmaX(i+1,j)))
C     &           /(sigmaR(i,j))
C     &           *_maskW(i,j,kLev,bi,bj)
C     &           *_maskW(i+1,j,kLev,bi,bj)
C     &           *maskC(i,j,kLev,bi,bj)
C          slopeX=0.0001
C      	  slopeY=0.0001
C      	  
C       drhoY=rhoNil*(
C     &  sBeta*(salt(i,j+1,kLev,bi,bj)-salt(i,j-1,kLev,bi,bj))
C     &  -tAlpha*(theta(i,j+1,kLev,bi,bj)-theta(i,j-1,kLev,bi,bj))
C     &  )
C       drhoX=rhoNil*(
C     &  sBeta*(salt(i+1,j,kLev,bi,bj)-salt(i-1,j,kLev,bi,bj))
C     &  -tAlpha*(theta(i+1,j,kLev,bi,bj)-theta(i-1,j,kLev,bi,bj))
C     &  )
C       drhoZ=rhoNil*(
C     &  sBeta*(salt(i,j,kLev-1,bi,bj)-salt(i,j,kLev+1,bi,bj))
C     &  -tAlpha*(theta(i,j,kLev-1,bi,bj)-theta(i,j,kLev+1,bi,bj))
C     &  ) 
C       slopeY=drhoY
C     &       *_recip_dyC(i,j,bi,bj)
C     &       /drhoZ
C     &       *drC(kLev)
C       slopeX=drhoX
C     &       *_recip_dxC(i,j,bi,bj)
C     &       /drhoZ
C     &       *drC(kLev)
       drhoY=rhoNil*(
     &  sBeta*(salt(i,j+1,kLev,bi,bj)-salt(i,j-1,kLev,bi,bj))
     &  -tAlpha*(theta(i,j+1,kLev,bi,bj)-theta(i,j-1,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i,j+1,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i,j-1,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoX=rhoNil*(
     &  sBeta*(salt(i+1,j,kLev,bi,bj)-salt(i-1,j,kLev,bi,bj))
     &  -tAlpha*(theta(i+1,j,kLev,bi,bj)-theta(i-1,j,kLev,bi,bj))
     &  +tCab/2.d0*((theta(i+1,j,kLev,bi,bj)-tRef(kLev))**2.d0 
     &  -(theta(i-1,j,kLev,bi,bj)-tRef(kLev))**2.d0 ))
       drhoZ=rhoNil*(
     &  sBeta*(salt(i,j,kLev-1,bi,bj)-salt(i,j,kLev+1,bi,bj))
     &  -tAlpha*(theta(i,j,kLev-1,bi,bj)-theta(i,j,kLev+1,bi,bj))
     &  +tCab/2.d0*((theta(i,j,kLev-1,bi,bj)-tRef(kLev-1))**2.d0 
     &  -(theta(i,j,kLev+1,bi,bj)-tRef(kLev+1))**2.d0 ))
       slopeY=drhoY
     &       /(dyC(i,j,bi,bj)+dyC(i,j-1,bi,bj))
     &       /drhoZ
     &       *(drC(kLev)+drC(kLev+1))
     &       *maskC(i,j-1,kLev,bi,bj)
     &       *maskC(i,j+1,kLev,bi,bj)
       slopeX=drhoX
     &       *_recip_dxC(i,j,bi,bj)
     &       /drhoZ
     &       *drC(kLev)
     &       *maskC(i-1,j,kLev,bi,bj)
     &       *maskC(i+1,j,kLev,bi,bj)
     
C      rhoLoc(i,j)=rhoNil*(
C     &     sBeta*(sFld(i,j)-refSalt)
C     &   -tAlpha*(tFld(i,j)-refTemp) 
C     &   +tCab/2.d0*(tFld(i,j)-refTemp)**2.d0 )
C     &        + dRho

C calculate du/dz    
         dSdz = 0.5d0* (recip_drC(kLev)*rkSign
     &    *( salt(i,j,kLev,bi,bj)-salt(i,j,kLev-1,bi,bj) )
     &    *maskC(i,j,kLev-1,bi,bj)*maskC(i,j,kLev,bi,bj)
     &    + recip_drC(kLev+1)*rkSign
     &    *( salt(i,j,kLev+1,bi,bj)-salt(i,j,kLev,bi,bj) )
     &    *maskC(i,j,kLev+1,bi,bj)*maskC(i,j,kLev,bi,bj)
     &     )
C Check addition is sensible drho/dz<0
       IF ( sigmaR(i,j) .LT. 0.d0 ) THEN            
         gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj)
     &      - (slopeX*U+slopeY*V)*dSdz
     &      *maskC(i,j,kLev,bi,bj)
       ENDIF
         
        ENDIF 
        ENDIF 
C End if not a boundary point         
       ENDIF
C ENDIF (isopycnalBFlow)


       ENDDO
      ENDDO
      ENDIF

#ifdef ALLOW_SHELFICE
      IF ( useShelfIce )
     &     CALL SHELFICE_FORCING_S(
     U                   gS(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
#endif /* ALLOW_SHELFICE */

#ifdef ALLOW_ICEFRONT
      IF ( useICEFRONT )
     &     CALL ICEFRONT_TENDENCY_APPLY_S(
     U                   gS(1-OLx,1-OLy,kLev,bi,bj),
     I                   kLev, bi, bj, myTime, 0, myThid )
#endif /* ALLOW_ICEFRONT */

#ifdef ALLOW_SALT_PLUME
      IF ( useSALT_PLUME )
     &     CALL SALT_PLUME_TENDENCY_APPLY_S(
     U                     gS(1-OLx,1-OLy,kLev,bi,bj),
     I                     iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                     myTime, 0, myThid )
#endif /* ALLOW_SALT_PLUME */

#ifdef ALLOW_RBCS
      IF (useRBCS) THEN
        CALL RBCS_ADD_TENDENCY(
     U                 gS(1-OLx,1-OLy,kLev,bi,bj),
     I                 kLev, bi, bj, 2,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_RBCS */

#ifdef ALLOW_OBCS
      IF (useOBCS) THEN
        CALL OBCS_SPONGE_S(
     U                   gS(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_OBCS */

#ifdef ALLOW_BBL
      IF ( useBBL ) CALL BBL_TENDENCY_APPLY_S(
     U                       gS(1-OLx,1-OLy,kLev,bi,bj),
     I                       iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                       myTime, 0, myThid )
#endif /* ALLOW_BBL */

#ifdef ALLOW_MYPACKAGE
      IF ( useMYPACKAGE ) THEN
        CALL MYPACKAGE_TENDENCY_APPLY_S(
     U                 gS(1-OLx,1-OLy,kLev,bi,bj),
     I                 iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                 myTime, 0, myThid )
      ENDIF
#endif /* ALLOW_MYPACKAGE */

#endif /* USE_OLD_EXTERNAL_FORCING */
      RETURN
      END
