# if defined TS_HADV_UP5 || defined TS_HADV_C6 || defined TS_HADV_WENO5
!----------------------------------------------------------
#  ifdef PREDICTOR
!
!----------------------------------------------------------
! Sixth order advection scheme [PREDICTOR]
!----------------------------------------------------------
!
#  define FLUX5 flux6
#  define FLUX3 flux4
#  define FLUX2 flux2
#  undef  UP5_MASKING
!
            cdif=0.
#  undef FLUX5
#  undef FLUX3
#  undef FLUX2
#  undef UP5_MASKING
#  else
!----------------------------------------------------------
! 5th or 6th order or WENO5 advection schemes [CORRECTOR]
!----------------------------------------------------------
!
#  ifdef TS_HADV_C6
#   define FLUX5 flux6
#   define FLUX3 flux4
#   define FLUX2 flux2
#   undef  UP5_MASKING
#  elif defined TS_HADV_WENO5
#   define FLUX5 flux5_weno
#   define FLUX3 flux3_weno
#   define FLUX2 flux1
#   define UP5_MASKING
#  else
#   define FLUX5 flux5
#   define FLUX3 flux3
#   define FLUX2 flux1
#   define UP5_MASKING
#  endif
!
            cdif=1. 
#  undef FLUX5
#  undef FLUX3
#  undef FLUX2
#  undef UP5_MASKING
#  endif  
!----------------------------------------------------------
# elif defined TS_HADV_C2
!
!----------------------------------------------------------
! Second order advection scheme
!----------------------------------------------------------
!
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=0.5*(t(i,j,k,itrc)+t(i-1,j,k,itrc))
     &                                                 *FlxU(i,j,k)
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=0.5*(t(i,j,k,itrc)+t(i,j-1,k,itrc))
     &                                                 *FlxV(i,j,k)
            enddo
          enddo  

# else   /*  --> UP3 (default) or C4 */
!
!----------------------------------------------------------
! Fourth or Third order advection scheme
!----------------------------------------------------------


!----------------------------------------------------------
#  ifdef PREDICTOR
!----------------------------------------------------------
! PREDICTOR [Fourth order advection scheme]
!
#    define grad WORK
!----------------------------------------------------------
#  else
!----------------------------------------------------------
!  CORRECTOR
!
#  if !(defined TS_HADV_C4 || defined TS_HADV_AKIMA)
#   define TS_HADV_UP3
#  endif
!
#  ifdef TS_HADV_UP3
#   define curv WORK
#  else
#   define grad WORK
#  endif
!------------------------------------------------------------
#  endif
!------------------------------------------------------------


!------------------------------------------------------------
#  ifdef EW_PERIODIC
#   define I_EXT_RANGE istr-1,iend+2
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=istr-1
          else
            imin=max(istr-1,1)
          endif
          if (EAST_INTER) then
            imax=iend+2
          else
            imax=min(iend+2,Lmmpi+1)
          endif
#    define I_EXT_RANGE imin,imax
#   else
#    define I_EXT_RANGE max(istr-1,1),min(iend+2,Lm+1)
#   endif
#  endif
#  ifdef NS_PERIODIC
#   define J_EXT_RANGE jstr-1,jend+2
#  else
#   ifdef MPI
          if (SOUTH_INTER) then
            jmin=jstr-1
          else
            jmin=max(jstr-1,1)
          endif
          if (NORTH_INTER) then
            jmax=jend+2
          else
            jmax=min(jend+2,Mmmpi+1)
          endif
#    define J_EXT_RANGE jmin,jmax
#   else
#    define J_EXT_RANGE max(jstr-1,1),min(jend+2,Mm+1)
#   endif
#  endif
!-------------------------------------------------------------------
          do j=jstr,jend
            do i=I_EXT_RANGE
              FX(i,j)=(t(i,j,k,itrc)-t(i-1,j,k,itrc))
#  ifdef MASKING
     &                                               *umask(i,j)
#  endif
            enddo
          enddo

#  undef I_EXT_RANGE
#  ifndef EW_PERIODIC 
          if (WESTERN_EDGE) then
            do j=jstr,jend
              FX(0,j)=FX(1,j)
            enddo
          endif
          if (EASTERN_EDGE) then
#   ifdef MPI
            do j=jstr,jend
              FX(Lmmpi+2,j)=FX(Lmmpi+1,j)
            enddo
#   else
             do j=jstr,jend
              FX(Lm+2,j)=FX(Lm+1,j)
            enddo
#   endif
          endif
#  endif
!---------------------------------------------------------------------
          do j=jstr,jend 
            do i=istr-1,iend+1
#  if (defined TS_HADV_C4 || defined PREDICTOR)
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
#  elif defined TS_HADV_AKIMA
              cff=2.*FX(i+1,j)*FX(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              else
                grad(i,j)=0.
              endif
#  elif defined TS_HADV_UP3
              curv(i,j)=FX(i+1,j)-FX(i,j)
#  endif
            enddo
          enddo             !--> discard FX
          do j=jstr,jend
            do i=istr,iend+1
#  if (defined TS_HADV_UP3 && !defined PREDICTOR)
              if (FlxU(i,j,k) .gt. 0.) then
                cff=curv(i-1,j)
              else
                cff=curv(i,j)
              endif
              FX(i,j)=0.5*( t(i,j,k,itrc)+t(i-1,j,k,itrc)
     &                           -0.333333333333*cff )*FlxU(i,j,k)
#   else
              FX(i,j)=0.5*( t(i,j,k,itrc)+t(i-1,j,k,itrc)
     &                     -0.333333333333*(grad(i,j)-grad(i-1,j))
     &                                                )*FlxU(i,j,k)
#   endif
            enddo
          enddo            !--> discard grad
!---------------------------------------------------------------------
          do j=J_EXT_RANGE
            do i=istr,iend
              FE(i,j)=(t(i,j,k,itrc)-t(i,j-1,k,itrc))
#  ifdef MASKING
     &                                               *vmask(i,j)
#  endif
            enddo
          enddo
#  undef J_EXT_RANGE
#  ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=istr,iend
              FE(i,0)=FE(i,1)
            enddo
          endif
          if (NORTHERN_EDGE) then
#   ifdef MPI
            do i=istr,iend
              FE(i,Mmmpi+2)=FE(i,Mmmpi+1)
            enddo
#   else
            do i=istr,iend
              FE(i,Mm+2)=FE(i,Mm+1)
            enddo
#   endif
          endif
#  endif
!---------------------------------------------------------------------
          do j=jstr-1,jend+1           !<-- C4 [only for pred]
            do i=istr,iend
#  if (defined TS_HADV_C4 || defined PREDICTOR)
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
#  elif defined TS_HADV_AKIMA
              cff=2.*FE(i,j+1)*FE(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              else
                grad(i,j)=0.
              endif
#  elif defined TS_HADV_UP3
              curv(i,j)=FE(i,j+1)-FE(i,j)
#  endif
            enddo
          enddo            !--> discard FE
          do j=jstr,jend+1
            do i=istr,iend
#  if (defined TS_HADV_UP3 && !defined PREDICTOR)
              if (FlxV(i,j,k) .gt. 0.) then
                cff=curv(i,j-1)
              else
                cff=curv(i,j)
              endif
              FE(i,j)=0.5*( t(i,j,k,itrc)+t(i,j-1,k,itrc)
     &                          -0.333333333333*cff )*FlxV(i,j,k)
#   undef curv
#  else
              FE(i,j)=0.5*( t(i,j,k,itrc)+t(i,j-1,k,itrc)
     &                     -0.333333333333*(grad(i,j)-grad(i,j-1))
     &                                               )*FlxV(i,j,k)
#   undef grad
#  endif
            enddo
          enddo            !--> discard grad
!---------------------------------------------------------------------
# endif /* TS_HADV_UP5 */


