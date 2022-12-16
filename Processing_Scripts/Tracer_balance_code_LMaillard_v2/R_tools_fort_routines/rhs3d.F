#include "cppdefs.h"
#ifdef SOLVE3D
# define UPSTREAM

# define NO_SLIP_BOTTOM


 
      subroutine rhs3d (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call rhs3d_tile (istr,iend,jstr,jend,  A3d(1,1), A3d(1,2),
     &                             A2d(1,1), A2d(1,2), A2d(1,3),
     &                             A2d(1,1), A2d(1,2), A2d(1,3),
     &                             A2d(1,4), A2d(1,5), A2d(1,6))
      return
      end
 
      subroutine rhs3d_tile (istr,iend,jstr,jend, ru,rv, CF,FC,DC,
     &                                wrk1,wrk2, UFx,UFe, VFx,VFe)
      implicit none
      integer istr,iend,jstr,jend, imin,imax,jmin,jmax, i,j,k
# include "param.h"
      real ru(PRIVATE_2D_SCRATCH_ARRAY,N),
     &     rv(PRIVATE_2D_SCRATCH_ARRAY,N),
     &     CF(PRIVATE_1D_SCRATCH_ARRAY,0:N),  cff,
     &     FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),  gamma,
     &     DC(PRIVATE_1D_SCRATCH_ARRAY,0:N), 
     &     wrk1(PRIVATE_2D_SCRATCH_ARRAY),
     &     wrk2(PRIVATE_2D_SCRATCH_ARRAY),
     &      UFx(PRIVATE_2D_SCRATCH_ARRAY),
     &      UFe(PRIVATE_2D_SCRATCH_ARRAY),
     &      VFx(PRIVATE_2D_SCRATCH_ARRAY),
     &      VFe(PRIVATE_2D_SCRATCH_ARRAY)
      parameter (gamma=0.25)
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "forces.h"
# include "scalars.h"
# include "mixing.h"
# ifdef DIAGNOSTICS_UV
#include "diags_uv.h"
# endif
# ifdef DIAGNOSTICS_VRT
#include "diags_vrt.h"
# endif

# ifdef DIAGNOSTICS_PV
#include "diags_pv.h"
# endif

!
# include "compute_auxiliary_bounds.h"
!
! Some parts of this code --- bottom stress term and vertical
! integration of r.h.s. to compute forcing terms for the barotropic
! mode --- needs to be computed only once per main time step, during
! predictor or corrector stage, whichever immediately preceeds the
! barotropic mode stepping.  This needs to be done differently for
! the precictor- and corrector-coupled versions.  
!
# ifdef CORR_COUPLED_MODE
#  define COUPLED_STAGE_ONLY nrhs.eq.3
# else
#  define COUPLED_STAGE_ONLY nrhs.eq.nstp
# endif

      do k=1,N
!
# if defined UV_COR || (defined CURVGRID && defined UV_ADV)


!
! Add in Coriolis and curvilinear transformation terms, if any.
!
        do j=jstrV-1,jend
          do i=istrU-1,iend
            cff=0.5*Hz(i,j,k)*(
#  ifdef UV_COR
     &              fomn(i,j)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     &             +0.5*( (v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)
     &                   -(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
#  endif
     &                                                             )
            UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
            VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
          enddo
        enddo
        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)+0.5*(UFx(i,j)+UFx(i-1,j))
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MCor(i,j,k,1) = 0.5*(UFx(i,j)+UFx(i-1,j))
            endif
# elif defined DIAGNOSTICS_VRT
            if (nnew.ne.3) then
              if (k.eq.1) then
                wrkCor(i,j,1) = 0.5*(UFx(i,j)+UFx(i-1,j))
              else
                wrkCor(i,j,1) = wrkCor(i,j,1) + 0.5*(UFx(i,j)+UFx(i-1,j))
              endif
            endif
#  endif
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-0.5*(VFe(i,j)+VFe(i,j-1))
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MCor(i,j,k,2) = -0.5*(VFe(i,j)+VFe(i,j-1))
            endif
# elif defined DIAGNOSTICS_VRT
            if (nnew.ne.3) then
              if (k.eq.1) then
                wrkCor(i,j,2) = -0.5*(VFe(i,j)+VFe(i,j-1))
              else
                wrkCor(i,j,2) = wrkCor(i,j,2) -0.5*(VFe(i,j)+VFe(i,j-1))
              endif
            endif
#  endif
          enddo
        enddo
# endif
# ifdef UV_ADV
!
! Add in horizontal advection of momentum: Compute diagonal [UFx,VFe]
! and off-diagonal [UFe,VFx] components of tensor of momentum flux
! due to horizontal advection; after that add divergence of these
! terms to r.h.s.
!
#  define uxx wrk1
#  define Huxx wrk2
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then        ! Sort out bounding indices of
          imin=istrU                  ! extended ranges: note that in
        else                          ! the vicinity of physical
          imin=istrU-1                ! boundaries values at the
        endif                         ! extremal points of stencil
        if (EASTERN_EDGE) then        ! are not available, so an
          imax=iend                   ! extrapolation rule needs to
        else                          ! be applied. Also note that
          imax=iend+1                 ! for this purpose periodic
        endif                         ! ghost points and MPI margins 
#  else
        imin=istrU-1                  ! are not considered as
        imax=iend+1                   ! physical boundaries.
#  endif
        do j=jstr,jend
          do i=imin,imax
            uxx(i,j)=u(i-1,j,k,nrhs)-2.*u(i,j,k,nrhs)
     &                                +u(i+1,j,k,nrhs)
            Huxx(i,j)=FlxU(i-1,j,k)-2.*FlxU(i,j,k)
     &                                  +FlxU(i+1,j,k)
          enddo
        enddo
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr,jend
            uxx(istrU-1,j) =uxx(istrU,j)
            Huxx(istrU-1,j)=Huxx(istrU,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
            uxx(iend+1,j) =uxx(iend,j)
            Huxx(iend+1,j)=Huxx(iend,j)
          enddo
        endif
#  endif
        do j=jstr,jend
          do i=istrU-1,iend
#  ifdef UPSTREAM
            cff=FlxU(i,j,k)+FlxU(i+1,j,k)-0.125*( Huxx(i  ,j)
     &                                           +Huxx(i+1,j))
            UFx(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*uxx(i  ,j)
     &                                  +min(cff,0.)*uxx(i+1,j)
     &                                                      ))
#  else
            UFx(i,j)=0.25*( u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
     &                         -0.125*(uxx(i,j)+uxx(i+1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i+1,j,k)
     &                      -0.125*(Huxx(i,j)+Huxx(i+1,j)))
#  endif
          enddo
        enddo
#  undef Huxx
#  undef uxx

#  define vee wrk1
#  define Hvee wrk2
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=jstrV
        else
          jmin=jstrV-1
        endif
        if (NORTHERN_EDGE) then
          jmax=jend
        else
          jmax=jend+1 
        endif
#  else
        jmin=jstrV-1
        jmax=jend+1
#  endif
        do j=jmin,jmax
          do i=istr,iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.*v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
            Hvee(i,j)=FlxV(i,j-1,k)-2.*FlxV(i,j,k)+FlxV(i,j+1,k)
          enddo
        enddo
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr,iend
            vee(i,jstrV-1)=vee(i,jstrV)
            Hvee(i,jstrV-1)=Hvee(i,jstrV)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
            vee(i,jend+1)=vee(i,jend)
            Hvee(i,jend+1)=Hvee(i,jend)
          enddo
        endif
#  endif
        do j=jstrV-1,jend
          do i=istr,iend
#  ifdef UPSTREAM
            cff=FlxV(i,j,k)+FlxV(i,j+1,k)-0.125*( Hvee(i,j  )
     &                                           +Hvee(i,j+1))
            VFe(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
     &                          -gamma*( max(cff,0.)*vee(i,j  )
     &                                  +min(cff,0.)*vee(i,j+1)
     &                                                      )) 
#  else
            VFe(i,j)=0.25*( v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
     &                        -0.125*(vee(i,j)+vee(i,j+1))
     &                  )*( FlxV(i,j,k)+FlxV(i,j+1,k)
     &                      -0.125*(Hvee(i,j)+Hvee(i,j+1)))
#  endif
          enddo
        enddo
#  undef Hvee
#  undef vee

#  define uee wrk1
#  define Hvxx wrk2
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=jstr
        else
          jmin=jstr-1
        endif
        if (NORTHERN_EDGE) then
          jmax=jend
        else
          jmax=jend+1
        endif
#  else
        jmin=jstr-1
        jmax=jend+1
#  endif
        do j=jmin,jmax
          do i=istrU,iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.*u(i,j,k,nrhs)
     &                                 +u(i,j+1,k,nrhs)
          enddo
        enddo
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istrU,iend
            uee(i,jstr-1)=uee(i,jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istrU,iend
            uee(i,jend+1)=uee(i,jend)
          enddo
        endif
#  endif
        do j=jstr,jend+1
          do i=istrU-1,iend
           Hvxx(i,j)=FlxV(i-1,j,k)-2.*FlxV(i,j,k)+FlxV(i+1,j,k)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istrU,iend
#  ifdef UPSTREAM
            cff=FlxV(i,j,k)+FlxV(i-1,j,k)-0.125*( Hvxx(i  ,j)
     &                                           +Hvxx(i-1,j))
            UFe(i,j)=0.25*( cff*(u(i,j,k,nrhs)+u(i,j-1,k,nrhs))
     &                          -gamma*( max(cff,0.)*uee(i,j-1)
     &                                  +min(cff,0.)*uee(i,j  )
     &                                                      ))
#  else
            UFe(i,j)=0.25*( u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
     &                        -0.125*(uee(i,j)+uee(i,j-1))
     &                  )*( FlxV(i,j,k)+FlxV(i-1,j,k)
     &                     -0.125*(Hvxx(i,j)+Hvxx(i-1,j)))
#  endif
          enddo
        enddo
#  undef Hvxx
#  undef uee

#  define vxx wrk1
#  define Huee wrk2
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          imin=istr
        else
          imin=istr-1
        endif
        if (EASTERN_EDGE) then
          imax=iend
        else
          imax=iend+1
        endif
#  else
        imin=istr-1
        imax=iend+1
#  endif
        do j=jstrV,jend
          do i=imin,imax
            vxx(i,j)=v(i-1,j,k,nrhs)-2.*v(i,j,k,nrhs)
     &                                 +v(i+1,j,k,nrhs)
          enddo
        enddo
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstrV,jend
            vxx(istr-1,j)=vxx(istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstrV,jend
            vxx(iend+1,j)=vxx(iend,j)
          enddo
        endif
#  endif
        do j=jstrV-1,jend
          do i=istr,iend+1
           Huee(i,j)=FlxU(i,j-1,k)-2.*FlxU(i,j,k)+FlxU(i,j+1,k)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend+1
#  ifdef UPSTREAM
            cff=FlxU(i,j,k)+FlxU(i,j-1,k)-0.125*( Huee(i,j  )
     &                                           +Huee(i,j-1))
            VFx(i,j)=0.25*( cff*(v(i,j,k,nrhs)+v(i-1,j,k,nrhs))
     &                          -gamma*( max(cff,0.)*vxx(i-1,j)
     &                                  +min(cff,0.)*vxx(i  ,j)
     &                                                      ))
#  else
            VFx(i,j)=0.25*( v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
     &                        -0.125*(vxx(i,j)+vxx(i-1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i,j-1,k)
     &                     -0.125*(Huee(i,j)+Huee(i,j-1)))
#  endif
          enddo
        enddo
#  undef Huee
#  undef vxx
        do j=jstr,jend
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)-UFx(i,j  )+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i  ,j)
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MXadv(i,j,k,1) = -UFx(i,j)+UFx(i-1,j)
              MYadv(i,j,k,1) = -UFe(i,j+1)+UFe(i,j)
            endif
# elif defined DIAGNOSTICS_VRT
            if (nnew.ne.3) then
              if (k.eq.1) then
                wrkXadv(i,j,1) = -UFx(i,j)+UFx(i-1,j)
                wrkYadv(i,j,1) = -UFe(i,j+1)+UFe(i,j)
              else
                wrkXadv(i,j,1) = wrkXadv(i,j,1) - UFx(i,j)+UFx(i-1,j)
                wrkYadv(i,j,1) = wrkYadv(i,j,1) - UFe(i,j+1)+UFe(i,j)
              endif
            endif
#  endif
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j  )
     &                         -VFe(i  ,j)+VFe(i,j-1)
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MXadv(i,j,k,2) = -VFx(i+1,j)+VFx(i,j)
              MYadv(i,j,k,2) = -VFe(i,j)+VFe(i,j-1)
            endif
# elif defined DIAGNOSTICS_VRT
            if (nnew.ne.3) then
              if (k.eq.1) then
                wrkXadv(i,j,2) = -VFx(i+1,j)+VFx(i,j)
                wrkYadv(i,j,2) = -VFe(i,j)+VFe(i,j-1)
              else
                wrkXadv(i,j,2) = wrkXadv(i,j,2) -VFx(i+1,j)+VFx(i,j)
                wrkYadv(i,j,2) = wrkYadv(i,j,2) -VFe(i,j)+VFe(i,j-1)
              endif
            endif
#  endif
          enddo
        enddo
# endif /* UV_ADV */
      enddo
































!
! Dynamic bottom drag coefficient
!
# define rd VFe

      if (COUPLED_STAGE_ONLY) then
        Zob=0.01
        do j=jstrV-1,jend
          do i=istrU-1,iend
# ifdef LINEAR_DRAG_ONLY
            rd(i,j)=rdrg
# else
            cff=sqrt( 0.333333333333*(
     &              u(i,j,1,nrhs)**2 +u(i+1,j,1,nrhs)**2
     &                    +u(i,j,1,nrhs)*u(i+1,j,1,nrhs)
     &              +v(i,j,1,nrhs)**2+v(i,j+1,1,nrhs)**2
     &                    +v(i,j,1,nrhs)*v(i,j+1,1,nrhs)
     &                                               ))

c**         rd(i,j)=rdrg + rdrg2*cff

            rd(i,j)=rdrg + cff*(vonKar/log(Hz(i,j,1)/Zob))**2


# endif
          enddo
        enddo
# if defined MPI  || defined EW_PERIODIC 
        if (WEST_INTER) then          ! Save "rd" into shared array
          imin=istr-1                 ! "r_D" for the subsequent use
        else                          ! in barotropic mode. Note that
          imin=istr                   ! in the case of MPI code the
        endif                         ! range of indices is extended
# else
        imin=istr                     ! into one row of points in MPI
# endif
# if defined MPI  || defined NS_PERIODIC
        if (SOUTH_INTER) then
          jmin=jstr-1                 ! exchange margines so that
        else                          ! "r_D" can be correctly
          jmin=jstr                   ! averaged to U- and V-points
        endif                         ! without exchanging its halo 
# else
        jmin=jstr                     ! points by MPI messages.

# endif
        do j=jmin,jend
          do i=imin,iend
            r_D(i,j)=rd(i,j)
          enddo
        enddo
      endif












      do j=jstr,jend
# ifdef UV_ADV
!
! Compute and add in vertical advection terms:
!
# define SPLINES
# define NEUMANN
 
#  ifdef SPLINES
        do i=istrU,iend
          DC(i,1)=0.5625*(Hz(i  ,j,1)+Hz(i-1,j,1))
     &           -0.0625*(Hz(i+1,j,1)+Hz(i-2,j,1))
#  if defined NEUMANN
          FC(i,0)=1.5*u(i,j,1,nrhs)
          CF(i,1)=0.5
#  elif defined LINEAR_CONTINUATION
          FC(i,0)=2.0*u(i,j,1,nrhs)
          CF(i,1)=1.
#  endif
        enddo
        do k=1,N-1,+1    !--> irreversible
          do i=istrU,iend
            DC(i,k+1)=0.5625*(Hz(i  ,j,k+1)+Hz(i-1,j,k+1))
     &               -0.0625*(Hz(i+1,j,k+1)+Hz(i-2,j,k+1))

            cff=1./(2.*DC(i,k)+DC(i,k+1)*(2.-CF(i,k)))
            CF(i,k+1)=cff*DC(i,k)
            FC(i,k)=cff*( 3.*( DC(i,k  )*u(i,j,k+1,nrhs)
     &                        +DC(i,k+1)*u(i,j,k  ,nrhs))
     &                              -DC(i,k+1)*FC(i,k-1))
          enddo
        enddo               !--> discard DC, keep CF,FC
        do i=istrU,iend
#  if defined NEUMANN
          FC(i,N)=(3.*u(i,j,N,nrhs)-FC(i,N-1))/(2.-CF(i,N))
#  elif defined LINEAR_CONTINUATION
          FC(i,N)=(2.*u(i,j,N,nrhs)-FC(i,N-1))/(1.-CF(i,N))
#  endif
          DC(i,N)=0.        !<-- uppermost W*U flux 
        enddo
        do k=N-1,1,-1       !--> irreversible
          do i=istrU,iend
            FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

            DC(i,k)=FC(i,k)*( 0.5625*(W(i  ,j,k)+W(i-1,j,k))
     &                       -0.0625*(W(i+1,j,k)+W(i-2,j,k)))

            ru(i,j,k+1)=ru(i,j,k+1) -DC(i,k+1)+DC(i,k)
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MVadv(i,j,k+1,1) = -DC(i,k+1)+DC(i,k)
            endif
#  endif
          enddo
        enddo                       !--> discard CF,FC
        do i=istrU,iend
          ru(i,j,1)=ru(i,j,1) -DC(i,1)

#  ifdef DIAGNOSTICS_UV
          if (nnew.ne.3) then
            MVadv(i,j,1,1) =  -DC(i,1)
          endif
#  endif
        enddo                          !--> discard DC

#  else
        do k=2,N-2
          do i=istrU,iend
            FC(i,k)=( 0.5625*(u(i,j,k  ,nrhs)+u(i,j,k+1,nrhs))
     &               -0.0625*(u(i,j,k-1,nrhs)+u(i,j,k+2,nrhs)))
     &                       *( 0.5625*(W(i  ,j,k)+W(i-1,j,k))
     &                         -0.0625*(W(i+1,j,k)+W(i-2,j,k)))
          enddo
        enddo
        do i=istrU,iend
          FC(i,N)=0.
          FC(i,N-1)=( 0.5625*(u(i,j,N-1,nrhs)+u(i,j,N,nrhs))
     &                 -0.0625*(u(i,j,N-2,nrhs)+u(i,j,N,nrhs)))
     &                   *( 0.5625*(W(i  ,j,N-1)+W(i-1,j,N-1))
     &                     -0.0625*(W(i+1,j,N-1)+W(i-2,j,N-1)))
 
          FC(i,  1)=( 0.5625*(u(i,j,  1,nrhs)+u(i,j,2,nrhs))
     &                 -0.0625*(u(i,j,  1,nrhs)+u(i,j,3,nrhs)))
     &                       *( 0.5625*(W(i  ,j,1)+W(i-1,j,1))
     &                         -0.0625*(W(i+1,j,1)+W(i-2,j,1)))
          FC(i,0)=0.
        enddo
c*      do k=1,N-1
c*        do i=istrU,iend
c*          FC(i,k)=0.25*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))
c*     &                            *(W(i,j,k)+W(i-1,j,k))
c*        enddo
c*      enddo
c*      do i=istrU,iend
c*        FC(i,0)=0.
c*        FC(i,N)=0.
c*      enddo
        do k=1,N
          do i=istrU,iend
            ru(i,j,k)=ru(i,j,k)-FC(i,k)+FC(i,k-1)
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MVadv(i,j,k,1) = -FC(i,k)+FC(i,k-1)
            endif
#  endif
          enddo
        enddo               !--> discard FC
#  endif
 
 
        if (j.ge.jstrV) then
#  ifdef SPLINES
          do i=istr,iend
            DC(i,1)=0.5625*(Hz(i  ,j,1)+Hz(i,j-1,1))
     &             -0.0625*(Hz(i,j+1,1)+Hz(i,j-2,1))
#  if defined NEUMANN
            FC(i,0)=1.5*v(i,j,1,nrhs)
            CF(i,1)=0.5
#  elif defined LINEAR_CONTINUATION
            FC(i,0)=2.0*v(i,j,1,nrhs)
            CF(i,1)=1.
#  endif
          enddo
          do k=1,N-1,+1       !--> irreversible
            do i=istr,iend
              DC(i,k+1)=0.5625*(Hz(i  ,j,k+1)+Hz(i,j-1,k+1))
     &                 -0.0625*(Hz(i,j+1,k+1)+Hz(i,j-2,k+1))

              cff=1./(2.*DC(i,k)+DC(i,k+1)*(2.-CF(i,k)))
              CF(i,k+1)=cff*DC(i,k)
              FC(i,k)=cff*( 3.*( DC(i,k  )*v(i,j,k+1,nrhs)
     &                          +DC(i,k+1)*v(i,j,k  ,nrhs))
     &                                -DC(i,k+1)*FC(i,k-1))
            enddo
          enddo               !--> discard DC, keep CF,FC
          do i=istr,iend
#  if defined NEUMANN
            FC(i,N)=(3.*v(i,j,N,nrhs)-FC(i,N-1))/(2.-CF(i,N))
#  elif defined LINEAR_CONTINUATION
            FC(i,N)=(2.*v(i,j,N,nrhs)-FC(i,N-1))/(1.-CF(i,N))
#  endif
            DC(i,N)=0.        !<-- uppermost W*V flux
          enddo
          do k=N-1,1,-1       !--> irreversible
            do i=istr,iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

              DC(i,k)=FC(i,k)*( 0.5625*(W(i,j  ,k)+W(i,j-1,k))
     &                         -0.0625*(W(i,j+1,k)+W(i,j-2,k)))

              rv(i,j,k+1)=rv(i,j,k+1) -DC(i,k+1)+DC(i,k)

#  ifdef DIAGNOSTICS_UV
              if (nnew.ne.3) then
                MVadv(i,j,k+1,2) = -DC(i,k+1)+DC(i,k)
              endif
#  endif

            enddo
          enddo               !--> discard CF,FC

          do i=istr,iend
            rv(i,j,1)=rv(i,j,1) -DC(i,1)
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MVadv(i,j,1,2) = -DC(i,1)
            endif
#  endif

          enddo                         !--> discard DC

#  else
          do k=2,N-2
            do i=istr,iend
              FC(i,k)=( 0.5625*(v(i,j,k ,nrhs)+v(i,j,k+1,nrhs))
     &                 -0.0625*(v(i,j,k-1,nrhs)+v(i,j,k+2,nrhs)))
     &                         *( 0.5625*(W(i,j  ,k)+W(i,j-1,k))
     &                           -0.0625*(W(i,j+1,k)+W(i,j-2,k)))
            enddo
          enddo
          do i=istr,iend
            FC(i,N)=0.
            FC(i,N-1)=(  0.5625*(v(i,j,N-1,nrhs)+v(i,j,N,nrhs))
     &                  -0.0625*(v(i,j,N-2,nrhs)+v(i,j,N,nrhs)))
     &                    *( 0.5625*(W(i,j  ,N-1)+W(i,j-1,N-1))
     &                      -0.0625*(W(i,j+1,N-1)+W(i,j-2,N-1)))
 
            FC(i,  1)=(  0.5625*(v(i,j,  1,nrhs)+v(i,j,2,nrhs))
     &                  -0.0625*(v(i,j,  1,nrhs)+v(i,j,3,nrhs)))
     &                        *( 0.5625*(W(i,j  ,1)+W(i,j-1,1))
     &                          -0.0625*(W(i,j+1,1)+W(i,j-2,1)))
            FC(i,0)=0.
          enddo
c*        do k=1,N-1
c*          do i=istr,iend
c*            FC(i,k)=0.25*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
c*     &                            *(W(i,j,k)+W(i,j-1,k))
c*          enddo
c*        enddo
c*        do i=istr,iend
c*          FC(i,0)=0.
c*          FC(i,N)=0.
c*        enddo
          do k=1,N
            do i=istr,iend
              rv(i,j,k)=rv(i,j,k)-FC(i,k)+FC(i,k-1)
#  ifdef DIAGNOSTICS_UV
            if (nnew.ne.3) then
              MVadv(i,j,k,2) = -FC(i,k)+FC(i,k-1)
            endif
#  endif
            enddo
          enddo
#  endif
        endif
# endif /* UV_ADV */
























!
! Start computation of the forcing terms for the 2D (barotropic mode)
! momentum equations: vertically integrate the just computed r.h.s
! "ru" and "rv". Also, if so prescribed, add in the difference
! between surface (wind) and bottom (drag) stresses. The computation
! of the 2D forcing terms will be finalized in "rhs2d" during the
! first barotropic time step, when the barotropically computed r.h.ss
! "rubar", "rvbar" will be subtracted from the vertically integrated
! (here) "rufrc", "rvfrc".
!
        if (COUPLED_STAGE_ONLY) then

          do i=istrU,iend
           DC(i,0)=dt*0.25*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))

            FC(i,N-1)=dt *2.*(Akv(i,j,N-1)+Akv(i-1,j,N-1))
     &                        /(  Hz(i,j,N  )+Hz(i-1,j,N  )
     &                          +Hz(i,j,N-1)+Hz(i-1,j,N-1))

            cff=1./(0.5*(Hz(i,j,N)+Hz(i-1,j,N))+FC(i,N-1))

            CF(i,N-1)=cff*FC(i,N-1)


            DC(i,N)=cff*( u(i,j,N,nnew) +DC(i,0)*ru(i,j,N)
     &                                      +dt*sustr(i,j))
          enddo
          do k=N-1,2,-1      !--> forward elimination
            do i=istrU,iend
              FC(i,k-1)= dt*2.*(Akv(i,j,k-1)+Akv(i-1,j,k-1))
     &                          /(  Hz(i,j,k  )+Hz(i-1,j,k  )
     &                            +Hz(i,j,k-1)+Hz(i-1,j,k-1))

              cff=1./( 0.5*(Hz(i,j,k)+Hz(i-1,j,k)) +FC(i,k-1)
     &                                   +FC(i,k)*(1.-CF(i,k))
     &                                                       )
              CF(i,k-1)=cff*FC(i,k-1)
              DC(i,k)=cff*( u(i,j,k,nnew) +DC(i,0)*ru(i,j,k)
     &                                    +FC(i,k)*DC(i,k+1))
            enddo
          enddo
          do i=istrU,iend
             DC(i,1)=(u(i,j,1,nnew) +DC(i,0)*ru(i,j,1)
     &                                +FC(i,1)*DC(i,2))
     &                   /( 0.5*(Hz(i,j,1)+Hz(i-1,j,1))
# ifdef NO_SLIP_BOTTOM
     &                    +dt * 0.5*(rd(i,j)+rd(i-1,j))
# endif
     &                         +FC(i,1)*(1.-CF(i,1)) )
          enddo
          cff=1./dt
          do k=2,N,+1          !--> backsubstitution
            do i=istrU,iend
              DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)
!>
              FC(i,k-1)=cff*FC(i,k-1)*(DC(i,k)-DC(i,k-1))
            enddo
          enddo       !--> now FC is visc. flux
          do i=istrU,iend
            DC(i,0)=dm_u(i,j)*dn_u(i,j)
            ru(i,j,N)=ru(i,j,N) +DC(i,0)*(sustr(i,j)-FC(i,N-1))
            ru(i,j,1)=ru(i,j,1) +DC(i,0)*( FC(i,1)
     &                        -0.5*(rd(i-1,j)+rd(i,j))*DC(i,1) )

            rufrc(i,j)=ru(i,j,N)+ru(i,j,1)


# ifdef DIAGNOSTICS_UV
                MVmix(i,j,1,1)= DC(i,0)*( FC(i,1)
     &                        -0.5*(rd(i-1,j)+rd(i,j))*DC(i,1) )
                MVmix(i,j,N,1)= DC(i,0)*(sustr(i,j)-FC(i,N-1))
# elif defined DIAGNOSTICS_VRT 
                wrkVmix(i,j,1)= DC(i,0)*( FC(i,1)
     &                        -0.5*(rd(i-1,j)+rd(i,j))*DC(i,1) )
     &                        + DC(i,0)*(sustr(i,j)-FC(i,N-1))
# endif /* DIAGNOSTICS_UV */


# if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_UV
               pvForc_bot(i,j,1) =  DC(i,0)*( FC(i,1)
     &                        -0.5*(rd(i-1,j)+rd(i,j))*DC(i,1) )

               pvForc_sfc(i,j,1) = DC(i,0)*(sustr(i,j)-FC(i,N-1))
# endif /* DIAGNOSTICS_PV */



# if defined DIAGNOSTICS_PV_SUP

                pvForc_0(i,j,1) =  DC(i,0)*(
     &                        -0.5*(rd(i-1,j)+rd(i,j))*DC(i,1) )

                pvForc_1(i,j,1) =  DC(i,0)* FC(i,1)

                pvForc_2(i,j,1) =  -FC(i,N-1)*DC(i,0)

                pvForc_3(i,j,1) =  DC(i,0)*(sustr(i,j))

# endif /* DIAGNOSTICS_PV_SUP */





# ifdef DIAGNOSTICS_VRT
                wrkDrag(i,j,1)= -0.5*(rd(i-1,j)+rd(i,j))*DC(i,1)
                wrkWind(i,j,1)= sustr(i,j)
# endif /* DIAGNOSTICS_VRT */


          enddo
          do k=N-1,2,-1
            do i=istrU,iend
              ru(i,j,k)=ru(i,j,k) +DC(i,0)*(FC(i,k)-FC(i,k-1))

              rufrc(i,j)=rufrc(i,j)+ru(i,j,k)

# ifdef DIAGNOSTICS_UV
                MVmix(i,j,k,1)=DC(i,0)*(FC(i,k)-FC(i,k-1))
# elif defined DIAGNOSTICS_VRT 
                wrkVmix(i,j,1)= wrkVmix(i,j,1)
     &                        + DC(i,0)*(FC(i,k)-FC(i,k-1))
# endif /* DIAGNOSTICS_UV */

            enddo
          enddo






















          if (j.ge.jstrV) then

            do i=istr,iend
              DC(i,0)=dt*0.25*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

              FC(i,N-1)=dt *2.*(Akv(i,j,N-1)+Akv(i,j-1,N-1))
     &                          /(  Hz(i,j,N  )+Hz(i,j-1,N  )
     &                            +Hz(i,j,N-1)+Hz(i,j-1,N-1))

              cff=1./(0.5*(Hz(i,j,N)+Hz(i,j-1,N))+FC(i,N-1))

              CF(i,N-1)=cff*FC(i,N-1)
              DC(i,N)=cff*( v(i,j,N,nnew) +DC(i,0)*rv(i,j,N)
     &                                        +dt*svstr(i,j))
            enddo
            do k=N-1,2,-1      !--> forward elimination
              do i=istr,iend
                FC(i,k-1)= dt*2.*(Akv(i,j,k-1)+Akv(i,j-1,k-1))
     &                            /(  Hz(i,j,k  )+Hz(i,j-1,k  )
     &                              +Hz(i,j,k-1)+Hz(i,j-1,k-1))

                cff=1./( 0.5*(Hz(i,j,k)+Hz(i,j-1,k)) +FC(i,k-1)
     &                                    +FC(i,k)*(1.-CF(i,k))
     &                                                        )
                CF(i,k-1)=cff*FC(i,k-1)
                DC(i,k)=cff*( v(i,j,k,nnew) +DC(i,0)*rv(i,j,k)
     &                                      +FC(i,k)*DC(i,k+1))
              enddo
            enddo
            do i=istr,iend
               DC(i,1)=( v(i,j,1,nnew) +DC(i,0)*rv(i,j,1)
     &                                   +FC(i,1)*DC(i,2))
     &                      /( 0.5*(Hz(i,j,1)+Hz(i,j-1,1))
# ifdef NO_SLIP_BOTTOM
     &                         +dt*0.5*(rd(i,j)+rd(i,j-1))
# endif
     &                            +FC(i,1)*(1.-CF(i,1)) )
            enddo
            cff=1./dt
            do k=2,N,+1          !<-- backsubstitution
              do i=istr,iend
                DC(i,k)=DC(i,k) +CF(i,k-1)*DC(i,k-1)
!>
                FC(i,k-1)=cff*FC(i,k-1)*(DC(i,k)-DC(i,k-1))
              enddo
            enddo        !--> now FC is visc. flux

            do i=istr,iend
              DC(i,0)=dm_v(i,j)*dn_v(i,j)
              rv(i,j,N)=rv(i,j,N) +DC(i,0)*(svstr(i,j)-FC(i,N-1))
              rv(i,j,1)=rv(i,j,1) +DC(i,0)*( FC(i,1)
     &                          -0.5*(rd(i,j-1)+rd(i,j))*DC(i,1) )
              rvfrc(i,j)=rv(i,j,N)+rv(i,j,1)

# ifdef DIAGNOSTICS_UV
                MVmix(i,j,1,2)= DC(i,0)*( FC(i,1)
     &                          -0.5*(rd(i,j-1)+rd(i,j))*DC(i,1) )
                MVmix(i,j,N,2)= DC(i,0)*(svstr(i,j)-FC(i,N-1))
# elif defined DIAGNOSTICS_VRT 
                wrkVmix(i,j,2)= DC(i,0)*( FC(i,1)
     &                        -0.5*(rd(i,j-1)+rd(i,j))*DC(i,1) )
     &                        + DC(i,0)*(svstr(i,j)-FC(i,N-1))
# endif /* DIAGNOSTICS_UV */

# if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_UV
              pvForc_bot(i,j,2) =  DC(i,0)*( FC(i,1)
     &                          -0.5*(rd(i,j-1)+rd(i,j))*DC(i,1) )

              pvForc_sfc(i,j,2) =  DC(i,0)*(svstr(i,j)-FC(i,N-1))
# endif /* DIAGNOSTICS_PV */

# if defined DIAGNOSTICS_PV_SUP

                pvForc_0(i,j,2) =  DC(i,0)*(
     &                          -0.5*(rd(i,j-1)+rd(i,j))*DC(i,1) )

                pvForc_1(i,j,2) =  DC(i,0)* FC(i,1)

                pvForc_2(i,j,2) =  -FC(i,N-1)*DC(i,0)

                pvForc_3(i,j,2) =  DC(i,0)*(svstr(i,j))

# endif /* DIAGNOSTICS_PV_SUP */

# ifdef DIAGNOSTICS_VRT
                wrkDrag(i,j,2)= -0.5*(rd(i,j-1)+rd(i,j))*DC(i,1)
                wrkWind(i,j,2)= svstr(i,j)

# endif /* DIAGNOSTICS_VRT */

            enddo
            do k=N-1,2,-1
              do i=istr,iend
                rv(i,j,k)=rv(i,j,k) +DC(i,0)*(FC(i,k)-FC(i,k-1))

                rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k)

# ifdef DIAGNOSTICS_UV
                MVmix(i,j,k,2)=DC(i,0)*(FC(i,k)-FC(i,k-1))
# elif defined DIAGNOSTICS_VRT
                wrkVmix(i,j,2)= wrkVmix(i,j,2)
     &                        + DC(i,0)*(FC(i,k)-FC(i,k-1))
# endif /* DIAGNOSTICS_UV */

              enddo
            enddo
          endif

        endif
# undef rd

      enddo      !<-- j




      return
      end


      subroutine check_rhs3d_switches (ierr)
!
! This routine keeps track of the status of local CPP-settings
! in "rhs3d35S.F". This code must be placed here rather than in a
! separate file in order to be exposed to the relevant CPP-settings.
! It does not affect any model results, other than signature in
! global attribute "CPPS" in output netCDF files.
!
      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2
      ie=is+30
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='<rhs3d35S.F>'
      MPI_master_only write(stdout,'(1x,A)') cpps(is:ie)
      is=ie+2
# ifdef UPSTREAM
      ie=is+7
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='UPSTREAM'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
      return
  99  MPI_master_only write(stdout,'(/1x,2A/12x,A/)') '### ERROR: ',
     &  'Unsufficient lenght of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      return
      end

#else
      subroutine rhs3d_empty
      end
#endif /* SOLVE3D */
 
 
