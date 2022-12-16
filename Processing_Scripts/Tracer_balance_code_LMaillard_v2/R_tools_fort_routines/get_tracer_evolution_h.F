
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_tracer_evolution_h (Lm,Mm,N,u,v, z_r,z_w
     &  ,pm,pn,dt
     &  ,t, W
     &  ,TXadv,TYadv,TVadv)    


      integer Lm,Mm,N,NT, i,j,k
     &      ,istr,iend,jstr,jend
     &      ,imin,imax,jmin,jmax
     &      itemp,isalt,dt,itrc

      parameter (NT=2)  
      parameter (itemp=1,isalt=2)  
      parameter (epsil=1.E-16)

      ! INPUTS
      real*8 t(0:Lm+1,0:Mm+1,N,NT)

      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N)
     &      ,z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N)
     &      ,pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1)
     &      ,W(0:Lm+1,0:Mm+1,0:N)

      ! OUTPUTS
      real*8 TXadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TYadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TVadv(0:Lm+1,0:Mm+1,N+1,NT)


      ! LOCAL
      real*8 wrk1(0:Lm+1,0:Mm+1), wrk2(0:Lm+1,0:Mm+1)
     &      ,FX(0:Lm+1,0:Mm+1), FE(0:Lm+1,0:Mm+1)
     &      ,WORK(0:Lm+1,0:Mm+1)
     &      ,FC(0:Lm+1,0:N), DC(0:Lm+1,0:N)
     &      ,CF(0:Lm+1,0:N)
     &      ,Hz(0:Lm+1,0:Mm+1,N)
     &      ,dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1)
     &      ,FlxU(1:Lm+1,0:Mm+1,N),  FlxV(0:Lm+1,1:Mm+1,N)
     &      ,tnew(0:Lm+1,0:Mm+1,N,NT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn,dt,t,stflx, srflx, ghat, swr_frac,W
Cf2py intent(out) TXadv,TYadv,TVadv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        istr=1
        iend=Lm
        jstr=1
        jend=Mm


        imin=0
        imax=Lm+1
        jmin=0
        jmax=Mm+1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j=jmin,jmax
        do i=imin,imax
          do k=0,N,+1
              W(i,j,k)=W(i,j,k)/(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo


      do j=jmin,jmax
        do i=imin,imax
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
        enddo
      enddo



      do j=jmin,jmax
        do i=imin+1,imax
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
            do k=1,N,+1
              FlxU(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
     &      * u(i,j,k) 
            enddo
          enddo
      enddo



      do j=jmin+1,jmax
        do i=imin,imax
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            do k=1,N,+1
              FlxV(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
     &      * v(i,j,k) 
            enddo
          enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! This part of the code is valid only if the ROMS simulations has been run with UPSTREAM scheme
! So it is not valid if ADV_ISONEUTRAL was used
!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! then using UPSTREAM scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! in the code
! Compute and add in contribution due to horizontal advection fluxes.
! Note that at this moment the content of array t(:,:,:,nnew,itrc) is
! concentration of tracer "itrc" taken at time corresponding to "nstp"
! and already MULTIPLIED by Hz(:,:,:) at "nstp" in "pre_step3d_tile".
! --> already multiplied 
      do itrc=1,NT
        do k=1,N
# include "compute_horiz_tracer_fluxes.h"
!# include "compute_horiz_tracer_fluxes_simplified.h"

          do j=jstr,jend
            do i=istr,iend
              TXadv(i,j,k,itrc) = FX(i,j)
              TYadv(i,j,k,itrc) = FE(i,j)
            enddo
          enddo                 !-->  discard FX,FE
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jstr,jend
        do itrc=1,NT
 

# include "compute_vert_tracer_fluxes.h"
 
          do k=0,N            ! Apply vertical advective fluxes.
            do i=istr,iend

              TVadv(i,j,k,itrc) =  FC(i,k)
     &                  * pm(i,j) * pn(i,j)

            enddo
          enddo               !--> discard FC
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      return
      end
