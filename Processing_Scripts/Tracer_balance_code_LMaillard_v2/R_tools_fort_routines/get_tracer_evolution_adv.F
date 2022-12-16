
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_tracer_evolution_adv (Lm,Mm,N,u,v, z_r,z_w
     &  ,z_r_bak,z_w_bak
     &  ,pm,pn,dt
     &  ,t, W
     &  ,rmask,umask,vmask
     &  ,TXadv,TYadv,TVadv)    


      integer Lm,Mm,N,NT, i,j,k
     &      ,istr,iend,jstr,jend
     &      ,imin,imax,jmin,jmax
     &      itemp,isalt,dt,itrc

      parameter (NT=1)  
      parameter (itemp=1,isalt=2)  
      parameter (epsil=1.E-16)

      ! INPUTS
      real*8 t(0:Lm+1,0:Mm+1,N,NT)

      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N)
     &      ,z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N)
     &      ,z_r_bak(0:Lm+1,0:Mm+1,N),  z_w_bak(0:Lm+1,0:Mm+1,0:N)
     &      ,pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1)
     &      ,W(0:Lm+1,0:Mm+1,0:N)
     &      ,umask(1:Lm+1,0:Mm+1),vmask(0:Lm+1,1:Mm+1)
     &      ,rmask(0:Lm+1,0:Mm+1)

      ! OUTPUTS
      real*8 TXadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TYadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TVadv(0:Lm+1,0:Mm+1,N,NT)


      ! LOCAL
      real*8 wrk1(0:Lm+1,0:Mm+1), wrk2(0:Lm+1,0:Mm+1)
     &      ,We(0:Lm+1,0:Mm+1,0:N)
     &      ,Cu_adv(0:Lm+1,0:N)
     &      ,wrk3(0:Lm+1)
     &      ,wrk4(0:Lm+1)
     &      ,FX(0:Lm+1,0:Mm+1), FE(0:Lm+1,0:Mm+1)
     &      ,WORK(0:Lm+1,0:Mm+1)
     &      ,FC(0:Lm+1,0:N), DC(0:Lm+1,0:N)
     &      ,CF(0:Lm+1,0:N)
     &      ,Hz(0:Lm+1,0:Mm+1,N)
     &      ,Hz_bak(0:Lm+1,0:Mm+1,N)
     &      ,on_u(0:Lm+1,0:Mm+1), om_v(0:Lm+1,0:Mm+1)
     &      ,on_r(0:Lm+1,0:Mm+1), om_r(0:Lm+1,0:Mm+1)
     &      ,FlxU(1:Lm+1,0:Mm+1,N), FlxV(0:Lm+1,1:Mm+1,N)
     &      ,tnew(0:Lm+1,0:Mm+1,N,NT)
     &      ,c_min,c_max,cu_cut,fcu
     &      ,cff,my_cu_adv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,z_r_bak,z_w_bak,pm,pn,dt,t,W
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
        
        c_min=0.6
        c_max=0.95D0
        cu_cut=2.D0*c_max-c_min
        fcu=4.D0*c_max*(c_max-c_min)


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
           Hz_bak(i,j,k) = z_w_bak(i,j,k) - z_w_bak(i,j,k-1)
          enddo
        enddo
      enddo

      do j=jmin,jmax
        do i=imin+1,imax
          on_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
          do k=1,N,+1
           FlxU(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)*u(i,j,k)
          enddo
        enddo
      enddo

      do j=jmin+1,jmax
        do i=imin,imax
          om_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
          do k=1,N,+1
           FlxV(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)*v(i,j,k)
          enddo
        enddo
      enddo

      do j=jstr,jend
        do i=istr,iend
          om_r(i,j)=1./pm(i,j)
          on_r(i,j)=1./pn(i,j)
        enddo
      enddo

      do j=jstr,jend
        do i=istr,iend
          We(i,j,0)=0.D0
          wrk3(i)=0.D0
        enddo
        do k=1,N,+1
          do i=istr,iend
            We(i,j,k)=We(i,j,k-1) -FlxU(i+1,j,k)+FlxU(i,j,k)
     &                            -FlxV(i,j+1,k)+FlxV(i,j,k)
     &                       -((Hz(i,j,k)-Hz_bak(i,j,k))/dt)
     &                                  *om_r(i,j)*on_r(i,j)
            wrk3(i)=wrk3(i)+Hz(i,j,k)
            Cu_adv(i,k)=max(FlxU(i+1,j,k),0.D0)-min(FlxU(i,j,k),0.D0)
     &                 +max(FlxV(i,j+1,k),0.D0)-min(FlxV(i,j,k),0.D0)
          enddo
        enddo
        do i=istr,iend
          wrk4(i)=0.D0
          wrk3(i)=We(i,j,N)/wrk3(i)
        enddo
        do k=1,N-1,+1
          do i=istr,iend
            wrk4(i)=wrk4(i)+Hz(i,j,k)
            We(i,j,k)=We(i,j,k)-wrk3(i)*wrk4(i)
          enddo
        enddo
        do i=istr,iend
           We(i,j,N)=0.D0
        enddo
        do k=1,N
          do i=istr,iend
            cff = dt*pm(i,j)*pn(i,j)/Hz(i,j,k)
            Cu_adv(i,k)=cff*(
     &                   Cu_adv(i,k)
     &                 + max(We(i,j,k),0.D0)-min(We(i,j,k-1),0.D0)
     &                       )
          enddo
        enddo
        do k=1,N-1
          do i=istr,iend
             my_cu_adv = max(Cu_adv(i,k+1),Cu_adv(i,k))
             if(my_cu_adv.lt.c_min) then
                 cff = 0.D0
             elseif(my_cu_adv.lt.cu_cut) then
                 cff = (my_cu_adv-c_min)**2
                 cff = cff / (fcu+cff)
             else
                 cff = (my_cu_adv-c_max)/my_cu_adv
             endif
             We(i,j,k) = (1.D0-cff)*We(i,j,k)
          enddo
        enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! then using UPSTREAM scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute and add in contribution due to horizontal advection fluxes.
      do itrc=1,NT
        do k=1,N
!# include "compute_horiz_tracer_fluxes.h"
# include "compute_horiz_tracer_fluxes_croco.h"
          do j=jstr,jend
            do i=istr,iend
              TXadv(i,j,k,itrc) = -(FX(i+1,j)-FX(i,j)) 
#  ifdef MASKING
     &         * rmask(i,j)
#  endif
              TYadv(i,j,k,itrc) = -(FE(i,j+1)-FE(i,j))
#  ifdef MASKING
     &         * rmask(i,j)
#  endif
            enddo
          enddo                 !-->  discard FX,FE
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jstr,jend
        do itrc=1,NT
 

!# include "compute_vert_tracer_fluxes.h"
# include "compute_vert_tracer_fluxes_croco.h"
 
          do k=1,N            ! Apply vertical advective fluxes.
            do i=istr,iend
              TVadv(i,j,k,itrc) =  -(FC(i,k)-FC(i,k-1))
#  ifdef MASKING
     &         * rmask(i,j)
#  endif

            enddo
          enddo               !--> discard FC
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! put in advection instead of flux : udx(tracer) instead of dx(u*tracer)
      do itrc=1,NT
        do k=1,N
          do j=jstr,jend
            do i=istr,iend
              TXadv(i,j,k,itrc)=TXadv(i,j,k,itrc)+
     &           t(i,j,k,itrc)*(FlxU(i+1,j,k)-FlxU(i,j,k))
              TYadv(i,j,k,itrc)=TYadv(i,j,k,itrc)+
     &           t(i,j,k,itrc)*(FlxV(i,j+1,k)-FlxV(i,j,k))
              TVadv(i,j,k,itrc)=TVadv(i,j,k,itrc)-
     &           t(i,j,k,itrc)*(FlxU(i+1,j,k)-FlxU(i,j,k)+
     &                            FlxV(i,j+1,k)-FlxV(i,j,k))
            enddo
          enddo           
        enddo
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set lateral boundary conditions; nudge toward tracer climatology;
! apply land-sea mask and exchange periodic boundary conditions.
!

!---------------------------------------------------------------
! Divide all diagnostic terms by the cell volume
! (Hz(i,j,k,itrc)/(pm(i,j).*pn(i,j)). There after the unit
! of diagnostic terms will be: (unit of tracers)* s-1.
!---------------------------------------------------------------
      do itrc=1,NT
       do k=1,N
         do j=jstr,jend
           do i=istr,iend
              cff=pm(i,j)*pn(i,j)/Hz(i,j,k)
              TXadv(i,j,k,itrc)=TXadv(i,j,k,itrc)*cff
              TYadv(i,j,k,itrc)=TYadv(i,j,k,itrc)*cff
              TVadv(i,j,k,itrc)=TVadv(i,j,k,itrc)*cff
           enddo
         enddo
        enddo
      enddo


      return
      end
