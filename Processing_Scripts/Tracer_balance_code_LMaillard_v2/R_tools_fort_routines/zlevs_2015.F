!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for 2015_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine zlevs_2015(Lm,Mm,N, h,zeta, hc, Cs_r, Cs_w,z_r,z_w)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), Cs_r(N),
     &       hc, 
     &       zeta(0:Lm+1,0:Mm+1),
     &       z_r(0:Lm+1,0:Mm+1,N),z_w(0:Lm+1,0:Mm+1,0:N),
     &       h(0:Lm+1,0:Mm+1),
     &       hinv_r(0:Lm+1,0:Mm+1), hinv_w(0:Lm+1,0:Mm+1)
      real*8 cff_r,cff1_r, cff_w,cff1_w, ds,Gcoord_w,Gcoord_r,s_r, s_w
      real, parameter :: eps=1.E-20

     
Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w, Cs_r
Cf2py intent(out) z_r,z_w


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      ds=1.D0/dble(N)

      do j=jmin,jmax
        do i=imin,imax
          z_w(i,j,0)=-h(i,j)
        enddo
        
        do k=1,N,+1  !--> irreversible because of recursion in Hz 
! first define s
          s_w = ds* dble(k-N)
          s_r = ds* (dble(k-N)-0.5)
! first define the Gcoord parameter
          Gcoord_w= 0.16 + abs(s_w)**0.8*(1+s_w)**0.2 + (0.3 *  exp(-1 /
     &     max(abs(s_w), eps)) / max(abs(s_w), eps)**1.7)


          Gcoord_r= 0.16 + abs(s_r)**0.8*(1+s_r)**0.2 + (0.3 *  exp(-1 /
     &     max(abs(s_r), eps)) / max(abs(s_r), eps)**1.7)

          cff_w=hc*s_w*Gcoord_w
          cff_r=hc*s_r*Gcoord_r

          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)


          do i=imin,imax

            hinv_w(i,j)=1./(h(i,j)+hc*Gcoord_w)
            z_w(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
     &                        *(cff_w+cff1_w*h(i,j))*hinv_w(i,j)

            hinv_r(i,j)=1./(h(i,j)+hc*Gcoord_r)
            z_r(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
     &                        *(cff_r+cff1_r*h(i,j))*hinv_r(i,j)

          enddo

        enddo
      enddo
      end


