! This is "compute_vert_tracer_fluxes.h" -- module which computes
! vertical fluxes for tracer equations. In the case of SPLINE_TS two
! versions of top and bottom boundary conditions are supported:
! Neumann (setting first derivative to zero at the top and bottom
! boundaries) and LINEAR_TS CONTINUATION (assumption that the tracer
! distributions are linear within the top-most and botom-most grid
! boxes).


c--#define SPLINE_TS
#define AKIMA_V
c--# define NEUMANN_TS
# define LINEAR_TS



# ifdef SPLINE_TS
          do i=istr,iend
#  if defined NEUMANN_TS
            FC(i,0)=1.5*t(i,j,1,itrc)
            CF(i,1)=0.5
#  elif defined LINEAR_TS
            FC(i,0)=2.0*t(i,j,1,itrc)
            CF(i,1)=1.
#  endif
          enddo
          do k=1,N-1,+1    !--> recursive
            do i=istr,iend
              cff=1./(2.*Hz(i,j,k)+Hz(i,j,k+1)*(2.-CF(i,k)))
              CF(i,k+1)=cff*Hz(i,j,k)
              FC(i,k)=cff*( 3.*( Hz(i,j,k  )*t(i,j,k+1,itrc)
     &                          +Hz(i,j,k+1)*t(i,j,k  ,itrc))
     &                                     -Hz(i,j,k+1)*FC(i,k-1))
            enddo
          enddo
          do i=istr,iend
#  if defined NEUMANN_TS
            FC(i,N)=(3.*t(i,j,N,itrc)-FC(i,N-1))/(2.-CF(i,N))
#  elif defined LINEAR_TS
            FC(i,N)=(2.*t(i,j,N,itrc)-FC(i,N-1))/(1.-CF(i,N))
#  endif
          enddo
          do k=N-1,0,-1    !<-- recursive
            do i=istr,iend
              FC(i,k)=FC(i,k)-CF(i,k+1)*FC(i,k+1)

              FC(i,k+1)=FC(i,k+1)*W(i,j,k+1)  !<-- Convert interface
            enddo                              !    value into vertical
          enddo              !--> discard CF   !    fluxe.
          do i=istr,iend
            FC(i,N)=0.                         ! Set top and bottom
            FC(i,0)=0.                         ! boundary conditions.
          enddo
# elif defined AKIMA_V
          do k=1,N-1
            do i=istr,iend
              FC(i,k)=t(i,j,k+1,itrc)-t(i,j,k,itrc)
            enddo
          enddo
          do i=istr,iend
            FC(i,0)=FC(i,1)
            FC(i,N)=FC(i,N-1)
          enddo
          do k=1,N
            do i=istr,iend
              cff=2.*FC(i,k)*FC(i,k-1)
              if (cff.gt.epsil) then
                CF(i,k)=cff/(FC(i,k)+FC(i,k-1))
              else
                CF(i,k)=0.
              endif
            enddo
          enddo            !--> discard FC
          do k=1,N-1
            do i=istr,iend
              FC(i,k)=0.5*( t(i,j,k,itrc)+t(i,j,k+1,itrc)
     &               -0.333333333333*(CF(i,k+1)-CF(i,k)) )*W(i,j,k)
            enddo
          enddo            !--> discard CF
          do i=istr,iend
            FC(i,0)=0.
            FC(i,N)=0.
          enddo
# else
          do k=2,N-2
            do i=istr,iend
              FC(i,k)=W(i,j,k)*(
     &                     0.58333333333333*( t(i,j,k  ,itrc)
     &                                       +t(i,j,k+1,itrc))
     &                    -0.08333333333333*( t(i,j,k-1,itrc)
     &                                       +t(i,j,k+2,itrc))
     &                                                            )
            enddo
          enddo
          do i=istr,iend
            FC(i, 0)=0.0
            FC(i,  1)=W(i,j,  1)*(     0.5*t(i,j,  1,itrc)
     &                       +0.58333333333333*t(i,j,  2,itrc)
     &                       -0.08333333333333*t(i,j,  3,itrc)
     &                                                            )
            FC(i,N-1)=W(i,j,N-1)*(     0.5*t(i,j,N  ,itrc)
     &                       +0.58333333333333*t(i,j,N-1,itrc)
     &                       -0.08333333333333*t(i,j,N-2,itrc)
     &                                                            )
            FC(i,N )=0.0
          enddo
# endif

c**       do k=1,N-1
c**         do i=istr,iend
c**           FC(i,k)=0.5*(t(i,j,k,itrc)+t(i,j,k+1,itrc))
c**     &                                                *W(i,j,k)
c**         enddo
c**       enddo
c**       do i=istr,iend
c**         FC(i, 0)=0.
c**         FC(i,N )=0.
c**       enddo




#ifdef BIO_1ST_USTREAM_TEST
        endif  !<-- itrc.gt.isalt, bio-components only.
#endif
