! Compute horizontal fluxes for tracers.  Essentially interpolate
! tracer values from their native locations on C grid to horizontal
! velocity points with simultaneous translation from grid-box-averages
! to instantaneous values at interface location.  Three options an be
! selected: 3-point upstream-biased parabolic interpolation (UPSTREAM_TS);
!  4-point symmetric fourth-order method (undefined state of both CPP
! switches); 4-point scheme where arithmetic averaging of elementary
! differences is replaced by harmonic averaging (AKIMA), resulting in
! mid-point values be bounded by two nearest values at native location,
! regardless of grid-scale roughness of the interpolated field, while
! still retaining asymptotic fourth-order behavior for smooth fields.
! This code is extracted into a special module because it is used
! twice, in predictor and corrector substeps for tracer variables.
#define UPSTREAM_TS 


#ifdef UPSTREAM_TS
# define curv wrk1
#else
# define grad wrk1
#endif
#ifndef EW_PERIODIC
          if (WESTERN_EDGE) then       ! Determine extended index
            imin=istr                  ! range for computation of
          else                         ! elementary differences: it
            imin=istr-1                ! needs to be restricted
          endif                        ! because in the vicinity of
          if (EASTERN_EDGE) then       ! physical boundary the extra
            imax=iend                  ! point may be not available,
          else                         ! and extrapolation of slope
            imax=iend+1                ! is used instead.
          endif
#else
          imin=istr-1
          imax=iend+1
#endif
          do j=jstr,jend
            do i=imin,imax+1
              FX(i,j)=(t(i,j,k,itrc)-t(i-1,j,k,itrc))
            enddo
          enddo
#ifndef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=jstr,jend
              FX(istr-1,j)=FX(istr,j)
            enddo
          endif
          if (EASTERN_EDGE) then
            do j=jstr,jend
              FX(iend+2,j)=FX(iend+1,j)
            enddo
          endif
#endif
          do j=jstr,jend
            do i=istr-1,iend+1
#if defined UPSTREAM_TS
              curv(i,j)=FX(i+1,j)-FX(i,j)
#elif defined AKIMA
              cff=2.*FX(i+1,j)*FX(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              else
                grad(i,j)=0.
              endif
#else
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
#endif
            enddo
          enddo             !--> discard FX
          do j=jstr,jend
            do i=istr,iend+1
#ifdef UPSTREAM_TS
              FX(i,j)=0.5*(t(i,j,k,itrc)+t(i-1,j,k,itrc))
     &                                                  *FlxU(i,j,k)
     &          -0.1666666666666666*( curv(i-1,j)*max(FlxU(i,j,k),0.)
     &                               +curv(i  ,j)*min(FlxU(i,j,k),0.))
# ifdef DIAGNOSTICS_TS
              TruncFX(i,j)=0.04166666666666667*(curv(i,j)-curv(i-1,j))
     &                                               *abs(FlxU(i,j,k))
# endif
#else
              FX(i,j)=0.5*( t(i,j,k,itrc)+t(i-1,j,k,itrc)
     &                  -0.3333333333333333*( grad(i,j)-grad(i-1,j))
     &                                                 )*FlxU(i,j,k)
#endif
            enddo           !--> discard curv,grad, keep FX
          enddo

#ifndef NS_PERIODIC
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
#else
          jmin=jstr-1
          jmax=jend+1
#endif
          do j=jmin,jmax+1
            do i=istr,iend
              FE(i,j)=(t(i,j,k,itrc)-t(i,j-1,k,itrc))
            enddo
          enddo
#ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=istr,iend
              FE(i,jstr-1)=FE(i,jstr)
            enddo
          endif
          if (NORTHERN_EDGE) then
            do i=istr,iend
              FE(i,jend+2)=FE(i,jend+1)
            enddo
          endif
#endif
          do j=jstr-1,jend+1
            do i=istr,iend
#if defined UPSTREAM_TS
              curv(i,j)=FE(i,j+1)-FE(i,j)
#elif defined AKIMA
              cff=2.*FE(i,j+1)*FE(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              else
                grad(i,j)=0.
              endif
#else
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
#endif
            enddo
          enddo            !--> discard FE

          do j=jstr,jend+1
            do i=istr,iend
#ifdef UPSTREAM_TS
              FE(i,j)=0.5*(t(i,j,k,itrc)+t(i,j-1,k,itrc))
     &                                                  *FlxV(i,j,k)
     &          -0.1666666666666666*( curv(i,j-1)*max(FlxV(i,j,k),0.)
     &                               +curv(i,j  )*min(FlxV(i,j,k),0.))
# ifdef DIAGNOSTICS_TS
              TruncFE(i,j)=0.04166666666666667*(curv(i,j)-curv(i-1,j))
     &                                               *abs(FlxV(i,j,k))
# endif
#else
              FE(i,j)=0.5*( t(i,j,k,itrc)+t(i,j-1,k,itrc)
     &                   -0.3333333333333333*(grad(i,j)-grad(i,j-1))
     &                                                 )*FlxV(i,j,k)
#endif
            enddo
          enddo             !--> discard curv,grad, keep FE

#ifdef PSOURCE
          do is=1,Nsrc             ! Set tracer fluxes due to point
            i=Isrc(is)             ! sources to simulate river run off
            j=Jsrc(is)
            if ( istr.le.i .and. i.le.iend+1 .and.
     &           jstr.le.j .and. j.le.jend+1 ) then
              if (Dsrc(is).eq.0) then
                FX(i,j)=FlxU(i,j,k)*Tsrc(is,k,itrc)
              else
                FE(i,j)=FlxV(i,j,k)*Tsrc(is,k,itrc)
              endif
            endif
          enddo
#endif
#ifdef BIO_1ST_USTREAM_TEST
        endif  !<-- itrc.gt.isalt, bio-components only.
#endif
