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

# define curv wrk1 
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

          do j=jstr,jend
            do i=imin,imax+1
              FX(i,j)=(t(i,j,k,itrc)-t(i-1,j,k,itrc))
            enddo
          enddo

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

          do j=jstr,jend
            do i=istr-1,iend+1
                curv(i,j)=FX(i+1,j)-FX(i,j)
            enddo
          enddo  

          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=0.5*(t(i,j,k,itrc)+t(i-1,j,k,itrc))
     &                                                  *FlxU(i,j,k)
     &          -0.1666666666666666*( curv(i-1,j)*max(FlxU(i,j,k),0.)
     &                               +curv(i  ,j)*min(FlxU(i,j,k),0.))
            enddo           !--> discard curv,grad, keep FX
          enddo

          do j=jmin,jmax+1
            do i=istr,iend
              FE(i,j)=(t(i,j,k,itrc)-t(i,j-1,k,itrc))
            enddo
          enddo
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
          do j=jstr-1,jend+1
            do i=istr,iend
                 curv(i,j)=FE(i,j+1)-FE(i,j)
            enddo
          enddo            !--> discard FE

          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=0.5*(t(i,j,k,itrc)+t(i,j-1,k,itrc))
     &                                                  *FlxV(i,j,k)
     &          -0.1666666666666666*( curv(i,j-1)*max(FlxV(i,j,k),0.)
     &                               +curv(i,j  )*min(FlxV(i,j,k),0.))
            enddo
          enddo             !--> discard curv,grad, keep FE

