module procedures
  implicit none
contains

  subroutine linked_list(pos, ngrid, gridmin, gridmax, ll, lirst, rgrid)
    implicit none
    integer*8 :: i, ng, ipx, ipy, ipz
    integer*8, intent(in) :: ngrid
    real*8, intent(out) :: rgrid
    real*8, intent(in) :: gridmin, gridmax
    real*8, dimension(:,:), intent(in) :: pos
    integer*8, dimension(:,:,:), intent(out) :: lirst
    integer*8, dimension(:), intent(out) :: ll

    rgrid = (gridmax - gridmin) / ngrid

    ng = size(pos, dim=2)
    lirst = 0
    ll = 0
    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.&
      ipz.gt.0.and.ipz.le.ngrid) lirst(ipx, ipy, ipz) = i

    end do

    do i = 1, ng
      ipx = int((pos(1, i) - gridmin) / rgrid + 1.)
      ipy = int((pos(2, i) - gridmin) / rgrid + 1.)
      ipz = int((pos(3, i) - gridmin) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.ipz&
      &.gt.0.and.ipz.le.ngrid) then
        ll(lirst(ipx, ipy, ipz)) = i
        lirst(ipx, ipy, ipz) = i
      endif
    end do

  end subroutine linked_list

end module procedures


program tophat_filter
  use procedures
  use OMP_LIB
  implicit none
  
  real*8 :: rgrid, rfilter, gridmax, gridmin
  real*8 :: disx, disy, disz, dis, dmax, dmin
  
  integer*8 :: ng, nc, nr
  integer*8 :: i, ii, ix, iy, iz
  integer*8 :: nrows, ncols
  integer*8 :: ipx, ipy, ipz, ndif
  integer*8 :: ngrid
  integer*4 :: nthreads, threadid
  integer*8, dimension(:, :, :), allocatable :: lirst_tracers, lirst_randoms

  integer*8, dimension(:), allocatable :: ll_tracers, ll_randoms
  
  real*8, dimension(3) :: r
  real*8, allocatable, dimension(:,:)  :: tracers, randoms, centres
  real*8, dimension(:), allocatable :: DD, RR, delta
  real*8, dimension(:), allocatable :: weights_tracers, weights_randoms

  logical :: debug = .false.
  
  character(20), external :: str
  character(len=500) :: input_tracers, input_centres, input_randoms, output_filter
  character(len=10) :: dmax_char, dmin_char, gridmin_char, gridmax_char
  character(len=10) :: ngrid_char, rfilter_char, nthreads_char
  
  if (iargc() .ne. 11) then
    write(*,*) 'Some arguments are missing.'
    write(*,*) '1) input_tracers'
    write(*,*) '2) input_centres'
    write(*,*) '3) input_randoms'
    write(*,*) '4) output_filter'
    write(*,*) '5) dmin'
    write(*,*) '6) dmax'
    write(*,*) '7) rfilter'
    write(*,*) '8) ngrid'
    write(*,*) '9) gridmin'
    write(*,*) '10) gridmax'
    write(*,*) '11) nthreads'
    write(*,*) ''
    stop
  end if

  ! read arguments from command line
  call getarg(1, input_tracers)
  call getarg(2, input_centres)
  call getarg(3, input_randoms)
  call getarg(4, output_filter)
  call getarg(5, dmin_char)
  call getarg(6, dmax_char)
  call getarg(7, rfilter_char)
  call getarg(8, ngrid_char)
  call getarg(9, gridmin_char)
  call getarg(10, gridmax_char)
  call getarg(11, nthreads_char)
  
  ! convert string arguments to corresponding data types
  read(dmin_char, *) dmin
  read(dmax_char, *) dmax
  read(rfilter_char, *) rfilter
  read(ngrid_char, *) ngrid
  read(gridmin_char, *) gridmin
  read(gridmax_char, *) gridmax
  read(nthreads_char, *) nthreads


  if (debug) then
    write(*,*) '-----------------------'
    write(*,*) 'Running tophat_filter.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'input_tracers: ', trim(input_tracers)
    write(*, *) 'input_centres: ', trim(input_centres)
    write(*, *) 'output_filter: ', trim(output_filter)
    write(*, *) 'dmin: ', trim(dmin_char), ' Mpc'
    write(*, *) 'dmax: ', trim(dmax_char), ' Mpc'
    write(*, *) 'rfilter: ', trim(rfilter_char), 'Mpc'
    write(*, *) 'ngrid: ', trim(ngrid_char)
    write(*, *) 'gridmin: ', trim(gridmin_char), 'Mpc'
    write(*, *) 'gridmax: ', trim(gridmax_char), 'Mpc'
    write(*, *) 'nthreads: ', trim(nthreads_char)
    write(*,*) ''
  end if

  ! read tracers catalogue
  open(10, file=input_tracers, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(tracers(ncols, nrows))
  allocate(weights_tracers(nrows))
  read(10) tracers
  close(10)
  ng = nrows
  if (ncols .eq. 4) then
    weights_tracers = tracers(4, :)
    if (debug) write(*,*) 'Tracer file has weight information.'
  else
    weights_tracers = 1.0
  end if
  if (debug) then
    write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
    write(*,*) 'tracers(min), tracers(max) = ', minval(tracers(:,:)), maxval(tracers(:,:))
  end if

  ! read centres catalogue
  open(11, file=input_centres, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(centres(ncols, nrows))
  read(11) centres
  close(11)
  nc = nrows
  if (debug) then
    write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)
    write(*,*) 'centres(min), tracers(max) = ', minval(centres(:,:)), maxval(centres(:,:))
  end if

  ! read random catalogue
  open(11, file=input_randoms, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(randoms(ncols, nrows))
  allocate(weights_randoms(nrows))
  read(11) randoms
  close(11)
  nr = nrows
  if (ncols .eq. 4) then
    weights_randoms = tracers(4, :)
    if (debug) write(*,*) 'Tracer file has weight information.'
  else
    weights_randoms = 1.0
  end if
  if (debug) then 
    write(*,*) 'nrandoms dim: ', size(randoms, dim=1), size(randoms, dim=2)
    write(*,*) 'randoms(min), randoms(max) = ', minval(randoms(:,:)), maxval(randoms(:,:))
  end if
  
  ! construct linked lists for tracers and randoms
  allocate(ll_tracers(ng))
  allocate(ll_randoms(nr))
  allocate(lirst_tracers(ngrid, ngrid, ngrid))
  allocate(lirst_randoms(ngrid, ngrid, ngrid))
  call linked_list(tracers, ngrid, gridmin, gridmax, ll_tracers, lirst_tracers, rgrid)
  call linked_list(randoms, ngrid, gridmin, gridmax, ll_randoms, lirst_randoms, rgrid)

  ! calculate number counts around each centre
  allocate(DD(nc))
  allocate(RR(nc))
  allocate(delta(nc))
  DD = 0
  RR = 0
  delta = 0
  ndif = int(dmax / rgrid + 1.)

  call OMP_SET_NUM_THREADS(nthreads)
  if (debug) then
    write(*, *) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  end if

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ipx, ipy, ipz, &
  !$OMP& ix, iy, iz, ii, disx, disy, disz, r, dis, threadid)
  do i = 1, nc
    ipx = int((centres(1, i) - gridmin) / rgrid + 1.)
    ipy = int((centres(2, i) - gridmin) / rgrid + 1.)
    ipz = int((centres(3, i) - gridmin) / rgrid + 1.)

    threadid = OMP_get_thread_num()
    write(*,*) i, threadid

    ! loop over cells around each centre
    do ix = ipx - ndif, ipx + ndif, 1
      do iy = ipy - ndif, ipy + ndif, 1
        do iz = ipz - ndif, ipz + ndif, 1

          ! loop over tracers in each cell
          ii = lirst_tracers(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_tracers(ii)
              disx = tracers(1, ii) - centres(1, i)
              disy = tracers(2, ii) - centres(2, i)
              disz = tracers(3, ii) - centres(3, i)
 
              r = (/ disx, disy, disz /)
              dis = norm2(r)

              if (dis .gt. dmin .and. dis .lt. dmax) then
                DD(i) = DD(i) + weights_tracers(ii)
              end if
  
              if(ii.eq.lirst_tracers(ix, iy, iz)) exit
  
            end do
          end if

          ! loop over randoms in each cell
          ii = lirst_randoms(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_randoms(ii)
              disx = randoms(1, ii) - centres(1, i)
              disy = randoms(2, ii) - centres(2, i)
              disz = randoms(3, ii) - centres(3, i)
 
              r = (/ disx, disy, disz /)
              dis = norm2(r)

              if (dis .gt. dmin .and. dis .lt. dmax) then
                RR(i) = RR(i) + weights_randoms(ii)
              end if
  
              if(ii .eq. lirst_randoms(ix, iy, iz)) exit
  
            end do
          end if

        end do
      end do
    end do
    ! normalize random counts
    RR = RR / (nr * 1./ng)

    ! calculate density
    delta(i) = DD(i) / RR(i) - 1

  end do
  !$OMP END PARALLEL DO
  
  ! write output  
  open(12, file=output_filter, status='replace', form='unformatted')
  write(12) nc
  write(12) delta

  end program tophat_filter
  