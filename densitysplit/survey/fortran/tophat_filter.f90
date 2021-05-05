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
  real*8 :: disx, disy, disz, dis2
  real*8 :: dim1_min, dim1_max, dim1_min2, dim1_max2
  
  integer*8 :: ng1, ng2, nr1, nr2
  integer*8 :: i, ii, ix, iy, iz
  integer*8 :: nrows, ncols
  integer*8 :: ipx, ipy, ipz, ndif
  integer*8 :: ngrid
  integer*4 :: nthreads
  integer*8 :: end, beginning, rate
  integer*8, dimension(:, :, :), allocatable :: lirst_data2, lirst_random2

  integer*8, dimension(:), allocatable :: ll_data2, ll_random2
  
  real*8, allocatable, dimension(:,:)  :: data1, data2, random1, random2
  real*8, dimension(:), allocatable :: D1D2, D1R2, R1D2, R1R2, xi_r
  real*8, dimension(:), allocatable :: weight_data1, weight_data2
  real*8, dimension(:), allocatable :: weight_random1, weight_random2

  logical :: debug = .true.
  
  character(20), external :: str
  character(len=500) :: data_filename2, data_filename1, output_filename
  character(len=500) :: random_filename1, random_filename2
  character(len=10) :: dmax_char, dmin_char, gridmin_char, gridmax_char
  character(len=10) :: ngrid_char, rfilter_char, nthreads_char
  character(len=2) :: estimator

  if (iargc() .ne. 13) then
    write(*,*) 'Some arguments are missing.'
    write(*,*) '1) data_filename1'
    write(*,*) '2) data_filename2'
    write(*,*) '3) random_filename1'
    write(*,*) '4) random_filename2'
    write(*,*) '5) output_filename'
    write(*,*) '6) dim1_min'
    write(*,*) '7) dim1_max'
    write(*,*) '8) rfilter'
    write(*,*) '9) ngrid'
    write(*,*) '10) gridmin'
    write(*,*) '11) gridmax'
    write(*,*) '12) estimator'
    write(*,*) '13) nthreads'
    write(*,*) ''
    stop
  end if

  call system_clock(beginning, rate)

  ! read arguments from command line
  call getarg(1, data_filename1)
  call getarg(2, data_filename2)
  call getarg(3, random_filename1)
  call getarg(4, random_filename2)
  call getarg(5, output_filename)
  call getarg(6, dmin_char)
  call getarg(7, dmax_char)
  call getarg(8, rfilter_char)
  call getarg(9, ngrid_char)
  call getarg(10, gridmin_char)
  call getarg(11, gridmax_char)
  call getarg(12, estimator)
  call getarg(13, nthreads_char)
  
  ! convert string arguments to corresponding data types
  read(dmin_char, *) dim1_min
  read(dmax_char, *) dim1_max
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
    write(*,*) 'data_filename1: ', trim(data_filename1)
    write(*,*) 'data_filename2: ', trim(data_filename2)
    write(*,*) 'random_filename1: ', trim(random_filename1)
    write(*,*) 'random_filename2: ', trim(random_filename2)
    write(*,*) 'output_filename: ', trim(output_filename)
    write(*,*) 'dmin: ', trim(dmin_char), ' Mpc'
    write(*,*) 'dmax: ', trim(dmax_char), ' Mpc'
    write(*,*) 'rfilter: ', trim(rfilter_char), 'Mpc'
    write(*,*) 'ngrid: ', trim(ngrid_char)
    write(*,*) 'gridmin: ', trim(gridmin_char), 'Mpc'
    write(*,*) 'gridmax: ', trim(gridmax_char), 'Mpc'
    write(*,*) 'estimator: ', trim(estimator)
    write(*,*) 'nthreads: ', trim(nthreads_char)
    write(*,*) ''
  end if

  
  ! read data catalogue #1
  open(11, file=data_filename1, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(data1(ncols, nrows))
  allocate(weight_data1(nrows))
  read(11) data1
  close(11)
  ng1 = nrows
  if (ncols .eq. 4) then
    weight_data1 = data1(4, :)
    if (debug) write(*,*) 'Data file 1 has weight information.'
  else
    weight_data1 = 1.0
  end if
  if (debug) then
    write(*,*) 'ndata1 dim: ', size(data1, dim=1), size(data1, dim=2)
    write(*,*) 'data1(min), data1(max) = ', minval(data1(:,:)), maxval(data1(:,:))
    write(*,*) 'weight_data1(min), weight_data1(max) = ', minval(weight_data1), maxval(weight_data1)
  end if

  ! read data catalogue #2
  open(10, file=data_filename2, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(data2(ncols, nrows))
  allocate(weight_data2(nrows))
  read(10) data2
  close(10)
  ng2 = nrows
  if (ncols .eq. 4) then
    weight_data2 = data2(4, :)
    if (debug) write(*,*) 'Data file 2 has weight information.'
  else
    weight_data2 = 1.0
  end if
  if (debug) then
    write(*,*) 'ndata2 dim: ', size(data2, dim=1), size(data2, dim=2)
    write(*,*) 'data2(min), data2(max) = ', minval(data2(:,:)), maxval(data2(:,:))
    write(*,*) 'weight_data2(min), weight_data2(max) = ', minval(weight_data2), maxval(weight_data2)
  end if

  if (estimator .eq. 'LS') then
    ! read random catalogue #1
    open(11, file=random_filename1, status='old', form='unformatted')
    read(11) nrows
    read(11) ncols
    allocate(random1(ncols, nrows))
    allocate(weight_random1(nrows))
    read(11) random1
    close(11)
    nr1 = nrows
    if (ncols .eq. 4) then
      weight_random1 = random1(4, :)
      if (debug) write(*,*) 'Random file 1 has weight information.'
    else
      weight_random1 = 1.0
    end if
    if (debug) then 
      write(*,*) 'nrandom1 dim: ', size(random1, dim=1), size(random1, dim=2)
      write(*,*) 'random1(min), random1(max) = ', minval(random1(:,:)), maxval(random1(:,:))
      write(*,*) 'weight_data1(min), weight_data1(max) = ', minval(weight_data1), maxval(weight_data1)
    end if
  end if

  ! read random2 catalogue
  open(11, file=random_filename2, status='old', form='unformatted')
  read(11) nrows
  read(11) ncols
  allocate(random2(ncols, nrows))
  allocate(weight_random2(nrows))
  read(11) random2
  close(11)
  nr2 = nrows
  if (ncols .eq. 4) then
    weight_random2 = random2(4, :)
    if (debug) write(*,*) 'Random file 2 has weight information.'
  else
    weight_random2 = 1.0
  end if
  if (debug) then 
    write(*,*) 'nrandom2 dim: ', size(random2, dim=1), size(random2, dim=2)
    write(*,*) 'random2(min), random2(max) = ', minval(random2(:,:)), maxval(random2(:,:))
    write(*,*) 'weight_random2(min), weight_random2(max) = ', minval(weight_random2), maxval(weight_random2)
  end if
  
  ! construct linked lists for data2 and random2
  allocate(ll_data2(ng2))
  allocate(ll_random2(nr2))
  allocate(lirst_data2(ngrid, ngrid, ngrid))
  allocate(lirst_random2(ngrid, ngrid, ngrid))
  call linked_list(data2, ngrid, gridmin, gridmax, ll_data2, lirst_data2, rgrid)
  call linked_list(random2, ngrid, gridmin, gridmax, ll_random2, lirst_random2, rgrid)

  ! calculate number counts around each centre
  allocate(D1D2(ng1))
  allocate(D1R2(ng1))
  allocate(xi_r(ng1))
  D1D2 = 0
  D1R2 = 0
  if (estimator .eq. 'LS') then
    allocate(R1R2(nr1))
    allocate(R1D2(nr1))
    R1R2 = 0
    R1D2 = 0
  end if
  
  ndif = int(dim1_max / rgrid + 1.)
  dim1_min2 = dim1_min ** 2
  dim1_max2 = dim1_max ** 2
  

  call OMP_SET_NUM_THREADS(nthreads)
  if (debug) then
    write(*,*) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  end if

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ipx, ipy, ipz, &
  !$OMP& ix, iy, iz, ii, disx, disy, disz, dis2)
  do i = 1, ng1
    ipx = int((data1(1, i) - gridmin) / rgrid + 1.)
    ipy = int((data1(2, i) - gridmin) / rgrid + 1.)
    ipz = int((data1(3, i) - gridmin) / rgrid + 1.)

    ! loop over cells around each centre
    do ix = ipx - ndif, ipx + ndif, 1
      do iy = ipy - ndif, ipy + ndif, 1
        do iz = ipz - ndif, ipz + ndif, 1
          if ((ix-ipx)**2 + (iy-ipy)**2 + (iz-ipz)**2 .gt. (ndif + 1)**2) cycle

          ! loop over data2 in each cell
          ii = lirst_data2(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_data2(ii)
              disx = data2(1, ii) - data1(1, i)
              disy = data2(2, ii) - data1(2, i)
              disz = data2(3, ii) - data1(3, i)

              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                D1D2(i) = D1D2(i) + weight_data1(i) * weight_data2(ii)
              end if
  
              if(ii.eq.lirst_data2(ix, iy, iz)) exit
  
            end do
          end if

          ! loop over random2 in each cell
          ii = lirst_random2(ix, iy, iz)
          if (ii .ne. 0) then
            do
              ii = ll_random2(ii)
              disx = random2(1, ii) - data1(1, i)
              disy = random2(2, ii) - data1(2, i)
              disz = random2(3, ii) - data1(3, i)
 
              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                D1R2(i) = D1R2(i) + weight_data1(i) * weight_random2(ii)
              end if
  
              if(ii .eq. lirst_random2(ix, iy, iz)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  if (estimator .eq. 'LS') then 
    ! Loop over randoms # 1
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, ipz, &
    !$OMP ix, iy, iz, disx, disy, disz, dis2)
    do i = 1, nr1
      ipx = int((random1(1, i) - gridmin) / rgrid + 1.)
      ipy = int((random1(2, i) - gridmin) / rgrid + 1.)
      ipz = int((random1(3, i) - gridmin) / rgrid + 1.)

      do ix = ipx - ndif, ipx + ndif, 1
        do iy = ipy - ndif, ipy + ndif, 1
          do iz = ipz - ndif, ipz + ndif, 1 
            if ((ix - ipx)**2 + (iy - ipy)**2 + (iz - ipz)**2 .gt. (ndif + 1)**2) cycle

            ii = lirst_random2(ix, iy, iz)
            if (ii .ne. 0) then
              do
                ii = ll_random2(ii)
                disx = random2(1, ii) - random1(1, i)
                disy = random2(2, ii) - random1(2, i)
                disz = random2(3, ii) - random1(3, i)

                dis2 = disx * disx + disy * disy + disz * disz

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  R1R2(i) = R1R2(i) + weight_random1(i) * weight_random2(ii)
                end if

                  if(ii .eq. lirst_random2(ix, iy, iz)) exit

              end do
            end if

            ii = lirst_data2(ix, iy, iz)
            if (ii .ne. 0) then
              do
                ii = ll_data2(ii)
                disx = data2(1, ii) - random1(1, i)
                disy = data2(2, ii) - random1(2, i)
                disz = data2(3, ii) - random1(3, i)

                dis2 = disx * disx + disy * disy + disz * disz

                if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                  R1D2(i) = R1D2(i) + weight_random1(i) * weight_data2(ii)
                end if

                  if(ii .eq. lirst_data2(ix, iy, iz)) exit

              end do
            end if
          end do
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end if

  ! Normalize pair counts
  D1D2 = D1D2 * 1. / (SUM(weight_data1) * SUM(weight_data2))
  D1R2 = D1R2 * 1. / (SUM(weight_data1) * SUM(weight_random2))
  if (estimator .eq. 'LS') then
    R1R2 = R1R2 * 1. / (SUM(weight_random1) * SUM(weight_random2))
    R1D2 = R1D2 * 1. / (SUM(weight_random1) * SUM(weight_data2))
  end if

  ! Calculate density contrast
  if (estimator .eq. 'DP') then
    xi_r = (D1D2 / D1R2) - 1
  else if (estimator .eq. 'LS') then
    xi_r = (D1D2 - D1R2 - R1D2 + R1R2) / R1R2
  else
    write(*,*) 'Estimator for the correlation function was not recognized.'
    stop
  end if
  
  ! write output  
  open(12, file=output_filename, status='replace', form='unformatted')
  write(12) ng1
  write(12) xi_r

  call system_clock(end)
  if (debug) then
    print *, "Elapsed time: ", real(end - beginning) / real(rate)
  end if

  end program tophat_filter
  
  
