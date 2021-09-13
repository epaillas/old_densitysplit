program tophat_filter
  use OMP_LIB
  use procedures
  implicit none
  
  real*8 :: bin_volume, mean_density
  real*8 :: boxsize, box2
  real*8 :: disx, disy, disz, dis, dis2
  real*8 :: dim1_min2, dim1_max2
  real*8 :: dim1_max, dim1_min, rfilter, rgrid
  real*8 :: pi = 4.*atan(1.)
  
  integer*8 :: ndata2, ndata1
  integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
  integer*8 :: ipx, ipy, ipz, ndif, ngrid
  integer*4 :: nthreads, use_weights
  
  integer*8, dimension(:, :, :), allocatable :: lirst
  integer*8, dimension(:), allocatable :: ll
  
  real*8, allocatable, dimension(:,:)  :: data2, data1
  real*8, dimension(:), allocatable :: D1D2, delta, weight1, weight2

  logical :: has_velocity1, has_velocity2

  character(20), external :: str
  character(len=500) :: data_filename2, data_filename1, output_filename
  character(len=10) :: dim1_max_char, dim1_min_char
  character(len=10) :: boxsize_char, rfilter_char, nthreads_char
  character(len=10) :: ngrid_char, use_weights_char, tracers_fileformat


  if (iargc() .lt. 11) then
      write(*,*) 'Some arguments are missing.'
      write(*,*) '1) data_filename1'
      write(*,*) '2) data_filename2'
      write(*,*) '3) output_filename'
      write(*,*) '4) boxsize'
      write(*,*) '5) dim1_min'
      write(*,*) '6) dim1_max'
      write(*,*) '7) rfilter'
      write(*,*) '8) ngrid'
      write(*,*) '9) nthreads'
      write(*,*) '10) use_weights'
      write(*,*) '11) tracers_fileformat'
      write(*,*) ''
      stop
    end if
    
    call getarg(1, data_filename1)
    call getarg(2, data_filename2)
    call getarg(3, output_filename)
    call getarg(4, boxsize_char)
    call getarg(5, dim1_min_char)
    call getarg(6, dim1_max_char)
    call getarg(7, rfilter_char)
    call getarg(8, ngrid_char)
    call getarg(9, nthreads_char)
    call getarg(10, use_weights_char)
    call getarg(11, tracers_fileformat)
  
    read(boxsize_char, *) boxsize
    read(dim1_min_char, *) dim1_min
    read(dim1_max_char, *) dim1_max
    read(rfilter_char, *) rfilter
    read(ngrid_char, *) ngrid
    read(nthreads_char, *) nthreads
    read(use_weights_char, *) use_weights

    ! read tracers file
    if (trim(tracers_fileformat) == 'ascii') then
        if (use_weights == 1) then
            call read_catalogue_type2(data_filename1, data1, weight1, ndata1, has_velocity1)
            call read_catalogue_type2(data_filename2, data2, weight2, ndata2, has_velocity2)
        else
            call read_catalogue_type1(data_filename1, data1, weight1, ndata1, has_velocity1)
            call read_catalogue_type1(data_filename2, data2, weight2, ndata2, has_velocity2)
        end if
    else
        if (use_weights == 1) then
            call read_catalogue_type6(data_filename1, data1, weight1, ndata1, has_velocity1)
            call read_catalogue_type6(data_filename2, data2, weight2, ndata2, has_velocity2)
        else
            call read_catalogue_type5(data_filename1, data1, weight1, ndata1, has_velocity1)
            call read_catalogue_type5(data_filename2, data2, weight2, ndata2, has_velocity2)
        end if
    end if

  call linked_list(data2, boxsize, ngrid, ll, lirst, rgrid)


  allocate(D1D2(ndata1))
  allocate(delta(ndata1))
  
  mean_density = ndata2 / (boxsize ** 3)
  D1D2 = 0
  delta = 0
  ndif = int(dim1_max / rgrid + 1.)
  dim1_min2 = dim1_min * dim1_min
  dim1_max2 = dim1_max * dim1_max
  box2 = boxsize / 2

  call OMP_SET_NUM_THREADS(nthreads)
  write(*,*) 'Maximum number of threads: ', OMP_GET_MAX_THREADS()
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, ii, ipx, ipy, &
  !$OMP ipz, ix, iy, iz, ix2, iy2, iz2, disx, disy, disz, dis2)
  do i = 1, ndata1
    ipx = int(data1(1, i) / rgrid + 1.)
    ipy = int(data1(2, i) / rgrid + 1.)
    ipz = int(data1(3, i) / rgrid + 1.)
  
    do ix = ipx - ndif, ipx + ndif
      do iy = ipy - ndif, ipy + ndif
        do iz = ipz - ndif, ipz + ndif
  
          ix2 = ix
          iy2 = iy
          iz2 = iz
  
          if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
          if (ix2 .lt. 1) ix2 = ix2 + ngrid
          if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
          if (iy2 .lt. 1) iy2 = iy2 + ngrid
          if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
          if (iz2 .lt. 1) iz2 = iz2 + ngrid
  
          ii = lirst(ix2,iy2,iz2)
          if(ii .ne. 0) then
            do
              ii = ll(ii)
              disx = data2(1, ii) - data1(1, i)
              disy = data2(2, ii) - data1(2, i)
              disz = data2(3, ii) - data1(3, i)

              if (disx .lt. -box2) disx = disx + boxsize
              if (disx .gt. box2) disx = disx - boxsize
              if (disy .lt. -box2) disy = disy + boxsize
              if (disy .gt. box2) disy = disy - boxsize
              if (disz .lt. -box2) disz = disz + boxsize
              if (disz .gt. box2) disz = disz - boxsize
  
              dis2 = disx * disx + disy * disy + disz * disz

              if (dis2 .gt. dim1_min2 .and. dis2 .lt. dim1_max2) then
                D1D2(i) = D1D2(i) + weight1(i) * weight2(ii)
              end if

              if(ii .eq. lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do
  !$OMP END PARALLEL DO

  bin_volume = 4./3 * pi * (dim1_max ** 3 - dim1_min ** 3)
  delta = D1D2 / (bin_volume * mean_density) - 1
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_filename, status='replace', form='unformatted')
  write(12) ndata1
  write(12) delta

end program tophat_filter
    
