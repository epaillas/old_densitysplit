module procedures
  implicit none
contains

  subroutine linked_list(data, boxsize, ngrid, ll, lirst, rgrid)
    implicit none
    integer*8 :: ndata
    integer*8 :: i, ipx, ipy, ipz
    integer*8, intent(in) :: ngrid
    real*8, intent(in) :: boxsize
    real*8, intent(out) :: rgrid
    real*8, dimension(:,:), intent(in) :: data
    integer*8, dimension(:,:,:), allocatable, intent(out) :: lirst
    integer*8, dimension(:), allocatable, intent(out) :: ll

    ndata = size(data, dim=2)
    allocate(lirst(ngrid, ngrid, ngrid))
    allocate(ll(ndata))
    rgrid = boxsize / real(ngrid)
    
    lirst = 0
    ll = 0
    do i = 1, ndata
      ipx = int(data(1, i) / rgrid + 1.)
      ipy = int(data(2, i) / rgrid + 1.)
      ipz = int(data(3, i) / rgrid + 1.)
      if(ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.&
      ipz.gt.0.and.ipz.le.ngrid) lirst(ipx, ipy, ipz) = i

    end do

    do i = 1, ndata
      ipx = int(data(1, i) / rgrid + 1.)
      ipy = int(data(2, i) / rgrid + 1.)
      ipz = int(data(3, i) / rgrid + 1.)

      if (ipx.gt.0.and.ipx.le.ngrid.and.ipy.gt.0.and.ipy.le.ngrid.and.ipz&
      &.gt.0.and.ipz.le.ngrid) then
        ll(lirst(ipx, ipy, ipz)) = i
        lirst(ipx, ipy, ipz) = i
      endif
    end do

  end subroutine linked_list

    subroutine read_catalogue_type1(input_filename, datarray,&
        weight, ndata)
        implicit none
        ! Read catalogue type 1: text file with 3 columns,
        ! corresponding to x, y, z positions.
        character(len=500), intent(in) :: input_filename
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight
        integer*8, intent(out) :: ndata
        integer*8 :: i


        ! figure out number of rows
        call count_rows_ascii(input_filename, ndata)

        ! read data
        allocate(datarray(3, ndata))
        allocate(weight(ndata))
        open(10, file=input_filename, status='old')
        do i = 1, ndata
            read(10, *) datarray(:, i) 
        end do
        close(10)

        ! equal weighting
        weight = 1.0

    end subroutine read_catalogue_type1

    subroutine read_catalogue_type2(input_filename, datarray,&
        weight, ndata)
        implicit none
        ! Read catalogue type 2: text file with 4 columns,
        ! corresponding to x, y, z positions plus one weight
        ! column.
        character(len=500), intent(in) :: input_filename
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight
        integer*8, intent(out) :: ndata
        integer*8 :: i


        ! figure out number of rows
        call count_rows_ascii(input_filename, ndata)

        ! read data
        allocate(datarray(4, ndata))
        allocate(weight(ndata))
        open(10, file=input_filename, status='old')
        do i = 1, ndata
            read(10, *) datarray(:, i) 
        end do
        close(10)

        weight = datarray(4, :)

    end subroutine read_catalogue_type2

    subroutine read_catalogue_type3(input_filename, datarray,&
        weight, ndata)
        implicit none
        ! Read catalogue type 3: text file with 6 columns,
        ! corresponding to x, y, z positions and velocities.
        character(len=500), intent(in) :: input_filename
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight
        integer*8, intent(out) :: ndata
        integer*8 :: i


        ! figure out number of rows
        call count_rows_ascii(input_filename, ndata)

        ! read data
        allocate(datarray(6, ndata))
        allocate(weight(ndata))
        open(10, file=input_filename, status='old')
        do i = 1, ndata
            read(10, *) datarray(:, i) 
        end do
        close(10)

        ! equal weighting
        weight = 1.0

    end subroutine read_catalogue_type3

    subroutine read_catalogue_type4(input_filename, datarray,&
        weight, ndata)
        implicit none
        ! Read catalogue type 4: text file with 7 columns,
        ! corresponding to x, y, z positions and velocities,
        ! plus one weight column.
        character(len=500), intent(in) :: input_filename
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight
        integer*8, intent(out) :: ndata
        integer*8 :: i


        ! figure out number of rows
        call count_rows_ascii(input_filename, ndata)

        ! read data
        allocate(datarray(7, ndata))
        allocate(weight(ndata))
        open(10, file=input_filename, status='old')
        do i = 1, ndata
            read(10, *) datarray(:, i) 
        end do
        close(10)

        weight = datarray(7, :)

    end subroutine read_catalogue_type4


    subroutine read_catalogue_type5(input_filename, datarray,&
         weight, ndata)
        implicit none
        ! Read catalogue type 5: binary file with 3 columns,
        ! corresponding to x, y, z positions
        integer*8 :: nrows, ncols
        character(len=500), intent(in) :: input_filename
        integer*8, intent(out) :: ndata
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight


        open(20, file=input_filename, status='old', form='unformatted')
        read(20) nrows
        read(20) ncols
        allocate(datarray(ncols, nrows))
        allocate(weight(nrows))
        read(20) datarray
        close(20)
        ndata = nrows

        ! equal weighting
        weight = 1.0

    end subroutine read_catalogue_type5


    subroutine read_catalogue_type6(input_filename, datarray,&
         weight, ndata)
        implicit none
        ! Read catalogue type 6: binary file with 4 columns,
        ! corresponding to x, y, z positions, plus one weight
        ! column
        integer*8 :: nrows, ncols
        character(len=500), intent(in) :: input_filename
        integer*8, intent(out) :: ndata
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight


        open(20, file=input_filename, status='old', form='unformatted')
        read(20) nrows
        read(20) ncols
        allocate(datarray(ncols, nrows))
        allocate(weight(nrows))
        read(20) datarray
        close(20)
        ndata = nrows

        weight = datarray(4, :)

    end subroutine read_catalogue_type6

    subroutine read_catalogue_type7(input_filename, datarray,&
         weight, ndata)
        implicit none
        ! Read catalogue type 7: binary file with 6 columns,
        ! corresponding to x, y, z positions and velocities.
        integer*8 :: nrows, ncols
        character(len=500), intent(in) :: input_filename
        integer*8, intent(out) :: ndata
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight


        open(20, file=input_filename, status='old', form='unformatted')
        read(20) nrows
        read(20) ncols
        allocate(datarray(ncols, nrows))
        allocate(weight(nrows))
        read(20) datarray
        close(20)
        ndata = nrows

        ! equal weighting
        weight = 1.0

    end subroutine read_catalogue_type7

    subroutine read_catalogue_type8(input_filename, datarray,&
         weight, ndata)
        implicit none
        ! Read catalogue type 8: binary file with 7 columns,
        ! corresponding to x, y, z positions and velocities,
        ! plus one weight column.
        integer*8 :: nrows, ncols
        character(len=500), intent(in) :: input_filename
        integer*8, intent(out) :: ndata
        real*8, allocatable, dimension(:,:), intent(out) :: datarray
        real*8, allocatable, dimension(:), intent(out) :: weight


        open(20, file=input_filename, status='old', form='unformatted')
        read(20) nrows
        read(20) ncols
        allocate(datarray(ncols, nrows))
        allocate(weight(nrows))
        read(20) datarray
        close(20)
        ndata = nrows

        weight = datarray(:, 7)

    end subroutine read_catalogue_type8


  subroutine read_unformatted(input_filename, data, weight, ndata, has_velocity)
    implicit none
    integer*8 :: nrows, ncols
    character(len=500), intent(in) :: input_filename
    integer*8, intent(out) :: ndata
    real*8, allocatable, dimension(:,:), intent(out) :: data
    real*8, allocatable, dimension(:), intent(out) :: weight
    logical, intent(out) :: has_velocity

    has_velocity = .false.

    open(20, file=input_filename, status='old', form='unformatted')
    read(20) nrows
    read(20) ncols
    allocate(data(ncols, nrows))
    allocate(weight(nrows))
    read(20) data
    close(20)
    ndata = nrows
    weight = 1.0 ! default weights
    if (ncols .ge. 4) then
      if (ncols .eq. 4) weight = data(4, :)
      if (ncols .eq. 6) has_velocity = .true.
      if (ncols .eq. 7) then
        weight = data(7, :)
        has_velocity = .true.
      end if
    end if
      
  end subroutine read_unformatted

  subroutine binning(rmin, rmax, nbin, bin, bin_edges, rwidth)
    implicit none

    integer*8 :: i
    integer*8, intent(in) :: nbin
    real*8, intent(in) :: rmin, rmax
    real*8, intent(out) :: rwidth
    real*8, allocatable, dimension(:), intent(out) :: bin, bin_edges

    allocate(bin(nbin))
    allocate(bin_edges(nbin + 1))

    rwidth = (rmax - rmin) / nbin
    do i = 1, nbin + 1
      bin_edges(i) = rmin + (i - 1) * rwidth
    end do
    do i = 1, nbin
      bin(i) = bin_edges(i + 1) - rwidth / 2.
    end do

  end subroutine binning

end module procedures
