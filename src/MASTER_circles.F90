program circles

  use alm_tools
  use healpix_types
  use fitstools
  use pix_tools
  use ran_tools
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
#include"read_archive_map.fh"

#include<fftw3.f>

#ifdef MPI
#include<mpif.h>
#endif

#ifdef DFFT
#define FP dp
#define FFT_COMPLEX double complex
#else
#define FP sp
#define FFT_COMPLEX complex
#endif

  type tcircle 
     real(dp)    :: phi, theta, opening_angle
     integer(dp) :: nphi_ring, pixel
  end type tcircle

  type tangle
     real(dp) :: theta, phi
  end type tangle

  type tpair
     type(tcircle):: circle1, circle2
     real :: smax
  end type tpair

  type tangle_range
     real    :: min, max, stepsize
     integer :: num_steps
  end type tangle_range

  integer, parameter        :: MAX_THREAD_NUM=32

  integer                   :: nside_map !512  !Healpix resolution of upgraded map
  integer                   :: nside_search
  integer                   :: m_min = 2
  integer                   :: iopening_angle, i, nphi_ring, npix,j, npix_search
  integer(dp)               :: analysed_circles
  logical                   :: use_pixel_average = .true.
  real(dp)                  :: mask_threshold = 0.5_dp
  real,         allocatable :: map(:)                     ! size: npix = 12*nside_map**2
  real,         allocatable :: obs_mask(:)
  real, allocatable         :: sw_filtered_map(:)
  real(dp), allocatable     :: sw_filter(:)

  ! these tables are used on the fine grid when computing the average in get_temperature_along_circle()
  integer,      allocatable :: nest2ring_table(:)         ! size: npix=12*nside_map**2
  integer,      allocatable :: ring2nest_table(:)         ! size: npix=12*nside_map**2


  ! this table is used to determine the angle to a given position in fill_circle_transforms()
  type(tangle), allocatable :: pix2ang_table(:)           ! size: npix = 12*nside_map**2

  real,         allocatable :: theta_rings(:)             ! size: nside_search/2:2*nside_search
  real(FP),     allocatable :: smax_angle(:)              ! size: nside_search/2:2*nside_search
  integer(FP)               :: smax_pos
  real(FP),     allocatable :: smax_result(:)             ! size: nside_search/2:2*nside_search
  integer(FP),  allocatable :: smax_circle1_pixel(:)      ! size: nside_search/2:2*nside_search
  integer(FP),  allocatable :: smax_circle2_pixel(:)      ! size: nside_search/2:2*nside_search
  integer(FP),  allocatable :: smax_circle1_pixel_result(:)   ! size: nside_search/2:2*nside_search
  integer(FP),  allocatable :: smax_circle2_pixel_result(:)   ! size: nside_search/2:2*nside_search
  integer(FP),  allocatable :: smax_beta_max_bin(:)
  integer(FP),  allocatable :: smax_beta_max_bin_result(:)

  integer,      allocatable :: opposite_list(:)           ! size: npix_search = 12*nside_search**2
  real,         allocatable :: circle_center_vectors(:,:) ! size: (3,npix_search)

  FFT_COMPLEX,  allocatable :: temperature_out(:)         ! size: nphi_ring/2+1
  FFT_COMPLEX,  allocatable :: sm_in(:)                   ! size: nphi_ring/2+1

  FFT_COMPLEX, allocatable  :: &
       & master_plan_in(:)

  FFT_COMPLEX, allocatable  :: &
       & master_plan_out(:)

  real(FP),     allocatable :: temperature_in(:)          ! size: nphi_ring
  real(FP),     allocatable :: sm_out(:)                  ! size: nphi_ring

  FFT_COMPLEX,  allocatable :: circle_transforms(:,:)     ! size: (nphi_ring/2+1, npix_search)
  FFT_COMPLEX,  allocatable :: masked_circle_transforms(:,:)
  FFT_COMPLEX,  allocatable :: mask_transforms(:,:)
  real(dp), allocatable :: master_auto_f_reg(:)
  real(dp), allocatable :: master_auto_denominator(:)
  logical, allocatable :: master_auto_valid(:)
  real(FP), allocatable :: master_auto_modes(:,:)
  
  integer(8) :: psd_fail_modes
  integer(8) :: expected_psd_modes
  real(dp) :: psd_fail_fraction
  integer, parameter :: n_rho_bins = 2000
  
  real(dp), parameter :: rho_min = -1.0_dp
  real(dp), parameter :: rho_max =  1.0_dp
  real(dp), parameter :: rho_bin_width = &
       (rho_max - rho_min) / real(n_rho_bins, dp)
  
  integer(8) :: rho_hist_thread(n_rho_bins, MAX_THREAD_NUM)
  integer(8) :: rho_hist(n_rho_bins)
  integer(8) :: rho_hist_result(n_rho_bins)

  integer(8)                :: plan_forward
  integer(8)                :: plan_backward
  integer(8)                :: plan_master_forward
  integer(8)                :: plan_master_backward
  real(dp)                  :: smax
  real(dp)                  :: selected_f_reg_ij
  real(dp)                  :: selected_f_reg_ii
  real(dp)                  :: selected_f_reg_jj
  real(dp)                  :: selected_f_lost
  real(dp)                  :: position_smax
  real(dp)                  :: opposite_vector(3)
  integer                   :: opposite_number
  integer                   :: verbosity = 0
  integer                   :: outfile = 101
  integer                   :: full_outfile = 110
  integer                   :: position_smax_file = 202
  integer                   :: rho_hist_file = 203
  logical                   :: orientable=.false.
  logical                   :: write_every_nth_smax=.false.
  logical                   :: write_position_smax=.false.
  integer                   :: every_nth_smax = 10000, every_nth_smax_file = 201
  logical                   :: scramble_sky = .false.
  logical                   :: use_sw_filter = .false.
  logical                   :: use_gcv_regularisation = .true.
  logical                   :: use_flost_regularisation = .false.
  real(dp)                  :: rho_tolerance = 1.0_dp
  logical                   :: write_rho_hist = .false.
  logical                   :: write_f_lost = .false.
  character(LEN=1000)       :: sw_filter_filename
  integer                   :: sw_filter_lmax = 300

  real(dp), parameter :: &
       & flost_A = 0.08912509381337459_dp

  integer, parameter :: n_flost_support = 5
  real(dp), parameter :: &
       & flost_support(n_flost_support) = (/ &
       & 0.05_dp, 0.15_dp, 0.25_dp, 0.35_dp, 0.45_dp /)

  type(tangle_range)        :: opening_angle_range
  type(tcircle)             :: circle1, circle2, bestcircle1, bestcircle2
  type(tpair)               :: bestpair(MAX_THREAD_NUM)
  type(tcircle)             :: position_circle2

  character(LEN=1000)       :: inifilename, outfilename, mapfilename
  character(LEN=1000)       :: maskfilename
  character(LEN=1000)       :: position_smax_filename
  character(LEN=1000)       :: rho_hist_filename
  real(dp)                  :: smax_thread(MAX_THREAD_NUM)
  real(dp)                  :: min_separation_angle
  logical                   :: random_s = .false.
  integer                   :: thread_id, beta_max_bin
  logical                   :: final_valid_pair
  integer                   :: position_beta_max_bin
  integer(8)                :: position_psd_fail_modes
  integer                   :: MPIid=0
  integer                   :: MPIsize=1

#ifdef MPI
  integer                   :: ierr
  character(LEN=1000)       :: mpiid_string
#endif
  real(dp) :: starttime, stoptime

#ifdef MPI
  call MPI_init(ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, MPIid, ierr)
  call MPI_comm_size(MPI_COMM_WORLD, MPIsize, ierr)
#endif

  if (iargc() .lt.1)then
     write(*,*) "Usage: circles <inifile>"
     stop
  end if

  call getarg(1,inifilename)
  call initialize(inifilename, nside_map, nphi_ring)
  circle1%nphi_ring=nphi_ring
  circle2%nphi_ring=nphi_ring

#ifdef MPI
  write(mpiid_string,*) mpiid
  outfilename=trim(adjustl(inifilename(1:len_trim(inifilename)-4)))//"_mpi_"//trim(adjustl(mpiid_string)) //".dat"
#endif

  call fill_theta_rings(nside_search, theta_rings)
  
#ifdef MPI
    open(unit=outfile, file=outfilename)
  
    if (write_position_smax) then
       position_smax_filename = &
            trim(adjustl(inifilename(1:len_trim(inifilename)-4))) // &
            "_position_smax_mpi_" // trim(adjustl(mpiid_string)) // ".dat"
  
       open(unit=position_smax_file, file=trim(position_smax_filename))
  
       write(position_smax_file,'(A)') &
            "# alpha_deg pixel_i pixel_j S_ij beta_max_bin beta_max_deg"
    end if
#endif

  smax_angle = 0.
  smax_result = 0.
  smax_circle1_pixel = 0
  smax_circle2_pixel = 0
  smax_beta_max_bin = 0
  rho_hist_thread = 0_8
  rho_hist = 0_8
  rho_hist_result = 0_8

  do iopening_angle = nside_search/2 + mpiid, 2*nside_search, mpisize
     circle1%opening_angle = theta_rings(iopening_angle)

     if ((rad_to_deg(circle1%opening_angle) .lt. opening_angle_range%min) .or. &
         (rad_to_deg(circle1%opening_angle) .gt. opening_angle_range%max)) cycle

     circle2%opening_angle = circle1%opening_angle
     smax_angle(iopening_angle) = 0.

#ifdef MPI
     if (verbosity > 0) write(*,*) "MPI process", mpiid, &
          & " filling circle transforms, opening angle = ", &
          & circle1%opening_angle/Pi*180.
#else
     if (verbosity > 0) write(*,*) &
          & "filling circle transforms, opening angle = ", &
          & circle1%opening_angle/Pi*180.
#endif

     if (use_sw_filter) then
        call fill_circle_transforms( &
             & npix_search, circle1, nside_map, &
             & sw_filtered_map, obs_mask, mask_threshold, &
             & use_pixel_average, nest2ring_table, &
             & circle_transforms, &
             & masked_circle_transforms, &
             & mask_transforms)
     else
        call fill_circle_transforms( &
             & npix_search, circle1, nside_map, &
             & map, obs_mask, mask_threshold, &
             & use_pixel_average, nest2ring_table, &
             & circle_transforms, &
             & masked_circle_transforms, &
             & mask_transforms)
     end if

     call fill_master_auto_cache( &
          & npix_search, circle1, &
          & masked_circle_transforms, &
          & mask_transforms, &
          & master_auto_f_reg, &
          & master_auto_denominator, &
          & master_auto_modes, &
          & master_auto_valid)

#ifdef MPI
     if (verbosity > 0) write(*,*) "MPI process", mpiid, " comparing circles"
#else
     if (verbosity > 0) write(*,*) "comparing circles"
#endif
     smax_thread = 0.
     analysed_circles = 0
     psd_fail_modes = 0_8
     bestpair%smax = 0.

     if (write_every_nth_smax) then
        call open_every_nth_smax_file(every_nth_smax_file, &
             & circle1%opening_angle, every_nth_smax, inifilename)
     end if

#ifdef MPI
     starttime = mpi_wtime()
#else
     call cpu_time(starttime)
#endif
     !$omp parallel private(i, smax, thread_id, opposite_number, &
     !$omp& opposite_list, opposite_vector, bestcircle1, bestcircle2, &
     !$omp& position_smax, position_circle2, position_beta_max_bin, &
     !$omp& position_psd_fail_modes) &
     !$omp& firstprivate(circle1, circle2) &
     !$omp& reduction(+:analysed_circles, psd_fail_modes) &
     !$omp& default(shared)

     smax = 0.

     allocate(opposite_list(0:npix_search-1))

     !$omp do schedule(dynamic)
     do i = 0, npix_search-1

        circle1%pixel = i
        opposite_vector(1:3) = -circle_center_vectors(1:3, circle1%pixel)

        call query_disc(nside_search, opposite_vector, &
             & dble(Pi/180.*min_separation_angle), opposite_list, &
             & opposite_number)

#ifdef MPI
        if (verbosity > 1) write(*,*) "MPI process ", mpiid, " found ", &
             & opposite_number, " opposing circles"
#else
        if (verbosity > 1) write(*,*) "found ", opposite_number, &
             & " opposing circles"
#endif
        call analyse_circle_pairs( &
             & circle1, circle2, circle_transforms, &
             & masked_circle_transforms, mask_transforms, &
             & opposite_list, opposite_number, &
             & smax, bestcircle1, bestcircle2, &
             & position_smax, position_circle2, &
             & position_beta_max_bin, analysed_circles, &
             & every_nth_smax_file, &
             & position_psd_fail_modes)

        psd_fail_modes = &
             psd_fail_modes + position_psd_fail_modes

        if (write_position_smax .and. &
            & position_beta_max_bin >= 0) then

           !$omp critical(position_smax_write)
           write(position_smax_file, *) &
                rad_to_deg(dble(circle1%opening_angle)), &
                circle1%pixel, &
                position_circle2%pixel, &
                position_smax, &
                position_beta_max_bin, &
                360.0_dp * dble(position_beta_max_bin) / &
                dble(nphi_ring)
           !$omp end critical(position_smax_write)
        end if

     end do
     !$omp end do

#ifdef _OPENMP
     thread_id = omp_get_thread_num() + 1
#else
     thread_id = 1
#endif

     if (thread_id .gt. MAX_THREAD_NUM) then
        write(*,*) "too many threads, please adjust MAX_THREAD_NUM and recompile"
        stop
     end if

     smax_thread(thread_id) = smax
     bestpair(thread_id)%smax = smax
     bestpair(thread_id)%circle1 = bestcircle1
     bestpair(thread_id)%circle2 = bestcircle2

     deallocate(opposite_list)
     !$omp end parallel

     expected_psd_modes = &
          int(analysed_circles, kind=8) * &
          int(nphi_ring/2 - m_min + 1, kind=8)

     if (psd_fail_modes > expected_psd_modes) then
        write(*,*) "ERROR: PSD failure count exceeds total checked modes."
        stop
     end if

     if (expected_psd_modes > 0_8) then
        psd_fail_fraction = &
             dble(psd_fail_modes) / &
             dble(expected_psd_modes)
     else
        psd_fail_fraction = 0.0_dp
     end if

#ifdef MPI
     stoptime = mpi_wtime()
#else
     call cpu_time(stoptime)
#endif

     write(*,*) &
          "PSD fail_fraction =", psd_fail_fraction
     if (write_every_nth_smax) close(every_nth_smax_file)

     smax_angle(iopening_angle) = maxval(smax_thread)
     smax_pos = maxloc(smax_thread, dim=1)

     smax_circle1_pixel(iopening_angle) = &
          & bestpair(smax_pos)%circle1%pixel + 1

     smax_circle2_pixel(iopening_angle) = &
          & bestpair(smax_pos)%circle2%pixel + 1

     call get_master_smax( &
          & bestpair(smax_pos)%circle1, &
          & bestpair(smax_pos)%circle2, &
          & masked_circle_transforms, &
          & mask_transforms, &
          & smax, beta_max_bin, final_valid_pair, &
          & selected_f_reg_ij, &
          & selected_f_reg_ii, &
          & selected_f_reg_jj, &
          & selected_f_lost)

     if (.not. final_valid_pair) then
        write(*,*) "ERROR: selected MASTER pair is invalid."
        stop
     end if

     smax_beta_max_bin(iopening_angle) = beta_max_bin

#ifdef MPI
     if (write_f_lost) then
        write(outfile,*) rad_to_deg(circle1%opening_angle), &
             & smax_angle(iopening_angle), &
             & selected_f_reg_ij, &
             & selected_f_reg_ii, &
             & selected_f_reg_jj, &
             & selected_f_lost
     else
        write(outfile,*) rad_to_deg(circle1%opening_angle), &
             & smax_angle(iopening_angle), &
             & selected_f_reg_ij, &
             & selected_f_reg_ii, &
             & selected_f_reg_jj
     end if
     
     call flush(outfile)
#endif

     if (verbosity .gt. 0) then
#ifdef MPI
        write(*,*) "MPI process", mpiid, " examined ", analysed_circles, &
             & " circles", " at a rate of ", &
             & (stoptime-starttime)/dble(analysed_circles), &
             & " seconds/circle"
#else
        write(*,*) "examined ", analysed_circles, " circles", &
             & " at a rate of ", &
             & (stoptime-starttime)/dble(analysed_circles), &
             & " seconds/circle"
#endif
        write(*,*) rad_to_deg(circle1%opening_angle), &
             & smax_angle(iopening_angle), mpiid
     end if
  end do

  if (write_rho_hist) then
     rho_hist = &
          sum(rho_hist_thread, dim=2)
  end if

! Gather the maxima from all MPI ranks.
#ifdef MPI
  close(outfile)

  if (write_position_smax) then
     close(position_smax_file)
  end if

  call MPI_reduce(smax_angle, smax_result, size(smax_angle), &
       & MPI_real, MPI_sum, 0, MPI_COMM_WORLD, ierr)

  call MPI_reduce(smax_circle1_pixel, smax_circle1_pixel_result, &
       & size(smax_circle1_pixel), MPI_integer, MPI_sum, 0, &
       & MPI_COMM_WORLD, ierr)

  call MPI_reduce(smax_circle2_pixel, smax_circle2_pixel_result, &
       & size(smax_circle2_pixel), MPI_integer, MPI_sum, 0, &
       & MPI_COMM_WORLD, ierr)

  call MPI_reduce(smax_beta_max_bin, smax_beta_max_bin_result, &
       & size(smax_beta_max_bin), MPI_integer, MPI_sum, 0, &
       & MPI_COMM_WORLD, ierr)

  if (write_rho_hist) then
     call MPI_reduce( &
          & rho_hist, rho_hist_result, n_rho_bins, &
          & MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  end if

  if (mpiid .eq. 0) then

     smax_circle1_pixel_result = smax_circle1_pixel_result - 1
     smax_circle2_pixel_result = smax_circle2_pixel_result - 1
  
       open(unit=full_outfile, &
            & file=trim(adjustl(inifilename(1:len_trim(inifilename)-4)))//".dat")
  
       do i = nside_search/2, 2*nside_search
  
          if ((rad_to_deg(dble(theta_rings(i))) .lt. &
               & opening_angle_range%min) .or. &
              (rad_to_deg(dble(theta_rings(i))) .gt. &
               & opening_angle_range%max)) cycle
  
          call write_pair_info(full_outfile, theta_rings(i), smax_result(i), &
               & smax_circle1_pixel_result(i), &
               & smax_circle2_pixel_result(i), &
               & smax_beta_max_bin_result(i), nside_search)
  
       end do
     close(full_outfile)

     if (write_rho_hist) then
        rho_hist_filename = &
             trim(adjustl( &
             inifilename(1:len_trim(inifilename)-4))) // &
             "_rho_hist.dat"
        open( &
             unit=rho_hist_file, &
             file=trim(rho_hist_filename))

        write(rho_hist_file,'(A)') &
             "# rho_low rho_high rho_center count"

        do i = 1, n_rho_bins
           write(rho_hist_file,*) &
                rho_min + &
                dble(i - 1) * rho_bin_width, &
                rho_min + &
                dble(i) * rho_bin_width, &
                rho_min + &
                (dble(i) - 0.5_dp) * rho_bin_width, &
                rho_hist_result(i)
        end do
        close(rho_hist_file)
     end if
  end if
#endif

  call finalize()

#ifdef MPI
  call MPI_finalize(ierr)
#endif

contains

  subroutine write_pair_info(file, opening_angle, smax, pixel1, pixel2, beta_max_bin, nside)
    implicit none
    integer, intent(in) :: file
    real,    intent(in) :: opening_angle, smax
    integer, intent(in) :: pixel1, pixel2, beta_max_bin, nside
    
    write(file, *) rad_to_deg(dble(opening_angle)), &
         & smax, &
         & pixel1, &
         & rad_to_deg(pix2ang_table(pixel1)%theta), &
         & rad_to_deg(pix2ang_table(pixel1)%phi), &
         & pixel2, &
         & rad_to_deg(pix2ang_table(pixel2)%theta), &
         & rad_to_deg(pix2ang_table(pixel2)%phi), &
         & rad_to_deg(angle_between_pixels(nside, pixel1, pixel2)), &
         & beta_max_bin, &
         & 360._dp * dble(beta_max_bin) / dble(nphi_ring)
  end subroutine write_pair_info

  real(dp) function angle_between_pixels(nside,pixel1, pixel2)
    implicit none
    integer, intent(in) :: nside
    integer, intent(in) :: pixel1, pixel2

    real(dp),  dimension(3) :: x, y
    call pix2vec_ring(nside, pixel1, x)
    call pix2vec_ring(nside, pixel2, y)
    angle_between_pixels = acos(dot_product(x,y)/(sqrt(dot_product(x,x)*dot_product(y,y))))
  end function angle_between_pixels

  subroutine finalize()
    implicit none

#ifdef DFFT
    call dfftw_destroy_plan(plan_forward)
    call dfftw_destroy_plan(plan_backward)    
    call dfftw_destroy_plan(plan_master_forward)
    call dfftw_destroy_plan(plan_master_backward)
#else
    call sfftw_destroy_plan(plan_forward)
    call sfftw_destroy_plan(plan_backward)
    call sfftw_destroy_plan(plan_master_forward)
    call sfftw_destroy_plan(plan_master_backward)
#endif

    deallocate(temperature_in)
    deallocate(obs_mask)
    deallocate(temperature_out)
    deallocate(sm_out)
    deallocate(sm_in)
    deallocate(master_plan_in)
    deallocate(master_plan_out)
    deallocate(smax_angle)
    deallocate(smax_result)
    deallocate(smax_circle1_pixel)
    deallocate(smax_circle2_pixel)
    deallocate(smax_circle1_pixel_result)
    deallocate(smax_circle2_pixel_result)
    deallocate(smax_beta_max_bin)
    deallocate(smax_beta_max_bin_result)

    if (allocated(circle_transforms)) then
       deallocate(circle_transforms)
    end if

    if (allocated(masked_circle_transforms)) then
       deallocate(masked_circle_transforms)
    end if

    if (allocated(mask_transforms)) then
       deallocate(mask_transforms)
    end if

    if (allocated(master_auto_f_reg)) then
       deallocate(master_auto_f_reg)
    end if

    if (allocated(master_auto_denominator)) then
       deallocate(master_auto_denominator)
    end if

    if (allocated(master_auto_valid)) then
       deallocate(master_auto_valid)
    end if

    if (allocated(master_auto_modes)) then
       deallocate(master_auto_modes)
    end if
    
    if (allocated(sw_filtered_map)) then
       deallocate(sw_filtered_map)
    end if
    
    if (allocated(sw_filter)) then
       deallocate(sw_filter)
    end if

    deallocate(nest2ring_table)
    deallocate(ring2nest_table)
    deallocate(theta_rings)
  end subroutine finalize

  subroutine initialize(inifilename, nside, nphi_ring)
    use inifile, only : ini_open, ini_close, ini_read_int, ini_read_string, ini_read_logical, ini_read_real
    implicit none
    character(LEN=*), intent(in) :: inifilename
    integer, intent(out) :: nside, nphi_ring

    call ini_open(inifilename, 100)
    verbosity = ini_read_int("verbosity", 0)
    nside = ini_read_int("nside", 128)
    opening_angle_range%min=ini_read_real("opening_angle_min")
    opening_angle_range%max=ini_read_real("opening_angle_max")
    opening_angle_range%num_steps=ini_read_int("n_opening_angles", 1)
    use_pixel_average=ini_read_logical("use_pixel_average", .true.)
    nphi_ring = ini_read_int("nphi_ring", 2*nside)
    m_min = ini_read_int("m_min", 1)

    if (m_min < 1 .or. m_min > nphi_ring/2) then
       write(*,*) "ERROR: m_min must be between 1 and nphi_ring/2"
       stop
    end if

    orientable = ini_read_logical("orientable", .false.)
    nside_search = ini_read_int("nside_search", 128)
    min_separation_angle = 180.-ini_read_real("min_separation_angle", 179.5)
    mapfilename=ini_read_string("mapfilename")
    maskfilename = ini_read_string("maskfilename")
    mask_threshold = ini_read_real( &
         & "mask_threshold", 0.5)
    write_every_nth_smax=ini_read_logical("write_every_nth_smax", .false.)
    every_nth_smax=ini_read_int("every_nth_smax", 10000)
    write_position_smax=ini_read_logical("write_position_smax", .false.)
    write_f_lost = ini_read_logical("write_f_lost", .false.)
    random_s = ini_read_logical("random_s", .false.)
    scramble_sky = ini_read_logical("scramble_sky", .false.)

    use_sw_filter = ini_read_logical( &
         & "use_sw_filter", .false.)

    use_gcv_regularisation = ini_read_logical( &
         & "use_gcv_regularisation", .true.)    

    use_flost_regularisation = ini_read_logical( &
         & "use_flost_regularisation", .false.)
    
    rho_tolerance = ini_read_real( &
         & "rho_tolerance", 1.0)
    
    write_rho_hist = ini_read_logical( &
         & "write_rho_hist", .false.)
    
    if (rho_tolerance < 0.0_dp .or. &
        & rho_tolerance > 1.0_dp) then
    
       write(*,*) &
            "ERROR: rho_tolerance must be between 0 and 1."
       stop
    
    end if

    if (use_gcv_regularisation .and. &
        & use_flost_regularisation) then
    
       write(*,*) &
            "ERROR: GCV and f_lost regularisation cannot both be enabled."
       stop
    end if

    if (use_sw_filter) then
       sw_filter_filename = ini_read_string( &
            & "sw_filter_filename")
       sw_filter_lmax = ini_read_int( &
            & "sw_filter_lmax", 300)
    else
       sw_filter_filename = ""
       sw_filter_lmax = 300
    end if
    
    if (use_sw_filter) then
       if (sw_filter_lmax < 0 .or. &
           & sw_filter_lmax > 3*nside - 1) then
    
          write(*,*) &
               "ERROR: sw_filter_lmax must be between 0 and ", &
               & 3*nside - 1
          stop
       end if
    end if
    
    call ini_close()

    npix = 12*nside**2
    npix_search = 12*nside_search**2
    allocate(map(0:npix-1))
    allocate(obs_mask(0:npix-1))

    if (use_sw_filter) then
       allocate(sw_filtered_map(0:npix-1))
       allocate(sw_filter(0:sw_filter_lmax))
       sw_filtered_map = 0.0
       sw_filter = 0.0_dp 
    end if

    allocate(nest2ring_table(0:npix-1))
    allocate(ring2nest_table(0:npix-1))
    allocate(pix2ang_table(0:npix-1))
    allocate(theta_rings(nside_search/2:2*nside_search))
    allocate(smax_angle(nside_search/2:2*nside_search))
    allocate(smax_result(nside_search/2:2*nside_search))
    allocate(smax_circle1_pixel(nside_search/2:2*nside_search))
    allocate(smax_circle2_pixel(nside_search/2:2*nside_search))
    allocate(smax_circle1_pixel_result(nside_search/2:2*nside_search))
    allocate(smax_circle2_pixel_result(nside_search/2:2*nside_search))
    allocate(smax_beta_max_bin(nside_search/2:2*nside_search))
    allocate(smax_beta_max_bin_result(nside_search/2:2*nside_search))

    allocate(circle_center_vectors(3,0:npix_search-1))
    allocate(temperature_out(nphi_ring/2+1))
    allocate(sm_in(nphi_ring/2+1))
    allocate(temperature_in(nphi_ring))    
    allocate(sm_out(nphi_ring))
    allocate(master_plan_in(nphi_ring))
    allocate(master_plan_out(nphi_ring))
    allocate(circle_transforms( &
         & 1:nphi_ring/2+1, 0:npix_search-1))
    allocate(masked_circle_transforms( &
         & 1:nphi_ring/2+1, 0:npix_search-1))
    allocate(mask_transforms( &
         & 1:nphi_ring/2+1, 0:npix_search-1))

    allocate(master_auto_f_reg(0:npix_search-1))
    allocate(master_auto_denominator(0:npix_search-1))
    allocate(master_auto_valid(0:npix_search-1))
    allocate(master_auto_modes( &
     & 0:nphi_ring/2, 0:npix_search-1))

    master_auto_f_reg = 0.0_dp
    master_auto_denominator = 0.0_dp
    master_auto_valid = .false.

#ifdef DFFT
    call dfftw_plan_dft_r2c_1d(plan_forward, nphi_ring, temperature_in, temperature_out, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(plan_backward, nphi_ring, sm_in, sm_out, FFTW_ESTIMATE)

    call dfftw_plan_dft_1d( &
         & plan_master_forward, &
         & nphi_ring, &
         & master_plan_in, &
         & master_plan_out, &
         & FFTW_FORWARD, &
         & FFTW_ESTIMATE)

    call dfftw_plan_dft_1d( &
         & plan_master_backward, &
         & nphi_ring, &
         & master_plan_in, &
         & master_plan_out, &
         & FFTW_BACKWARD, &
         & FFTW_ESTIMATE)
#else
    call sfftw_plan_dft_r2c_1d(plan_forward, nphi_ring, temperature_in, temperature_out, FFTW_ESTIMATE)
    call sfftw_plan_dft_c2r_1d(plan_backward, nphi_ring, sm_in, sm_out, FFTW_ESTIMATE)

    call sfftw_plan_dft_1d( &
         & plan_master_forward, &
         & nphi_ring, &
         & master_plan_in, &
         & master_plan_out, &
         & FFTW_FORWARD, &
         & FFTW_ESTIMATE)

    call sfftw_plan_dft_1d( &
         & plan_master_backward, &
         & nphi_ring, &
         & master_plan_in, &
         & master_plan_out, &
         & FFTW_BACKWARD, &
         & FFTW_ESTIMATE)
#endif

    if (trim(adjustl(mapfilename)) .eq. "simsky.dat") then
       call read_simulated_sky(nside, scramble_sky, map)
    else
!      call read_map(nside, trim(adjustl(mapfilename)), map)
       call read_map_fits(nside, trim(adjustl(mapfilename)), map)
    end if

    if (use_sw_filter) then
       call read_d_filter(sw_filter_filename, sw_filter_lmax, sw_filter)
       call apply_d_filter(nside, sw_filter_lmax, map, sw_filter, sw_filtered_map)
    end if

    call read_map_fits( &
         & nside, &
         & trim(adjustl(maskfilename)), &
         & obs_mask)

    call initialize_nest2ring_table(nside, npix, nest2ring_table)
    call initialize_ring2nest_table(nside, npix, ring2nest_table)
    call initialize_pix2ang_table(nside_search, pix2ang_table, circle_center_vectors)

    opening_angle_range%stepsize=(opening_angle_range%max - opening_angle_range%min) / (1.*opening_angle_range%num_steps)

    if (min_separation_angle .lt. 0.5) then
       write(*,*) "Are you sure that you want to set min_separation_angle to such a small value (=", min_separation_angle, ")?"
       stop
    end if
  end subroutine initialize

  subroutine initialize_pix2ang_table(nside, pix2ang_table, circle_center_vectors)
    implicit none
    integer,      intent(in)  :: nside
    type(tangle), intent(out) :: pix2ang_table(0:)
    real,         intent(out) :: circle_center_vectors(1:,0:)
    real :: theta, phi
    integer i
    do i=0, 12*nside**2-1
       call pix2ang_ring(nside, i, pix2ang_table(i)%theta, pix2ang_table(i)%phi)
       theta=pix2ang_table(i)%theta
       phi=pix2ang_table(i)%phi
       circle_center_vectors(1, i)=sin(theta)*cos(phi)
       circle_center_vectors(2, i)=sin(theta)*sin(phi)
       circle_center_vectors(3, i)=cos(theta)
    end do
  end subroutine initialize_pix2ang_table

  subroutine initialize_ring2nest_table(nside, npix, ring2nest_table)
    implicit none
    integer,  intent(in) :: nside, npix
    integer,  intent(out):: ring2nest_table(0:)
    integer :: i

    !$omp parallel do
    do i=0, npix-1
       call ring2nest(nside, i, ring2nest_table(i))
    end do
    !$omp end parallel do
  end subroutine initialize_ring2nest_table

  real(dp) function rad_to_deg(theta)
    implicit none
    real(dp), intent(in) :: theta
    rad_to_deg=theta/pi * 180.
  end function rad_to_deg

  subroutine get_smax(circle1, circle2, circle_transforms, smax, random_s, beta_max_bin)
    implicit none
    type(tcircle), intent(in) :: circle1, circle2
    FFT_COMPLEX,   intent(in) :: circle_transforms(1:,0:)
    real(dp),      intent(out):: smax
    logical, optional, intent(in) :: random_s
    integer, optional, intent(out) :: beta_max_bin
    real(dp) denominator
    integer n
    FFT_COMPLEX, dimension(circle1%nphi_ring/2+1) :: smax_tmp_in, t_i, t_j
    real(FP), dimension(circle1%nphi_ring) :: smax_tmp_out

    t_i = circle_transforms(:,circle1%pixel)
    if (orientable) then
       t_j = circle_transforms(:,circle2%pixel)
    else
       t_j = conjg(circle_transforms(:,circle2%pixel))
    end if

    denominator = 0.0_dp
    smax_tmp_in = 0.0

    do n = m_min + 1, circle1%nphi_ring/2 + 1
       denominator = denominator + (n-1.) * &
            (abs(t_i(n))**2 + abs(t_j(n))**2)

       smax_tmp_in(n) = (n-1.) * t_i(n) * t_j(n)
    end do

    if (denominator <= 0.0_dp) then
       smax = -huge(1.0_dp)
       return
    end if

#ifdef DFFT
    call dfftw_execute_dft_c2r(plan_backward, smax_tmp_in, smax_tmp_out)
#else
    call sfftw_execute_dft_c2r(plan_backward, smax_tmp_in, smax_tmp_out)
#endif
    if (present(beta_max_bin)) then
       beta_max_bin = maxloc(smax_tmp_out, dim=1) - 1
    end if

    smax = maxval(smax_tmp_out) / denominator

    if (present(random_s)) then
       if (random_s) smax = smax_tmp_out(1)/denominator
    end if
  end subroutine get_smax

  subroutine get_f_lost_by_beta( &
       & circle1, circle2, mask_transforms, &
       & f_lost_by_beta, has_eligible_beta, &
       & joint_modes)

    implicit none

    type(tcircle), intent(in) :: circle1
    type(tcircle), intent(in) :: circle2

    FFT_COMPLEX, intent(in) :: mask_transforms(1:,0:)

    real(dp), intent(out) :: &
         & f_lost_by_beta(circle1%nphi_ring)

    logical, intent(out) :: has_eligible_beta

    integer :: nphi
    integer :: beta_bin
    integer :: n_joint

    FFT_COMPLEX, intent(out) :: &
         & joint_modes(circle1%nphi_ring/2 + 1)

    real(FP) :: &
         & joint_by_beta(circle1%nphi_ring)

    FFT_COMPLEX :: &
         & joint_modes_work(circle1%nphi_ring/2 + 1)

    nphi = circle1%nphi_ring

    ! Fourier product of the two binary mask rings.
    if (orientable) then
       ! Matching convention:
       ! partner = beta - k
       joint_modes = &
            mask_transforms(:,circle1%pixel) * &
            mask_transforms(:,circle2%pixel)
    else
       ! Matching convention:
       ! partner = k + beta
       joint_modes = &
            mask_transforms(:,circle1%pixel) * &
            conjg(mask_transforms(:,circle2%pixel))
    end if

    joint_modes_work = joint_modes

    ! Recover N_joint(beta) for every beta.
#ifdef DFFT
    call dfftw_execute_dft_c2r( &
         & plan_backward, joint_modes_work, joint_by_beta)
#else
    call sfftw_execute_dft_c2r( &
         & plan_backward, joint_modes_work, joint_by_beta)
#endif

    do beta_bin = 0, nphi - 1

       ! mask_transforms were divided by Nphi.
       ! Therefore the inverse FFT result must be
       ! multiplied by Nphi to recover the sample count.
       n_joint = nint( &
            real(joint_by_beta(beta_bin + 1), dp) * &
            dble(nphi))

       ! Protect against tiny FFT round-off errors.
       n_joint = max(0, min(nphi, n_joint))

       f_lost_by_beta(beta_bin + 1) = &
            1.0_dp - &
            dble(n_joint) / dble(nphi)

    end do

    has_eligible_beta = any( &
         f_lost_by_beta <= 0.5_dp)

  end subroutine get_f_lost_by_beta

  subroutine expand_r2c_spectrum( &
       & r2c_modes, nphi, full_modes)

    implicit none

    integer, intent(in) :: nphi

    FFT_COMPLEX, intent(in) :: &
         & r2c_modes(nphi/2 + 1)

    FFT_COMPLEX, intent(out) :: &
         & full_modes(nphi)

    integer :: m

    ! FFTW r2c stores the unique non-negative-frequency
    ! entries m = 0, ..., Nphi/2.
    full_modes(1:nphi/2 + 1) = &
         r2c_modes

    ! Explicitly store the negative-frequency entries
    ! implied by the real-valued input:
    ! X(-m) = conjugate[X(m)].
    do m = 1, nphi/2 - 1

       full_modes(nphi - m + 1) = &
            conjg(r2c_modes(m + 1))
    end do
  end subroutine expand_r2c_spectrum

  subroutine calculate_coupling_eigenvalues( &
       & nphi, &
       & M_column_ij, M_column_ii, M_column_jj, &
       & mu_ij, mu_ii, mu_jj)

    implicit none

    integer, intent(in) :: nphi

    ! First columns of the coupling matrices.
    FFT_COMPLEX, intent(in) :: &
         & M_column_ij(nphi)

    FFT_COMPLEX, intent(in) :: &
         & M_column_ii(nphi)

    FFT_COMPLEX, intent(in) :: &
         & M_column_jj(nphi)

    ! Coupling eigenvalues in the q basis.
    FFT_COMPLEX, intent(out) :: &
         & mu_ij(nphi)

    FFT_COMPLEX, intent(out) :: &
         & mu_ii(nphi)

    FFT_COMPLEX, intent(out) :: &
         & mu_jj(nphi)

    ! For a circulant matrix:
    !
    ! mu_q = DFT[M(:,0)]_q
#ifdef DFFT
    call dfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_ij, mu_ij)

    call dfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_ii, mu_ii)

    call dfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_jj, mu_jj)
#else
    call sfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_ij, mu_ij)

    call sfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_ii, mu_ii)

    call sfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_jj, mu_jj)
#endif
  end subroutine calculate_coupling_eigenvalues

  subroutine fill_master_auto_cache( &
       & npix, circle, &
       & masked_circle_transforms, &
       & mask_transforms, &
       & auto_f_reg, &
       & auto_denominator, &
       & auto_modes, &
       & auto_valid)

    implicit none
    integer, intent(in) :: npix
    type(tcircle), intent(in) :: circle

    FFT_COMPLEX, intent(in) :: &
         & masked_circle_transforms(1:,0:)

    FFT_COMPLEX, intent(in) :: &
         & mask_transforms(1:,0:)

    real(dp), intent(out) :: &
         & auto_f_reg(0:npix-1)

    real(dp), intent(out) :: &
         & auto_denominator(0:npix-1)

    real(FP), intent(out) :: &
         & auto_modes(0:,0:)

    logical, intent(out) :: &
         & auto_valid(0:npix-1)

    integer :: pixel
    integer :: nphi
    integer :: m

    real(dp) :: f_reg_auto
    real(dp) :: lambda_auto
    real(dp) :: denominator_auto
    real(dp) :: f_lost_auto

    logical :: boundary_auto

    FFT_COMPLEX :: &
         & masked_modes(circle%nphi_ring)

    FFT_COMPLEX :: &
         & window_modes(circle%nphi_ring)

    FFT_COMPLEX :: &
         & C_hat_auto(circle%nphi_ring)

    FFT_COMPLEX :: &
         & C_hat_F_auto(circle%nphi_ring)

    FFT_COMPLEX :: &
         & M_column_auto(circle%nphi_ring)

    FFT_COMPLEX :: &
         & mu_auto(circle%nphi_ring)

    FFT_COMPLEX :: &
         & C_corrected_auto(circle%nphi_ring)

    nphi = circle%nphi_ring

    ! Reset the complete cache for the current
    ! opening angle.
    auto_f_reg = 0.0_dp
    auto_denominator = 0.0_dp
    auto_valid = .false.

    ! Every pixel is an independent circle-centre
    ! auto problem.
    !$omp parallel do schedule(static) &
    !$omp& default(shared) &
    !$omp& private( &
    !$omp& pixel, m, &
    !$omp& masked_modes, window_modes, &
    !$omp& C_hat_auto, C_hat_F_auto, &
    !$omp& M_column_auto, &
    !$omp& mu_auto, C_corrected_auto, &
    !$omp& f_reg_auto, lambda_auto, &
    !$omp& f_lost_auto, &
    !$omp& denominator_auto, boundary_auto)

    do pixel = 0, npix - 1

       ! Complete the cyclic masked-temperature
       ! spectrum for this circle.
       call expand_r2c_spectrum( &
            & masked_circle_transforms(:,pixel), &
            & nphi, masked_modes)

       ! Complete the cyclic mask/window spectrum.
       call expand_r2c_spectrum( &
            & mask_transforms(:,pixel), &
            & nphi, window_modes)

       ! A completely masked circle has no usable
       ! auto inverse problem.
       if (maxval(abs(window_modes)) <= 0.0) cycle

       ! Masked auto pseudo-mode product:
       !
       ! C_hat_aa(m) =
       ! T_hat_a(m) conjugate[T_hat_a(m)].
       C_hat_auto = &
            masked_modes * conjg(masked_modes)

       ! First column of the auto-coupling matrix:
       !
       ! M_aa(m,0) =
       ! W_a(m) conjugate[W_a(m)].
       M_column_auto = &
            window_modes * conjg(window_modes)

       ! For the circulant auto-coupling matrix:
       !
       ! mu_aa(q) = DFT[M_aa(:,0)](q).
#ifdef DFFT
       call dfftw_execute_dft( &
            & plan_master_forward, &
            & M_column_auto, mu_auto)
#else
       call sfftw_execute_dft( &
            & plan_master_forward, &
            & M_column_auto, mu_auto)
#endif
       if (use_flost_regularisation) then
          ! Zero-frequency mask mode is the observed
          ! fraction of this individual circle.
          f_lost_auto = &
               1.0_dp - &
               real(mask_transforms(1,pixel), dp)

#ifdef DFFT
          call dfftw_execute_dft( &
               & plan_master_forward, &
               & C_hat_auto, C_hat_F_auto)
#else
          call sfftw_execute_dft( &
               & plan_master_forward, &
               & C_hat_auto, C_hat_F_auto)
#endif

          call apply_flost_regularisation( &
               & nphi, C_hat_F_auto, mu_auto, &
               & f_lost_auto, flost_A, &
               & C_corrected_auto)
          f_reg_auto = 0.0_dp
       else
          call select_regularisation( &
               & nphi, C_hat_auto, mu_auto, &
               & C_hat_F_auto, &
               & f_reg_auto, lambda_auto, &
               & boundary_auto)

          call apply_tikhonov_correction( &
               & nphi, C_hat_F_auto, mu_auto, &
               & lambda_auto, C_corrected_auto)
       end if

       ! Cache exactly the auto contribution used
       ! in the S-statistic denominator:
       !
       ! D_aa =
       ! sum_{m=m_min}^{Nphi/2}
       ! m Re[C_corrected_aa(m)].
       denominator_auto = 0.0_dp
       
       do m = m_min, nphi/2
       
          auto_modes(m,pixel) = &
               real(C_corrected_auto(m + 1), FP)
       
          denominator_auto = &
               denominator_auto + &
               dble(m) * &
               real( &
                    C_corrected_auto(m + 1), &
                    dp &
               )
       end do

       ! Store by circle-centre pixel.
       auto_f_reg(pixel) = f_reg_auto

       auto_denominator(pixel) = &
            denominator_auto

       auto_valid(pixel) = .true.
    end do

    !$omp end parallel do

  end subroutine fill_master_auto_cache

  subroutine get_master_smax( &
       & circle1, circle2, &
       & masked_circle_transforms, &
       & mask_transforms, &
       & smax, beta_max_bin, valid_pair, &
       & selected_f_reg_ij_out, &
       & selected_f_reg_ii_out, &
       & selected_f_reg_jj_out, &
       & selected_f_lost_out, &
       & psd_fail_modes_out)

    implicit none

    type(tcircle), intent(in) :: circle1
    type(tcircle), intent(in) :: circle2

    FFT_COMPLEX, intent(in) :: &
         & masked_circle_transforms(1:,0:)

    FFT_COMPLEX, intent(in) :: &
         & mask_transforms(1:,0:)

    real(dp), intent(out) :: smax

    integer, intent(out) :: beta_max_bin

    logical, intent(out) :: valid_pair

    real(dp), optional, intent(out) :: &
         & selected_f_reg_ij_out

    real(dp), optional, intent(out) :: &
         & selected_f_reg_ii_out

    real(dp), optional, intent(out) :: &
         & selected_f_reg_jj_out

     real(dp), optional, intent(out) :: &
          & selected_f_lost_out

     integer(8), optional, intent(out) :: &
          & psd_fail_modes_out

     integer :: nphi
     integer :: m
     integer :: support_index
     integer :: beta_index
     integer :: bin_index
     integer :: rho_bin_index
     integer :: rho_thread_id
     integer(8) :: psd_fail_local
     
     real(dp) :: auto_i
     real(dp) :: auto_j
     real(dp) :: cross_real
     real(dp) :: cross_imag
     real(dp) :: trace_g
     real(dp) :: determinant_g
     real(dp) :: rho_m

    real(dp) :: &
         & f_lost_by_beta(circle1%nphi_ring)

    FFT_COMPLEX :: &
         & joint_modes(circle1%nphi_ring/2 + 1)

    logical :: has_eligible_beta
    logical :: valid_s

    ! Stored non-negative-frequency cross product.
    FFT_COMPLEX :: &
         & C_hat_half_ij(circle1%nphi_ring/2 + 1)

    ! Complete cyclic cross pseudo-mode product.
    FFT_COMPLEX :: &
         & C_hat_ij(circle1%nphi_ring)

    FFT_COMPLEX :: &
         & C_hat_F_ij(circle1%nphi_ring)

    FFT_COMPLEX :: &
         & M_column_ij(circle1%nphi_ring)

    FFT_COMPLEX :: &
         & mu_ij(circle1%nphi_ring)

    FFT_COMPLEX :: &
         & C_corrected_ij(circle1%nphi_ring)

    real(dp) :: f_reg_ij
    real(dp) :: f_reg_ii
    real(dp) :: f_reg_jj

    real(dp) :: lambda_ij

    real(dp) :: denominator_ii
    real(dp) :: denominator_jj
    real(dp) :: support_f_lost

    logical :: boundary_ij
    real(dp) :: &
         & S_by_beta(circle1%nphi_ring)
    real(dp) :: &
         & S_support( &
         & circle1%nphi_ring, n_flost_support)
    
    FFT_COMPLEX :: &
         & C_corrected_support( &
         & circle1%nphi_ring, n_flost_support)

    real(dp) :: &
         & rho_modes(0:circle1%nphi_ring/2)
    
    real(dp) :: &
         & rho_support( &
         & 0:circle1%nphi_ring/2, n_flost_support)

    nphi = circle1%nphi_ring

    ! Default result for skipped or invalid pairs.
    smax = -huge(1.0_dp)
    beta_max_bin = -1
    valid_pair = .false.
    psd_fail_local = 0_8
    
    if (present(psd_fail_modes_out)) then
       psd_fail_modes_out = 0_8
    end if

    ! Calculate the lost matched-pair fraction
    ! separately for every beta.
    call get_f_lost_by_beta( &
         & circle1, circle2, &
         & mask_transforms, &
         & f_lost_by_beta, &
         & has_eligible_beta, &
         & joint_modes)

    ! Skip the circle pair if no beta satisfies
    ! f_lost(beta) <= 0.5.
    if (.not. has_eligible_beta) return

    ! The cached auto quantities are indexed by the
    ! circle-centre pixel. They were rebuilt for the
    ! current opening angle before the pair search.
    if ( &
         .not. master_auto_valid(circle1%pixel) .or. &
         .not. master_auto_valid(circle2%pixel) &
       ) return

    ! Read the correct cached auto quantities.
    f_reg_ii = &
         master_auto_f_reg(circle1%pixel)

    f_reg_jj = &
         master_auto_f_reg(circle2%pixel)

    denominator_ii = &
         master_auto_denominator(circle1%pixel)

    denominator_jj = &
         master_auto_denominator(circle2%pixel)

    ! Reuse the cross-window product already calculated
    ! for f_lost(beta).
    call expand_r2c_spectrum( &
         & joint_modes, nphi, M_column_ij)
    
    if (orientable) then
    
       ! Orientable cross product on the stored
       ! non-negative-frequency modes.
       C_hat_half_ij = &
            masked_circle_transforms(:,circle1%pixel) * &
            masked_circle_transforms(:,circle2%pixel)
    else
       ! Non-orientable cross product on the stored
       ! non-negative-frequency modes.
       C_hat_half_ij = &
            masked_circle_transforms(:,circle1%pixel) * &
            conjg(masked_circle_transforms(:,circle2%pixel))
    end if
    
    ! Complete the cyclic cross pseudo-mode product.
    call expand_r2c_spectrum( &
         & C_hat_half_ij, nphi, C_hat_ij)

    ! Because the coupling matrix is circulant,
    ! its eigenvalues are the DFT of its first column:
    !
    ! mu_ij(q) = DFT[M_ij(:,0)](q).
#ifdef DFFT
    call dfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_ij, mu_ij)
#else
    call sfftw_execute_dft( &
         & plan_master_forward, &
         & M_column_ij, mu_ij)
#endif

    if (use_flost_regularisation) then
       f_reg_ij = 0.0_dp

       ! Transform the cross pseudo-product into the
       ! coupling eigenbasis once for all five bins.
#ifdef DFFT
       call dfftw_execute_dft( &
            & plan_master_forward, &
            & C_hat_ij, C_hat_F_ij)
#else
       call sfftw_execute_dft( &
            & plan_master_forward, &
            & C_hat_ij, C_hat_F_ij)
#endif

       ! One reconstruction for each f_lost bin.
       do support_index = 1, n_flost_support
          support_f_lost = &
               flost_support(support_index)

          call apply_flost_regularisation( &
               & nphi, C_hat_F_ij, mu_ij, &
               & support_f_lost, flost_A, &
               & C_corrected_support(:,support_index))

          ! Calculate rho_m once for this reconstructed support,
          ! only when rho filtering is actually requested.
          if (rho_tolerance < 1.0_dp .or. &
              & write_rho_hist) then
             call calculate_rho_modes( &
                  & nphi, &
                  & circle1%pixel, &
                  & circle2%pixel, &
                  & C_corrected_support(:,support_index), &
                  & rho_support(:,support_index))
          end if

          call calculate_master_s_by_beta( &
               & nphi, &
               & C_corrected_support(:,support_index), &
               & denominator_ii, &
               & denominator_jj, &
               & circle1%pixel, &
               & circle2%pixel, &
               & rho_support(:,support_index), &
               & S_support(:,support_index), &
               & valid_s)
          if (.not. valid_s) return
       end do

       ! For each beta, use the S(beta) curve
       ! corresponding to its own f_lost bin.
       do beta_index = 1, nphi
          if (f_lost_by_beta(beta_index) > 0.5_dp) then
             S_by_beta(beta_index) = &
                  -huge(1.0_dp)
          else
             bin_index = &
                  min( &
                  n_flost_support, &
                  int( &
                  f_lost_by_beta(beta_index) / &
                  0.1_dp) + 1)

             S_by_beta(beta_index) = &
                  S_support(beta_index,bin_index)
          end if
       end do
    else
       ! Existing GCV or no-regularisation route.
       call select_regularisation( &
            & nphi, C_hat_ij, mu_ij, &
            & C_hat_F_ij, &
            & f_reg_ij, lambda_ij, boundary_ij)

       call apply_tikhonov_correction( &
            & nphi, C_hat_F_ij, mu_ij, &
            & lambda_ij, C_corrected_ij)
       
       ! Calculate rho_m once only when S filtering needs it.
       if (rho_tolerance < 1.0_dp .or. &
           & write_rho_hist) then
          call calculate_rho_modes( &
               & nphi, &
               & circle1%pixel, &
               & circle2%pixel, &
               & C_corrected_ij, &
               & rho_modes)
       end if
       
       call calculate_master_s_by_beta( &
            & nphi, &
            & C_corrected_ij, &
            & denominator_ii, &
            & denominator_jj, &
            & circle1%pixel, &
            & circle2%pixel, &
            & rho_modes, &
            & S_by_beta, valid_s)

       if (.not. valid_s) return
       where (f_lost_by_beta > 0.5_dp)
          S_by_beta = -huge(1.0_dp)
       end where
    end if

    ! Return all three selected fractions when requested.
    ! For f_lost regularisation there is no scalar f_reg,
    ! so these values are zero.
    if (present(selected_f_reg_ij_out)) then
       selected_f_reg_ij_out = f_reg_ij
    end if

    if (present(selected_f_reg_ii_out)) then
       selected_f_reg_ii_out = f_reg_ii
    end if

    if (present(selected_f_reg_jj_out)) then
       selected_f_reg_jj_out = f_reg_jj
    end if

    ! Same beta maximisation as the original
    ! matched-circle S-statistic.
    beta_max_bin = &
         maxloc(S_by_beta, dim=1) - 1
    smax = &
         S_by_beta(beta_max_bin + 1)

    if (use_flost_regularisation) then
       bin_index = &
            min( &
            n_flost_support, &
            int( &
            f_lost_by_beta(beta_max_bin + 1) / &
            0.1_dp) + 1)
    
       C_corrected_ij = &
            C_corrected_support(:,bin_index)
    
       if (rho_tolerance < 1.0_dp .or. &
           & write_rho_hist) then
    
          rho_modes = &
               rho_support(:,bin_index)
       end if
    end if

    ! ==================================================
    ! PSD check and optional rho histogram
    ! ==================================================
    if (present(psd_fail_modes_out)) then
       if (write_rho_hist) then
#ifdef _OPENMP
          rho_thread_id = omp_get_thread_num() + 1
#else
          rho_thread_id = 1
#endif
          if (rho_thread_id > MAX_THREAD_NUM) then
             write(*,*) &
                  "too many threads, please adjust MAX_THREAD_NUM and recompile"
             stop
          end if
       end if
    
       do m = m_min, nphi/2
          auto_i = &
               real( &
                    master_auto_modes(m,circle1%pixel), &
                    dp &
               )
    
          auto_j = &
               real( &
                    master_auto_modes(m,circle2%pixel), &
                    dp &
               )
    
          cross_real = &
               real(C_corrected_ij(m + 1), dp)
    
          cross_imag = &
               real(aimag(C_corrected_ij(m + 1)), dp)
    
          trace_g = &
               auto_i + auto_j
    
          determinant_g = &
               auto_i * auto_j - &
               cross_real * cross_real - &
               cross_imag * cross_imag
    
          if (write_rho_hist) then
             rho_m = rho_modes(m)
    
             rho_bin_index = &
                  int((rho_m - rho_min) / rho_bin_width) + 1
    
             rho_bin_index = &
                  max(1, min(n_rho_bins, rho_bin_index))
    
             rho_hist_thread(rho_bin_index,rho_thread_id) = &
                  rho_hist_thread(rho_bin_index,rho_thread_id) + 1_8
          end if
    
          if (trace_g < 0.0_dp .or. &
              & determinant_g < 0.0_dp) then
    
             psd_fail_local = &
                  psd_fail_local + 1_8
          end if
       end do
       psd_fail_modes_out = &
            psd_fail_local
    end if

    if (present(selected_f_lost_out)) then
       selected_f_lost_out = &
            f_lost_by_beta(beta_max_bin + 1)
    end if
    valid_pair = .true.
  end subroutine get_master_smax

  subroutine select_regularisation( &
       & nphi, C_hat, mu, C_hat_F, &
       & selected_f_reg, &
       & selected_lambda, &
       & selected_at_boundary)

    implicit none

    integer, intent(in) :: nphi
    FFT_COMPLEX, intent(in) :: C_hat(nphi)
    FFT_COMPLEX, intent(in) :: mu(nphi)
    FFT_COMPLEX, intent(out) :: C_hat_F(nphi)
    real(dp), intent(out) :: selected_f_reg
    real(dp), intent(out) :: selected_lambda
    logical, intent(out) :: selected_at_boundary

    if (use_gcv_regularisation) then
       call select_gcv_regularisation( &
            & nphi, C_hat, mu, C_hat_F, &
            & selected_f_reg, &
            & selected_lambda, &
            & selected_at_boundary)
    else

#ifdef DFFT
       call dfftw_execute_dft( &
            & plan_master_forward, &
            & C_hat, C_hat_F)
#else
       call sfftw_execute_dft( &
            & plan_master_forward, &
            & C_hat, C_hat_F)
#endif
       selected_f_reg = 0.0_dp
       selected_lambda = 0.0_dp
       selected_at_boundary = .false.
    end if

  end subroutine select_regularisation

  subroutine select_gcv_regularisation( &
       & nphi, C_hat, mu, C_hat_F, &
       & selected_f_reg, &
       & selected_lambda, &
       & selected_at_boundary)
  
    implicit none
    integer, intent(in) :: nphi

    ! Masked pseudo-mode product:
    ! C_hat(m), m = 0, ..., Nphi - 1.
    FFT_COMPLEX, intent(in) :: &
         & C_hat(nphi)

    ! Coupling eigenvalues:
    ! mu(q), q = 0, ..., Nphi - 1.
    FFT_COMPLEX, intent(in) :: &
         & mu(nphi)

    FFT_COMPLEX, intent(out) :: &
         & C_hat_F(nphi)

    real(dp), intent(out) :: selected_f_reg
    real(dp), intent(out) :: selected_lambda
    logical, intent(out)  :: selected_at_boundary

    ! Original grid:
    ! 21 points from 10^-9 to 10^-7,
    ! separated by 0.1 dex.
    integer, parameter :: n_gcv_grid = 21

    ! Each extension adds ten new points,
    ! corresponding to one complete decade.
    integer, parameter :: n_extension_points = 10

    ! Lower side: at most three one-decade extensions.
    integer, parameter :: &
         & max_lower_gcv_extensions = 3

    ! Upper side: at most six one-decade extensions.
    integer, parameter :: &
         & max_upper_gcv_extensions = 6

    real(dp), parameter :: &
         & log10_f_reg_min = -9.0_dp

    real(dp), parameter :: &
         & log10_f_reg_max = -7.0_dp

    real(dp), parameter :: &
         & log10_grid_step = 0.1_dp

    integer :: grid_index
    integer :: q
    integer :: extension_number
    integer :: points_in_block
    integer :: selected_point_in_block
    integer :: extension_direction

    real(dp) :: log10_f_reg
    real(dp) :: f_reg

    real(dp) :: block_start
    real(dp) :: block_step
    real(dp) :: current_boundary

    real(dp) :: lambda
    real(dp) :: lambda_squared
    real(dp) :: maximum_mu

    real(dp) :: mu_squared(nphi)
    real(dp) :: C_hat_F_power(nphi)

    real(dp) :: filter_denominator
    real(dp) :: residual_filter
    real(dp) :: residual_power
    real(dp) :: mean_residual_filter
    real(dp) :: gcv_value
    real(dp) :: best_gcv_value

    ! Transform from the circle-mode-product index m
    ! into the coupling eigenbasis q:
    !
    ! C_hat_F(q) = DFT[C_hat](q).
#ifdef DFFT
    call dfftw_execute_dft( &
         & plan_master_forward, &
         & C_hat, C_hat_F)
#else
    call sfftw_execute_dft( &
         & plan_master_forward, &
         & C_hat, C_hat_F)
#endif

    ! These quantities depend on the current geometry,
    ! but not on the tested regularisation value.
    ! Calculate them only once for this GCV call.
    mu_squared = &
         dble(abs(mu))**2
    C_hat_F_power = &
         dble(abs(C_hat_F))**2

    ! Define:
    !
    ! lambda = f_reg max_q |mu(q)|.
    maximum_mu = &
         dble(maxval(abs(mu)))

    if (maximum_mu <= 0.0_dp) then
       stop "All coupling eigenvalues are zero."
    end if

    ! Initialise the global minimum over the original
    ! grid and every subsequently evaluated extension.
    best_gcv_value = huge(1.0_dp)

    selected_f_reg = &
         10.0_dp**log10_f_reg_min
    selected_lambda = &
         selected_f_reg * maximum_mu
    selected_at_boundary = .false.
    extension_direction = 0
    current_boundary = 0.0_dp

    ! extension_number = 0:
    !     original 21-point grid.
    !
    ! extension_number = 1, ..., 6:
    !     optional ten-point, one-decade extensions.
    
    do extension_number = 0, max_upper_gcv_extensions
       ! The loop permits six upper extensions,
       ! but only three lower extensions.
       if ( &
            extension_number > 0 .and. &
            extension_direction == -1 .and. &
            extension_number > &
            max_lower_gcv_extensions &
       ) then
          exit
       end if

       if (extension_number == 0) then
          ! Original grid:
          ! -9.0, -8.9, ..., -7.0.
          points_in_block = n_gcv_grid
          block_start = log10_f_reg_min
          block_step = log10_grid_step

          ! This preserves the original initial
          ! lower-boundary selection.
          selected_point_in_block = 1
       else
          points_in_block = n_extension_points

          ! Begin one step beyond the current boundary,
          ! so the old boundary is not recalculated.
          block_start = &
               current_boundary + &
               dble(extension_direction) * &
               log10_grid_step

          block_step = &
               dble(extension_direction) * &
               log10_grid_step

          ! Zero means no point in the new block has
          ! improved the accumulated global minimum.
          selected_point_in_block = 0
       end if

       do grid_index = 1, points_in_block
          log10_f_reg = &
               block_start + &
               dble(grid_index - 1) * &
               block_step

          f_reg = &
               10.0_dp**log10_f_reg

          lambda = &
               f_reg * maximum_mu

          lambda_squared = &
               lambda**2

          residual_power = 0.0_dp
          mean_residual_filter = 0.0_dp
          do q = 1, nphi
             filter_denominator = &
                  mu_squared(q) + &
                  lambda_squared

             ! Residual filter:
             !
             ! r_q = 1 - h_q
             !     = lambda^2 /
             !       (|mu_q|^2 + lambda^2).
             residual_filter = &
                  lambda_squared / &
                  filter_denominator

             ! Accumulate:
             !
             ! sum_q r_q^2 |C_hat_F(q)|^2.
             residual_power = &
                  residual_power + &
                  residual_filter**2 * &
                  C_hat_F_power(q)

             ! Accumulate:
             !
             ! sum_q r_q.
             mean_residual_filter = &
                  mean_residual_filter + &
                  residual_filter
          end do

          ! Parseval-normalised residual power:
          !
          ! (1/Nphi^2)
          ! sum_q r_q^2 |C_hat_F(q)|^2.
          residual_power = &
               residual_power / &
               dble(nphi)**2

          ! This equals:
          !
          ! 1 - tr(H)/Nphi.
          mean_residual_filter = &
               mean_residual_filter / &
               dble(nphi)
          if (mean_residual_filter <= 0.0_dp) then
             gcv_value = huge(1.0_dp)
          else
             ! Ordinary GCV:
             !
             ! V_GCV =
             ! residual_power /
             ! [1 - tr(H)/Nphi]^2.
             gcv_value = &
                  residual_power / &
                  mean_residual_filter**2
          end if
          ! Compare with the accumulated minimum over
          ! the original grid and all earlier extensions.
          if (gcv_value < best_gcv_value) then
             best_gcv_value = gcv_value
             selected_point_in_block = grid_index
             selected_f_reg = f_reg
             selected_lambda = lambda
          end if
       end do

       if (extension_number == 0) then
          if (selected_point_in_block == 1) then
             ! Original lower boundary 10^-9 selected.
             ! Extend toward smaller f_reg.
             extension_direction = -1
             current_boundary = log10_f_reg_min
             selected_at_boundary = .true.

          else if ( &
               selected_point_in_block == &
               n_gcv_grid) then
             ! Original upper boundary 10^-7 selected.
             ! Extend toward larger f_reg.
             extension_direction = 1
             current_boundary = log10_f_reg_max
             selected_at_boundary = .true.
          else
             ! The original minimum is interior,
             ! so no extension is required.
             selected_at_boundary = .false.
             exit
          end if
       else
          if ( &
               selected_point_in_block == &
               points_in_block) then
             ! The newest outermost point is the global
             ! minimum. Update the current boundary and
             ! continue to the next extension decade.
             current_boundary = &
                  block_start + &
                  dble(points_in_block - 1) * &
                  block_step
             selected_at_boundary = .true.
          else
             ! Either no new point improved the old
             ! minimum, or an interior extension point
             ! became the minimum. Stop extending.
             selected_at_boundary = .false.
             exit
          end if
       end if
    end do
  end subroutine select_gcv_regularisation

  subroutine apply_tikhonov_correction( &
       & nphi, C_hat_F, mu, lambda, C_corrected)

    implicit none
    integer, intent(in) :: nphi

  ! Masked pseudo-mode product already transformed
  ! into the coupling eigenbasis by GCV.
  FFT_COMPLEX, intent(in) :: &
       & C_hat_F(nphi)

    ! Coupling eigenvalues in the q basis.
    FFT_COMPLEX, intent(in) :: &
         & mu(nphi)

    ! Regularisation scale selected by GCV.
    real(dp), intent(in) :: lambda

    ! Corrected mode-product estimator,
    ! returned in the circle-mode index m.
    FFT_COMPLEX, intent(out) :: &
         & C_corrected(nphi)

    integer  :: q
    real(dp) :: lambda_squared
    real(dp) :: filter_denominator
    FFT_COMPLEX :: C_corrected_F(nphi)

    lambda_squared = lambda**2

    ! Apply the Tikhonov inverse response
    ! independently to every coupling eigenmode.

    do q = 1, nphi
       filter_denominator = &
            dble(abs(mu(q)))**2 + &
            lambda_squared

       ! C_corrected_F(q) =
       !
       ! conjugate[mu(q)] /
       ! (|mu(q)|^2 + lambda^2)
       ! times C_hat_F(q).
       if (filter_denominator > 0.0_dp) then
          C_corrected_F(q) = &
               conjg(mu(q)) * C_hat_F(q) / &
               filter_denominator
       else
          C_corrected_F(q) = 0.0
       end if
    end do

    ! Transform from the coupling eigenbasis q
    ! back to the circle-mode-product index m.
#ifdef DFFT
    call dfftw_execute_dft( &
         & plan_master_backward, &
         & C_corrected_F, C_corrected)
#else
    call sfftw_execute_dft( &
         & plan_master_backward, &
         & C_corrected_F, C_corrected)
#endif
    ! FFTW's backward complex transform does not
    ! include the inverse-DFT factor 1/Nphi.
    C_corrected = &
         C_corrected / dble(nphi)
  end subroutine apply_tikhonov_correction

  subroutine apply_flost_regularisation( &
       & nphi, C_hat_F, mu, f_lost, A, C_corrected)

    implicit none

    integer, intent(in) :: nphi
    FFT_COMPLEX, intent(in) :: &
         & C_hat_F(nphi)
    FFT_COMPLEX, intent(in) :: &
         & mu(nphi)
    real(dp), intent(in) :: f_lost
    real(dp), intent(in) :: A
    FFT_COMPLEX, intent(out) :: &
         & C_corrected(nphi)
    integer :: q
    real(dp) :: maximum_mu
    real(dp) :: mu_squared
    real(dp) :: spectral_weight
    real(dp) :: lambda_squared
    real(dp) :: filter_denominator

    FFT_COMPLEX :: &
         & C_corrected_F(nphi)

    if (f_lost < 0.0_dp .or. &
        & f_lost >= 1.0_dp) then

       stop "f_lost must satisfy 0 <= f_lost < 1."
    end if

    maximum_mu = &
         dble(maxval(abs(mu)))

    if (maximum_mu <= 0.0_dp) then
       stop "All coupling eigenvalues are zero."
    end if

    do q = 1, nphi
       mu_squared = &
            dble(abs(mu(q)))**2

       spectral_weight = &
            1.0_dp - &
            mu_squared / maximum_mu**2

       lambda_squared = &
            A * maximum_mu**2 * &
            f_lost / (1.0_dp - f_lost) * &
            spectral_weight

       filter_denominator = &
            mu_squared + lambda_squared

       if (filter_denominator > 0.0_dp) then
          C_corrected_F(q) = &
               conjg(mu(q)) * C_hat_F(q) / &
               filter_denominator
       else
          C_corrected_F(q) = 0.0
       end if
    end do

    ! Transform back to the circle-mode index m.
#ifdef DFFT
    call dfftw_execute_dft( &
         & plan_master_backward, &
         & C_corrected_F, C_corrected)
#else
    call sfftw_execute_dft( &
         & plan_master_backward, &
         & C_corrected_F, C_corrected)
#endif
    C_corrected = &
         C_corrected / dble(nphi)

  end subroutine apply_flost_regularisation

  subroutine calculate_rho_modes( &
       & nphi, pixel_i, pixel_j, &
       & C_corrected_ij, rho_out)
  
    implicit none
    integer, intent(in) :: nphi
    integer(dp), intent(in) :: pixel_i
    integer(dp), intent(in) :: pixel_j
    FFT_COMPLEX, intent(in) :: &
         & C_corrected_ij(nphi)
    real(dp), intent(out) :: &
         & rho_out(0:nphi/2)
  
    integer :: m
    real(dp) :: auto_i
    real(dp) :: auto_j
    real(dp) :: cross_real
    real(dp) :: cross_imag
    real(dp) :: trace_g
    real(dp) :: eigen_gap
    real(dp) :: lambda_minus
    real(dp) :: lambda_plus
    real(dp) :: eigen_scale

    rho_out = 0.0_dp
    do m = m_min, nphi/2
  
       auto_i = &
            real(master_auto_modes(m,pixel_i), dp)
  
       auto_j = &
            real(master_auto_modes(m,pixel_j), dp)
  
       cross_real = &
            real(C_corrected_ij(m + 1), dp)
  
       cross_imag = &
            real(aimag(C_corrected_ij(m + 1)), dp)
  
       trace_g = &
            auto_i + auto_j
  
       eigen_gap = &
            sqrt( &
            (auto_i - auto_j) * (auto_i - auto_j) + &
            4.0_dp * ( &
            cross_real * cross_real + &
            cross_imag * cross_imag &
            ) &
            )
  
       lambda_minus = &
            0.5_dp * (trace_g - eigen_gap)
  
       lambda_plus = &
            0.5_dp * (trace_g + eigen_gap)
  
       eigen_scale = &
            max(abs(lambda_minus), abs(lambda_plus))
  
       if (eigen_scale > 0.0_dp) then
          rho_out(m) = &
               lambda_minus / eigen_scale
       else
          rho_out(m) = 0.0_dp
       end if
    end do
  
  end subroutine calculate_rho_modes

  subroutine calculate_master_s_by_beta( &
       & nphi, &
       & C_corrected_ij, &
       & denominator_ii, &
       & denominator_jj, &
       & pixel_i, &
       & pixel_j, &
       & rho_in, &
       & S_by_beta, valid_s)
  
    implicit none
  
    integer, intent(in) :: nphi
    integer(dp), intent(in) :: pixel_i
    integer(dp), intent(in) :: pixel_j
  
    FFT_COMPLEX, intent(in) :: &
         & C_corrected_ij(nphi)
  
    real(dp), intent(in) :: &
         & rho_in(0:nphi/2)

    ! Cached weighted auto contributions:
    !
    ! denominator_ii =
    ! sum_m m Re[C_corrected_ii(m)]
    !
    ! denominator_jj =
    ! sum_m m Re[C_corrected_jj(m)]
    real(dp), intent(in) :: denominator_ii
    real(dp), intent(in) :: denominator_jj
    real(dp), intent(out) :: &
         & S_by_beta(nphi)
    logical, intent(out) :: valid_s

    integer :: m
    real(dp) :: denominator
    real(dp) :: auto_i
    real(dp) :: auto_j
    FFT_COMPLEX :: &
         & numerator_modes(nphi/2 + 1)

    real(FP) :: &
         & numerator_by_beta(nphi)

    numerator_modes = 0.0
    ! ==================================================
    ! No rho filtering: use the original calculation
    ! ==================================================
    if (rho_tolerance >= 1.0_dp) then
       denominator = &
            denominator_ii + denominator_jj
       do m = m_min, nphi/2
          numerator_modes(m + 1) = &
               real(m, FP) * &
               C_corrected_ij(m + 1)
       end do
    
    ! ==================================================
    ! Rho filtering
    ! ==================================================
    else
       denominator = 0.0_dp
       do m = m_min, nphi/2
          ! Reject the whole recovered mode m.
          if (abs(rho_in(m)) > rho_tolerance) cycle
    
          ! Accepted mode: numerator contribution.
          numerator_modes(m + 1) = &
               real(m, FP) * &
               C_corrected_ij(m + 1)
    
          ! Accepted mode: corresponding denominator
          ! contribution from the same m.
          auto_i = &
               real(master_auto_modes(m,pixel_i), dp)
    
          auto_j = &
               real(master_auto_modes(m,pixel_j), dp)
    
          denominator = &
               denominator + &
               dble(m) * (auto_i + auto_j)
       end do
    end if
    
    if (denominator <= 0.0_dp) then
       S_by_beta = -huge(1.0_dp)
       valid_s = .false.
       return
    end if

    ! Obtain the cross-correlation numerator
    ! simultaneously for every beta.
#ifdef DFFT
    call dfftw_execute_dft_c2r( &
         & plan_backward, &
         & numerator_modes, &
         & numerator_by_beta)
#else
    call sfftw_execute_dft_c2r( &
         & plan_backward, &
         & numerator_modes, &
         & numerator_by_beta)
#endif
    S_by_beta = &
         real(numerator_by_beta, dp) / denominator
    valid_s = .true.

  end subroutine calculate_master_s_by_beta



  subroutine fill_circle_transforms( &
       & npix, circle, nside_map, map, &
       & obs_mask, mask_threshold, &
       & use_pixel_average, nest2ring_table, &
       & circle_transforms, &
       & masked_circle_transforms, &
       & mask_transforms)

    implicit none

    integer, intent(in) :: npix

    type(tcircle), intent(inout) :: circle

    integer, intent(in) :: nside_map

    real, intent(in) :: map(0:)
    real, intent(in) :: obs_mask(0:)

    real(dp), intent(in) :: mask_threshold

    logical, intent(in) :: use_pixel_average

    integer, intent(in) :: &
         & nest2ring_table(0:)

    FFT_COMPLEX, intent(out) :: &
         & circle_transforms(1:,0:)

    FFT_COMPLEX, intent(out) :: &
         & masked_circle_transforms(1:,0:)

    FFT_COMPLEX, intent(out) :: &
         & mask_transforms(1:,0:)

    FFT_COMPLEX, dimension( &
         & circle%nphi_ring/2 + 1) :: &
         & plan_out, &
         & temp_out, &
         & masked_temp_out, &
         & mask_out

    real(FP), dimension( &
         & circle%nphi_ring) :: &
         & plan_in, &
         & temp_in, &
         & masked_temp_in, &
         & mask_in

    integer(8) :: plan

    integer :: i

    ! Create one real-to-complex FFT plan.
#ifdef DFFT
    call dfftw_plan_dft_r2c_1d( &
         & plan, &
         & int(circle%nphi_ring), &
         & plan_in, &
         & plan_out, &
         & FFTW_ESTIMATE)
#else
    call sfftw_plan_dft_r2c_1d( &
         & plan, &
         & int(circle%nphi_ring), &
         & plan_in, &
         & plan_out, &
         & FFTW_ESTIMATE)
#endif

    !$omp parallel &
    !$omp& private( &
    !$omp& i, &
    !$omp& temp_in, temp_out, &
    !$omp& masked_temp_in, masked_temp_out, &
    !$omp& mask_in, mask_out) &
    !$omp& firstprivate(circle)

    !$omp do
    do i = 0, npix - 1

       circle%pixel = i

       circle%theta = &
            pix2ang_table(circle%pixel)%theta

       circle%phi = &
            pix2ang_table(circle%pixel)%phi

    ! Sample the temperature and mask using the
    ! same circle geometry.
    call get_temperature_along_circle( &
         & circle, &
         & nside_map, &
         & map, &
         & use_pixel_average, &
         & nest2ring_table, &
         & temp_in, &
         & obs_mask, &
         & mask_in)

       ! Convert the sampled mask to a binary window:
       ! 1 = observed, 0 = masked.
       where ( &
            real(mask_in, dp) >= mask_threshold)
          mask_in = real(1, FP)

       elsewhere
          mask_in = real(0, FP)

       end where

       ! Apply the binary window to the temperature ring.
       masked_temp_in = &
            mask_in * temp_in

       ! Fourier transform:
       !   T
       !   wT
       !   w
#ifdef DFFT
       call dfftw_execute_dft_r2c( &
            & plan, temp_in, temp_out)

       call dfftw_execute_dft_r2c( &
            & plan, &
            & masked_temp_in, &
            & masked_temp_out)

       call dfftw_execute_dft_r2c( &
            & plan, mask_in, mask_out)
#else
       call sfftw_execute_dft_r2c( &
            & plan, temp_in, temp_out)

       call sfftw_execute_dft_r2c( &
            & plan, &
            & masked_temp_in, &
            & masked_temp_out)

       call sfftw_execute_dft_r2c( &
            & plan, mask_in, mask_out)
#endif

       ! Original unnormalised temperature FFT.
       ! This is unchanged from the original S pipeline.
       circle_transforms(:,i) = &
            temp_out

       ! Masked-temperature FFT uses the same
       ! unnormalised convention as the original S pipeline.
       masked_circle_transforms(:,i) = &
            masked_temp_out

       ! The mask/window coefficients contain 1/Nphi
       ! so that they are the convolution coefficients W_m.
       mask_transforms(:,i) = &
            mask_out / dble(circle%nphi_ring)

    end do
    !$omp end do

    !$omp end parallel

    ! Destroy the FFT plan after every circle is processed.
#ifdef DFFT
    call dfftw_destroy_plan(plan)
#else
    call sfftw_destroy_plan(plan)
#endif

  end subroutine fill_circle_transforms
  
  subroutine initialize_nest2ring_table(nside, npix, nest2ring_table)
    implicit none
    integer, intent(in)  :: nside, npix
    integer, intent(out) :: nest2ring_table(0:)

    integer :: i

    !$omp parallel do
    do i=0, npix-1
       call nest2ring(nside, i, nest2ring_table(i))
    end do
    !$omp end parallel do

  end subroutine initialize_nest2ring_table

  subroutine get_rotation_matrix(phi, theta, rotation_matrix)
    implicit none
    real(dp), intent(in) :: phi, theta
    real(dp), intent(out):: rotation_matrix(1:3,1:3)
    !<PMV>
    !original code
!!$    rotation_matrix(1,1) = cos(phi) ! rotates from new to old
!!$    rotation_matrix(1,2) = cos(theta)*sin(phi)
!!$    rotation_matrix(1,3) = sin(phi)*sin(theta)
!!$    rotation_matrix(2,1) = -sin(phi)
!!$    rotation_matrix(2,2) = cos(theta)*cos(phi)
!!$    rotation_matrix(2,3) = cos(phi)*sin(theta)
!!$    rotation_matrix(3,1) = 0.0
!!$    rotation_matrix(3,2) = -sin(theta)
!!$    rotation_matrix(3,3) = cos(theta)
    !</PMV>

    !<PM> 
    !rotate around y by theta and then around z by phi
    !</PMV>
    rotation_matrix(1,1) = cos(theta)*cos(phi)
    rotation_matrix(1,2) = -sin(phi)
    rotation_matrix(1,3) = sin(theta)*cos(phi)
    rotation_matrix(2,1) = cos(theta)*sin(phi)
    rotation_matrix(2,2) = cos(phi)
    rotation_matrix(2,3) = sin(theta)*sin(phi)
    rotation_matrix(3,1) = -sin(theta)
    rotation_matrix(3,2) = 0.
    rotation_matrix(3,3) = cos(theta)

  end subroutine get_rotation_matrix

  subroutine get_temperature_along_circle( &
       & circle, nside_map, map, pixel_average, &
       & nest2ring_table, t, obs_mask, mask_t, &
       & mask_threshold, valid)

    implicit none
    type(tcircle), intent(in) :: circle
    integer,  intent(in) :: nside_map
    real(sp), intent(in) :: map(0:)
    logical,  intent(in) :: pixel_average
    integer,  intent(in) :: nest2ring_table(0:)
    real(FP), intent(out) :: t(:)
    real(sp), intent(in) :: obs_mask(0:)
    real(FP), intent(out) :: mask_t(:)
    
    real(dp), intent(in), optional :: mask_threshold
    logical,  intent(out), optional :: valid(:)
    real(dp) :: rotation_matrix(1:3,1:3), n_ipix(1:3), n_kpix(1:3)
    real(dp) :: theta_k, phi_k, phi_ring, sinang
    integer jpix, iphi,m, numnb, jnest

    integer :: nbnest(8), nb(8)

    call get_rotation_matrix(circle%phi, circle%theta,rotation_matrix)

    sinang=sin(circle%opening_angle)
    n_ipix(3) = cos(circle%opening_angle)
    do iphi=1, circle%nphi_ring
       phi_ring=2.*PI*DBLE(iphi)/DBLE(circle%nphi_ring)
       n_ipix(1) = sinang*cos(phi_ring)
       n_ipix(2) = sinang*sin(phi_ring)
       n_kpix = matmul(rotation_matrix,n_ipix)
       theta_k = acos(n_kpix(3))
       if(theta_k.ge.pi) theta_k = theta_k -0.00000001

       phi_k = atan2(n_kpix(2),n_kpix(1))+Pi

       if( (theta_k .lt. 0) .or. (theta_k.gt.PI).or. (phi_k.gt.2.*PI) .or. (phi_k.lt.0)) then
          write(*,*) theta_k,phi_k
          write(*,*) n_kpix(1:3)          
          stop 'error in get_temperature_along_circle'
       endif

       call ang2pix_ring(nside_map,theta_k,phi_k,jpix)

       if (present(valid)) then
          valid(iphi) = dble(obs_mask(jpix)) >= mask_threshold
       end if

       if (pixel_average) then
          jnest=ring2nest_table(jpix)
          call neighbours_nest(nside_map,jnest,nbnest,numnb)
          do m = 1,numnb
             if(nbnest(m).gt.12*nside_map*nside_map) then
                write(*,*) 'error here,nb:',nbnest(m)
                stop
             endif
             nb(m)=nest2ring_table(nbnest(m))
          enddo

          if(numnb.eq.8) then
             t(iphi) = DBLE( &
                  (map(jpix)*2.6137104 + &
                  (map(nb(1))+map(nb(3))+map(nb(5))+map(nb(7)))*0.9583863776 + &
                  (map(nb(2))+map(nb(4))+map(nb(6))+map(nb(8)))*0.692777175) / &
                  9.21836461)

             mask_t(iphi) = DBLE( &
                  (obs_mask(jpix)*2.6137104 + &
                  (obs_mask(nb(1))+obs_mask(nb(3))+obs_mask(nb(5))+obs_mask(nb(7)))*0.9583863776 + &
                  (obs_mask(nb(2))+obs_mask(nb(4))+obs_mask(nb(6))+obs_mask(nb(8)))*0.692777175) / &
                  9.21836461)
          else
             t(iphi) = DBLE( &
                  (map(jpix)*2.6137104 + &
                  (map(nb(1))+map(nb(3))+map(nb(5))+map(nb(7)))*0.9583863776 + &
                  (map(nb(2))+map(nb(4))+map(nb(6)))*0.692777175) / &
                  8.525587431)

             mask_t(iphi) = DBLE( &
                  (obs_mask(jpix)*2.6137104 + &
                  (obs_mask(nb(1))+obs_mask(nb(3))+obs_mask(nb(5))+obs_mask(nb(7)))*0.9583863776 + &
                  (obs_mask(nb(2))+obs_mask(nb(4))+obs_mask(nb(6)))*0.692777175) / &
                  8.525587431)
          endif
       else
          t(iphi) = DBLE(map(jpix))
          mask_t(iphi) = DBLE(obs_mask(jpix))
       end if

    end do
  end subroutine get_temperature_along_circle

  subroutine read_d_filter(filename, lmax, filter)
  
    implicit none
  
    character(len=*), intent(in)  :: filename
    integer,          intent(in)  :: lmax
    real(dp),         intent(out) :: filter(0:lmax)
  
    integer             :: filter_file
    integer             :: ios
    integer             :: ell
    real(dp)            :: filter_value
    character(len=1000) :: line
  
    filter = 0.0_dp
  
    open(newunit=filter_file, file=trim(adjustl(filename)), &
         & status="old", action="read", iostat=ios)
  
    if (ios /= 0) then
       write(*,*) "ERROR: could not open D filter file:"
       write(*,*) trim(adjustl(filename))
       stop
    end if
  
    do
  
       read(filter_file, "(A)", iostat=ios) line
  
       if (ios < 0) exit
  
       if (ios > 0) then
          write(*,*) "ERROR while reading D filter file"
          stop
       end if
  
       line = adjustl(line)
  
       if (len_trim(line) == 0) cycle
       if (line(1:1) == "#") cycle
  
       read(line, *, iostat=ios) ell, filter_value
  
       if (ios /= 0) then
          write(*,*) "ERROR: invalid row in D filter file:"
          write(*,*) trim(line)
          stop
       end if
  
       if (ell >= 0 .and. ell <= lmax) then
          filter(ell) = filter_value
       end if
  
    end do
  
    close(filter_file)
  
  end subroutine read_d_filter

  subroutine apply_d_filter(nside, lmax, input_map, filter, filtered_map)
  
    implicit none
  
    integer,  intent(in)  :: nside, lmax
    real(sp), intent(in)  :: input_map(0:)
    real(dp), intent(in)  :: filter(0:lmax)
    real(sp), intent(out) :: filtered_map(0:)
  
    complex(spc), allocatable :: alm(:,:,:)
    integer :: ell
  
    allocate(alm(1, 0:lmax, 0:lmax))
  
    call map2alm(nside, lmax, lmax, input_map, alm)
  
    do ell = 0, lmax
       alm(1,ell,0:ell) = &
            & alm(1,ell,0:ell) * real(filter(ell), sp)
    end do
  
    call alm2map(nside, lmax, lmax, alm, filtered_map)
  
    deallocate(alm)
  
  end subroutine apply_d_filter

  subroutine read_simulated_sky(nside_upgraded, scramble_sky, map)
    use ran_tools
    implicit none
    integer,  intent(in) :: nside_upgraded
    logical,  intent(in) :: scramble_sky
    real(sp), intent(out):: map(0:)


    integer :: simskyfile = 7
    integer(i4b), parameter :: nside_simskymap = 128
    integer(i4b), parameter :: ellmax = 256
    ! smoothing with gaussian beam of FWHM sqrt(8*ln2)/lsmooth
    integer(i4b), parameter :: lsmooth = 337 

    integer :: nside_input
    real(sp), allocatable :: simskymap(:)
    complex(spc), allocatable :: alm(:,:,:)
    real(dp), allocatable :: w8ring(:,:)

    real(dp) :: cos_theta_cut(2)=(/0,0/)
    character(LEN=1000) filename
    integer l, m
    
    ! for scrambling the sky
    COMPLEX(SPC) alm_big(1:1,0:1024,0:1024)
    REAL(DP), DIMENSION(1:2*512,1) :: w8ring_big 
    complex(spc) :: ahold
    integer :: n
    integer :: idum=-2234

    nside_input=nside_simskymap


    allocate(alm(1,0:ellmax, 0:ellmax))

    allocate(simskymap(0:12*nside_input**2-1))
    allocate(w8ring(2*nside_input, 1))

    !    write(*,*) ' reading in  test data'
    filename = 'simsky.dat'
    open(unit=simskyfile,file=filename,status='old', form='formatted')
    do l=0,12*nside_input*nside_input-1,1
       read(simskyfile,*) simskymap(l)
    end do
    close(simskyfile)

    ! ********************** Upgrade the pixelation and filter  **************************

    w8ring = 1.d0
    call map2alm(nside_input, ellmax, ellmax, simskymap, alm, cos_theta_cut, w8ring) 

    !  A spatial filter can be implemented here lsmooth=337 is a 0.4 degree beam

    do l=0,ellmax,1
       do m=0,l,1
          alm(1,l,m)=alm(1,l,m)*exp(-dble(l*(l+1))/(dble(lsmooth**2))/2.0)
       end do
    end do
    call alm2map(nside_upgraded,ellmax,ellmax,alm,map)

    if (scramble_sky) then
       w8ring_big = 1.
       call map2alm(512, 1024, 1024, map, alm_big, cos_theta_cut, w8ring_big) !Turns input map into alm's
       do l=1,1024,1
          do m=1,l,1     !can't scramble the m=0's without messing up the C_l's
             ahold = alm_big(1,l,m)
             n = int(ran_mwc(idum)*dble(l))+1
             alm_big(1,l,m) = 1.*alm_big(1, l, n)
             alm_big(1,l,n) = ahold
          enddo
       enddo
       write(*,*) 'SCRAMBLED SKY'
       call alm2map(512,1024,1024,alm_big,map)
    end if

    deallocate(w8ring, simskymap, alm)

  end subroutine read_simulated_sky

  subroutine read_map(nside_map, mapfilename, map)
    implicit none
    integer,          intent(in) :: nside_map
    character(LEN=*), intent(in) :: mapfilename
    real(sp),         intent(out):: map(0:)

    integer, parameter    :: nside_input = 512
    real, allocatable     :: nobs(:)
    integer               :: npix, status, npix_input

    if (nside_map .ne. nside_input) then
       write(*,*) "ERROR in read_map(): nside_map != ", nside_input
       stop
    end if
    npix= 12*nside_map**2
    allocate(nobs(0:npix-1))
    call read_archive_map(trim(adjustl(mapfilename)), map, nobs, npix_input, status)
    if (status .ne. 0) then
       write(*,*) "ERROR in read_map(): could not read file ", trim(adjustl(mapfilename))
    end if
    if (npix_input .ne. npix) then 
       write(*,*) "ERROR in read_map(): not enough pixels in input map"
       stop
    end if
    !call convert_nest2ring(nside_map, map)
    deallocate(nobs)
  end subroutine read_map

  subroutine read_map_fits(nside, mapfilename, map)
    implicit none
    integer,          intent(in)  :: nside
    character(len=*), intent(in)  :: mapfilename
    real,             intent(out) :: map(:)

    integer           :: npix
    real, allocatable :: in_map(:,:)
    npix = nside2npix(nside)
    allocate(in_map(npix, 1))

    call input_map(trim(adjustl(mapfilename)), in_map, 12*nside**2, 1)

    map = in_map(:,1)
    !call convert_nest2ring(nside, map)

    deallocate(in_map)
  end subroutine read_map_fits

  subroutine fill_theta_rings(nside_search, theta_rings)
    implicit none
    integer, intent(in) :: nside_search
    real,    intent(out):: theta_rings(nside_search/2:2*nside_search)

    real :: dth1, dth2, cth
    integer :: ith

    dth1 = 1.d0 / (3.d0*DBLE(nside_search)**2)
    dth2 = 2.d0 / (3.d0*DBLE(nside_search))

    do ith = nside_search/2,2*nside_search      ! pi/8 <= cos theta < 1
       !        cos(theta) in the pixelisation scheme
       if (ith.lt.nside_search) then  ! polar cap (north)
          cth = 1.d0  - DBLE(ith)**2 * dth1
       else                   ! tropical band (north) + equator
          cth = DBLE(2*nside_search-ith) * dth2
       endif

       theta_rings(ith) = acos(cth)
    enddo

  end subroutine fill_theta_rings

  subroutine keep_best_pairs(circle1, circle2, smax, pairlist)
    implicit none
    type(tcircle),   intent(in)  :: circle1, circle2
    real,            intent(in)  :: smax
    type(tpair), intent(out) :: pairlist(:, :)

  end subroutine keep_best_pairs

  real function smax_min(alpha)
    implicit none
    real(dp), intent(in) :: alpha

    real, parameter  :: theta_c=0.5
    smax_min = 1./sqrt(2.*Pi*sin(alpha)/theta_c)
  end function smax_min

  subroutine analyse_circle_pairs( &
       & circle1, circle2, circle_transforms, &
       & masked_circle_transforms, mask_transforms, &
       & opposite_list, opposite_number, &
       & smax, bestcircle1, bestcircle2, &
       & position_smax, position_circle2, &
       & position_beta_max_bin, &
       & analysed_circles, every_nth_smax_file, &
       & psd_fail_modes_out)

    implicit none

    type(tcircle), intent(in)    :: circle1
    type(tcircle), intent(inout) :: circle2

    FFT_COMPLEX, intent(in) :: circle_transforms(1:,0:)
    FFT_COMPLEX, intent(in) :: masked_circle_transforms(1:,0:)
    FFT_COMPLEX, intent(in) :: mask_transforms(1:,0:)

    integer, intent(in) :: opposite_list(0:)
    integer, intent(in) :: opposite_number
    real(dp), intent(inout) :: smax
    type(tcircle), intent(inout) :: bestcircle1, bestcircle2
    real(dp),      intent(out) :: position_smax
    type(tcircle), intent(out) :: position_circle2
    integer,       intent(out) :: position_beta_max_bin
    integer(dp), intent(inout) :: analysed_circles
    integer, intent(in) :: every_nth_smax_file
    integer(8), intent(out) :: psd_fail_modes_out
    integer  :: j
    integer  :: beta_try
    real(dp) :: smax_try
    logical :: valid_pair
    integer(8) :: pair_psd_fail
    position_smax = -huge(1.0_dp)
    position_beta_max_bin = -1
    position_circle2%pixel = -1
    psd_fail_modes_out = 0_8
    do j = 0, opposite_number - 1

       circle2%pixel = opposite_list(j)

       if (.not. write_position_smax) then
          if (circle2%pixel .lt. circle1%pixel) cycle
       end if

       call get_master_smax( &
            & circle1, circle2, &
            & masked_circle_transforms, &
            & mask_transforms, &
            & smax_try, beta_try, valid_pair, &
            & psd_fail_modes_out=pair_psd_fail)

       if (.not. valid_pair) cycle
       psd_fail_modes_out = &
            psd_fail_modes_out + pair_psd_fail

       if (smax_try .gt. position_smax) then
          position_smax = smax_try
          position_circle2 = circle2
          position_beta_max_bin = beta_try
       end if

       analysed_circles = analysed_circles + 1

       if (write_every_nth_smax .and. &
           mod(analysed_circles, every_nth_smax) .eq. 0) then
          !$omp critical
          write(every_nth_smax_file, *) smax_try
          !$omp end critical
       end if

       if (smax_try .gt. smax) then
          smax = smax_try
          bestcircle1 = circle1
          bestcircle2 = circle2
       end if
    end do

  end subroutine analyse_circle_pairs

  subroutine open_every_nth_smax_file(every_nth_smax_file, angle, every_nth_smax, inifilename)
    implicit none
    integer,          intent(in) :: every_nth_smax_file, every_nth_smax
    real(8),          intent(in) :: angle
    character(LEN=*), intent(in) :: inifilename

    character(LEN=1000) :: filename, string_angle

    write(string_angle, '(g16.5)') angle/Pi * 180.
    filename = trim(adjustl(inifilename(1:len_trim(inifilename)-4)))//"_smax_statistics_"//trim(adjustl(string_angle))//".dat"
    open(unit=every_nth_smax_file, file=trim(adjustl(filename)))
    write(every_nth_smax_file, *) "# writing every ", every_nth_smax, "th circle of angle ", angle
  end subroutine open_every_nth_smax_file

end program circles
