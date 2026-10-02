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
  integer                   :: d_diff_max_separation = 2

  integer                   :: nside_map !512  !Healpix resolution of upgraded map
  integer                   :: nside_search
  integer                   :: m_min = 2
  integer                   :: iopening_angle, i, nphi_ring, npix,j, npix_search
  integer(dp)               :: analysed_circles
  integer(dp)               :: d_diff_analysed_pairs
  logical                   :: use_pixel_average = .true.
  logical                   :: use_d_statistics = .false.
  logical                   :: use_s_statistics = .true.
  logical                   :: use_d_diff_statistics = .false.
  logical                   :: write_d_diff_null_histogram = .false.
  logical                   :: use_d_diff_bias_correction = .false.
  logical                   :: use_d_diff_q_statistics = .false.
  integer                   :: d_diff_hist_nbins = 2000
  integer                   :: d_diff_hist_n_joint
  integer                   :: d_diff_hist_bin_out
  integer                   :: d_min_joint_samples = 1
  real(dp)                  :: mask_threshold = 0.5_dp
  integer                   :: d_filter_lmax = 300

  real,         allocatable :: map(:)                     ! size: npix = 12*nside_map**2
  real,         allocatable :: obs_mask(:)
  real,         allocatable :: d_map(:)
  real(dp),     allocatable :: d_filter(:)

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

  real(dp),     allocatable :: dmax_angle(:)
  real(dp),     allocatable :: dmax_result(:)
  integer,      allocatable :: dmax_circle1_pixel(:)
  integer,      allocatable :: dmax_circle2_pixel(:)
  integer,      allocatable :: dmax_circle1_pixel_result(:)
  integer,      allocatable :: dmax_circle2_pixel_result(:)
  integer,      allocatable :: dmax_beta_max_bin(:)
  integer,      allocatable :: dmax_beta_max_bin_result(:)
  integer,      allocatable :: dmax_n_joint_max(:)
  integer,      allocatable :: dmax_n_joint_max_result(:)

! Independent global D_diff search
  real(dp),     allocatable :: d_diff_max_angle(:)
  real(dp),     allocatable :: d_diff_max_result(:)
  integer,      allocatable :: d_diff_max_circle1_pixel(:)
  integer,      allocatable :: d_diff_max_circle2_pixel(:)
  integer,      allocatable :: d_diff_max_circle1_pixel_result(:)
  integer,      allocatable :: d_diff_max_circle2_pixel_result(:)
  integer,      allocatable :: d_diff_max_beta_bin(:)
  integer,      allocatable :: d_diff_max_beta_bin_result(:)
  integer,      allocatable :: d_diff_max_n_joint(:)
  integer,      allocatable :: d_diff_max_n_joint_result(:)

  ! Gaussian null lookup for Q standardisation.
  ! Indexed by (iopening_angle, N_joint).
  real(dp), allocatable :: d_diff_gaussian_mu(:,:)
  real(dp), allocatable :: d_diff_gaussian_sigma(:,:)
  logical,  allocatable :: d_diff_gaussian_valid(:,:)
  integer,  allocatable :: d_diff_gaussian_fit_method(:,:)

  integer,      allocatable :: opposite_list(:)           ! size: npix_search = 12*nside_search**2
  real,         allocatable :: circle_center_vectors(:,:) ! size: (3,npix_search)

  FFT_COMPLEX,  allocatable :: temperature_out(:)         ! size: nphi_ring/2+1
  FFT_COMPLEX,  allocatable :: sm_in(:)                   ! size: nphi_ring/2+1
  real(FP),     allocatable :: temperature_in(:)          ! size: nphi_ring
  real(FP),     allocatable :: sm_out(:)                  ! size: nphi_ring

  FFT_COMPLEX,  allocatable :: circle_transforms(:,:)     ! size: (nphi_ring/2+1, npix_search)

  ! Transforms for independent D search.
  FFT_COMPLEX, allocatable :: d_w_t_transforms(:,:)
  FFT_COMPLEX, allocatable :: d_w_transforms(:,:)
  FFT_COMPLEX, allocatable :: d_w_t_sq_transforms(:,:)

  ! Third index is the physical separation s.
  FFT_COMPLEX, allocatable :: d_diff_w_delta_t_transforms(:,:,:)
  FFT_COMPLEX, allocatable :: d_diff_w_transforms(:,:,:)
  FFT_COMPLEX, allocatable :: d_diff_w_delta_t_sq_transforms(:,:,:)

  ! Null counts indexed by:
  ! (N_joint^diff, D_ij^diff histogram bin)
  integer(kind=8), allocatable :: d_diff_null_hist_counts(:,:)

  integer(8)                :: plan_forward, plan_backward
  real(dp)                  :: smax
  real(dp)                  :: position_smax
  real(dp)                  :: position_dij
  integer                   :: position_n_joint
  real(dp)                  :: opposite_vector(3)
  integer                   :: opposite_number
  integer                   :: verbosity = 0
  integer                   :: outfile = 101
  integer                   :: full_outfile = 110
  integer                   :: d_full_outfile = 111
  integer                   :: position_smax_file = 202
  integer                   :: d_diff_null_hist_file
  logical                   :: orientable=.false.
  logical                   :: write_every_nth_smax=.false.
  logical                   :: write_position_smax=.false.
  integer                   :: every_nth_smax = 10000, every_nth_smax_file = 201
  logical                   :: scramble_sky = .false.
  logical                   :: run_symmetry_test = .false.
  real(dp)                  :: symmetry_center_theta_deg = 0.0_dp
  real(dp)                  :: symmetry_center_phi_deg = 0.0_dp
  logical                   :: symmetry_input_has_missing_pixels = .true.

  type(tangle_range)        :: opening_angle_range
  type(tcircle)             :: circle1, circle2, bestcircle1, bestcircle2
  type(tpair)               :: bestpair(MAX_THREAD_NUM)
  type(tcircle)             :: position_circle2

  ! Thread-local results for the independent D search.
  type(tcircle) :: d_best_circle1_thread(MAX_THREAD_NUM)
  type(tcircle) :: d_best_circle2_thread(MAX_THREAD_NUM)

  type(tcircle) :: d_best_circle1_local
  type(tcircle) :: d_best_circle2_local

  real(dp) :: d_max_thread(MAX_THREAD_NUM)
  real(dp) :: d_best_local

  integer :: d_beta_thread(MAX_THREAD_NUM)
  integer :: d_n_joint_thread(MAX_THREAD_NUM)

  integer :: d_beta_local
  integer :: d_n_joint_local
  integer :: d_best_thread

  logical :: d_found_thread(MAX_THREAD_NUM)
  logical :: d_found_local

  ! Thread-local results for the independent D_diff search.
  type(tcircle) :: d_diff_best_circle1_thread(MAX_THREAD_NUM)
  type(tcircle) :: d_diff_best_circle2_thread(MAX_THREAD_NUM)
  
  type(tcircle) :: d_diff_best_circle1_local
  type(tcircle) :: d_diff_best_circle2_local
  
  real(dp) :: d_diff_max_thread(MAX_THREAD_NUM)
  real(dp) :: d_diff_best_local
  
  integer :: d_diff_beta_thread(MAX_THREAD_NUM)
  integer :: d_diff_n_joint_thread(MAX_THREAD_NUM)
  
  integer :: d_diff_beta_local
  integer :: d_diff_n_joint_local
  integer :: d_diff_best_thread
  
  logical :: d_diff_found_thread(MAX_THREAD_NUM)
  logical :: d_diff_found_local

  character(LEN=1000)       :: inifilename, outfilename, mapfilename, maskfilename
  character(LEN=1000)       :: d_filter_filename
  character(LEN=1000)       :: d_diff_gaussian_lookup_filename
  character(LEN=1000)       :: position_smax_filename
  character(LEN=1000)       :: d_diff_null_hist_filename
  real(dp)                  :: smax_thread(MAX_THREAD_NUM)
  real(dp)                  :: min_separation_angle
  logical                   :: random_s = .false.
  integer                   :: thread_id, beta_max_bin
  integer                   :: position_beta_max_bin
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
  if (run_symmetry_test) then
   call run_symmetry_test_mode()
   call finalize()
#ifdef MPI
   call MPI_finalize(ierr)
#endif
   stop
  end if

#ifdef MPI
  write(mpiid_string,*) mpiid
  outfilename=trim(adjustl(inifilename(1:len_trim(inifilename)-4)))//"_mpi_"//trim(adjustl(mpiid_string)) //".dat"
#endif

  if (write_d_diff_null_histogram) then
  
#ifdef MPI
     d_diff_null_hist_filename = &
          trim(adjustl( &
          inifilename(1:len_trim(inifilename)-4))) // &
          "_d_diff_null_hist_mpi_" // &
          trim(adjustl(mpiid_string)) // ".dat"
#else
     d_diff_null_hist_filename = &
          trim(adjustl( &
          inifilename(1:len_trim(inifilename)-4))) // &
          "_d_diff_null_hist.dat"
#endif
  
     open( &
          newunit=d_diff_null_hist_file, &
          file=trim(d_diff_null_hist_filename), &
          status="replace", &
          action="write")
  
     write(d_diff_null_hist_file, *) &
          "# d_diff_hist_nbins = ", d_diff_hist_nbins
  
     write(d_diff_null_hist_file, '(A)') &
          "# alpha_deg n_joint_diff d_diff_bin count"
  end if

  call fill_theta_rings(nside_search, theta_rings)

  if (use_d_diff_q_statistics) then
     call load_d_diff_gaussian_lookup()
  end if
  
#ifdef MPI
    open(unit=outfile, file=outfilename)
  
    if (write_position_smax) then
       position_smax_filename = &
            trim(adjustl(inifilename(1:len_trim(inifilename)-4))) // &
            "_position_smax_mpi_" // trim(adjustl(mpiid_string)) // ".dat"
  
       open(unit=position_smax_file, file=trim(position_smax_filename))
  
       write(position_smax_file,'(A)') &
            "# alpha_deg pixel_i pixel_j Smax beta_max_bin beta_max_deg Dij_at_betaSij n_joint"
    end if
#endif

  smax_angle = 0.
  smax_result = 0.
  smax_circle1_pixel = 0
  smax_circle2_pixel = 0
  smax_beta_max_bin = 0

  if (use_d_statistics) then
     dmax_angle = 0.0_dp
     dmax_result = 0.0_dp
     dmax_circle1_pixel = 0
     dmax_circle2_pixel = 0
     dmax_beta_max_bin = 0
     dmax_n_joint_max = 0
  end if

  if (use_d_diff_statistics) then
     d_diff_max_angle = 0.0_dp
     d_diff_max_result = 0.0_dp
     d_diff_max_circle1_pixel = 0
     d_diff_max_circle2_pixel = 0
     d_diff_max_beta_bin = 0
     d_diff_max_n_joint = 0
  end if

  do iopening_angle = nside_search/2 + mpiid, 2*nside_search, mpisize

     circle1%opening_angle = theta_rings(iopening_angle)

     if ((rad_to_deg(circle1%opening_angle) .lt. opening_angle_range%min) .or. &
         (rad_to_deg(circle1%opening_angle) .gt. opening_angle_range%max)) cycle

     circle2%opening_angle = circle1%opening_angle
     smax_angle(iopening_angle) = 0.

     ! Start a new null histogram for this alpha.
     if (write_d_diff_null_histogram) then
        d_diff_null_hist_counts = 0_8
     end if

     if (use_s_statistics) then

#ifdef MPI
        if (verbosity > 0) write(*,*) "MPI process", mpiid, &
             & " filling circle transforms, opening angle = ", &
             & circle1%opening_angle/Pi*180.
#else
        if (verbosity > 0) write(*,*) &
             & "filling circle transforms, opening angle = ", &
             & circle1%opening_angle/Pi*180.
#endif
        call fill_circle_transforms( &
             & npix_search, circle1, nside_map, map, &
             & use_pixel_average, nest2ring_table, &
             & circle_transforms)
     end if

     if (use_d_statistics .and. .not. use_s_statistics) then
        call fill_d_transforms( &
             & npix_search, circle1, nside_map, d_map, &
             & use_pixel_average, nest2ring_table, &
             & obs_mask, mask_threshold, &
             & d_w_t_transforms, &
             & d_w_transforms, &
             & d_w_t_sq_transforms)
     end if

     if (use_d_diff_statistics) then
        call fill_d_diff_transforms( &
             & npix_search, circle1, nside_map, d_map, &
             & use_pixel_average, nest2ring_table, &
             & obs_mask, mask_threshold, &
             & d_diff_w_delta_t_transforms, &
             & d_diff_w_transforms, &
             & d_diff_w_delta_t_sq_transforms)
     end if

#ifdef MPI
     if (verbosity > 0) write(*,*) "MPI process", mpiid, " comparing circles"
#else
     if (verbosity > 0) write(*,*) "comparing circles"
#endif
     smax_thread = 0.
     analysed_circles = 0
     d_diff_analysed_pairs = 0
     bestpair%smax = 0.

     if (use_d_statistics .and. .not. use_s_statistics) then
        d_max_thread = -huge(1.0_dp)
        d_found_thread = .false.
        d_beta_thread = 0
        d_n_joint_thread = 0
     end if

     if (use_d_diff_statistics) then
        d_diff_max_thread = -huge(1.0_dp)
        d_diff_found_thread = .false.
        d_diff_beta_thread = 0
        d_diff_n_joint_thread = 0
     end if

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
     !$omp& position_dij, position_n_joint, &
     !$omp& d_best_local, d_best_circle1_local, &
     !$omp& d_best_circle2_local, d_beta_local, &
     !$omp& d_n_joint_local, d_found_local, &
     !$omp& d_diff_best_local, d_diff_best_circle1_local, &
     !$omp& d_diff_best_circle2_local, d_diff_beta_local, &
     !$omp& d_diff_n_joint_local, d_diff_found_local) &
     !$omp& firstprivate(circle1, circle2) &
     !$omp& reduction(+:analysed_circles, d_diff_analysed_pairs) &
     !$omp& default(shared)

     smax = 0.
     d_best_local = -huge(1.0_dp)
     d_beta_local = 0
     d_n_joint_local = 0
     d_found_local = .false.

     d_diff_best_local = -huge(1.0_dp)
     d_diff_beta_local = 0
     d_diff_n_joint_local = 0
     d_diff_found_local = .false.

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
        if (use_s_statistics) then

           call analyse_circle_pairs( &
                & circle1, circle2, circle_transforms, &
                & opposite_list, opposite_number, &
                & smax, bestcircle1, bestcircle2, &
                & position_smax, position_circle2, &
                & position_beta_max_bin, analysed_circles, &
                & every_nth_smax_file)

           if (write_position_smax .and. &
               & position_beta_max_bin >= 0) then

              position_dij = 0.0_dp
              position_n_joint = 0

              if (use_d_statistics) then
                 call get_d_at_beta( &
                      & circle1, &
                      & position_circle2, &
                      & position_beta_max_bin, &
                      & position_dij, &
                      & position_n_joint)
              end if

              !$omp critical(position_smax_write)
              write(position_smax_file, *) &
                   rad_to_deg(dble(circle1%opening_angle)), &
                   circle1%pixel, &
                   position_circle2%pixel, &
                   position_smax, &
                   position_beta_max_bin, &
                   360.0_dp * dble(position_beta_max_bin) / &
                   dble(nphi_ring), &
                   position_dij, &
                   position_n_joint
              !$omp end critical(position_smax_write)
           end if  
        end if

        if (use_d_statistics .and. .not. use_s_statistics) then
           call analyse_circle_pairs_d( &
                & circle1, circle2, &
                & opposite_list, opposite_number, &
                & d_best_local, &
                & d_best_circle1_local, &
                & d_best_circle2_local, &
                & d_beta_local, &
                & d_n_joint_local, &
                & d_found_local)
        end if

        if (opposite_number > 0) then
           d_diff_analysed_pairs = d_diff_analysed_pairs + int( &
                count( &
                opposite_list(0:opposite_number-1) >= circle1%pixel &
                ), kind=dp)
        end if

        if (use_d_diff_statistics) then

           call analyse_circle_pairs_d_diff( &
                & circle1, circle2, &
                & opposite_list, opposite_number, &
                & d_diff_best_local, &
                & d_diff_best_circle1_local, &
                & d_diff_best_circle2_local, &
                & d_diff_beta_local, &
                & d_diff_n_joint_local, &
                & d_diff_found_local)
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

     if (use_s_statistics) then
        smax_thread(thread_id) = smax
        bestpair(thread_id)%smax = smax
        bestpair(thread_id)%circle1 = bestcircle1
        bestpair(thread_id)%circle2 = bestcircle2
     end if

     if (use_d_statistics .and. .not. use_s_statistics &
         .and. d_found_local) then

        d_max_thread(thread_id) = d_best_local

        d_best_circle1_thread(thread_id) = &
             & d_best_circle1_local

        d_best_circle2_thread(thread_id) = &
             & d_best_circle2_local

        d_beta_thread(thread_id) = d_beta_local

        d_n_joint_thread(thread_id) = d_n_joint_local

        d_found_thread(thread_id) = .true.
     end if

     if (use_d_diff_statistics .and. d_diff_found_local) then
        d_diff_max_thread(thread_id) = d_diff_best_local

        d_diff_best_circle1_thread(thread_id) = &
             & d_diff_best_circle1_local

        d_diff_best_circle2_thread(thread_id) = &
             & d_diff_best_circle2_local

        d_diff_beta_thread(thread_id) = &
             & d_diff_beta_local

        d_diff_n_joint_thread(thread_id) = &
             & d_diff_n_joint_local

        d_diff_found_thread(thread_id) = .true.
     end if

     deallocate(opposite_list)
     !$omp end parallel
     
     if (use_d_statistics .and. .not. use_s_statistics) then

        if (any(d_found_thread)) then
           d_best_thread = maxloc( &
                & d_max_thread, &
                & dim=1, &
                & mask=d_found_thread)

           dmax_angle(iopening_angle) = &
                & d_max_thread(d_best_thread)

           dmax_circle1_pixel(iopening_angle) = &
                & d_best_circle1_thread(d_best_thread)%pixel + 1

           dmax_circle2_pixel(iopening_angle) = &
                & d_best_circle2_thread(d_best_thread)%pixel + 1

           dmax_beta_max_bin(iopening_angle) = &
                & d_beta_thread(d_best_thread)

           dmax_n_joint_max(iopening_angle) = &
                & d_n_joint_thread(d_best_thread)
        else
           dmax_angle(iopening_angle) = 0.0_dp
           dmax_circle1_pixel(iopening_angle) = 0
           dmax_circle2_pixel(iopening_angle) = 0
           dmax_beta_max_bin(iopening_angle) = 0
           dmax_n_joint_max(iopening_angle) = 0
        end if
     end if

     if (use_d_diff_statistics .and. verbosity > 0) then
        write(*,*) &
             "D_diff examined ", d_diff_analysed_pairs, &
             " unique (i,j) pairs at alpha = ", &
             rad_to_deg(dble(circle1%opening_angle))
     end if

#ifdef MPI
     stoptime = mpi_wtime()
#else
     call cpu_time(stoptime)
#endif

     if (use_d_diff_statistics .and. verbosity > 0) then
        write(*,*) &
             "D_diff pair-search elapsed = ", &
             stoptime - starttime, " seconds"

        write(*,*) &
             "D_diff seconds/pair = ", &
             (stoptime - starttime) / &
             dble(d_diff_analysed_pairs)
     end if

     if (write_every_nth_smax) close(every_nth_smax_file)

     if (use_s_statistics) then
        smax_angle(iopening_angle) = maxval(smax_thread)
        smax_pos = maxloc(smax_thread, dim=1)

        smax_circle1_pixel(iopening_angle) = &
             & bestpair(smax_pos)%circle1%pixel + 1

        smax_circle2_pixel(iopening_angle) = &
             & bestpair(smax_pos)%circle2%pixel + 1

        call get_smax( &
             & bestpair(smax_pos)%circle1, &
             & bestpair(smax_pos)%circle2, &
             & circle_transforms, smax, random_s, beta_max_bin)

        smax_beta_max_bin(iopening_angle) = beta_max_bin

        if (use_d_statistics) then
           call get_d_at_beta( &
                & bestpair(smax_pos)%circle1, &
                & bestpair(smax_pos)%circle2, &
                & beta_max_bin, &
                & dmax_angle(iopening_angle), &
                & dmax_n_joint_max(iopening_angle))

           dmax_circle1_pixel(iopening_angle) = &
                & smax_circle1_pixel(iopening_angle)

           dmax_circle2_pixel(iopening_angle) = &
                & smax_circle2_pixel(iopening_angle)

           dmax_beta_max_bin(iopening_angle) = beta_max_bin

        end if
     end if

     if (use_d_diff_statistics) then
        if (any(d_diff_found_thread)) then

           d_diff_best_thread = maxloc( &
                & d_diff_max_thread, &
                & dim=1, &
                & mask=d_diff_found_thread)

           d_diff_max_angle(iopening_angle) = &
                & d_diff_max_thread(d_diff_best_thread)

           d_diff_max_circle1_pixel(iopening_angle) = &
                & d_diff_best_circle1_thread( &
                & d_diff_best_thread)%pixel + 1

           d_diff_max_circle2_pixel(iopening_angle) = &
                & d_diff_best_circle2_thread( &
                & d_diff_best_thread)%pixel + 1

           d_diff_max_beta_bin(iopening_angle) = &
                & d_diff_beta_thread(d_diff_best_thread)

           d_diff_max_n_joint(iopening_angle) = &
                & d_diff_n_joint_thread(d_diff_best_thread)

        else

           d_diff_max_angle(iopening_angle) = 0.0_dp
           d_diff_max_circle1_pixel(iopening_angle) = 0
           d_diff_max_circle2_pixel(iopening_angle) = 0
           d_diff_max_beta_bin(iopening_angle) = 0
           d_diff_max_n_joint(iopening_angle) = 0

        end if

     end if

! Write the non-zero null-histogram counts for this alpha.
     if (write_d_diff_null_histogram) then
     
        do d_diff_hist_n_joint = d_min_joint_samples, nphi_ring
     
           do d_diff_hist_bin_out = 0, d_diff_hist_nbins - 1
     
              if (d_diff_null_hist_counts( &
                   & d_diff_hist_n_joint, &
                   & d_diff_hist_bin_out) == 0_8) cycle
     
              write(d_diff_null_hist_file, *) &
                   rad_to_deg(dble(circle1%opening_angle)), &
                   d_diff_hist_n_joint, &
                   d_diff_hist_bin_out, &
                   d_diff_null_hist_counts( &
                   & d_diff_hist_n_joint, &
                   & d_diff_hist_bin_out)
           end do
        end do
        call flush(d_diff_null_hist_file)
     end if

#ifdef MPI
     write(outfile,*) rad_to_deg(circle1%opening_angle), &
          & smax_angle(iopening_angle)
     call flush(outfile)
#endif

     if (verbosity .gt. 0 .and. use_s_statistics) then
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

        if (use_d_statistics) then
           write(*,*) "Dmax = ", dmax_angle(iopening_angle), &
                & " beta_bin = ", dmax_beta_max_bin(iopening_angle), &
                & " n_joint = ", dmax_n_joint_max(iopening_angle)
        end if
     end if

  end do

if (write_d_diff_null_histogram) then
   close(d_diff_null_hist_file)
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

  if (use_d_statistics) then

     call MPI_reduce(dmax_angle, dmax_result, size(dmax_angle), &
          & MPI_double_precision, MPI_sum, 0, MPI_COMM_WORLD, ierr)

     call MPI_reduce(dmax_circle1_pixel, dmax_circle1_pixel_result, &
          & size(dmax_circle1_pixel), MPI_integer, MPI_sum, 0, &
          & MPI_COMM_WORLD, ierr)

     call MPI_reduce(dmax_circle2_pixel, dmax_circle2_pixel_result, &
          & size(dmax_circle2_pixel), MPI_integer, MPI_sum, 0, &
          & MPI_COMM_WORLD, ierr)

     call MPI_reduce(dmax_beta_max_bin, dmax_beta_max_bin_result, &
          & size(dmax_beta_max_bin), MPI_integer, MPI_sum, 0, &
          & MPI_COMM_WORLD, ierr)

     call MPI_reduce(dmax_n_joint_max, dmax_n_joint_max_result, &
          & size(dmax_n_joint_max), MPI_integer, MPI_sum, 0, &
          & MPI_COMM_WORLD, ierr)

  end if

  if (use_d_diff_statistics) then
     call MPI_reduce( &
          & d_diff_max_angle, &
          & d_diff_max_result, &
          & size(d_diff_max_angle), &
          & MPI_double_precision, MPI_sum, &
          & 0, MPI_COMM_WORLD, ierr)

     call MPI_reduce( &
          & d_diff_max_circle1_pixel, &
          & d_diff_max_circle1_pixel_result, &
          & size(d_diff_max_circle1_pixel), &
          & MPI_integer, MPI_sum, &
          & 0, MPI_COMM_WORLD, ierr)

     call MPI_reduce( &
          & d_diff_max_circle2_pixel, &
          & d_diff_max_circle2_pixel_result, &
          & size(d_diff_max_circle2_pixel), &
          & MPI_integer, MPI_sum, &
          & 0, MPI_COMM_WORLD, ierr)

     call MPI_reduce( &
          & d_diff_max_beta_bin, &
          & d_diff_max_beta_bin_result, &
          & size(d_diff_max_beta_bin), &
          & MPI_integer, MPI_sum, &
          & 0, MPI_COMM_WORLD, ierr)

     call MPI_reduce( &
          & d_diff_max_n_joint, &
          & d_diff_max_n_joint_result, &
          & size(d_diff_max_n_joint), &
          & MPI_integer, MPI_sum, &
          & 0, MPI_COMM_WORLD, ierr)

  end if

  if (mpiid .eq. 0) then

     if (use_s_statistics) then
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

     end if

     if (use_d_statistics) then

        dmax_circle1_pixel_result = dmax_circle1_pixel_result - 1
        dmax_circle2_pixel_result = dmax_circle2_pixel_result - 1

        open(unit=d_full_outfile, &
             & file=trim(adjustl(inifilename(1:len_trim(inifilename)-4)))// &
             & "_d.dat")

        do i = nside_search/2, 2*nside_search

           if ((rad_to_deg(dble(theta_rings(i))) .lt. &
                & opening_angle_range%min) .or. &
               (rad_to_deg(dble(theta_rings(i))) .gt. &
                & opening_angle_range%max)) cycle

           call write_d_pair_info(d_full_outfile, theta_rings(i), &
                & dmax_result(i), dmax_circle1_pixel_result(i), &
                & dmax_circle2_pixel_result(i), &
                & dmax_beta_max_bin_result(i), &
                & dmax_n_joint_max_result(i), nside_search)
        end do

        close(d_full_outfile)
     end if

     if (use_d_diff_statistics) then

        d_diff_max_circle1_pixel_result = &
             & d_diff_max_circle1_pixel_result - 1

        d_diff_max_circle2_pixel_result = &
             & d_diff_max_circle2_pixel_result - 1

     if (use_d_diff_q_statistics) then
        open(unit=d_full_outfile, &
             & file=trim(adjustl( &
             & inifilename(1:len_trim(inifilename)-4))) // &
             & "_q_diff.dat")
     else
        open(unit=d_full_outfile, &
             & file=trim(adjustl( &
             & inifilename(1:len_trim(inifilename)-4))) // &
             & "_d_diff.dat")
     end if

        do i = nside_search/2, 2*nside_search

           if ((rad_to_deg(dble(theta_rings(i))) .lt. &
                & opening_angle_range%min) .or. &
               (rad_to_deg(dble(theta_rings(i))) .gt. &
                & opening_angle_range%max)) cycle

           call write_d_pair_info( &
                & d_full_outfile, &
                & theta_rings(i), &
                & d_diff_max_result(i), &
                & d_diff_max_circle1_pixel_result(i), &
                & d_diff_max_circle2_pixel_result(i), &
                & d_diff_max_beta_bin_result(i), &
                & d_diff_max_n_joint_result(i), &
                & nside_search)
        end do
        close(d_full_outfile)
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

  subroutine write_d_pair_info(file, opening_angle, dmax, pixel1, pixel2, &
       & beta_max_bin, n_joint_max, nside)

    implicit none

    integer,  intent(in) :: file
    real,     intent(in) :: opening_angle
    real(dp), intent(in) :: dmax
    integer,  intent(in) :: pixel1, pixel2
    integer,  intent(in) :: beta_max_bin, n_joint_max, nside

    if (pixel1 < 0 .or. pixel2 < 0 .or. beta_max_bin < 0) then
       write(file, *) rad_to_deg(dble(opening_angle)), dmax, &
            & pixel1, 0.0_dp, 0.0_dp, pixel2, 0.0_dp, 0.0_dp, &
            & 0.0_dp, beta_max_bin, 0.0_dp, n_joint_max
       return
    end if

    write(file, *) rad_to_deg(dble(opening_angle)), &
         & dmax, &
         & pixel1, &
         & rad_to_deg(pix2ang_table(pixel1)%theta), &
         & rad_to_deg(pix2ang_table(pixel1)%phi), &
         & pixel2, &
         & rad_to_deg(pix2ang_table(pixel2)%theta), &
         & rad_to_deg(pix2ang_table(pixel2)%phi), &
         & rad_to_deg(angle_between_pixels(nside, pixel1, pixel2)), &
         & beta_max_bin, &
         & 360.0_dp * dble(beta_max_bin) / dble(nphi_ring), &
         & n_joint_max

  end subroutine write_d_pair_info

  real(dp) function angle_between_pixels(nside,pixel1, pixel2)
    implicit none
    integer, intent(in) :: nside
    integer, intent(in) :: pixel1, pixel2

    real(dp),  dimension(3) :: x, y
    call pix2vec_ring(nside, pixel1, x)
    call pix2vec_ring(nside, pixel2, y)
    angle_between_pixels = acos(max( &
         & -1.0_dp, min(1.0_dp, &
         & dot_product(x,y) / &
         & sqrt(dot_product(x,x)*dot_product(y,y)))))

  end function angle_between_pixels

  subroutine finalize()
    implicit none

#ifdef DFFT
    call dfftw_destroy_plan(plan_forward)
    call dfftw_destroy_plan(plan_backward)    
#else
    call sfftw_destroy_plan(plan_forward)
    call sfftw_destroy_plan(plan_backward)
#endif
    deallocate(temperature_in)
    deallocate(obs_mask)
    deallocate(temperature_out)
    deallocate(sm_out)
    deallocate(sm_in)
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

    if (allocated(d_w_t_transforms)) then
       deallocate(d_w_t_transforms)
    end if

    if (allocated(d_w_transforms)) then
       deallocate(d_w_transforms)
    end if

    if (allocated(d_w_t_sq_transforms)) then
       deallocate(d_w_t_sq_transforms)
    end if

    if (allocated(d_diff_w_delta_t_transforms)) then
       deallocate(d_diff_w_delta_t_transforms)
    end if

    if (allocated(d_diff_w_transforms)) then
       deallocate(d_diff_w_transforms)
    end if

    if (allocated(d_diff_w_delta_t_sq_transforms)) then
       deallocate(d_diff_w_delta_t_sq_transforms)
    end if

    if (allocated(d_diff_null_hist_counts)) then
       deallocate(d_diff_null_hist_counts)
    end if

    if (allocated(d_map)) deallocate(d_map)
    if (allocated(d_filter)) deallocate(d_filter)
    if (allocated(dmax_angle)) deallocate(dmax_angle)
    if (allocated(dmax_result)) deallocate(dmax_result)
    if (allocated(dmax_circle1_pixel)) deallocate(dmax_circle1_pixel)
    if (allocated(dmax_circle2_pixel)) deallocate(dmax_circle2_pixel)
    if (allocated(dmax_circle1_pixel_result)) then
       deallocate(dmax_circle1_pixel_result)
    end if
    if (allocated(dmax_circle2_pixel_result)) then
       deallocate(dmax_circle2_pixel_result)
    end if
    if (allocated(dmax_beta_max_bin)) deallocate(dmax_beta_max_bin)
    if (allocated(dmax_beta_max_bin_result)) then
       deallocate(dmax_beta_max_bin_result)
    end if
    if (allocated(dmax_n_joint_max)) deallocate(dmax_n_joint_max)
    if (allocated(dmax_n_joint_max_result)) then
       deallocate(dmax_n_joint_max_result)
    end if
    
    if (allocated(d_diff_max_angle)) then
       deallocate(d_diff_max_angle)
    end if
    
    if (allocated(d_diff_max_result)) then
       deallocate(d_diff_max_result)
    end if
    
    if (allocated(d_diff_max_circle1_pixel)) then
       deallocate(d_diff_max_circle1_pixel)
    end if
    
    if (allocated(d_diff_max_circle2_pixel)) then
       deallocate(d_diff_max_circle2_pixel)
    end if
    
    if (allocated(d_diff_max_circle1_pixel_result)) then
       deallocate(d_diff_max_circle1_pixel_result)
    end if
    
    if (allocated(d_diff_max_circle2_pixel_result)) then
       deallocate(d_diff_max_circle2_pixel_result)
    end if
    
    if (allocated(d_diff_max_beta_bin)) then
       deallocate(d_diff_max_beta_bin)
    end if
    
    if (allocated(d_diff_max_beta_bin_result)) then
       deallocate(d_diff_max_beta_bin_result)
    end if
    
    if (allocated(d_diff_max_n_joint)) then
       deallocate(d_diff_max_n_joint)
    end if
    
    if (allocated(d_diff_max_n_joint_result)) then
       deallocate(d_diff_max_n_joint_result)
    end if

    if (allocated(d_diff_gaussian_mu)) then
       deallocate(d_diff_gaussian_mu)
    end if
    
    if (allocated(d_diff_gaussian_sigma)) then
       deallocate(d_diff_gaussian_sigma)
    end if
    
    if (allocated(d_diff_gaussian_valid)) then
       deallocate(d_diff_gaussian_valid)
    end if
    
    if (allocated(d_diff_gaussian_fit_method)) then
       deallocate(d_diff_gaussian_fit_method)
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

    use_s_statistics = ini_read_logical( &
         & "use_s_statistics", .true.)
    
    use_d_statistics = ini_read_logical( &
         & "use_d_statistics", .false.)
    
    use_d_diff_statistics = ini_read_logical( &
         & "use_d_diff_statistics", .false.)

    d_diff_max_separation = ini_read_int( &
         & "d_diff_max_separation", 2)

    use_d_diff_bias_correction = ini_read_logical( &
         & "use_d_diff_bias_correction", .false.)

    use_d_diff_q_statistics = ini_read_logical( &
         & "use_d_diff_q_statistics", .false.)

    if (use_d_diff_bias_correction .and. &
        & .not. use_d_diff_statistics) then
    
       write(*,*) &
            "ERROR: use_d_diff_bias_correction requires ", &
            "use_d_diff_statistics = T"
       stop
    end if

    if (use_d_diff_q_statistics .and. &
        & .not. use_d_diff_statistics) then
    
       write(*,*) &
            "ERROR: use_d_diff_q_statistics requires ", &
            "use_d_diff_statistics = T"
       stop
    end if

    write_d_diff_null_histogram = ini_read_logical( &
         & "write_d_diff_null_histogram", .false.)
    
    d_diff_hist_nbins = ini_read_int( &
         & "d_diff_hist_nbins", 2000)
    
    if (write_d_diff_null_histogram .and. &
        & .not. use_d_diff_statistics) then
    
       write(*,*) &
            "ERROR: write_d_diff_null_histogram requires ", &
            "use_d_diff_statistics = T"
       stop
    end if
    
    if (d_diff_hist_nbins < 1) then
       write(*,*) "ERROR: d_diff_hist_nbins must be positive"
       stop
    end if

    maskfilename = ""
    d_filter_filename = ""
    d_diff_gaussian_lookup_filename = ""

    if (use_d_statistics .or. use_d_diff_statistics) then
       maskfilename = ini_read_string("maskfilename")
       d_filter_filename = ini_read_string("d_filter_filename")
       d_filter_lmax = ini_read_int("d_filter_lmax", 300)
       mask_threshold = ini_read_real("mask_threshold", 0.5)
       d_min_joint_samples = ini_read_int("d_min_joint_samples")

       if (use_d_diff_q_statistics) then
          d_diff_gaussian_lookup_filename = ini_read_string( &
               & "d_diff_gaussian_lookup_filename")
       end if
    end if

    if (m_min < 1 .or. m_min > nphi_ring/2) then
       write(*,*) "ERROR: m_min must be between 1 and nphi_ring/2"
       stop
    end if

    if (use_d_diff_statistics) then
       if (d_diff_max_separation < 1 .or. &
           & d_diff_max_separation > nphi_ring/2) then
          write(*,*) &
               "ERROR: d_diff_max_separation must be between 1 and nphi_ring/2"
          stop
       end if
    end if

    if (use_d_statistics .or. use_d_diff_statistics) then
       if (d_min_joint_samples < 1 .or. &
           d_min_joint_samples > nphi_ring) then
          write(*,*) "ERROR: d_min_joint_samples must be between 1 and nphi_ring"
          stop
       end if

       if (d_filter_lmax < 0 .or. d_filter_lmax > 3*nside - 1) then
          write(*,*) "ERROR: d_filter_lmax must be between 0 and ", &
               & 3*nside - 1
          stop
       end if
    end if

    orientable = ini_read_logical("orientable", .false.)
    nside_search = ini_read_int("nside_search", 128)
    min_separation_angle = 180.-ini_read_real("min_separation_angle", 179.5)
    mapfilename=ini_read_string("mapfilename")
    write_every_nth_smax=ini_read_logical("write_every_nth_smax", .false.)
    every_nth_smax=ini_read_int("every_nth_smax", 10000)
    write_position_smax=ini_read_logical("write_position_smax", .false.)
    random_s = ini_read_logical("random_s", .false.)
    scramble_sky = ini_read_logical("scramble_sky", .false.)
    run_symmetry_test = ini_read_logical("run_symmetry_test", .false.)
    symmetry_input_has_missing_pixels = ini_read_logical( &
     & "symmetry_input_has_missing_pixels", .true.)

    if (run_symmetry_test) then
       symmetry_center_theta_deg = ini_read_real("symmetry_center_theta_deg")
       symmetry_center_phi_deg   = ini_read_real("symmetry_center_phi_deg")
    end if
    
    call ini_close()

    npix = 12*nside**2
    npix_search = 12*nside_search**2
    allocate(map(0:npix-1))
    allocate(obs_mask(0:npix-1))
    obs_mask = 1.0

    if (use_d_statistics .or. use_d_diff_statistics) then
       allocate(d_map(0:npix-1))
       allocate(d_filter(0:d_filter_lmax))
    
       d_map = 0.0
       d_filter = 0.0_dp
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

    if (use_d_statistics) then
       allocate(dmax_angle(nside_search/2:2*nside_search))
       allocate(dmax_result(nside_search/2:2*nside_search))
       allocate(dmax_circle1_pixel(nside_search/2:2*nside_search))
       allocate(dmax_circle2_pixel(nside_search/2:2*nside_search))
       allocate(dmax_circle1_pixel_result(nside_search/2:2*nside_search))
       allocate(dmax_circle2_pixel_result(nside_search/2:2*nside_search))
       allocate(dmax_beta_max_bin(nside_search/2:2*nside_search))
       allocate(dmax_beta_max_bin_result(nside_search/2:2*nside_search))
       allocate(dmax_n_joint_max(nside_search/2:2*nside_search))
       allocate(dmax_n_joint_max_result(nside_search/2:2*nside_search))
    end if

    if (use_d_diff_statistics) then
       allocate(d_diff_max_angle(nside_search/2:2*nside_search))
       allocate(d_diff_max_result(nside_search/2:2*nside_search))
    
       allocate(d_diff_max_circle1_pixel( &
            nside_search/2:2*nside_search))
       allocate(d_diff_max_circle2_pixel( &
            nside_search/2:2*nside_search))
    
       allocate(d_diff_max_circle1_pixel_result( &
            nside_search/2:2*nside_search))
       allocate(d_diff_max_circle2_pixel_result( &
            nside_search/2:2*nside_search))
    
       allocate(d_diff_max_beta_bin( &
            nside_search/2:2*nside_search))
       allocate(d_diff_max_beta_bin_result( &
            nside_search/2:2*nside_search))
    
       allocate(d_diff_max_n_joint( &
            nside_search/2:2*nside_search))
       allocate(d_diff_max_n_joint_result( &
            nside_search/2:2*nside_search))
    end if

    allocate(circle_center_vectors(3,0:npix_search-1))
    allocate(temperature_out(nphi_ring/2+1))
    allocate(sm_in(nphi_ring/2+1))
    allocate(temperature_in(nphi_ring))    
    allocate(sm_out(nphi_ring))
    if (use_s_statistics) then
       allocate(circle_transforms( &
            & 1:nphi_ring/2+1, 0:npix_search-1))
    end if

    if (use_d_statistics .and. .not. use_s_statistics) then
       allocate(d_w_t_transforms( &
            & 1:nphi_ring/2+1, 0:npix_search-1))

       allocate(d_w_transforms( &
            & 1:nphi_ring/2+1, 0:npix_search-1))

       allocate(d_w_t_sq_transforms( &
            & 1:nphi_ring/2+1, 0:npix_search-1))
    end if

    if (use_d_diff_statistics) then
       allocate(d_diff_w_delta_t_transforms( &
            & 1:nphi_ring/2+1, &
            & 0:npix_search-1, &
            & 1:d_diff_max_separation))
    
       allocate(d_diff_w_transforms( &
            & 1:nphi_ring/2+1, &
            & 0:npix_search-1, &
            & 1:d_diff_max_separation))
    
       allocate(d_diff_w_delta_t_sq_transforms( &
            & 1:nphi_ring/2+1, &
            & 0:npix_search-1, &
            & 1:d_diff_max_separation))
    end if

    if (write_d_diff_null_histogram) then
    
       allocate(d_diff_null_hist_counts( &
            & 0:nphi_ring, &
            & 0:d_diff_hist_nbins-1))
    
       d_diff_null_hist_counts = 0_8
    
    end if

#ifdef DFFT
    call dfftw_plan_dft_r2c_1d(plan_forward, nphi_ring, temperature_in, temperature_out, FFTW_ESTIMATE)
    call dfftw_plan_dft_c2r_1d(plan_backward, nphi_ring, sm_in, sm_out, FFTW_ESTIMATE)
#else
    call sfftw_plan_dft_r2c_1d(plan_forward, nphi_ring, temperature_in, temperature_out, FFTW_ESTIMATE)
    call sfftw_plan_dft_c2r_1d(plan_backward, nphi_ring, sm_in, sm_out, FFTW_ESTIMATE)
#endif

    if (trim(adjustl(mapfilename)) .eq. "simsky.dat") then
       call read_simulated_sky(nside, scramble_sky, map)
    else
!      call read_map(nside, trim(adjustl(mapfilename)), map)
       call read_map_fits(nside, trim(adjustl(mapfilename)), map)
    end if

    if (use_d_statistics .or. use_d_diff_statistics) then
       call read_d_filter(d_filter_filename, d_filter_lmax, d_filter)
       call apply_d_filter( &
            & nside, d_filter_lmax, map, d_filter, d_map)
       call read_map_fits(nside, trim(adjustl(maskfilename)), obs_mask)
    end if

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
  
  subroutine get_d_at_beta(circle1, circle2, beta_bin, dvalue, n_joint)
  
    implicit none
  
    type(tcircle), intent(in)  :: circle1, circle2
    integer,       intent(in)  :: beta_bin
    real(dp),      intent(out) :: dvalue
    integer,       intent(out) :: n_joint
  
    type(tcircle) :: d_circle1, d_circle2
  
    integer :: nphi
    integer :: iphi
    integer :: partner
  
    real(dp) :: numerator
    real(dp) :: denominator
    real(dp) :: d_diff_try
  
    real(FP) :: ring_i(circle1%nphi_ring)
    real(FP) :: ring_j(circle1%nphi_ring)
  
    FFT_COMPLEX :: modes_i(circle1%nphi_ring/2 + 1)
    FFT_COMPLEX :: modes_j(circle1%nphi_ring/2 + 1)
  
    logical :: valid_i(circle1%nphi_ring)
    logical :: valid_j(circle1%nphi_ring)
  
    nphi = circle1%nphi_ring
  
    d_circle1 = circle1
    d_circle2 = circle2
  
    d_circle1%theta = pix2ang_table(d_circle1%pixel)%theta
    d_circle1%phi   = pix2ang_table(d_circle1%pixel)%phi
  
    d_circle2%theta = pix2ang_table(d_circle2%pixel)%theta
    d_circle2%phi   = pix2ang_table(d_circle2%pixel)%phi
  
    call get_temperature_along_circle(d_circle1, nside_map, d_map, &
         & use_pixel_average, nest2ring_table, ring_i, obs_mask, &
         & mask_threshold, valid_i)
  
    call get_temperature_along_circle(d_circle2, nside_map, d_map, &
         & use_pixel_average, nest2ring_table, ring_j, obs_mask, &
         & mask_threshold, valid_j)
  
#ifdef DFFT
    call dfftw_execute_dft_r2c(plan_forward, ring_i, modes_i)
    call dfftw_execute_dft_r2c(plan_forward, ring_j, modes_j)
#else
    call sfftw_execute_dft_r2c(plan_forward, ring_i, modes_i)
    call sfftw_execute_dft_r2c(plan_forward, ring_j, modes_j)
#endif
  
    modes_i(1:m_min) = 0.0
    modes_j(1:m_min) = 0.0
  
#ifdef DFFT
    call dfftw_execute_dft_c2r(plan_backward, modes_i, ring_i)
    call dfftw_execute_dft_c2r(plan_backward, modes_j, ring_j)
#else
    call sfftw_execute_dft_c2r(plan_backward, modes_i, ring_i)
    call sfftw_execute_dft_c2r(plan_backward, modes_j, ring_j)
#endif
  
    ring_i = ring_i / dble(nphi)
    ring_j = ring_j / dble(nphi)
  
    numerator = 0.0_dp
    denominator = 0.0_dp
    n_joint = 0
  
    do iphi = 1, nphi
  
       if (orientable) then
          partner = modulo(beta_bin - (iphi - 1), nphi) + 1
       else
          partner = modulo((iphi - 1) + beta_bin, nphi) + 1
       end if
  
       if (.not. (valid_i(iphi) .and. valid_j(partner))) cycle
  
       numerator = numerator + 2.0_dp * &
            & real(ring_i(iphi), dp) * real(ring_j(partner), dp)
  
       denominator = denominator + real(ring_i(iphi), dp)**2 + &
            & real(ring_j(partner), dp)**2
  
       n_joint = n_joint + 1
    end do
  
    if (n_joint < d_min_joint_samples .or. &
        denominator <= tiny(1.0_dp)) then
       dvalue = 0.0_dp
       return
    end if
  
    dvalue = numerator / denominator
    dvalue = max(-1.0_dp, min(1.0_dp, dvalue))
  end subroutine get_d_at_beta

  subroutine get_d_max(circle1, circle2, d_max, &
       & beta_max_bin, n_joint_max)

    implicit none

    type(tcircle), intent(in)  :: circle1, circle2
    real(dp),      intent(out) :: d_max
    integer,       intent(out) :: beta_max_bin
    integer,       intent(out) :: n_joint_max

    integer :: nphi
    integer :: beta_bin
    integer :: beta_index
    integer :: n_joint

    real(dp) :: numerator
    real(dp) :: denominator
    real(dp) :: d_try

    FFT_COMPLEX :: numerator_modes( &
         & circle1%nphi_ring/2 + 1)

    FFT_COMPLEX :: denominator_modes( &
         & circle1%nphi_ring/2 + 1)

    FFT_COMPLEX :: n_joint_modes( &
         & circle1%nphi_ring/2 + 1)

    real(FP) :: numerator_by_beta(circle1%nphi_ring)
    real(FP) :: denominator_by_beta(circle1%nphi_ring)
    real(FP) :: n_joint_by_beta(circle1%nphi_ring)

    logical :: found_valid_beta

    nphi = circle1%nphi_ring


    ! ==================================================
    ! Construct numerator, denominator and N_joint
    ! for every beta at once.
    ! ==================================================

    if (orientable) then

       ! partner = beta - k
       ! Circular convolution.

       numerator_modes = &
            d_w_t_transforms(:,circle1%pixel) * &
            d_w_t_transforms(:,circle2%pixel)

       denominator_modes = &
            d_w_t_sq_transforms(:,circle1%pixel) * &
            d_w_transforms(:,circle2%pixel) + &
            d_w_transforms(:,circle1%pixel) * &
            d_w_t_sq_transforms(:,circle2%pixel)

       n_joint_modes = &
            d_w_transforms(:,circle1%pixel) * &
            d_w_transforms(:,circle2%pixel)
    else
       ! partner = k + beta
       ! Circular correlation.

       numerator_modes = &
            conjg(d_w_t_transforms(:,circle1%pixel)) * &
            d_w_t_transforms(:,circle2%pixel)

       denominator_modes = &
            conjg(d_w_t_sq_transforms(:,circle1%pixel)) * &
            d_w_transforms(:,circle2%pixel) + &
            conjg(d_w_transforms(:,circle1%pixel)) * &
            d_w_t_sq_transforms(:,circle2%pixel)

       n_joint_modes = &
            conjg(d_w_transforms(:,circle1%pixel)) * &
            d_w_transforms(:,circle2%pixel)
    end if

#ifdef DFFT
    call dfftw_execute_dft_c2r( &
         & plan_backward, numerator_modes, &
         & numerator_by_beta)

    call dfftw_execute_dft_c2r( &
         & plan_backward, denominator_modes, &
         & denominator_by_beta)

    call dfftw_execute_dft_c2r( &
         & plan_backward, n_joint_modes, &
         & n_joint_by_beta)
#else
    call sfftw_execute_dft_c2r( &
         & plan_backward, numerator_modes, &
         & numerator_by_beta)

    call sfftw_execute_dft_c2r( &
         & plan_backward, denominator_modes, &
         & denominator_by_beta)

    call sfftw_execute_dft_c2r( &
         & plan_backward, n_joint_modes, &
         & n_joint_by_beta)
#endif
    d_max = -huge(1.0_dp)
    beta_max_bin = 0
    n_joint_max = 0
    found_valid_beta = .false.


    do beta_bin = 0, nphi - 1

       beta_index = beta_bin + 1

       denominator = real( &
            denominator_by_beta(beta_index), dp)

       n_joint = nint( &
            real(n_joint_by_beta(beta_index), dp) / &
            dble(nphi))

       n_joint = max(0, min(nphi, n_joint))

       if (n_joint < d_min_joint_samples .or. &
           denominator <= tiny(1.0_dp)) cycle

       numerator = 2.0_dp * real( &
            numerator_by_beta(beta_index), dp)

       d_try = numerator / denominator
       d_try = max(-1.0_dp, min(1.0_dp, d_try))

       if (.not. found_valid_beta .or. &
           d_try > d_max) then

          d_max = d_try
          beta_max_bin = beta_bin
          n_joint_max = n_joint
          found_valid_beta = .true.
       end if

    end do

    if (.not. found_valid_beta) then
       d_max = 0.0_dp
       beta_max_bin = 0
       n_joint_max = 0
    end if

  end subroutine get_d_max

  real(dp) function d_diff_bias_mean(alpha_deg, n_joint)
  
    implicit none
  
    real(dp), intent(in) :: alpha_deg
    integer,  intent(in) :: n_joint
  
    real(dp) :: x
    real(dp) :: A_alpha
    real(dp) :: gamma_alpha
  
    x = (alpha_deg - 55.0_dp) / 30.0_dp
  
    A_alpha = &
         & -0.0854624_dp * x**3 &
         & -0.15252581_dp * x**2 &
         & +0.49238932_dp * x &
         & +1.61226065_dp
  
    gamma_alpha = &
         & -0.00348046_dp * x**3 &
         & -0.03711699_dp * x**2 &
         & +0.07374485_dp * x &
         & +0.36047432_dp
  
    d_diff_bias_mean = &
         & A_alpha * dble(n_joint)**(-gamma_alpha)
  
  end function d_diff_bias_mean

  subroutine get_d_diff_max(circle1, circle2, d_diff_max, &
       & beta_max_bin, n_joint_max)

    implicit none

    type(tcircle), intent(in)  :: circle1, circle2
    real(dp),      intent(out) :: d_diff_max
    integer,       intent(out) :: beta_max_bin
    integer,       intent(out) :: n_joint_max

    integer :: nphi
    integer :: beta_bin
    integer :: beta_array_index
    integer :: separation_bin
    integer :: n_joint
    integer :: n_joint_separation

    real(dp) :: numerator
    real(dp) :: denominator
    real(dp) :: d_diff_try
    real(dp) :: score_try

    FFT_COMPLEX :: numerator_modes( &
         & circle1%nphi_ring/2 + 1)

    FFT_COMPLEX :: denominator_modes( &
         & circle1%nphi_ring/2 + 1)

    FFT_COMPLEX :: n_joint_modes( &
         & circle1%nphi_ring/2 + 1)

    real(FP) :: numerator_by_beta( &
         & circle1%nphi_ring, d_diff_max_separation)

    real(FP) :: denominator_by_beta( &
         & circle1%nphi_ring, d_diff_max_separation)

    real(FP) :: n_joint_by_beta( &
         & circle1%nphi_ring, d_diff_max_separation)

    logical :: found_valid_beta

    nphi = circle1%nphi_ring

    ! ==================================================
    ! Calculate numerator, denominator and N_joint
    ! for every beta and every separation.
    ! ==================================================
    do separation_bin = 1, d_diff_max_separation

       if (orientable) then
          ! Circular convolution.
          numerator_modes = &
               d_diff_w_delta_t_transforms( &
               :, circle1%pixel, separation_bin) * &
               d_diff_w_delta_t_transforms( &
               :, circle2%pixel, separation_bin)

          denominator_modes = &
               d_diff_w_delta_t_sq_transforms( &
               :, circle1%pixel, separation_bin) * &
               d_diff_w_transforms( &
               :, circle2%pixel, separation_bin) + &
               d_diff_w_transforms( &
               :, circle1%pixel, separation_bin) * &
               d_diff_w_delta_t_sq_transforms( &
               :, circle2%pixel, separation_bin)

          n_joint_modes = &
               d_diff_w_transforms( &
               :, circle1%pixel, separation_bin) * &
               d_diff_w_transforms( &
               :, circle2%pixel, separation_bin)
       else
          ! Circular correlation.
          numerator_modes = &
               conjg(d_diff_w_delta_t_transforms( &
               :, circle1%pixel, separation_bin)) * &
               d_diff_w_delta_t_transforms( &
               :, circle2%pixel, separation_bin)

          denominator_modes = &
               conjg(d_diff_w_delta_t_sq_transforms( &
               :, circle1%pixel, separation_bin)) * &
               d_diff_w_transforms( &
               :, circle2%pixel, separation_bin) + &
               conjg(d_diff_w_transforms( &
               :, circle1%pixel, separation_bin)) * &
               d_diff_w_delta_t_sq_transforms( &
               :, circle2%pixel, separation_bin)

          n_joint_modes = &
               conjg(d_diff_w_transforms( &
               :, circle1%pixel, separation_bin)) * &
               d_diff_w_transforms( &
               :, circle2%pixel, separation_bin)
       end if


       ! Inverse FFTs give all beta values for this
       ! separation.
#ifdef DFFT
       call dfftw_execute_dft_c2r( &
            & plan_backward, numerator_modes, &
            & numerator_by_beta(:,separation_bin))

       call dfftw_execute_dft_c2r( &
            & plan_backward, denominator_modes, &
            & denominator_by_beta(:,separation_bin))

       call dfftw_execute_dft_c2r( &
            & plan_backward, n_joint_modes, &
            & n_joint_by_beta(:,separation_bin))
#else
       call sfftw_execute_dft_c2r( &
            & plan_backward, numerator_modes, &
            & numerator_by_beta(:,separation_bin))

       call sfftw_execute_dft_c2r( &
            & plan_backward, denominator_modes, &
            & denominator_by_beta(:,separation_bin))

       call sfftw_execute_dft_c2r( &
            & plan_backward, n_joint_modes, &
            & n_joint_by_beta(:,separation_bin))
#endif

    end do


    ! ==================================================
    ! Search for the best beta.
    ! ==================================================
    d_diff_max = -huge(1.0_dp)
    beta_max_bin = 0
    n_joint_max = 0
    found_valid_beta = .false.

    do beta_bin = 0, nphi - 1

       numerator = 0.0_dp
       denominator = 0.0_dp

       ! This will become the smallest N_joint among
       ! all requested separations.
       n_joint = nphi

       do separation_bin = 1, d_diff_max_separation

          if (orientable) then
             ! For separation s, the cached forward
             ! difference on circle j must be read at
             ! beta-s and multiplied by -1.
             beta_array_index = modulo( &
                  & beta_bin - separation_bin, nphi) + 1

             numerator = numerator - 2.0_dp * real( &
                  numerator_by_beta( &
                  beta_array_index, separation_bin), dp)

          else

             beta_array_index = beta_bin + 1

             numerator = numerator + 2.0_dp * real( &
                  numerator_by_beta( &
                  beta_array_index, separation_bin), dp)

          end if

          denominator = denominator + real( &
               denominator_by_beta( &
               beta_array_index, separation_bin), dp)

          ! FFTW's backward transform contains an
          ! N_phi factor. Remove it from the count.
          n_joint_separation = nint( &
               real(n_joint_by_beta( &
               beta_array_index, separation_bin), dp) / &
               dble(nphi))

          n_joint_separation = max( &
               0, min(nphi, n_joint_separation))

          ! Require all separations to have adequate
          ! jointly valid support.
          n_joint = min( &
               n_joint, n_joint_separation)

       end do

       if (n_joint < d_min_joint_samples .or. &
           & denominator <= tiny(1.0_dp)) cycle

       ! Raw numerator and denominator contributions
       ! from all requested separations have already
       ! been summed.
       d_diff_try = numerator / denominator

       d_diff_try = max( &
            -1.0_dp, min(1.0_dp, d_diff_try))

       if (use_d_diff_q_statistics) then
       
          score_try = -log10( &
               0.5_dp * erfc( &
               (d_diff_try - &
                d_diff_gaussian_mu(iopening_angle, n_joint)) / &
               (sqrt(2.0_dp) * &
                d_diff_gaussian_sigma(iopening_angle, n_joint)) &
               ) &
               )
       else
          score_try = d_diff_try
       end if

       if (.not. found_valid_beta .or. &
           & score_try > d_diff_max) then
          d_diff_max = score_try
          beta_max_bin = beta_bin
          n_joint_max = n_joint
          found_valid_beta = .true.
       end if

    end do

    if (.not. found_valid_beta) then
       d_diff_max = 0.0_dp
       beta_max_bin = 0
       n_joint_max = 0
    end if

  end subroutine get_d_diff_max

  subroutine load_d_diff_gaussian_lookup()
  
    implicit none
  
    integer :: lookup_unit
    integer :: ios
  
    integer :: alpha_index
    integer :: n_joint
    integer :: n_geometry
    integer :: n_pooled_values
    integer :: fit_valid
    integer :: fit_method
  
    integer :: first_iopening_angle
    integer :: n_active_alpha
    integer :: iopening_lookup
    integer :: ith
    integer :: n
  
    real(dp) :: alpha_deg
    real(dp) :: mu_value
    real(dp) :: sigma_value
    
    character(len=1000) :: line
  
    ! ==================================================
    ! This Gaussian table was calibrated specifically for:
    !   nphi = 512
    !   N_joint = 100,...,512
    ! ==================================================
  
    if (nphi_ring /= 512) then
       write(*,*) &
            "ERROR: D_diff Gaussian lookup requires nphi_ring = 512"
       stop
    end if
  
    if (d_min_joint_samples < 100 .or. &
        & d_min_joint_samples > 512) then
       write(*,*) &
            "ERROR: D_diff Gaussian lookup requires ", &
            "100 <= d_min_joint_samples <= 512"
       stop
    end if

    if (.not. orientable) then
       write(*,*) &
            "ERROR: D_diff Gaussian lookup requires orientable = T"
       stop
    end if
    
    if (d_diff_max_separation /= 2) then
       write(*,*) &
            "ERROR: D_diff Gaussian lookup requires ", &
            "d_diff_max_separation = 2"
       stop
    end if
    
    ! ==================================================
    ! Find the first active opening-angle grid point.
    !
    ! For the current 25--45 degree search this gives
    ! 52 consecutive theta_rings.
    !
    ! Table alpha_index = 0 corresponds to the first one,
    ! alpha_index = 1 to the next one, etc.
    ! ==================================================
  
    first_iopening_angle = -1
    n_active_alpha = 0
  
    do ith = nside_search/2, 2*nside_search
  
       if ( &
            rad_to_deg(dble(theta_rings(ith))) >= &
            dble(opening_angle_range%min) .and. &
            rad_to_deg(dble(theta_rings(ith))) <= &
            dble(opening_angle_range%max) &
          ) then
  
          if (first_iopening_angle < 0) then
             first_iopening_angle = ith
          end if
  
          n_active_alpha = n_active_alpha + 1
       end if
  
    end do
  
    if (n_active_alpha /= 52) then
       write(*,*) &
            "ERROR: Gaussian lookup contains 52 alpha values, ", &
            "but current search contains ", n_active_alpha
       stop
    end if
  
    ! ==================================================
    ! Allocate lookup:
    !
    ! first index  = iopening_angle
    ! second index = N_joint
    ! ==================================================
  
    allocate( &
         d_diff_gaussian_mu( &
         nside_search/2:2*nside_search, &
         0:nphi_ring) &
    )
  
    allocate( &
         d_diff_gaussian_sigma( &
         nside_search/2:2*nside_search, &
         0:nphi_ring) &
    )
  
    d_diff_gaussian_mu = 0.0_dp
    d_diff_gaussian_sigma = -1.0_dp
  
    ! ==================================================
    ! Open lookup table.
    ! Same basic I/O pattern as read_d_filter().
    ! ==================================================
  
    open( &
         newunit=lookup_unit, &
         file=trim(adjustl( &
         d_diff_gaussian_lookup_filename)), &
         status="old", &
         action="read", &
         iostat=ios &
    )
  
    if (ios /= 0) then
       write(*,*) &
            "ERROR: could not open D_diff Gaussian lookup:"
       write(*,*) &
            trim(adjustl( &
            d_diff_gaussian_lookup_filename))
       stop
    end if
  
    ! ==================================================
    ! File columns:
    !
    ! 1 alpha_index
    ! 2 alpha_deg
    ! 3 N_joint
    ! 4 mu
    ! 5 sigma
    ! 6 n_geometry
    ! 7 n_pooled_values
    ! 8 fit_valid
    ! 9 fit_method
    ! ==================================================
    do
       read(lookup_unit, "(A)", iostat=ios) line
  
       ! End of file.
       if (ios < 0) exit
  
       if (ios > 0) then
          write(*,*) &
               "ERROR while reading D_diff Gaussian lookup"
          stop
       end if
  
       line = adjustl(line)
  
       ! Ignore blank/comment lines.
       if (len_trim(line) == 0) cycle
       if (line(1:1) == "#") cycle
  
       read( &
            line, *, iostat=ios &
       ) &
            alpha_index, &
            alpha_deg, &
            n_joint, &
            mu_value, &
            sigma_value, &
            n_geometry, &
            n_pooled_values, &
            fit_valid, &
            fit_method
  
       if (ios /= 0) then
          write(*,*) &
               "ERROR: invalid row in D_diff Gaussian lookup:"
          write(*,*) trim(line)
          stop
       end if
  
       ! Table has alpha_index = 0,...,51.
       if (alpha_index < 0 .or. &
           alpha_index >= n_active_alpha) then
          write(*,*) &
               "ERROR: invalid alpha_index in Gaussian lookup: ", &
               alpha_index
          stop
       end if
  
       if (n_joint < 100 .or. &
           n_joint > nphi_ring) then
          write(*,*) &
               "ERROR: invalid N_joint in Gaussian lookup: ", &
               n_joint
          stop
       end if
  
       if (fit_valid /= 1) then
          write(*,*) &
               "ERROR: invalid Gaussian lookup cell: ", &
               alpha_index, n_joint
          stop
       end if
  
       if (sigma_value <= 0.0_dp) then
          write(*,*) &
               "ERROR: non-positive sigma in Gaussian lookup: ", &
               alpha_index, n_joint, sigma_value
          stop
       end if
  
       ! Direct mapping:
       !
       ! alpha_index=0  -> first active theta_ring
       ! alpha_index=1  -> next theta_ring
       ! ...
       ! alpha_index=51 -> last active theta_ring
       iopening_lookup = &
            first_iopening_angle + alpha_index
  
       d_diff_gaussian_mu( &
            iopening_lookup, n_joint) = &
            mu_value
  
       d_diff_gaussian_sigma( &
            iopening_lookup, n_joint) = &
            sigma_value
    end do
  
    close(lookup_unit)
  
    ! ==================================================
    ! Final completeness check.
    !
    ! Every active alpha must have mu/sigma for
    ! every N_joint = 100,...,512.
    ! ==================================================
  
    do ith = first_iopening_angle, &
         first_iopening_angle + n_active_alpha - 1
  
       do n = 100, nphi_ring
          if (d_diff_gaussian_sigma(ith,n) <= 0.0_dp) then
  
             write(*,*) &
                  "ERROR: missing Gaussian lookup cell:"
             write(*,*) &
                  "iopening_angle = ", ith, &
                  " N_joint = ", n
             stop
          end if
       end do
  
    end do
  
    if (verbosity > 0) then
       write(*,*) &
            "Loaded complete D_diff Gaussian lookup"
  
       write(*,*) &
            "alpha cells = ", n_active_alpha, &
            ", N_joint = 100...", nphi_ring
    end if
  
  end subroutine load_d_diff_gaussian_lookup

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

  subroutine get_temperature_along_circle(circle, nside_map, map, pixel_average, nest2ring_table, t, obs_mask, mask_threshold, valid)
    implicit none
    type(tcircle), intent(in) :: circle
    integer,  intent(in) :: nside_map
    real(sp), intent(in) :: map(0:)
    logical,  intent(in) :: pixel_average
    integer,  intent(in) :: nest2ring_table(0:)
    real(FP), intent(out):: t(:)
    real(sp), intent(in), optional :: obs_mask(0:)
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
             t(iphi) = DBLE((map(jpix)*2.6137104+(map(nb(1))+map(nb(3))+map(nb(5))+map(nb(7)))*0.9583863776+(map(nb(2))+map(nb(4))+map(nb(6))+map(nb(8)))*0.692777175)/9.21836461)
          else
             t(iphi) = DBLE((map(jpix)*2.6137104+(map(nb(1))+map(nb(3))+map(nb(5))+map(nb(7)))*0.9583863776+(map(nb(2))+map(nb(4))+map(nb(6)))*0.692777175)/8.525587431)
          endif
       else
          t(iphi) = DBLE(map(jpix))
       end if
    end do
  end subroutine get_temperature_along_circle

  subroutine get_symmetry_ring_and_validity(circle, nside_map, map, t, valid, &
     & sample_pixel, sample_theta_deg, sample_phi_deg)

    implicit none

    type(tcircle), intent(in) :: circle
    integer,       intent(in) :: nside_map
    real(sp),      intent(in) :: map(0:)
    real(FP),      intent(out) :: t(:)
    logical,       intent(out) :: valid(:)

    integer,  intent(out) :: sample_pixel(:)
    real(dp), intent(out) :: sample_theta_deg(:)
    real(dp), intent(out) :: sample_phi_deg(:)

    real(dp) :: rotation_matrix(1:3,1:3)
    real(dp) :: n_ipix(1:3), n_kpix(1:3)
    real(dp) :: theta_k, phi_k, phi_ring, sinang
    integer  :: jpix, iphi

    call get_rotation_matrix(circle%phi, circle%theta, rotation_matrix)

    sinang = sin(circle%opening_angle)
    n_ipix(3) = cos(circle%opening_angle)

    do iphi = 1, circle%nphi_ring

       phi_ring = 2.0_dp * Pi * dble(iphi) / dble(circle%nphi_ring)

       n_ipix(1) = sinang * cos(phi_ring)
       n_ipix(2) = sinang * sin(phi_ring)

       n_kpix = matmul(rotation_matrix, n_ipix)

       theta_k = acos(n_kpix(3))
       if (theta_k >= Pi) theta_k = theta_k - 0.00000001_dp

       phi_k = atan2(n_kpix(2), n_kpix(1)) + Pi

       call ang2pix_ring(nside_map, theta_k, phi_k, jpix)

       sample_pixel(iphi) = jpix
       sample_theta_deg(iphi) = rad_to_deg(theta_k)
       sample_phi_deg(iphi) = rad_to_deg(phi_k)

       valid(iphi) = dble(map(jpix)) > -1.0e29_dp

       if (valid(iphi)) then
          t(iphi) = dble(map(jpix))
       else
          t(iphi) = 0.0
       end if

    end do
    
  end subroutine get_symmetry_ring_and_validity

  subroutine fill_circle_transforms(npix, circle, nside_map, map, use_pixel_average, nest2ring_table, circle_transforms)
    implicit none
    integer,       intent(in) :: npix
    type(tcircle), intent(inout) :: circle
    integer,       intent(in) :: nside_map
    real,          intent(in) :: map(0:)
    logical,       intent(in) :: use_pixel_average
    integer,       intent(in) :: nest2ring_table(0:)
    FFT_COMPLEX,   intent(out):: circle_transforms(1:,0:)

    FFT_COMPLEX, dimension(circle%nphi_ring/2+1) :: plan_out, temp_out
    real(FP), dimension(circle%nphi_ring) :: plan_in, temp_in
    integer(8) :: plan
    integer :: i

#ifdef DFFT
    call dfftw_plan_dft_r2c_1d(plan, int(circle%nphi_ring), plan_in, plan_out, FFTW_ESTIMATE)
#else
    call sfftw_plan_dft_r2c_1d(plan, int(circle%nphi_ring), plan_in, plan_out, FFTW_ESTIMATE)
#endif

    !$omp parallel private(i, temp_in, temp_out) firstprivate(circle)
    !$omp do
    do i=0, npix-1
       circle%pixel=i
       circle%theta=pix2ang_table(circle%pixel)%theta
       circle%phi=pix2ang_table(circle%pixel)%phi
       call get_temperature_along_circle(circle, nside_map, map, use_pixel_average, nest2ring_table, temp_in)
#ifdef DFFT
       call dfftw_execute_dft_r2c(plan, temp_in, temp_out)
#else
       call sfftw_execute_dft_r2c(plan, temp_in, temp_out)
#endif
       circle_transforms(1:circle%nphi_ring/2+1,i)=temp_out(1:circle%nphi_ring/2+1)
    end do
    !$omp end do
    !$omp end parallel

#ifdef DFFT
    call dfftw_destroy_plan(plan)
#else
    call sfftw_destroy_plan(plan)
#endif

  end subroutine fill_circle_transforms

  subroutine fill_d_transforms( &
       & npix, circle, nside_map, map, pixel_average, &
       & nest2ring_table, obs_mask, mask_threshold, &
       & w_t_transforms, w_transforms, w_t_sq_transforms)

    implicit none

    integer,       intent(in)    :: npix
    type(tcircle), intent(inout) :: circle
    integer,       intent(in)    :: nside_map
    real(sp),      intent(in)    :: map(0:)
    logical,       intent(in)    :: pixel_average
    integer,       intent(in)    :: nest2ring_table(0:)
    real(sp),      intent(in)    :: obs_mask(0:)
    real(dp),      intent(in)    :: mask_threshold

    FFT_COMPLEX, intent(out) :: w_t_transforms(1:,0:)
    FFT_COMPLEX, intent(out) :: w_transforms(1:,0:)
    FFT_COMPLEX, intent(out) :: w_t_sq_transforms(1:,0:)

    integer :: i
    integer :: iphi
    integer :: nphi

    real(FP) :: ring(circle%nphi_ring)
    real(FP) :: w_t(circle%nphi_ring)
    real(FP) :: w(circle%nphi_ring)
    real(FP) :: w_t_sq(circle%nphi_ring)

    FFT_COMPLEX :: modes(circle%nphi_ring/2 + 1)

    logical :: valid(circle%nphi_ring)

    nphi = circle%nphi_ring


    !$omp parallel private( &
    !$omp& i, iphi, ring, w_t, w, w_t_sq, modes, valid) &
    !$omp& firstprivate(circle)

    !$omp do schedule(dynamic)
    do i = 0, npix - 1

       circle%pixel = i
       circle%theta = pix2ang_table(i)%theta
       circle%phi   = pix2ang_table(i)%phi

       call get_temperature_along_circle( &
            & circle, nside_map, map, pixel_average, &
            & nest2ring_table, ring, obs_mask, &
            & mask_threshold, valid)


#ifdef DFFT
       call dfftw_execute_dft_r2c( &
            & plan_forward, ring, modes)
#else
       call sfftw_execute_dft_r2c( &
            & plan_forward, ring, modes)
#endif

       modes(1:m_min) = 0.0

#ifdef DFFT
       call dfftw_execute_dft_c2r( &
            & plan_backward, modes, ring)
#else
       call sfftw_execute_dft_c2r( &
            & plan_backward, modes, ring)
#endif

       ring = ring / dble(nphi)


       w      = 0.0
       w_t    = 0.0
       w_t_sq = 0.0

       do iphi = 1, nphi

          if (.not. valid(iphi)) cycle

          w(iphi)      = 1.0
          w_t(iphi)    = ring(iphi)
          w_t_sq(iphi) = ring(iphi)**2

       end do


#ifdef DFFT
       call dfftw_execute_dft_r2c( &
            & plan_forward, w_t, &
            & w_t_transforms(:,i))

       call dfftw_execute_dft_r2c( &
            & plan_forward, w, &
            & w_transforms(:,i))

       call dfftw_execute_dft_r2c( &
            & plan_forward, w_t_sq, &
            & w_t_sq_transforms(:,i))
#else
       call sfftw_execute_dft_r2c( &
            & plan_forward, w_t, &
            & w_t_transforms(:,i))

       call sfftw_execute_dft_r2c( &
            & plan_forward, w, &
            & w_transforms(:,i))

       call sfftw_execute_dft_r2c( &
            & plan_forward, w_t_sq, &
            & w_t_sq_transforms(:,i))
#endif

    end do
    !$omp end do
    !$omp end parallel

  end subroutine fill_d_transforms

  subroutine fill_d_diff_transforms( &
       & npix, circle, nside_map, map, pixel_average, &
       & nest2ring_table, obs_mask, mask_threshold, &
       & w_delta_t_transforms, w_transforms, &
       & w_delta_t_sq_transforms)

    implicit none

    integer,       intent(in)    :: npix
    type(tcircle), intent(inout) :: circle
    integer,       intent(in)    :: nside_map
    real(sp),      intent(in)    :: map(0:)
    logical,       intent(in)    :: pixel_average
    integer,       intent(in)    :: nest2ring_table(0:)
    real(sp),      intent(in)    :: obs_mask(0:)
    real(dp),      intent(in)    :: mask_threshold

    FFT_COMPLEX, intent(out) :: w_delta_t_transforms(1:,0:,1:)
    FFT_COMPLEX, intent(out) :: w_transforms(1:,0:,1:)
    FFT_COMPLEX, intent(out) :: w_delta_t_sq_transforms(1:,0:,1:)

    integer :: i
    integer :: iphi
    integer :: next_iphi
    integer :: nphi
    integer :: separation_bin

    real(FP) :: delta_t

    real(FP) :: ring(circle%nphi_ring)
    real(FP) :: w_delta_t(circle%nphi_ring)
    real(FP) :: w(circle%nphi_ring)
    real(FP) :: w_delta_t_sq(circle%nphi_ring)

    FFT_COMPLEX :: modes(circle%nphi_ring/2 + 1)

    logical :: valid(circle%nphi_ring)

    nphi = circle%nphi_ring

    !$omp parallel private( &
    !$omp& i, iphi, next_iphi, separation_bin, delta_t, ring, &
    !$omp& w_delta_t, w, w_delta_t_sq, modes, valid) &
    !$omp& firstprivate(circle)

    !$omp do schedule(dynamic)
    do i = 0, npix - 1

       circle%pixel = i
       circle%theta = pix2ang_table(i)%theta
       circle%phi   = pix2ang_table(i)%phi

       ! Sample the filtered temperature and the mask.
       call get_temperature_along_circle( &
            & circle, nside_map, map, pixel_average, &
            & nest2ring_table, ring, obs_mask, &
            & mask_threshold, valid)

       ! Apply the same m_min cut used by D_ij and D_diff.
#ifdef DFFT
       call dfftw_execute_dft_r2c( &
            & plan_forward, ring, modes)
#else
       call sfftw_execute_dft_r2c( &
            & plan_forward, ring, modes)
#endif

       modes(1:m_min) = 0.0

#ifdef DFFT
       call dfftw_execute_dft_c2r( &
            & plan_backward, modes, ring)
#else
       call sfftw_execute_dft_c2r( &
            & plan_backward, modes, ring)
#endif

       ring = ring / dble(nphi)

       ! Construct the required quantities separately for
       ! real-space separations s = 1,...,d_diff_max_separation.
       do separation_bin = 1, d_diff_max_separation

          w = 0.0
          w_delta_t = 0.0
          w_delta_t_sq = 0.0
       
          do iphi = 1, nphi
       
             next_iphi = modulo( &
                  & iphi - 1 + separation_bin, nphi) + 1
       
             ! Both endpoints of this difference must be observed.
             if (.not. ( &
                  & valid(iphi) .and. valid(next_iphi) &
                  & )) cycle
       
             delta_t = ring(next_iphi) - ring(iphi)
       
             w(iphi) = 1.0
             w_delta_t(iphi) = delta_t
             w_delta_t_sq(iphi) = delta_t**2
       
          end do
       
          ! Store the three Fourier transforms for this
          ! circle and this separation.
#ifdef DFFT
          call dfftw_execute_dft_r2c( &
               & plan_forward, w_delta_t, &
               & w_delta_t_transforms(:,i,separation_bin))
       
          call dfftw_execute_dft_r2c( &
               & plan_forward, w, &
               & w_transforms(:,i,separation_bin))
       
          call dfftw_execute_dft_r2c( &
               & plan_forward, w_delta_t_sq, &
               & w_delta_t_sq_transforms(:,i,separation_bin))
#else
          call sfftw_execute_dft_r2c( &
               & plan_forward, w_delta_t, &
               & w_delta_t_transforms(:,i,separation_bin))
       
          call sfftw_execute_dft_r2c( &
               & plan_forward, w, &
               & w_transforms(:,i,separation_bin))
       
          call sfftw_execute_dft_r2c( &
               & plan_forward, w_delta_t_sq, &
               & w_delta_t_sq_transforms(:,i,separation_bin))
#endif
       end do

    end do
    !$omp end do
    !$omp end parallel

  end subroutine fill_d_diff_transforms

  subroutine get_validity_along_circle(circle, nside_map, obs_mask, &
       & mask_threshold, valid)

    implicit none

    type(tcircle), intent(in) :: circle
    integer,       intent(in) :: nside_map
    real(sp),      intent(in) :: obs_mask(0:)
    real(dp),      intent(in) :: mask_threshold
    logical,       intent(out) :: valid(:)

    real(dp) :: rotation_matrix(1:3,1:3)
    real(dp) :: n_ipix(1:3), n_kpix(1:3)
    real(dp) :: theta_k, phi_k, phi_ring, sinang
    integer  :: iphi, jpix

    call get_rotation_matrix(circle%phi, circle%theta, rotation_matrix)

    sinang = sin(circle%opening_angle)
    n_ipix(3) = cos(circle%opening_angle)

    do iphi = 1, circle%nphi_ring

       phi_ring = 2.0_dp * Pi * dble(iphi) / dble(circle%nphi_ring)

       n_ipix(1) = sinang * cos(phi_ring)
       n_ipix(2) = sinang * sin(phi_ring)

       n_kpix = matmul(rotation_matrix, n_ipix)

       theta_k = acos(n_kpix(3))
       if (theta_k >= Pi) theta_k = theta_k - 0.00000001_dp

       phi_k = atan2(n_kpix(2), n_kpix(1)) + Pi

       call ang2pix_ring(nside_map, theta_k, phi_k, jpix)

       valid(iphi) = dble(obs_mask(jpix)) >= mask_threshold

    end do

  end subroutine get_validity_along_circle

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

    if (run_symmetry_test .and. symmetry_input_has_missing_pixels) then
       call input_map(trim(adjustl(mapfilename)), in_map, 12*nside**2, 1, &
            & fmissval=-1.0e30_sp)
    else
       call input_map(trim(adjustl(mapfilename)), in_map, 12*nside**2, 1)
    end if

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

  subroutine analyse_circle_pairs(circle1, circle2, circle_transforms, &
       & opposite_list, opposite_number, smax, bestcircle1, bestcircle2, &
       & position_smax, position_circle2, position_beta_max_bin, &
       & analysed_circles, every_nth_smax_file)

    implicit none

    type(tcircle), intent(in)    :: circle1
    type(tcircle), intent(inout) :: circle2

    FFT_COMPLEX, intent(in) :: circle_transforms(1:,0:)

    integer, intent(in) :: opposite_list(0:)
    integer, intent(in) :: opposite_number
    real(dp), intent(inout) :: smax
    type(tcircle), intent(inout) :: bestcircle1, bestcircle2
    real(dp),      intent(out) :: position_smax
    type(tcircle), intent(out) :: position_circle2
    integer,       intent(out) :: position_beta_max_bin
    integer(dp), intent(inout) :: analysed_circles
    integer, intent(in) :: every_nth_smax_file
    integer  :: j
    integer  :: beta_try
    real(dp) :: smax_try
    position_smax = -huge(1.0_dp)
    position_beta_max_bin = -1
    position_circle2%pixel = -1

    do j = 0, opposite_number - 1

       circle2%pixel = opposite_list(j)

       if (.not. write_position_smax) then
          if (circle2%pixel .lt. circle1%pixel) cycle
       end if

       call get_smax(circle1, circle2, circle_transforms, &
            & smax_try, random_s, beta_try)

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

  subroutine analyse_circle_pairs_d( &
       & circle1, circle2, opposite_list, opposite_number, &
       & d_best_at_alpha, &
       & best_circle1_at_alpha, best_circle2_at_alpha, &
       & best_beta_at_alpha, best_n_joint_at_alpha, &
       & found_pair_at_alpha)

    implicit none

    type(tcircle), intent(in)    :: circle1
    type(tcircle), intent(inout) :: circle2

    integer, intent(in) :: opposite_list(0:)
    integer, intent(in) :: opposite_number

    real(dp), intent(inout) :: d_best_at_alpha

    type(tcircle), intent(inout) :: best_circle1_at_alpha
    type(tcircle), intent(inout) :: best_circle2_at_alpha

    integer, intent(inout) :: best_beta_at_alpha
    integer, intent(inout) :: best_n_joint_at_alpha

    logical, intent(inout) :: found_pair_at_alpha

    integer :: j
    integer :: beta_try
    integer :: n_joint_try

    real(dp) :: d_try

    do j = 0, opposite_number - 1

       circle2%pixel = opposite_list(j)

       ! Avoid evaluating both (i,j) and (j,i).
       if (circle2%pixel < circle1%pixel) cycle

       call get_d_max( &
            & circle1, circle2, d_try, &
            & beta_try, n_joint_try)

       if (n_joint_try < d_min_joint_samples) cycle

       if (.not. found_pair_at_alpha .or. &
           d_try > d_best_at_alpha) then

          d_best_at_alpha = d_try

          best_circle1_at_alpha = circle1
          best_circle2_at_alpha = circle2

          best_beta_at_alpha = beta_try
          best_n_joint_at_alpha = n_joint_try

          found_pair_at_alpha = .true.
       end if
    end do

  end subroutine analyse_circle_pairs_d

  subroutine analyse_circle_pairs_d_diff( &
       & circle1, circle2, opposite_list, opposite_number, &
       & d_diff_best_at_alpha, &
       & best_circle1_at_alpha, best_circle2_at_alpha, &
       & best_beta_at_alpha, best_n_joint_at_alpha, &
       & found_pair_at_alpha)
  
    implicit none
  
    ! Fixed circle i and fixed radius alpha.
    type(tcircle), intent(in)    :: circle1
  
    ! Candidate circle j, with the same fixed radius alpha.
    type(tcircle), intent(inout) :: circle2
  
    ! Allowed j-centres for the current circle i.
    integer, intent(in) :: opposite_list(0:)
    integer, intent(in) :: opposite_number
  
    ! Best result found so far at this fixed alpha.
    real(dp), intent(inout) :: d_diff_best_at_alpha
  
    type(tcircle), intent(inout) :: best_circle1_at_alpha
    type(tcircle), intent(inout) :: best_circle2_at_alpha
  
    integer, intent(inout) :: best_beta_at_alpha
    integer, intent(inout) :: best_n_joint_at_alpha
  
    logical, intent(inout) :: found_pair_at_alpha
  
    ! Result for the current circle pair (i,j).
    integer  :: j
    integer  :: beta_try
    integer  :: n_joint_try
    integer  :: d_diff_hist_bin
    real(dp) :: d_diff_try
  
    ! Loop only over the allowed circle centres j.
    do j = 0, opposite_number - 1
  
       circle2%pixel = opposite_list(j)
  
       ! Avoid evaluating both (i,j) and (j,i).
       if (circle2%pixel < circle1%pixel) cycle
  
       ! At this fixed (i,j,alpha), maximise only over beta.
       call get_d_diff_max( &
            & circle1, circle2, d_diff_try, &
            & beta_try, n_joint_try)
  
       ! No beta had enough jointly valid differences.
       if (n_joint_try < d_min_joint_samples) cycle

       if (use_d_diff_bias_correction) then
          d_diff_try = d_diff_try - &
               & d_diff_bias_mean( &
               & rad_to_deg(dble(circle1%opening_angle)), &
               & n_joint_try)
       end if

       ! Record D_ij,max-beta^diff using its corresponding N_joint^diff.
       if (write_d_diff_null_histogram) then
       
          d_diff_hist_bin = int( &
               & 0.5_dp * (d_diff_try + 1.0_dp) * &
               & dble(d_diff_hist_nbins))
       
          d_diff_hist_bin = max( &
               & 0, min(d_diff_hist_nbins - 1, &
               & d_diff_hist_bin))
       
          !$omp atomic update
          d_diff_null_hist_counts( &
               & n_joint_try, d_diff_hist_bin) = &
               & d_diff_null_hist_counts( &
               & n_joint_try, d_diff_hist_bin) + 1_8
       
       end if
  
       ! Keep the largest result found so far at this alpha.
       if (.not. found_pair_at_alpha .or. &
           & d_diff_try > d_diff_best_at_alpha) then
  
          d_diff_best_at_alpha = d_diff_try
  
          best_circle1_at_alpha = circle1
          best_circle2_at_alpha = circle2
  
          best_beta_at_alpha = beta_try
          best_n_joint_at_alpha = n_joint_try
  
          found_pair_at_alpha = .true.
  
       end if
  
    end do  
  end subroutine analyse_circle_pairs_d_diff

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

  subroutine run_symmetry_test_mode()
    implicit none
    
    integer  :: iphi, partner, beta_bin
    integer  :: beta_max_bin, beta_min_bin
    integer  :: n_valid, n_masked, ialpha, n_alpha
    real(dp) :: numerator, denominator
    real(dp) :: dvalue, dmax, dmin, alpha_deg
    logical, allocatable :: ring_valid(:)
    integer :: sample_out
    integer, allocatable :: sample_pixel(:)
    real(dp), allocatable :: sample_theta_deg(:)
    real(dp), allocatable :: sample_phi_deg(:)

    if (symmetry_center_theta_deg < 0.0_dp .or. &
        symmetry_center_theta_deg > 180.0_dp) then
       write(*,*) "ERROR: symmetry_center_theta_deg must lie between 0 and 180."
       stop
    end if

    circle1%theta = symmetry_center_theta_deg * Pi / 180.0_dp
    circle1%phi   = modulo(symmetry_center_phi_deg, 360.0_dp) * Pi / 180.0_dp
    allocate(ring_valid(nphi_ring))
    allocate(sample_pixel(nphi_ring))
    allocate(sample_theta_deg(nphi_ring))
    allocate(sample_phi_deg(nphi_ring))
    
    open(newunit=sample_out, &
         file='symmetry_samples.dat', &
         status='replace', action='write')
    
    write(sample_out,'(A)') &
         '# alpha_deg iphi phi_ring_deg pixel theta_deg phi_deg b_deg valid T'

    write(*,'(A)') '# alpha_deg n_valid n_masked Dmax minus_Dmin beta_max_deg beta_min_deg'

    n_alpha = opening_angle_range%num_steps

    do ialpha = 0, n_alpha - 1

       alpha_deg = dble(opening_angle_range%min) + &
            & dble(ialpha) * &
            & (dble(opening_angle_range%max) - &
            &  dble(opening_angle_range%min)) / &
            & dble(n_alpha - 1)

       circle1%opening_angle = alpha_deg * Pi / 180.0_dp
        
       call get_symmetry_ring_and_validity(circle1, nside_map, &
            & map, temperature_in, ring_valid, sample_pixel, &
            & sample_theta_deg, sample_phi_deg)

    do iphi = 1, nphi_ring
    
       write(sample_out,*) &
            alpha_deg, &
            iphi, &
            360.0_dp * dble(iphi) / dble(nphi_ring), &
            sample_pixel(iphi), &
            sample_theta_deg(iphi), &
            sample_phi_deg(iphi), &
            90.0_dp - sample_theta_deg(iphi), &
            ring_valid(iphi), &
            temperature_in(iphi)  
    end do

    n_valid = count(ring_valid)
    n_masked = nphi_ring - n_valid

    ! Transform the sampled real-space ring T_k into Fourier coefficients.
#ifdef DFFT
    call dfftw_execute_dft_r2c(plan_forward, temperature_in, temperature_out)
#else
    call sfftw_execute_dft_r2c(plan_forward, temperature_in, temperature_out)
#endif

    ! Copy the Fourier coefficients, then remove low ring modes.
    ! With m_min = 2, this removes m = 0 and m = +/- 1.
    sm_in = temperature_out
    sm_in(1:m_min) = 0.0

    ! Transform back to the filtered real-space ring.
#ifdef DFFT
    call dfftw_execute_dft_c2r(plan_backward, sm_in, sm_out)
#else
    call sfftw_execute_dft_c2r(plan_backward, sm_in, sm_out)
#endif

    ! FFTW's inverse transform is not normalised.
    sm_out = sm_out / dble(nphi_ring)

    ! Scan every front-to-front phase shift beta.
    dmax = -huge(1.0_dp)
    dmin =  huge(1.0_dp)

    beta_max_bin = -1
    beta_min_bin = -1

    do beta_bin = 0, nphi_ring - 1

       numerator = 0.0_dp
       denominator = 0.0_dp

       do iphi = 1, nphi_ring

          partner = modulo(beta_bin - (iphi - 1), nphi_ring) + 1

          if (.not. (ring_valid(iphi) .and. ring_valid(partner))) cycle

          numerator = numerator + 2.0_dp * dble(sm_out(iphi)) * &
               & dble(sm_out(partner))

          denominator = denominator + dble(sm_out(iphi))**2 + &
               & dble(sm_out(partner))**2

       end do

       if (denominator <= tiny(1.0_dp)) cycle

       dvalue = numerator / denominator

       if (dvalue > dmax) then
          dmax = dvalue
          beta_max_bin = beta_bin
       end if

       if (dvalue < dmin) then
          dmin = dvalue
          beta_min_bin = beta_bin
       end if

    end do

    if (beta_max_bin < 0 .or. beta_min_bin < 0) then
       write(*,*) "ERROR: filtered ring has zero power."
       stop
    end if

       write(*,'(F8.3,1X,I6,1X,I6,4(1X,F14.8))') &
            & alpha_deg, n_valid, n_masked, dmax, -dmin, &
            & 360.0_dp * dble(beta_max_bin) / dble(nphi_ring), &
            & 360.0_dp * dble(beta_min_bin) / dble(nphi_ring)
    end do
    close(sample_out)

    deallocate(sample_pixel)
    deallocate(sample_theta_deg)
    deallocate(sample_phi_deg)
    deallocate(ring_valid)

  end subroutine run_symmetry_test_mode

end program circles
