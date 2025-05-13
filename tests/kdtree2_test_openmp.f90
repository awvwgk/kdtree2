! Test kdtree2 module with OpenMP
program test_openmp
  use iso_fortran_env, only: int64
  use kdtree2_module

  implicit none

  type(kdtree2)                     :: tree
  real(kdkind), allocatable         :: points(:, :)
  real(kdkind), allocatable         :: query_points(:, :)
  type(kdtree2_result), allocatable :: results(:)
  integer, allocatable              :: nearest_neighbors_seq(:, :)
  integer, allocatable              :: nearest_neighbors_omp(:, :)
  integer                           :: n, n_dim, n_points
  integer                           :: n_queries, n_neighbors
  character(len=80)                 :: arg_str
  integer(int64)                    :: t_start, t_end, count_rate
  real(kdkind)                      :: wallclock_time

  if (command_argument_count() == 4) then
     call get_command_argument(1, arg_str)
     read(arg_str, *) n_dim
     call get_command_argument(2, arg_str)
     read(arg_str, *) n_points
     call get_command_argument(3, arg_str)
     read(arg_str, *) n_queries
     call get_command_argument(4, arg_str)
     read(arg_str, *) n_neighbors
  else
     print *, "Usage: ./kdtree2_test_openmp n_dim n_points n_queries n_neighbors"
     n_dim = 3
     n_points = 100*1000
     n_queries = 10*1000
     n_neighbors = 1
  end if

  write(*, "(A20,I12)") "n_dim               ", n_dim
  write(*, "(A20,I12)") "n_points            ", n_points
  write(*, "(A20,I12)") "n_queries           ", n_queries
  write(*, "(A20,I12)") "n_neighbors         ", n_neighbors

  allocate(points(n_dim, n_points))
  allocate(query_points(n_dim, n_queries))
  allocate(nearest_neighbors_omp(n_neighbors, n_queries))
  allocate(nearest_neighbors_seq(n_neighbors, n_queries))
  allocate(results(n_neighbors))

  call random_number(points)
  call random_number(query_points)

  tree = kdtree2_create(points, sort=.true., rearrange=.true.)

  call system_clock(t_start, count_rate)
  do n = 1, n_queries
     call kdtree2_n_nearest(tree, query_points(:, n), n_neighbors, results)
     nearest_neighbors_seq(:, n) = results(:)%idx
  end do
  call system_clock(t_end, count_rate)
  wallclock_time = (t_end-t_start) / real(count_rate, kdkind)

  write(*, "(A20,E12.3,A)") "Wallclock-seq.      ", wallclock_time
  write(*, "(A20,E12.3,A)") "Queries/s-seq.      ", n_queries/wallclock_time

  call system_clock(t_start, count_rate)
  !$omp parallel do private(results)
  do n = 1, n_queries
     call kdtree2_n_nearest(tree, query_points(:, n), n_neighbors, results)
     nearest_neighbors_omp(:, n) = results(:)%idx
  end do
  !$omp end parallel do
  call system_clock(t_end, count_rate)
  wallclock_time = (t_end-t_start) / real(count_rate, kdkind)

  write(*, "(A20,E12.3,A)") "Wallclock-par.      ", wallclock_time
  write(*, "(A20,E12.3,A)") "Queries/s-par.      ", n_queries/wallclock_time

  if (any(nearest_neighbors_omp /= nearest_neighbors_seq)) &
       error stop "Parallel and sequential results do not agree"

  call kdtree2_destroy(tree)
  deallocate(points, query_points, nearest_neighbors_seq, &
       nearest_neighbors_omp, results)

end program test_openmp
