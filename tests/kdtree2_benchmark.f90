program kdtree_benchmark
  use kdtree2_module
  implicit none

  type(kdtree2) :: tree
  real(kdkind), dimension(:,:), allocatable :: data, queries
  type(kdtree2_result) :: results(1)

  real :: t0, t1, t2
  integer, parameter :: m(3) = [3, 8, 16]
  integer :: n = 10000, r = 1000
  integer :: idx
  integer :: i, s
  integer, allocatable :: seed(:), idxs(:)
  character(len=80) :: arg_str

  if (command_argument_count() == 2) then
     call get_command_argument(1, arg_str)
     read(arg_str, *) n
     call get_command_argument(2, arg_str)
     read(arg_str, *) r
  end if

  write(*, "(A,I0,A,I0)") " Using n_points = ", n, ", n_samples = ", r

  allocate(idxs(r))

  call random_seed(size=s)
  allocate(seed(s))
  seed = 1234
  call random_seed(put=seed)

  do i = 1, 3

    write(*,*) "Dimension = ", m(i)

    allocate(data(m(i),n))
    call random_number(data)

    allocate(queries(m(i),r))
    call random_number(queries)

    ! Populate tree
    call cpu_time(t0)
    tree = kdtree2_create(data,sort=.true.,rearrange=.true.)
    call cpu_time(t1)

    ! Query random vectors
    do idx = 1, r
      call kdtree2_n_nearest(tp=tree,qv=queries(:,idx),nn=1,results=results)
      idxs(idx) = results(1)%idx
    end do
    call cpu_time(t2)

    write(*,*) "Tree build (s): ", t1 - t0
    write(*,*) "Tree query (s): ", t2 - t1

    call kdtree2_destroy(tree)  
    deallocate(data)
    deallocate(queries)

  end do

end program
