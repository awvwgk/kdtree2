!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
! K-D tree routines in Fortran 90 by Matt Kennel.
! Original program was written in Sather by Steve Omohundro and
! Matt Kennel.  Only the Euclidean metric is supported.
!
!
! This module is identical to 'kd_tree', except that the order
! of subscripts is reversed in the data file.
! In otherwords for an embedding of N D-dimensional vectors, the
! data file is here, in natural Fortran order  data(1:D, 1:N)
! because Fortran lays out columns first,
! whereas conventionally (C-style) it is data(1:N,1:D)
! as in the original kd_tree module.
!
! ** Documentation for priority queue **
! A heap-based priority queue lets one efficiently implement the following
! operations, each in log(N) time, as opposed to linear time.
!
! 1)  add a datum (push a datum onto the queue, increasing its length)
! 2)  return the priority value of the maximum priority element
! 3)  pop-off (and delete) the element with the maximum priority, decreasing
!     the size of the queue.
! 4)  replace the datum with the maximum priority with a supplied datum
!     (of either higher or lower priority), maintaining the size of the
!     queue.
!
! In the k-d tree case, the 'priority' is the square distance of a point in
! the data set to a reference point.   The goal is to keep the smallest M
! distances to a reference point.  The tree algorithm searches terminal
! nodes to decide whether to add points under consideration.
!
! A priority queue is useful here because it lets one quickly return the
! largest distance currently existing in the list.  If a new candidate
! distance is smaller than this, then the new candidate ought to replace
! the old candidate.  In priority queue terms, this means removing the
! highest priority element, and inserting the new one.
!
! Algorithms based on Cormen, Leiserson, Rivest, _Introduction
! to Algorithms_, 1990, with further optimization by the author.
!
! Originally informed by a C implementation by Sriranga Veeraraghavan.
!
! This module is not written in the most clear way, but is implemented such
! for speed, as it its operations will be called many times during searches
! of large numbers of neighbors.
module kdtree2_module
  use iso_fortran_env, only: error_unit

  implicit none
  private

  integer, parameter :: kdkind = kind(0.0d0)

  ! Warn user if dimension exceeds number of points and this threshold
  integer, parameter :: warning_ndim_threshold = 20

  !-------------DATA TYPE, CREATION, DELETION---------------------
  public :: kdkind
  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create, kdtree2_destroy

  !-------------------SEARCH ROUTINES-----------------------------
  ! Return fixed number of nearest neighbors around arbitrary vector,
  ! or extant point in dataset, with decorrelation window.
  public :: kdtree2_n_nearest, kdtree2_n_nearest_around_point

  ! Return points within a fixed ball of arb vector/extant point
  public :: kdtree2_r_nearest, kdtree2_r_nearest_around_point

  ! Sort, in order of increasing distance, rseults from above.
  public :: kdtree2_sort_results

  ! Count points within a fixed ball of arb vector/extant point
  public :: kdtree2_r_count, kdtree2_r_count_around_point

  ! brute force of kdtree2_[n|r]_nearest
  public :: kdtree2_n_nearest_brute_force, kdtree2_r_nearest_brute_force

  ! The priority queue consists of elements priority(1:heap_size), with
  ! associated payload(:).
  type pq
    ! There are heap_size active elements. Assumes the allocation is always
    ! sufficient. Will NOT increase it to match.
    integer :: heap_size = 0
  end type pq

  type interval
    real(kdkind) :: lower, upper
  end type interval

  ! An internal tree node
  type :: tree_node

    private
    ! The dimension to cut
    integer :: cut_dim
    ! Where to cut the dimension
    real(kdkind) :: cut_val
    ! Improved cutoffs knowing the spread in child boxes.
    real(kdkind) :: cut_val_left, cut_val_right

    ! Child pointers
    ! Points included in this node are indexes[k] with k \in [l,u]
    integer :: l, u
    type(tree_node), pointer :: left, right
    type(interval), allocatable :: box(:)
  end type tree_node

  ! Global information about the tree, one per tree
  type :: kdtree2
    ! Dimensionality
    integer :: dimen = 0

    ! Total # of points
    integer :: n = 0

    ! Copy of the input data
    real(kdkind), allocatable :: input_data(:, :)

    !  IMPORTANT NOTE:  IT IS DIMENSIONED   input_data(1:d,1:N)
    !  which may be opposite of what may be conventional.
    !  This is, because in Fortran, the memory layout is such that
    !  the first dimension is in sequential order.  Hence, with
    !  (1:d,1:N), all components of the vector will be in consecutive
    !  memory locations.  The search time is dominated by the
    !  evaluation of distances in the terminal nodes.  Putting all
    !  vector components in consecutive memory location improves
    !  memory cache locality, and hence search speed, and may enable
    !  vectorization on some processors and compilers.

    ! Permuted index into the data, so that indexes[l..u] of some
    ! bucket represent the indexes of the actual points in that
    ! bucket.
    integer, allocatable :: ind(:)

    ! The maximum number of points to keep in a terminal node.
    integer :: bucket_size = -1

    ! do we always sort output results?
    logical       :: sort

    ! if (rearrange .eqv. .true.) then rearranged data has been stored
    logical       :: rearrange

    ! Rearranged input data
    real(kdkind), allocatable :: rearranged_data(:, :)

    ! Root pointer of the tree
    type(tree_node), pointer :: root => null()

    ! If .true. print warnings to stderr
    logical :: verbose
  end type kdtree2

  ! One of these is created for each search.
  type :: tree_search_record
    private
    integer           :: nn, nfound
    real(kdkind)      :: ballsize
    integer           :: centeridx = 999, correltime = 9999
    ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
    ! did the # of points found overflow the storage provided?
    logical           :: overflow
    real(kdkind), allocatable :: qv(:)  ! query vector
    type(pq) :: pq
  end type tree_search_record

  ! A pair of distances, indexes
  type kdtree2_result
     real(kdkind) :: dis
     integer      :: idx
  end type kdtree2_result

contains

  ! Create the actual tree structure, given an input array of data.
  !
  ! Note, input data is input_data(1:d,1:N), NOT the other way around.
  ! THIS IS THE REVERSE OF THE PREVIOUS VERSION OF THIS MODULE.
  ! The reason for it is cache friendliness, improving performance.
  !
  ! Optional arguments:  If 'dim' is specified, then the tree
  !                      will only search the first 'dim' components
  !                      of input_data, otherwise, dim is inferred
  !                      from SIZE(input_data,1).
  !
  !                      if sort .eqv. .true. then output results
  !                      will be sorted by increasing distance.
  !                      default=.false., as it is faster to not sort.
  !
  !                      if rearrange .eqv. .true. then an internal
  !                      copy of the data, rearranged by terminal node,
  !                      will be made for cache friendliness.
  !                      default=.true., as it speeds searches, but
  !                      building takes longer, and extra memory is used.
  function kdtree2_create(input_data, dim, sort, rearrange, &
       bucket_size, verbose) result(mr)
    type(kdtree2)                 :: mr
    integer, intent(in), optional :: dim
    logical, intent(in), optional :: sort
    logical, intent(in), optional :: rearrange
    integer, intent(in), optional :: bucket_size
    logical, intent(in), optional :: verbose
    real(kdkind)                  :: input_data(:, :)
    integer                       :: i

    if (present(dim)) then
      mr%dimen = dim
    else
      mr%dimen = size(input_data, 1)
    end if
    mr%n = size(input_data, 2)

    if (present(verbose)) then
       mr%verbose = verbose
    else
       mr%verbose = .true.
    end if

    if (mr%dimen > max(mr%n, warning_ndim_threshold) .and. mr%verbose) then
      write (error_unit, '(A,A,I0,A,I0,A)') 'kdtree2_create warning: ', &
           'n_dim = ', mr%dimen, ', n_points = ', mr%n, &
           ' - typically n_dim < n_points'
    end if

    allocate(mr%input_data(mr%dimen, mr%n))
    mr%input_data(:, :) = input_data(1:mr%dimen, :)

    if (present(bucket_size)) then
       mr%bucket_size = bucket_size
    else
       mr%bucket_size = 12
    end if

    call build_tree(mr)

    if (present(sort)) then
      mr%sort = sort
    else
      mr%sort = .false.
    end if

    if (present(rearrange)) then
      mr%rearrange = rearrange
    else
      mr%rearrange = .true.
    end if

    if (mr%rearrange) then
      allocate(mr%rearranged_data(mr%dimen, mr%n))
      do i = 1, mr%n
        mr%rearranged_data(:, i) = input_data(1:mr%dimen, mr%ind(i))
      end do
    end if

  end function kdtree2_create

  subroutine build_tree(tp)
    type(kdtree2), intent(inout) :: tp
    integer                      :: j
    type(tree_node), pointer     :: dummy => null()

    allocate(tp%ind(tp%n))
    do concurrent (j = 1:tp%n)
      tp%ind(j) = j
    end do
    tp%root => build_tree_for_range(tp, 1, tp%n, dummy)
  end subroutine build_tree

  recursive function build_tree_for_range(tp, l, u, parent) result(res)
    type(tree_node), pointer             :: res
    type(kdtree2), intent(inout)         :: tp
    type(tree_node), pointer, intent(in) :: parent
    integer, intent(In)                  :: l, u
    integer                              :: i, c, m, dimen, n_below_average
    logical                              :: recompute
    real(kdkind)                         :: average, balance

    ! first compute min and max
    dimen = tp%dimen
    allocate (res)
    allocate (res%box(dimen))

    ! First, compute an APPROXIMATE bounding box of all points associated with this node.
    if (u < l) then
      ! no points in this box
      nullify (res)
      return
    end if

    if ((u - l) <= tp%bucket_size) then
      !
      ! always compute true bounding box for terminal nodes.
      !
      do i = 1, dimen
        call spread_in_coordinate(tp, i, l, u, res%box(i))
      end do
      res%cut_dim = 0
      res%cut_val = 0.0
      res%l = l
      res%u = u
      res%left => null()
      res%right => null()
    else
      !
      ! modify approximate bounding box.  This will be an
      ! overestimate of the true bounding box, as we are only recomputing
      ! the bounding box for the dimension that the parent split on.
      !
      ! Going to a true bounding box computation would significantly
      ! increase the time necessary to build the tree, and usually
      ! has only a very small difference.  This box is not used
      ! for searching but only for deciding which coordinate to split on.
      !
      do i = 1, dimen
        recompute = .true.
        if (associated(parent)) then
          if (i .ne. parent%cut_dim) then
            recompute = .false.
          end if
        end if
        if (recompute) then
          call spread_in_coordinate(tp, i, l, u, res%box(i))
        else
          res%box(i) = parent%box(i)
        end if
      end do

      ! c is the identity of which coordinate has the greatest spread.
      c = maxloc(res%box(1:dimen)%upper - res%box(1:dimen)%lower, 1)

      ! Determine arithmetic average
      average = sum(tp%input_data(c, tp%ind(l:u)))/real(u - l + 1, kdkind)

      ! Determine how balanced a split on the average would be
      n_below_average = count(tp%input_data(c, tp%ind(l:u)) < average)
      balance = n_below_average / real(u - l + 1, kdkind)

      if (balance < 0.25_kdkind .or. balance > 0.75_kdkind) then
        ! Bad balance, so select exact median to have 'perfect' balance
        m = (l + u)/2
        call select_on_coordinate(tp%input_data, tp%ind, c, m, l, u)
      else
        m = select_on_coordinate_value(tp%input_data, tp%ind, c, average, l, u)
      end if

      ! moves indexes around
      res%cut_dim = c
      res%l = l
      res%u = u

      res%left => build_tree_for_range(tp, l, m, res)
      res%right => build_tree_for_range(tp, m + 1, u, res)

      if (associated(res%right) .eqv. .false.) then
        res%box = res%left%box
        res%cut_val_left = res%left%box(c)%upper
        res%cut_val = res%cut_val_left
      elseif (associated(res%left) .eqv. .false.) then
        res%box = res%right%box
        res%cut_val_right = res%right%box(c)%lower
        res%cut_val = res%cut_val_right
      else
        res%cut_val_right = res%right%box(c)%lower
        res%cut_val_left = res%left%box(c)%upper
        res%cut_val = (res%cut_val_left + res%cut_val_right)/2

        ! now remake the true bounding box for self.
        ! Since we are taking unions (in effect) of a tree structure,
        ! this is much faster than doing an exhaustive
        ! search over all points
        res%box%upper = max(res%left%box%upper, res%right%box%upper)
        res%box%lower = min(res%left%box%lower, res%right%box%lower)
      end if
    end if
  end function build_tree_for_range

  ! Move elts of ind around between l and u, so that all points
  ! <= than alpha (in c cooordinate) are first, and then
  ! all points > alpha are second.
  !
  ! Algorithm (matt kennel).
  !
  ! Consider the list as having three parts: on the left,
  ! the points known to be <= alpha.  On the right, the points
  ! known to be > alpha, and in the middle, the currently unknown
  ! points.   The algorithm is to scan the unknown points, starting
  ! from the left, and swapping them so that they are added to
  ! the left stack or the right stack, as appropriate.
  !
  ! The algorithm finishes when the unknown stack is empty.
  integer function select_on_coordinate_value(v, ind, c, alpha, li, ui) result(res)
    integer, intent(In)       :: c, li, ui
    real(kdkind), intent(in)  :: alpha
    real(kdkind) , intent(in) :: v(1:,1:)
    integer, intent(inout)    :: ind(1:)

    integer :: tmp
    integer :: lb, rb

    ! The points known to be <= alpha are in
    ! [l,lb-1]
    !
    ! The points known to be > alpha are in
    ! [rb+1,u].
    !
    ! Therefore we add new points into lb or
    ! rb as appropriate.  When lb=rb
    ! we are done.  We return the location of the last point <= alpha.

    lb = li
    rb = ui

    do while (lb < rb)
      if (v(c, ind(lb)) <= alpha) then
        ! it is good where it is.
        lb = lb + 1
      else
        ! swap it with rb.
        tmp = ind(lb)
        ind(lb) = ind(rb)
        ind(rb) = tmp
        rb = rb - 1
      end if
    end do

    ! now lb .eq. ub
    if (v(c, ind(lb)) <= alpha) then
      res = lb
    else
      res = lb - 1
    end if

  end function select_on_coordinate_value

  ! Move elts of ind around between l and u, so that the kth element is >=
  ! those below, <= those above, in the coordinate c.
  subroutine select_on_coordinate(v, ind, c, k, li, ui)
    integer, intent(In) :: c, k, li, ui
    integer             :: i, l, m, s, t, u
    real(kdkind)        :: v(:, :)
    integer             :: ind(:)

    l = li
    u = ui
    do while (l < u)
      t = ind(l)
      m = l
      do i = l + 1, u
        if (v(c, ind(i)) < v(c, t)) then
          m = m + 1
          s = ind(m)
          ind(m) = ind(i)
          ind(i) = s
        end if
      end do
      s = ind(l)
      ind(l) = ind(m)
      ind(m) = s
      if (m <= k) l = m + 1
      if (m >= k) u = m - 1
    end do
  end subroutine select_on_coordinate

  ! The spread in coordinate 'c', between l and u.
  !
  ! Return lower bound in 'smin', and upper in 'smax',
  subroutine spread_in_coordinate(tp, c, l, u, interv)
    type(kdtree2), intent(in)   :: tp
    type(interval), intent(out) :: interv
    integer, intent(In)         :: c, l, u
    real(kdkind)                :: last, lmax, lmin, t, smin, smax
    integer                     :: i, ulocal

    smin = tp%input_data(c, tp%ind(l))
    smax = smin

    ulocal = u

    do i = l + 2, ulocal, 2
      lmin = tp%input_data(c, tp%ind(i - 1))
      lmax = tp%input_data(c, tp%ind(i))
      if (lmin > lmax) then
        t = lmin
        lmin = lmax
        lmax = t
      end if
      if (smin > lmin) smin = lmin
      if (smax < lmax) smax = lmax
    end do
    if (i == ulocal + 1) then
      last = tp%input_data(c, tp%ind(ulocal))
      if (smin > last) smin = last
      if (smax < last) smax = last
    end if

    interv%lower = smin
    interv%upper = smax

  end subroutine spread_in_coordinate

  ! Deallocates all memory for the tree, except input data matrix
  subroutine kdtree2_destroy(tp)
    type(kdtree2), intent(inout) :: tp

    call destroy_node(tp%root)

    deallocate (tp%ind)

  contains

    recursive subroutine destroy_node(np)
      type(tree_node), pointer :: np

      if (associated(np%left)) then
        call destroy_node(np%left)
        nullify (np%left)
      end if
      if (associated(np%right)) then
        call destroy_node(np%right)
        nullify (np%right)
      end if
      if (allocated(np%box)) deallocate (np%box)
      deallocate (np)
    end subroutine destroy_node

  end subroutine kdtree2_destroy

  ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
  ! returning their indexes and distances in 'indexes' and 'distances'
  ! arrays already allocated passed to this subroutine.
  subroutine kdtree2_n_nearest(tp, qv, nn, results)
    type(kdtree2), intent(in)                   :: tp
    real(kdkind), intent(In)                    :: qv(:)
    integer, intent(In)                         :: nn
    type(kdtree2_result), intent(inout), target :: results(nn)
    type(tree_search_record)                    :: sr

    sr%ballsize = huge(1.0)
    sr%qv = qv
    sr%nn = nn
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false.
    sr%pq = pq_create()

    call search(tp, sr, tp%root, nn, results)
    if (tp%sort) call kdtree2_sort_results(nn, results)
  end subroutine kdtree2_n_nearest

  ! Find the 'nn' vectors in the tree nearest to point 'idxin',
  ! with correlation window 'correltime', returing results in
  ! results(:), which must be pre-allocated upon entry.
  subroutine kdtree2_n_nearest_around_point(tp, idxin, correltime, nn, results)
    type(kdtree2), intent(in)                   :: tp
    integer, intent(In)                         :: idxin, correltime, nn
    type(kdtree2_result), intent(inout), target :: results(nn)
    type(tree_search_record)                    :: sr

    allocate (sr%qv(tp%dimen))
    sr%qv = tp%input_data(:, idxin) ! copy the vector
    sr%ballsize = huge(1.0)       ! the largest real(kdkind) number
    sr%centeridx = idxin
    sr%correltime = correltime
    sr%nn = nn
    sr%nfound = 0
    sr%pq = pq_create()

    call search(tp, sr, tp%root, nn, results)
    if (tp%sort) call kdtree2_sort_results(nn, results)
  end subroutine kdtree2_n_nearest_around_point

  ! Find the nearest neighbors to point 'idxin', within SQUARED
  ! Euclidean distance 'r2'.   Upon ENTRY, nalloc must be the
  ! size of memory allocated for results(1:nalloc).  Upon
  ! EXIT, nfound is the number actually found within the ball.
  !
  !  Note that if nfound .gt. nalloc then more neighbors were found
  !  than there were storage to store.  The resulting list is NOT
  !  the smallest ball inside norm r^2
  !
  ! Results are NOT sorted unless tree was created with sort option.
  subroutine kdtree2_r_nearest(tp, qv, r2, nfound, nalloc, results)
    type(kdtree2), intent(in)                   :: tp
    real(kdkind), intent(In)                    :: qv(:)
    real(kdkind), intent(in)                    :: r2
    integer, intent(out)                        :: nfound
    integer, intent(In)                         :: nalloc
    type(kdtree2_result), intent(inout), target :: results(nalloc)
    type(tree_search_record)                    :: sr

    sr%qv = qv
    sr%ballsize = r2
    sr%nn = 0      ! flag for fixed ball search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false.

    call search(tp, sr, tp%root, nalloc, results)
    nfound = sr%nfound
    if (tp%sort) call kdtree2_sort_results(nfound, results)

    if (sr%overflow .and. tp%verbose) then
      write (error_unit, '(A)') 'kdtree2_r_nearest warning: nfound > nalloc'
      write (error_unit, '(A)') 'Answer is not smallest ball (thus wrong)'
    end if

  end subroutine kdtree2_r_nearest

  ! Like kdtree2_r_nearest, but around a point 'idxin' already existing
  ! in the data set.
  !
  ! Results are NOT sorted unless tree was created with sort option.
  subroutine kdtree2_r_nearest_around_point(tp, idxin, correltime, r2, &
                                            nfound, nalloc, results)
    type(kdtree2), intent(in)                   :: tp
    integer, intent(In)                         :: idxin, correltime, nalloc
    real(kdkind), intent(in)                    :: r2
    integer, intent(out)                        :: nfound
    type(kdtree2_result), intent(inout), target :: results(nalloc)
    type(tree_search_record)                    :: sr

    allocate (sr%qv(tp%dimen))
    sr%qv = tp%input_data(:, idxin) ! copy the vector
    sr%ballsize = r2
    sr%nn = 0    ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    sr%overflow = .false.

    call search(tp, sr, tp%root, nalloc, results)
    nfound = sr%nfound
    if (tp%sort) call kdtree2_sort_results(nfound, results)

    if (sr%overflow .and. tp%verbose) then
      write (error_unit, '(A)') &
           'kdtree2_r_nearest_around_point warning: nfound > nalloc'
      write (error_unit, '(A)') 'Answer is not smallest ball (thus wrong)'
    end if

  end subroutine kdtree2_r_nearest_around_point

  ! Count the number of neighbors within square distance 'r2'.
  function kdtree2_r_count(tp, qv, r2) result(nfound)
    type(kdtree2), intent(in) :: tp
    real(kdkind), intent(In)  :: qv(:)
    real(kdkind), intent(in)  :: r2
    integer                   :: nfound
    type(tree_search_record)  :: sr
    type(kdtree2_result)      :: dummy_results(0)

    sr%qv       = qv
    sr%ballsize = r2
    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false.

    call search(tp, sr, tp%root, 0, dummy_results)
    nfound = sr%nfound

  end function kdtree2_r_count

  ! Count the number of neighbors within square distance 'r2' around
  ! point 'idxin' with decorrelation time 'correltime'.
  function kdtree2_r_count_around_point(tp, idxin, correltime, r2) &
    result(nfound)
    type(kdtree2), intent(in) :: tp
    integer, intent(In)       :: correltime, idxin
    real(kdkind), intent(in)  :: r2
    integer                   :: nfound
    type(tree_search_record)  :: sr
    type(kdtree2_result)      :: dummy_results(0)

    allocate (sr%qv(tp%dimen))
    sr%qv = tp%input_data(:, idxin)
    sr%ballsize = r2
    sr%nn = 0       ! flag for fixed r search
    sr%nfound = 0
    sr%centeridx = idxin
    sr%correltime = correltime
    sr%overflow = .false.

    call search(tp, sr, tp%root, 0, dummy_results)
    nfound = sr%nfound

  end function kdtree2_r_count_around_point

  ! Distance between iv[1:n] and qv[1:n]
  pure function square_distance(d, iv, qv) result(res)
    real(kdkind)             :: res
    integer, intent(in)      :: d
    real(kdkind), intent(in) :: iv(:), qv(:)

    res = sum((iv(1:d) - qv(1:d))**2)
  end function square_distance

  ! This is the innermost core routine of the kd-tree search.  Along
  ! with "process_terminal_node", it is the performance bottleneck.
  !
  ! This version uses a logically complete secondary search of
  ! "box in bounds", whether the sear
  recursive subroutine search(tp, sr, node, n_max, results)
    type(kdtree2), intent(in)               :: tp
    type(tree_search_record), intent(inout) :: sr
    type(Tree_node), intent(in)             :: node
    integer, intent(in)                     :: n_max
    type(kdtree2_result), intent(inout)     :: results(n_max)
    type(tree_node), pointer                :: ncloser, nfarther
    integer                                 :: cut_dim, i
    real(kdkind)                            :: qval, dis
    real(kdkind)                            :: ballsize

    if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
      ! we are on a terminal node
      if (sr%nn .eq. 0) then
        call process_terminal_node_fixedball(tp, sr, node, &
             n_max, results)
      else
        call process_terminal_node(tp, sr, node, n_max, results)
      end if
    else
      ! we are not on a terminal node
      cut_dim = node%cut_dim
      qval = sr%qv(cut_dim)

      if (qval < node%cut_val) then
        ncloser => node%left
        nfarther => node%right
        dis = (node%cut_val_right - qval)**2
      else
        ncloser => node%right
        nfarther => node%left
        dis = (node%cut_val_left - qval)**2
!          extra = qval- node%cut_val_left
      end if

      if (associated(ncloser)) call search(tp, sr, ncloser, n_max, results)

      ! we may need to search the second node.
      if (associated(nfarther)) then
        ballsize = sr%ballsize
        if (dis <= ballsize) then
          !
          ! we do this separately as going on the first cut dimen is often
          ! a good idea.
          ! note that if extra**2 < sr%ballsize, then the next
          ! check will also be false.
          !
          do i = 1, tp%dimen
            if (i .ne. cut_dim) then
              dis = dis + dis2_from_bnd(sr%qv(i), node%box(i)%lower, node%box(i)%upper)
              if (dis > ballsize) then
                return
              end if
            end if
          end do

          !
          ! if we are still here then we need to search mroe.
          !
          call search(tp, sr, nfarther, n_max, results)
        end if
      end if
    end if
  end subroutine search

  pure function dis2_from_bnd(x, amin, amax) result(res)
    real(kdkind), intent(in) :: x, amin, amax
    real(kdkind) :: res

    if (x > amax) then
      res = (x - amax)**2
    else
      if (x < amin) then
        res = (amin - x)**2
      else
        res = 0.0_kdkind
      end if
    end if
  end function dis2_from_bnd

  ! Look for actual near neighbors in 'node', and update
  ! the search results on the sr data structure.
  subroutine process_terminal_node(tp, sr, node, n_max, results)
    type(kdtree2), intent(in)               :: tp
    type(tree_search_record), intent(inout) :: sr
    type(tree_node), intent(in)             :: node
    integer, intent(in)                     :: n_max
    type(kdtree2_result), intent(inout)     :: results(n_max)
    integer                                 :: i, indexofi, k
    real(kdkind)                            :: sd, newpri

    mainloop: do i = node%l, node%u
      if (tp%rearrange) then
        sd = 0.0
        do k = 1, tp%dimen
          sd = sd + (tp%rearranged_data(k, i) - sr%qv(k))**2
          if (sd > sr%ballsize) cycle mainloop
        end do
        indexofi = tp%ind(i)  ! only read it if we have not broken out
      else
        indexofi = tp%ind(i)
        sd = 0.0
        do k = 1, tp%dimen
          sd = sd + (tp%input_data(k, indexofi) - sr%qv(k))**2
          if (sd > sr%ballsize) cycle mainloop
        end do
      end if

      if (sr%centeridx > 0) then ! doing correlation interval?
        if (abs(indexofi - sr%centeridx) < sr%correltime) cycle mainloop
      end if

      !
      ! two choices for any point.  The list so far is either undersized,
      ! or it is not.
      !
      ! If it is undersized, then add the point and its distance
      ! unconditionally.  If the point added fills up the working
      ! list then set the sr%ballsize, maximum distance bound (largest distance on
      ! list) to be that distance, instead of the initialized +infinity.
      !
      ! If the running list is full size, then compute the
      ! distance but break out immediately if it is larger
      ! than sr%ballsize, "best squared distance" (of the largest element),
      ! as it cannot be a good neighbor.
      !
      ! Once computed, compare to best_square distance.
      ! if it is smaller, then delete the previous largest
      ! element and add the new one.

      if (sr%nfound .lt. sr%nn) then
        !
        ! add this point unconditionally to fill list.
        !
        sr%nfound = sr%nfound + 1
        newpri = pq_insert(sr%pq, results, sd, indexofi)
        if (sr%nfound .eq. sr%nn) sr%ballsize = newpri
        ! we have just filled the working list.
        ! put the best square distance to the maximum value
        ! on the list, which is extractable from the PQ.

      else
        !
        ! now, if we get here,
        ! we know that the current node has a squared
        ! distance smaller than the largest one on the list, and
        ! belongs on the list.
        ! Hence we replace that with the current one.
        !
        sr%ballsize = pq_replace_max(sr%pq, results, sd, indexofi)
      end if
    end do mainloop

  end subroutine process_terminal_node

  ! Look for actual near neighbors in 'node', and update
  ! the search results on the sr data structure, i.e.
  ! save all within a fixed ball.
  subroutine process_terminal_node_fixedball(tp, sr, node, n_max, results)
    type(kdtree2), intent(in)               :: tp
    type(tree_search_record), intent(inout) :: sr
    type(tree_node), intent(in)             :: node
    integer, intent(in)                     :: n_max
    type(kdtree2_result), intent(inout)     :: results(n_max)
    integer                                 :: i, indexofi, k
    real(kdkind)                            :: sd

    ! search through terminal bucket.
    mainloop: do i = node%l, node%u

      !
      ! two choices for any point.  The list so far is either undersized,
      ! or it is not.
      !
      ! If it is undersized, then add the point and its distance
      ! unconditionally.  If the point added fills up the working
      ! list then set the sr%ballsize, maximum distance bound (largest distance on
      ! list) to be that distance, instead of the initialized +infinity.
      !
      ! If the running list is full size, then compute the
      ! distance but break out immediately if it is larger
      ! than sr%ballsize, "best squared distance" (of the largest element),
      ! as it cannot be a good neighbor.
      !
      ! Once computed, compare to best_square distance.
      ! if it is smaller, then delete the previous largest
      ! element and add the new one.

      ! which index to the point do we use?

      if (tp%rearrange) then
        sd = 0.0
        do k = 1, tp%dimen
          sd = sd + (tp%rearranged_data(k, i) - sr%qv(k))**2
          if (sd > sr%ballsize) cycle mainloop
        end do
        indexofi = tp%ind(i)  ! only read it if we have not broken out
      else
        indexofi = tp%ind(i)
        sd = 0.0
        do k = 1, tp%dimen
          sd = sd + (tp%input_data(k, indexofi) - sr%qv(k))**2
          if (sd > sr%ballsize) cycle mainloop
        end do
      end if

      if (sr%centeridx > 0) then ! doing correlation interval?
        if (abs(indexofi - sr%centeridx) < sr%correltime) cycle mainloop
      end if

      sr%nfound = sr%nfound + 1
      if (sr%nfound .gt. n_max) then
        ! oh nuts, we have to add another one to the tree but
        ! there isn't enough room.
        sr%overflow = .true.
      else
        results(sr%nfound)%dis = sd
        results(sr%nfound)%idx = indexofi
      end if
    end do mainloop

  end subroutine process_terminal_node_fixedball

  ! Find the 'n' nearest neighbors to 'qv' by exhaustive search. only use this
  ! subroutine for testing, as it is SLOW! The whole point of a k-d tree is to
  ! avoid doing what this subroutine does.
  subroutine kdtree2_n_nearest_brute_force(tp, qv, nn, results)
    type(kdtree2), intent(in) :: tp
    real(kdkind), intent(In)  :: qv(:)
    integer, intent(In)       :: nn
    type(kdtree2_result)      :: results(nn)
    integer                   :: i, j, k
    real(kdkind), allocatable :: all_distances(:)

    allocate (all_distances(tp%n))
    do i = 1, tp%n
      all_distances(i) = square_distance(tp%dimen, qv, tp%input_data(:, i))
    end do
    ! now find 'n' smallest distances
    do i = 1, nn
      results(i)%dis = huge(1.0)
      results(i)%idx = -1
    end do
    do i = 1, tp%n
      if (all_distances(i) < results(nn)%dis) then
        ! insert it somewhere on the list
        do j = 1, nn
          if (all_distances(i) < results(j)%dis) exit
        end do
        ! now we know 'j'
        do k = nn - 1, j, -1
          results(k + 1) = results(k)
        end do
        results(j)%dis = all_distances(i)
        results(j)%idx = i
      end if
    end do
    deallocate (all_distances)
  end subroutine kdtree2_n_nearest_brute_force

  ! find the nearest neighbors to 'qv' with distance**2 <= r2 by exhaustive
  ! search. only use this subroutine for testing, as it is SLOW! The whole
  ! point of a k-d tree is to avoid doing what this subroutine does.
  subroutine kdtree2_r_nearest_brute_force(tp, qv, r2, nfound, results)
    type(kdtree2), intent(in) :: tp
    real(kdkind), intent(In)  :: qv(:)
    real(kdkind), intent(In)  :: r2
    integer, intent(out)      :: nfound
    type(kdtree2_result)      :: results(:)
    integer                   :: i, nalloc
    real(kdkind), allocatable :: all_distances(:)

    allocate (all_distances(tp%n))
    do i = 1, tp%n
      all_distances(i) = square_distance(tp%dimen, qv, tp%input_data(:, i))
    end do

    nfound = 0
    nalloc = size(results, 1)

    do i = 1, tp%n
      if (all_distances(i) < r2) then
        ! insert it somewhere on the list
        if (nfound .lt. nalloc) then
          nfound = nfound + 1
          results(nfound)%dis = all_distances(i)
          results(nfound)%idx = i
        end if
      end if
    end do
    deallocate (all_distances)

    call kdtree2_sort_results(nfound, results)

  end subroutine kdtree2_r_nearest_brute_force

  ! Use after search to sort results(1:nfound) in order of increasing
  ! distance.
  subroutine kdtree2_sort_results(nfound, results)
    integer, intent(in)  :: nfound
    type(kdtree2_result) :: results(nfound)

    if (nfound .gt. 1) call heapsort_struct(results, nfound)
  end subroutine kdtree2_sort_results

  ! Sort a(1:n) in ascending order
  subroutine heapsort_struct(a, n)
    integer, intent(in)                 :: n
    type(kdtree2_result), intent(inout) :: a(:)
    type(kdtree2_result)                :: tmpval ! temporary value

    integer :: i, j
    integer :: ileft, iright

    ileft = n/2 + 1
    iright = n

    if (n .eq. 1) return

    do
      if (ileft > 1) then
        ileft = ileft - 1
        tmpval = a(ileft)
      else
        tmpval = a(iright)
        a(iright) = a(1)
        iright = iright - 1
        if (iright == 1) then
          a(1) = tmpval
          return
        end if
      end if
      i = ileft
      j = 2*ileft
      do while (j <= iright)
        if (j < iright) then
          if (a(j)%dis < a(j + 1)%dis) j = j + 1
        end if
        if (tmpval%dis < a(j)%dis) then
          a(i) = a(j);
          i = j
          j = j + j
        else
          j = iright + 1
        end if
      end do
      a(i) = tmpval
    end do
  end subroutine heapsort_struct

  ! Create an empty priority queue
  function pq_create() result(res)
    type(pq) :: res
    res%heap_size = 0
  end function pq_create

  ! Insert a new element and return the new maximum priority, which may or may
  ! not be the same as the old maximum priority.
  real(kdkind) function pq_insert(a, elems, dis, idx)
    type(pq), intent(inout)             :: a
    type(kdtree2_result), intent(inout) :: elems(*)
    real(kdkind), intent(in)            :: dis
    integer, intent(in)                 :: idx
    integer                             :: i, isparent
    real(kdkind)                        :: parentdis

    a%heap_size = a%heap_size + 1
    i = a%heap_size

    do while (i .gt. 1)
      isparent = int(i/2)
      parentdis = elems(isparent)%dis
      if (dis .gt. parentdis) then
        ! move what was in i's parent into i.
        elems(i)%dis = parentdis
        elems(i)%idx = elems(isparent)%idx
        i = isparent
      else
        exit
      end if
    end do

    ! insert the element at the determined position
    elems(i)%dis = dis
    elems(i)%idx = idx

    pq_insert = elems(1)%dis

  end function pq_insert

  ! Replace the extant maximum priority element in the PQ with (dis,idx).
  ! Return the new maximum priority, which may be larger or smaller than the
  ! old one.
  real(kdkind) function pq_replace_max(a, elems, dis, idx)
    type(pq), intent(inout)             :: a
    type(kdtree2_result), intent(inout) :: elems(*)
    real(kdkind), intent(in)            :: dis
    integer, intent(in)                 :: idx
    integer                             :: parent, child, N
    real(kdkind)                        :: prichild, prichildp1

    N = a%heap_size
    if (N .ge. 1) then
      parent = 1
      child = 2

      loop: do while (child .le. N)
        prichild = elems(child)%dis

        !
        ! posibly child+1 has higher priority, and if
        ! so, get it, and increment child.
        !

        if (child .lt. N) then
          prichildp1 = elems(child + 1)%dis
          if (prichild .lt. prichildp1) then
            child = child + 1
            prichild = prichildp1
          end if
        end if

        if (dis .ge. prichild) then
          exit loop
          ! we have a proper place for our new element,
          ! bigger than either children's priority.
        else
          ! move child into parent.
          elems(parent) = elems(child)
          parent = child
          child = 2*parent
        end if
      end do loop
      elems(parent)%dis = dis
      elems(parent)%idx = idx
      pq_replace_max = elems(1)%dis
    else
      elems(1)%dis = dis
      elems(1)%idx = idx
      pq_replace_max = dis
    end if

  end function pq_replace_max

end module kdtree2_module
