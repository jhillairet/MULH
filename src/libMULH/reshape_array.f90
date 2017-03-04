module reshape_array

contains
subroutine add21Darray(A,extra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lengthens rank-one array A by extra
!
! ARGUMENTS
!  A : real allocatable rank-one array
! extra : integer indicating how many more spaces are required

use precision_def
! Subroutine arguments
real(long), allocatable, target, intent(inout) :: A(:)
integer, intent(in) :: extra
! Dummy stuff
real(long), allocatable, target :: temp (:)
real(long), pointer :: p(:)
integer :: lb,i
  

p => A  
lb = size(A)
  
allocate (temp(lb))  
temp = p 
deallocate(A)  
allocate (A(lb+extra))  

  
A(1:lb) = temp(:)

do i = lb+1,lb+extra
  A(i) = 0.
enddo

deallocate(temp) 


end subroutine add21Darray

subroutine add21DarrayI(A,extra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lengthens rank-one array A by extra
!
! ARGUMENTS
!  A : real allocatable rank-one array
! extra : integer indicating how many more spaces are required

! Subroutine arguments
integer, allocatable, target, intent(inout) :: A(:)
integer, intent(in) :: extra
! Dummy stuff
integer, allocatable, target :: temp (:)
integer, pointer :: p(:)
integer :: lb,i
  

p => A  
lb = size(A)
  
allocate (temp(lb))  
temp = p 
deallocate(A)  
allocate (A(lb+extra))  

  
A(1:lb) = temp(:)

do i = lb+1,lb+extra
  A(i) = 0
enddo

deallocate(temp) 

end subroutine add21DarrayI

subroutine add23Darray(A,extra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lengthens rank-3 array A by extra
!
! ARGUMENTS
!  A : real allocatable rank-three array
! extra : integer array indicating how many more spaces are required
!	  in each dimension


!!!!!! Need to change this to make it truly general. Right now only accepts adding 1 at a time along 1 dimension

use precision_def
! Subroutine arguments
real(long), allocatable, intent(inout) :: A(:,:,:)
integer, intent(in) :: extra(3)
! Dummy stuff
real(long), allocatable :: temp (:,:,:)
integer :: i,j,k
integer :: lb(3)
  

lb = shape(A)
  
allocate (temp(lb(1),lb(2),lb(3)))  
temp = A 
deallocate(A)  
allocate (A(lb(1)+extra(1),lb(2)+extra(2),lb(3)+extra(3)))  

A(1:lb(1),1:lb(2),1:lb(3)) = temp(:,:,:)

!if (extra(1)/=0 .AND. extra(2)/=0 .AND. extra(3)/=0) then
!  do i = lb(1)+1,lb(1)+extra
!    do j = lb(2)+1,lb(2)+extra
!      do k = lb(3)+1,lb(3)+extra

!	A(i,j,k) = 0.

!      enddo
!    enddo
!  enddo
!elseif (extra(1)/=0 .AND. extra(2)==0 .AND. extra(3)==0) then
  do i = 1,lb(1)
    do j = 1,lb(2)
      do k = lb(3)+1,lb(3)+extra(3)
	A(i,j,k) = 0.
      enddo
    enddo
  enddo
!i = lb(1)
!do
!  if (extra(1) /= 0) i = i + 1
!  if (extra(1)==0) exit
!  j = lb(2)
!  do
!    if (extra(2) /= 0) j = j + 1
!    if (extra(2)==0) exit
!    k = lb(3)
!    do
!      if (extra(3) /= 0) k = k + 1
!      if (extra(3)==0) exit

!      A(i,j,k) = 0.

!      if (k == lb(3)+extra(3)) exit
!    enddo
!    if (j == lb(2)+extra(2)) exit
!  enddo
!  if (i == lb(1)+extra(1)) exit
!enddo

deallocate(temp) 

end subroutine add23Darray

subroutine add22Darray(A,extra)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lengthens rank-2 array A by extra
!
! ARGUMENTS
!  A : real allocatable rank-three array
! extra : integer array indicating how many more spaces are required
!	  in each dimension


!!!!!! Need to change this to make it truly general. Right now only accepts adding 1 at a time along 1 dimension

use precision_def
! Subroutine arguments
real(long), allocatable, intent(inout) :: A(:,:)
integer, intent(in) :: extra(2)
! Dummy stuff
real(long), allocatable :: temp (:,:)
integer :: i,j
integer :: lb(2)
  

lb = shape(A)
  
allocate (temp(lb(1),lb(2)))  
temp = A 

deallocate(A)  
allocate (A(lb(1)+extra(1),lb(2)+extra(2)))  

A(1:lb(1),1:lb(2)) = temp(:,:)

do i = lb(1)+1,lb(1)+extra(1)
  do j = 1,lb(2)
    A(i,j) = 0.
  enddo
enddo

deallocate(temp) 

end subroutine add22Darray

subroutine removeI(A,ir)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Removes item ir from allocatable integer array A and
! resizes array A
!
! ARGUMENTS
! A: integer allocatable rank-one array
! ir: index of the item to be removed

! Subroutine arguments
integer, allocatable, intent(inout) :: A(:)
integer, intent(in) :: ir
! Dummy stuff
integer, allocatable, target :: temp (:)
integer :: lb,i,j
  

lb = size(A)
  
allocate (temp(lb-1))

if (ir == 1) then
  
  do i = 2,lb
    temp(i-1) = A(i)
  enddo

else
  
  j = 1
  do i = 1,lb
    if (i /= ir) then
      temp(j) = A(i)
      j = j + 1
    endif
  enddo

endif

deallocate(A)  
allocate (A(lb-1))

A = temp

deallocate(temp)

end subroutine removeI

real(long) function mean(A)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!		MEAN
!
! Find the arithmetic mean of elements in array A
!

use precision_def
real(long) :: A(:)
integer :: i

mean = 0.
do i = 1,size(A)
  mean = mean + A(i)
enddo
mean = mean/size(A)

return
end function mean

subroutine file_row_count ( file_name, line_num )

!*****************************************************************************80
!
!! FILE_ROW_COUNT counts the number of rows in a file.
!
!  Discussion:
!
!    The file is assumed to be a simple text file.
!
!    Blank lines and comment lines, which begin with '#', are not counted.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) FILE_NAME, the name of the file.
!
!    Output, integer ( kind = 4 ) LINE_NUM, the number of lines found in the 
!    file.  If the file could not be opened, then LINE_NUM is returned as -1.
!
  implicit none

  character ( len = * ) file_name
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  character ( len = 255 ) line
  integer ( kind = 4 ) line_num
  logical, parameter :: verbose = .false.

  line_num = 0
!
!  Open the file.
!
  call get_unit ( iunit )

  open ( unit = iunit, file = file_name, status = 'old', form = 'formatted', &
    access = 'sequential', iostat = ios )

  if ( ios /= 0 ) then

    line_num = -1

    if ( verbose ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file:'
      write ( *, '(4x,a)' ) '"' // trim ( file_name ) // '".'
    end if

    return

  end if
!
!  Count the lines.
!
  do

    read ( iunit, '(a)', iostat = ios ) line

    if ( ios /= 0 ) then
      exit
    end if

    if ( len_trim ( line ) == 0 ) then
      cycle
    end if

    if ( line(1:1) == '#' ) then
      cycle
    end if

    line_num = line_num + 1

  end do

  close ( unit = iunit )

  return
end subroutine file_row_count

subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen
 
  iunit = 0
 
  do i = 1, 99
 
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )
 
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if
 
  end do

  return
end subroutine get_unit

end module reshape_array
