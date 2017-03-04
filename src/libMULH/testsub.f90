module testsub

contains
subroutine mallocate(A)

implicit none
integer, allocatable :: A(:,:)

allocate (A(1,2))

A(1,1) = 13
A(1,2) = 27


end subroutine mallocate

end module testsub
