program gs

 implicit none

	REAL, Dimension(3,3) :: a
	REAL, Dimension(3) :: b, x_old,x_ref, x, erra
        INTEGER :: i,j,k,maxiter,num
	real :: err
	real, parameter :: tol = 1E-6
! 	a= (/1, 2, 3, &
!		 4, 5, 6, &
!		 7, 8, 9 /)
        num = 1
        a = transpose(reshape((/ 10, 3, 1,&
                                 2,-10, 3,&
                                 1, 3, 10 /), shape(a)))

        b = reshape((/ 14, -5, 14 /), shape(b))

        x_old = reshape((/ 0, 0, 0 /), shape(x_old))
        x_ref = reshape((/ 0, 0, 0 /), shape(x_ref))

!	do i = 1, 3
!		do j = 1, 3				
!		end do
!		write(*,*), a(i,:)	
!	end do
! 	do while (err .ge. tol)
      maxiter = 25

  do k = 1,maxiter 
          !write(*,*) 'k,x_old',k, x_old(:)
        do i =  1, 3 
         do j = 1, 3
          !write(*,*) 'after j', x_old
          if (j .ne. i) then
           x(i) = a(i,j)*x_old(j) + x_ref(i)
           x_ref(i) = x(i)
          end if!else
!          write(*,*) 'i,j,x_ref', i,j,x_ref
         end do 
           x(i) = (b(i) - x_ref(i)) / a(i,i)! - a(i,j)*x(j))
           !x(i) = ((b(i) - a(i,i+1)*x_old(i+1) - a(i,i+2)*x_old(i+2))) / a(i,i)! - a(i,j)*x(j))
        
         !end do
         erra(i) = x(i) - x_old(i)
         err = abs(erra(i))
         if (err <= tol) then
          write(*,*) 'iters,x,err:',k, x_old(:),err
          exit !stop
         end if 
! 		write(*,*), err
!         if (num==1)
 !x_old = x    ! if uncommented=point-Jacobi, else Gauss-Seidel

        end do
        
        x_old = x
        x_ref = reshape((/ 0, 0, 0 /), shape(x_ref))
         !write(*,*) 'x_old', x_old(:) 
      end do

      write(*,*) 'Computation complete...'
      write(*,*) 'Number of iters, error', k,err
 !h	write(*,*), err 	
!	if ( ios /= 0 ) stop "Write error in file "

!	if (allocated(a)) deallocate(a, stat=err)
!	if (err /= 0) print *, "a: Deallocation request denied"

end program gs
