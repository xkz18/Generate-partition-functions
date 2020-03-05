program bipartite_deriv
implicit none
    ! INPUT:
    ! n1=Nsim(iset,1),n2=Nsim(iset,2)
    ! ntot=Nsim(iset,1)+Nsim(iset,2)
    ! mat_b=qmatrix
    ! facinv(0:ntot)
    ! Note: In upper level, initialize qmatrix(0:n1+1,0:n2+1)=0
    ! and then give qmatrix input values.
    ! (+1 is to protect code from overflow of data.)
    !callbipartite_deriv(n1,n2,ntot,facinv(0:ntot),qmatrix(0:n1+1,0:n2+1),nk(0:n1+1,0:n2+1))
    integer:: n1,n2,ntot!,Nmax?
    real*8,dimension(0:40) :: facinv
real*8,dimension(0:40) :: factorial
    real*8,dimension(0:21,0:21):: mat_b
    real*8,dimension(0:21,0:21):: nk
    real*8 :: Qn
    integer :: i,j,k1,k2,p,q
    real*8,dimension(0:21,0:21) :: mat_a,mat_g,mat_k,mat_c
    real*8,dimension(0:21,0:21) :: mat_x,mat_y,mat_m
    real*8,dimension(0:21,0:21) :: mat_o,C
    real*8,dimension(0:21,0:21) :: deriv_b,deriv_d,deriv_x
facinv =0
facinv(0)=1.0d0
factorial(0)=1.0d0
  do i=1,40
    factorial(i)=factorial(i-1)*i
    facinv(i)=1.0d0/factorial(i)
  enddo
    C=0
    Qn=0
    nk=0
    mat_a=0
    mat_c=0
    mat_g=0
    mat_k=0
    mat_x=0
    mat_y=0
    mat_m=0
    mat_o=0
    deriv_x=0
    deriv_b=0
    deriv_d=0

    n1 = 20
    n2 = 20
    mat_b = 0

    do i=0,n1

        read(*,*) (mat_b(i,j),j=0,n2)
 !   write(*,*) i,j,mat_b(i,j)

enddo
!write(*,*) mat_b

    do i=0,n1
        do j=0,n2
            deriv_b(i,j)=mat_b(i+1,j)*(i+1)
!write(*,*) "mat_b",i,j,mat_b(i,j)
!write(*,*) "deriv_b_test",i,j,deriv_b(i,j)
 
       enddo
    enddo

    do i=0,n1
        do j=0,n2
            deriv_d(i,j)=mat_b(i,j+1)*(j+1)
        enddo
    enddo
 !write(*,*) deriv_b
    ! Extract Q(N1,N2,V,T) from Xi(lambda1,lambda2,V,T)


    mat_a(0,0)=mat_b(0,0)
    do j=0,n2
        deriv_x(0,j)=mat_b(0,j+1)*(j+1)
    enddo
    mat_a(0,1)=deriv_x(0,0)

    mat_x=0
    mat_x=deriv_x
    do j=2,n2
        call calc2(n1,n2,deriv_d(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1),mat_y(0:n1+1,0:n2+1))
        mat_a(0,j)=mat_y(0,0)*facinv(j)
        mat_x=0
        mat_x=mat_y
    enddo
    !write(*,*) mat_x
    mat_m=0
    do j=0,n2
        do i=0,n1
            mat_m(i,j)=mat_b(i+1,j)*(i+1)
        enddo
    enddo
    mat_a(1,0)=mat_m(0,0)
    mat_x=0
    mat_x=mat_m
    do j=1,n2
        call calc2(n1,n2,deriv_d(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1),mat_y(0:n1+1,0:n2+1))
        mat_a(1,j)=mat_y(0,0)*facinv(j)
        mat_x=0
        mat_x=mat_y
    enddo
    !write(*,*) mat_y
    do i=2,n1

       mat_x=0

        !write(*,*) mat_m
        call calc1(n1,n2,deriv_b(0:n1+1,0:n2+1),mat_m(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1))
        !if (i .EQ. 2) then
	!write(*,*) "mat_x",mat_x
        !write(*,*) "deriv_b",deriv_b
        !end if
        mat_a(i,0)=mat_x(0,0)*facinv(i) !! Corrected error here 2017.02.10
        mat_c=mat_x
        do j=1,n2
            mat_y=0
            call calc2(n1,n2,deriv_d(0:n1+1,0:n2+1),mat_x(0:n1+1,0:n2+1),mat_y(0:n1+1,0:n2+1))
            do p=0,n1
              do q=0,n2
if (mat_y(p,q) .NE.0)	then	
write(*,*) "i",i,"j",j,p,q,mat_y(p,q)
endif           
 enddo
enddo
            mat_a(i,j)=mat_y(0,0)*facinv(i)*facinv(j)
            mat_x=0
            mat_x=mat_y
        enddo
        mat_m=mat_c
        mat_c=0
    enddo
!write(*,*) mat_y
    Qn=mat_y(0,0)*facinv(n1)*facinv(n2)
    write (*,*) Qn
    !Calculate <n_i,j>
    do i=0,n1
        do j=0,n2
            k1=n1-i
            k2=n2-j
            C(k1,k2)=mat_b(k1,k2)*mat_a(i,j)
            nk(k1,k2)=C(k1,k2)/Qn
        enddo
    enddo
end


subroutine calc1(n1,n2,deriv_b,mat_m,mat_x)
    implicit none
    integer,intent(in):: n1,n2
    real*8,dimension(0:n1+1,0:n2+1),intent(in) :: mat_m,deriv_b
    real*8,dimension(0:n1+1,0:n2+1),intent(out)::mat_x
    integer :: r,s,m,j,p,l,k,i
    real*8,dimension(0:n1+1,0:n2+1) :: deriv_u,mat_u
    real*8,dimension(0:n1*2,0:n2*2) :: mat_v
    real*8,dimension(0:n1+1,0:n2+1) :: mat_w

    mat_x=0
    mat_u=mat_m
    deriv_u=0
    mat_v=0
    mat_w=0
    r=0
    s=0
i = 0
    do m=0,n2
       do j=0,n1
i = i+1
            deriv_u(j,m)=mat_u(j+1,m)*(j+1)
        enddo
    enddo
!write(*,*) i,"mat_u",mat_u
   do j=0,n1
       do p=0,n2
!write(*,*) "deriv_b",j,p,deriv_b(j,p)
enddo
enddo
   do j=0,n1
        do p=0,n1
            do l=0,n2
                do k=0,n2
                    r=j+p
                    s=l+k
                    mat_v(r,s)=mat_v(r,s)+mat_u(j,l)*deriv_b(p,k)
!               if ((r .EQ. 2) .AND. (s .EQ.2)) then
!write(*,*) r,s,j,p,l,k,mat_v(r,s),mat_u(j,l),deriv_b(p,k)
!end if
!if ((r .EQ. 1) .AND. (s .EQ.1)) then
!write(*,*) r,s,j,p,l,k,mat_v(r,s),mat_u(j,l),deriv_b(p,k)
!end if
               enddo
            enddo
        enddo
    enddo
!write(*,*) mat_v
   do l=0,n1
        do m=0,n2
            mat_w(l,m)=deriv_u(l,m)+mat_v(l,m)
 	    if (mat_w(l,m) .NE. 0) then
!write(*,*) "l,m,mat_w",l,m,mat_w(l,m)
end if
       enddo
    enddo
    mat_x=mat_w
write(*,*) "mat_w",mat_w(3,:)
end

subroutine calc2(n1,n2,deriv_d,mat_n,mat_y)
    implicit none
    integer,intent(in):: n1,n2
    real*8,dimension(0:n1+1,0:n2+1),intent(in) :: deriv_d,mat_n
    real*8,dimension(0:n1+1,0:n2+1),intent(out)::mat_y
    integer :: s,j,l,k
    real*8,dimension(0:n1+1,0:n2+1) :: deriv_o
    real*8,dimension(0:n1*2,0:n2*2) :: mat_p
    real*8,dimension(0:n1+1,0:n2+1) :: mat_q,mat_o
    mat_y=0
    mat_o=mat_n
    deriv_o=0
    mat_p=0
    mat_q=0
    s=0

    do j=0,n2
        deriv_o(0,j)=mat_o(0,j+1)*(j+1)
    enddo
    do l=0,n2
        do k=0,n2
            s=l+k
            mat_p(0,s)=mat_p(0,s)+mat_o(0,l)*deriv_d(0,k)
        enddo
    enddo
    do l=0,n2
        mat_q(0,l)=deriv_o(0,l)+mat_p(0,l)
    enddo
    mat_y=mat_q
end
