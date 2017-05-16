subroutine lgwt(N,a,b,gpt,gwt)
    implicit none
    integer(8)::N
    real(8)::a,b
    real(8)::gpt(100),gwt(100)
    !*****************************
    real(8)::L(100,100),Lp(100),y(100),y0(100),xu(100)
    integer(8)::N1,N2,i,j,k
    real(8)::pi,temp

    N=N-1
    N1=N+1
    N2=N+2

    do i=1,N1
        xu(i)=-1.d0+2.d0*(i-1)/(N1-1)
    enddo
    !
    pi=4.d0*datan(1.d0)

    do i=1,N1
        y(i)=dcos((2.d0*(i-1)+1.0)*pi/(2.d0*N+2.d0))+(0.27d0/N1)*dsin(pi*xu(i)*N/N2)
    enddo

    ! Legendre-Gauss Vandermonde Matrix

    do i=1,N1
        do j=1,N2
            L(i,j)=0.d0
        enddo
        Lp(i)=0.d0 !Derivative of LGVM
    enddo
    !
    ! Compute the zeros of the N+1 Legendre Polynomial
    ! using the recursion relation and the Newton-Raphson method

    do i=1,N1
        y0(i)=2.d0
    enddo

    ! Iterative until new points are uniformly within epsilon of old points

    do while (1)
        temp=dabs(y(1)-y0(1))
        do i=2,N1
            if(dabs(y(i)-y0(i))>temp) then
                temp=abs(y(i)-y0(i))
            endif
        enddo
        if(temp<1.e-15) then
            goto 10
        endif
        do i=1,N1
            L(i,1)=1.d0
            L(i,2)=y(i)
        enddo

        do i=1,N1
            do k=2,N1
                L(i,k+1)=((2.d0*k-1.d0)*y(i)*L(i,k)-(k-1)*L(i,k-1))/k
            enddo
        enddo

        do i=1,N1
            Lp(i)=N2*(L(i,N1)-y(i)*L(i,N2))/(1.d0-y(i)*y(i))
        enddo

        do i=1,N1
            y0(i)=y(i)
        enddo
        do i=1,N1
            y(i)=y0(i)-L(i,N2)/Lp(i)
        enddo
    enddo
10  temp=10.0
    ! Linear map from [-1,1] to [a,b]

    do i=1,N1
        gpt(i)=(a*(1.d0-y(i))+b*(1.d0+y(i)))/2.d0
    enddo

    ! Compute the weights

    do i=1,N1
        gwt(i)=(b-a)/((1.d0-y(i)**2)*(Lp(i)**2))*((1.d0*N2/N1)**2)
    enddo

    N=N1

    end