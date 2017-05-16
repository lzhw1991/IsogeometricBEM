Module Solver
    contains
    
    subroutine AGAUSS(A,B,N,X,L)
    !全主元Guass消去法求解方程AX=B
    !A(N,N),B(N),X(N),方程的解为X
    !如果L=0则方程奇异，求解失败
    implicit none
    integer(8),intent(in)::n
    real(8),intent(in out)::A(:,:),B(:)
    real(8),intent(out)::X(:)
    integer(8),intent(in out)::L
    !*********************
    integer(8)::i,j,k
    real(8)::d,t
    integer::js(n),is
    !*********************
    l=1
    do k=1,n-1
        d=0.d0
        do i=k,n
            do j=k,n
                if(abs(A(i,j))>d) then
                    d=abs(A(i,j))
                    js(k)=j
                    is=i
                endif
            enddo
        enddo
        if((d+1.d0)==1.d0) then
            l=0
        else
            if(js(k).ne.k) then
                do i=1,n
                    t=A(i,k)
                    A(i,k)=A(i,js(k))
                    A(i,js(k))=t
                enddo
            endif
            if(is.ne.k) then
                do j=k,n
                    t=A(k,j)
                    A(k,j)=A(is,j)
                    A(is,j)=t
                enddo
                t=B(k)
                B(k)=B(is)
                B(is)=t
            endif
        endif
        if(l==0) then
            write(*,100)
            return
        endif
        do j=k+1,n
            A(k,j)=A(k,j)/A(k,k)
        enddo
        B(k)=B(k)/A(k,k)
        do i=k+1,n
            do j=k+1,n
                A(i,j)=A(i,j)-A(i,k)*A(k,j)
            enddo
            B(i)=B(i)-A(i,k)*B(k)
        enddo
    enddo
    if((abs(A(n,n))+1.d0)==1.d0) then
        l=0
        write(*,100)
        return
    endif
    X(n)=B(n)/A(n,n)
    do i=n-1,1,-1
        t=0.d0
        do j=i+1,n
            t=t+A(i,j)*X(j)
        enddo
        X(i)=B(i)-t
    enddo
100 format(1X,'Fail!')
    js(N)=N
    do k=n,1,-1
        if(js(k).ne.k) then
            t=X(k)
            X(k)=X(js(K))
            X(js(k))=t
        endif
    enddo
    return
    end subroutine AGAUSS
    
    
    subroutine AGGJE(A,B,N,X,L)
    !全主元高斯-若当消去法解大型系数矩阵
    !AX=B，X(N)为解，A(N,N),B(N)
    !L=0则矩阵奇异，求解失败
    implicit none
    real(8),intent(in out)::A(:,:),X(:)
    real(8),intent(in)::B(:)
    integer(8),intent(in)::N
    integer(8),intent(in out)::L
    !*********************************
    integer(8)::i,j,k,is,js(N)
    real(8)::d,t
    X=B
    L=1
    do k=1,n
        d=0.d0
        do i=k,n
            do j=k,n
                if(abs(a(i,j))>d) then
                    d=abs(a(i,j))
                    js(k)=j
                    is=i
                endif
            enddo
        enddo
        if(d+1.d0==1.d0) then
            write(*,101)
            L=0
            return
        endif
101     format(1X,'Failed!')
        do j=k,n
            t=a(k,j)
            a(k,j)=a(is,j)
            a(is,j)=t
        enddo
        t=x(k)
        x(k)=x(is)
        x(is)=t
        do i=1,n
            t=a(i,k)
            a(i,k)=a(i,js(k))
            a(i,js(k))=t
        enddo
        t=a(k,k)
        do j=k+1,n
            if(a(k,j).ne.0.d0) then
                a(k,j)=a(k,j)/t
            endif
        enddo
        x(k)=x(k)/t
        do j=k+1,n
            if(a(k,j).ne.0.d0) then
                do i=1,n
                    if((i.ne.k).and.(a(i,k).ne.0.d0)) then
                        a(i,j)=a(i,j)-a(i,k)*a(k,j)
                    endif
                enddo
            endif
        enddo
        do i=1,n
            if((i.ne.k).and.(a(i,k).ne.0.d0)) then
                x(i)=x(i)-a(i,k)*x(k)
            endif
        enddo
    enddo
    do k=n,1,-1
        if(k.ne.js(k)) then
            t=x(k)
            x(k)=x(js(k))
            x(js(k))=t
        endif
    enddo
    return
    end subroutine AGGJE
    
    subroutine BEM_Solver(A,B,X,N,d)
    implicit none
    integer(8)::n
    real(8)::A(n,n),B(n),X(n),d
    !*********************
    integer(8)::i,j,k,k1,n1,L
    real(8)::c,tol
    !*********************
    tol=1.e-20

    n1=n-1
    do 100 k=1,n1
        k1=k+1
        c=a(k,k)
        if(Abs(c)-tol) 1,1,3
1       do 7 j=k1,n

            if(abs((A(j,k)))-tol) 7,7,5
5           do 6 L=k,n
                c=A(k,L)
                A(k,L)=A(j,L)
6           A(j,L)=c
            c=B(k)
            B(k)=B(j)
            B(j)=c
            c=A(k,k)
            goto 3
7       continue
        goto 8

3       c=A(k,k)
        do 4 j=k1,n
4       A(k,j)=A(k,j)/c
        B(k)=B(k)/c


        do 10 i=k1,n
            c=A(i,k)
            do 9 j=k1,n
9           A(i,j)=A(i,j)-c*A(k,j)
10      B(i)=B(i)-c*B(k)
100 continue


    if(abs((A(n,n)))-tol) 8,8,101
101 B(n)=B(n)/A(n,n)


    do 200 L=1,n1
        k=n-L
        k1=k+1
        do 200 j=k1,n
200 B(k)=B(k)-A(k,j)*B(j)


    d=1.
    do 250 i=1,n
250 d=d*A(i,i)
    goto 300
8   write(*,2) k
2   format(' **** singularity in row',i5)
    d=0.
    do i=1,n
        X(i)=B(i)
    enddo
300 return
    end subroutine BEM_Solver
    
        
        
    
    end Module Solver