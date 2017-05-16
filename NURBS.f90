Module NURBS
    contains
    
    subroutine FindSpan(n,p,u,knot,span)
    implicit none
    integer(8),intent(in)::n,p
    real(8),intent(in)::u,knot(0:)
    integer(8),intent(out)::span
    !**********************************
    integer(8)::low,high,mid
    !**********************************
    if(u>=knot(n+1)) then
        span=n
        return
    endif
    
    if(u<=knot(p)) then
        span=p
        return
    endif
    
    low=p;high=n+1;mid=(low+high)/2
    do while((u<knot(mid)).or.(u>=knot(mid+1)))
        if(u<knot(mid)) then
            high=mid
        else
            low=mid
        endif
        mid=(low+high)/2
    enddo
    span=mid
    return
    end subroutine FindSpan
    
    subroutine BasisFuns(i,u,p,knot,N)
    implicit none
    integer(8),intent(in)::i,p
    real(8),intent(in)::u,knot(0:)
    real(8),intent(in out)::N(0:)
    !*********************************
    integer(8)::j,r
    real(8),allocatable::left(:),right(:)
    real(8)::saved,temp
    !*********************************
    allocate(left(0:p),right(0:p))
    
    N(0)=1.d0
    
    do j=1,p
        left(j)=u-knot(i+1-j)
        right(j)=knot(i+j)-u
        saved=0.d0
        do r=0,j-1
            temp=N(r)/(right(r+1)+left(j-r))
            N(r)=saved+right(r+1)*temp
            saved=left(j-r)*temp
        enddo
        N(j)=saved
    enddo
    
    deallocate(left,right)
    
    end subroutine BasisFuns
    
    subroutine DersBasisFuns(i,u,p,order,knot,ders)
    implicit none
    integer(8),intent(in)::i,p,order
    real(8),intent(in)::u,knot(0:)
    real(8),intent(in out)::ders(0:,0:)
    !**********************************************
    real(8)::saved,temp,d
    integer(8)::j,k,j1,j2,r
    integer(8)::s1,s2,rk,pk
    real(8),allocatable::left(:),right(:),ndu(:,:),a(:,:)
    !**********************************************
    allocate(left(0:p),right(0:p),ndu(0:p,0:p),a(0:p,0:p))
    
    ndu(0,0)=1.d0
    do j=1,p
        left(j)=u-knot(i+1-j)
        right(j)=knot(i+j)-u
        saved=0.d0
        do r=0,j-1
            ndu(j,r)=right(r+1)+left(j-r)
            temp=ndu(r,j-1)/ndu(j,r)
            
            ndu(r,j)=saved+right(r+1)*temp
            saved=left(j-r)*temp
        enddo
        ndu(j,j)=saved
    enddo
    
    do j=0,p
        ders(0,j)=ndu(j,p)
    enddo
    
    if(order==0) then
        return
    endif
    
    do r=0,p
        s1=0;s2=1
        a(0,0)=1.d0
        do k=1,order
            d=0.d0
            rk=r-k;pk=p-k
            if(r>=k) then
                a(s2,0)=a(s1,0)/ndu(pk+1,rk)
                d=a(s2,0)*ndu(rk,pk)
            endif
            
            if(rk>=-1) then
                j1=1
            else
                j1=-rk
            endif
            
            if(r-1<=pk) then
                j2=k-1
            else
                j2=p-r
            endif
            
            do j=j1,j2
                a(s2,j)=(a(s1,j)-a(s1,j-1))/ndu(pk+1,rk+j)
                d=d+a(s2,j)*ndu(rk+j,pk)
            enddo
            
            if(r<=pk) then
                a(s2,k)=-a(s1,k-1)/ndu(pk+1,r)
                d=d+a(s2,k)*ndu(r,pk)
            endif
            ders(k,r)=d
            j=s1;s1=s2;s2=j
        enddo
    enddo
    r=p
    do k=1,order
        do j=0,p
            ders(k,j)=ders(k,j)*r
        enddo
        r=r*(p-k)
    enddo
    
    deallocate(left,right,ndu,a)
    end subroutine DersBasisFuns
    
    subroutine OneBasisFun(p,m,knot,i,u,Nip)
    implicit none
    integer(8),intent(in)::p,m,i
    real(8),intent(in)::knot(0:),u
    real(8),intent(out)::Nip
    !***************************************
    integer(8)::j,k
    real(8)::saved,Uleft,Uright,temp
    real(8),allocatable::N(:)
    !***************************************
    allocate(N(0:p))
    
    if((i==0.and.u==knot(0)).or.(i==m-p-1.and.u==knot(m))) then
        Nip=1.d0
        return
    endif
    
    if((u<knot(i)).or.(u>=knot(i+p+1))) then
        Nip=0.d0
        return
    endif
    
    do j=0,p
        if((u>=knot(i+j)).and.(u<knot(i+j+1))) then
            N(j)=1.d0
        else
            N(j)=0.d0
        endif
    enddo
    
    do k=1,p
        if(N(0)==0.d0) then
            saved=0.d0
        else
            saved=((u-knot(i))*N(0))/(knot(i+k)-knot(i))
        endif
        do j=0,p-k
            Uleft=knot(i+j+1)
            Uright=knot(i+j+k+1)
            if(N(j+1)==0.d0) then
                N(j)=saved
                saved=0.d0
            else
                temp=N(j+1)/(Uright-Uleft)
                N(j)=saved+(Uright-u)*temp
                saved=(u-Uleft)*temp
            endif
        enddo
    enddo
    
    Nip=N(0)
    
    deallocate(N)
    
    end subroutine OneBasisFun
    
    subroutine DersOneBasisFuns(p,m,knot,i,u,order,ders)
    implicit none
    integer(8),intent(in)::p,m,i,order
    real(8),intent(in)::knot(0:),u
    real(8),intent(in out)::ders(0:)
    !***************************************************
    integer(8)::k,j,jj
    real(8)::Uleft,Uright,saved,temp
    real(8),allocatable::N(:,:),ND(:)
    !***************************************************
    allocate(N(0:order,0:order),ND(0:order))
    
    if((u<knot(i)).or.(u>=knot(i+p+1))) then
        do k=0,order
            ders(k)=0.d0
        enddo
        return
    endif
    
    do j=0,p
        if((u>=knot(i+j)).and.(u<knot(i+j+1))) then
            N(j,0)=1.d0
        else
            N(j,0)=0.d0
        endif
    enddo
    
    do k=1,p
        if(N(0,k-1)==0.d0) then
            saved=0.d0
        else
            saved=((u-knot(i))*N(0,k-1))/(knot(i+k)-knot(i))
        endif
        do j=0,p-k
            Uleft=knot(i+j+1)
            Uright=knot(i+j+k+1)
            if(N(j+1,k-1)==0.d0) then
                N(j,k)=saved
                saved=0.d0
            else
                temp=N(j+1,k-1)/(Uright-Uleft)
                N(j,k)=saved+(Uright-u)*temp
                saved=(u-Uleft)*temp
            endif
        enddo
    enddo
    
    ders(0)=N(0,p)
    
    do k=1,order
        do j=0,k
            ND(j)=N(j,p-k)
        enddo
        do jj=1,k
            if(ND(0)==0.d0) then
                saved=0.d0
            else
                saved=ND(0)/(knot(i+p-k+jj)-knot(i))
            endif

            do j=0,k-jj
                Uleft=knot(i+j+1)
                Uright=knot(i+j+p+jj) !***
                if(ND(j+1)==0.d0) then
                    ND(j)=(p-k+jj)*saved
                    saved=0.d0
                else
                    temp=ND(j+1)/(Uright-Uleft)
                    ND(j)=(p-k+jj)*(saved-temp)
                    saved=temp
                endif
            enddo
        enddo
        ders(k)=ND(0)
    enddo
    
    deallocate(N,ND)
    
    end subroutine DersOneBasisFuns
    
    subroutine NURBSbasis(ii,p,xi,mm,knot,weight,Rip,dRip)
    implicit none
    integer(8),intent(in)::ii,p,mm
    real(8),intent(in)::xi,knot(0:),weight(0:)
    real(8),intent(out)::Rip,dRip
    !*************************************************
    integer(8)::i,m,c,k,n,span,numKnot,numWeights
    real(8)::Nip,w_interp,dw_interp_dxi,eps,srcXi
    real(8),allocatable::NN(:),dN(:),ders(:,:)
    !*************************************************
    i=ii-1
    m=mm-1
    n=m-p-1
    eps=1.e-12
    
    allocate(NN(0:p),dN(0:p),ders(0:p,0:p))
    if(dabs(xi-knot(m))<eps) then
        srcXi=knot(m)-eps
    else
        srcXi=xi
    endif
        
    call FindSpan(n,p,srcXi,knot,span)
    call BasisFuns(span,srcXi,p,knot,NN)
    call DersBasisFuns(span,srcXi,p,p,knot,ders)
    
    w_interp=0.d0;dw_interp_dxi=0.d0
    
    do c=0,p
        w_interp=w_interp+NN(c)*weight(span-p+c)
        dw_interp_dxi=dw_interp_dxi+ders(1,c)*weight(span-p+c)
    enddo
    
    call DersOneBasisFuns(p,m,knot,i,srcXi,p,dN)
    call OneBasisFun(p,m,knot,i,srcXi,Nip)
    
    Rip=Nip*weight(i)/w_interp
    dRip=weight(i)*(w_interp*dN(1)-dw_interp_dxi*Nip)/(w_interp**2)
    
    deallocate(ders,NN,dN)
    
    end subroutine NURBSbasis
    !************************************************************
    !*******                                           **********
    !************************************************************
    subroutine NURBSinterpolation(src,p,length_knotVec,knot,points,weight,interp,interp_deriv)
    implicit none
    real(8),intent(in)::src,knot(0:),points(:),weight(0:)
    integer(8),intent(in)::p,length_knotVec
    real(8),intent(in out)::interp,interp_deriv
    !******************************************************
    integer(8)::m,n,span,k,c
    real(8)::xi,w_interp,dw_interp_dxi,eps
    real(8),allocatable::NN(:),NURBS(:),NURBS_deriv(:),ders(:,:),tpoints(:)
    !*******************************************************
    m=length_knotVec-1
    n=m-p-1
    eps=1.e-16
    allocate(NN(0:p),NURBS(0:p),NURBS_deriv(0:p),ders(0:n,0:p),tpoints(0:n+2))
    do k=1,n
        tpoints(k-1)=points(k)
    enddo
    tpoints(n)=points(1)
    tpoints(n+1)=points(n)
    if(abs(knot(m)-src)<eps) then
        xi=knot(m)-eps
    else
        xi=src
    endif
    call FindSpan(n,p,xi,knot,span)
    call BasisFuns(span,xi,p,knot,NN)
    call DersBasisFuns(span,xi,p,n,knot,ders)
    do k=0,p
        w_interp=0.d0;dw_interp_dxi=0.d0
        do c=0,p
            w_interp=w_interp+NN(c)*weight(span-p+c)
            dw_interp_dxi=dw_interp_dxi+ders(1,c)*weight(span-p+c)
        enddo
        NURBS(k)=NN(k)*weight(span-p+k)/w_interp
        NURBS_deriv(k)=weight(span-p+k)*(w_interp*ders(1,k)-dw_interp_dxi*NN(k))/(w_interp**2)
    enddo
    interp=0.d0;interp_deriv=0.d0
    do c=0,p
        interp=interp+NURBS(c)*tpoints(span-p+c)
        interp_deriv=interp_deriv+NURBS_deriv(c)*tpoints(span-p+c)
    enddo
    
    deallocate(NN,NURBS,NURBS_deriv,ders,tpoints)
    
    end subroutine NURBSinterpolation
    
    end Module NURBS