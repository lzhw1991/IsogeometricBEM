subroutine CalculateJumpTerm(n1,n2,Cij)
    use MaterialPara
    implicit none
    real(8)::n1(2),n2(2),Cij(2,2)
    !**************************
    real(8)::t(2,2),xVector(2),elmTheta(2),dotProd
    real(8)::theta,pi,thetaBar
    real(8)::term1,term2,term3,term4
    integer(8)::i,j
    
    do i=1,2
        do j=1,2
            cij(i,j)=0.d0
            t(i,j)=0.d0
        enddo
    enddo
    
    t(1,1)=-n1(2);t(1,2)=n1(1)
    t(2,1)=-n2(2);t(2,2)=n2(1)
    
    dotProd=n1(1)*n2(1)+n1(2)*n2(2)
    if(dotProd-1.d0>1.e-12) then
        theta=0.d0
    else
        theta=dacos(dotProd)
    endif
    
    pi=4.d0*datan(1.d0)
    if(n1(1)*n2(2)>n1(2)*n2(1)) then
        thetaBar=pi-theta
    else
        thetaBar=pi+theta
    endif
    
    xVector(1)=1.d0;xVector(2)=0.d0
    elmTheta(1)=0.d0;elmTheta(2)=0.d0
    
    do i=1,2
        theta=dacos(xVector(1)*t(i,1)+xVector(2)*t(i,2))
        if(xVector(1)*t(i,2)>xVector(2)*t(i,1)) then
            elmTheta(i)=theta
        else
            elmTheta(i)=2.d0*pi-theta
        endif
    enddo
    
    term1=1.d0/(8.d0*pi*(1.d0-mu))
    term2=4.d0*(1.d0-mu)*thetaBar
    term3=dsin(2.d0*elmTheta(1))-dsin(2.d0*elmTheta(2))
    term4=dcos(2.d0*elmTheta(2))-dcos(2.d0*elmTheta(1))
    
    Cij(1,1)=(term2+term3)*term1
    Cij(1,2)=term1*term4
    Cij(2,1)=Cij(1,2)
    Cij(2,2)=(term2-term3)*term1
    
    end 