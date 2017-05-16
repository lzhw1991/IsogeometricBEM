subroutine DBIEkernels(i,j,r,dr,drdn,normal,Tij,Uij)
    use GlobalPara
    implicit none
    integer(8)::i,j
    real(8)::ij,r,dr(2),drdn,normal(2),Tij,Uij
    
    if(i==j) then
        ij=1.d0
    else
        ij=0.d0
    endif
    
    Tij=(1.d0/r)*const3*(-drdn*(const4*ij+2.d0*(dr(i)*dr(j)))+&
        const4*(dr(i)*normal(j)-dr(j)*normal(i)))
    
    Uij=const1*(dr(i)*dr(j)+ij*(const2*log(1.d0/r)))
    end