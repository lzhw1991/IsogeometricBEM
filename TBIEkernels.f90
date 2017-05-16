subroutine TBIEkernels(i,j,k,r,dr,drdn,normal,Skij,Dkij)
    use GlobalPara
    use MaterialPara
    implicit none
    integer(8)::i,j,k
    real(8)::r,dr(2),drdn,normal(2),Skij,Dkij
    real(8)::ij,ki,jk,pi
    
    pi=4.d0*datan(1.d0)
    if(i==j) then
        ij=1.d0
    else
        ij=0.d0
    endif
    
    if(k==i) then
        ki=1.d0
    else
        ki=0.d0
    endif
    
    if(j==k) then
        jk=1.d0
    else
        jk=0.d0
    endif
    
    Skij=shearMod/(2.d0*pi*(1.d0-mu))*1.d0/(r**2)
    Skij=Skij*(2.d0*drdn*(const4*ij*dr(k)+mu*(dr(j)*ki+dr(i)*jk)-&
        4.d0*dr(i)*dr(j)*dr(k))+2.d0*mu*(normal(i)*dr(j)*dr(k)+&
        normal(j)*dr(i)*dr(k))+const4*(2.d0*normal(k)*dr(i)*dr(j)+&
        normal(j)*ki+normal(i)*jk)-(1.d0-4.d0*mu)*normal(k)*ij)

    Dkij=const3*(1.d0/r)*(const4*(-dr(k)*ij+dr(j)*ki+dr(i)*jk)+&
        2.d0*dr(i)*dr(j)*dr(k))
    
end