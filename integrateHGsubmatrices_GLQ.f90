subroutine integrateHGsubmatrices_GLQ(ngp,elcoords,&
        glbBsFnConn,glbCollocCoords,srange,Hsub,Gsub)
    use GlobalPara
    use NURBS
    implicit none
    integer(8)::ngp,glbBsFnConn(p+1)
    real(8)::elCoords(p+1,2),glbCollocCoords(2),srange(2)
    real(8)::Hsub(2,6),Gsub(2,6)
    !*********************************************
    real(8)::gpt(100),gwt(100)
    integer(8)::numBasisFns,i,j,k,pt
    real(8)::Rip(p+1),dRip(p+1),jacob_param,xi_param
    real(8)::jacob_xi,normals(2),r,dr(2),drdn,jacob,Ttemp(2,2),Utemp(2,2)
    
    call lgwt(ngp,-1.d0,1.d0,gpt,gwt)
    
    numBasisFns=p+1
    
    do i=1,2
        do j=1,6
            Hsub(i,j)=0.d0
            Gsub(i,j)=0.d0
        enddo
    enddo
    
    jacob_param=(srange(2)-srange(1))/2.d0
    
    do pt=1,ngp
        call convertToParamSpace(gpt(pt),srange,xi_param)
        
        do j=1,numBasisFns
            i=glbBsFnConn(j)
            call NURBSbasis(i,p,xi_param,length_knotVec,knotVec,weights,Rip(j),dRip(j))
        enddo
        call getKernelParameters(elcoords,glbCollocCoords,Rip,dRip,&
            jacob_xi,normals,r,dr,drdn)
        
        jacob=jacob_xi*jacob_param
        
        do i=1,2
            do j=1,2
                call DBIEkernels(i,j,r,dr,drdn,normals,Ttemp(i,j),Utemp(i,j))
            enddo
        enddo
        
        do k=1,3
            do i=1,2
                do j=1,2
                    Hsub(i,j+2*(k-1))=Hsub(i,j+2*(k-1))+Rip(k)*Ttemp(i,j)*jacob*gwt(pt)
                    Gsub(i,j+2*(k-1))=Gsub(i,j+2*(k-1))+Rip(k)*Utemp(i,j)*jacob*gwt(pt)
                enddo
            enddo
        enddo
    enddo
    
    end