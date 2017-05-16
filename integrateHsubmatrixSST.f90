subroutine integrateHsubmatrixSST(ngp,elCoords,glbBsFnConn,&
        glbCollocCoords,srcXi_param,srange,jumpTerm,Hsub)
    use GlobalPara
    use NURBS
    implicit none
    integer(8)::ngp,glbBsFnConn(3)
    real(8)::elCoords(3,2),glbCollocCoords(2)
    real(8)::srcXi_param,srange(2),jumpTerm(2,2)
    real(8)::Hsub(2,6)
    !***************************************
    integer(8)::i,j,k,pt
    real(8)::gpt(100),gwt(100)
    integer(8)::numBasisFns
    real(8),allocatable::srcN(:),srcdN(:),Rip(:),dRip(:)
    real(8)::nudgedXi,eps,srcXi
    integer(8)::c
    real(8)::hterm(2,2),htermMatrix(2,6),jacob_param
    real(8)::xi,xi_param,jacob_xi,normals(2),r,dr(2),drdn
    real(8)::jacob,Ttemp(2,2),tempMatrix(2,6),jacob_s,beta_m,sign
    real(8)::Utemp(2,2),jumpMatrix(2,6)
    
    eps=1.e-15
    call lgwt(ngp,-1.d0,1.d0,gpt,gwt)
    
    do i=1,2
        do j=1,6
            Hsub(i,j)=0.d0
        enddo
    enddo
    
    numBasisFns=p+1
    allocate(srcN(p+1),srcdN(p+1),Rip(p+1),dRip(p+1))
    do i=1,numBasisFns
        srcN(i)=0.d0;srcdN(i)=0.d0
        Rip(i)=0.d0;dRip(i)=0.d0
    enddo
    
    if(srcXi_param==srange(1)) then
        nudgedXi=srcXi_param+eps
    elseif(srcXi_param==srange(2)) then
        nudgedXi=srcXi_param-eps
    else
        nudgedXi=srcXi_param
    endif
    
    call convertToParentCoordSpace(srcXi_param,srange,srcXi)
    
    do c=1,numBasisFns
        i=glbBsFnConn(c)
        call NURBSbasis(i,p,nudgedXi,length_knotVec,knotVec,weights,srcN(c),srcdN(c))
    enddo
    
    hterm(1,1)=0.d0;
    hterm(1,2)=-const4*const3
    hterm(2,1)=const4*const3
    hterm(2,2)=0.d0
    do k=1,3
        do i=1,2
            do j=1,2
                htermMatrix(i,2*(k-1)+j)=hterm(i,j)*srcN(k)
            enddo
        enddo
    enddo
    
    jacob_param=(srange(2)-srange(1))/2.d0
    
    do pt=1,ngp
        
        xi=gpt(pt)
        call convertToParamSpace(xi,srange,xi_param)
        
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
                    tempMatrix(i,j+2*(k-1))=Rip(k)*Ttemp(i,j)*(xi-srcXi)*jacob*gwt(pt)
                enddo
            enddo
        enddo
        
       
        do i=1,2
            do j=1,6
                tempMatrix(i,j)=tempMatrix(i,j)-htermMatrix(i,j)
            enddo
        enddo
        
        do i=1,2
            do j=1,6
                tempMatrix(i,j)=tempMatrix(i,j)/(xi-srcXi)
            enddo
        enddo
        
        do i=1,2
            do j=1,6
                Hsub(i,j)=Hsub(i,j)+tempMatrix(i,j)
            enddo
        enddo
    enddo
    
    if(abs(srcXi)<1.d0-100.d0*eps) then
        do i=1,2
            do j=1,6
                Hsub(i,j)=Hsub(i,j)+htermMatrix(i,j)*log(abs((1.d0-srcXi)/(1.d0+srcXi)))
            enddo
        enddo
    else
        call getKernelParameters(elcoords,glbCollocCoords,srcN,srcdN,&
            jacob_xi,normals,r,dr,drdn)
        jacob_s=jacob_xi*jacob_param
        beta_m=1.d0/jacob_s
        if(xi>srcXi) then
            sign=1.d0
        elseif(xi<srcXi) then
            sign=-1.d0
        else
            sign=0.d0
        endif
        do i=1,2
            do j=1,6
                Hsub(i,j)=Hsub(i,j)+htermMatrix(i,j)*log(abs(2.d0/beta_m))*sign
            enddo
        enddo
    endif
    
    if(abs(srcXi_param-srange(2))>100.d0*eps) then
        do k=1,3
            do i=1,2
                do j=1,2
                    jumpMatrix(i,j+2*(k-1))=jumpTerm(i,j)*srcN(k)
                enddo
            enddo
        enddo
        do i=1,2
            do j=1,6
                Hsub(i,j)=Hsub(i,j)+jumpMatrix(i,j)
            enddo
        enddo
    endif
    
    end
        
        
    
    
    