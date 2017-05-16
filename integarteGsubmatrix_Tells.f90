subroutine integrateGsubmatrix_Tells(ngp,elcoords,&
        glbBsFnConn,glbCollocCoords,srcXi_param,srange,Gsub)
    use GlobalPara
    use NURBS
    implicit none
    integer(8)::ngp,glbBsFnConn(p+1)
    real(8)::elCoords(3,2),glbCollocCoords(2),srcXi_param,srange(2)
    real(8)::Gsub(2,6)
    !***************************************
    integer(8)::numBasisFns,pt,i,j,k
    real(8)::gpt(100),gwt(100)
    real(8)::Rip(p+1),dRip(p+1),srcXi,xiStar
    real(8)::gamBar,jacob_param,gXi,xi,xi_param
    real(8)::jacob_xi,normals(2),r,dr(2),drdn,jacob
    real(8)::Ttemp(2,2),Utemp(2,2),jacobTelles
    real(8)::nthroot
    
    call lgwt(ngp,-1.d0,1.d0,gpt,gwt)
    
    
    numBasisFns=p+1
    do i=1,numBasisFns
        Rip(i)=0.d0
        dRip(i)=0.d0
    enddo
    
    do i=1,2
        do j=1,6
            Gsub(i,j)=0.d0
        enddo
    enddo
    
    call convertToParentCoordSpace(srcXi_param,srange,srcXi)
    xiStar=srcXi**2-1
    gamBar=nthroot(srcXi*xiStar+abs(xiStar),3)+nthroot(srcXi*xiStar-abs(xiStar),3)+srcXi
    jacob_param=(srange(2)-srange(1))/2.d0
    
    do pt=1,ngp
        gXi=gpt(pt)
        xi=((gXi-gamBar)**3+gamBar*(gamBar**2+3.d0))/(1.d0+3.d0*gamBar**2)
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
        jacobTelles=3.d0*((gXi-gamBar)**2)/(1.d0+3.d0*gamBar**2)
        
        do k=1,3
            do i=1,2
                do j=1,2
                    Gsub(i,j+2*(k-1))=Gsub(i,j+2*(k-1))+&
                        Rip(k)*Utemp(i,j)*jacob*gwt(pt)*jacobTelles
                enddo
            enddo
        enddo
    enddo
        
    end
    
    real(8) function nthroot(x,n)
    implicit none
    real(8)::x
    integer(8)::n
    real(8)::root
    !*******************
    real(8)::temp
    temp=dabs(x)
    root=temp**(1.d0/n)
    if(x<0.d0) then
        root=-1.d0*root
    endif
    nthroot=root
    end function nthroot