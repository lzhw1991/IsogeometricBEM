subroutine findNormals(collocNormals)
    use GlobalPara
    implicit none
    real(8)::collocNormals(10000,2,2)
    !*******************************
    integer(8)::c,element
    real(8)::eps,srcXi_param
    integer(8)::i,j,k
    real(8)::range(2),elCoords(3,2),n(2)
    integer(8)::glbBsFnConn(3),dElConn(3)
    
    do k=1,nPts
        do i=1,2
            do j=1,2
                collocNormals(k,i,j)=0.d0
            enddo
        enddo
    enddo
    
                
    eps=1.e-12
    
    do c=1,nPts
        srcXi_param=collocPts(c)
        do element=1,ne
            do j=1,2
                range(j)=elRange(element,j)
            enddo
            do j=1,p+1
                glbBsFnConn(j)=bsFnConn(element,j)
                dElConn(j)=dispConn(element,j)
            enddo
            do j=1,p+1
                elCoords(j,1)=controlPts(dElConn(j),1)
                elCoords(j,2)=controlPts(dElConn(j),2)
            enddo
            
            if((srcXi_param<=range(2)).and.((srcXi_param>=range(1)).or.((element==ne).and.(srcXi_param==0.d0)))) then
                if(abs(srcXi_param-range(2))<eps) then
                    call getNormal(elCoords,srcXi_param-eps,glbBsFnConn,n)
                    collocNormals(c,1,1)=n(1)
                    collocNormals(c,2,1)=n(2)
                elseif(abs(srcXi_param-range(1))<eps) then
                    call getNormal(elCoords,srcXi_param,glbBsFnConn,n)
                    collocNormals(c,1,2)=n(1)
                    collocNormals(c,2,2)=n(2)
                elseif((element==ne).and.(srcXi_param==0.d0)) then
                    srcXi_param=1.d0-eps
                    call getNormal(elCoords,1.d0-eps,glbBsFnConn,n)
                    collocNormals(c,1,1)=n(1)
                    collocNormals(c,2,1)=n(2)
                else
                    call getNormal(elCoords,srcXi_param,glbBsFnConn,n)
                    collocNormals(c,1,1)=n(1)
                    collocNormals(c,2,1)=n(2)
                    
                    collocNormals(c,1,2)=n(1)
                    collocNormals(c,2,2)=n(2)
                endif
            endif
        enddo 
    enddo
    
    end
    
    subroutine getNormal(elCoords,xi_param,glbBsFnConn,collocNormal)
    use GlobalPara
    use NURBS
    implicit none
    real(8)::elCoords(3,2),xi_param
    integer(8)::glbBsFnConn(3)
    real(8)::collocNormal(2)
    !************************
    real(8)::Rip(3),dRip(3)
    integer(8)::i,j
    real(8)::dxydxi(2),jacob
    
    collocNormal(1)=0.d0;collocNormal(2)=0.d0
    do i=1,p+1
        j=glbBsFnConn(i)
        call NURBSbasis(j,p,xi_param,length_knotVec,knotVec,weights,Rip(i),dRip(i))
    enddo
    
    dxydxi(1)=0.d0;dxydxi(2)=0.d0
    do i=1,p+1
        dxydxi(1)=dxydxi(1)+dRip(i)*elCoords(i,1)
        dxydxi(2)=dxydxi(2)+dRip(i)*elCoords(i,2)
    enddo
    jacob=sqrt(dxydxi(1)**2+dxydxi(2)**2)
    collocNormal(1)=dxydxi(2)/jacob
    collocNormal(2)=-dxydxi(1)/jacob
    
    end