subroutine getKernelParameters(elcoords,collocCoords,N,dN,&
        jacob, normals, r, dr, drdn)
    implicit none
    real(8)::elCoords(3,2),collocCoords(2),N(3),dN(3)
    real(8)::jacob,normals(2),r,dr(2),drdn
    !***************************************
    real(8)::dxydxi(2),fieldPt(2),relDist(2)
    integer(8)::i,j
    
    dxydxi(1)=0.d0;dxydxi(2)=0.d0
    do i=1,3
        dxydxi(1)=dxydxi(1)+dN(i)*elCoords(i,1)
        dxydxi(2)=dxydxi(2)+dN(i)*elCoords(i,2)
    enddo
    jacob=sqrt(dxydxi(1)**2+dxydxi(2)**2)
    normals(1)=dxydxi(2)/jacob
    normals(2)=-dxydxi(1)/jacob
    
    fieldPt(1)=0.d0;fieldPt(2)=0.d0
    do i=1,3
        fieldPt(1)=fieldPt(1)+N(i)*elCoords(i,1)
        fieldPt(2)=fieldPt(2)+N(i)*elCoords(i,2)
    enddo
    relDist(1)=fieldPt(1)-collocCoords(1)
    relDist(2)=fieldPt(2)-collocCoords(2)
    r=sqrt(relDist(1)**2+relDist(2)**2)
    dr(1)=relDist(1)/r
    dr(2)=relDist(2)/r
    drdn=dr(1)*normals(1)+dr(2)*normals(2)
    
    end
    