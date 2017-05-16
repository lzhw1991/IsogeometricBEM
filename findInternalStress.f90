subroutine findInternalStress(displacements,tractions,XY,stress)
    use GlobalPara
    use NURBS
    implicit none
    real(8)::displacements(nDof),tractions(tracNdof)
    real(8)::XY(2),stress(2,2)
    !**************************************
    integer(8)::ngp,element,pt,i,j,k
    real(8)::gpt(100),gwt(100),range(2),elCoords(3,2)
    integer(8)::glbBsFnConn(p+1),dElConn(p+1),tElConn(p+1)
    integer(8)::dispSctrX(p+1),dispSctrY(p+1),tracSctrX(p+1),tracSctrY(p+1)
    real(8)::jacob_param,xi_param,Rip(p+1),dRip(p+1)
    real(8)::trac(2),disp(2)
    real(8)::jacob_xi,normal(2),r,dr(2),drdn,jacob
    real(8)::tempS(2,2,2),tempD(2,2,2)
    
    ngp=40
    do i=1,2
        do j=1,2
            stress(i,j)=0.d0
        enddo
    enddo
    call lgwt(ngp,-1.d0,1.d0,gpt,gwt)
    
    do element=1,ne
        do i=1,2
            range(i)=elRange(element,i)
        enddo
        do i=1,p+1
            glbBsFnConn(i)=bsFnConn(element,i)
            dElConn(i)=dispConn(element,i)
            tElConn(i)=tracConn(element,i)
        enddo
        do i=1,p+1
            elCoords(i,1)=controlPts(dElConn(i),1)
            elCoords(i,2)=controlPts(dElConn(i),2)
        enddo
        do i=1,p+1
            dispSctrX(i)=2*dElConn(i)-1
            dispSctrY(i)=2*dElConn(i)
            tracSctrX(i)=2*tElConn(i)-1
            tracSctrY(i)=2*tElConn(i)
        enddo
        
        jacob_param=(range(2)-range(1))/2.d0
        
        do pt=1,ngp
            call convertToParamSpace(gpt(pt),range,xi_param)
            do j=1,p+1
                i=glbBsFnConn(j)
                call NURBSbasis(i,p,xi_param,length_knotVec,knotVec,weights,Rip(j),dRip(j))
            enddo
            trac(1)=0.d0;trac(2)=0.d0
            disp(1)=0.d0;disp(2)=0.d0
            do j=1,p+1
                trac(1)=trac(1)+Rip(j)*tractions(tracSctrX(j))
                trac(2)=trac(2)+Rip(j)*tractions(tracSctrY(j))
                disp(1)=disp(1)+Rip(j)*displacements(dispSctrX(j))
                disp(2)=disp(2)+Rip(j)*displacements(dispSctrY(j))
            enddo
            call getKernelParameters(elCoords,XY,Rip,dRip,&
            jacob_xi,normal,r,dr,drdn)
            
            jacob=jacob_xi*jacob_param
            
            do i=1,2
                do j=1,2
                    do k=1,2
                        tempS(i,j,k)=0.d0
                        tempD(i,j,k)=0.d0
                    enddo
                enddo
            enddo
            
            do i=1,2
                do j=1,2
                    do k=1,2
                        call TBIEkernels(i,j,k,r,dr,drdn,normal,&
                            tempS(i,j,k),tempD(i,j,k))
                    enddo
                enddo
            enddo
            
            do i=1,2
                do j=1,2
                    stress(i,j)=stress(i,j)+((tempD(i,j,1)*trac(1)+&
                        tempD(i,j,2)*trac(2))-(tempS(i,j,1)*disp(1)+&
                        tempS(i,j,2)*disp(2)))*gwt(pt)*jacob
                enddo
            enddo
        enddo
    enddo
    
    end
        
            
            
            