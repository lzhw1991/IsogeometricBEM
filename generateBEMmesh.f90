    subroutine generateBEMmesh(refinement)
    use GlobalPara
    use NURBS
    implicit none
    integer::refinement
    !*******************************
    integer(8)::m,n,i,j,e
    real(8),allocatable::uniqueKnots(:),weightedPts(:,:),oldKnotVec(:)
    integer(8)::length_uniqueKnots
    real(8),allocatable::newKnots(:,:),distributed(:)
    integer(8)::length_newKnots
    real(8),allocatable::newControlPts(:,:)
    integer(8)::length_newControlPts,knot,kval
    real(8)::knot_bar
    real(8),allocatable::temp(:,:)
    real(8)::alpha
    real(8),allocatable::newPoints(:,:)
    integer(8),allocatable::elConn(:,:),elKnotIndices(:,:)
    integer(8)::element
    real(8),allocatable::previousKnotVals(:),currentKnotVals(:)
    integer(8)::numRepeatedKnots
    integer(8),allocatable::indices(:)
    integer(8)::isequal,length_nonzeros_previousKnotVals
    integer(8)::totalTractionDOF
    real(8),allocatable::xi(:)
    integer(8)::length_xi,numBasisFns
    real(8),allocatable::NURBSvalues(:,:),NURBSderivVals(:,:),r(:,:)
    real(8)::Rip,dRip
    integer(8)::point
    real(8),allocatable::NURBSCoords(:,:)
    
    integer(8)::ipr  !文件输出
    
    do i=1,length_knotVec
        knotVec(i)=knotVec(i)/knotVec(length_knotVec)
    enddo
    
    do i=1,length_ControlPts
        controlPts(i,3)=weights(i)
    enddo
    allocate(uniqueKnots(length_knotVec))
    uniqueKnots(1)=knotVec(1)
    j=1
    do i=2,length_knotVec
        if(knotVec(i).ne.knotVec(i-1)) then
            j=j+1
            uniqueKnots(j)=knotVec(i)
        endif
    enddo
    length_uniqueKnots=j
    n=length_knotVec-1-p
    allocate(weightedPts(length_controlPts,3))
    do i=1,length_controlPts
        weightedPts(i,1)=controlPts(i,1)*controlPts(i,3)
        weightedPts(i,2)=controlPts(i,2)*controlPts(i,3)
        weightedPts(i,3)=controlPts(i,3)
    enddo
    
    if(refinement>1) then
        allocate(oldKnotVec(length_knotVec))
        do i=1,length_knotVec
            oldKnotVec(i)=knotVec(i)
        enddo
        m=length_uniqueKnots-1;n=refinement-1
        allocate(newKnots(m,n))
        allocate(distributed(refinement+1))
        do i=1,m
            do j=1,refinement+1
                distributed(j)=uniqueKnots(i)+&
                    (j-1)*(uniqueKnots(i+1)-uniqueKnots(i))/refinement
            enddo
            do j=2,refinement
                newKnots(i,j-1)=distributed(j)
            enddo
        enddo
        deallocate(distributed)
        deallocate(uniqueKnots)
        allocate(temp(m,n))
        do i=1,m
            do j=1,n
                temp(i,j)=newKnots(i,j)
            enddo
        enddo
        deallocate(newKnots)
        allocate(newKnots(m*n,1))
        do i=1,m
            do j=1,n
                newKnots((i-1)*n+j,1)=temp(i,j)
            enddo
        enddo
        length_newKnots=m*n
        deallocate(temp)
        allocate(temp(1,1))
        do i=1,length_newKnots-1
            do j=1,length_newKnots-i
                if(newKnots(j,1)>=newKnots(j+1,1)) then
                    temp(1,1)=newKnots(j,1)
                    newKnots(j,1)=newKnots(j+1,1)
                    newKnots(j+1,1)=temp(1,1)
                endif
            enddo
        enddo
        deallocate(temp)
        
        do knot=1,length_newKnots
            allocate(newControlPts(length_controlPts+1,3))
            length_newControlPts=length_ControlPts+1
            knot_bar=newKnots(knot,1)
            do i=1,length_knotVec
                if(knotVec(i)>knot_bar) then
                    kval=i-1  
                    exit
                endif
            enddo
            allocate(newPoints(1,3))
            do i=1,length_newControlPts
                if(i<=(kval-p)) then
                    alpha=1.d0
                elseif((i>=kval-p+1).and.(i<=kval)) then
                    alpha=(newKnots(knot,1)-oldKnotVec(i))/(oldKnotVec(i+p)-oldKnotVec(i))
                else
                    alpha=0.d0
                endif
                do j=1,3
                    newPoints(1,j)=0.d0
                enddo
                if(i.ne.1) then
                    do j=1,3
                        newPoints(1,j)=(1-alpha)*weightedPts(i-1,j)
                    enddo
                endif
                if(i.ne.length_newControlPts) then
                    do j=1,3
                        newPoints(1,j)=newPoints(1,j)+alpha*weightedPts(i,j)
                    enddo
                endif
                do j=1,3
                    newControlPts(i,j)=newPoints(1,j)
                enddo
            enddo
            deallocate(newPoints)
            deallocate(weightedPts)
            length_ControlPts=length_newControlPts
            allocate(weightedPts(length_ControlPts,3))
            do i=1,length_newControlPts
                do j=1,3
                    weightedPts(i,j)=newControlPts(i,j)
                enddo
            enddo
            length_knotVec=length_knotVec+1
            knotVec(length_knotVec)=knot_bar
            allocate(temp(1,1))
            do i=1,length_knotVec-1
                do j=1,length_knotVec-i
                    if(knotVec(j)>=knotVec(j+1)) then
                        temp(1,1)=knotVec(j)
                        knotVec(j)=knotVec(j+1)
                        knotVec(j+1)=temp(1,1)
                    endif
                enddo
            enddo
            deallocate(oldKnotVec)
            allocate(oldKnotVec(length_knotVec))
            do i=1,length_knotVec
                oldKnotVec(i)=knotVec(i)
            enddo
            deallocate(temp)
            deallocate(newControlPts)
        enddo
        
        do i=1,length_controlPts
            controlPts(i,1)=weightedPts(i,1)/weightedPts(i,3)
            controlPts(i,2)=weightedPts(i,2)/weightedPts(i,3)
            controlPts(i,3)=weightedPts(i,3)
        enddo
        deallocate(weightedPts)
        
        allocate(uniqueKnots(length_KnotVec))
        uniqueKnots(1)=knotVec(1)
        j=1
        do i=2,length_KnotVec
            if(knotVec(i).ne.knotVec(i-1)) then
                j=j+1
                uniqueKnots(j)=knotVec(i)
            endif
        enddo
        length_uniqueKnots=j
    endif
    
    ! --------------------------------------------------------
    ! ------------- Define element connectivities ------------
    ! --------------------------------------------------------
    ne=length_uniqueKnots-1
    allocate(elConn(ne,p+1))
    allocate(elKnotIndices(length_KnotVec,2))
    do i=1,ne
        elRange(i,:)=0.d0
        elConn(i,:)=0
        elKnotIndices(i,:)=0
        tracConn(i,:)=0
    enddo
    allocate(previousKnotVals(1),currentKnotVals(1))
    element=1
    previousKnotVals(1)=0.d0
    do i=1,length_knotVec
        currentKnotVals(1)=knotVec(i)
        if(knotVec(i).ne.previousKnotVals(1)) then
            elRange(element,1)=previousKnotVals(1)
            elRange(element,2)=currentKnotVals(1)
            elKnotIndices(element,1)=i-1
            elKnotIndices(element,2)=i
            element=element+1
        endif
        previousKnotVals(1)=currentKnotVals(1)
    enddo
    numRepeatedKnots=0
    allocate(indices(p))
    deallocate(previousKnotVals,currentKnotVals)
    allocate(previousKnotVals(p),currentKnotVals(p))
    do e=1,ne
        do i=elKnotIndices(e,1)-p+1,elKnotIndices(e,1)
            indices(i-elKnotIndices(e,1)+p)=i
        enddo
        do i=1,p
            previousKnotVals(i)=knotVec(indices(i))
            currentKnotVals(i)=knotVec(elKnotIndices(e,1))
        enddo
        isequal=1
        do i=1,p
            if(previousKnotVals(i).ne.currentKnotVals(i)) then
                isequal=0
                exit
            endif
        enddo
        length_nonzeros_previousKnotVals=0
        do i=1,p
            if(previousKnotVals(i).ne.0.d0) then
                length_nonzeros_previousKnotVals=length_nonzeros_previousKnotVals+1
            endif
        enddo
        
        if((isequal.eq.1).and.(length_nonzeros_previousKnotVals>1)) then
            numRepeatedKnots=numRepeatedKnots+1
        endif
        do i=elKnotIndices(e,1)-p,elKnotIndices(e,1)
            elConn(e,i+1-elKnotIndices(e,1)+p)=i
        enddo
        do i=1,p+1
            tracConn(e,i)=elConn(e,i)+numRepeatedKnots
        enddo
    enddo
    do e=1,ne
        do j=1,p+1
            bsFnConn(e,j)=elConn(e,j)
            dispConn(e,j)=elConn(e,j)
        enddo
    enddo
    dispConn(ne,p+1)=1
    elConn(ne,p+1)=1
    totalTractionDOF=length_knotVec-p-1+numRepeatedKnots
    do e=1,ne
        do j=1,p+1
            tracDispConn(tracConn(e,j))=elConn(e,j)
        enddo
    enddo
    length_tracDispConn=totalTractionDOF
    
    length_xi=1000
    allocate(xi(length_xi))
    do i=1,length_xi
        xi(i)=(i-1)*knotVec(length_knotVec)/(length_xi-1)
    enddo
    n=length_xi
    numBasisFns=length_knotVec-1-p
    allocate(NURBSvalues(n,numBasisFns),NURBSderivVals(n,numBasisFns))
    do i=1,numBasisFns
        CollocPts(i)=0.d0
        do j=i+1,i+p
            CollocPts(i)=CollocPts(i)+knotVec(j)
        enddo
        CollocPts(i)=CollocPts(i)/p
    enddo
    allocate(r(n,numBasisFns))
    !*************************
    !***注意！在NURBSbasis和NURBSinterpolation函数中
    !***KnotVec和Weights数据下边均是从0开始，所以需处理一下
    do i=1,length_knotVec
        knotVec(i-1)=knotVec(i)
    enddo
    knotVec(length_knotVec+1)=knotVec(length_knotVec)
    do i=1,length_ControlPts
        weights(i)=controlPts(i,3)
    enddo
    do i=1,length_ControlPts
        weights(i-1)=weights(i)
    enddo
    weights(length_ControlPts)=weights(0)
    
    do i=1,numBasisFns
        do point=1,n
            call NURBSbasis(i,p,xi(point),length_knotVec,knotVec,weights,Rip,dRip)
            NURBSvalues(point,i)=Rip
            NURBSderivVals(point,i)=dRip
            r(point,i)=dabs(CollocPts(i)-xi(point))
        enddo
    enddo
    deallocate(xi)
    
    ipr=9
    open(ipr,file='NURBSvalue.txt')
    do i=1,n
        do j=1,numBasisFns-1
            write(ipr,1008,advance='no') NURBSvalues(i,j)
        enddo
        write(ipr,1008) NURBSvalues(i,numBasisFns)
    enddo
1008 format(E14.5,1X)
    close(ipr)
    
    do point=1,numBasisFns
        call NURBSinterpolation(CollocPts(point),p,length_knotVec,knotVec,&
            ControlPts(:,1),weights,Rip,dRip)
        CollocCoords(point,1)=Rip
        call NURBSinterpolation(CollocPts(point),p,length_knotVec,knotVec,&
            controlPts(:,2),weights,Rip,dRip)
        CollocCoords(point,2)=Rip
    enddo
    
    allocate(NURBSCoords(n,2))
    do i=1,n
        NURBSCoords(i,1)=0.0
        NURBSCoords(i,2)=0.0
        do j=1,numBasisFns
            NURBSCoords(i,1)=NURBSCoords(i,1)+NURBSvalues(i,j)*ControlPts(j,1)
            NURBSCoords(i,2)=NURBSCoords(i,2)+NURBSvalues(i,j)*ControlPts(j,2)
        enddo
    enddo
    
    length_ControlPts=length_ControlPts-1
    nPts=length_ControlPts
    nDof=2*nPts
    length_dispDofs=0
    do i=1,2*DispConn(ne,p)
        dispDofs(i)=i
        length_dispDofs=length_dispDofs+1
    enddo
    length_tracDofs=0
    do i=1,2*tracConn(ne,p+1)
        tracDofs(i)=i
        length_tracDofs=length_tracDofs+1
    enddo
    tracNdof=2*tracConn(ne,p+1)
    
    !****************
    !***输出数据以便检查
    ipr=6
    open(ipr,file='Output.txt')
    write(ipr,*)
    write(ipr,*) "******************************************"
    write(ipr,*) "-------Information of ControlPoints"
    write(ipr,*) "******************************************"
    write(ipr,100) 
100 format('Order',4X,'X',5X,'Y',5X,'Weight')    
    do i=1,nPts
        write(ipr,101) i,ControlPts(i,1),ControlPts(i,2),weights(i-1)
    enddo

    
    write(ipr,*)
    write(ipr,*) "******************************************"
    write(ipr,*) "-------Information of CollocCoords and CollocPts"
    write(ipr,*) "******************************************"
    write(ipr,1001)
1001 format('Order',4X,'CollocX',5X,'CollocY',5X,'CollocPts')    
    do i=1,nPts
        write(ipr,101) i,CollocCoords(i,1),CollocCoords(i,2),CollocPts(i)
    enddo
101 format(I5,1X,3(E12.3,1X))
    
    write(ipr,*)
    write(ipr,*) "******************************************"
    write(ipr,*) "-------Information of elRange"
    write(ipr,*) "******************************************"
    write(ipr,1002)
1002 format('Element',2X,'Begin',9X,'End')
     do e=1,ne
         write(ipr,102) e,elRange(e,1),elRange(e,2)
     enddo
102 format(I5,2X,2(E12.3,1X))     
    
    write(ipr,*) "******************************************"
    write(ipr,*) "-------Information of DispConn"
    write(ipr,*) "******************************************"
    do i=1,ne
        write(ipr,103) (DispConn(i,j),j=1,p+1)
    enddo
    
    write(ipr,*)
    write(ipr,*) "******************************************"
    write(ipr,*) "-------Information of BsFnConn"
    write(ipr,*) "******************************************"
    do i=1,ne
        write(ipr,103) (BsFnConn(i,j),j=1,p+1)
    enddo
    
    write(ipr,*)
    write(ipr,*) "******************************************"
    write(ipr,*) "-------Information of TracConn"
    write(ipr,*) "******************************************"
    do i=1,ne
        write(ipr,103) (TracConn(i,j),j=1,p+1)
    enddo
103 format(3(I5,1X))  
    
    write(ipr,*)
    write(ipr,*) "******************************************"
    write(ipr,*) "-------Information of TracConn"
    write(ipr,*) "******************************************"
    
    do i=0,length_knotVec-1
        write(ipr,104) i+1,knotVec(i)
    enddo
104 format('i=',I5,3X,'Knot=',E14.3)    
    close(ipr)
    
    open(ipr,file='NURBS.txt')
    do i=1,n
        write(ipr,105) NURBSCoords(i,1),NURBSCoords(i,2)
    enddo
105 format(2(E14.5,1X))   
    close(ipr)
    
    end
    
    
    
    
    
            
        
        
    
        
                
            
            
            
                
                    
                
                
            
            
            
        
        
            
    
    
    
    