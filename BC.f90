subroutine BC(n,presDispDOFs,length_presDispDOFs,&
        presTracDOFs,length_presTracDOFs,&
        nonZeroXTracDOFs,length_nonZeroXTracDOFs,&
        nonZeroYTracDOFs,length_nonZeroYTracDOFs)
    use GlobalPara
    implicit none
    integer(8)::n
    integer(8)::presDispDOFs(10000), presTracDOFs(10000)
    integer(8)::length_presDispDOFs,length_presTracDOFs
    integer(8)::nonZeroXTracDOFs(10000),nonZeroYTracDOFs(10000)
    integer(8)::length_nonZeroXTracDOFs,length_nonZeroYtracDOFs
    !**********************************************
    integer(8)::elms,i,j,str
    integer(8),allocatable::sctrX(:,:),sctrY(:,:),temp(:)
    integer(8)::length_sctrX,length_sctrY
    
    !*********************************************
    !****二维平板问题
    !****底边固定，顶端沿Y轴方向均布拉伸载荷
    !*********************************************
    
    !***第一条边的Ux,Uy均固定为0
    !取出Ux=0的节点编号
    allocate(sctrX(n,p+1))
    do elms=1,n
        do i=1,p+1
            sctrX(elms,i)=dispConn(elms,i)
        enddo
    enddo
    !进行重排
    allocate(temp(n*(p+1)))
    do i=1,n
        do j=1,p+1
            temp((i-1)*(p+1)+j)=sctrX(i,j)
        enddo
    enddo
    !排序
    do i=1,n*(p+1)-1
        do j=1,n*(p+1)-i
            if(temp(j)>temp(j+1)) then
                str=temp(j)
                temp(j)=temp(j+1)
                temp(j+1)=str
            endif
        enddo
    enddo
    !取出unique元素
    deallocate(sctrX)
    allocate(sctrX(n*(p+1),1))
    sctrX(1,1)=temp(1)
    j=1
    do i=2,n*(p+1)
        if(temp(i).ne.temp(i-1)) then
            j=j+1
            sctrX(j,1)=temp(i)
        endif
    enddo
    length_sctrX=j
    deallocate(temp)
    !******************
    !取出Uy=0的节点编号
    allocate(sctrY(n,p+1))
    do elms=1,n
        do j=1,p+1
            sctrY(elms,j)=dispConn(elms,j)
        enddo
    enddo
    allocate(temp(n*(p+1)))
    do i=1,n
        do j=1,p+1
            temp((i-1)*(p+1)+j)=sctrY(i,j)
        enddo
    enddo
    do i=1,n*(p+1)-1
        do j=1,n*(p+1)-i
            if(temp(j)>temp(j+1)) then
                str=temp(j)
                temp(j)=temp(j+1)
                temp(j+1)=str
            endif
        enddo
    enddo
    deallocate(sctrY)
    allocate(sctrY(n*(p+1),1))
    sctrY(1,1)=temp(1)
    j=1
    do i=2,n*(p+1)
        if(temp(i).ne.temp(i-1)) then
            j=j+1
            sctrY(j,1)=temp(i)
        endif
    enddo
    length_sctrY=j
    deallocate(temp)
    !*************************
    do i=1,length_sctrX
        presDispDOFs(i)=2*sctrX(i,1)-1
    enddo
    do i=1,length_sctrY
        presDispDOFs(i+length_sctrX)=2*sctrY(i,1)
    enddo
    length_presDispDOFs=length_sctrX+length_sctrY
    !排序
    do i=1,length_presDispDOFs-1
        do j=1,length_presDispDOFs-i
            if(presDispDOFs(j)>presDispDOFs(j+1)) then
                str=presDispDOFs(j)
                presDispDOFs(j)=presDispDOFs(j+1)
                presDispDOFs(j+1)=str
            endif
        enddo
    enddo
    deallocate(sctrX,sctrY)
    !**************************
    !***提取面力Tx自由度编号
    allocate(sctrX(3*n,p+1))
    do elms=(n+1),4*n
        do j=1,p+1
            sctrX(elms-n,j)=tracConn(elms,j)
        enddo
    enddo
    allocate(temp(3*n*(p+1)))
    do i=1,3*n
        do j=1,p+1
            temp((i-1)*(p+1)+j)=sctrX(i,j)
        enddo
    enddo
    !排序
    do i=1,3*n*(p+1)-1
        do j=1,3*n*(p+1)-i
            if(temp(j)>temp(j+1)) then
                str=temp(j)
                temp(j)=temp(j+1)
                temp(j+1)=str
            endif
        enddo
    enddo
    deallocate(sctrX)
    allocate(sctrX(3*n*(p+1),1))
    sctrX(1,1)=temp(1)
    j=1
    do i=2,3*n*(p+1)
        if(temp(i).ne.temp(i-1)) then
            j=j+1
            sctrX(j,1)=temp(i)
        endif
    enddo
    length_sctrX=j
    deallocate(temp)
    !**********************
    !***提取面力Ty的自由度编号
    allocate(sctrY(3*n,p+1))
    do elms=n+1,4*n
        do j=1,p+1
            sctrY(elms-n,j)=tracConn(elms,j)
        enddo
    enddo
    allocate(temp(3*n*(p+1)))
    do i=1,3*n
        do j=1,p+1
            temp((i-1)*(p+1)+j)=sctrY(i,j)
        enddo
    enddo
    
    !排序
    do i=1,3*n*(p+1)-1
        do j=1,3*n*(p+1)-i
            if(temp(j)>temp(j+1)) then
                str=temp(j)
                temp(j)=temp(j+1)
                temp(j+1)=str
            endif
        enddo
    enddo
    deallocate(sctrY)
    allocate(sctrY(3*n*(p+1),1))
    sctrY(1,1)=temp(1)
    j=1
    do i=2,3*n*(p+1)
        if(temp(i).ne.temp(i-1)) then
            j=j+1
            sctrY(j,1)=temp(i)
        endif
    enddo
    length_sctrY=j
    deallocate(temp)
    
    length_presTracDOFs=length_sctrX+length_sctrY
    do i=1,length_sctrX
        presTracDOFs(i)=2*sctrX(i,1)-1
    enddo
    do i=1,length_sctrY
        presTracDOFs(i+length_sctrX)=2*sctrY(i,1)
    enddo
    !排序
    do i=1,length_presTracDOFs-1
        do j=1,length_presTracDOFs-i
            if(presTracDOFs(j)>presTracDOFs(j+1)) then
                str=presTracDOFs(j)
                presTracDOFs(j)=presTracDOFs(j+1)
                presTracDOFs(j+1)=str
            endif
        enddo
    enddo
    deallocate(sctrX,sctrY)
    !********************************
    !提取出非零的面力自由度
    !Tx非零的编号没有
    length_nonZeroXTracDOFs=0
    
    !************************
    !***Ty非零的面力自由度编号
    allocate(sctrY(n,(p+1)))
    do elms=2*n+1,3*n
        do j=1,p+1
            sctrY(elms-2*n,j)=tracConn(elms,j)
        enddo
    enddo
    allocate(temp(n*(p+1)))
    do i=1,n
        do j=1,p+1
            temp((i-1)*(p+1)+j)=sctrY(i,j)
        enddo
    enddo
    do i=1,n*(p+1)-1
        do j=1,n*(p+1)-i
            if(temp(j)>temp(j+1)) then
                str=temp(j)
                temp(j)=temp(j+1)
                temp(j+1)=str
            endif
        enddo
    enddo
    deallocate(sctrY)
    allocate(sctrY(n*(p+1),1))
    sctrY(1,1)=temp(1)
    j=1
    do i=2,n*(p+1)
        if(temp(i).ne.temp(i-1)) then
            j=j+1
            sctrY(j,1)=temp(i)
        endif
    enddo
    length_nonZeroYTracDOFs=j
    do j=1,length_nonZeroYTracDOFs
        nonZeroYTracDOFs(j)=2*sctrY(j,1)
    enddo
    
    deallocate(sctrY,temp)
        
     
    !*********************************************
    !****平板左右拉伸问题
    !****所有位移位置，只给了面力，这种情况
    !****需要在Main函数中处理刚体位移
    !**********************************************
     
    end