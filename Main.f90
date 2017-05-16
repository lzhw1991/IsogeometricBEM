! 等几何边界元-线弹性力学问题
! 作者：飞狐
! Bug report: walkandthinker@gmail.com
! 讨论:计算力学QQ 467597866
Program IsoBEM
    use MaterialPara
    use GlobalPara
    use Solver
    use NURBS
    implicit none
    
    integer(8)::i,j,flag
    
    integer(8)::refinement
    integer(8)::presDispDOFs(10000),length_presDispDOFs
    integer(8)::presTracDOFs(10000),length_presTracDOFs
    integer(8)::nonZeroXTracDOFs(10000),length_nonZeroXTracDOFs
    integer(8)::nonZeroYTracDOFs(10000),length_nonZeroYTracDOFs
    
    integer(8)::unknownDispDofs(10000),unknownTracDofs(10000)
    integer(8)::length_unknownDispDofs,length_unknownTracDofs
    
    real(8)::pi
    
    real(8),allocatable::A(:,:),H(:,:),G(:,:),Z(:)
    real(8)::collocNormals(10000,2,2)
    integer(8)::ngp_s,ngp_r     !高斯积分点
    
    integer(8)::c,element,glbBsFnConn(3),dElConn(3),tElConn(3)
    real(8)::srcXi_param,range(2),elCoords(3,2),collocGlbPt(2)
    real(8)::n1(2),n2(2),jumpTerm(2,2),Hsub(2,6),Gsub(2,6)
    integer(8)::sctrVec(6),rowSctrVec(2)
    
    integer(8)::tracDispConnDOF(10000),mappedTractionDofs(10000)
    integer(8)::length_mappedTractionDofs,L
    real(8),allocatable::displacement(:),soln(:),traction(:)
    real(8),allocatable::disp(:),trac(:)
    
    integer(8)::inp,ipr        !输出结果
    integer(8)::Pa,Pb
    
    real(8)::start,finish
    
    real(8)::XY(2),stress(2,2),Rip(3),dRip(3)
    
    call CPU_TIME(start)
    call About()
    pi=4.d0*atan(1.d0)
    
    refinement=1  !单元细化数，越大自由度数越多
    ngp_s=50
    ngp_r=50 !ngp为高斯点积分的高斯点个数
    
    !读入控制点的坐标和权值
    !读入材料参数以及内点坐标
    call InputData()
    
    !生成等几何的配置点和自由度编号等信息
    call generateBEMmesh(refinement)
    
    !施加边界条件
    !这里只计算了下端固定，上端均布拉伸的均质平板算例
    !施加的是沿Y轴正方向的单位1的均布载荷
    call BC(refinement,presDispDOFs,length_presDispDOFs,presTracDOFs,&
        length_presTracDOFs,nonZeroXTracDOFs,length_nonZeroXTracDOFs,&
        nonZeroYTracDOFs,length_nonZeroYTracDOFs)

    write(*,*) "nDof=",nDof
    
    !得到未知的位移自由度编号
    length_unknownDispDofs=0
    do i=1,length_dispDofs
        flag=0
        do j=1,length_presDispDOFs
            if(dispDofs(i).eq.presDispDOFs(j)) then
                flag=flag+1
            endif
        enddo
        if(flag.eq.0) then
            length_unknownDispDofs=length_unknownDispDofs+1
            unknownDispDofs(length_unknownDispDofs)=dispDofs(i)
        endif 
    enddo
    
    !得到未知的面力自由度编号
    length_unknownTracDofs=0
    do i=1,length_tracDofs
        flag=0
        do j=1,length_presTracDOFs
            if(tracDofs(i).eq.presTracDOFs(j)) then
                flag=flag+1
            endif
        enddo
        if(flag.eq.0) then
            length_unknownTracDofs=length_unknownTracDofs+1
            unknownTracDofs(length_unknownTracDofs)=tracDofs(i)
        endif
    enddo
   
    allocate(H(nDof,nDof),A(nDof,nDof),G(nDof,tracNdof),Z(nDof))
    do i=1,nDof
        do j=1,nDof
            H(i,j)=0.d0
            A(i,j)=0.d0
        enddo
        do j=1,tracNdof
            G(i,j)=0.d0
        enddo
        Z(i)=0.d0
    enddo
    
    call findNormals(collocNormals) !得到配置点的前后法向量
    
    do c=1,nPts  !对源点进行循环
        srcXi_param=collocPts(c) !节点在参数空间knotVec中对应的值
        do element=1,ne     !对单元循环
            
            range(1)=elRange(element,1) !单元的knotVec值的范围
            range(2)=elRange(element,2)
            
            do j=1,3
                glbBsFnConn(j)=bsFnConn(element,j) !单元的插值函数编号
                dElConn(j)=dispConn(element,j)     !单元的位移自由度编号
                tElConn(j)=tracConn(element,j)     !单元的面力自由度编号
            enddo
            do j=1,3
                elCoords(j,1)=controlPts(dElConn(j),1) !单元上的p+1个控制点的坐标
                elCoords(j,2)=controlPts(dElConn(j),2)
            enddo
            do j=1,2
                collocGlbPt(j)=collocCoords(c,j)  !源点的坐标
            enddo
            
            if((srcXi_param<=range(2)).and.((srcXi_param>=range(1)).or.((element==ne).and.(srcXi_param==0.d0)))) then
                !计算奇异积分
                
                !当第1个点处在最后一个单元的末端是，值要赋为1.0
                if((element==ne).and.(srcXi_param==0.d0)) then
                    srcXi_param=1.d0
                endif
                
                n1(1)=collocNormals(c,1,1);n1(2)=collocNormals(c,2,1)
                n2(1)=collocNormals(c,1,2);n2(2)=collocNormals(c,2,2)
                
                call CalculateJumpTerm(n1,n2,jumpTerm) !计算CU+HU=GT中的C
                
                
                call integrateHsubmatrixSST(ngp_s,elcoords,&
                    glbBsFnConn,collocGlbPt,srcXi_param,range,jumpTerm,Hsub)
                
                call integrateGsubmatrix_Tells(ngp_s,elcoords,&
                    glbBsFnConn,collocGlbPt,srcXi_param,range,Gsub)
                
                
            else
                
                !计算非奇异积分
                call integrateHGsubmatrices_GLQ(ngp_r,elcoords,&
                    glbBsFnConn,collocGlbPt,range,Hsub,Gsub)
                    
            endif
            
            !装配系数矩阵
            do j=1,2
                rowSctrVec(j)=2*(c-1)+j
            enddo
            do j=1,p+1
                sctrVec(2*j-1)=2*dElConn(j)-1
                sctrVec(2*j)=2*dElConn(j)
            enddo
            do i=1,2
                do j=1,6
                    H(rowSctrVec(i),sctrVec(j))=H(rowSctrVec(i),sctrVec(j))+Hsub(i,j)
                enddo
            enddo
            !*************************
            do j=1,p+1
                sctrVec(2*j-1)=2*tElConn(j)-1
                sctrVec(2*j)=2*tElConn(j)
            enddo
            do i=1,2
                do j=1,6
                    G(rowSctrVec(i),sctrVec(j))=G(rowSctrVec(i),sctrVec(j))+Gsub(i,j)
                enddo
            enddo
        enddo
    enddo
    !**************************************
    !*****装配系数矩阵
    
    !将所有位置的位移自由度对应的H矩阵的系数装配到A矩阵中
    do i=1,length_unknownDispDofs
        do j=1,nDof
            A(j,unknownDispDofs(i))=H(j,unknownDispDofs(i))
        enddo
    enddo
    
    do i=1,length_tracDispConn
        tracDispConnDOF(2*i-1)=2*tracDispConn(i)-1
        tracDispConnDOF(2*i)=2*tracDispConn(i)
    enddo
    !得到位置面力自由度对应的与H矩阵关联的编号
    do i=1,length_unknownTracDofs
        mappedTractionDofs(i)=tracDispConnDOF(unknownTracDofs(i))
    enddo
    length_mappedTractionDofs=length_unknownTracDofs
    
    !将位置的面力自由度对应的G矩阵的系数装配到A矩阵中
    do i=1,length_mappedTractionDofs
        do j=1,nDof
            A(j,mappedTractionDofs(i))=A(j,mappedTractionDofs(i))-G(j,unknownTracDofs(i))
        enddo
    enddo
    !************************
    !*****形成右端项
    do i=1,nDof
        do j=1,length_nonZeroXTracDOFs
            Z(i)=Z(i)+G(i,nonZeroXTracDOFs(j))*(0.d0)
        enddo
    enddo
    do i=1,nDof
        do j=1,length_nonZeroYTracDOFs
            Z(i)=Z(i)+G(i,nonZeroYTracDOFs(j))*(1.d0) !施加单位1的均布载荷
        enddo
    enddo
    
    ipr=5
    open(ipr,file='check.txt')
    write(ipr,*) "The Matrix A is:"
    do i=1,nDof
        do j=1,nDof-1
            write(ipr,1001,advance='no') A(i,j)
        enddo
        write(ipr,1001) A(i,nDof)
    enddo
    
    write(ipr,*) "The Matrix H is:"
    do i=1,nDof
        do j=1,nDof-1
            write(ipr,1001,advance='no') H(i,j)
        enddo
        write(ipr,1001) H(i,nDof)
    enddo
    
    write(ipr,*) "The Matrix G is:"
    do i=1,nDof
        do j=1,tracNdof-1
            write(ipr,1001,advance='no') G(i,j)
        enddo
        write(ipr,1001) G(i,tracNdof)
    enddo
    
    write(ipr,*) "The Matrix Z is:"
    do i=1,nDof
        write(ipr,1001) Z(i)
    enddo
1001 format(E14.6,1X)   
     close(ipr)
    
     
    allocate(displacement(nDof),soln(nDof),traction(tracNdof))
    allocate(disp(nDof),trac(tracNdof))
    do i=1,nDof
        displacement(i)=0.d0
        disp(i)=0.d0
        soln(i)=0.d0
    enddo
    do i=1,tracNdof
        traction(i)=0.d0
        trac(i)=0.d0
    enddo
    
    call AGAUSS(A,Z,nDof,soln,L)
    !call AGGJE(A,Z,nDof,soln,L)
    !call BEM_Solver(A,Z,soln,nDof,1.0*L)
    if(L==0) then
        write(*,*) "方程奇异，求解失败!"
    else
        do i=1,length_presDispDOFs
            displacement(presDispDOFs(i))=0.d0
        enddo
        do i=1,length_unknownDispDofs
            displacement(unknownDispDofs(i))=soln(unknownDispDofs(i))
        enddo
        
        do i=1,length_nonZeroXTracDOFs
            traction(nonZeroXTracDOFs(i))=0.d0
        enddo
        do i=1,length_nonZeroYTracDOFs
            traction(nonZeroYTracDOFs(i))=1.d0
        enddo
        do i=1,length_unknownTracDofs
            traction(unknownTracDofs(i))=soln(mappedTractionDofs(i))
        enddo
        
        inp=6
        open(inp,file='result.txt')
        write(inp,100)
100     format('Order',5X,'X',12X,'Y',12X,'Ux',12X,'Uy')   
        do i=1,nPts
            write(inp,101) i,CollocCoords(i,1),CollocCoords(i,2),&
                displacement(2*i-1),displacement(2*i)
101             format(I5,1X,2(E10.3,1X),2(E13.6,1X))
        enddo
        
        write(inp,102) 
102     format('Order',4X,'X',12X,'Y',10X,'Tax',10X,'Tbx',12X,'Tay',10X,'Tby')
        do i=1,nPts
            do element=1,ne
                do j=1,p+1
                    if(i==dispConn(element,j)) then
                        if(j==1) then
                            Pb=tracConn(element,j)
                            if(element==1) then
                                Pa=tracConn(ne,p+1)
                            else
                                Pa=tracConn(element-1,p+1)
                            endif
                        elseif(j==2) then
                            Pa=tracConn(element,j)
                            Pb=tracConn(element,j)
                        elseif(j==p+1) then
                            Pa=tracConn(element,j)
                            if(element==ne) then
                                Pb=tracConn(1,1)
                            else
                                Pb=tracConn(element+1,1)
                            endif
                        endif
                    endif
                enddo
            enddo
            
            write(inp,103) i,CollocCoords(i,1),CollocCoords(i,2),&
                traction(2*Pa-1),traction(2*Pb-1),traction(2*Pa),traction(2*Pb)
103             format(I4,1X,2(E10.3,1X),4(E13.6,1X))
        enddo
        
        do element=1,ne
            range(1)=elRange(element,1) !单元的knotVec值的范围
            range(2)=elRange(element,2)

            do j=1,3
                glbBsFnConn(j)=bsFnConn(element,j) !单元的插值函数编号
                dElConn(j)=dispConn(element,j)     !单元的位移自由度编号
                tElConn(j)=tracConn(element,j)     !单元的面力自由度编号
            enddo
            do c=1,p+1
                if(c==1) then
                    srcXi_param=range(1)+2.226e-16
                elseif(c==2) then
                    srcXi_param=collocPts(dElConn(c))
                elseif(c==p+1) then
                    srcXi_param=range(2)-2.226e-16
                endif                    
                do j=1,p+1
                    i=glbBsFnConn(j)
                    call NURBSbasis(i,p,srcXi_param,length_knotVec,knotVec,weights,Rip(j),dRip(j))
                    trac(2*tElConn(c)-1)=trac(2*tElConn(c)-1)+Rip(j)*traction(2*tElConn(j)-1)
                    trac(2*tElConn(c))=trac(2*tElConn(c))+Rip(j)*traction(2*tElConn(j))
                    disp(2*dElConn(c)-1)=disp(2*dElConn(c)-1)+Rip(j)*displacement(2*dElConn(j)-1)
                    disp(2*dElConn(c))=disp(2*dElConn(c))+Rip(j)*displacement(2*dElConn(j))
                enddo
            enddo
        enddo
        !******************************
        do c=1,nPts
            do element=1,ne
                if(c==dispConn(element,1)) then
                    disp(2*c-1)=disp(2*c-1)*0.5d0
                    disp(2*c)=disp(2*c)*0.5d0
                endif
            enddo
        enddo
        write(inp,*) "********** The Real Displacements **********"
        write(inp,100)
        do i=1,nPts
            write(inp,101) i,CollocCoords(i,1),CollocCoords(i,2),&
                disp(2*i-1),disp(2*i)
        enddo
        write(inp,*) "********** The Real Tractions **********"
        write(inp,102) 
        do i=1,nPts
            do element=1,ne
                do j=1,p+1
                    if(i==dispConn(element,j)) then
                        if(j==1) then
                            Pb=tracConn(element,j)
                            if(element==1) then
                                Pa=tracConn(ne,p+1)
                            else
                                Pa=tracConn(element-1,p+1)
                            endif
                        elseif(j==2) then
                            Pa=tracConn(element,j)
                            Pb=tracConn(element,j)
                        elseif(j==p+1) then
                            Pa=tracConn(element,j)
                            if(element==ne) then
                                Pb=tracConn(1,1)
                            else
                                Pb=tracConn(element+1,1)
                            endif
                        endif
                    endif
                enddo
            enddo
            
            write(inp,103) i,CollocCoords(i,1),CollocCoords(i,2),&
                trac(2*Pa-1),trac(2*Pb-1),trac(2*Pa),trac(2*Pb)
        enddo
        
        
        
        !*******************************
        !***计算内点应力
        write(inp,*) "********** The Stress of the Inner points**********"
        write(inp,104)
104     format('Order',5X,'X',9X,'Y',9X,'Sigma_X',8X,'Sigma_Y',6X,'Sigma_XY')
        do i=1,Num_InnerPts
            XY(1)=cx(i);XY(2)=cy(i)
            call findInternalStress(displacement,traction,XY,stress)
            write(inp,105) i,cx(i),cy(i),stress(1,1),stress(2,2),stress(1,2)
105         format(I4,1X,2(E10.3,1X),3(E13.6,1X))
        enddo  
        
        
        write(inp,*)
        write(inp,*) "**********Finished!**********"
        close(inp)
    endif
    call CPU_TIME(finish)
    write(*,*) "****************************"
    write(*,*) "*****计算用时：",finish-start
    write(*,*) "****************************"
    stop
    end