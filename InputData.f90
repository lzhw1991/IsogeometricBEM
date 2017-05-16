subroutine InputData()
    use GlobalPara
    use MaterialPara
    implicit none
    !*******************
    !***读取材料参数并对相关系数进行初始化
    !***读取KnotVec和ControlPts的信息
    character(len=80)::str
    integer(8)::inp,m,n
    integer(8)::i,j
    real(8)::pi
    inp=2
    open(inp,file='input.txt')
    
    read(inp,*) str
    read(inp,*) E,mu,iE,imu
    
    read(inp,*) str
    read(inp,*) length_knotVec,length_ControlPts,p,num_InnerPts
    
    read(inp,*) str
    do i=1,length_knotVec
        read(inp,*) knotVec(i)
    enddo
    
    read(inp,*) str
    do i=1,length_ControlPts
        read(inp,*) controlPts(i,1),controlPts(i,2),weights(i)
    enddo
    
    read(inp,*) str
    do i=1,Num_InnerPts
        read(inp,*) cx(i),cy(i)
    enddo
    
    
    !平面应力问题
    E=(1.d0+mu)*E/((1.d0+mu)**2)
    mu=mu/(1.d0+mu)
    
    !初始化材料参数和相关变量
    pi=4.d0*datan(1.d0)
    shearMod=E/(2.d0*(1.d0+mu))
    const4=(1.d0-2.d0*mu)
    const3=1.d0/(4.d0*pi*(1.d0-mu))
    const2=(3.d0-4.d0*mu)
    const1=1.d0/(8.d0*pi*shearMod*(1.d0-mu))
    
    end    