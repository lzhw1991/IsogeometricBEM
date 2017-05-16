Module GlobalPara
    implicit none
    integer(8)::p
    real(8)::const1,const2,const3,const4,kappa
    real(8)::iconst1,iconst2,iconst3,iconst4,ikappa
    
    real(8)::knotVec(0:10000),iknotVec(0:10000),incknotVec(0:10000)
    real(8)::weights(0:10000),iweights(0:10000),incweights(0:10000)
    integer(8)::length_knotVec,length_iknotVec,length_incknotVec
    integer(8)::nDof,ne,nPts
    integer(8)::dispDofs(20000),tracDofs(30000),tracNdof
    integer(8)::length_dispDofs,length_tracDofs
    integer(8)::length_controlPts
    integer(8)::mDof,mne,mPts
    integer(8)::inDof,ine,iPts
    integer(8)::bsFnConn(10000,4),dispConn(10000,4),tracConn(10000,4),tracDispConn(10000)
    integer(8)::length_tracDispConn
    real(8)::controlPts(10000,4),collocPts(10000),collocCoords(10000,4),elRange(10000,4)
    
    !***********************
    !***内点的坐标信息
    real(8)::cx(5000),cy(5000)
    integer(8)::Num_InnerPts
    
    
    
    end Module GlobalPara