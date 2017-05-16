    subroutine convertToParamSpace(xi,range,xi_param)
    implicit none
    real(8)::xi,range(2)
    real(8)::xi_param

    xi_param=((range(2)-range(1))*xi+(range(2)+range(1)))/2.0
    end

    subroutine convertToParentCoordSpace(xi_param,range,xi)
    implicit none
    real(8)::xi_param,range(2)
    real(8)::xi
    xi=(2*xi_param-(range(2)+range(1)))/(range(2)-range(1))
    
    end