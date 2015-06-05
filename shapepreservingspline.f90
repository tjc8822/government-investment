                                             !********************************************************************
! MODULE DSPSP (version 2.0, October 25, 2003)
! Quadratic shape preserving spline proposed by L. Schumaker (1983).
! The algorithm is from "Numerical Method in Economics" by K.
! Judd (1998) pp231-235. If the data implies monotonicity or
! concavity, these shape is passed on to the interpolated function.
! F90 code.
! Double precision.
! DSPSPLINE outputs the necessary data for interpolation, using
! grids and value of the function at the grdis points as inputs.
! The output is the additional breaking points, one for each
! interval.
! DSPSPLINT evaluates the interpolated function at a given point,
! using the output of DSPSPLINE.
! DSPSPLEASY does both of them in one subroutine.
!********************************************************************
MODULE DSPSP
  implicit none
  private
  public:: dspspline, dspsplint, dspspleasy, dspsplintd, dspspleasyd
  public:: dspspline2, dspspleasy2, dspspleasyd2
  public:: dspspline3, dspspleasy3, dspspleasyd3

CONTAINS

  !********************************************************************
  ! SUBROUTINE DSPSPLINE
  !********************************************************************
  SUBROUTINE DSPSPLINE(sp_np, sp_x, sp_y, sp_work, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    !number of grid points (>2)
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    !sp_x(1:sp_np) is the array of grid points
    !sp_y(1:sp_np) is the array of the value of functions

    !******* output variables *******
    real(sp_prec), intent(out), dimension(2*sp_np-1):: sp_work
    integer, intent(out):: sp_error
    !error code

    !******* local variables *******
    real(sp_prec), dimension(1:sp_np-1):: sp_zeta
    !array of additional breaking points
    real(sp_prec), dimension(1:sp_np):: sp_d1
    !sp_d1(1:sp_np) is the array to store 1st derivatives
    real(sp_prec), dimension(1:sp_np-1):: sp_delta
    !internally used

    !******* check inputs *******
    call sp_check_inputs1(sp_np, sp_x, sp_y, sp_error)
    if (sp_error>0) return

    !******* compute d1 and delta *******
    call sp_compt_d1(sp_np, sp_x, sp_y, sp_d1, sp_delta)

    !******* compute additional break points *******
    call sp_compt_zeta(sp_np, sp_x, sp_y, sp_d1, sp_delta, sp_zeta)

    !******* construct output array *******
    sp_work(1:sp_np-1)=sp_zeta(1:sp_np-1)
    sp_work(sp_np:2*sp_np-1)=sp_d1(1:sp_np)

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLINE

  !********************************************************************
  ! SUBROUTINE DSPSPLINE2
  !********************************************************************
  SUBROUTINE DSPSPLINE2(sp_np, sp_x, sp_y, flag_1, sp_d1_1, &
       flag_n, sp_d1_n, sp_work, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np, flag_1, flag_n
    !number of grid points (>2)
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    !sp_x(1:sp_np) is the array of grid points
    !sp_y(1:sp_np) is the array of the value of functions

    !******* input variables (special for DSPSPLINE2) *******
    real(sp_prec), intent(in):: sp_d1_1, sp_d1_n
    ! sp_d1 is the derivative at point 1
    ! sp_dn is the derivative at point sp_np

    !******* output variables *******
    real(sp_prec), intent(out), dimension(2*sp_np-1):: sp_work
    integer, intent(out):: sp_error
    !error code

    !******* local variables *******
    real(sp_prec), dimension(1:sp_np-1):: sp_zeta
    !array of additional breaking points
    real(sp_prec), dimension(1:sp_np):: sp_d1
    !sp_d1(1:sp_np) is the array to store 1st derivatives
    real(sp_prec), dimension(1:sp_np-1):: sp_delta
    !internally used

    !******* check inputs *******
    call sp_check_inputs1(sp_np, sp_x, sp_y, sp_error)
    if (sp_error>0) return

    !******* compute d1 and delta *******
    call sp_compt_d1(sp_np, sp_x, sp_y, sp_d1, sp_delta)

    !******* replace sp_d1(1) and sp_d1(sp_np) *******
    if (flag_1/=0) sp_d1(1)=sp_d1_1
    if (flag_n/=0) sp_d1(sp_np)=sp_d1_n

    !******* compute additional break points *******
    call sp_compt_zeta(sp_np, sp_x, sp_y, sp_d1, sp_delta, sp_zeta)

    !******* construct output array *******
    sp_work(1:sp_np-1)=sp_zeta(1:sp_np-1)
    sp_work(sp_np:2*sp_np-1)=sp_d1(1:sp_np)

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLINE2

  !********************************************************************
  ! SUBROUTINE DSPSPLINE3
  !********************************************************************
  SUBROUTINE DSPSPLINE3(sp_np, sp_x, sp_y, sp_d, sp_work, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    !number of grid points (>2)
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    !sp_x(1:sp_np) is the array of grid points
    !sp_y(1:sp_np) is the array of the value of functions

    !******* input variables (special for DSPSPLINE3) *******
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_d
    !sp_d(1:sp_np) is the array of derivatives at grid points

    !******* output variables *******
    real(sp_prec), intent(out), dimension(2*sp_np-1):: sp_work
    integer, intent(out):: sp_error
    !error code

    !******* local variables *******
    real(sp_prec), dimension(1:sp_np-1):: sp_zeta
    !array of additional breaking points
    real(sp_prec), dimension(1:sp_np):: sp_d1
    !sp_d1(1:sp_np) is the array to store 1st derivatives
    real(sp_prec), dimension(1:sp_np-1):: sp_delta
    !internally used

    !******* check inputs *******
    call sp_check_inputs1(sp_np, sp_x, sp_y, sp_error)
    if (sp_error>0) return

    !******* compute d1 and delta *******
    call sp_compt_d1(sp_np, sp_x, sp_y, sp_d1, sp_delta)

    !******* replace sp_d1 with sp_d *******
    sp_d1=sp_d

    !******* compute additional break points *******
    call sp_compt_zeta(sp_np, sp_x, sp_y, sp_d1, sp_delta, sp_zeta)

    !******* construct output array *******
    sp_work(1:sp_np-1)=sp_zeta(1:sp_np-1)
    sp_work(sp_np:2*sp_np-1)=sp_d1(1:sp_np)

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLINE3

  !********************************************************************
  ! SUBROUTINE SP_CHECK_INPUTS1
  !********************************************************************
  SUBROUTINE SP_CHECK_INPUTS1(sp_np, sp_x, sp_y, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y

    !******* output variables *******
    integer, intent(out):: sp_error

    !******* local variables *******
    integer:: p0

    !******* set default error code *******
    sp_error=0

    !******* check number of sp_np *******
    if (sp_np<3) then
       sp_error=1
       return
    end if

    !******* check order of grids *******
    do p0=1, sp_np-1
       if (sp_x(p0)>=sp_x(p0+1)) then
          sp_error=2
          return
       end if
    end do

    !******* end of subroutine *******
    return

  END SUBROUTINE SP_CHECK_INPUTS1

  !********************************************************************
  ! SUBROUTINE SP_COMPT_D1
  !********************************************************************
  SUBROUTINE SP_COMPT_D1(sp_np, sp_x, sp_y, sp_d1, sp_delta)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y

    !******* output variables *******
    real(sp_prec), intent(out), dimension(1:sp_np):: sp_d1
    real(sp_prec), intent(out), dimension(1:sp_np-1):: sp_delta

    !******* local variables *******
    integer:: p0
    real(sp_prec), dimension(1:sp_np-1):: sp_l

    !******* calculate l and delta *******
    do p0=1, sp_np-1
       sp_l(p0)=((sp_x(p0+1)-sp_x(p0))**2._sp_prec+&
            (sp_y(p0+1)-sp_y(p0))**2._sp_prec)**0.5_sp_prec
       sp_delta(p0)=(sp_y(p0+1)-sp_y(p0))/(sp_x(p0+1)-sp_x(p0))
    end do

    !******* compute sp_d1(2:sp_np-1) *******
    do p0=2, sp_np-1
       if ((sp_delta(p0-1)*sp_delta(p0))<=0._sp_prec) then
          sp_d1(p0)=0._sp_prec
       else
          sp_d1(p0)=((sp_l(p0-1)*sp_delta(p0-1))+(sp_l(p0)*sp_delta(p0)))&
               /(sp_l(p0-1)+sp_l(p0))
       end if
    end do

    !******* compute sp_d1(1) and sp_d1(sp_np) *******
    sp_d1(1)=(3._sp_prec*sp_delta(1)-sp_d1(2))/2._sp_prec
    sp_d1(sp_np)=(3._sp_prec*sp_delta(sp_np-1)-sp_d1(sp_np-1))&
         /2._sp_prec

    !******* end of subroutine *******
    return

  END SUBROUTINE SP_COMPT_D1

  !********************************************************************
  ! SUBROUTINE SP_COMPT_ZETA
  !********************************************************************
  SUBROUTINE sp_compt_zeta&
       (sp_np, sp_x, sp_y, sp_d1, sp_delta, sp_zeta)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y, sp_d1
    real(sp_prec), intent(in), dimension(1:sp_np-1):: sp_delta

    !******* output variables *******
    real(sp_prec), intent(out), dimension(1:sp_np-1):: sp_zeta

    !******* local variables *******
    integer:: p0
    real(sp_prec):: zetatemp, delta1, delta2

    !******* compute sp_zeta *******
    do p0=1, sp_np-1
       delta1=sp_d1(p0)-sp_delta(p0)
       delta2=sp_d1(p0+1)-sp_delta(p0)
       if (delta1*delta2>=0._sp_prec) then
          sp_zeta(p0)=(sp_x(p0)+sp_x(p0+1))/2._sp_prec
       else if (abs(delta1)>abs(delta2)) then
          zetatemp=sp_x(p0)+2._sp_prec*(sp_x(p0+1)-sp_x(p0))*delta2/&
               (sp_d1(p0+1)-sp_d1(p0))
          sp_zeta(p0)=(sp_x(p0)+zetatemp)/2._sp_prec
       else
          zetatemp=sp_x(p0+1)+2._sp_prec*(sp_x(p0+1)-sp_x(p0))*delta1/&
               (sp_d1(p0+1)-sp_d1(p0))
          sp_zeta(p0)=(zetatemp+sp_x(p0+1))/2._sp_prec
       end if
    end do

    !******* end of subroutine *******
    return

  END SUBROUTINE SP_COMPT_ZETA

  !********************************************************************
  ! SUBROUTINE DSPSPLINT
  !********************************************************************
  SUBROUTINE DSPSPLINT(sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    !number of grid points
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    !sp_x(1:sp_np) is the array of grid points
    !sp_y(1:sp_np) is the array of the value of functions
    real(sp_prec), intent(in), dimension(1:2*sp_np-1):: sp_work
    !array of additional breaking points
    real(sp_prec), intent(in):: sp_xv
    !point where the interpolated function is evaluated

    !******* output variables *******
    real(sp_prec), intent(out):: sp_yv
    !value of the interpolated function at sp_xv
    integer, intent(out):: sp_error
    !error code

    !******* local variables *******
    real(sp_prec), dimension(1:sp_np-1):: sp_zeta
    !array of additional breaking points
    real(sp_prec), dimension(1:sp_np):: sp_d1
    !1st derivatives
    real(sp_prec), dimension(1:sp_np-1):: sp_delta
    !internally used

    !******* construct zeta and d1 *******
    sp_zeta(1:sp_np-1)=sp_work(1:sp_np-1)
    sp_d1(1:sp_np)=sp_work(sp_np:2*sp_np-1)

    !******* check inputs *******
    sp_error=0
    call sp_check_inputs2(sp_np, sp_x, sp_y, sp_zeta, sp_xv, sp_error)
    if (sp_error>0) return

    !******* evaluate the interpolated function *******
    call sp_eval(sp_np, sp_x, sp_y, sp_d1, sp_zeta, sp_xv, sp_yv)

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLINT

  !********************************************************************
  ! SUBROUTINE DSPSPLINTD
  ! DSPSPLINT plus calculating derivative
  !********************************************************************
  SUBROUTINE DSPSPLINTD(sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_yd,&
       sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    !number of grid points
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    !sp_x(1:sp_np) is the array of grid points
    !sp_y(1:sp_np) is the array of the value of functions
    real(sp_prec), intent(in), dimension(1:2*sp_np-1):: sp_work
    !array of additional breaking points
    real(sp_prec), intent(in):: sp_xv
    !point where the interpolated function is evaluated

    !******* output variables *******
    real(sp_prec), intent(out):: sp_yv
    !value of the interpolated function at sp_xv
    real(sp_prec), intent(out):: sp_yd
    !derivative of the interpolated function at sp_xv
    integer, intent(out):: sp_error
    !error code

    !******* local variables *******
    real(sp_prec), dimension(1:sp_np-1):: sp_zeta
    !array of additional breaking points
    real(sp_prec), dimension(1:sp_np):: sp_d1
    !1st derivatives
    real(sp_prec), dimension(1:sp_np-1):: sp_delta
    !internally used

    !******* construct zeta and d1 *******
    sp_zeta(1:sp_np-1)=sp_work(1:sp_np-1)
    sp_d1(1:sp_np)=sp_work(sp_np:2*sp_np-1)

    !******* check inputs *******
    sp_error=0
    !call sp_check_inputs2(sp_np, sp_x, sp_y, sp_zeta, sp_xv, sp_error)
    !if (sp_error>0) return

    !******* evaluate the value and the derivative *******
    call sp_evald(sp_np, sp_x, sp_y, sp_d1, sp_zeta, sp_xv, sp_yv,&
         sp_yd)

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLINTD

  !********************************************************************
  ! SUBROUTINE SP_CHECK_INPUTS2
  !********************************************************************
  SUBROUTINE SP_CHECK_INPUTS2&
       (sp_np, sp_x, sp_y, sp_zeta, sp_xv, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    real(sp_prec), intent(in), dimension(1:sp_np-1):: sp_zeta
    real(sp_prec), intent(in):: sp_xv

    !******* output variables *******
    integer, intent(out):: sp_error

    !******* local variables *******
    integer:: p0

    !******* set default error code *******
    sp_error=0

    !******* check number of sp_np *******
    if (sp_np<3) then
       sp_error=1
       return
    end if

    !******* check order of grids *******
    do p0=1, sp_np-1
       if (sp_x(p0)>=sp_x(p0+1)) then
          sp_error=2
          return
       end if
    end do

    !******* check position of zeta *******
    do p0=1, sp_np-1
       if ((sp_zeta(p0)<sp_x(p0)) .or. (sp_zeta(p0)>sp_x(p0+1))) then
          sp_error=3
          return
       end if
    end do

    !******* check position of sp_xv *******
    if ((sp_xv<sp_x(1)) .or. (sp_xv>sp_x(sp_np))) then
       sp_error=4
       return
    end if

    !******* end of subroutine *******
    return

  END SUBROUTINE SP_CHECK_INPUTS2

  !********************************************************************
  ! SUBROUTINE SP_EVAL
  !********************************************************************
  SUBROUTINE SP_EVAL(sp_np,sp_x,sp_y,sp_d1,sp_zeta,sp_xv,sp_yv)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y, sp_d1
    real(sp_prec), intent(in), dimension(1:sp_np-1):: sp_zeta
    real(sp_prec), intent(in):: sp_xv

    !******* output variables *******
    real(sp_prec):: sp_yv

    !******* local variables *******
    integer:: p0, p00
    real(sp_prec):: a1, a2, b1, b2, c1, c2, d1bar, a0, b0

    !	!******* find in which interval: p0 *******
    !	do p0=1, sp_np-1
    !		if (sp_x(p0+1)>=sp_xv) exit
    !	end do
    !
    !	!******* find whether left or right of zeta *******
    !	if (sp_xv>=sp_x(p0) .and. sp_zeta(p0)>=sp_xv) then
    !		p00=1
    !	else if (sp_xv>sp_zeta(p0) .and. sp_x(p0+1)>=sp_xv) then
    !		p00=2
    !	else
    !		print *,'Something wrong...'
    !		print *,sp_np,sp_xv
    !		print *,sp_x(:)
    !		print *,sp_y(:)
    !		print *,sp_zeta(:)
    !		stop
    !	end if

    !******* find p0 and p00 *******
    !******* updated from above on 09/15/2003 *******
    do p0=1, sp_np-1
       if (sp_zeta(p0)>=sp_xv) then
          p00=1
          exit
       else if (sp_x(p0+1)>=sp_xv) then
          p00=2
          exit
       else if (p0==sp_np-1) then
          if (sp_zeta(p0)>=sp_xv) then
             p00=1
             exit
          else
             p00=2
             exit
          end if
       end if
    end do

    !******* evaluate the function *******
    a1=sp_y(p0)
    b1=sp_d1(p0)
    a0=sp_zeta(p0)-sp_x(p0)
    b0=sp_x(p0+1)-sp_zeta(p0)
    d1bar=(2._sp_prec*(sp_y(p0+1)-sp_y(p0))-&
         (a0*sp_d1(p0)+b0*sp_d1(p0+1)))/(sp_x(p0+1)-sp_x(p0))
    c1=(d1bar-sp_d1(p0))/(2._sp_prec*a0)
    if (p00==1) then
       sp_yv=a1+b1*(sp_xv-sp_x(p0))+c1*(sp_xv-sp_x(p0))**2._sp_prec
    else
       a2=a1+a0*b1+a0**2._sp_prec*c1
       b2=d1bar
       c2=(sp_d1(p0+1)-d1bar)/(2._sp_prec*b0)
       sp_yv=a2+b2*(sp_xv-sp_zeta(p0))+c2*(sp_xv-sp_zeta(p0))**2._sp_prec
    end if

    !******* end of subroutine *******
    return

  END SUBROUTINE SP_EVAL

  !********************************************************************
  ! SUBROUTINE SP_EVALD
  !********************************************************************
  SUBROUTINE SP_EVALD(sp_np,sp_x,sp_y,sp_d1,sp_zeta,sp_xv,sp_yv,sp_yd)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8
    !double precision

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y, sp_d1
    real(sp_prec), intent(in), dimension(1:sp_np-1):: sp_zeta
    real(sp_prec), intent(in):: sp_xv

    !******* output variables *******
    real(sp_prec):: sp_yv, sp_yd

    !******* local variables *******
    integer:: p0, p00
    real(sp_prec):: a1, a2, b1, b2, c1, c2, d1bar, a0, b0

    !	!******* find in which interval: p0 *******
    !	do p0=1, sp_np-1
    !		if (sp_x(p0+1)>=sp_xv) exit
    !	end do
    !
    !	!******* find whether left or right of zeta *******
    !	if (sp_xv>=sp_x(p0) .and. sp_zeta(p0)>=sp_xv) then
    !		p00=1
    !	else if (sp_xv>sp_zeta(p0) .and. sp_x(p0+1)>=sp_xv) then
    !		p00=2
    !	else
    !		print *,'Something wrong...'
    !		print *,sp_np,sp_xv
    !		print *,sp_x(:)
    !		print *,sp_y(:)
    !		print *,sp_zeta(:)
    !		stop
    !	end if

    !******* find p0 and p00 *******
    !******* updated from above on 09/15/2003 *******
    do p0=1, sp_np-1
       if (sp_zeta(p0)>=sp_xv) then
          p00=1
          exit
       else if (sp_x(p0+1)>=sp_xv) then
          p00=2
          exit
       else if (p0==sp_np-1) then
          if (sp_zeta(p0)>=sp_xv) then
             p00=1
             exit
          else
             p00=2
             exit
          end if
       end if
    end do

    !******* evaluate the function *******
    a1=sp_y(p0)
    b1=sp_d1(p0)
    a0=sp_zeta(p0)-sp_x(p0)
    b0=sp_x(p0+1)-sp_zeta(p0)
    d1bar=(2._sp_prec*(sp_y(p0+1)-sp_y(p0))-&
         (a0*sp_d1(p0)+b0*sp_d1(p0+1)))/(sp_x(p0+1)-sp_x(p0))
    c1=(d1bar-sp_d1(p0))/(2._sp_prec*a0)
    if (p00==1) then
       sp_yv=a1+b1*(sp_xv-sp_x(p0))+c1*(sp_xv-sp_x(p0))**2._sp_prec
       sp_yd=b1+2._sp_prec*c1*(sp_xv-sp_x(p0))
    else
       a2=a1+a0*b1+a0**2._sp_prec*c1
       b2=d1bar
       c2=(sp_d1(p0+1)-d1bar)/(2._sp_prec*b0)
       sp_yv=a2+b2*(sp_xv-sp_zeta(p0))+c2*(sp_xv-sp_zeta(p0))**2._sp_prec
       sp_yd=b2+2._sp_prec*c2*(sp_xv-sp_zeta(p0))
    end if

    !******* end of subroutine *******
    return

  END SUBROUTINE SP_EVALD

  !********************************************************************
  ! SUBROUTINE DSPSPLEASY
  !********************************************************************
  SUBROUTINE DSPSPLEASY (sp_np, sp_x, sp_y, sp_xv, sp_yv, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    real(sp_prec), intent(in):: sp_xv

    !******* output variables *******
    real(sp_prec), intent(out):: sp_yv
    integer, intent(out):: sp_error

    !******* local variables *******
    real(sp_prec), dimension(2*sp_np-1):: sp_work

    !******* construct interpolated function *******
    call dspspline(sp_np, sp_x, sp_y, sp_work, sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+100
       return
    end if

    !******* compute the value of interpolated function at xv *******
    call dspsplint(sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+200
       return
    end if

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLEASY

  !********************************************************************
  ! SUBROUTINE DSPSPLEASYD
  !********************************************************************
  SUBROUTINE DSPSPLEASYD(sp_np,sp_x,sp_y,sp_xv,sp_yv,sp_yd,sp_error)
    implicit none

    !******* precision controller *******
    integer,parameter:: sp_prec=8

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec),intent(in),dimension(1:sp_np):: sp_x,sp_y
    real(sp_prec),intent(in):: sp_xv

    !******* output variables *******
    real(sp_prec),intent(out):: sp_yv,sp_yd
    integer,intent(out):: sp_error

    !******* local variables *******
    real(sp_prec),dimension(2*sp_np-1):: sp_work

    !******* construct interpolated function *******
    call dspspline(sp_np,sp_x,sp_y,sp_work,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+100
       return
    end if

    !******* compute the value of interpolated function at xv *******
    call dspsplintd&
         (sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_yd,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+200
       return
    end if

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLEASYD

  !********************************************************************
  ! SUBROUTINE DSPSPLEASY2
  !********************************************************************
  SUBROUTINE DSPSPLEASY2 (sp_np, sp_x, sp_y, sp_xv, flag_1, sp_d1_1, &
       sp_d1_n, flag_n, sp_yv, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8

    !******* input variables *******
    integer, intent(in):: sp_np, flag_1, flag_n
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y
    real(sp_prec), intent(in):: sp_xv, sp_d1_1, sp_d1_n

    !******* output variables *******
    real(sp_prec), intent(out):: sp_yv
    integer, intent(out):: sp_error

    !******* local variables *******
    real(sp_prec), dimension(2*sp_np-1):: sp_work

    !******* construct interpolated function *******
    call dspspline2(sp_np, sp_x, sp_y, flag_1, sp_d1_1, flag_n, &
         sp_d1_n, sp_work, sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+100
       return
    end if

    !******* compute the value of interpolated function at xv *******
    call dspsplint(sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+200
       return
    end if

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLEASY2

  !********************************************************************
  ! SUBROUTINE DSPSPLEASY3
  !********************************************************************
  SUBROUTINE DSPSPLEASY3 (sp_np, sp_x, sp_y, sp_xv, sp_d, &
       sp_yv, sp_error)
    implicit none

    !******* precision controller *******
    integer, parameter:: sp_prec=8

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec), intent(in), dimension(1:sp_np):: sp_x, sp_y, sp_d
    real(sp_prec), intent(in):: sp_xv

    !******* output variables *******
    real(sp_prec), intent(out):: sp_yv
    integer, intent(out):: sp_error

    !******* local variables *******
    real(sp_prec), dimension(2*sp_np-1):: sp_work

    !******* construct interpolated function *******
    call dspspline3(sp_np, sp_x, sp_y, sp_d, sp_work, sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+100
       return
    end if

    !******* compute the value of interpolated function at xv *******
    call dspsplint(sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+200
       return
    end if

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLEASY3

  !********************************************************************
  ! SUBROUTINE DSPSPLEASYD2
  !********************************************************************
  SUBROUTINE DSPSPLEASYD2(sp_np,sp_x,sp_y,sp_xv,flag_1,sp_d1_1,&
       flag_n,sp_d1_n,sp_yv,sp_yd,sp_error)
    implicit none

    !******* precision controller *******
    integer,parameter:: sp_prec=8

    !******* input variables *******
    integer, intent(in):: sp_np,flag_1,flag_n
    real(sp_prec),intent(in),dimension(1:sp_np):: sp_x,sp_y
    real(sp_prec),intent(in):: sp_xv,sp_d1_1,sp_d1_n

    !******* output variables *******
    real(sp_prec),intent(out):: sp_yv,sp_yd
    integer,intent(out):: sp_error

    !******* local variables *******
    real(sp_prec),dimension(2*sp_np-1):: sp_work

    !******* construct interpolated function *******
    call dspspline2(sp_np,sp_x,sp_y,flag_1,sp_d1_1,flag_n,sp_d1_n,&
         sp_work,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+100
       return
    end if

    !******* compute the value of interpolated function at xv *******
    call dspsplintd&
         (sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_yd,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+200
       return
    end if

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLEASYD2

  !********************************************************************
  ! SUBROUTINE DSPSPLEASYD3
  !********************************************************************
  SUBROUTINE DSPSPLEASYD3(sp_np,sp_x,sp_y,sp_xv,sp_d,&
       sp_yv,sp_yd,sp_error)
    implicit none

    !******* precision controller *******
    integer,parameter:: sp_prec=8

    !******* input variables *******
    integer, intent(in):: sp_np
    real(sp_prec),intent(in),dimension(1:sp_np):: sp_x,sp_y,sp_d
    real(sp_prec),intent(in):: sp_xv

    !******* output variables *******
    real(sp_prec),intent(out):: sp_yv,sp_yd
    integer,intent(out):: sp_error

    !******* local variables *******
    real(sp_prec),dimension(2*sp_np-1):: sp_work

    !******* construct interpolated function *******
    call dspspline3(sp_np,sp_x,sp_y,sp_d,sp_work,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+100
       return
    end if

    !******* compute the value of interpolated function at xv *******
    call dspsplintd&
         (sp_np,sp_x,sp_y,sp_work,sp_xv,sp_yv,sp_yd,sp_error)

    !******* if error occurs, return with the error code+100 *******
    if (sp_error>0) then
       sp_error=sp_error+200
       return
    end if

    !******* end  of the subroutine *******
    return

  END SUBROUTINE DSPSPLEASYD3

END MODULE DSPSP