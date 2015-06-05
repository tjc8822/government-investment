              module supp

    
    use parameters
    implicit none
    save
    contains
    
    subroutine getEfromIKZ()
        real(8)  ::  tempR1,tempR2,tempR3,tempR4, k, y, z, c, efrt
        integer  ::  ik, iz, iEFRT
        
      
        do ik = 1, nk
            do iz = 1, nz
            K=Kgrid(iK)
            Z=Zgrid(iZ)
            
            do iEFRT=1,nEFRT
                EFRT=0.0001+0.9998*real(iEFRT-1)/real(nEFRT-1)
                Y=Z*(K**alphaK)*(EFRT**(1.0-alphaN))
                if((ik.eq.1).and.(iz.eq.1)) EIgrid1(iEFRT)=EFRT
                if(phi0.gt.0) then
                    C=((1.0-phi0)/phi0)*(1.0-EFRT)*Z*(1.0-alphaN)*(K**alphaK)*(EFRT**(-alphaN))
                else
                    C=.1*Y
                end if
                EIgrid2(ik,iz,iEFRT)=Y-C
            end do
            
            enddo
        enddo
        
    
    
    end subroutine getEfromIKZ
    
    SUBROUTINE findspot(Knext,iKnext)
        implicit none
        integer, intent(out)    ::    iKnext
        integer                 ::    i
        integer                 ::    ik  
        real(8), intent(in)   ::    Knext

      if(nK.eq.1) then
       iKnext=1
       goto 103
      end if

      if(Knext.lt.Kgrid(iK)) then
       do i=0,iK-1
        if(Knext.lt.Kgrid(iK-i)) then
         iKnext=iK-i-1
        else
         goto 101
        end if
       end do
 101   continue
       if(iKnext.lt.1) iKnext=1
      else
       if(iK.eq.nK) then
        iKnext=nK-1
        goto 102
       end if
       do i=iK,nK-1
        if(Knext.ge.Kgrid(i)) then
         iKnext=i
        else
         goto 102
        end if
       end do
 102   continue
      end if

 103  continue

    end subroutine findspot
    
    SUBROUTINE writedata(arrout,filename,nlength,nwidth,mlength, mwidth)
       implicit none
       integer nlength,nwidth,iunit,i,j,mlength,mwidth
       character(len=10) filename
       real(8) arrout(mlength,mwidth)
       data iunit /11/

       open(unit=iunit,file=filename,status='old',form='FORMATTED')

       do i=1,nlength
        write(iunit,*) (arrout(i,j), j=1,nwidth)
       end do

       close(iunit)

       return
      end subroutine writedata

      subroutine linspace(x,N, x_start,x_end,x_len)
            integer, intent(in)                ::   N
            real*8, dimension(N), intent(out)  ::   x
            real*8, intent(in)                 ::   x_start, x_end
            real*8                             ::   dx
            integer, intent(in)                ::   x_len
            integer                            ::   i

            dx  =   (x_end - x_start)/(x_len - 1)
            x(1:x_len) = [(x_start + ((i-1)*dx), i=1,x_len )]


     end subroutine linspace
    
    
    
    
end module supp