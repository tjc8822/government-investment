! Jincheng: This program extend Jack Favilikus's Code in the following way:
 
! 1. Howard Policy Iteration 
! 2. MacQueen-Porteus Error Bounds
! 3. Cubic spline to approximate value function

! Conslusion: Due to the avoidance of rott finding and all the extensions, the speed has decreased by 25 times!!!!!!!! 160 seconds to 4-5 seconds. WOW
program main ! The main program to do Value Function Iteration
      
      !include 'link_fnl_shared.h'
      use parameters
      use dspsp
      use supp
      USE RNUN_INT
      !USE UMACH_INT
      !USE RNSET_INT       
      !use CSIEZ_INT
      !use UMACH_INT
      implicit none
      integer,parameter         ::           Iseed = 86456
      real(8)                   ::           multK,C,EFRT,INV,K,Knext,Y,Z,INVcost,multEI,minival,Util,EVnext,Vnext,R(1)
      integer                   ::           iIter, ierr    
      integer                   ::           t, th 
      real(8)                   ::           tempR1
      integer                   ::           iI,j,iZnext, iHoward
      integer, parameter        ::           nI0=100,nI1=100,nI2=100, nHoward = 15 ! nI0, nI1, nI2 are critical for convergence
      real(8), parameter        ::           toler = 0.000000001
      real(8)                   ::           Igrid2(2),Igrid1(2),Igrid(nI0),x(nk),fx(nk),xx,fxx, xxx(1), fxxx(1)
      real(8)                   ::           Vgrid(nK,nZ), Vmax(nK,nZ)
      real(8)                   ::           perError(nK,nZ)

      
      real(8)                   ::           g(nk), temp1, temp
      real(8)                   ::           start, finish
      real(8)                   ::           mpl, mpu
      real(8)                   ::           xtemp(nefrt), fxtemp(nefrt)
      integer                   ::           iK,iZ,iAgg,nAgg,iEFRT, i, ikp


      !! For simulations 
      integer, parameter        ::           years = 2500
      real(8)                   ::           SDF, Q, Jincheng, TrProbCum(nz,nz), wage, Clast, EFRTlast
      real(8)                   ::           output1(years,15) 
                                                               !  


      mpl = 0.
      mpu = 5.
      call CPU_TIME(start)
      call assign_value
      !do i=1,nK
       !Kgrid(i)=.01+1.0*real(i-1)/real(nK-1)    !Grid for problem that can be solved analytically: delta=1, alpha=.36, beta=.96, v=0
       !Kgrid(i)=0.01+(3.0-0.01)*real(i-1)/real(nK-1)  !annual
       !Kgrid(i)=0.01+(20.0-0.01)*real(i-1)/real(nK-1)  !quarterly
      !end do
      
         perError = 5.       
      

      ! Try Unequally Spaced Point Here
        temp   =   (3.0-0.01)**(1.0/gamma0)
        temp1  =   0.
        call linspace(g,nk,temp1,temp,nk)
        do i=1, nk
           kgrid(i) = g(i) ** (gamma0) + 0.01     
        enddo
              
      call getEfromIKZ
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!                                   VFI
      !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      minival=-10.0**30.0
      nAgg=nK*nZ
      do iI=1,nI0
       Igrid(iI)=-2.0+4.0*real(iI-1)/real(nI0-1)
      end do

      do ik = 1, nk
        do iz = 1, nz
           K=Kgrid(iK)
           if(abs(gamma0-1.0).lt..0001) then   !set initial value function
            Vgrid(ik, iz)=(alphaK/(1-alphaK*beta0))*log(K) ! Preserve Concavity and Monotonicity of Value Function
           else
            Vgrid(ik, iz)=(K**(1.0-gamma0))/(1.0-gamma0)
           end if
        enddo
      end do

do while ((mpu-mpl) .ge. toler)
   do ik = 1, nk
      do iz = 1, nz  
       
       
       Vmax(ik,iz)=minival

       K=Kgrid(iK)
       Z=Zgrid(iZ) 
        do iI=1,nI0+nI1+nI2
            if(iI.le.nI0) then
             INV=Igrid(iI)
            else
             if(iI.le.nI0+nI1) then
              INV=Igrid1(1)+(Igrid1(2)-Igrid1(1))*real(iI-nI0)/real(nI1)
             else
              INV=Igrid2(1)+(Igrid2(2)-Igrid2(1))*real(iI-nI0-nI1)/real(nI2)
             end if
            end if
            InvCost=INV+v0*(((INV/K)-delta0)**2.0)*K   !allows for quadratic adjustment cost

           

            ! JT: Use Cubic Spline to Pin Down Effort Level Given Investment Grid
             do j=1, nefrt
                  xtemp(j)      =   EIgrid2(ik,iz,j)
                  fxtemp(j)     =   EIgrid1(j)
             enddo
             xx         =   INVCost
             
             if (xx.lt.xtemp(1) ) then
                xx = xtemp(1)
             elseif (xx.gt.xtemp(nefrt)) then
                xx = xtemp(nefrt)
             endif
                call dspspleasy(nefrt,xtemp,fxtemp,xx,fxx,ierr)
                !call CSIEZ(xtemp, fxtemp, xxx, fxxx)
             if (ierr>0) then
                print*, "Error in Interpolation EFRT GRID: Error Code", ierr
             stop
             end if
             EFRT = fxx
             !EFRT  =  fxxx(1)
             if (EFRT.ge.EIgrid1(nefrt)) then
                 EFRT =   EIgrid1(nefrt)
             endif

             if (EFRT.le.EIgrid1(1)) then
                 EFRT =   EIgrid1(1)
             endif

            Y=Z*(K**alphaK)*(EFRT**(1-alphaN))
            C=Y-InvCost
            if(C.gt.0.0001) then
             if(abs(gamma0-1.0).lt..00001) then
              Util=(1.0-phi0)*log(C)+phi0*log(1.0-EFRT)
              if(abs(phi0).lt.0.00001) Util=log(C)
             else
              Util=(((C**(1.0-phi0))*((1.0-EFRT)**phi0))**(1.0-gamma0))/(1.0-gamma0)
              if(abs(phi0).lt.0.00001) Util=(C**(1.0-gamma0))/(1.0-gamma0)
             end if
            else
             Util=minival
            end if
            Knext=(1.0-delta0)*K+INV
            if(Knext.lt.0.0) then
             Knext=0.0
             Util=minival
            end if
           
            
            EVnext=0.0
            
            do iZnext=1,nZ
                 x(:) = kgrid(:)
                 do  ikp = 1 , nK                
                     fx(ikp) = Vgrid(ikp,iZnext)
                 enddo
                 xx = knext
                 !xxx(1) = knext
                 if (xx.lt.kgrid(1)) then
                    xx = kgrid(1)
                 elseif (xx.gt.kgrid(nk)) then 
                    xx = kgrid(nk)
                 endif

                  
                 call dspspleasy(nK,x,fx,xx,fxx,ierr)
                 !call CSIEZ(x,fx,xxx,fxxx)
                 if (ierr>0) then
                    print*, "Error in Interpolation knext grid: Error Code", IERR
                    stop
                 endif
                 Vnext=fxx
                 !Vnext = fxxx(1)
                 EVnext=EVnext+TrProb(iZ,iZnext)*Vnext
            end do

            
            tempR1=Util+beta0*EVnext
            
            if(tempR1.gt.Vmax(ik,iz)) then
                 Vmax(ik,iz)=tempR1
                 policyINV(ik,iz)=INV
                 policyEFRT(ik,iz)=EFRT
                 policyC(ik,iz) = C
                 if(iI.le.nI0) then
                  Igrid1(1)=INV-2.0*(Igrid(2)-Igrid(1))
                  Igrid1(2)=INV+2.0*(Igrid(2)-Igrid(1))
                  Igrid2(1)=INV-0.25*(Igrid(2)-Igrid(1))
                  Igrid2(2)=INV+0.25*(Igrid(2)-Igrid(1))
                  if(iI.le.2) then
                   Igrid1(1)=INV-2.0
                   Igrid2(1)=INV-2.0
                  end if
                  if(iI.ge.nI0-1) then
                   Igrid1(2)=INV+2.0
                   Igrid2(2)=INV+2.0
                  end if
                end if
                 if(iI.gt.nI0.and.iI.le.nI0+nI1) then
                  Igrid2(1)=INV-.25*(Igrid(2)-Igrid(1))
                  Igrid2(2)=INV+.25*(Igrid(2)-Igrid(1))
                  if(iI.le.nI0+2) then
                   Igrid2(1)=INV-2.0
                  end if
                  if(iI.ge.nI0+nI1-1) then
                   Igrid2(2)=INV+2.0
                 end if
                 end if
            end if
        
            
        enddo   ! iI

            do i = 1, nK
                do j = 1, nz
                perError(i,j) = abs(Vmax(i,j)-Vgrid(i,j)) 
                enddo
            enddo
            
     enddo        
   enddo ! iAgg
   mpl = -1* beta0/(1-beta0) * maxval(-(Vmax-Vgrid))
   mpu = beta0/(1-beta0) * maxval(Vmax-Vgrid)
   print*, 'Upper Bound is ', mpu, 'Lower BOUND IS ', mpl
   print*, 'Max Difference is ', maxval(abs(Vmax-Vgrid))
   
   Vgrid=Vmax
   
   
   t = t + 1
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!   Howard Policy Iteration
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (maxval(perError).le. 2. ) then
       th =th+1       
       do iHoward = 1, nHoward
          do ik = 1, nk
             do iz = 1, nz   
             
             
             K=Kgrid(iK)
             Z=Zgrid(iZ)  
             Efrt = policyEFRT(ik, iz)
             INV=policyINV(ik, iz)
             InvCost=INV+v0*(((INV/K)-delta0)**2.0)*K
             Y=Z*(K**alphaK)*(EFRT**(1-alphaN))
             C=Y-InvCost
             if(C.gt.0.0001) then
                if(abs(gamma0-1.0).lt..00001) then
                    Util=(1.0-phi0)*log(C)+phi0*log(1.0-EFRT)
                    if(abs(phi0).lt.0.00001) Util=log(C)
                else
                    Util=(((C**(1.0-phi0))*((1.0-EFRT)**phi0))**(1.0-gamma0))/(1.0-gamma0)
                    if(abs(phi0).lt.0.00001) Util=(C**(1.0-gamma0))/(1.0-gamma0)
                end if
                else
                    Util=minival
                end if
                
                Knext=(1.0-delta0)*K+INV
                
                
            
                EVnext=0.0
            
                do iZnext=1,nZ
                     x(:) = kgrid(:)
                     do  ikp = 1 , nK                
                     fx(ikp) = Vgrid(ikp, iznext)
                     enddo
                     xx = knext
                     !xxx(1) = knext
                     if (xx.lt.kgrid(1)) then
                        xx = kgrid(1)
                        elseif (xx.gt.kgrid(nk)) then 
                        xx = kgrid(nk)
                     endif 
                 
                     call dspspleasy(nk,x,fx,xx,fxx,ierr)
                     !call CSIEZ(x,fx,xxx,fxxx)
                     Vnext=fxx  
                     EVnext=EVnext+TrProb(iZ,iZnext)*Vnext

                end do
                tempR1=Util+beta0*EVnext
                Vmax(ik, iz) = tempR1         
          enddo   !iAgg
       enddo

         
                Vgrid = Vmax
          
       enddo
       
       
   endif
   
   
   
enddo    ! End do while Loop


open(unit=10,file='polyINV_1.txt')
open(unit=20,file='polyEFRT_1.txt')
open(unit=25,file='polyC_1.txt')
open(unit=30,file='k.txt')
do i=1,nK
        write(10,*) (policyINV(i,j), j=1,nz)
        write(20,*) (policyefrt(i,j), j=1,nz)
        write(25,*) (policyc(i,j), j=1,nz)
        write(30,*) kgrid(i)
end do

close(10)
close(20)
close(30)

call CPU_TIME(finish)
print '("Time = ",f10.3," seconds.")',finish-start


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!                   Simulation
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! Set transition probability

       do iz = 1, nZ
          TrProbCum(iz, 1) = TrProb(iz, 1)  
          do iZnext = 2, nZ
            TrProbCum(iz, iznext) = TrProbCum(iz, iznext-1) + TrProb(iz, iznext)
          end do     
       end do
       write(6,*) TrProbCum(1,1),TrProbCum(1,2),TrProbCum(1,3)
       write(6,*) TrProbCum(2,1),TrProbCum(2,2),TrProbCum(2,3)
       write(6,*) TrProbCum(3,1),TrProbCum(3,2),TrProbCum(3,3)
       !! Initialize Simulation
       
       Knext=0.5*Kgrid(nK)
       C=0.1*Knext
       EFRT=0.5
       iZnext=nint(0.51*real(nZ))

       !! Simulation 
       do t = 1, years
          K = Knext
          Clast = C
          EFRTlast = EFRT
          iZ=iZnext

          xx        =     K
          call dspspleasy(nk,Kgrid, policyINV(:, iZ), xx, fxx,ierr)
          INV       =     fxx
          call dspspleasy(nk,Kgrid, policyEFRT(:, iZ), xx, fxx,ierr)
          EFRT      =     fxx

          Y=Zgrid(iZ)*(K**alphaK)*(EFRT**(1.0-alphaN))
          InvCost=INV+v0*(((INV/K)-delta0)**2.0)*K
          C=Y-InvCost
          Wage=(1.0-alphaN)*Zgrid(iZ)*(K**alphaK)*(EFRT**(-alphaN))
          Knext=K*(1-delta0)+INV

          SDF=beta0*((C/Clast)**((1.0-phi0)*(1.0-gamma0)-1.0))*(((1.0-EFRT)/(1.0-EFRTlast))**(phi0*(1.0-gamma0)))

            if(v0.gt.0) then
                tempR1=(Zgrid(iZ)**(1.0/alphaN))*(Wage**((alphaN-1.0)/alphaN))* alphaN*(1.0-alphaN)**((1.0-alphaN)/alphaN)
                Q=tempR1-(INV/K)-v0*(((INV/K)-delta0)**2.0)+(1.0-delta0+(INV/K))*(1.0+2*v0*((INV/K)-delta0))
           else
                tempR1=alphaK*Zgrid(iZ)*(K/EFRT)**(alphaK-1.0)
                Q=tempR1+1.0-delta0
           end if

           output1(t,1)=real(iZ)
           output1(t,2)=K
           output1(t,3)=Y
           output1(t,4)=INV
           output1(t,5)=C
           output1(t,6)=EFRT
           output1(t,7)=K/Y
           output1(t,8)=alphaK*Zgrid(iZ)*((K/EFRT)**(alphaK-1))-delta0 !Z*(df/dK)-delta=return on capital
           output1(t,9)=Y/EFRT
           output1(t,10)=Zgrid(iZ)
           output1(t,11)=Wage
           output1(t,12)=C-Wage*EFRT !dividend
           output1(t,13)=SDF
           output1(t,14)=Q
           output1(t,15)=Q*K
           
           
           call RANDOM_NUMBER(tempR1)
           iZnext=1
           do i=1,nZ-1
              if(tempR1.gt.TrProbCum(iZ,i)) iZnext=i+1
           end do
       end do
       1000 format (15F10.4)
       open(unit=40,file='output1.txt')
        do i=1,years
                write(40,1000) (output1(i,j), j=1,15)                
        end do

        close(40)











end program main