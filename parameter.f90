module parameters 


    
    implicit none  
    integer, parameter       ::     nK=100,nZ=3,nM=15,nEFRT=30
    real(8)                 ::     delta0, alphaK, alphaN, v0, beta0, phi0,gamma0
    real(8), dimension(3,3) ::     TrProb
    real(8)                 ::     Zgrid(3), Kgrid(nK)
    
    real(8)                 ::     EIgrid1(nEFRT),EIgrid2(nK,nZ,nEFRT)
    real(8)                 ::     policyINV(nK,nZ),policyEFRT(nK,nZ), policyC(nK,nZ)  
                                                      
    
    save
     !Annual
      !1.0 for log utility  (note that FOC for leisure is the same for log and power utility so that getEfromIKZ does not need gamma0)
      
    contains 
    subroutine assign_value()  ! Can use Rowenhurst Later
    integer   :: i
    
      Zgrid(1)=0.97
      Zgrid(2)=1.0
      Zgrid(3)=1.03
      TrProb(1,1)=0.8 !0.75
      TrProb(1,2)=0.2 !0.25
      TrProb(1,3)=0.0 !0.0
      TrProb(2,1)=0.1 !0.25
      TrProb(2,2)=0.8 !0.5
      TrProb(2,3)=0.1 !0.25
      TrProb(3,1)=0.0 !0.0
      TrProb(3,2)=0.2 !0.25
      TrProb(3,3)=0.8 !0.75
      
!Quarterly
      delta0=0.025
      alphaK=0.36 !capital share
      alphaN=0.36 !1 - labor share
      v0=0  !adjustment cost
      beta0=0.99 !time preference
      phi0=0.667 !lisure share in utility, set to 0 to make labor inelastic
      gamma0=1.0 !risk aversion, set to 0 for log utility
!Annual
      delta0=0.10
      alphaK=0.36
      alphaN=0.36 !labor share
      v0=0.0  !adjustment cost
      beta0=0.96 !time preference
      phi0=0.667 !lisure share in utility, set to 0 to make labor inelastic
      gamma0=8.!0.0 !risk aversion, set to  
      
      

    
    end subroutine assign_value
    
    
end module parameters