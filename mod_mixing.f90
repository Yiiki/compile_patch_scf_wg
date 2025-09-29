
module mod_mixing
   implicit none
   integer npulay_max
   parameter (npulay_max=90)
   real*8 AA(npulay_max,npulay_max)        
   integer nreset
   real*8,allocatable, dimension (:) :: R0,w_in0
   real*8,allocatable, dimension (:,:) :: dw,dR

end module mod_mixing
