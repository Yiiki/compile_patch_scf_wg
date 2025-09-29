      subroutine mch_pulay(w_in,w_out,iter,nr)

      use mod_mixing

      implicit double precision (a-h,o-z)

      real*8 w_in(nr),w_out(nr)

      real*8, allocatable, dimension(:,:)  :: AA1
      real*8, allocatable, dimension(:)  :: B

      allocate(AA1(npulay_max,npulay_max))
      allocate(B(npulay_max))

      alpha2=1.d0
cccc alpha2 controls how many recent charge densities
cccc to be used. If alpha2 > 1, then, the very old
cccc charge density is not used effectivly. 
cccc We find that alpha2=1 is O.K.
******************************************************
      if(iter.eq.1) nreset=0
      nint=iter-nreset

      if(nint.gt.npulay_max) then
      write(6,*) "restart pulay, nint0,npulay_max",
     &    iter,npulay_max
      nreset=iter-1
      nint=1
      endif

      if(iter.eq.1) then
      allocate(R0(nr))
      allocate(w_in0(nr))
      allocate(dw(nr,npulay_max))
      allocate(dR(nr,npulay_max))
      endif

      if(nint.eq.1) then
      R0=w_out-w_in
      w_in0=w_in
      endif
******************************************************
      if(nint.gt.1) then

      do i=1,nr
      dw(i,nint-1)=w_in(i)-w_in0(i)
      dR(i,nint-1)=w_out(i)-w_in(i)-R0(i)
      R0(i)=w_out(i)-w_in(i)
      w_in0(i)=w_in(i)
      enddo


      do m=1,nint-1
      s=0.d0
      do i=1,nr
      s=s+dR(i,m)*R0(i)
      enddo
      B(m)=-s
      enddo

******************************************************
      do m=1,nint-1
      s1=0.d0
      do i=1,nr
      s1=s1+dR(i,m)*dR(i,nint-1)
      enddo

      AA(m,nint-1)=s1
      AA(nint-1,m)=s1
      enddo

cccccccccc pulay optimization
**********************************************************
      do m1=1,nint-1
      do m2=1,nint-1
      AA1(m1,m2)=AA(m1,m2)
      enddo
      enddo

      w=1.0d0
      do m=nint-1,1,-1
      AA1(m,m)=AA1(m,m)*w
      w=w*alpha2
      enddo

**********************************************************
      call gaussj(AA1,nint-1,npulay_max,B,1,1)

      w_out=R0

      
      do m=1,nint-1
      do i=1,nr
      w_in(i)=w_in(i)+B(m)*dw(i,m)
      w_out(i)=w_out(i)+B(m)*dR(i,m)
      enddo
      enddo

      w_out=w_out+w_in

      endif

      deallocate(AA1)
      deallocate(B)

      return
      end
      

