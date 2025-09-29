      subroutine getpot2(rho,vtot,n1,n2,n3,AL)
!*****************************************
!c     Written by Lin-Wang Wang, March 30, 2001.  
!************************************************************************
!*  copyright (c) 2003, The Regents of the University of California,
!*  through Lawrence Berkeley National Laboratory (subject to receipt of any
!*  required approvals from the U.S. Dept. of Energy).  All rights reserved.
!************************************************************************

!*****************************************



      implicit double precision (a-h,o-z)

      real*8 rho(n1,n2,n3),vtot(n1,n2,n3)
      real*8 AL(3,3),ALI(3,3),ALI2p(3,3)
      real*8 vr(n1,n2,n3),vi(n1,n2,n3)
      real*8, allocatable, dimension(:) :: wrk


!*************************************************
!cccc generate the potential vtot from rho(i)
      pi=4*datan(1.d0)
      call get_ALI(AL,ALI)

      ALI2p=ALI*2*pi

      lwrk=6*(n1+n2+n3)+15
      allocate(wrk(lwrk))

      vr=rho
      vi=0.d0
   
      call cfft(n1,n2,n3,vr,vi,wrk,lwrk,1)

            do k=1,n3
               ak=0.d0
               if(k.le.(n3+1)/2) ak=k-1
               if(k.ge.n3/2+1) ak=k-n3-1
            do j=1,n2
               aj=0.d0
               if(j.le.(n2+1)/2) aj=j-1
               if(j.ge.n2/2+1) aj=j-n2-1
            do i=1,n1
               ai=0.d0
               if(i.le.(n1+1)/2) ai=i-1
               if(i.ge.n1/2+1) ai=i-n1-1

        akx=ALI2p(1,1)*ai+ALI2p(1,2)*aj+ALI2p(1,3)*ak
        aky=ALI2p(2,1)*ai+ALI2p(2,2)*aj+ALI2p(2,3)*ak
        akz=ALI2p(3,1)*ai+ALI2p(3,2)*aj+ALI2p(3,3)*ak

        akk=akx**2+aky**2+akz**2

        if(akk.gt.1.D-6) then
        vr(i,j,k)=vr(i,j,k)*4*pi/akk
        vi(i,j,k)=vi(i,j,k)*4*pi/akk
        else
        vr(i,j,k)=0.d0
        vi(i,j,k)=0.d0
        endif
            enddo
            enddo
            enddo
      
      call cfft(n1,n2,n3,vr,vi,wrk,lwrk,-1)

!  Assume there is no core charge !
      do k=1,n3
      do j=1,n2
      do i=1,n1
      vtot(i,j,k)=vr(i,j,k)+UxcCA(rho(i,j,k),uxc2)
      enddo
      enddo
      enddo

      return
      end
      


