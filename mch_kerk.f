      subroutine mch_kerk(w_in,w_out,n1,n2,n3,ALI)

      implicit double precision (a-h,o-z)

      real*8 w_in(n1,n2,n3),w_out(n1,n2,n3)

      real*8 ALI(3,3),ALI2p(3,3)

      real*8, allocatable,dimension (:,:,:)  :: vi,vr

      real*8, allocatable, dimension(:) :: wrk

      pi=4*datan(1.d0)
      ALI2p=ALI*2*pi

       open(77,file="alpha.input",action="read")
       rewind(77)
       read(77,*) alpha1
       read(77,*) alpha2
       read(77,*) alpha3
       close(77)



      lwrk=6*(n1+n2+n3)+15
      allocate(wrk(lwrk))
      allocate(vr(n1,n2,n3))
      allocate(vi(n1,n2,n3))

      w_out=w_out-w_in
      vr=w_out
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

         if(akk.lt.1.D-10) then
         vr(i,j,k)=(alpha3-alpha1)*vr(i,j,k)
         vi(i,j,k)=(alpha3-alpha1)*vi(i,j,k)
         else
         vr(i,j,k)=((alpha1*akk+alpha2)/(akk+0.5)-alpha1)*vr(i,j,k)
         vi(i,j,k)=((alpha1*akk+alpha2)/(akk+0.5)-alpha1)*vi(i,j,k)
         endif
             enddo
             enddo
             enddo

      call cfft(n1,n2,n3,vr,vi,wrk,lwrk,-1)

      w_in=w_in+vr+alpha1*w_out

      deallocate(wrk)
      deallocate(vr)
      deallocate(vi)

**************************************************
      return
      end
      

