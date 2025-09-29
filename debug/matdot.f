      program main
      implicit none
      real*8 :: A(3,3),B(3,3),C(3,3),D(3,3)
      open(11,file="matrix")
      read(11,*) A(:,1)
      read(11,*) A(:,2)
      read(11,*) A(:,3)
      read(11,*) B(:,1)
      read(11,*) B(:,2)
      read(11,*) B(:,3)
      close(11)
      C=matmul(A,B) 
      write(6,*) C(:,1)
      write(6,*) C(:,2)
      write(6,*) C(:,3)

      call get_ALI(C,D)
      D=TRANSPOSE(D)
      write(6,*) D(:,1)
      write(6,*) D(:,2)
      write(6,*) D(:,3)
      stop
      end
