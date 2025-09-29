      program patch_new_mpi2
ccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc ind(4(i),iatom): the index of ith neighbore of iatom, used in iiatom(ind), ind=[1,natom]
cccc ind_motif(4(i),iatom): the motif leg index of ith neighbore of iatom, in dxyz_motif(3,ind)
cccc                        ind_motif=[1,4]
cccc dxyz_motif(3,4(j)): the x,y,z direction of the jth motif leg
cccc iat_motif(4(j)): the iiatom of the jth motif leg  
cccc                  iiatom(ind(i,iatom))=iat_motif(ind_motif(i,iatom))
cccc dxyz_indZ(3,4(j),iatom): the x,y,z direction of the jth neighbore of atom iatom

      implicit double precision (a-h,o-z)

      include 'mpif.h'

c      parameter (m=80,rad_box=8.d0)
      parameter (nmap=500000)

      integer status(MPI_STATUS_SIZE)

      real*8, allocatable, dimension(:,:,:) ::vr,rho0,vr_ext,
     &    rho_in,rho_out,vr0

      real*8, allocatable, dimension(:) :: vr_tmp

      real*8 x(3,5000),AL(3,3),ALt(3,3),v1(3),v2(3),v3(3)
      real*8 AL_au(3,3)
      real*8 xZ(3,5000),ALZ(3,3)
      integer ipatch_flag(5000)
      real*8 xv1(3),xv2(3),xv3(3),axis(3,3)
      real*8 xaxis(3,3),xaxis_tmp1(3,3),xaxis_tmp2(3,3)
      real*8 dxyz_ind(3,50,5000),w(5000)
      real*8 dxyz_indZ(3,50,5000)
      real*8 dxyz_motif(3,50,500)
      real*8 op(3,3),Efield(3,5000)
      integer ind(50,5000),iiatom(5000),nH_neigh(5000)
      integer iatype(5000),iatype_el(10)
      integer num_neigh_ALL(5000)
      integer ind_motif0(5000),err_motif(5000)
      real*8 op_m(3,3,5000)
      integer ind_config(5000)

      real*8, allocatable, dimension(:,:,:) :: dens,dens_tmp

      real*4, allocatable, dimension(:,:,:) ::  motif_tmp
      real*4, allocatable, dimension(:,:,:,:,:) :: dens_mEM

      real*8 x_neighM(3,50,500)
      integer iat_neighM(50,500),num_neighM(500)
      integer iat_neigh(50,5000)

      character*200 file_Emotif(500)

      real*8 ALI(3,3),amatrix(3,3)
      real*8 AL_mbox(3,3),AL_tmp(3,3)
      integer ijk_map(3,nmap)
      real*4 dens_map(nmap)
      character*30 filename,filenameZ,file_tmp
      real*8 size_Efield(10),EfieldM(3,10),Efieldtmp(3)
      real*8 fact_motif(1000)
      integer ind_motif_fact(1000)
      real*8,allocatable,dimension(:,:,:) :: dipmot
      real*8 :: dipall(3,3),pi,mot_x,mot_y,mot_z,vol_fac
      integer*4 :: ichi,jchi,ixyz,jxyz,kxyz
      real*8 :: volc,vole,mean,mot_u,mot_v,mot_w,oq(3,3),mat_tmp(3,3)

      real*8,allocatable,dimension(:,:,:) :: matset,d1,d2,d3
      integer*4 :: num_mat_fac
      logical :: res
!!!!! initialized as identical matrix
      dipall=0.d0
      dipall(1,1)=1.d0      
      dipall(2,2)=1.d0      
      dipall(3,3)=1.d0      

      open(11,file="patchE.input")
      rewind(11)
      read(11,*) filename,filenameZ
      read(11,*) n1,n2,n3
      read(11,*) ipatch_selective
      read(11,*) num_motif_fact
      do ii=1,num_motif_fact
      read(11,*) ind_motif_fact(ii),fact_motif(ii)
      enddo
      close(11)

      write(6,*) "patchE.input ... read finished."

      open(11,file="MOTIF_list_all")
      rewind(11)
!      read(11,*)
      read(11,*) num_Emotif
      allocate(matset(1:3,1:3,num_Emotif))
      matset=0.d0

      do ii=1,num_Emotif

      matset(1,1,ii)=1.d0
      matset(2,2,ii)=1.d0
      matset(3,3,ii)=1.d0


      read(11,*) ii0,jj0,file_Emotif(ii)
      read(11,*)
      enddo
      close(11)


      write(6,*) "MOTIF_list_all ... read finished."


!     new: we need a matrix factor
      inquire(file="MOTIF_matrix_factor",exist=res)
      if(res) then
              write(6,*) "MOTIF_matrix_factor ... enabled"
              open(11,file="MOTIF_matrix_factor")
              rewind(11)
              read(11,*) num_mat_fac
                  if(num_mat_fac.ne.num_Emotif) then
                       write(6,*) "wrong num_mat_fac", num_mat_fac
                       stop
                  end if
                  do ii=1,num_mat_fac
                     read(11,*)
                     read(11,*)
                     read(11,*) matset(1:3,1,ii) 
                     read(11,*) matset(1:3,2,ii) 
                     read(11,*) matset(1:3,3,ii)
                     matset(1,1,ii)=matset(1,1,ii)+1.d0 
                     matset(2,2,ii)=matset(2,2,ii)+1.d0 
                     matset(3,3,ii)=matset(3,3,ii)+1.d0 
                  end do 
              close(11)
      end if


      open(11,file="motif.input")
      rewind(11)
      read(11,*) m,rad_box
      read(11,*) AL_mbox(1,1),AL_mbox(2,1),AL_mbox(3,1)! [A]
      read(11,*) AL_mbox(1,2),AL_mbox(2,2),AL_mbox(3,2)
      read(11,*) AL_mbox(1,3),AL_mbox(2,3),AL_mbox(3,3)
      read(11,*) dcut
      read(11,*) ntype
      do ii=1,ntype
      read(11,*) iatype_el(ii),file_tmp,size_Efield(ii)       ! Ga  2
      enddo
      close(11)

      write(6,*) "motif.input ... read finished."


      rad_box2=rad_box**2

      allocate(dens_mEM(-m:m,-m:m,-m:m,3,num_Emotif))
      allocate(d1(-m:m,-m:m,-m:m))
      allocate(d2(-m:m,-m:m,-m:m))
      allocate(d3(-m:m,-m:m,-m:m))


      do ii=1,num_Emotif
      open(12,file="MISSING_MOTIFS/E"//file_Emotif(ii),
     &        form="unformatted")
      rewind(12)
      read(12) num_neighM(ii),rad_tmp
      read(12) AL_tmp                                 ! [A]

      if(abs(rad_tmp-rad_box).gt.1E-8) then
      write(6,*) "rad_box changed,stop",rad_box,rad_tmp
      stop
      endif
      if(sum(abs(AL_mbox-AL_tmp)).gt.1.E-6) then
      write(6,*) "AL_mbox changed, stop"
      write(6,*) AL_mbox
      write(6,*) AL_tmp
      stop
      endif

      do ia=1,num_neighM(ii)
       read(12) iat_neighM(ia,ii),x_neighM(1,ia,ii),x_neighM(2,ia,ii),
     &      x_neighM(3,ia,ii)
      enddo
      read(12) mb,EfieldM(1,ii),EfieldM(2,ii),EfieldM(3,ii)
      if(mb.ne.m) then
      write(6,*) "mb changed, stop",mb,m
      stop
      endif

      read(12) dens_mEM(:,:,:,1,ii)
      read(12) dens_mEM(:,:,:,2,ii)
      read(12) dens_mEM(:,:,:,3,ii)

      d1=dens_mEM(:,:,:,1,ii)
      d2=dens_mEM(:,:,:,2,ii)
      d3=dens_mEM(:,:,:,3,ii)


! new introduce M matrix
      dens_mEM(:,:,:,1,ii)=matset(1,1,ii)*d1+matset(1,2,ii)*d2+
     &                     matset(1,3,ii)*d3
      dens_mEM(:,:,:,2,ii)=matset(2,1,ii)*d1+matset(2,2,ii)*d2+
     &                     matset(2,3,ii)*d3
      dens_mEM(:,:,:,3,ii)=matset(3,1,ii)*d1+matset(3,2,ii)*d2+
     &                     matset(3,3,ii)*d3
!

      write(6,"(A,3x,I,3x,A)") "EMOTIF",ii,"read finished"
      close(12)
      enddo

      do ii=1,num_motif_fact
      jj=ind_motif_fact(ii)
      dens_mEM(:,:,:,:,jj)=dens_mEM(:,:,:,:,jj)*fact_motif(ii)
      enddo

!cccccccccccccccccccccccccccccccccccccccccccccccc
      open(10,file=filename)
      rewind(10)
      read(10,*) natom
      read(10,*) 
      read(10,*) AL(1,1),AL(2,1),AL(3,1)         ! [A]
      read(10,*) AL(1,2),AL(2,2),AL(3,2)
      read(10,*) AL(1,3),AL(2,3),AL(3,3)
      read(10,*) 
      do i=1,natom
      if(ipatch_selective.eq.0) then
      read(10,*) iiatom(i),x(1,i),x(2,i),x(3,i)
      ipatch_flag(i)=1
      else
      read(10,*) iiatom(i),x(1,i),x(2,i),x(3,i),im1,im2,im3,
     &  ipatch_flag(i)
      endif
      x(1,i)=mod(x(1,i)+2.d0,1.d0)
      x(2,i)=mod(x(2,i)+2.d0,1.d0)
      x(3,i)=mod(x(3,i)+2.d0,1.d0)
      enddo
      close(10)


!!!!! summ up the motif with the weight they are used in the list
      open(12,file="OUT.atom_motif.list")
      rewind(12)
      read(12,*) natom_patch
      do ii=1,natom_patch
      read(12,*) ind_config(ii),ind_motif0(ii),err_motif(ii)
      read(12,*)
      read(12,*) op_m(1,1,ii),op_m(2,1,ii),op_m(3,1,ii)
      read(12,*) op_m(1,2,ii),op_m(2,2,ii),op_m(3,2,ii)
      read(12,*) op_m(1,3,ii),op_m(2,3,ii),op_m(3,3,ii)
      enddo
      close(12)

      write(6,*) "natom_patch=",natom_patch
      if(inode.eq.1) then
              if(natom_patch.le.natom) then
                      write(6,*) "natom_patch < natom, eps might bad"
              endif
      endif

      allocate(dipmot(3,3,natom_patch))
      dipmot=0.d0
      num=0
      do ii=1,natom
      if(ipatch_flag(ii).eq.1) then
      num=num+1
        if(ind_config(num).ne.ii) then
        write(6,*) "OUT.atom_motif.list,atom.config not match"
        stop
        endif
      endif
      enddo

      if(num.ne.natom_patch) then
       write(6,*) "natom_patch not right",num,natom_patch
       stop
      endif

!!!!  config box volume
      volc=AL(1,1)*(AL(2,2)*AL(3,3)-AL(3,2)*AL(2,3))+
     &        AL(2,1)*(AL(3,2)*AL(1,3)-AL(1,2)*AL(3,3))+
     &        AL(3,1)*(AL(1,2)*AL(2,3)-AL(2,2)*AL(1,3))
      volc=dabs(volc)
!!!!  emotif box volume
      AL=AL_mbox
      vole=AL(1,1)*(AL(2,2)*AL(3,3)-AL(3,2)*AL(2,3))+
     &        AL(2,1)*(AL(3,2)*AL(1,3)-AL(1,2)*AL(3,3))+
     &        AL(3,1)*(AL(1,2)*AL(2,3)-AL(2,2)*AL(1,3))
      vole=dabs(vole)

      pi=4.d0*datan(1.d0)
      vol_fac=-4.d0*pi/((2*mb)**3) ! this "-" sign comes from Gauss law
     &*(vole/volc)

      do 100 iatom_patch=1,natom_patch

      iatom=ind_config(iatom_patch)
      itype=ind_motif0(iatom_patch)
!*    one needs this to rotate the EMOTIF to its direct counterpart      
!*    originally, one should have one-to-one mapping of EMOTIF to atoms
!*    however, only limited EMOTIF will be stored, thus op is necessary      
      op(:,:)=op_m(:,:,iatom_patch)
      call get_ALI(op,oq)
      oq=TRANSPOSE(oq)

      do jchi=1,3

        do ixyz=-mb,mb-1
          mot_x=ixyz*AL_mbox(1,1)/(2*mb)/0.529177d0

        do jxyz=-mb,mb-1
          mot_y=jxyz*AL_mbox(2,2)/(2*mb)/0.529177d0

        do kxyz=-mb,mb-1
          mot_z=kxyz*AL_mbox(3,3)/(2*mb)/0.529177d0

          dipmot(1,jchi,iatom)=dipmot(1,jchi,iatom)+
     &    mot_x*dens_mEM(ixyz,jxyz,kxyz,jchi,itype)*vol_fac

          dipmot(2,jchi,iatom)=dipmot(2,jchi,iatom)+
     &    mot_y*dens_mEM(ixyz,jxyz,kxyz,jchi,itype)*vol_fac

          dipmot(3,jchi,iatom)=dipmot(3,jchi,iatom)+
     &    mot_z*dens_mEM(ixyz,jxyz,kxyz,jchi,itype)*vol_fac
        enddo
        enddo 
        enddo

      enddo
!     tensor transform as P^-1 R P
!     where R = r ⊗ ∇⋅J 

      mat_tmp=matmul(oq,dipmot(1:3,1:3,iatom))
      mat_tmp=matmul(mat_tmp,op)
      dipmot(1:3,1:3,iatom)=mat_tmp
      ! sum atom
      dipall(1:3,1:3)=dipall(1:3,1:3)+dipmot(1:3,1:3,iatom)

!                                                                 !
      write(6,"(A23,I4,A21)") "------4πχ--of--atom--",iatom,
     &"--------------------"
      write(6,"(3F16.9)") dipmot(:,1,iatom)
      write(6,"(3F16.9)") dipmot(:,2,iatom)
      write(6,"(3F16.9)") dipmot(:,3,iatom)

cccccccccccccccccccccccccccccccccccccccccccccccc
100   continue         ! do all the atoms for this node 

      mean=(dipall(1,1)+dipall(2,2)+dipall(3,3))/3.d0
!                                                                 !
!     write(6,"(A48)") "-------1-+-sum--of--4πχ--to--atoms--------"
      if(res) then
      write(6,"(A21,A16,F16.9)") "-ε(∞)+4πχ(ion)--",
     &"---trace-avg-->",mean
      else
      write(6,"(A19,A16,F16.9)") "------ε(∞)------",
     &"---trace-avg-->",mean
      end if
      write(6,"(3F16.9)") dipall(:,1)
      write(6,"(3F16.9)") dipall(:,2)
      write(6,"(3F16.9)") dipall(:,3)

      stop
      end
