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

      real*8 :: chi(3,3),epsr(3,3),pi

      integer*4 :: edim

      chi=0.d0
      epsr=0.d0

      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,inode,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,nnodes,ierr)

      inode=inode+1

!     do i=1,nnodes
!     if(inode.eq.i) then
!             write(6,"(A,3x,I3,3x,A)") "procs",i,"ready"
!     endif
!             call mpi_barrier(MPI_COMM_WORLD,ierr)
!     enddo

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

      if(inode.eq.1) then
              write(6,*) "patchE.input ... read finished."
      endif

      open(11,file="MOTIF_list_all")
      rewind(11)
!      read(11,*)
      read(11,*) num_Emotif
      do ii=1,num_Emotif
      read(11,*) ii0,jj0,file_Emotif(ii)
      read(11,*)
      enddo
      close(11)


      if(inode.eq.1) then
              write(6,*) "MOTIF_list_all ... read finished."
      endif


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

      if(inode.eq.1) then
              write(6,*) "motif.input ... read finished."
      endif


      rad_box2=rad_box**2

      allocate(dens_mEM(-m:m,-m:m,-m:m,3,num_Emotif))


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

      if(inode.eq.1) then
      read(12) dens_mEM(:,:,:,1,ii)
      read(12) dens_mEM(:,:,:,2,ii)
      read(12) dens_mEM(:,:,:,3,ii)
      write(6,"(A,3x,I,3x,A)") "EMOTIF",ii,"read finished"
      endif
      close(12)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo


      if(inode.eq.1) then
              write(6,*) "start bcast dens_mEM"
      endif
      
      call mpi_bcast(dens_mEM,(2*m+1)**3*3*num_Emotif,MPI_REAL4,0,
     &   MPI_COMM_WORLD,ierr)    ! send vr to other processors

      if(inode.eq.1) then
              write(6,*) "finished bcast dens_mEM"
      endif
 

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

      w(i)=1.d0
        iatype(i)=0
        do itype=1,ntype
        if(iiatom(i).eq.iatype_el(itype)) iatype(i)=itype
        enddo

        if(iatype(i).eq.0) then
        write(6,*) "something wrong",iiatom(i),iatype(i)
        stop
        endif
cccccccc  iatype(i)=0 for H atoms
c        if(iatype(i).eq.0) then
c        write(6,*) "atom type not found in motif.input"
c        endif
      enddo
      close(10)

      open(10,file=filenameZ)
      rewind(10)
      read(10,*) natomZ
      read(10,*) 
      read(10,*) ALZ(1,1),ALZ(2,1),ALZ(3,1)
      read(10,*) ALZ(1,2),ALZ(2,2),ALZ(3,2)
      read(10,*) ALZ(1,3),ALZ(2,3),ALZ(3,3)
       if(natom.ne.natomZ) then
       write(6,*) "natom.ne.natomZ, stop"
      call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
      read(10,*) 
      do i=1,natomZ
      read(10,*) iat_tmp2,xZ(1,i),xZ(2,i),xZ(3,i)
      if(iat_tmp2.ne.iiatom(i)) then
      write(6,*) "iat_tmp2.ne.iiatom(i), stop", 
     &   i, iat_tmp2, iiatom(i)
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif
      enddo
      close(10)


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





      call get_ALI(AL,ALI)

ccccccc inside neighbore, x, xZ cannot be shifted by 1. !
      call neighbore(dcut)

      allocate(dens(n1,n2,n3))

cccccccccccccccccccccccccccccccccccccccccccccc
!cccc now, calculate the electric field and the average electric field around each atom
!ccccccccccccccccccccccccccccccccccccccccccccccc

      vol=AL(1,1)*(AL(2,2)*AL(3,3)-AL(3,2)*AL(2,3))+
     &    AL(2,1)*(AL(3,2)*AL(1,3)-AL(1,2)*AL(3,3))+
     &    AL(3,1)*(AL(1,2)*AL(2,3)-AL(2,2)*AL(1,3))
      vol=dabs(vol)
      vol_fac=vol/(n1*n2*n3)
      AL_au=AL/0.529177d0


!      call getpot2(rho0,vr0,n1,n2,n3,AL_au)

      rrcut=(rad_box)**2
      rcutx1=rad_box*dsqrt(ALI(1,1)**2+ALI(2,1)**2+ALI(3,1)**2)
      rcutx2=rad_box*dsqrt(ALI(1,2)**2+ALI(2,2)**2+ALI(3,2)**2)
      rcutx3=rad_box*dsqrt(ALI(1,3)**2+ALI(2,3)**2+ALI(3,3)**2)
cccc rcutx1,2,3 are the cutoff in terms of x1,x2,x3 (in the AL coord) 
cccc for the patched sphere of radius: rad_box
      ircut1=n1*rcutx1+1
      ircut2=n2*rcutx2+1
      ircut3=n3*rcutx3+1


      do edim=1,3

      allocate(motif_tmp(-m:m,-m:m,-m:m))
      sum_tot_all=0.d0
      dens=0.d0

      num_atom=0
      do 100 iatom_patch=1,natom_patch

      iatom=ind_config(iatom_patch)

      if(x(1,iatom).ge.1.d0) x(1,iatom)=x(1,iatom)-1.d0
      if(x(1,iatom).lt.0.d0) x(1,iatom)=x(1,iatom)+1.d0
c      xtmp=x(1,iatom)*nnodes
c      if(xtmp.lt.(inode-1).or.xtmp.ge.inode) goto 100
      if(mod(iatom_patch-1,nnodes).ne.inode-1) goto 100
      num_atom=num_atom+1    ! selected atom for this node

ccccccccccccccccccccccccccccccccccccccc
      if(iatype(iatom).eq.0) then
      write(6,*) "something wrong for iatype(iatom)",
     &     iatype(iatom),iatom
      call mpi_abort(MPI_COMM_WORLD,ierr)
      endif

      itype=ind_motif0(iatom_patch)
      op(:,:)=op_m(:,:,iatom_patch)
     

      Efield(1:3,iatom)=0.d0
      Efield(edim,iatom)=-1.d0
       Efieldtmp(1)=op(1,1)*Efield(1,iatom)+op(1,2)*Efield(2,iatom)+
     &        op(1,3)*Efield(3,iatom)
       Efieldtmp(2)=op(2,1)*Efield(1,iatom)+op(2,2)*Efield(2,iatom)+
     &        op(2,3)*Efield(3,iatom)
       Efieldtmp(3)=op(3,1)*Efield(1,iatom)+op(3,2)*Efield(2,iatom)+
     &        op(3,3)*Efield(3,iatom)
       write(6,*) "Efieldtmp=",Efieldtmp
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      motif_tmp(:,:,:)=Efieldtmp(1)*dens_mEM(:,:,:,1,itype)+
     &                 Efieldtmp(2)*dens_mEM(:,:,:,2,itype)+
     &                 Efieldtmp(3)*dens_mEM(:,:,:,3,itype)

cccccccccccccccccccccccccccccccccccccccccccccccc

       num_map=0
       sum_dens=0.d0
       sum_dens0=0.d0

       ix1=x(1,iatom)*n1
       ix2=x(2,iatom)*n2
       ix3=x(3,iatom)*n3
       do 300 ib=-ircut1,ircut1
       do 300 jb=-ircut2,ircut2
       do 300 kb=-ircut3,ircut3
       i=mod(ib+ix1+10*n1,n1)+1
       j=mod(jb+ix2+10*n2,n2)+1
       k=mod(kb+ix3+10*n3,n3)+1

       dx1=(ib+ix1)*1.d0/n1-x(1,iatom)
       dx2=(jb+ix2)*1.d0/n2-x(2,iatom)
       dx3=(kb+ix3)*1.d0/n3-x(3,iatom)

       dxt=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
       dyt=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
       dzt=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
       rr=dxt**2+dyt**2+dzt**2
       if(rr.lt.rrcut) then
       num_map=num_map+1
       if(num_map.gt.nmap) then
       write(6,*) "num_map.gt.nmap, stop"
       call mpi_abort(MPI_COMM_WORLD,ierr)
       endif
       ijk_map(1,num_map)=i  ! this is the i,j,k in vr(n1,n2,n3)
       ijk_map(2,num_map)=j
       ijk_map(3,num_map)=k

cccccccccccccccccccccccccccccccccccccc
cccc   convert the coordinates into motif coordinates
cccccccccccccccccccccccccccccccccccccc
       dx=op(1,1)*dxt+op(1,2)*dyt+op(1,3)*dzt
       dy=op(2,1)*dxt+op(2,2)*dyt+op(2,3)*dzt
       dz=op(3,1)*dxt+op(3,2)*dyt+op(3,3)*dzt
cccccccccccccccccccccccccccccccccccccc

       x1=dx*2*m/AL_mbox(1,1)
       x2=dy*2*m/AL_mbox(1,1)
       x3=dz*2*m/AL_mbox(1,1)

cccccccc x1,x2,x3 are the positions in the small box coordinates. 

       i0=x1+m*1.d0          ! bring x1 to positive for float->integer truncation
       j0=x2+m*1.d0
       k0=x3+m*1.d0
       i0=i0-m
       j0=j0-m
       k0=k0-m
       xi0=x1-i0
       xj0=x2-j0
       xk0=x3-k0

       call interpolate(vr00,i0,j0,k0,xi0,xj0,xk0,motif_tmp,m)
       dens_map(num_map)=vr00

       sum_dens=sum_dens+dens_map(num_map)
       sum_dens0=sum_dens0+abs(dens_map(num_map))
       endif        ! rrcut

300    continue        ! sweep all the i,j,k for the big box
       fact=sum_dens/sum_dens0


       do ii=1,num_map
       i=ijk_map(1,ii)
       j=ijk_map(2,ii)
       k=ijk_map(3,ii)
!      dens(i,j,k)=dens(i,j,k)+dens_map(ii)-abs(dens_map(ii))*fact
       dens(i,j,k)=dens(i,j,k)+dens_map(ii)
       enddo

100   continue         ! do all the atoms for this node 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      deallocate(motif_tmp)


      allocate(dens_tmp(n1,n2,n3))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call mpi_allreduce(dens,dens_tmp,n1*n2*n3,
     &  MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       dens=dens_tmp
       deallocate(dens_tmp)

      
!!!!! since it has been made 0-means for dens
!!!!! one has
!!!!! we assume AL is orthogonal box to avoid complicate math     
      if(inode.eq.1) then
             write(6,*) "sum dens =",sum(dens)
             write(6,*) "avg dens =",sum(dabs(dens))/(n1*n2*n3)
             write(6,*) "start to sum up" 
             !dens=dens-sum(dens)/(n1*n2*n3)
      do i=1,n1
      do j=1,n2
      do k=1,n3
      chi(1,edim)=chi(1,edim)+dens(i,j,k)*(i-1)*1.d0/n1*AL_au(1,1)
      chi(2,edim)=chi(2,edim)+dens(i,j,k)*(j-1)*1.d0/n2*AL_au(2,2)
      chi(3,edim)=chi(3,edim)+dens(i,j,k)*(k-1)*1.d0/n3*AL_au(3,3)
      enddo
      enddo
      enddo

      endif


      call  mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo ! end of edim

      if(inode.eq.1) then
      chi=chi*1.d0/(n1*n2*n3)
      pi=4*datan(1.d0)
      epsr(1,1)=1.d0
      epsr(2,2)=1.d0
      epsr(3,3)=1.d0
      epsr=epsr+4.d0*pi*chi

      write(6,*) "epsr"
      write(6,"(3F16.9)") epsr(:,1)
      write(6,"(3F16.9)") epsr(:,2)
      write(6,"(3F16.9)") epsr(:,3)
      endif

      call  mpi_barrier(MPI_COMM_WORLD,ierr)
      call  mpi_finalize(ierr)

      contains

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! what data does this subroutine manipulate ?
! iHflag
! ddcut
! dx10,dx20,dx30
! dx10Z,dx20Z,dx30Z
! dx1,dx2,dx3
! dx,dy,dz
! dx1Z,dx2Z,dx3Z
! dxZ,dyZ,dzZ
! dd              
      subroutine neighbore(dcut)
      implicit double precision (a-h,o-z)

cccccccccccccccccccccccccccccccccccccccccc
cccc The x, xZ relative position can not be shifted by 1 !

      do 500 i=1,natom
      iHflag=0
       if(iiatom(i).eq.105.or.iiatom(i).eq.115.
     & or.iiatom(i).eq.102.or.iiatom(i).eq.101.
     & or.iiatom(i).eq.1) iHflag=1       ! Hydrogen atom

      if(iHflag.eq.1) then
      ddcut=(dcut/2)**2
      else
      ddcut=dcut**2
      endif

      num=0 
      do j=1,natom
      if(i.ne.j) then
      dx10=x(1,j)-x(1,i)
      dx20=x(2,j)-x(2,i)
      dx30=x(3,j)-x(3,i)
      dx10Z=xZ(1,j)-xZ(1,i)
      dx20Z=xZ(2,j)-xZ(2,i)
      dx30Z=xZ(3,j)-xZ(3,i)
      do ii1=-1,1
      do ii2=-1,1
      do ii3=-1,1
      dx1=dx10+ii1
      dx2=dx20+ii2
      dx3=dx30+ii3
      dx=AL(1,1)*dx1+AL(1,2)*dx2+AL(1,3)*dx3
      dy=AL(2,1)*dx1+AL(2,2)*dx2+AL(2,3)*dx3
      dz=AL(3,1)*dx1+AL(3,2)*dx2+AL(3,3)*dx3
      dx1Z=dx10Z+ii1
      dx2Z=dx20Z+ii2
      dx3Z=dx30Z+ii3
      dxZ=ALZ(1,1)*dx1Z+ALZ(1,2)*dx2Z+ALZ(1,3)*dx3Z
      dyZ=ALZ(2,1)*dx1Z+ALZ(2,2)*dx2Z+ALZ(2,3)*dx3Z
      dzZ=ALZ(3,1)*dx1Z+ALZ(3,2)*dx2Z+ALZ(3,3)*dx3Z

      dd=dxZ**2+dyZ**2+dzZ**2

      if(dd.lt.ddcut) then
      num=num+1
      ind(num,i)=j
      dxyz_ind(1,num,i)=dx
      dxyz_ind(2,num,i)=dy
      dxyz_ind(3,num,i)=dz

      dxyz_indZ(1,num,i)=dxZ
      dxyz_indZ(2,num,i)=dyZ
      dxyz_indZ(3,num,i)=dzZ

      iat_neigh(num,i)= iiatom(j)

      endif

      enddo
      enddo
      enddo
      endif
      enddo

      num_neigh_All(i)=num
      
500   continue

      return
      end  subroutine neighbore


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      end

ccccccccccccccccccccccccccccccccccccccccccccccccC
cccccccccccccccccccccccccccccccccccccccccccccccccc
cccc The following subroutines are not contained subroutine, 
cccc so the parameters have to be passed in completely. 
ccccccccccccccccccccccccccccccccccccccccccccccccC
ccccccccccccccccccccccccccccccccccccccccccccccccc
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine interpolate(vr00,i0,j0,k0,xi0,xj0,xk0,
     &  motif_tmp,m)

       implicit double precision (a-h,o-z)
       real*8 vt(-1:2,-1:2,-1:2),vt1(-1:2,-1:2),vt2(-1:2)
       real*4 motif_tmp(-m:m,-m:m,-m:m)

       do it=-1,2
       itm=it+i0
       if(itm.gt.m) itm=m
       if(itm.lt.-m) itm=-m
       do jt=-1,2
       jtm=jt+j0
       if(jtm.gt.m) jtm=m
       if(jtm.lt.-m) jtm=-m
       do kt=-1,2
       ktm=kt+k0
       if(ktm.gt.m) ktm=m
       if(ktm.lt.-m) ktm=-m

       vt(it,jt,kt)=motif_tmp(itm,jtm,ktm)
       enddo
       enddo
       enddo

      do it=-1,2
      do jt=-1,2
      vt1(it,jt)=vt(it,jt,0)+(6*vt(it,jt,1)-2*vt(it,jt,-1)-
     &     vt(it,jt,2)-3*vt(it,jt,0))/6*xk0
     &  +(vt(it,jt,1)+vt(it,jt,-1)-2*vt(it,jt,0))/2*xk0**2
     &  +(vt(it,jt,2)+3*vt(it,jt,0)-3*vt(it,jt,1)
     &  -vt(it,jt,-1))/6*xk0**3
      enddo
      enddo

      do it=-1,2
      vt2(it)=vt1(it,0)+(6*vt1(it,1)-2*vt1(it,-1)-
     &      vt1(it,2)-3*vt1(it,0))/6*xj0
     &  +(vt1(it,1)+vt1(it,-1)-2*vt1(it,0))/2*xj0**2
     &  +(vt1(it,2)+3*vt1(it,0)-3*vt1(it,1)
     &  -vt1(it,-1))/6*xj0**3
      enddo

      vt3=vt2(0)+(6*vt2(1)-2*vt2(-1)-vt2(2)-3*vt2(0))/6*xi0
     &  +(vt2(1)+vt2(-1)-2*vt2(0))/2*xi0**2
     &  +(vt2(2)+3*vt2(0)-3*vt2(1)-vt2(-1))/6*xi0**3

       vr00=vt3

       return
       end 




!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       subroutine find_avg_E(AL,vr,n1,n2,n3,x,natom,
     &  size_Efield,Efield,iatype,inode,nnodes)

       implicit double precision (a-h,o-z)
       include 'mpif.h'

       real*8 vr(n1,n2,n3)
       real*8 AL(3,3),ALt(3,3),ALIt(3,3)
       real*8 x(3,natom),size_Efield(10),size_Efieldt(10)
       real*8 Efield(3,natom),Efieldt(3,natom)
       real*8 pi 
       integer iatype(natom)

!ccccccccccc want to get the atomic unit for electric field
       ALt=AL/0.529177 ! convert to Bohr
       size_Efieldt=size_Efield/0.529177
       pi=4*datan(1.d0)

       call get_ALI(ALt,ALIt)

       Efield=0.d0

       do 100 ii=1,natom

       if(mod(ii-1,nnodes).ne.inode-1) goto 99 

       dc=size_Efieldt(iatype(ii))

       rcutx1=dc*dsqrt(ALIt(1,1)**2+ALIt(2,1)**2+ALIt(3,1)**2)
       rcutx2=dc*dsqrt(ALIt(1,2)**2+ALIt(2,2)**2+ALIt(3,2)**2)
       rcutx3=dc*dsqrt(ALIt(1,3)**2+ALIt(2,3)**2+ALIt(3,3)**2)
       ircut1=n1*rcutx1+1
       ircut2=n2*rcutx2+1
       ircut3=n3*rcutx3+1


       sum1x=0.d0
       sum2x=0.d0
       sum1y=0.d0
       sum2y=0.d0
       sum1z=0.d0
       sum2z=0.d0

       xc1=mod(x(1,ii)+10,1.d0)
       xc2=mod(x(2,ii)+10,1.d0)
       xc3=mod(x(3,ii)+10,1.d0)

       ix1=xc1*n1
       ix2=xc2*n2
       ix3=xc3*n3

       num_map=0
       do 300 ib=-ircut1,ircut1
       do 300 jb=-ircut2,ircut2
       do 300 kb=-ircut3,ircut3

       i=mod(ib+ix1+10*n1,n1)+1
       j=mod(jb+ix2+10*n2,n2)+1
       k=mod(kb+ix3+10*n3,n3)+1

       dx1=(ib+ix1-xc1*n1)/n1
       dx2=(jb+ix2-xc2*n2)/n2
       dx3=(kb+ix3-xc3*n3)/n3

       dx=ALt(1,1)*dx1+ALt(1,2)*dx2+ALt(1,3)*dx3
       dy=ALt(2,1)*dx1+ALt(2,2)*dx2+ALt(2,3)*dx3
       dz=ALt(3,1)*dx1+ALt(3,2)*dx2+ALt(3,3)*dx3
       d=dsqrt(dx**2+dy**2+dz**2)
       if(d.lt.dc) then

       num_map=num_map+1

!       fact=cos(d/dc*pi/2)**2
!       sum1x=sum1x+fact*dx**2
!       sum2x=sum2x+vr(i,j,k)*fact*dx
!       sum1y=sum1y+fact*dy**2
!       sum2y=sum2y+vr(i,j,k)*fact*dy
!       sum1z=sum1z+fact*dz**2
!       sum2z=sum2z+vr(i,j,k)*fact*dz

        fact=cos(d/dc*pi/2)**2               
        sum1x=sum1x+fact*dx**2
        sum2x=sum2x+vr(i,j,k)*fact*dx
        sum1y=sum1y+fact*dy**2
        sum2y=sum2y+vr(i,j,k)*fact*dy
        sum1z=sum1z+fact*dz**2
        sum2z=sum2z+vr(i,j,k)*fact*dz

!! start revise by YZY
!       fact=2*cos((d/dc)**2*pi/2)*sin((d/dc)**2*pi/2)
!       ptx=dx
!       pty=dy
!       ptz=dz
!       sum1x=sum1x+fact*ptx*dx
!       sum2x=sum2x+vr(i,j,k)*fact*ptx
!       sum1y=sum1y+fact*pty*dy
!       sum2y=sum2y+vr(i,j,k)*fact*pty
!       sum1z=sum1z+fact*ptz*dz
!       sum2z=sum2z+vr(i,j,k)*fact*ptz
!! end revise by YZY
       endif

300    continue

       
       Efield(1,ii)=sum2x/sum1x   ! Hartree/Bohr, atomic unit
       Efield(2,ii)=sum2y/sum1y   ! Hartree/Bohr, atomic unit
       Efield(3,ii)=sum2z/sum1z   ! Hartree/Bohr, atomic unit

!       write(6,*) "Efield",ii,x(1,ii),Efield(1,ii)

99     continue
100    continue


      call mpi_allreduce(Efield,Efieldt,3*natom,
     $     MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
       Efield=Efieldt

       return
       end subroutine find_avg_E


      

