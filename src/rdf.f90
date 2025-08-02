        program rdf 
        implicit none
        integer:: natom,n_atomtypes 
        integer:: dimn,kAdimn,kBdimn
        integer:: i,j,k,m,l,nbin,r_dimn,norm_dimn
        real*8:: bx,by,bz,x,y,z,rmin,rmax,dr,ri,rf
        real*8:: rhoA_B,rhoA_A,rhoB_B
        real*8,allocatable,dimension(:,:):: xt
        real*8,allocatable,dimension(:):: frac_atomtypes,b_atomtypes
        real*8,allocatable,dimension(:,:):: xt_sub_A,xt_sub_B
        real*8,allocatable,dimension(:):: r_atom,distFunc,r_mid
        integer,allocatable,dimension(:):: kdimn
        character(len=5),allocatable,dimension(:):: ATOM,atom_types

        ! Reading input parameters 
        read(*,*) natom
        read(*,*) n_atomtypes
        allocate(atom_types(n_atomtypes))
        !allocate(frac_atomtypes(n_atomtypes))
        !allocate(b_atomtypes(n_atomtypes))
        !allocate(kdimn(n_atomtypes))
        read(*,*) (atom_types(i),i=1,n_atomtypes)
        read(*,*) rmin,rmax,dr
        !read(*,*) (frac_atomtypes(i),i=1,n_atomtypes)
        !read(*,*) (b_atomtypes(i),i=1,n_atomtypes) 
        
        ! *******************************************************
        allocate(ATOM(natom))
        dimn=3
        allocate(xt(natom,dimn))
        !********************************************************

        open(unit=20,file='Structure.xyz')
        read(20,*) natom
        read(20,*) bx, by, bz 
         
        ! Reading the positions
        do i=1,natom 
           read(20,*) ATOM(i),(xt(i,m),m=1,3)
        end do
        close(20)
        allocate(xt_sub_A(natom,dimn))
        allocate(xt_sub_B(natom,dimn))
        k=0
        do i=1,natom 
           if(ATOM(i).eq.atom_types(1)) then
              k=k+1   
              do m=1,3  
                 xt_sub_A(k,m)=xt(i,m)
              enddo   
           endif        
        enddo
        kAdimn=k
        k=0
        do i=1,natom 
           if(ATOM(i).eq.atom_types(2)) then
              k=k+1   
              do m=1,3  
                 xt_sub_B(k,m)=xt(i,m)
              enddo   
           endif        
        enddo
        kBdimn=k
        !write(*,*) "number of A and B", kAdimn,kBdimn 

        nbin=((rmax-rmin)/dr)+1 
        allocate(distFunc(nbin)) 
        allocate(r_mid(nbin))
        rhoA_B=(real(KBdimn))/(bx*by*bz)
        rhoA_A=(real(KAdimn-1))/(bx*by*bz)
        rhoB_B=(real(KBdimn-1))/(bx*by*bz) 
        ! Distance calculations A and B  
        r_dimn=KAdimn*KBdimn
        allocate(r_atom(r_dimn))
        l=0;
        do i=1,KAdimn
           do j=1,KBdimn 
              x=xt_sub_B(j,1)-xt_sub_A(i,1)
              y=xt_sub_B(j,2)-xt_sub_A(i,2)
              z=xt_sub_B(j,3)-xt_sub_A(i,3)
              x=x-(bx*Anint(x/bx))
              y=y-(by*Anint(y/by))
              z=z-(bz*Anint(z/bz))
              l=l+1
              r_atom(l)=sqrt(x**2+y**2+z**2)
           enddo
        end do
        !!!cal rdf !!!
        do i=1,nbin
              distFunc(i)=0.0
        enddo
        call cal_rdf(r_mid,rmin,rmax,dr,nbin,r_atom,r_dimn,rhoA_B,KAdimn,distFunc)
        open(unit=40,file="gofr-AB.txt")
        do i=1,nbin
           write(40,"(2(F12.6,1X))") r_mid(i),distFunc(i)
        enddo
        close(40)
        deallocate(r_atom)

        ! Distance calculations A and other A  
        r_dimn=KAdimn*(KAdimn-1)
        allocate(r_atom(r_dimn))
        l=0;
        do i=1,KAdimn
           do j=1,KAdimn
            if(i.ne.j) then  
               x=xt_sub_A(j,1)-xt_sub_A(i,1)
               y=xt_sub_A(j,2)-xt_sub_A(i,2)
               z=xt_sub_A(j,3)-xt_sub_A(i,3)
               x=x-(bx*Anint(x/bx))
               y=y-(by*Anint(y/by))
               z=z-(bz*Anint(z/bz))
               l=l+1
               r_atom(l)=sqrt(x**2+y**2+z**2)
            endif 
           enddo
        end do
        !!!cal rdf !!!
        do i=1,nbin
              distFunc(i)=0.0
        enddo
        norm_dimn=KAdimn-1
        call cal_rdf(r_mid,rmin,rmax,dr,nbin,r_atom,r_dimn,rhoA_A,norm_dimn,distFunc)
        open(unit=40,file="gofr-AA.txt")
        do i=1,nbin
           write(40,"(2(F12.6,1X))") r_mid(i),distFunc(i)
        enddo
        close(40)
        deallocate(r_atom)

        ! Distance calculations B and other B  
        r_dimn=KBdimn*(KBdimn-1)
        allocate(r_atom(r_dimn))
        l=0;
        do i=1,KBdimn
           do j=1,KBdimn
            if(i.ne.j) then  
               x=xt_sub_B(j,1)-xt_sub_B(i,1)
               y=xt_sub_B(j,2)-xt_sub_B(i,2)
               z=xt_sub_B(j,3)-xt_sub_B(i,3)
               x=x-(bx*Anint(x/bx))
               y=y-(by*Anint(y/by))
               z=z-(bz*Anint(z/bz))
               l=l+1
               r_atom(l)=sqrt(x**2+y**2+z**2)
            endif 
           enddo
        end do
        !!!cal rdf !!!
        do i=1,nbin
              distFunc(i)=0.0
        enddo
        norm_dimn=KBdimn-1
        call cal_rdf(r_mid,rmin,rmax,dr,nbin,r_atom,r_dimn,rhoB_B,norm_dimn,distFunc)
        open(unit=40,file="gofr-BB.txt")
        do i=1,nbin
           write(40,"(2(F12.6,1X))") r_mid(i),distFunc(i)
        enddo
        close(40) 

        deallocate(distFunc,r_mid,r_atom,xt_sub_A,xt_sub_B,xt,ATOM,atom_types) 
        end program rdf

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine cal_rdf(r_mid,rmin,rmax,dr,nbin,r_atom,r_dimn,rho,norm_dimn,distFunc)
        implicit none
        integer::i,j
        real*8::ri,rf
        integer,intent(in)::nbin,r_dimn,norm_dimn
        real*8,intent(in)::rmin,rmax,dr,rho
        real*8,dimension(:),intent(in)::r_atom(r_dimn)
        real*8,dimension(:),intent(out)::r_mid(nbin),distFunc(nbin)

        ri=rmin
        rf=ri+dr
        do i=1,nbin
           do j=1,r_dimn
             if(((ri<r_atom(j)).and.(r_atom(j)<rf)))then 
                  distFunc(i)=distFunc(i)+1.0
             endif
           enddo
          ri=rf
          rf=ri+dr 
        enddo   
        
        ri=rmin
        rf=ri+dr
        do i=1,nbin
           distFunc(i)=distFunc(i)/((real(norm_dimn))*(4.0*3.141592)*((ri+0.5*dr)*(ri+0.5*dr)*dr)*rho)
           r_mid(i)=ri+(0.5*dr)
           ri=rf
           rf=ri+dr 
        enddo
        end subroutine cal_rdf
        
        


