!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
! Hogil Kim Memorial Building #501 POSTECH,
! San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!E-mail: hoangdanhtai@gmail.com  Personal site: https://sites.google.com/site/hoangdanhtai
!-----------------------------------------------------------------------------------------!
!!! 2014.06.04: DYNAMIC IN PANCREATIC ISLET
!!! alpha: 11, beta: 12, delta: 13
!!! 2015.04.07: Sua lai chi dung cho Human, Mouse, Human with delta
!!! 2015.04.21: Bo sung tinh delta R,...
!!! 2015.04.24: Sua lai cach tinh R-->(R2): R*exp(2*i*Phi)= sum(exp(2*i*phi_i)
!!! 2015.04.27: Bo sung HCP structure
!!! 2015.05.04: Bo sung cau truc Human sau khi da remove 50% beta cells
!!! 06.05.2015: Bo sung cau truc Human sau khi da remove 70% beta cells
!!! 07.05.2015: Bo sung cau truc Human sau khi da remove 30%, 90% beta cells
!!! 12.05.2015: Bo sung cau truc moi HCP voi 725 cells: N501,N502,N503,....
!!! 15.05.2015: Bo sung 10 structure ensembles N503, N603, N703, N803, N903
!!! 16.05.2015: Bo sung 10 structure ensembles N502, N602, N702, N802, N902
!!! 20.05.2015: Bo sung HCP725_R (725 cells, random)
!!! 07.09.2015: Khong chia cho so luong Neighbor (o phan model) !!!!
!!! 14.09.2015: Bo sung structure sau khi remove 50% beta va reagregation
!!!                         (H01_50_R1 (nearest neighbor) va H01_50_R2 (random))
!!!  29.09.2015: Bo sung cau truc HCP1357 va HCP2493 voi cac Pb=0.6, 0.7, 0.8, 0.9
!!!  02.10.2015: Bo sung HCP725 voi Pb=0.6, 0.7, 0.8, 0.9 va Jaa=1, Jbb=1, Jab=2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
   PROGRAM main_BIO
   IMPLICIT NONE

   CHARACTER (LEN=150):: LOAD_PHI,NAME_DATA
   CHARACTER (LEN=3)  :: SAt
   REAL      (KIND=8),PARAMETER :: nul=0.,delt=0.01,aaa=1.0,Kaa=1.,Kbb=1.,Kdd=1.

   CHARACTER (256)    :: Ligne20,Ligne21,Ligne22,Ligne23,Ligne40

   INTEGER (KIND=8) :: i,j,i_n,nn_total,k,i_n2,tmp2
   INTEGER (KIND=8):: i_loop,i_loop1,i_av,n_equi1,n_equi2,n_average,i_times,n_times
   INTEGER (KIND=8) :: time_2,time_3,time_5,time_6,time_7,natom,na,nb
   
   REAL (KIND=8) :: n_atom,number_a,number_b,number_d,rdn_s,Kab,Kba,Kad,Kda,Kbd,Kdb,r
   REAL (KIND=8) :: pi,delt1,Zax,Zay,Za,Zbx,Zby,Zb,Za_av,Zb_av,Zdx,Zdy,Zd,Zd_av,Sabd_load
   REAL (KIND=8) :: Zx,Zy,Z_total,Z_total_av,run_time,S_tmp,S_load,WS,WS_av
   REAL (KIND=8) :: phi_a,phi_b,phi_d,delta_ab,delta_ad,Pb,rdn_config
   REAL (KIND=8) :: r0,tmp1,r02,alpha0,alpha,OA2,OB2,AB2

   INTEGER (KIND=8),DIMENSION(:),ALLOCATABLE :: nn,nn1
   INTEGER (KIND=8),DIMENSION(:,:),ALLOCATABLE :: name_in
   REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: S,delS,S1,x,y,z,Sabd
   REAL (KIND=8),DIMENSION(:,:),ALLOCATABLE :: Ks,Ks1,r3,a

   !!! Khai bao phan Histogram
   INTEGER (KIND=8),PARAMETER :: n_histo=100
   INTEGER (KIND=8) :: i_histo
   REAL    (KIND=8) :: del_histo
   REAL (KIND=8),DIMENSION(0:n_histo) :: S_histo,H1_histo,H2_histo,H3_histo
     
!!!=======================================================================================   
!!!=======================================================================================
   CALL system('rm config_ini_3D.pdb')
   CALL system('rm *.dat*')

   CALL ini_rdm_number()
   CALL read_input_parameter()
   
   !!! Set Kba, Kda: 
   Kba=-1./Kab ; Kda=-1./Kad ; Kbd=Kad/Kab ; Kdb=-1./Kbd

   CALL open_data()

   n_atom=real(natom)

   ALLOCATE(S(natom),delS(natom),S1(natom))
   ALLOCATE(nn(natom))
   ALLOCATE(x(natom),y(natom),z(natom),Sabd(natom))
   ALLOCATE(name_in(natom,0:natom))
   ALLOCATE(Ks(natom,natom),Ks1(natom,natom))
   ALLOCATE(r3(natom,0:natom))
   ALLOCATE(a(natom,0:natom))
   ALLOCATE(nn1(0:natom))

   IF ((NAME_DATA=='HCP').or.(NAME_DATA=='HCP725_R')) THEN
      CALL read_data_random()
   ELSE
      CALL read_data()
   END IF
   !!WRITE(*,*)number_a,number_b
   
   r0=sqrt(r02)
   
   delt1=delt+(1./2.)*delt**2. ! +(1./6.)*delt**3.+(1./24.)*delt**4. !! 2nd order
   
   pi=acos(-1.)

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   OPEN(unit=15,file='nn_i_cell.dat') 
   OPEN(unit=16,file='nn_i_final.dat')
   
   OPEN(unit=20,file='phi_equi.dat')
   OPEN(unit=21,file='rtime.dat')
   OPEN(unit=22,file='average.dat')
   OPEN(unit=23,file='phi_time.dat')

   IF (LOAD_PHI=='YES') THEN
      WRITE(*,*)'Load phi equi'
      CALL load_phi_equi()
   ELSE
      CALL generate_ini_s()
   END IF

   CALL write_config_ini_3D()
  
   !!!---------- Cho experiment data : ------------
   CALL nearest_neighbor1()
   CALL nearest_neighbor2()
   !!!----------------------------------------------
   CALL Ks_value()   

   DO i=1,natom
      DO j=1,natom
      !!! Ks1(i,j)=Ks(i,j)/nn(i)
      !! 07.09.2015: khong chia cho so neighbor
      Ks1(i,j)=Ks(i,j)
      
      END DO
   END DO
   
   CALL equi_lattice1()
   WRITE(*,*)'Finish Equilibre'

   CALL average_thermal()

   CALL save_phi_equi()

   CALL computation_time()
   
   CLOSE(15)
   CLOSE(16)

   CONTAINS
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE init_rdm_number()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE ini_rdm_number()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time
   INTEGER (KIND=8),DIMENSION(50) :: seed

   CALL DATE_AND_TIME(values=time)     ! Get the current time
   seed(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
   CALL RANDOM_SEED(PUT=seed)

   time_2=time(2) ; time_3=time(3) ; time_5=time(5) ; time_6=time(6) ; time_7=time(7)

   END SUBROUTINE ini_rdm_number
   
!!!=======================================================================================  
   SUBROUTINE computation_time()
   IMPLICIT NONE

   INTEGER (KIND=8),DIMENSION(8) :: time1
   INTEGER (KIND=8) :: run_date,run_hour,run_minute,run_second

   OPEN (90,file='time_run.dat')

   CALL DATE_AND_TIME(values=time1)     ! Get the current time
      
   run_time = (time1(2)-time_2)*1296000.+(time1(3)-time_3)*86400.&
             +(time1(5)-time_5)*3600.+(time1(6)-time_6)*60.+(time1(7)-time_7) !! second (s)
   run_date=int(run_time/86400.)
   run_hour=int(run_time/3600.-run_date*24.)
   run_minute=int(run_time/60.-run_date*1440.-run_hour*60.)
   run_second=int(run_time-run_date*86400.-run_hour*3600.-run_minute*60.)

   WRITE(90,*)'run_date  :',run_date
   WRITE(90,*)'run_hour  :',run_hour
   WRITE(90,*)'run_minute:',run_minute
   WRITE(90,*)'run_second:',run_second

   CLOSE(90)

   END SUBROUTINE computation_time   

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE read_input_parameter() 
!!! OPEN the parameter from file "parameter.in"
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE read_input_parameter()
   IMPLICIT NONE

   CHARACTER (LEN=150) :: tamp
   OPEN(11,file='1parameter.in')   
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(A15))')   tamp, NAME_DATA
   READ(11, '(A30,(F12.6))') tamp, Pb
   READ(11, '(A30,(A10))')   tamp, LOAD_PHI
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, Kab
   READ(11, '(A30,(F12.6))') tamp, Kad
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(I12))')   tamp,n_equi1
   READ(11, '(A30,(I12))')   tamp,n_equi2
   READ(11, '(A30,(I12))')   tamp,n_average
   READ(11, '(A30,(I12))')   tamp,n_times
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp

   CLOSE(11) 

   END SUBROUTINE read_input_parameter

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE equi_lattice1()
   IMPLICIT NONE

   DO i_loop=1,n_equi1
      CALL equi_lattice()
   END DO

   END SUBROUTINE equi_lattice1
    
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE generate_ini_s() !!! Generate initial value for S:
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE generate_ini_s()
   IMPLICIT NONE 
   
   DO i=1,natom

      CALL random_number(rdn_s)
      !S(i)=2.*pi*rdn_s       !! random
      
      S(i)=(rdn_s-0.5)*pi      !! around pi/2
      !S(i)=(rdn_s-0.5)*pi/2.  !! around pi/4
      !S(i)=(rdn_s-0.5)*pi/4.  !! around pi/8
      
      !S(i)=0.

   ENDDO

   END SUBROUTINE generate_ini_s

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE save_phi_equi()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE save_phi_equi()
   IMPLICIT NONE 
   
   DO i=1,natom
      WRITE(Ligne20,*)i,Sabd(i),S(i)
      WRITE(20,'(a)') trim(Ligne20)
   ENDDO

   CLOSE(20)

   END SUBROUTINE save_phi_equi

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE load_phi_equi()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE load_phi_equi()
   IMPLICIT NONE 

   OPEN(unit=30,file='phi_equi.txt')

   DO i_loop=1,natom
      READ(30,*)i,Sabd_load,S_load
      Sabd(i)=Sabd_load
      S(i)=S_load
 
   END DO

   END SUBROUTINE load_phi_equi

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! WRITE initial position configuration in 3D : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE write_config_ini_3D()
   IMPLICIT NONE
 
   OPEN(unit=12,file='config_ini_3D.pdb')

   DO i=1,natom
      IF (Sabd(i)==11) THEN 
         SAt='C'  !! white , alpha
         WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul
      ELSE
                                        
      IF (Sabd(i)==12) THEN 
         SAt='S'  !! yellow , beta
         WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul
      ELSE
           SAt='Cd'  !! pink, delta
         WRITE(12,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            SAt,x(i)*aaa,y(i)*aaa,z(i)*aaa,nul

      END IF
      END IF

   END DO

   CLOSE(12)

   END SUBROUTINE write_config_ini_3D

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE equi_lattice1()
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE equi_lattice()
   IMPLICIT NONE
     
   DO i=1,natom
      S_tmp=S(i)
      delS(i)=0.

      DO j=1,nn(i)
      delS(i)=delS(i)+Ks1(i,j)*sin(S(name_in(i,j))-S_tmp)
      END DO

   END DO

   DO i=1,natom
      S(i)=S(i)+delS(i)*delt1
   END DO

   END SUBROUTINE equi_lattice
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE value_thermal() : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE value_thermal()
   IMPLICIT NONE

!!!===========================================================
   DO i=1,natom
      S_tmp=S(i)
      delS(i)=0.
      DO j=1,nn(i)
      delS(i)=delS(i)+Ks1(i,j)*sin(S(name_in(i,j))-S_tmp)
      END DO

   END DO

   DO i=1,natom
      S(i)=S(i)+delS(i)*delt1
   END DO

!!!===========================================================
   Zax=0. ; Zay=0. ; Zbx=0. ; Zby=0. ;  Zdx=0. ; Zdy=0. ; WS=0.
   DO i=1,natom
         S_tmp=S(i)
         !!! Order parameter:     
         
         !!!! 2014.04.24: Luu y da thay doi lai cach tinh R2=sin(2*phi)
         IF (Sabd(i)<11.5) THEN
            Zax=Zax+cos(2.*S_tmp) 
            Zay=Zay+sin(2.*S_tmp)
         ELSE
            IF (Sabd(i)<12.5) THEN
            Zbx=Zbx+cos(2.*S_tmp) 
            Zby=Zby+sin(2.*S_tmp)
            ELSE
            Zdx=Zdx+cos(2.*S_tmp) 
            Zdy=Zdy+sin(2.*S_tmp)
            END IF
         END IF

         !!!----------------------
         !!! Wave speed:
         WS=WS+delS(i)

   END DO
   Zx=(Zax+Zbx+Zdx)/n_atom ; Zy=(Zay+Zby+Zdy)/n_atom
   Zax=Zax/number_a ; Zay=Zay/number_a
   Zbx=Zbx/number_b ; Zby=Zby/number_b
   Zdx=Zdx/number_d ; Zdy=Zdy/number_d

   Za=sqrt(Zax**2.+Zay**2.)
   Zb=sqrt(Zbx**2.+Zby**2.)
   Zd=sqrt(Zdx**2.+Zdy**2.)
   Z_total=sqrt(Zx**2.+Zy**2.)
   WS=WS/n_atom
     
   END SUBROUTINE value_thermal

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!   SUBROUTINE histogram() : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE histogram()
   IMPLICIT NONE

   DO i=1,natom

      S1(i)=S(i)-int(S(i)/2./pi)*2.*pi
   
      IF (S1(i)<0.) THEN
         S1(i)=S1(i)+2.*pi
      END IF

      i_histo=int(S1(i)/del_histo)+1

      IF (Sabd(i)<11.5) THEN
         H1_histo(i_histo)=H1_histo(i_histo)+1.
      ELSE
         IF (Sabd(i)<12.5) THEN    
            H2_histo(i_histo)=H2_histo(i_histo)+1.
         ELSE
            H3_histo(i_histo)=H3_histo(i_histo)+1.
         END IF
      END IF

   END DO

   END SUBROUTINE histogram

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE average_thermal() : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE average_thermal()
   IMPLICIT NONE   

   Za_av=0. ; Zb_av=0. ; Zd_av=0. ; Z_total_av=0. ; WS_av=0.

   DO i_times=1,n_times
      DO i_loop1=1,n_equi2
         CALL equi_lattice()
      END DO

      DO i_av=1,n_average

         CALL value_thermal()
         CALL average_phase()

         IF (mod(i_av,100)==0) THEN
         WRITE(Ligne21,*)i_av,Za,Zb,Zd,Z_total
         WRITE(21,'(a)') trim(Ligne21)

         WRITE(Ligne23,*)i_av,phi_a,phi_b,phi_d,delta_ab,delta_ad
         WRITE(23,'(a)') trim(Ligne23)

         END IF     
                  
         Za_av=Za_av+Za
         Zb_av=Zb_av+Zb
         Zd_av=Zd_av+Zd
         Z_total_av=Z_total_av+Z_total
         WS_av=WS_av+WS                     

      END DO
      
   END DO

   Za_av=Za_av/real(n_times*n_average)
   Zb_av=Zb_av/real(n_times*n_average)
   Zd_av=Zd_av/real(n_times*n_average)
   Z_total_av=Z_total_av/real(n_times*n_average)
   WS_av=WS_av/real(n_times*n_average)

   WRITE(Ligne22,*)Kad,Kab,Za_av,Zb_av,Zd_av,Z_total_av,WS_av
   WRITE(22,'(a)') trim(Ligne22)

   CLOSE(21)
   CLOSE(22)
   CLOSE(23)

   !!! Histogram
   !!del_histo = 2.*pi/(n_histo-1)
   del_histo = 2.*pi/n_histo
   H1_histo(:)=0. ; H2_histo(:)=0. ; H3_histo(:)=0.

   WRITE(*,*)'Running histogram'
   CALL histogram()
   OPEN(unit=40,file='histogram.dat')
   DO i_histo=1,n_histo
      H1_histo(i_histo)=H1_histo(i_histo)/number_a
      H2_histo(i_histo)=H2_histo(i_histo)/number_b
      H3_histo(i_histo)=H3_histo(i_histo)/number_d     
      S_histo(i_histo)=(i_histo-0.5)*del_histo
      WRITE(Ligne40,*)S_histo(i_histo),H1_histo(i_histo),H2_histo(i_histo),H3_histo(i_histo)
      WRITE(40,'(a)') trim(Ligne40)
   ENDDO
   CLOSE(40)
      
   END SUBROUTINE average_thermal


!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE nearest_neighbor1() : find nearest neighbor j of every i and number of nn
!!!! Su dung cho cau truc chuan (HCP,...), khong sung dung duoc cho experiment data
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE nearest_neighbor1()
   IMPLICIT NONE
     
   !!!--------------------------------------------------------------
   nn_total=0 ;  nn(:)=0
      
   DO i=1,natom
   DO j=1,natom
      r=sqrt((x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.)
      IF ((0.1<r).and.(r<r0)) THEN
         nn(i)=nn(i)+1
      END IF

   END DO
                                  !!! nn(i): number neigbors of i cell
   nn_total=nn_total+nn(i)        !!! number neigbors of system
   WRITE(15,*)i,nn(i)

   END DO

  !! WRITE(*,*)'nn_total=',nn_total/2

   !!!--------------------------------------------------------
   !!! Assign name for nearest neighbor of cell i
           
   name_in(:,:)=0
      
   DO i=1,natom
      i_n=0           
      DO j=1,natom
      r=sqrt((x(i)-x(j))**2.+(y(i)-y(j))**2.+(z(i)-z(j))**2.)
         IF ((0.1<r).and.(r<r0)) THEN
            i_n=i_n+1
            name_in(i,i_n)=j         

         END IF
      END DO
   END DO

   END SUBROUTINE nearest_neighbor1

!!!=======================================================================================   
   SUBROUTINE nearest_neighbor2()
   IMPLICIT NONE

   !!!---------------------------------------------------------
   !!! Tim khoang cach r3(i,i_n)
   DO i=1,natom
      DO i_n=1,nn(i) 
      r3(i,i_n)=(x(i)-x(name_in(i,i_n)))**2.+(y(i)-y(name_in(i,i_n)))**2.+(z(i)-z(name_in(i,i_n)))**2.
      END DO
   END DO

   !!!---------------------------------------------------------
   !!! Sap xep theo thu tu tang dan khoang cach r3(i,i_n)
   DO i=1,natom
      DO i_n=1,nn(i)-1
         k=i_n

         DO j=i_n+1,nn(i)
            IF (r3(i,j)<r3(i,k)) THEN
               k=j
            END IF
         END DO

         IF (k/=i_n) THEN
            tmp1=r3(i,i_n)
            r3(i,i_n)=r3(i,k)
            r3(i,k)=tmp1

            tmp2=name_in(i,i_n)
            name_in(i,i_n)=name_in(i,k)
            name_in(i,k)=tmp2

         END IF

      END DO
   END DO

   !!!---------------------------------------------------------
   !!! Kiem tra goc de tim cac nearest neighbor
   a(:,:)=1
   DO i=1,natom
   DO i_n=1,nn(i)-1

   IF (a(i,i_n)==1) THEN
   DO i_n2=i_n+1,nn(i)

   OA2=(x(i)-x(name_in(i,i_n)))**2.+(y(i)-y(name_in(i,i_n)))**2.+(z(i)-z(name_in(i,i_n)))**2.
   OB2=(x(i)-x(name_in(i,i_n2)))**2.+(y(i)-y(name_in(i,i_n2)))**2.+(z(i)-z(name_in(i,i_n2)))**2.
   AB2=(x(name_in(i,i_n))-x(name_in(i,i_n2)))**2.+(y(name_in(i,i_n))-y(name_in(i,i_n2)))**2.&
        +(z(name_in(i,i_n))-z(name_in(i,i_n2)))**2.

   alpha=acos((OA2+OB2-AB2)/(2.*sqrt(OA2)*sqrt(OB2)))*180./pi

   IF (alpha<alpha0) THEN
      a(i,i_n2)=0
   END IF

   END DO
   END IF

   END DO
   END DO

   !!!---------------------------------------------------------
   !!! Sap xep lai (theo thu tu khoang cach tang dan) sau khi loai bo
   DO i=1,natom
      k=0
      DO i_n=1,nn(i)

      IF (a(i,i_n)==1) THEN
         k=k+1
         r3(i,k)=r3(i,i_n)
         name_in(i,k)=name_in(i,i_n)
      END IF
      END DO
            
      !! Update nn(i):
      nn(i)=k            

   END DO
                                 
   !!! Number neighbors
   nn1(:)=0

   DO j=1,MAXVAL(nn)
      DO i=1,natom
      IF (nn(i)==j) THEN
         nn1(j)=nn1(j)+1            
      END IF
      END DO
   WRITE(16,*)j,nn1(j)
   !WRITE(*,*)j,nn1(j)

   END DO

   END SUBROUTINE nearest_neighbor2
   
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE Ks_value : Assign Ks=Kab,Kba,Kaa,Kbb,... : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE Ks_value()
   IMPLICIT NONE
     
   Ks(:,:)=0.

   DO i=1,natom
   DO j=1,nn(i)

      !!! For a:
      IF (Sabd(i)==11) THEN
         IF (Sabd(name_in(i,j))==11) THEN
            Ks(i,j)=Kaa
         END IF

         IF (Sabd(name_in(i,j))==12) THEN
            Ks(i,j)=Kba
         END IF

         IF (Sabd(name_in(i,j))==13) THEN
            Ks(i,j)=Kda
         END IF

      END IF

      !!! For b:
      IF (Sabd(i)==12) THEN
         IF (Sabd(name_in(i,j))==11) THEN
            Ks(i,j)=Kab
         END IF

         IF (Sabd(name_in(i,j))==12) THEN
            Ks(i,j)=Kbb
         END IF

         IF (Sabd(name_in(i,j))==13) THEN
            Ks(i,j)=Kdb
         END IF
      
      END IF

      !!! For d:
      IF (Sabd(i)==13) THEN
         IF (Sabd(name_in(i,j))==11) THEN
            Ks(i,j)=Kad
         END IF

         IF (Sabd(name_in(i,j))==12) THEN
            Ks(i,j)=Kbd
         END IF

         IF (Sabd(name_in(i,j))==13) THEN
            Ks(i,j)=Kdd
         END IF
      
      END IF

   END DO
   END DO

   END SUBROUTINE Ks_value

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!  SUBROUTINE Open data and Read_data() Read experimental data :
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE open_data()
   IMPLICIT NONE

!!!==================================================================  
!!! 02.10.2015:
   !!! HCP725_60M mixing, N=725, pb=0.6, Jaa=1, Jbb=1, Jab=2 :
   IF (name_data=='HCP725_60M') THEN     
      OPEN(unit=12,file='0HCP725_60M.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP725_70M mixing, N=725, pb=0.7, Jaa=1, Jbb=1, Jab=2 :
   IF (name_data=='HCP725_70M') THEN     
      OPEN(unit=12,file='0HCP725_70M.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP725_80M mixing, N=725, pb=0.8, Jaa=1, Jbb=1, Jab=2 :
   IF (name_data=='HCP725_80M') THEN     
      OPEN(unit=12,file='0HCP725_80M.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP725_90M mixing, N=725, pb=0.9, Jaa=1, Jbb=1, Jab=2 :
   IF (name_data=='HCP725_90M') THEN     
      OPEN(unit=12,file='0HCP725_90M.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 29.09.2015: Bo sung HCP1357 va HCP2493 voi cac Pb khac nhau
   !!! HCP1357_60S shell-core, N=1357, pb=0.6, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP1357_60S') THEN     
      OPEN(unit=12,file='0HCP1357_60S.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP1357_60M partial mixing, N=1357, pb=0.6, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP1357_60M') THEN     
      OPEN(unit=12,file='0HCP1357_60M.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP1357_70S shell-core, N=1357, pb=0.7, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP1357_70S') THEN     
      OPEN(unit=12,file='0HCP1357_70S.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP1357_70M partial mixing, N=1357, pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP1357_70M') THEN     
      OPEN(unit=12,file='0HCP1357_70M.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! HCP1357_80S shell-core, N=1357, pb=0.8, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP1357_80S') THEN     
      OPEN(unit=12,file='0HCP1357_80S.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP1357_80M partial mixing, N=1357, pb=0.8, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP1357_80M') THEN     
      OPEN(unit=12,file='0HCP1357_80M.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF   

   !!! HCP1357_90S shell-core, N=1357, pb=0.9, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP1357_90S') THEN     
      OPEN(unit=12,file='0HCP1357_90S.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP1357_90M partial mixing, N=1357, pb=0.9, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP1357_90M') THEN     
      OPEN(unit=12,file='0HCP1357_90M.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF  
   
!!!------------------------------------------------------------------------------------
   !!! HCP2493_60S shell-core, N=2493, pb=0.6, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP2493_60S') THEN     
      OPEN(unit=12,file='0HCP2493_60S.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP2493_60M partial mixing, N=2493, pb=0.6, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_60M') THEN     
      OPEN(unit=12,file='0HCP2493_60M.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP2493_70S shell-core, N=2493, pb=0.7, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP2493_70S') THEN     
      OPEN(unit=12,file='0HCP2493_70S.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP2493_70M partial mixing, N=2493, pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_70M') THEN     
      OPEN(unit=12,file='0HCP2493_70M.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! HCP2493_80S shell-core, N=2493, pb=0.8, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP2493_80S') THEN     
      OPEN(unit=12,file='0HCP2493_80S.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP2493_80M partial mixing, N=2493, pb=0.8, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_80M') THEN     
      OPEN(unit=12,file='0HCP2493_80M.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF   

   !!! HCP2493_90S shell-core, N=2493, pb=0.9, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP2493_90S') THEN     
      OPEN(unit=12,file='0HCP2493_90S.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP2493_90M partial mixing, N=2493, pb=0.9, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_90M') THEN     
      OPEN(unit=12,file='0HCP2493_90M.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF  

!!!!!=================================================================================
!!! == 14.09.2015: Human go bo 50% beta cells, sau do reagregation nearest neighbor====
      IF (name_data=='H01_50_R1') THEN     
            OPEN(unit=12,file='0H01_50_R1.txt')
            natom=1392 ; r02= 350.0; alpha0= 44.2  
      END IF   
      IF (name_data=='H02_50_R1') THEN     
            OPEN(unit=12,file='0H02_50_R1.txt')
            natom=690 ; r02= 350.0; alpha0= 43.7 
      END IF    
      IF (name_data=='H03_50_R1') THEN     
            OPEN(unit=12,file='0H03_50_R1.txt')
            natom=1010 ; r02= 350.0; alpha0= 42.5 
      END IF  
      IF (name_data=='H04_50_R1') THEN     
            OPEN(unit=12,file='0H04_50_R1.txt')
            natom=382 ; r02= 350.0 ; alpha0= 41.9
      END IF  
      IF (name_data=='H05_50_R1') THEN     
            OPEN(unit=12,file='0H05_50_R1.txt')
            natom=798 ; r02= 360.0; alpha0= 40.5
      END IF
      IF (name_data=='H06_50_R1') THEN     
            OPEN(unit=12,file='0H06_50_R1.txt')
            natom=212 ; r02= 350.0; alpha0= 40.0 
      END IF      
      IF (name_data=='H07_50_R1') THEN     
            OPEN(unit=12,file='0H07_50_R1.txt')
            natom=1817 ; r02= 350.0; alpha0= 42.0 
      END IF 
      IF (name_data=='H08_50_R1') THEN     
            OPEN(unit=12,file='0H08_50_R1.txt')
            natom=1082 ; r02= 350.0; alpha0= 43.0 
      END IF   
      IF (name_data=='H09_50_R1') THEN     
            OPEN(unit=12,file='0H09_50_R1.txt')
            natom=512 ; r02= 350.0; alpha0= 42.0
      END IF
      IF (name_data=='H10_50_R1') THEN     
            OPEN(unit=12,file='0H10_50_R1.txt')
            natom=1410 ; r02= 350.0; alpha0= 42.21 
      END IF           
      IF (name_data=='H11_50_R1') THEN     
            OPEN(unit=12,file='0H11_50_R1.txt')
            natom=544 ; r02= 350.0; alpha0= 41.2
      END IF     
      IF (name_data=='H12_50_R1') THEN     
            OPEN(unit=12,file='0H12_50_R1.txt')
            natom=271 ; r02= 350.0; alpha0= 42.54
      END IF      
      IF (name_data=='H13_50_R1') THEN     
            OPEN(unit=12,file='0H13_50_R1.txt')
            natom=1291 ; r02= 350.0; alpha0= 42.6
      END IF
      IF (name_data=='H14_50_R1') THEN     
            OPEN(unit=12,file='0H14_50_R1.txt')
            natom=347 ; r02= 360.0; alpha0= 40.2
      END IF     
      IF (name_data=='H15_50_R1') THEN     
            OPEN(unit=12,file='0H15_50_R1.txt')
            natom=1090 ; r02= 350.0; alpha0= 41.9
      END IF 
      IF (name_data=='H16_50_R1') THEN     
            OPEN(unit=12,file='0H16_50_R1.txt')
            natom=473 ; r02= 350.0; alpha0= 42.1
      END IF       
      IF (name_data=='H17_50_R1') THEN     
            OPEN(unit=12,file='0H17_50_R1.txt')
            natom=295 ; r02= 350.0; alpha0= 42.1
      END IF 
      IF (name_data=='H18_50_R1') THEN     
            OPEN(unit=12,file='0H18_50_R1.txt')
            natom=889 ; r02= 350.0; alpha0= 43.9
      END IF       
      IF (name_data=='H19_50_R1') THEN     
            OPEN(unit=12,file='0H19_50_R1.txt')
            natom=1134 ; r02= 350.0; alpha0= 43.5
      END IF       
      IF (name_data=='H20_50_R1') THEN     
            OPEN(unit=12,file='0H20_50_R1.txt')
            natom=911 ; r02= 350.0; alpha0= 42.61
      END IF
      IF (name_data=='H21_50_R1') THEN     
            OPEN(unit=12,file='0H21_50_R1.txt')
            natom=1259 ; r02= 350.0; alpha0= 41.1
      END IF     
      IF (name_data=='H22_50_R1') THEN     
            OPEN(unit=12,file='0H22_50_R1.txt')
            natom=604 ; r02= 350.0; alpha0= 40.8
      END IF      
      IF (name_data=='H23_50_R1') THEN     
            OPEN(unit=12,file='0H23_50_R1.txt')
            natom=407 ; r02= 350.0; alpha0= 41.97
      END IF
      IF (name_data=='H24_50_R1') THEN     
            OPEN(unit=12,file='0H24_50_R1.txt')
            natom=4000 ; r02= 350.0; alpha0= 43.01
      END IF    
      IF (name_data=='H25_50_R1') THEN     
            OPEN(unit=12,file='0H25_50_R1.txt')
            natom=437 ; r02= 350.0; alpha0= 41.8
      END IF 
      IF (name_data=='H26_50_R1') THEN     
            OPEN(unit=12,file='0H26_50_R1.txt')
            natom=981 ; r02= 350.0; alpha0= 42.5
      END IF       
      IF (name_data=='H27_50_R1') THEN     
            OPEN(unit=12,file='0H27_50_R1.txt')
            natom=1452 ; r02= 350.0; alpha0= 42.71
      END IF 
      IF (name_data=='H28_50_R1') THEN     
            OPEN(unit=12,file='0H28_50_R1.txt')
            natom=1825 ; r02= 350.0; alpha0= 42.71
      END IF       
      IF (name_data=='H29_50_R1') THEN     
            OPEN(unit=12,file='0H29_50_R1.txt')
            natom=1453 ; r02= 350.0; alpha0= 42.82
      END IF       
      IF (name_data=='H30_50_R1') THEN     
            OPEN(unit=12,file='0H30_50_R1.txt')
            natom=792 ; r02= 350.0; alpha0= 42.5
      END IF  

!!!!!=================================================================================
!!! == 14.09.2015: Human sau khi da go bo 50% beta cells, sau do reagregation random ====

      IF (name_data=='H01_50_R2') THEN     
            OPEN(unit=12,file='0H01_50_R2.txt')
            natom=1392 ; r02= 350.0; alpha0= 44.2      
      END IF   
      IF (name_data=='H02_50_R2') THEN     
            OPEN(unit=12,file='0H02_50_R2.txt')
            natom=690 ; r02= 350.0; alpha0= 43.7 
      END IF    
      IF (name_data=='H03_50_R2') THEN     
            OPEN(unit=12,file='0H03_50_R2.txt')
            natom=1010 ; r02= 350.0; alpha0= 42.5 
      END IF  
      IF (name_data=='H04_50_R2') THEN     
            OPEN(unit=12,file='0H04_50_R2.txt')
            natom=382 ; r02= 350.0 ; alpha0= 41.9
      END IF  
      IF (name_data=='H05_50_R2') THEN     
            OPEN(unit=12,file='0H05_50_R2.txt')
            natom=798 ; r02= 360.0; alpha0= 40.5
      END IF
      IF (name_data=='H06_50_R2') THEN     
            OPEN(unit=12,file='0H06_50_R2.txt')
            natom=212 ; r02= 350.0; alpha0= 40.0 
      END IF      
      IF (name_data=='H07_50_R2') THEN     
            OPEN(unit=12,file='0H07_50_R2.txt')
            natom=1817 ; r02= 350.0; alpha0= 42.0 
      END IF 
      IF (name_data=='H08_50_R2') THEN     
            OPEN(unit=12,file='0H08_50_R2.txt')
            natom=1082 ; r02= 350.0; alpha0= 43.0 
      END IF   
      IF (name_data=='H09_50_R2') THEN     
            OPEN(unit=12,file='0H09_50_R2.txt')
            natom=512 ; r02= 350.0; alpha0= 42.0
      END IF
      IF (name_data=='H10_50_R2') THEN     
            OPEN(unit=12,file='0H10_50_R2.txt')
            natom=1410 ; r02= 350.0; alpha0= 42.21 
      END IF           
      IF (name_data=='H11_50_R2') THEN     
            OPEN(unit=12,file='0H11_50_R2.txt')
            natom=544 ; r02= 350.0; alpha0= 41.2
      END IF     
      IF (name_data=='H12_50_R2') THEN     
            OPEN(unit=12,file='0H12_50_R2.txt')
            natom=271 ; r02= 350.0; alpha0= 42.54
      END IF      
      IF (name_data=='H13_50_R2') THEN     
            OPEN(unit=12,file='0H13_50_R2.txt')
            natom=1291 ; r02= 350.0; alpha0= 42.6
      END IF
      IF (name_data=='H14_50_R2') THEN     
            OPEN(unit=12,file='0H14_50_R2.txt')
            natom=347 ; r02= 360.0; alpha0= 40.2
      END IF     
      IF (name_data=='H15_50_R2') THEN     
            OPEN(unit=12,file='0H15_50_R2.txt')
            natom=1090 ; r02= 350.0; alpha0= 41.9
      END IF 
      IF (name_data=='H16_50_R2') THEN     
            OPEN(unit=12,file='0H16_50_R2.txt')
            natom=473 ; r02= 350.0; alpha0= 42.1
      END IF       
      IF (name_data=='H17_50_R2') THEN     
            OPEN(unit=12,file='0H17_50_R2.txt')
            natom=295 ; r02= 350.0; alpha0= 42.1
      END IF 
      IF (name_data=='H18_50_R2') THEN     
            OPEN(unit=12,file='0H18_50_R2.txt')
            natom=889 ; r02= 350.0; alpha0= 43.9
      END IF       
      IF (name_data=='H19_50_R2') THEN     
            OPEN(unit=12,file='0H19_50_R2.txt')
            natom=1134 ; r02= 350.0; alpha0= 43.5
      END IF       
      IF (name_data=='H20_50_R2') THEN     
            OPEN(unit=12,file='0H20_50_R2.txt')
            natom=911 ; r02= 350.0; alpha0= 42.61
      END IF
      IF (name_data=='H21_50_R2') THEN     
            OPEN(unit=12,file='0H21_50_R2.txt')
            natom=1259 ; r02= 350.0; alpha0= 41.1
      END IF     
      IF (name_data=='H22_50_R2') THEN     
            OPEN(unit=12,file='0H22_50_R2.txt')
            natom=604 ; r02= 350.0; alpha0= 40.8
      END IF      
      IF (name_data=='H23_50_R2') THEN     
            OPEN(unit=12,file='0H23_50_R2.txt')
            natom=407 ; r02= 350.0; alpha0= 41.97
      END IF
      IF (name_data=='H24_50_R2') THEN     
            OPEN(unit=12,file='0H24_50_R2.txt')
            natom=4000 ; r02= 350.0; alpha0= 43.01
      END IF    
      IF (name_data=='H25_50_R2') THEN     
            OPEN(unit=12,file='0H25_50_R2.txt')
            natom=437 ; r02= 350.0; alpha0= 41.8
      END IF 
      IF (name_data=='H26_50_R2') THEN     
            OPEN(unit=12,file='0H26_50_R2.txt')
            natom=981 ; r02= 350.0; alpha0= 42.5
      END IF       
      IF (name_data=='H27_50_R2') THEN     
            OPEN(unit=12,file='0H27_50_R2.txt')
            natom=1452 ; r02= 350.0; alpha0= 42.71
      END IF 
      IF (name_data=='H28_50_R2') THEN     
            OPEN(unit=12,file='0H28_50_R2.txt')
            natom=1825 ; r02= 350.0; alpha0= 42.71
      END IF       
      IF (name_data=='H29_50_R2') THEN     
            OPEN(unit=12,file='0H29_50_R2.txt')
            natom=1453 ; r02= 350.0; alpha0= 42.82
      END IF       
      IF (name_data=='H30_50_R2') THEN     
            OPEN(unit=12,file='0H30_50_R2.txt')
            natom=792 ; r02= 350.0; alpha0= 42.5
      END IF  

!!!==================================================================  
!!! 20.05.2015: HCP725_R: 725 cells, random distribution
   IF (name_data=='HCP725_R') THEN     
      OPEN(unit=12,file='0HCP725_R.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   
!!!==================================================================  
!!! 16.05.2015: N902 voi 10 ensembles  
!!! N902: shell-core: pb=0.9, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N902_01') THEN     
      OPEN(unit=12,file='0N902_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N902_02') THEN     
      OPEN(unit=12,file='0N902_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N902_03') THEN     
      OPEN(unit=12,file='0N902_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N902_04') THEN     
      OPEN(unit=12,file='0N902_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N902_05') THEN     
      OPEN(unit=12,file='0N902_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N902_06') THEN     
      OPEN(unit=12,file='0N902_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N902_07') THEN     
      OPEN(unit=12,file='0N902_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N902_08') THEN     
      OPEN(unit=12,file='0N902_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N902_09') THEN     
      OPEN(unit=12,file='0N902_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N902_10') THEN     
      OPEN(unit=12,file='0N902_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 16.05.2015: N802 voi 10 ensembles  
!!! N802: shell-core: pb=0.8, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N802_01') THEN     
      OPEN(unit=12,file='0N802_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N802_02') THEN     
      OPEN(unit=12,file='0N802_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N802_03') THEN     
      OPEN(unit=12,file='0N802_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N802_04') THEN     
      OPEN(unit=12,file='0N802_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N802_05') THEN     
      OPEN(unit=12,file='0N802_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N802_06') THEN     
      OPEN(unit=12,file='0N802_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N802_07') THEN     
      OPEN(unit=12,file='0N802_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N802_08') THEN     
      OPEN(unit=12,file='0N802_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N802_09') THEN     
      OPEN(unit=12,file='0N802_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N802_10') THEN     
      OPEN(unit=12,file='0N802_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 16.05.2015: N702 voi 10 ensembles  
!!! N702: shell-core: pb=0.7, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N702_01') THEN     
      OPEN(unit=12,file='0N702_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N702_02') THEN     
      OPEN(unit=12,file='0N702_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N702_03') THEN     
      OPEN(unit=12,file='0N702_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N702_04') THEN     
      OPEN(unit=12,file='0N702_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N702_05') THEN     
      OPEN(unit=12,file='0N702_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N702_06') THEN     
      OPEN(unit=12,file='0N702_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N702_07') THEN     
      OPEN(unit=12,file='0N702_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N702_08') THEN     
      OPEN(unit=12,file='0N702_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N702_09') THEN     
      OPEN(unit=12,file='0N702_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N702_10') THEN     
      OPEN(unit=12,file='0N702_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 16.05.2015: N602 voi 10 ensembles  
!!! N602: shell-core: pb=0.6, Jaa=1, Jbb=3, Jab=1.8 :

   IF (name_data=='N602_01') THEN     
      OPEN(unit=12,file='0N602_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N602_02') THEN     
      OPEN(unit=12,file='0N602_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N602_03') THEN     
      OPEN(unit=12,file='0N602_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N602_04') THEN     
      OPEN(unit=12,file='0N602_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N602_05') THEN     
      OPEN(unit=12,file='0N602_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N602_06') THEN     
      OPEN(unit=12,file='0N602_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N602_07') THEN     
      OPEN(unit=12,file='0N602_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N602_08') THEN     
      OPEN(unit=12,file='0N602_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N602_09') THEN     
      OPEN(unit=12,file='0N602_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N602_10') THEN     
      OPEN(unit=12,file='0N602_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 16.05.2015: N502 voi 10 ensembles  
!!! N502: shell-core: pb=0.5, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N502_01') THEN     
      OPEN(unit=12,file='0N502_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N502_02') THEN     
      OPEN(unit=12,file='0N502_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N502_03') THEN     
      OPEN(unit=12,file='0N502_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N502_04') THEN     
      OPEN(unit=12,file='0N502_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N502_05') THEN     
      OPEN(unit=12,file='0N502_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N502_06') THEN     
      OPEN(unit=12,file='0N502_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N502_07') THEN     
      OPEN(unit=12,file='0N502_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N502_08') THEN     
      OPEN(unit=12,file='0N502_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N502_09') THEN     
      OPEN(unit=12,file='0N502_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N502_10') THEN     
      OPEN(unit=12,file='0N502_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 15.05.2015: HCP903 voi 10 ensembles  
!!! N903: partial mixing: pb=0.9, Jaa=1, Jbb=1, Jab=0.98 :

   IF (name_data=='N903_01') THEN     
      OPEN(unit=12,file='0N903_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N903_02') THEN     
      OPEN(unit=12,file='0N903_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N903_03') THEN     
      OPEN(unit=12,file='0N903_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N903_04') THEN     
      OPEN(unit=12,file='0N903_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N903_05') THEN     
      OPEN(unit=12,file='0N903_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N903_06') THEN     
      OPEN(unit=12,file='0N903_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N903_07') THEN     
      OPEN(unit=12,file='0N903_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N903_08') THEN     
      OPEN(unit=12,file='0N903_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N903_09') THEN     
      OPEN(unit=12,file='0N903_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N903_10') THEN     
      OPEN(unit=12,file='0N903_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 15.05.2015: HCP803 voi 10 ensembles  
!!! N803: partial mixing: pb=0.8, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N803_01') THEN     
      OPEN(unit=12,file='0N803_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N803_02') THEN     
      OPEN(unit=12,file='0N803_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N803_03') THEN     
      OPEN(unit=12,file='0N803_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N803_04') THEN     
      OPEN(unit=12,file='0N803_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N803_05') THEN     
      OPEN(unit=12,file='0N803_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N803_06') THEN     
      OPEN(unit=12,file='0N803_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N803_07') THEN     
      OPEN(unit=12,file='0N803_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N803_08') THEN     
      OPEN(unit=12,file='0N803_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N803_09') THEN     
      OPEN(unit=12,file='0N803_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N803_10') THEN     
      OPEN(unit=12,file='0N803_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 15.05.2015: HCP703 voi 10 ensembles  
!!! N703: partial mixing: pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N703_01') THEN     
      OPEN(unit=12,file='0N703_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N703_02') THEN     
      OPEN(unit=12,file='0N703_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N703_03') THEN     
      OPEN(unit=12,file='0N703_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N703_04') THEN     
      OPEN(unit=12,file='0N703_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N703_05') THEN     
      OPEN(unit=12,file='0N703_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N703_06') THEN     
      OPEN(unit=12,file='0N703_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N703_07') THEN     
      OPEN(unit=12,file='0N703_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N703_08') THEN     
      OPEN(unit=12,file='0N703_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N703_09') THEN     
      OPEN(unit=12,file='0N703_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N703_10') THEN     
      OPEN(unit=12,file='0N703_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!==================================================================  
!!! 15.05.2015: HCP603 voi 10 ensembles  
!!! N603: partial mixing: pb=0.6, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N603_01') THEN     
      OPEN(unit=12,file='0N603_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N603_02') THEN     
      OPEN(unit=12,file='0N603_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N603_03') THEN     
      OPEN(unit=12,file='0N603_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N603_04') THEN     
      OPEN(unit=12,file='0N603_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N603_05') THEN     
      OPEN(unit=12,file='0N603_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N603_06') THEN     
      OPEN(unit=12,file='0N603_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N603_07') THEN     
      OPEN(unit=12,file='0N603_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N603_08') THEN     
      OPEN(unit=12,file='0N603_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N603_09') THEN     
      OPEN(unit=12,file='0N603_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N603_10') THEN     
      OPEN(unit=12,file='0N603_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF    
   

!!!==================================================================  
!!! 15.05.2015: HCP503 voi 10 ensembles  
!!! N503: partial mixing: pb=0.5, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N503_01') THEN     
      OPEN(unit=12,file='0N503_01.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N503_02') THEN     
      OPEN(unit=12,file='0N503_02.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N503_03') THEN     
      OPEN(unit=12,file='0N503_03.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N503_04') THEN     
      OPEN(unit=12,file='0N503_04.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N503_05') THEN     
      OPEN(unit=12,file='0N503_05.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N503_06') THEN     
      OPEN(unit=12,file='0N503_06.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N503_07') THEN     
      OPEN(unit=12,file='0N503_07.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N503_08') THEN     
      OPEN(unit=12,file='0N503_08.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
   IF (name_data=='N503_09') THEN     
      OPEN(unit=12,file='0N503_09.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   IF (name_data=='N503_10') THEN     
      OPEN(unit=12,file='0N503_10.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF    
   

!!!==================================================================  
!!! 12.05.2015: Bo sung cau truc HCP voi N=725
   !!! N501: complete mixing: pb=0.5, Jaa=1, Jbb=1, Jab=0.2 :
   IF (name_data=='N501') THEN     
      OPEN(unit=12,file='0N501.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! N502: shell-core: pb=0.5, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N502') THEN     
      OPEN(unit=12,file='0N502.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF  
   !!! N503: partial mixing: pb=0.5, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N503') THEN     
      OPEN(unit=12,file='0N503.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

   !!! N601: complete mixing: pb=0.6, Jaa=1, Jbb=1, Jab=0.2 :
   IF (name_data=='N601') THEN     
      OPEN(unit=12,file='0N601.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! N602: shell-core: pb=0.6, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N602') THEN     
      OPEN(unit=12,file='0N602.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF  
   !!! N603: partial mixing: pb=0.6, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N603') THEN     
      OPEN(unit=12,file='0N603.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   
   !!! N701: complete mixing: pb=0.7, Jaa=1, Jbb=1, Jab=0.2 :
   IF (name_data=='N701') THEN     
      OPEN(unit=12,file='0N701.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! N702: shell-core: pb=0.7, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N702') THEN     
      OPEN(unit=12,file='0N702.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF  
   !!! N703: partial mixing: pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N703') THEN     
      OPEN(unit=12,file='0N703.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

   !!! N801: complete mixing: pb=0.8, Jaa=1, Jbb=1, Jab=0.2 :
   IF (name_data=='N801') THEN     
      OPEN(unit=12,file='0N801.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! N802: shell-core: pb=0.8, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N802') THEN     
      OPEN(unit=12,file='0N802.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF  
   !!! N803: partial mixing: pb=0.8, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N803') THEN     
      OPEN(unit=12,file='0N803.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

   !!! N901: complete mixing: pb=0.9, Jaa=1, Jbb=1, Jab=0.2 :
   IF (name_data=='N901') THEN     
      OPEN(unit=12,file='0N901.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! N902: shell-core: pb=0.9, Jaa=1, Jbb=3, Jab=1.8 :
   IF (name_data=='N902') THEN     
      OPEN(unit=12,file='0N902.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF  
   !!! N903: partial mixing: pb=0.9, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='N903') THEN     
      OPEN(unit=12,file='0N903.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 

!!!!!=================================================================================
 !!! == 07.05.2015: Human sau khi da go bo 30% beta cells =====
      IF (name_data=='H01_30') THEN     
            OPEN(unit=12,file='0H01_30.txt')
            natom=1745 ; r02= 350.0; alpha0= 44.2  
      END IF   
      IF (name_data=='H02_30') THEN     
            OPEN(unit=12,file='0H02_30.txt')
            natom=904 ; r02= 350.0; alpha0= 43.7 
      END IF    
      IF (name_data=='H04_30') THEN     
            OPEN(unit=12,file='0H04_30.txt')
            natom=502 ; r02= 350.0 ; alpha0= 41.9
      END IF  
      IF (name_data=='H05_30') THEN     
            OPEN(unit=12,file='0H05_30.txt')
            natom=950 ; r02= 360.0; alpha0= 40.5
      END IF
      IF (name_data=='H06_30') THEN     
            OPEN(unit=12,file='0H06_30.txt')
            natom=263 ; r02= 350.0; alpha0= 40.0 
      END IF      
      IF (name_data=='H07_30') THEN     
            OPEN(unit=12,file='0H07_30.txt')
            natom=2079 ; r02= 350.0; alpha0= 42.0 
      END IF 
      IF (name_data=='H08_30') THEN     
            OPEN(unit=12,file='0H08_30.txt')
            natom=1382 ; r02= 350.0; alpha0= 43.0 
      END IF   
      IF (name_data=='H09_30') THEN     
            OPEN(unit=12,file='0H09_30.txt')
            natom=674 ; r02= 350.0; alpha0= 42.0
      END IF
      IF (name_data=='H10_30') THEN     
            OPEN(unit=12,file='0H10_30.txt')
            natom=1797 ; r02= 350.0; alpha0= 42.21 
      END IF           
      IF (name_data=='H11_30') THEN     
            OPEN(unit=12,file='0H11_30.txt')
            natom=654 ; r02= 350.0; alpha0= 41.2
      END IF     
      IF (name_data=='H12_30') THEN     
            OPEN(unit=12,file='0H12_30.txt')
            natom=346 ; r02= 350.0; alpha0= 42.54
      END IF      
      IF (name_data=='H13_30') THEN     
            OPEN(unit=12,file='0H13_30.txt')
            natom=1565 ; r02= 350.0; alpha0= 42.6
      END IF
      IF (name_data=='H14_30') THEN     
            OPEN(unit=12,file='0H14_30.txt')
            natom=423 ; r02= 360.0; alpha0= 40.2
      END IF     
      IF (name_data=='H15_30') THEN     
            OPEN(unit=12,file='0H15_30.txt')
            natom=1394 ; r02= 350.0; alpha0= 41.9
      END IF 
      IF (name_data=='H16_30') THEN     
            OPEN(unit=12,file='0H16_30.txt')
            natom=623 ; r02= 350.0; alpha0= 42.1
      END IF       
      IF (name_data=='H17_30') THEN     
            OPEN(unit=12,file='0H17_30.txt')
            natom=388 ; r02= 350.0; alpha0= 42.1
      END IF 
      IF (name_data=='H18_30') THEN     
            OPEN(unit=12,file='0H18_30.txt')
            natom=1127 ; r02= 350.0; alpha0= 43.9
      END IF       
      IF (name_data=='H19_30') THEN     
            OPEN(unit=12,file='0H19_30.txt')
            natom=1396 ; r02= 350.0; alpha0= 43.5
      END IF       
      IF (name_data=='H20_30') THEN     
            OPEN(unit=12,file='0H20_30.txt')
            natom=1148 ; r02= 350.0; alpha0= 42.61
      END IF
      IF (name_data=='H21_30') THEN     
            OPEN(unit=12,file='0H21_30.txt')
            natom=1522 ; r02= 350.0; alpha0= 41.1
      END IF     
      IF (name_data=='H22_30') THEN     
            OPEN(unit=12,file='0H22_30.txt')
            natom=729 ; r02= 350.0; alpha0= 40.8
      END IF      
      IF (name_data=='H23_30') THEN     
            OPEN(unit=12,file='0H23_30.txt')
            natom=523 ; r02= 350.0; alpha0= 41.97
      END IF
      IF (name_data=='H25_30') THEN     
            OPEN(unit=12,file='0H25_30.txt')
            natom=508 ; r02= 350.0; alpha0= 41.8
      END IF 
      IF (name_data=='H26_30') THEN     
            OPEN(unit=12,file='0H26_30.txt')
            natom=1201 ; r02= 350.0; alpha0= 42.5
      END IF       
      IF (name_data=='H27_30') THEN     
            OPEN(unit=12,file='0H27_30.txt')
            natom=1669 ; r02= 350.0; alpha0= 42.71
      END IF 
      IF (name_data=='H28_30') THEN     
            OPEN(unit=12,file='0H28_30.txt')
            natom=2143 ; r02= 350.0; alpha0= 42.71
      END IF       
      IF (name_data=='H29_30') THEN     
            OPEN(unit=12,file='0H29_30.txt')
            natom=1685 ; r02= 350.0; alpha0= 42.82
      END IF       
      IF (name_data=='H30_30') THEN     
            OPEN(unit=12,file='0H30_30.txt')
            natom=941 ; r02= 350.0; alpha0= 42.5
      END IF   


!!!!!=================================================================================
 !!! == 07.05.2015: Human sau khi da go bo 90% beta cells =====
      IF (name_data=='H01_90') THEN     
            OPEN(unit=12,file='0H01_90.txt')
            natom=687 ; r02= 350.0; alpha0= 44.2  
      END IF   
      IF (name_data=='H02_90') THEN     
            OPEN(unit=12,file='0H02_90.txt')
            natom=262 ; r02= 350.0; alpha0= 43.7 
      END IF    
      IF (name_data=='H04_90') THEN     
            OPEN(unit=12,file='0H04_90.txt')
            natom=142 ; r02= 350.0 ; alpha0= 41.9
      END IF  
      IF (name_data=='H05_90') THEN     
            OPEN(unit=12,file='0H05_90.txt')
            natom=494 ; r02= 360.0; alpha0= 40.5
      END IF
      IF (name_data=='H06_90') THEN     
            OPEN(unit=12,file='0H06_90.txt')
            natom=109 ; r02= 350.0; alpha0= 40.0 
      END IF      
      IF (name_data=='H07_90') THEN     
            OPEN(unit=12,file='0H07_90.txt')
            natom=1293 ; r02= 350.0; alpha0= 42.0 
      END IF 
      IF (name_data=='H08_90') THEN     
            OPEN(unit=12,file='0H08_90.txt')
            natom=483 ; r02= 350.0; alpha0= 43.0 
      END IF   
      IF (name_data=='H09_90') THEN     
            OPEN(unit=12,file='0H09_90.txt')
            natom=189 ; r02= 350.0; alpha0= 42.0
      END IF
      IF (name_data=='H10_90') THEN     
            OPEN(unit=12,file='0H10_90.txt')
            natom=638 ; r02= 350.0; alpha0= 42.21 
      END IF           
      IF (name_data=='H11_90') THEN     
            OPEN(unit=12,file='0H11_90.txt')
            natom=324 ; r02= 350.0; alpha0= 41.2
      END IF     
      IF (name_data=='H12_90') THEN     
            OPEN(unit=12,file='0H12_90.txt')
            natom=122 ; r02= 350.0; alpha0= 42.54
      END IF      
      IF (name_data=='H13_90') THEN     
            OPEN(unit=12,file='0H13_90.txt')
            natom=743 ; r02= 350.0; alpha0= 42.6
      END IF
      IF (name_data=='H14_90') THEN     
            OPEN(unit=12,file='0H14_90.txt')
            natom=195 ; r02= 360.0; alpha0= 40.2
      END IF     
      IF (name_data=='H15_90') THEN     
            OPEN(unit=12,file='0H15_90.txt')
            natom=482 ; r02= 350.0; alpha0= 41.9
      END IF 
      IF (name_data=='H16_90') THEN     
            OPEN(unit=12,file='0H16_90.txt')
            natom=174 ; r02= 350.0; alpha0= 42.1
      END IF       
      IF (name_data=='H17_90') THEN     
            OPEN(unit=12,file='0H17_90.txt')
            natom=110 ; r02= 350.0; alpha0= 42.1
      END IF 
      IF (name_data=='H18_90') THEN     
            OPEN(unit=12,file='0H18_90.txt')
            natom=415 ; r02= 350.0; alpha0= 43.9
      END IF       
      IF (name_data=='H19_90') THEN     
            OPEN(unit=12,file='0H19_90.txt')
            natom=611 ; r02= 350.0; alpha0= 43.5
      END IF       
      IF (name_data=='H20_90') THEN     
            OPEN(unit=12,file='0H20_90.txt')
            natom=437 ; r02= 350.0; alpha0= 42.61
      END IF
      IF (name_data=='H21_90') THEN     
            OPEN(unit=12,file='0H21_90.txt')
            natom=732 ; r02= 350.0; alpha0= 41.1
      END IF     
      IF (name_data=='H22_90') THEN     
            OPEN(unit=12,file='0H22_90.txt')
            natom=355 ; r02= 350.0; alpha0= 40.8
      END IF      
      IF (name_data=='H23_90') THEN     
            OPEN(unit=12,file='0H23_90.txt')
            natom=177 ; r02= 350.0; alpha0= 41.97
      END IF
      IF (name_data=='H25_90') THEN     
            OPEN(unit=12,file='0H25_90.txt')
            natom=297 ; r02= 350.0; alpha0= 41.8
      END IF 
      IF (name_data=='H26_90') THEN     
            OPEN(unit=12,file='0H26_90.txt')
            natom=541 ; r02= 350.0; alpha0= 42.5
      END IF       
      IF (name_data=='H27_90') THEN     
            OPEN(unit=12,file='0H27_90.txt')
            natom=1019 ; r02= 350.0; alpha0= 42.71
      END IF 
      IF (name_data=='H28_90') THEN     
            OPEN(unit=12,file='0H28_90.txt')
            natom=1189 ; r02= 350.0; alpha0= 42.71
      END IF       
      IF (name_data=='H29_90') THEN     
            OPEN(unit=12,file='0H29_90.txt')
            natom=989 ; r02= 350.0; alpha0= 42.82
      END IF       
      IF (name_data=='H30_90') THEN     
            OPEN(unit=12,file='0H30_90.txt')
            natom=493 ; r02= 350.0; alpha0= 42.5
      END IF   

!!!!!=================================================================================
 !!! == 06.05.2015: Human sau khi da go bo 70% beta cells =====
      IF (name_data=='H01_70') THEN     
            OPEN(unit=12,file='0H01_70.txt')
            natom=1039 ; r02= 350.0; alpha0= 44.2  
      END IF   
      IF (name_data=='H02_70') THEN     
            OPEN(unit=12,file='0H02_70.txt')
            natom=476 ; r02= 350.0; alpha0= 43.7 
      END IF    
      IF (name_data=='H04_70') THEN     
            OPEN(unit=12,file='0H04_70.txt')
            natom=262 ; r02= 350.0 ; alpha0= 41.9
      END IF  
      IF (name_data=='H05_70') THEN     
            OPEN(unit=12,file='0H05_70.txt')
            natom=646 ; r02= 360.0; alpha0= 40.5
      END IF
      IF (name_data=='H06_70') THEN     
            OPEN(unit=12,file='0H06_70.txt')
            natom=161 ; r02= 350.0; alpha0= 40.0 
      END IF      
      IF (name_data=='H07_70') THEN     
            OPEN(unit=12,file='0H07_70.txt')
            natom=1556 ; r02= 350.0; alpha0= 42.0 
      END IF 
      IF (name_data=='H08_70') THEN     
            OPEN(unit=12,file='0H08_70.txt')
            natom=782 ; r02= 350.0; alpha0= 43.0 
      END IF   
      IF (name_data=='H09_70') THEN     
            OPEN(unit=12,file='0H09_70.txt')
            natom=351 ; r02= 350.0; alpha0= 42.0
      END IF
      IF (name_data=='H10_70') THEN     
            OPEN(unit=12,file='0H10_70.txt')
            natom=1024 ; r02= 350.0; alpha0= 42.21 
      END IF           
      IF (name_data=='H11_70') THEN     
            OPEN(unit=12,file='0H11_70.txt')
            natom=434 ; r02= 350.0; alpha0= 41.2
      END IF     
      IF (name_data=='H12_70') THEN     
            OPEN(unit=12,file='0H12_70.txt')
            natom=197 ; r02= 350.0; alpha0= 42.54
      END IF      
      IF (name_data=='H13_70') THEN     
            OPEN(unit=12,file='0H13_70.txt')
            natom=1017 ; r02= 350.0; alpha0= 42.6
      END IF
      IF (name_data=='H14_70') THEN     
            OPEN(unit=12,file='0H14_70.txt')
            natom=271 ; r02= 360.0; alpha0= 40.2
      END IF     
      IF (name_data=='H15_70') THEN     
            OPEN(unit=12,file='0H15_70.txt')
            natom=786 ; r02= 350.0; alpha0= 41.9
      END IF 
      IF (name_data=='H16_70') THEN     
            OPEN(unit=12,file='0H16_70.txt')
            natom=324 ; r02= 350.0; alpha0= 42.1
      END IF       
      IF (name_data=='H17_70') THEN     
            OPEN(unit=12,file='0H17_70.txt')
            natom=202 ; r02= 350.0; alpha0= 42.1
      END IF 
      IF (name_data=='H18_70') THEN     
            OPEN(unit=12,file='0H18_70.txt')
            natom=652 ; r02= 350.0; alpha0= 43.9
      END IF       
      IF (name_data=='H19_70') THEN     
            OPEN(unit=12,file='0H19_70.txt')
            natom=873 ; r02= 350.0; alpha0= 43.5
      END IF       
      IF (name_data=='H20_70') THEN     
            OPEN(unit=12,file='0H20_70.txt')
            natom=674 ; r02= 350.0; alpha0= 42.61
      END IF
      IF (name_data=='H21_70') THEN     
            OPEN(unit=12,file='0H21_70.txt')
            natom=996 ; r02= 350.0; alpha0= 41.1
      END IF     
      IF (name_data=='H22_70') THEN     
            OPEN(unit=12,file='0H22_70.txt')
            natom=479 ; r02= 350.0; alpha0= 40.8
      END IF      
      IF (name_data=='H23_70') THEN     
            OPEN(unit=12,file='0H23_70.txt')
            natom=292 ; r02= 350.0; alpha0= 41.97
      END IF
      IF (name_data=='H25_70') THEN     
            OPEN(unit=12,file='0H25_70.txt')
            natom=367 ; r02= 350.0; alpha0= 41.8
      END IF 
      IF (name_data=='H26_70') THEN     
            OPEN(unit=12,file='0H26_70.txt')
            natom=761 ; r02= 350.0; alpha0= 42.5
      END IF       
      IF (name_data=='H27_70') THEN     
            OPEN(unit=12,file='0H27_70.txt')
            natom=1235 ; r02= 350.0; alpha0= 42.71
      END IF 
      IF (name_data=='H28_70') THEN     
            OPEN(unit=12,file='0H28_70.txt')
            natom=1507 ; r02= 350.0; alpha0= 42.71
      END IF       
      IF (name_data=='H29_70') THEN     
            OPEN(unit=12,file='0H29_70.txt')
            natom=1221 ; r02= 350.0; alpha0= 42.82
      END IF       
      IF (name_data=='H30_70') THEN     
            OPEN(unit=12,file='0H30_70.txt')
            natom=643 ; r02= 350.0; alpha0= 42.5
      END IF   
      

!!!!!=================================================================================
 !!! == 2015.05.04: Human sau khi da go bo 50% beta cells =====
      IF (name_data=='H01_50') THEN     
            OPEN(unit=12,file='0H01_50.txt')
            natom=1392 ; r02= 350.0; alpha0= 44.2  
      END IF   
      IF (name_data=='H02_50') THEN     
            OPEN(unit=12,file='0H02_50.txt')
            natom=690 ; r02= 350.0; alpha0= 43.7 
      END IF    
      IF (name_data=='H03_50') THEN     
            OPEN(unit=12,file='0H03_50.txt')
            natom=1010 ; r02= 350.0; alpha0= 42.5 
      END IF  
      IF (name_data=='H04_50') THEN     
            OPEN(unit=12,file='0H04_50.txt')
            natom=382 ; r02= 350.0 ; alpha0= 41.9
      END IF  
      IF (name_data=='H05_50') THEN     
            OPEN(unit=12,file='0H05_50.txt')
            natom=798 ; r02= 360.0; alpha0= 40.5
      END IF
      IF (name_data=='H06_50') THEN     
            OPEN(unit=12,file='0H06_50.txt')
            natom=212 ; r02= 350.0; alpha0= 40.0 
      END IF      
      IF (name_data=='H07_50') THEN     
            OPEN(unit=12,file='0H07_50.txt')
            natom=1817 ; r02= 350.0; alpha0= 42.0 
      END IF 
      IF (name_data=='H08_50') THEN     
            OPEN(unit=12,file='0H08_50.txt')
            natom=1082 ; r02= 350.0; alpha0= 43.0 
      END IF   
      IF (name_data=='H09_50') THEN     
            OPEN(unit=12,file='0H09_50.txt')
            natom=512 ; r02= 350.0; alpha0= 42.0
      END IF
      IF (name_data=='H10_50') THEN     
            OPEN(unit=12,file='0H10_50.txt')
            natom=1410 ; r02= 350.0; alpha0= 42.21 
      END IF           
      IF (name_data=='H11_50') THEN     
            OPEN(unit=12,file='0H11_50.txt')
            natom=544 ; r02= 350.0; alpha0= 41.2
      END IF     
      IF (name_data=='H12_50') THEN     
            OPEN(unit=12,file='0H12_50.txt')
            natom=271 ; r02= 350.0; alpha0= 42.54
      END IF      
      IF (name_data=='H13_50') THEN     
            OPEN(unit=12,file='0H13_50.txt')
            natom=1291 ; r02= 350.0; alpha0= 42.6
      END IF
      IF (name_data=='H14_50') THEN     
            OPEN(unit=12,file='0H14_50.txt')
            natom=347 ; r02= 360.0; alpha0= 40.2
      END IF     
      IF (name_data=='H15_50') THEN     
            OPEN(unit=12,file='0H15_50.txt')
            natom=1090 ; r02= 350.0; alpha0= 41.9
      END IF 
      IF (name_data=='H16_50') THEN     
            OPEN(unit=12,file='0H16_50.txt')
            natom=473 ; r02= 350.0; alpha0= 42.1
      END IF       
      IF (name_data=='H17_50') THEN     
            OPEN(unit=12,file='0H17_50.txt')
            natom=295 ; r02= 350.0; alpha0= 42.1
      END IF 
      IF (name_data=='H18_50') THEN     
            OPEN(unit=12,file='0H18_50.txt')
            natom=889 ; r02= 350.0; alpha0= 43.9
      END IF       
      IF (name_data=='H19_50') THEN     
            OPEN(unit=12,file='0H19_50.txt')
            natom=1134 ; r02= 350.0; alpha0= 43.5
      END IF       
      IF (name_data=='H20_50') THEN     
            OPEN(unit=12,file='0H20_50.txt')
            natom=911 ; r02= 350.0; alpha0= 42.61
      END IF
      IF (name_data=='H21_50') THEN     
            OPEN(unit=12,file='0H21_50.txt')
            natom=1259 ; r02= 350.0; alpha0= 41.1
      END IF     
      IF (name_data=='H22_50') THEN     
            OPEN(unit=12,file='0H22_50.txt')
            natom=604 ; r02= 350.0; alpha0= 40.8
      END IF      
      IF (name_data=='H23_50') THEN     
            OPEN(unit=12,file='0H23_50.txt')
            natom=407 ; r02= 350.0; alpha0= 41.97
      END IF
      IF (name_data=='H24_50') THEN     
            OPEN(unit=12,file='0H24_50.txt')
            natom=4000 ; r02= 350.0; alpha0= 43.01
      END IF    
      IF (name_data=='H25_50') THEN     
            OPEN(unit=12,file='0H25_50.txt')
            natom=437 ; r02= 350.0; alpha0= 41.8
      END IF 
      IF (name_data=='H26_50') THEN     
            OPEN(unit=12,file='0H26_50.txt')
            natom=981 ; r02= 350.0; alpha0= 42.5
      END IF       
      IF (name_data=='H27_50') THEN     
            OPEN(unit=12,file='0H27_50.txt')
            natom=1452 ; r02= 350.0; alpha0= 42.71
      END IF 
      IF (name_data=='H28_50') THEN     
            OPEN(unit=12,file='0H28_50.txt')
            natom=1825 ; r02= 350.0; alpha0= 42.71
      END IF       
      IF (name_data=='H29_50') THEN     
            OPEN(unit=12,file='0H29_50.txt')
            natom=1453 ; r02= 350.0; alpha0= 42.82
      END IF       
      IF (name_data=='H30_50') THEN     
            OPEN(unit=12,file='0H30_50.txt')
            natom=792 ; r02= 350.0; alpha0= 42.5
      END IF   

!!!!!=================================================================================
   !!! ---------- For Human --------------
      IF (name_data=='H01') THEN     
            OPEN(unit=12,file='0H01.txt')
            natom=2273 ; r02= 350.0; alpha0= 44.2  
      END IF   
      IF (name_data=='H02') THEN     
            OPEN(unit=12,file='0H02.txt')
            natom=1225 ; r02= 350.0; alpha0= 43.7 
      END IF    
      IF (name_data=='H03') THEN     
            OPEN(unit=12,file='0H03.txt')
            natom=1114 ; r02= 350.0; alpha0= 42.5 
      END IF  
      IF (name_data=='H04') THEN     
            OPEN(unit=12,file='0H04.txt')
            natom=681 ; r02= 350.0 ; alpha0= 41.9
      END IF  
      IF (name_data=='H05') THEN     
            OPEN(unit=12,file='0H05.txt')
            natom=1178 ; r02= 360.0; alpha0= 40.5
      END IF
      IF (name_data=='H06') THEN     
            OPEN(unit=12,file='0H06.txt')
            natom=340 ; r02= 350.0; alpha0= 40.0 
      END IF      
      IF (name_data=='H07') THEN     
            OPEN(unit=12,file='0H07.txt')
            natom=2472 ; r02= 350.0; alpha0= 42.0 
      END IF 
      IF (name_data=='H08') THEN     
            OPEN(unit=12,file='0H08.txt')
            natom=1831 ; r02= 350.0; alpha0= 43.0 
      END IF   
      IF (name_data=='H09') THEN     
            OPEN(unit=12,file='0H09.txt')
            natom=916 ; r02= 350.0; alpha0= 42.0
      END IF
      IF (name_data=='H10') THEN     
            OPEN(unit=12,file='0H10.txt')
            natom=2376 ; r02= 350.0; alpha0= 42.21 
      END IF           
      IF (name_data=='H11') THEN     
            OPEN(unit=12,file='0H11.txt')
            natom=818 ; r02= 350.0; alpha0= 41.2
      END IF     
      IF (name_data=='H12') THEN     
            OPEN(unit=12,file='0H12.txt')
            natom=458 ; r02= 350.0; alpha0= 42.54
      END IF      
      IF (name_data=='H13') THEN     
            OPEN(unit=12,file='0H13.txt')
            natom=1976 ; r02= 350.0; alpha0= 42.6
      END IF
      IF (name_data=='H14') THEN     
            OPEN(unit=12,file='0H14.txt')
            natom=536 ; r02= 360.0; alpha0= 40.2
      END IF     
      IF (name_data=='H15') THEN     
            OPEN(unit=12,file='0H15.txt')
            natom=1850 ; r02= 350.0; alpha0= 41.9
      END IF 
      IF (name_data=='H16') THEN     
            OPEN(unit=12,file='0H16.txt')
            natom=847 ; r02= 350.0; alpha0= 42.1
      END IF       
      IF (name_data=='H17') THEN     
            OPEN(unit=12,file='0H17.txt')
            natom=526 ; r02= 350.0; alpha0= 42.1
      END IF 
      IF (name_data=='H18') THEN     
            OPEN(unit=12,file='0H18.txt')
            natom=1482 ; r02= 350.0; alpha0= 43.9
      END IF       
      IF (name_data=='H19') THEN     
            OPEN(unit=12,file='0H19.txt')
            natom=1788 ; r02= 350.0; alpha0= 43.5
      END IF       
      IF (name_data=='H20') THEN     
            OPEN(unit=12,file='0H20.txt')
            natom=1503 ; r02= 350.0; alpha0= 42.61
      END IF
      IF (name_data=='H21') THEN     
            OPEN(unit=12,file='0H21.txt')
            natom=1917 ; r02= 350.0; alpha0= 41.1
      END IF     
      IF (name_data=='H22') THEN     
            OPEN(unit=12,file='0H22.txt')
            natom=915 ; r02= 350.0; alpha0= 40.8
      END IF      
      IF (name_data=='H23') THEN     
            OPEN(unit=12,file='0H23.txt')
            natom=695 ; r02= 350.0; alpha0= 41.97
      END IF
      IF (name_data=='H24') THEN     
            OPEN(unit=12,file='0H24.txt')
            natom=5005 ; r02= 350.0; alpha0= 43.01
      END IF    
      IF (name_data=='H25') THEN     
            OPEN(unit=12,file='0H25.txt')
            natom=613 ; r02= 350.0; alpha0= 41.8
      END IF 
      IF (name_data=='H26') THEN     
            OPEN(unit=12,file='0H26.txt')
            natom=1531 ; r02= 350.0; alpha0= 42.5
      END IF       
      IF (name_data=='H27') THEN     
            OPEN(unit=12,file='0H27.txt')
            natom=1993 ; r02= 350.0; alpha0= 42.71
      END IF 
      IF (name_data=='H28') THEN     
            OPEN(unit=12,file='0H28.txt')
            natom=2620 ; r02= 350.0; alpha0= 42.71
      END IF       
      IF (name_data=='H29') THEN     
            OPEN(unit=12,file='0H29.txt')
            natom=2033 ; r02= 350.0; alpha0= 42.82
      END IF       
      IF (name_data=='H30') THEN     
            OPEN(unit=12,file='0H30.txt')
            natom=1165 ; r02= 350.0; alpha0= 42.5
      END IF   
      
      !!! ---------- For MOUSE --------------
      IF (name_data=='M01') THEN     
            OPEN(unit=12,file='0M01.txt')
            natom=893 ; r02= 450.0; alpha0= 41.3
      END IF 
      IF (name_data=='M02') THEN     
            OPEN(unit=12,file='0M02.txt')
            natom=2598 ; r02= 450.0; alpha0= 42.0
      END IF       
      IF (name_data=='M03') THEN     
            OPEN(unit=12,file='0M03.txt')
            natom=1870 ; r02= 450.0; alpha0= 41.75
      END IF      
      IF (name_data=='M04') THEN     
            OPEN(unit=12,file='0M04.txt')
            natom=1680 ; r02= 550.0; alpha0= 39.3
      END IF      
      IF (name_data=='M05') THEN     
            OPEN(unit=12,file='0M05.txt')
            natom=518 ; r02= 450.0; alpha0= 38.5
      END IF      
      IF (name_data=='M06') THEN     
            OPEN(unit=12,file='0M06.txt')
            natom=1895 ; r02= 450.0; alpha0= 40.6
      END IF 
      IF (name_data=='M07') THEN     
            OPEN(unit=12,file='0M07.txt')
            natom=571 ; r02= 450.0; alpha0= 40.87
      END IF 
      IF (name_data=='M08') THEN     
            OPEN(unit=12,file='0M08.txt')
            natom=1192 ; r02= 500.0; alpha0= 38.905
      END IF
      IF (name_data=='M09') THEN     
            OPEN(unit=12,file='0M09.txt')
            natom=1153 ; r02= 500.0; alpha0= 39.0
      END IF 
      IF (name_data=='M10') THEN     
            OPEN(unit=12,file='0M10.txt')
            natom=1063 ; r02= 450.0; alpha0= 42.02
      END IF 
     IF (name_data=='M11') THEN     
            OPEN(unit=12,file='0M11.txt')
            natom=2143 ; r02= 500.0; alpha0= 40.8
      END IF
      IF (name_data=='M12') THEN     
            OPEN(unit=12,file='0M12.txt')
            natom=1576 ; r02= 450.0; alpha0= 49.0
      END IF
      IF (name_data=='M13') THEN     
            OPEN(unit=12,file='0M13.txt')
            natom=2471 ; r02= 450.0; alpha0= 42.0
      END IF
      IF (name_data=='M14') THEN     
            OPEN(unit=12,file='0M14.txt')
            natom=1113 ; r02= 450.0; alpha0= 42.0
      END IF
      IF (name_data=='M15') THEN     
            OPEN(unit=12,file='0M15.txt')
            natom=1230 ; r02= 450.0; alpha0= 40.3
      END IF
      IF (name_data=='M16') THEN     
            OPEN(unit=12,file='0M16.txt')
            natom=3294 ; r02= 450.0; alpha0= 39.01
      END IF
      IF (name_data=='M17') THEN     
            OPEN(unit=12,file='0M17.txt')
            natom=1645 ; r02= 450.0; alpha0= 41.0
      END IF
      IF (name_data=='M18') THEN     
            OPEN(unit=12,file='0M18.txt')
            natom=4159 ; r02= 450.0; alpha0= 41.43
      END IF
      IF (name_data=='M19') THEN     
            OPEN(unit=12,file='0M19.txt')
            natom=2248 ; r02= 450.0; alpha0= 43.6
      END IF
      IF (name_data=='M20') THEN     
            OPEN(unit=12,file='0M20.txt')
            natom=1182 ; r02= 450.0; alpha0= 43.86
      END IF
      IF (name_data=='M21') THEN     
            OPEN(unit=12,file='0M21.txt')
            natom=4010 ; r02= 450.0; alpha0= 41.0
      END IF
      IF (name_data=='M22') THEN     
            OPEN(unit=12,file='0M22.txt')
            natom=3884 ; r02= 450.0; alpha0= 41.6
      END IF
      IF (name_data=='M23') THEN     
            OPEN(unit=12,file='0M23.txt')
            natom=3956 ; r02= 450.0; alpha0= 43.0
      END IF
      IF (name_data=='M24') THEN     
            OPEN(unit=12,file='0M24.txt')
            natom=2097 ; r02= 450.0; alpha0= 42.65
      END IF
      IF (name_data=='M25') THEN     
            OPEN(unit=12,file='0M25.txt')
            natom=3209 ; r02= 450.0; alpha0= 43.01
      END IF
      IF (name_data=='M26') THEN     
            OPEN(unit=12,file='0M26.txt')
            natom=9314 ; r02= 450.0; alpha0= 43.0
      END IF
      IF (name_data=='M27') THEN     
            OPEN(unit=12,file='0M27.txt')
            natom=1736 ; r02= 450.0; alpha0= 43.0
      END IF
      IF (name_data=='M28') THEN     
            OPEN(unit=12,file='0M28.txt')
            natom=1430 ; r02= 450.0; alpha0= 44.0
      END IF
      IF (name_data=='M29') THEN     
            OPEN(unit=12,file='0M29.txt')
            natom=2224 ; r02= 450.0; alpha0= 44.61
      END IF
      IF (name_data=='M30') THEN     
            OPEN(unit=12,file='0M30.txt')
            natom=1631 ; r02= 450.0; alpha0= 44.60
      END IF
      
!!! ------- For Human da luoc bo delta, chi con lai alpha va beta-----
      IF (name_data=='H41') THEN     
            OPEN(unit=12,file='0H41.txt')
            natom=467 ; r02= 400.0 ; alpha0= 42.0
      END IF 
      IF (name_data=='H42') THEN     
            OPEN(unit=12,file='0H42.txt')
            natom=1898 ; r02= 400.0 ; alpha0= 42.0
      END IF
      IF (name_data=='H43') THEN     
            OPEN(unit=12,file='0H43.txt')
            natom=2635 ; r02= 400.0 ; alpha0= 42.0
      END IF 
      IF (name_data=='H44') THEN     
            OPEN(unit=12,file='0H44.txt')
            natom=3194 ; r02= 400.0 ; alpha0= 42.0
      END IF 
      IF (name_data=='H45') THEN     
            OPEN(unit=12,file='0H45.txt')
            natom=1822 ; r02= 400.0 ; alpha0= 42.0
      END IF
      IF (name_data=='H46') THEN     
            OPEN(unit=12,file='0H46.txt')
            natom=2198 ; r02= 400.0 ; alpha0= 42.0
      END IF 
      
!!! ---------- For Human alpha, beta, delta--------------
      IF (name_data=='H51') THEN     
            OPEN(unit=12,file='0H51.txt')
            natom=588 ; r02= 400.0 ; alpha0= 42.0
      END IF   
      IF (name_data=='H52') THEN     
            OPEN(unit=12,file='0H52.txt')
            natom=2264 ; r02= 400.0 ; alpha0= 42.0
      END IF
      IF (name_data=='H53') THEN     
            OPEN(unit=12,file='0H53.txt')
            natom=3253 ; r02= 400.0 ; alpha0= 42.0
      END IF
      IF (name_data=='H54') THEN     
            OPEN(unit=12,file='0H54.txt')
            natom=3544 ; r02= 400.0 ; alpha0= 42.0
      END IF
      IF (name_data=='H55') THEN     
            OPEN(unit=12,file='0H55.txt')
            natom=2096 ; r02= 400.0 ; alpha0= 42.5
      END IF
      IF (name_data=='H56') THEN     
            OPEN(unit=12,file='0H56.txt')
            natom=2858 ; r02= 400.0 ; alpha0= 42.0
      END IF 
      
!!! 2015.04.27:       
!!!==================================================================
!!!   SC
   !!! SC complete sorting, pb=0.7, Jab=0.7 :
   IF (name_data=='SC701') THEN     
      OPEN(unit=12,file='0SC701.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! SC shell-core sorting, pb=0.7, Jab=0.98 :
   IF (name_data=='SC702') THEN     
      OPEN(unit=12,file='0SC702.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! SC partial mixing, pb=0.7, Jab=1.2 :
   IF (name_data=='SC703') THEN     
      OPEN(unit=12,file='0SC703.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! SC partial mixing, pb=0.8 :
   IF (name_data=='SC80') THEN     
      OPEN(unit=12,file='0SC80.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! SC partial mixing, pb=0.9 :
   IF (name_data=='SC90') THEN     
      OPEN(unit=12,file='0SC90.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
!!!==================================================================
!!!   Pb=50
   !!! Complete sorting, pb=0.5, Jaa=1, Jbb=1, Jab=0.2 :
   IF (name_data=='HCP501') THEN     
      OPEN(unit=12,file='0HCP501.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF  
   
   !!! HCP shell-core sorting, pb=0.5, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP502') THEN     
      OPEN(unit=12,file='0HCP502.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
  
   !!! HCP503 partial mixing, pb=0.5, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP503') THEN     
      OPEN(unit=12,file='0HCP503.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   
   !!! HCP504 complete mixing, pb=0.5, Jaa=1, Jbb=1, Jab=3:
   IF (name_data=='HCP504') THEN     
      OPEN(unit=12,file='0HCP504.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
!!!==================================================================
   !!! HCP sorting, pb=6 :
   IF (name_data=='HCP601') THEN     
      OPEN(unit=12,file='0HCP601.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP partial mixing, pb=0.6, Jab=0.98 :
   IF (name_data=='HCP602') THEN     
      OPEN(unit=12,file='0HCP602.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP shell-core sorting, pb=0.6, Jaa=1, Jbb=3, Jab=1.4 :
   IF (name_data=='HCP602S') THEN     
      OPEN(unit=12,file='0HCP602S.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP complete sorting, pb=0.7, Jab=0.7 :
   IF (name_data=='HCP701') THEN     
      OPEN(unit=12,file='0HCP701.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! HCP, pb=0.7, Jab=0.98 :
   IF (name_data=='HCP702') THEN     
      OPEN(unit=12,file='0HCP702.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP shell-core sorting, pb=0.7, Jaa=1, Jbb=3, Jab=1.81 :
   IF (name_data=='HCP702S') THEN     
      OPEN(unit=12,file='0HCP702S.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! HCP partial mixing, pb=0.7, Jab=1.2 :
   IF (name_data=='HCP703') THEN     
      OPEN(unit=12,file='0HCP703.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! HCP partial mixing, pb=0.8, Jab=0.98 :
   IF (name_data=='HCP802') THEN     
      OPEN(unit=12,file='0HCP802.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP shell-core, pb=0.8, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP802S') THEN     
      OPEN(unit=12,file='0HCP802S.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP803 partial mixing, pb=0.8, Jab=1.2 :
   IF (name_data=='HCP803') THEN     
      OPEN(unit=12,file='0HCP803.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
!!!Pb=0.93
   !!! HCP complete sorting, pb=0.93, Jab=0.7 :
   IF (name_data=='HCP901') THEN     
      OPEN(unit=12,file='0HCP901.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! HCP shell-core sorting, pb=0.93, Jab=0.91 :
   IF (name_data=='HCP902') THEN     
      OPEN(unit=12,file='0HCP902.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

   !!! HCP partial mixing, pb=0.93, Jab=1.2 :
   IF (name_data=='HCP903') THEN     
      OPEN(unit=12,file='0HCP903.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
!!!==================================================================   
!!! HCP random
   IF (name_data=='HCP') THEN     
      OPEN(unit=12,file='0HCP.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
!!!==================================================================  
!!! N=725
   !!! HCP725_70 partial mixing, N=725, pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP725_702') THEN     
      OPEN(unit=12,file='0HCP725_702.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   
   !!! HCP725_70S partial mixing, N=725, pb=0.7, Jaa=1, Jbb=3, Jab=1.5:
   IF (name_data=='HCP725_70S') THEN     
      OPEN(unit=12,file='0HCP725_70S.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   
   !!! HCP725_701 partial mixing, N=725, pb=0.7, Jaa=1, Jbb=1, Jab=0.5 :
   IF (name_data=='HCP725_701') THEN     
      OPEN(unit=12,file='0HCP725_701.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   
   !!! HCP725_80 partial mixing, N=725, pb=0.8, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP725_802') THEN     
      OPEN(unit=12,file='0HCP725_802.txt')
      natom=725 ; r02= 16.5 ; alpha0= 55.0
   END IF
!!!==================================================================  
!!! N=2493
   !!! HCP2493_60 partial mixing, N=2493, pb=0.6, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_60') THEN     
      OPEN(unit=12,file='0HCP2493_60.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF 
   !!! HCP2493_70 partial mixing, N=2493, pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_70') THEN     
      OPEN(unit=12,file='0HCP2493_70.txt')
      natom=2493 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
!!!==================================================================
!!! delta cells
   !!! HCP with delta cells, Pb=0.6, Pa=0.3, Pd=0.1, Jab=0.98, Jad=0.98, Jbd=0.98, N=1357 :
   IF (name_data=='HCP3') THEN     
      OPEN(unit=12,file='0HCP3.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP with delta cells, Pb=0.5, Pa=0.4, Pd=0.1, Jab=0.98, Jad=0.98, Jbd=0.98, N=1357 :
   IF (name_data=='HCP541') THEN     
      OPEN(unit=12,file='0HCP541.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF
   
   !!! HCP with delta cells, Pb=0.5, Pa=0.3, Pd=0.2, Jab=0.98, Jad=0.98, Jbd=0.98, N=1357 :
   IF (name_data=='HCP532') THEN     
      OPEN(unit=12,file='0HCP532.txt')
      natom=1357 ; r02= 16.5 ; alpha0= 55.0
   END IF

      
   END SUBROUTINE open_data
!!!======================================================================================
!!!======================================================================================
   SUBROUTINE read_data()
   IMPLICIT NONE
  
   number_a=0.; number_b=0. ; number_d=0.
   DO i=1,natom
      READ(12,*)Sabd(i),x(i),y(i),z(i)

      IF (Sabd(i)==11) THEN
         number_a=number_a+1.
      END IF

      IF (Sabd(i)==12) THEN
         number_b=number_b+1.
      END IF

      IF (Sabd(i)==13) THEN
         number_d=number_d+1.
      END IF
   
   END DO
      
   CLOSE(12)
       
   WRITE(*,*)number_a,number_b,number_d
   END SUBROUTINE read_data

!!!=======================================================================================
   SUBROUTINE read_data_random()
   IMPLICIT NONE
  
   DO i=1,natom
      READ(12,*)Sabd(i),x(i),y(i),z(i)  
   END DO
   CLOSE(12)
       
   number_b=real(nint(n_atom*pb))  ; number_a=n_atom-number_b ; number_d=0.
   Sabd(:)=0. ; na=0 ; nb=0
   
         DO WHILE (na<number_a)
            i=0

            DO WHILE(i==0)
                  CALL random_number(rdn_config)
                  i=int(rdn_config*n_atom)
            ENDDO

            IF (Sabd(i)==0.) THEN                      
                  Sabd(i)=11.
                  na=na+1
            END IF

      END DO

      DO i=1,natom

         IF (Sabd(i)==0.) THEN   
            Sabd(i)=12.
            nb=nb+1
         END IF

      ENDDO

   WRITE(*,*)'HCP random, na=',na,'nb=',nb
   
   !!WRITE(*,*)number_a,number_b,number_d
   
   END SUBROUTINE read_data_random

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! SUBROUTINE average phase phi_a, phi_b, phi_total: OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE average_phase()
   IMPLICIT NONE
     
   !!! Tinh phi_a
   IF (Zay >=0.) THEN      
      phi_a=acos(Zax/Za)
   ELSE
      phi_a=2.*pi-acos(Zax/Za)
   END IF

   !!! Tinh phi_b
   IF (Zby >=0.) THEN      
      phi_b=acos(Zbx/Zb)
   ELSE
      phi_b=2.*pi-acos(Zbx/Zb)
   END IF

   !!! Tinh phi_d
   IF (Zdy >=0.) THEN      
      phi_d=acos(Zdx/Zd)
   ELSE
      phi_d=2.*pi-acos(Zdx/Zd)
   END IF

   delta_ab=phi_a-phi_b
   delta_ad=phi_a-phi_d

   IF (delta_ab < 0.) THEN
      delta_ab=delta_ab+2.*pi
   END IF
   IF (delta_ad < 0.) THEN
      delta_ad=delta_ad+2.*pi
   END IF

   END SUBROUTINE average_phase

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   END PROGRAM

