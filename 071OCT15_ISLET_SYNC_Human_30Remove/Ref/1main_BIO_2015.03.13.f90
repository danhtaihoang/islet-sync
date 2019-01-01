!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  HOANG Danh Tai - Asia Pacific Center for Theoretical Physics, 
!  Hogil Kim Memorial Building #501 POSTECH,
!  San 31, Hyoja-dong, Namgum, Pohang, Gyeongbuk 790-784, Korea.
!  E-mail: hoangdanhtai@gmail.com    Personal site: hoangdanhtai.com 
!-----------------------------------------------------------------------------------------!
!!! 2014.06.04: DYNAMIC IN PANCREATIC ISLET
!!! alpha: 11, beta: 12, delta: 13
!!! 2015.02.17: Bo sung structure, chuyen HCP503 va HCP504, ...
!!! 2015.03.13: Sua lai histogram n=100
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%! 
   PROGRAM main_BIO
   IMPLICIT NONE

   CHARACTER (LEN=150):: LOAD_PHI,NAME_DATA
   CHARACTER (LEN=3)  :: SAt
   REAL      (KIND=8),PARAMETER :: nul=0.,delt=0.01,aaa=1.0,r0=4.5,Kaa=1.,Kbb=1.,Kdd=1.

   CHARACTER (256)    :: Ligne20,Ligne21,Ligne22,Ligne23,Ligne40

   INTEGER (KIND=8) :: i,j,i_n,nn_total
   INTEGER (KIND=8):: i_loop,i_loop1,i_av,n_equi1,n_equi2,n_average,i_times,n_times
   INTEGER (KIND=8) :: time_2,time_3,time_5,time_6,time_7,natom,na,nb
   
   REAL (KIND=8) :: n_atom,number_a,number_b,number_d,rdn_s,Kab,Kba,Kad,Kda,Kbd,Kdb,r
   REAL (KIND=8) :: pi,delt1,Zax,Zay,Za,Zbx,Zby,Zb,Za_av,Zb_av,Zdx,Zdy,Zd,Zd_av,Sabd_load
   REAL (KIND=8) :: Zx,Zy,Z_total,Z_total_av,run_time,S_tmp,S_load,WS,WS_av
   REAL (KIND=8) :: phi_a,phi_b,phi_d,delta_ab,delta_ad,Pb,rdn_config

   INTEGER (KIND=8),DIMENSION(:),ALLOCATABLE :: nn
   INTEGER (KIND=8),DIMENSION(:,:),ALLOCATABLE :: name_in
   REAL (KIND=8),DIMENSION(:),ALLOCATABLE :: S,delS,S1,x,y,z,Sabd
   REAL (KIND=8),DIMENSION(:,:),ALLOCATABLE :: Ks,Ks1

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


   IF (NAME_DATA=='HCP') THEN
      CALL read_data_random()
   ELSE
      CALL read_data()
   END IF
   !!WRITE(*,*)number_a,number_b
   

   delt1=delt+(1./2.)*delt**2. ! +(1./6.)*delt**3.+(1./24.)*delt**4. !! 2nd order
   
   pi=acos(-1.)

!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!! ====== MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM === MAIN PROGRAM ======
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
   CALL nearest_neighbor()
   
   CALL Ks_value()   

   DO i=1,natom
      DO j=1,natom
      Ks1(i,j)=Ks(i,j)/nn(i)
      END DO
   END DO
   
   CALL equi_lattice1()
   WRITE(*,*)'Finish Equilibre'

   CALL average_thermal()

   CALL save_phi_equi()

   CALL computation_time()

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
   READ(11, '(A30,(A10))')   tamp, NAME_DATA
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
         IF (Sabd(i)<11.5) THEN
            Zax=Zax+cos(S_tmp) 
            Zay=Zay+sin(S_tmp)
         ELSE
            IF (Sabd(i)<12.5) THEN
            Zbx=Zbx+cos(S_tmp) 
            Zby=Zby+sin(S_tmp)
            ELSE
            Zdx=Zdx+cos(S_tmp) 
            Zdy=Zdy+sin(S_tmp)
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

   WRITE(Ligne22,*)Za_av,Zb_av,Zd_av,Z_total_av,WS_av
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
!!! SUBROUTINE nearest_neighbor() : find nearest neighbor j of every i and number of nn
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE nearest_neighbor()
   IMPLICIT NONE

   OPEN(unit=15,file='nn_i_cell.dat') 
      
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

   CLOSE(15)

   END SUBROUTINE nearest_neighbor
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
!!!  SUBROUTINE Open data and Read_data() Read experimental data : OK
!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   SUBROUTINE open_data()
   IMPLICIT NONE

!!!==================================================================
   !!! TEST : only 1 alpha + 1 beta :
   IF (name_data=='TEST') THEN     
      OPEN(unit=12,file='0TEST.txt')
      natom=2
   END IF   
   !!! TEST1 : only 1 alpha + 12 beta around :
   IF (name_data=='TEST1') THEN     
      OPEN(unit=12,file='0TEST1.txt')
      natom=13
   END IF
   !!! TEST11 : only 1 alpha on the peripery+ 12 beta :
   IF (name_data=='TEST11') THEN     
      OPEN(unit=12,file='0TEST11.txt')
      natom=13
   END IF
   !!! TEST2 : only 2 alpha on the peripery:
   IF (name_data=='TEST2') THEN     
      OPEN(unit=12,file='0TEST2.txt')
      natom=13
   END IF
   
   !!! TEST12:
   IF (name_data=='TEST12') THEN     
      OPEN(unit=12,file='0TEST12.txt')
      natom=12
   END IF
   
!!!==================================================================
!!!   SC
   !!! SC complete sorting, pb=0.7, Jab=0.7 :
   IF (name_data=='SC701') THEN     
      OPEN(unit=12,file='0SC701.txt')
      natom=1357
   END IF

   !!! SC shell-core sorting, pb=0.7, Jab=0.98 :
   IF (name_data=='SC702') THEN     
      OPEN(unit=12,file='0SC702.txt')
      natom=1357
   END IF

   !!! SC partial mixing, pb=0.7, Jab=1.2 :
   IF (name_data=='SC703') THEN     
      OPEN(unit=12,file='0SC703.txt')
      natom=1357
   END IF

   !!! SC partial mixing, pb=0.8 :
   IF (name_data=='SC80') THEN     
      OPEN(unit=12,file='0SC80.txt')
      natom=1357
   END IF

   !!! SC partial mixing, pb=0.9 :
   IF (name_data=='SC90') THEN     
      OPEN(unit=12,file='0SC90.txt')
      natom=1357
   END IF
!!!==================================================================
!!!   Pb=50
   !!! Complete sorting, pb=0.5, Jaa=1, Jbb=1, Jab=0.2 :
   IF (name_data=='HCP501') THEN     
      OPEN(unit=12,file='0HCP501.txt')
      natom=1357
   END IF  
   
   !!! HCP shell-core sorting, pb=0.5, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP502') THEN     
      OPEN(unit=12,file='0HCP502.txt')
      natom=1357
   END IF
  
   !!! HCP503 partial mixing, pb=0.5, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP503') THEN     
      OPEN(unit=12,file='0HCP503.txt')
      natom=1357
   END IF 
   
   !!! HCP504 complete mixing, pb=0.5, Jaa=1, Jbb=1, Jab=3:
   IF (name_data=='HCP504') THEN     
      OPEN(unit=12,file='0HCP504.txt')
      natom=1357
   END IF
!!!==================================================================
   !!! HCP sorting, pb=6 :
   IF (name_data=='HCP601') THEN     
      OPEN(unit=12,file='0HCP601.txt')
      natom=1357
   END IF
   
   !!! HCP partial mixing, pb=0.6, Jab=0.98 :
   IF (name_data=='HCP602') THEN     
      OPEN(unit=12,file='0HCP602.txt')
      natom=1357
   END IF
   
   !!! HCP shell-core sorting, pb=0.6, Jaa=1, Jbb=3, Jab=1.4 :
   IF (name_data=='HCP602S') THEN     
      OPEN(unit=12,file='0HCP602S.txt')
      natom=1357
   END IF
   
   !!! HCP complete sorting, pb=0.7, Jab=0.7 :
   IF (name_data=='HCP701') THEN     
      OPEN(unit=12,file='0HCP701.txt')
      natom=1357
   END IF

   !!! HCP, pb=0.7, Jab=0.98 :
   IF (name_data=='HCP702') THEN     
      OPEN(unit=12,file='0HCP702.txt')
      natom=1357
   END IF
   
   !!! HCP shell-core sorting, pb=0.7, Jaa=1, Jbb=3, Jab=1.81 :
   IF (name_data=='HCP702S') THEN     
      OPEN(unit=12,file='0HCP702S.txt')
      natom=1357
   END IF

   !!! HCP partial mixing, pb=0.7, Jab=1.2 :
   IF (name_data=='HCP703') THEN     
      OPEN(unit=12,file='0HCP703.txt')
      natom=1357
   END IF

   !!! HCP partial mixing, pb=0.8, Jab=0.98 :
   IF (name_data=='HCP802') THEN     
      OPEN(unit=12,file='0HCP802.txt')
      natom=1357
   END IF
   
   !!! HCP shell-core, pb=0.8, Jaa=1, Jbb=3, Jab=1.5 :
   IF (name_data=='HCP802S') THEN     
      OPEN(unit=12,file='0HCP802S.txt')
      natom=1357
   END IF
   
   !!! HCP803 partial mixing, pb=0.8, Jab=1.2 :
   IF (name_data=='HCP803') THEN     
      OPEN(unit=12,file='0HCP803.txt')
      natom=1357
   END IF
   
!!!Pb=0.93
   !!! HCP complete sorting, pb=0.93, Jab=0.7 :
   IF (name_data=='HCP901') THEN     
      OPEN(unit=12,file='0HCP901.txt')
      natom=1357
   END IF

   !!! HCP shell-core sorting, pb=0.93, Jab=0.91 :
   IF (name_data=='HCP902') THEN     
      OPEN(unit=12,file='0HCP902.txt')
      natom=1357
   END IF

   !!! HCP partial mixing, pb=0.93, Jab=1.2 :
   IF (name_data=='HCP903') THEN     
      OPEN(unit=12,file='0HCP903.txt')
      natom=1357
   END IF
   
!!!==================================================================   
!!! HCP random
   IF (name_data=='HCP') THEN     
      OPEN(unit=12,file='0HCP.txt')
      natom=1357
   END IF
   
!!!==================================================================  
!!! N=725
   !!! HCP725_70 partial mixing, N=725, pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP725_70') THEN     
      OPEN(unit=12,file='0HCP725_70.txt')
      natom=725
   END IF 
   !!! HCP725_80 partial mixing, N=725, pb=0.8, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP725_80') THEN     
      OPEN(unit=12,file='0HCP725_80.txt')
      natom=725
   END IF
!!!==================================================================  
!!! N=2493
   !!! HCP2493_60 partial mixing, N=2493, pb=0.6, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_60') THEN     
      OPEN(unit=12,file='0HCP2493_60.txt')
      natom=2493
   END IF 
   !!! HCP2493_70 partial mixing, N=2493, pb=0.7, Jaa=1, Jbb=1, Jab=0.98 :
   IF (name_data=='HCP2493_70') THEN     
      OPEN(unit=12,file='0HCP2493_70.txt')
      natom=2493
   END IF
   
!!!==================================================================
!!! delta cells
   !!! HCP with delta cells, Pb=0.6, Pa=0.3, Pd=0.1, Jab=0.98, Jad=0.98, Jbd=0.98, N=1357 :
   IF (name_data=='HCP3') THEN     
      OPEN(unit=12,file='0HCP3.txt')
      natom=1357
   END IF
   
   !!! HCP with delta cells, Pb=0.5, Pa=0.4, Pd=0.1, Jab=0.98, Jad=0.98, Jbd=0.98, N=1357 :
   IF (name_data=='HCP541') THEN     
      OPEN(unit=12,file='0HCP541.txt')
      natom=1357
   END IF
   
   !!! HCP with delta cells, Pb=0.5, Pa=0.3, Pd=0.2, Jab=0.98, Jad=0.98, Jbd=0.98, N=1357 :
   IF (name_data=='HCP532') THEN     
      OPEN(unit=12,file='0HCP532.txt')
      natom=1357
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

