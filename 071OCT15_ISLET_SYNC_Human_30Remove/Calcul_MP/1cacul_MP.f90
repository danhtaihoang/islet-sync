!!! ======================================================================================
!!! ======================================================================================
!!! 2015.04.14: Chuong trinh tong hop ket qua tu cach tinh MP
!!! 2015.05.01: Sua lai cho tinh Ra, Rb, Rd, R 
!!! 11.05.2015: Sua lai dung standard error, khong phai la standard deviation 
!!! ======================================================================================
!!! ======================================================================================      
   PROGRAM calcul_MP
   IMPLICIT NONE
      
   INTEGER (KIND=8), PARAMETER :: nF=28 , nKab=17
   INTEGER (KIND=8) :: i,j
   REAL    (KIND=8) :: Kad,Rmin,Rmax,WS

   REAL    (KIND=8),DIMENSION(nF,nKab)  :: Kab,Ra,Rb,Rd,R
   REAL    (KIND=8),DIMENSION(nKab) :: Ra_av,Rb_av,Rd_av,R_av,Ra_dev,Rb_dev,Rd_dev,R_dev
      
    
!!!=======================================================================================
!!!=======================================================================================
!!! Doc gia tri tu cac file       
   DO i=1,nF
      IF (i==1) THEN
         OPEN(unit=12,file='average_H01.txt')
      END IF

      IF (i==2) THEN
         OPEN(unit=12,file='average_H02.txt')
      END IF

      IF (i==3) THEN
         OPEN(unit=12,file='average_H30.txt')
      END IF
      
      IF (i==4) THEN
         OPEN(unit=12,file='average_H04.txt')
      END IF
      
      IF (i==5) THEN
         OPEN(unit=12,file='average_H05.txt')
      END IF
      
      IF (i==6) THEN
         OPEN(unit=12,file='average_H06.txt')
      END IF
      
      IF (i==7) THEN
         OPEN(unit=12,file='average_H07.txt')
      END IF
      
      IF (i==8) THEN
         OPEN(unit=12,file='average_H08.txt')
      END IF
      
      IF (i==9) THEN
         OPEN(unit=12,file='average_H09.txt')
      END IF
      
      IF (i==10) THEN
         OPEN(unit=12,file='average_H10.txt')
      END IF
      
      IF (i==11) THEN
         OPEN(unit=12,file='average_H11.txt')
      END IF

      IF (i==12) THEN
         OPEN(unit=12,file='average_H12.txt')
      END IF

      IF (i==13) THEN
         OPEN(unit=12,file='average_H13.txt')
      END IF
      
      IF (i==14) THEN
         OPEN(unit=12,file='average_H14.txt')
      END IF
      
      IF (i==15) THEN
         OPEN(unit=12,file='average_H15.txt')
      END IF
      
      IF (i==16) THEN
         OPEN(unit=12,file='average_H16.txt')
      END IF
      
      IF (i==17) THEN
         OPEN(unit=12,file='average_H17.txt')
      END IF
      
      IF (i==18) THEN
         OPEN(unit=12,file='average_H18.txt')
      END IF
      
      IF (i==19) THEN
         OPEN(unit=12,file='average_H19.txt')
      END IF     
      
      IF (i==20) THEN
         OPEN(unit=12,file='average_H20.txt')
      END IF         
      
      IF (i==21) THEN
         OPEN(unit=12,file='average_H21.txt')
      END IF

      IF (i==22) THEN
         OPEN(unit=12,file='average_H22.txt')
      END IF

      IF (i==23) THEN
         OPEN(unit=12,file='average_H23.txt')
      END IF
      
      IF (i==24) THEN
         OPEN(unit=12,file='average_H29.txt')
      END IF
      
      IF (i==25) THEN
         OPEN(unit=12,file='average_H25.txt')
      END IF
      
      IF (i==26) THEN
         OPEN(unit=12,file='average_H26.txt')
      END IF
      
      IF (i==27) THEN
         OPEN(unit=12,file='average_H27.txt')
      END IF
      
      IF (i==28) THEN
         OPEN(unit=12,file='average_H28.txt')
      END IF
      
      
      DO j=1,nKab
         READ(12,*)Kad,Kab(i,j),Ra(i,j),Rb(i,j),Rd(i,j),R(i,j),Rmin,Rmax,WS
      END DO     
                  
   END DO 
                 
!!!=======================================================================================
!!! Tinh gia tri trung binh

   OPEN(unit=21,file='average_MP.dat')  
   OPEN(unit=22,file='Ra_MP.dat')     
   OPEN(unit=23,file='Rd_MP.dat') 
       
   Ra_av(:)=0. ; Rb_av(:)=0. ; Rd_av(:)=0. ; R_av(:)=0.
 
   DO j=1,nKab
      DO i=1,nF
         Ra_av(j)=Ra_av(j)+Ra(i,j)
         Rb_av(j)=Rb_av(j)+Rb(i,j)
         Rd_av(j)=Rd_av(j)+Rd(i,j)
         R_av(j)=R_av(j)+R(i,j)
                
      END DO
      
      Ra_av(j)=Ra_av(j)/real(nF)
      Rb_av(j)=Rb_av(j)/real(nF)
      Rd_av(j)=Rd_av(j)/real(nF)
      R_av(j)=R_av(j)/real(nF)
      
   END DO

!!! Tinh deviation        
   Ra_dev(:)=0. ; Rb_dev(:)=0. ; Rd_dev(:)=0. ; R_dev(:)=0.
   DO j=1,nKab
      DO i=1,nF
         Ra_dev(j)=Ra_dev(j)+(Ra(i,j)-Ra_av(j))**2.
         Rb_dev(j)=Rb_dev(j)+(Rb(i,j)-Rb_av(j))**2.
         Rd_dev(j)=Rd_dev(j)+(Rd(i,j)-Rd_av(j))**2.
         R_dev(j)=R_dev(j)+(R(i,j)-R_av(j))**2.
      END DO
      !Ra_dev(j)=sqrt(Ra_dev(j)/nF)
      !Rb_dev(j)=sqrt(Rb_dev(j)/nF)
      !Rd_dev(j)=sqrt(Rd_dev(j)/nF)
      !R_dev(j)=sqrt(R_dev(j)/nF)
      
      !!! Standard error:
      Ra_dev(j)=sqrt(Ra_dev(j))/nF
      Rb_dev(j)=sqrt(Rb_dev(j))/nF
      Rd_dev(j)=sqrt(Rd_dev(j))/nF
      R_dev(j)=sqrt(R_dev(j))/nF
      
      
   END DO        

!!! WRITE         
   DO j=1,nKab
      WRITE(21,*)Kad,Kab(1,j),Ra_av(j),Rb_av(j),Rd_av(j),R_av(j),R_av(j)-R_dev(j),R_av(j)+R_dev(j)
      WRITE(22,*)Kad,Kab(1,j),Ra_av(j),Ra_av(j)-Ra_dev(j),Ra_av(j)+Ra_dev(j),Rb_av(j) &
      ,Rb_av(j)-Rb_dev(j),Rb_av(j)+Rb_dev(j)
      WRITE(23,*)Kad,Kab(1,j),Rd_av(j),Rd_av(j)-Rd_dev(j),Rd_av(j)+Rd_dev(j)
      
   END DO
         
         
   CLOSE(12)
   CLOSE(21) 
   CLOSE(22)
   CLOSE(23)

!!! ================================================================================================
      
   END PROGRAM calcul_MP


