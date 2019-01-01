   PROGRAM entropy_value
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::n=100

   INTEGER (KIND=8):: i
   REAL    (KIND=8):: S1,S2,S3,P1_av,P2_av,P3_av,Smax,Kab,Kad,n1
   REAL    (KIND=8),DIMENSION(n):: P1,P2,P3,phi

!!!=======================================================================================
!!!======================================================================================= 
   OPEN(unit=12,file='histogram.dat')
   OPEN(unit=13,file='entropy.dat')

   n1=real(n)

   CALL read_input_parameter()
   CALL entropy()

   CONTAINS
!!!=======================================================================================
!!!======================================================================================= 
   SUBROUTINE read_input_parameter()
   IMPLICIT NONE

   CHARACTER (LEN=150) :: tamp
   OPEN(11,file='1parameter.in')   
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A30,(F12.6))') tamp, Kab
   READ(11, '(A30,(F12.6))') tamp, Kad
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp
   READ(11, '(A50)')         tamp

   CLOSE(11) 

   END SUBROUTINE read_input_parameter

!!!=======================================================================================
!!!======================================================================================= 
   SUBROUTINE entropy()
   IMPLICIT NONE
 
!!! Doc gia tri P
   P1_av=0. ; P2_av=0. ; P3_av=0.
   DO i=1,n
      READ(12,*)phi(i),P1(i),P2(i),P3(i)
      P1_av=P1_av+P1(i)
      P2_av=P2_av+P2(i)
      P3_av=P3_av+P3(i)
   END DO

!!! Smax

  ! Smax=log(100.)/log(2.)

!!! Tinh entropy S
   S1=0. ; S2=0. ; S3=0.  
   
   DO i=0,n  
     
      !IF (P1(i)/=0.) THEN
      IF (P1(i)>0.001) THEN         
         S1=S1-P1(i)*log(P1(i))
      END IF
      
      !IF (P2(i)/=0.) THEN   
      IF (P2(i)>0.001) THEN      
         S2=S2-P2(i)*log(P2(i))
      END IF
      
      !IF (P3(i)/=0.) THEN  
      IF (P3(i)>0.001) THEN       
         S3=S3-P3(i)*log(P3(i))
      END IF

   END DO
      
   !S1=S1/log(2.)/Smax ; S2=S2/log(2.)/Smax ; S3=S3/log(2.)/Smax
   S1=S1/log(n1) ; S2=S2/log(n1) ; S3=S3/log(n1)

   WRITE(13,*)Kad,Kab,S1,S2,S3
   WRITE(*,*)S1,S2,S3,P1_av,P2_av,P3_av,Smax

   CLOSE(12)
   CLOSE(13)

   END SUBROUTINE entropy
    
   END PROGRAM entropy_value

