   PROGRAM mean_entropy
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::n=5

   INTEGER (KIND=8):: i
   REAL    (KIND=8):: Ra_av,Rb_av,Rd_av,Kad,Kab
   REAL    (KIND=8):: Ra_dev,Rb_dev,Rd_dev
   REAL    (KIND=8),DIMENSION(n):: Ra,Rb,Rd

!!!=======================================================================================
!!!======================================================================================= 
   OPEN(unit=11,file='entropy.dat')
   OPEN(unit=12,file='entropy_av.dat')
   OPEN(unit=13,file='S_av.dat')

!!! Tinh gia tri trung binh
   Ra_av=0. ; Rb_av=0. ; Rd_av=0.

   DO i=1,n
      READ(11,*)Kad,Kab,Ra(i),Rb(i),Rd(i)
      Ra_av=Ra_av+Ra(i)
      Rb_av=Rb_av+Rb(i)
      Rd_av=Rd_av+Rd(i)
   END DO
 
   Ra_av=Ra_av/n ; Rb_av=Rb_av/n ; Rd_av=Rd_av/n

!!! Tinh standard deviation
   Ra_dev=0.   
   DO i=1,n
      Ra_dev=Ra_dev+(Ra(i)-Ra_av)**2.
   END DO
   Ra_dev=sqrt(Ra_dev/n)

   Rb_dev=0.   
   DO i=1,n
      Rb_dev=Rb_dev+(Rb(i)-Rb_av)**2.
   END DO
   Rb_dev=sqrt(Rb_dev/n)

   Rd_dev=0.   
   DO i=1,n
      Rd_dev=Rd_dev+(Rd(i)-Rd_av)**2.
   END DO
   Rd_dev=sqrt(Rd_dev/n)

   WRITE(12,*)Kad,Kab,Ra_av,Rb_av,Rd_av
   WRITE(13,*)Ra_av,Ra_dev,Rb_av,Rb_dev,Rd_av,Rd_dev

   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
    
   END PROGRAM mean_entropy

