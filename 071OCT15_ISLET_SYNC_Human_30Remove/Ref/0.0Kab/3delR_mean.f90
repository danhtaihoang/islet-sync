   PROGRAM delR_mean
   IMPLICIT NONE

   INTEGER (KIND=8),PARAMETER::n=5

   INTEGER (KIND=8):: i
   REAL    (KIND=8):: Ra_av,Rb_av,Rd_av,R_av,Kad,Kab
   REAL    (KIND=8):: Ra_dev,Rb_dev,Rd_dev,R_dev
   REAL    (KIND=8),DIMENSION(n):: Ra,Rb,Rd,R

!!!=======================================================================================
!!!======================================================================================= 
   OPEN(unit=11,file='delR.dat')
   OPEN(unit=12,file='delR_av.dat')
   OPEN(unit=13,file='delRa_av.dat')

!!! Tinh gia tri trung binh
   Ra_av=0. ; Rb_av=0. ; Rd_av=0. ; R_av=0.

   DO i=1,n
      READ(11,*)Kad,Kab,Ra(i),Rb(i),Rd(i),R(i)
      Ra_av=Ra_av+Ra(i)
      Rb_av=Rb_av+Rb(i)
      Rd_av=Rd_av+Rd(i)
      R_av=R_av+R(i)
   END DO
 
   Ra_av=Ra_av/n ; Rb_av=Rb_av/n ; Rd_av=Rd_av/n ; R_av=R_av/n

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

   R_dev=0.   
   DO i=1,n
      R_dev=R_dev+(R(i)-R_av)**2.
   END DO
   R_dev=sqrt(R_dev/n)

   WRITE(12,*)Kad,Kab,Ra_av,Rb_av,Rd_av,R_av
   WRITE(13,*)Ra_av,Ra_dev,Rb_av,Rb_dev,Rd_av,Rd_dev,R_av,R_dev
   
   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
    
   END PROGRAM delR_mean

