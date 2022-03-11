      program binodal
      implicit real*8 (A-H, O-Z)
      logical*1 text

      COMMON Q(3), R(3), PAR(3,3), TAU(3,3), G(3,3), ALP(3,3), TCK
      dimension Y(6),YOLD(4), DY(4), DYOLD(4), DMAT(4,4), text(3,40),
     &          P(3,3),Ys(6) 

      open(5, file='sample-2.txt', status='unknown')
      open(8, file='cal-tie-2.dat',status='unknown')
      open(9, file='cal-spn-2.dat',status='unknown')
      read(5,799) model, tc, step, y1, y3
 799  format(I1, 9X, 4F10.1)
      read(5,800) (R(I), Q(I),(TEXT(I,J),J=1,40), I=1,3)
 800  format(2F10.4,40A1)
      read(5,801) P(1,2),P(2,1),ALP(1,2)
      read(5,801) P(1,3),P(3,1),ALP(1,3)
      read(5,801) P(2,3),P(3,2),ALP(2,3)
 801  format(3F10.4)
      if(dabs(alp(1,2)).lt.1.d-10) alp(1,2)=.2
      if(dabs(alp(1,3)).lt.1.d-10) alp(1,3)=.2
      if(dabs(alp(2,3)).lt.1.d-10) alp(2,3)=.2

      write(*,799) model, tc, step, y1, y3
      write(*,800) (R(I), Q(I),(TEXT(I,J),J=1,40), I=1,3)
      write(*,801) P(1,2),P(2,1),ALP(1,2)
      write(*,801) P(1,3),P(3,1),ALP(1,3)
      write(*,801) P(2,3),P(3,2),ALP(2,3)

      write(*,909) 
 909  format(/)
      if(model.eq.1) write(6,910)
 910  format(' binodal curve construction using uniquac',/)
      if(model.eq.2) write(6,911)
 911  format(' binodal curve construnction using nrtl',/)
      if(model.eq.1) write(6,912)(i,(text(i,j),j=1,40),R(i),Q(i),i=1,3)
 912  format(/,49x,'R',7x,'Q',/,3(/,' (',I1,')',40A1,F8.4,F7.3))      
      if(model.eq.2) write(6,913)(i,(text(i,j),j=1,40),i=1,3)
 913  format(/,3(/,' (',I1,')',40A1))
      write(6,914)TC
 914  format(/,' Temperature = ', F4.1,'DEG C')
      TCK = TC + 273.15D0
      write(6,915)
 915  format(/,' Interaction Parameters: ')
      I1 = 1
      I2 = 2
      I3 = 3
      if(model.eq.1)write(6,916)I1,I2,P(1,2),P(2,1),I1,I3,P(1,3),
     & P(3,1),I2,I3,P(2,3),P(3,2)
 916  format(/,' I J',5x,' A(I,J)',6x,' A(J,I)',/,10x,'Kelvin',6x,
     & 'Kelvin',/,3(/,1x,2I2,3x,2G12.5))
      if(model.eq.2)write(6,917)I1,I2,P(1,2),P(2,1),ALP(1,2),I1,I3,
     & P(1,3),P(3,1),ALP(1,3),I2,I3,P(2,3),P(3,2),ALP(2,3)
 917  format(/,' I J',5x,' A(I,J)',6x,' A(J,I)',7x,'ALPHA',/,10x,
     & 'Kelvin',6x,'Kelvin',/,3(/,1x,2I2,3x,2G12.5,F9.4))
      ALP(2,1) = ALP(1,2)
      ALP(3,1) = ALP(1,3)
      ALP(3,2) = ALP(2,3)
      do 120 i=1,3
      P(i,i)=0.D0
      ALP(i,i)=0.D0
      do 120 j=1,3
      TAU(i,j)=P(i,j)/TCK
      G(i,j)   = DEXP(-ALP(i,j)*TAU(i,j))
 120  PAR(i,j) = DEXP(-TAU(i,j))      
      write(*,'("TAU is",/,3(3(F8.5,2x),/))')((TAU(i,j),j=1,3),i=1,3)
      if(step.lt.1.D-10) step=.02D0
      write(6,918)step
 918  format(/,' step length used in the binodal curve construction ='
     & 2PF4.1,' Mole PCT')
      write(6,919)Y1,Y3
 919  format(/,' initial Guess at mutual solubilities of components
     & (1) and (3):',//,' (1) in (3):',F5.1,' Mole PCT',/,
     &                  ' (3) in (1):',F5.1,' Mole PCT')
      write(6,909)
      write(6,920)
 920  format(' computed binodal curve:',//,9x,'Left Phase',17x,
     & 'Right Phase',/,2(5x,'(1)      (2)      (3)     '),/)
      do 1 i = 1,4
      do 1 j = 1,4
 1    DMAT(i,j) = 0.D0
      IRUND = 200
      IC = 0
      NOLD = 2
      N = 0
      Y(1) = 1.D0 - Y3/100.D0
      Y(2) = 0.D0
      Y(3) = Y1/100.D0
      Y(4) = 0.D0


 12   call solve(model,Y,DY,NOLD,NEW,NITER,N)

      if(NITER.le.10) goto 16
      if(N.GT.0) goto 19
      write(6,902)
 902  format('The base line calculation did not converge in 10
     & iterations')
      write(6,901)
 901  format('Try an initial guess at base line concentration Y1',/,
     & 'and Y3.   If this does not help, perhaps the base line',/,
     & 'parameters correspond to only one liquid phase')
      write(6,925)
 925  format('Are the components positioned as on figure 1 of
     &the introduction?')
      goto 3            
 19   if(IHALF.lt.5) goto 20
      write(6,900)
 900  format('The next tie line calculation did not converge in 10  
     &iterations')
      goto 3                ! end
 20   IHALF = IHALF + 1
      ST = ST/2.D0
      goto 17

 16   if(DABS(Y(1)-Y(3))+DABS(Y(2)-Y(4)).gt.1.D-5) goto 21
      if(N.gt.0) goto 19
      write(6,903)
 903  format('The calculated concentrations on the base line are',/,
     & 'identical in the two phases')
      write(6,901)
      write(6,925)
      goto 3
 21   N=N+1
      IHALF=0
      Y(5)=1.D0-Y(1)-Y(2)
      Y(6)=1.D0-Y(3)-Y(4)

C     stone for spinodal
      call SPINODAL(model,Y,Ys)
      write(8,899) Y(2),Y(5),Y(4),Y(6) 
      write(9,899) Ys(2),Ys(3),Ys(5),Ys(6)
 899  format(4(F5.3,1x))

      if(IC.eq.1) goto 3
      if(IC.eq.2.and.Y(1).lt.1.D-10) goto 3

C     FIND coefficients of third degress polynomium to extrapolate
C     binodal curve
      
      DYMAX = DABS(DY(NEW))

      do 4 I=1,4
 4    DY(I) = DY(I)/DYMAX    
      if(N.eq.1) goto 5
      STAP = DABS(Y(NEW)-YOLD(NEW))

      if(DY(NEW)*DYOLD(NEW).gt.0.D0) goto 6
      do 7 i=1,4
 7    DY(i) = -DY(i)
 6    if(NEW.eq.NOLD) goto 8
      RR = DY(NEW)/DYOLD(NEW)
      do 9 I=1,4
 9    DYOLD(I) = DYOLD(I)*RR
 8    do 10 I=1,4
      Z = (YOLD(I)-Y(I))/STAP
      DMAT(I,3)=(3.D0*Z+2.D0*DY(I)+DYOLD(I))/STAP
 10   DMAT(I,4)=(2.D0*Z+DY(I)+DYOLD(I))/STAP**2
 5    ST=RUND(Y(NEW),DY(NEW),step,IRUND)

      do 18 i=1,4
      DMAT(i,1)=Y(i)
      DMAT(i,2)=DY(i)
      YOLD(i) = Y(i)
 18   DYOLD(i) = DY(i)

C     Initial Guess at Next TieLine
 17   do 11 i=1,4
      Y(i) = DMAT(i,4)
      do 11 j=1,3
 11   Y(i)=ST*Y(i)+DMAT(i,4-j)

      if(IHALF.gt.0) goto 12

C     Check for End of Binodal Curve
      call TERM(Y,DMAT,IC,NEW)
      NOLD = NEW

      if(IC.eq.0.OR.IC.eq.2) goto 12
      if(IC.eq.1 ) goto 21
      if(IC.eq.-2) write(6,898)
 898  format(' Components (1) and (2) are predicted to be
     & immiscible')
      write(6,925)
      if(IC.eq.-1) write(6,897)
 897  format(' Plait point calculation did not converge')
 3    write(6,909)
      stop
      close(5)
      close(8)
      close(9)
      end



      subroutine solve(model,Y,DY,NOLD,NEW,NITER,N)
      implicit real*8 (A-H,O-Z)
      dimension Y1(3),Y2(3),ACT1(3),ACT2(3),DACT1(3,3),DACT2(3,3),
     & Y(6),DY(4),AMAT(3,5),INO(3)
      NITER = 0
 
C     CONVERGE The tie line corresponding to Y(NOLD) and find the
C     derivatives of the concentrations with respect to Y(NOLD)


 11   NITER = NITER+1
      if(NITER.gt.10) return
      do 2 i=1,4
 2    if(Y(i).lt.0.D0)Y(i)=0.D0
      do 3 i=1,2
      Y1(i) = Y(i)
 3    Y2(i) = Y(i+2)
      if(model.ne.1) goto 4
      call UNIQ(Y1,ACT1,DACT1)
      call UNIQ(Y2,ACT2,DACT2)

      goto 5
 4    call NRTL(Y1,ACT1,DACT1)
      call NRTL(Y2,ACT2,DACT2)
 5    J=0
      do 6 i=1,4
      if(i.eq.NOLD) goto 6
      J=J+1
      INO(J)=I
 6    continue

      do 7 i=1,3
      do 7 j=1,2
      AMAT(i,j)  = DACT1(i,j)-DACT1(i,3)
 7    AMAT(i,j+2)= DACT2(i,3)-DACT2(i,j)
      do 8 i=1,3
      AMAT(i,5) = AMAT(i,NOLD)
      do 9 j=1,3
 9    AMAT(i,j) = AMAT(i,INO(j))
 8    AMAT(i,4) = ACT1(i)-ACT2(i)      


      call GAUSL(3,5,3,2,AMAT)

      RES=0.D0
      do 10 i=1,3
      Y(INO(i)) =Y(INO(i))-AMAT(i,4)
      DY(INO(i))=-AMAT(i,5)
 10   RES=RES+AMAT(i,4)**2

      if(RES.gt.1.D-10) goto 11

C     Check for unstable Phases
      call GMIX(Y1,ACT1,DACT1,ICVEX)
      if(ICVEX.eq.-1)write(6,900)
 900  format('  Left  Phase Unstable in the Following Tie Line')
      call GMIX(Y2,ACT2,DACT2,ICVEX)
      if(ICVEX.eq.-1)write(6,901)
 901  format('  Right Phase Unstable in the Following Tie Line')

C     Find new,the number of the concentration with greatest derivative

      DY(NOLD) = 1.0
      NEW = NOLD
      DYMAX = 1.D0
      do 12 i=1,4
      if(DABS(DY(i)).le.DYMAX) goto 12
      NEW = i
      DYMAX=DABS(DY(i))
 12   continue
      return
      end

      subroutine GMIX(X,ACT,DACT,ICVEX)
      implicit real*8(A-H,O-Z)
      dimension X(3), DG(2), DDG(2,2), ACT(3), DACT(3,3)
C     Check for Stability of Each Phase
      NK = 3
      ICVEX = 1
      X(3)=1.0D0-X(1)-X(2)
      IRETUR=0
      do 2 i=1,NK
 2    if(X(I).lt.1.D-10) IRETUR=1
      if(IRETUR.eq.1) return
      do 5 i=1,NK
      do 5 j=1,NK
 5    DACT(i,j) = DACT(i,j)/ACT(i)
      do 20 i=2,NK
      II = i-1
      do 20 j=2,NK
      JJ = j-1
 20   DDG(II,JJ)=DACT(i,j)-DACT(1,j)-DACT(i,1)+DACT(1,1)
      DET=DDG(1,1)*DDG(2,2)-DDG(2,1)*DDG(1,2)
      if(DET.le.0.D0.OR.DDG(1,1).le.0.D0.OR.DDG(2,2).le.0.D0)ICVEX=-1
      return 
      end

C     stone 
      subroutine SPINODAL(model,Y,Ys)
      implicit real*8(A-H,O-Z)
      dimension Y(6),Ys(6),x1(3),x2(3)
      dimension Z(3),DG(2),DDG(2,2),ACT(3),DACT(3,3)

      NK = 3
      write(*,'(4(F5.3,1x))') Y(2),Y(5),Y(4),Y(6)
      x1(1)=Y(1);x1(2)=Y(2);x1(3)=Y(5)
      x2(1)=Y(3);x2(2)=Y(4);x2(3)=Y(6)  
      IRETUR=0
      do 2 i=1,NK
      if(x1(i).lt.1.D-10) IRETUR=1
 2    if(x2(i).lt.1.D-10) IRETUR=1
      if(IRETUR.eq.1) return

      step = 1.d-5 
      do 22 km=1,2
      ratio=0.0D0
      IBACK=0      
      do 21 nstep=1,50000
         do 18 i=1,3
            if(km.eq.1) z(i)=x1(i)+ratio*(x2(i)-x1(i))
            if(km.eq.2) z(i)=x2(i)+ratio*(x1(i)-x2(i))
18       continue
         if(modle.eq.1) then
            call UNIQ(z,ACT,DACT)
         else
            call NRTL(z,ACT,DACT)
         end if 
         do 5 i=1,NK
         do 5 j=1,NK
  5      DACT(i,j) = DACT(i,j)/ACT(i)
         do 20 i=2,NK
         II = i-1
         do 20 j=2,NK
         JJ = j-1
 20      DDG(II,JJ)=DACT(i,j)-DACT(1,j)-DACT(i,1)+DACT(1,1)
         det=DDG(1,1)*DDG(2,2)-DDG(2,1)*DDG(1,2)
         if(det.gt.0.0D0) then
           if(IBACK.eq.1) exit
           ratio=ratio+step
         endif
         if(det.lt.0.0D0) then
           IBACK=1
           ratio=ratio-step/5000
         endif
 21   continue
      do 19 i=1,3
         if(km.eq.1) then
           Ys(i)=z(i)
         else
           Ys(i+3)=z(i)
         endif
 19   continue 
      write(*,*)'km is',km,' ratio is ', ratio,' det is ',det
 22   continue
      write(*,'("Ys is",/,4(F5.3,1x),/)')Ys(2),Ys(3),Ys(5),Ys(6) 
      return
      end

      function RUND(Y,DY,S,IRUND)
      implicit real*8(A-H,O-Z)
C     find round value for concentration step

      X = Y+S*DY+1.D-8*DY**2
      IX = IRUND*X
      Z = DFLOAT(IX)/IRUND-Y
      RUND=DABS(Z)
      return
      end

      subroutine TERM(Y,DMAT,ICOND,NEW)
      implicit real*8(A-H,O-Z)
      dimension Y(6),DMAT(4,4),A(4)
C     Check for end of binodal curve
      if(Y(1).lt.1.D-10.OR.Y(3).lt.1.D-10)       goto 1
      if(Y(1)+Y(2).gt.1.D0.OR.Y(3)+Y(4).gt.1.D0) goto 2
      if(Y(1)+Y(2)-.01D0.lt.Y(3)+Y(4).AND.Y(1)-.01D0.lt.Y(3)) goto 3
      return
 1    ICOND=2
      DS=DMAT(1,1)/(DMAT(1,1)-Y(1))
      DO 5 I=1,4
 5    Y(I)=DMAT(I,1)+DS*(Y(I)-DMAT(I,1))
      Y(1)=0.D0
      NEW=1
      RETURN
 2    ICOND=-2
      RETURN
 3    ICOND=1
      ND=2+NEW
      IF(ND.GT.4) ND=ND-4
      DO 6 I=1,4
 6    A(I)=DMAT(NEW,I)-DMAT(ND,I)
      DS=0.D0
      NITER=0
 7    NITER=NITER+1
      IF(NITER.LE.10) GOTO 8
      ICOND=-1
      RETURN
 8    F=((A(4)*DS+A(3))*DS+A(2))*DS+A(1)
      DF=(3.D0*A(4)*DS+2.D0*A(3))*DS+A(2)
      DF=-F/DF
      DS=DS+DF
      IF(DABS(DF).GT.1.D-6) GOTO 7
      DO 9 I=1,4
      Y(I) = DMAT(I,4)
      DO 9 J=1,3
 9    Y(I)=Y(I)*DS+DMAT(I,4-J)
      RETURN
      END

      SUBROUTINE UNIQ(X,ACT,DACT)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON Q(3),R(3),PAR(3,3),TAU(3,3),G(3,3),ALP(3,3),TCK
      DIMENSION X(3),ACT(3),DACT(3,3),THETA(3),PHI(3),THS(3),
     & QI(3),QIX(3),RI(3),PARA(3,3),PARB(3,3),GAM(3),QID(3)

C      CALCULATION OF ACTIVITIES AND DERIVATIVES OF ACTIVITIES WITH
C      RESPECT TO CONCENTRATIONS USING THE UNIQUAC EQUATION
      NK=3
      X(3)=1.D0-X(1)-X(2)

      IF(X(3).LT.0.D0) X(3)=0.D0
      NCOR=5
      THETS=0.D0
      PHS=0.D0
      DO 10 I=1,NK
      THETA(I)=Q(I)*X(I)
      THETS=THETS+THETA(I)
      PHI(I)=R(I)*X(I)
 10   PHS=PHS+PHI(I)
      DO 20 I=1,NK
      THETA(I)=THETA(I)/THETS
      PHI(I)=PHI(I)/PHS
      RI(I)=R(I)/PHS
      QIX(I)=Q(I)/THETS
      QI(I)=RI(I)/QIX(I)
 20   QID(I)=1.D0-QI(I)
      DO 30 I=1,NK
      THS(I)=0.D0
      DO 30 J=1,NK
 30   THS(I)=THS(I)+PAR(J,I)*THETA(J)
      DO 40 I=1,NK
      GA=1.0D0-RI(I)
      VAL=DLOG(QI(I)**NCOR*THS(I))
      GB=NCOR*QID(I)+VAL-1.D0
      DO 45 J=1,NK
      PARA(I,J)=PAR(I,J)/THS(J)
      PARB(I,J)=PARA(I,J)*THETA(J)
 45   GB=GB+PARB(I,J)
      GAM(I)=DEXP(GA-Q(I)*GB)*RI(I)
 40   ACT(I)=X(I)*GAM(I)
      DO 50 I=1,NK
      DO 50 J=I,NK
      PSUM=1.D0-PARA(I,J)-PARA(J,I)
      DO 55 K=1,NK
 55   PSUM=PSUM+PARA(I,K)*PARB(J,K)
      PSUM=PSUM-NCOR*QID(I)*QID(J)
      DACT(I,J)=Q(I)*QIX(J)*PSUM+(1.D0-RI(I))*(1.D0-RI(J))
 50   DACT(J,I)=DACT(I,J)
      DO 60 I=1,NK
      DO 60 J=1,NK
      DACT(I,J)=ACT(I)*DACT(I,J)
      IF (J.EQ.I) DACT(I,J) = DACT(I,J)+GAM(I)
 60   CONTINUE
      RETURN
      END

      SUBROUTINE NRTL(X,ACT,DACT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON Q(3),R(3),PAR(3,3),TAU(3,3),G(3,3),ALP(3,3),TCK
      DIMENSION X(3),ACT(3),DACT(3,3),GAM(3),G1(3,3),TAU1(3,3),
     & S(3,3), G2(3,3),A1(3),B1(3)

      NK=3
      X(3)=1.D0-X(1)-X(2)
      IF(X(3).LT.0.D0)X(3)=0.D0
      DO 20 I=1,NK
      AA=0.D0
      BB=0.D0
      DO 30 J=1,NK
      Z=G(J,I)*X(J)
      AA=AA+Z
 30   BB=BB+Z*TAU(J,I)
      A1(I)=AA
      B1(I)=BB
      GAM(I) = BB/AA
      DO 20 J=1,NK
      G1(J,I)=G(J,I)/AA
      TAU1(J,I)=TAU(J,I)-GAM(I)
      G2(J,I)=G1(J,I)*TAU1(J,I)
 20   S(J,I)=X(I)*G2(J,I)
      DO 40 I=1,NK
      DO 50 J=1,NK
 50   GAM(I)=GAM(I)+S(I,J)
      GAM(I)=DEXP(GAM(I))
 40   ACT(I)=X(I)*GAM(I)
      DO 60 I=1,NK
      DO 60 J=I,NK
      PSUM=G2(J,I)+G2(I,J)
      DO 65 K=1,NK
 65   PSUM=PSUM-G1(I,K)*S(J,K)-G1(J,K)*S(I,K)
      DACT(I,J)=PSUM
 60   DACT(J,I)=PSUM
      DO 70 I=1,NK
      DO 70 J=1,NK
      DACT(I,J)=DACT(I,J)*ACT(I)
      IF(J.EQ.I)DACT(I,I)=DACT(I,I)+GAM(I)
 70   CONTINUE
      RETURN
      END

      SUBROUTINE GAUSL(ND,NCOL,N,NS,A)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(ND,NCOL)

C     GAUSL solves A*X=B, where A is N*N and B is N*NS, by gaussian
C     elimination with partial pivoting, the matrix(or vector) B is
C     placed adjacent to a in columns N+1 to N+NS
C     A is destroyed, and the resulting Matrix X replaces B

      N1=N+1
      NT=N+NS
      IF(N.EQ.1) GOTO 50

      DO 10 I=2,N
      IP=I-1
      I1=IP
      X=DABS(A(I1,I1))
      DO 11 J=I,N
      IF(DABS(A(J,I1)).LT.X) GOTO 11
      X=DABS(A(J,I1))
      IP=J
 11   CONTINUE
      IF(IP.EQ.I1) GOTO 13

C     ROW INTERCHANGE
      DO 12 J=I1,NT
      X=A(I1,J)
      A(I1,J)=A(IP,J)
 12   A(IP,J)=X

 13   DO 10 J=I,N
      IF(DABS(A(I1,I1)).LT.1.D-10)A(I1,I1)=1.D0
      X=A(J,I1)/A(I1,I1)

      DO 17 K=I1,NT
 17   A(J,K)=A(J,K)-X*A(I1,K)
 10   continue



C     ELIMINATION FINISHED, NOW BACKSUBSTITUTION
 50   DO 20 IP=1,N
      I=N1-IP
      DO 20 K=N1,NT
      IF(DABS(A(I,I)).LT.1.D-10)A(I,I)=1.D0
      A(I,K)=A(I,K)/A(I,I)
      IF(I.EQ.1) GOTO 20
      I1=I-1
      DO 25 J=1,I1
 25   A(J,K)=A(J,K)-A(I,K)*A(J,I)
 20   CONTINUE
      RETURN
      END

