       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION X(1000),Y(1000),Z(1000)


        PI = DACOS(-1.D0)
        

        ETA = 0.4D0
        
        NTOT = 124
C       NTOT = 230


C       AL = ((NTOT+10)*PI/ETA/6.D0)**(1.D0/3.D0)
C       AL = (NTOT*PI/ETA/6.D0)**(1.D0/3.D0)
        
        AL = 6.2D0

        ALH = 0.5D0*AL


        ND = INT(AL) + 1


        DR = 1.D0

         
        PRINT*,"ETA,AL,DR",ETA,AL,DR

          IC = 0

        DO I = 1 , ND-1

         DO J = 1, ND-1

           DO K = 1, ND-1

           
             XTRY = -ALH + (I-0.5)*DR

             YTRY = -ALH + (J-0.5)*DR

             ZTRY = -ALH + (K-0.5)*DR

             IF(XT.LT.-ALH .OR. XT.GT.ALH) GOTO 10
             IF(YT.LT.-ALH .OR. YT.GT.ALH) GOTO 10
             IF(ZT.LT.-ALH .OR. ZT.GT.ALH) GOTO 10

             IF (IC.GT.1) THEN
                DO MN = 1, IC-1

                 TX = X(MN)-XTRY
                 TY = Y(MN)-YTRY
                 TZ = Z(MN)-ZTRY

C                TX = TX - AL*DNINT(TX/AL)
C                TY = TY - AL*DNINT(TY/AL)
C                TZ = TZ - AL*DNINT(TZ/AL)

                 R2 = TX*TX + TY*TY + TZ*TZ

                 IF(R2+1.D-6 .LT. 1.D0) GO TO 10

                 ENDDO


             ENDIF

               IC = IC + 1

               X(IC) = XTRY
               Y(IC) = YTRY
               Z(IC) = ZTRY

               print*,"IC",IC,XTRY,YTRY,ZTRY

               IF(IC.EQ.NTOT) GO TO 20

10        CONTINUE

           ENDDO

          ENDDO        

         ENDDO

          PRINT*,"NOT ENOUGH PARTICLES"
          GOTO 30

20        CONTINUE

          AREA = (AL+1.D0)**2
          AREA = 6.D0*AREA

           OPEN(14,FILE="INPUT", STATUS="UNKNOWN")

           WRITE(14,*) 1000000
           WRITE(14,*) AL,AL,AL
           WRITE(14,*) NTOT
C          WRITE(14,*) 64
           WRITE(14,*) AREA

           CLOSE(14)



           OPEN(14,FILE="POS.TXT", STATUS="UNKNOWN")

            DO I = 1 , NTOT

             WRITE(14,*) X(I),Y(I),Z(I)

            ENDDO

           CLOSE(14)

30        CONTINUE

          END
