       PROGRAM HS
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       PARAMETER (PI=3.141592653589793D0)
       DIMENSION X(1000),Y(1000),Z(1000)
       COMMON/POS/X,Y,Z
       COMMON/TOT/NTOT
       COMMON/BOX/ALX,ALY,ALZ
       COMMON/HBOX/HALX,HALY,HALZ
       COMMON/SIDE/AREA
       COMMON/PRESSURE/PS
       COMMON/P123/PG(250),PC(250)
       COMMON/ATT/ASP,AW
       DATA X/1000*0.D0/
       DATA Y/1000*0.D0/
       DATA Z/1000*0.D0/
       DATA PG/250*0.D0/
       DATA PC/250*0.D0/

          NSEED = TIME()
          CALL SRAND(NSEED)

          OPEN(14,FILE="INPUT",STATUS = "UNKNOWN")

            READ(14,*) NCON
            READ(14,*) ALX,ALY,ALZ
            READ(14,*) NTOT
            READ(14,*) AREA
           CLOSE(14)

           HALX = 0.5D0*ALX
           HALY = 0.5D0*ALY
           HALZ = 0.5D0*ALZ

          CALL SRAND(NSEED)


C         PS = 20.D0
C         PS = 12.D0
C         PS = 1.D0
C         PS = 9.D0
          
          NSKIP = 50

C         AREA = (ALX+1.D0)*(ALY+1.D0) 
C         AREA = AREA*(ALX+1.D0)*(ALZ+1.D0)
C         AREA = AREA*2.D0*(ALY+1.D0)*(ALZ+1.D0)
          
          VBOX = (ALX+1.D0)*(ALY+1.D0)*(ALZ+1.D0)
          print*,"packing fraction", NTOT*PI/6.D0/VBOX
C         print*, "AREA", AREA

    
           DR = RS/100.D0

            OPEN(12, FILE="POS.TXT",STATUS="UNKNOWN")

              DO K = 1 , NTOT

               READ(12,*)X(K),Y(K),Z(K)
   
              ENDDO

            CLOSE(12)

            DO K1 = 1 , NTOT

             IF (X(K1).LT.-HALX .OR. X(K1).GT.HALX) THEN
                  print*,"x out-of-bound", K1, X(K1)
                  GOTO 999
             ENDIF
             IF (Y(K1).LT.-HALY .OR. Y(K1).GT.HALY) THEN
                  print*,"y out-of-bound", K1, Y(K1)
                  GOTO 999
             ENDIF
             IF (Z(K1).LT.-HALZ .OR. Z(K1).GT.HALZ) THEN
                  print*,"z out-of-bound", K1, Z(K1)
                  GOTO 999
             ENDIF

            ENDDO


           DO K1 = 1 , NTOT-1

             DO K2 = K1+1, NTOT

                TX = X(K1)-X(K2) 
                TY = Y(K1)-Y(K2) 
                TZ = Z(K1)-Z(K2) 

C               TX = TX - ALX*DNINT(TX/ALX)
C               TY = TY - ALY*DNINT(TY/ALY)
C               TZ = TZ - ALZ*DNINT(TZ/ALZ)

             DIST2 = TX*TX + TY*TY + TZ*TZ


             IF(DIST2 + 1.D-6 .LT. 1.D0) THEN

                print*,"err-Overlap", K1, K2, DIST2
                print*,"err-Overlap", X(K1), Y(K1), Z(K1)
                print*,"err-Overlap", X(K2), Y(K2), Z(K2)
                print*," "

               STOP

             ENDIF
C             ENDIF

5            CONTINUE

              ENDDO

            ENDDO
            
              
             print*,"simulation begins"
C            PS = 3.D0
C            PS = 15.D0
             PS = 3.D0
C            ASP = 6.D0
             ASP = 0.D0
             AW = 0.D0

C            STEP = 0.003D0*AL
C            STEP = 0.0033D0*(ALX+ALY+ALZ)/3.D0
             STEP = 0.005D0*(ALX+ALY+ALZ)/3.D0
C            STEP = 0.020D0*(ALX+ALY+ALZ)/3.D0
C            STEPV = 0.01*AL
C            STEPV = 0.00058*(ALX+ALY+ALZ)/3.D0
C            STEPV = 0.00058*(ALX+ALY+ALZ)/3.D0
             STEPV = 0.0002*(ALX+ALY+ALZ)/3.D0
C            STEPV = 0.006*(ALX+ALY+ALZ)/3.D0
             TSUMO = 0.D0
             TSUMV = 0.D0
             ACCMO = 0.D0
             ACCMV = 0.D0
             TAV = 0.D0
             TALX = 0.D0
             TALY = 0.D0
             TALZ = 0.D0
             STALX = 0.D0
             STALY = 0.D0
             STALZ = 0.D0

             TCAL = 0.D0


            DO 1000  ICON = 1 , NCON

                PSEL=RAND()

                IF(PSEL .LT. 1.) THEN

                    TSUMO = TSUMO + 1.D0

                  CALL MOVE(ISUC,STEP)

                  IF(ISUC.EQ.0) ACCMO = ACCMO + 1.

                ELSE

                   TSUMV =  TSUMV + 1.D0

                   CALL MVOL(ISUC,STEPV)

                  IF(ISUC.EQ.0) ACCMV = ACCMV + 1.

                ENDIF


                IF(MOD(ICON,NSKIP).EQ.0) THEN
                    CALL OMEGA()
                    TCAL = TCAL + 1.D0
                    TAV = TAV + ALX*ALY*ALZ
                    STALX = STALX + ALX
                    STALY = STALY + ALY
                    STALZ = STALZ + ALZ
C                   NSEED = TIME()
C                   CALL SRAND(NSEED)
C                   DO IR = 1, 200
C                     XX = RAND()
C                   ENDDO
                ENDIF


1000      continue

           print*,"AL",ALX,ALY,ALZ
           print*,ACCMO/TSUMO,ACCMV/TSUMV
           print*, "AVG-Vol", TAV/TCAL
           print*, "AVG-AL", STALX/TCAL,STALY/TCAL,STALZ/TCAL
           print*, "AREA",AREA


          OPEN(14,FILE="INPUT",STATUS = "UNKNOWN")

            WRITE(14,*) NCON
            WRITE(14,*) ALX,ALY,ALZ
            WRITE(14,*) NTOT
            WRITE(14,*) AREA

           CLOSE(14)

           OPEN(12,FILE = "POS.TXT", STATUS = "UNKNOWN")

             DO K = 1 , NTOT

                WRITE(12,*) X(K),Y(K),Z(K)

             ENDDO

            CLOSE(12)

            VOL = ALX*ALY*ALZ

            OPEN(14,FILE = "GR.TXT", STATUS = "UNKNOWN")

              DO I = 1, 250

                 R1 = (I-1)*0.04D0
                 R2 = I*0.04D0

                 V12 = 4.D0*PI*(R2**3-R1**3)/3.D0
                 BID = V12*NTOT*NTOT/VOL
                 PG(I) = PG(I)/BID/TCAL

           WRITE(14,*) 0.5D0*(R1+R2),PG(I),PC(I)/V12/TCAL

              ENDDO

            CLOSE(14)


111    FORMAT(1X,4F15.7)
112    FORMAT(1X,3F12.5)
999    CONTINUE

       END

       SUBROUTINE MOVE(ISUC,STEP)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION X(1000),Y(1000),Z(1000)
       COMMON/POS/X,Y,Z
       COMMON/TOT/NTOT
       COMMON/BOX/ALX,ALY,ALZ
       COMMON/ATT/ASP,AW
       COMMON/HBOX/HALX,HALY,HALZ
        
       ISUC = 1 


       J = INT(NTOT*RAND()) +1

          XTRY = X(J) + STEP*(2.*RAND()-1.)       
          YTRY = Y(J) + STEP*(2.*RAND()-1.)       
          ZTRY = Z(J) + STEP*(2.*RAND()-1.)       

C         print*,"TRY",XTRY,YTRY,ZTRY
C         print*,"TRY",XTRY-X(J),YTRY-Y(J),ZTRY-Z(J)
C         print*," "

          IF(XTRY.LT.-HALX .OR. XTRY.GT.HALX) RETURN
          IF(YTRY.LT.-HALY .OR. YTRY.GT.HALY) RETURN
          IF(ZTRY.LT.-HALZ .OR. ZTRY.GT.HALZ) RETURN


C excluded volume interaction

         DO K = 1, NTOT 

            IF(J.NE.K) THEN

               TX = XTRY - X(K)
               TY = YTRY - Y(K)
               TZ = ZTRY - Z(K)

C              TX = TX1 - ALX*DNINT(TX1/ALX)
C              TY = TY1 - ALY*DNINT(TY1/ALY)
C              TZ = TZ1 - ALZ*DNINT(TZ1/ALZ)

               DIST2 = TX*TX + TY*TY + TZ*TZ

                 IF(DIST2+1.D-6 .LT. 1.D0) RETURN 


             ENDIF        

          ENDDO


          EOLD = 0.D0
          ENEW = 0.D0

          DO K1 = 1, NTOT

           IF(J.NE.K1) THEN

              TX = X(K1)-X(J)
              TY = Y(K1)-Y(J)
              TZ = Z(K1)-Z(J)

              D2 = TX*TX + TY*TY + TZ*TZ
              D6 = D2*D2*D2

              EOLD = EOLD - 1.D0/D6

              TX = X(K1)-XTRY
              TY = Y(K1)-YTRY
              TZ = Z(K1)-ZTRY

              D2 = TX*TX + TY*TY + TZ*TZ
              D6 = D2*D2*D2

              ENEW = ENEW - 1.D0/D6

           ENDIF

          ENDDO

            EOLD = ASP*EOLD
            ENEW = ASP*ENEW
            TT1 = 0.D0
            TT2 = 0.D0
            EWOLD=0.D0
            EWNEW = 0.D0

C           print*,"ASP",ENEW,EOLD

            TX1 = (X(J)-HALX-0.5D0)**6
            TX2 = (X(J)+HALX+0.5D0)**6

            TY1 = (Y(J)-HALY-0.5D0)**6
            TY2 = (Y(J)+HALY+0.5D0)**6

            TZ1 = (Z(J)-HALZ-0.5D0)**6
            TZ2 = (Z(J)+HALZ+0.5D0)**6

            TT1 = 1.D0/TX1 +1.D0/TX2 +1.D0/TY1+1.D0/TY2
     +             + 1.D0/TZ1 + 1.D0/TZ2
         
            EWOLD = EWOLD - TT1

C           print*,"OLD",TX1,TX2
C           print*,"OLD",TY1,TY2

            TX1 = (XTRY-HALX-0.5D0)**6
            TX2 = (XTRY+HALX+0.5D0)**6

            TY1 = (YTRY-HALY-0.5D0)**6
            TY2 = (YTRY+HALY+0.5D0)**6

            TZ1 = (ZTRY-HALZ-0.5D0)**6
            TZ2 = (ZTRY+HALZ+0.5D0)**6

            TT2 =1.D0/TX1 +1.D0/TX2 +1.D0/TY1+1.D0/TY2
     +           + 1.D0/TZ1 + 1.D0/TZ2
         
            EWNEW = EWNEW - TT2

C           print*,"NEW",TX1,TX2
C           print*,"NEW",TY1,TY2
C           print*," "

            EOLD = EOLD + EWOLD*AW
            ENEW = ENEW + EWNEW*AW

C           print*,"AW",EWNEW*AW,EWOLD*AW
C           print*,"ENE",ENEW,EOLD
C           print*," "


            IF(RAND().GT.DEXP(EOLD-ENEW)) RETURN

C


         ISUC = 0

            X(J) = XTRY
            Y(J) = YTRY
            Z(J) = ZTRY
            
 
       RETURN
       END 

       SUBROUTINE MVOL(ISUC,STEPV)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION X(1000),Y(1000),Z(1000)
       DIMENSION XN(1000),YN(1000),ZN(1000)
       COMMON/POS/X,Y,Z
       COMMON/TOT/NTOT
       COMMON/BOX/ALX,ALY,ALZ
       COMMON/ATT/ASP,AW
       COMMON/HBOX/HALX,HALY,HALZ
       COMMON/SIDE/AREA
       COMMON/PRESSURE/PS
        
       ISUC = 1 

       AL_NEW_X = ALX 
       AL_NEW_Y = ALY 
       AL_NEW_Z = ALZ

       IOP = INT(3*RAND()) + 1 

       IF(IOP .EQ. 1) THEN

        
       AL_NEW_X = ALX + STEPV*(2.*RAND()-1.)

       IF (AL_NEW_X .LE. 0.D0) RETURN

C      AL_NEW_Y = ALY + STEPV*(2.*RAND()-1.)
C      AL_NEW_Z = ALZ + STEPV*(2.*RAND()-1.)

C      IF(AL_NEW_Y .LE. 0.D0) RETURN
C      IF(AL_NEW_Z .LE. 0.D0) RETURN

C      AL_NEW_X = 0.5D0*AREA-(AL_NEW_Y+1.D0)*(AL_NEW_Z+1.D0)
C      AL_NEW_X = (AL_NEW_X/(AL_NEW_Y+AL_NEW_Z+2.D0))-1.D0

       ELSEIF(IOP .EQ. 2) THEN

       AL_NEW_Y = ALY + STEPV*(2.*RAND()-1.)

       IF(AL_NEW_Y .LE. 0.D0) RETURN

C      AL_NEW_X = ALX + STEPV*(2.*RAND()-1.)
C      AL_NEW_Z = ALZ + STEPV*(2.*RAND()-1.)

C      IF(AL_NEW_X .LE. 0.) RETURN
C      IF(AL_NEW_Z .LE. 0.) RETURN

C      AL_NEW_Y = 0.5D0*AREA-(AL_NEW_X+1.D0)*(AL_NEW_Z+1.D0)
C      AL_NEW_Y = (AL_NEW_Y/(AL_NEW_X+AL_NEW_Z+2.D0))-1.D0

       ELSEIF(IOP .EQ. 3) THEN

       AL_NEW_Z = ALZ + STEPV*(2.*RAND()-1.)

       IF(AL_NEW_Z .LE. 0.D0) RETURN

C      AL_NEW_X = ALX + STEPV*(2.*RAND()-1.)
C      AL_NEW_Y = ALY + STEPV*(2.*RAND()-1.)

C      IF(AL_NEW_X .LE. 0.) RETURN
C      IF(AL_NEW_Y .LE. 0.) RETURN

C      AL_NEW_Z = 0.5D0*AREA-(AL_NEW_X+1.D0)*(AL_NEW_Y+1.D0)
C      AL_NEW_Z = (AL_NEW_Z/(AL_NEW_X+AL_NEW_Y+2.D0))-1.D0

       ENDIF

C      print*,"NEW",AL_NEW_X,AL_NEW_Y,AL_NEW_Z

       Ratio_X = (AL_NEW_X+1.D0)/(ALX+1.D0)
       Ratio_Y = (AL_NEW_Y+1.D0)/(ALY+1.D0)
       Ratio_Z = (AL_NEW_Z+1.D0)/(ALZ+1.D0)

C      print*,"IOP",IOP
C      print*,Ratio_X,Ratio_Y,Ratio_Z

              HALX_NEW = 0.5D0*AL_NEW_X
              HALY_NEW = 0.5D0*AL_NEW_Y
              HALZ_NEW = 0.5D0*AL_NEW_Z

        DO J = 1 , NTOT
          XN(J) = X(J) * Ratio_X      
        IF(XN(J).GT.HALX_NEW .OR. XN(J).LT.-HALX_NEW) RETURN
          YN(J) = Y(J) * Ratio_Y     
        IF(YN(J).GT.HALY_NEW .OR. YN(J).LT.-HALY_NEW) RETURN
          ZN(J) = Z(J) * Ratio_Z     
        IF(ZN(J).GT.HALZ_NEW .OR. ZN(J).LT.-HALZ_NEW) RETURN

        ENDDO

C excluded volume interaction

         DO K1 = 1, NTOT - 1

            DO K2 = K1+1, NTOT

               TX = XN(K1) - XN(K2)
               TY = YN(K1) - YN(K2)
               TZ = ZN(K1) - ZN(K2)

C              TX = TX1 - AL_NEW_X*DNINT(TX1/AL_NEW_X)
C              TY = TY1 - AL_NEW_Y*DNINT(TY1/AL_NEW_Y)
C              TZ = TZ1 - AL_NEW_Z*DNINT(TZ1/AL_NEW_Z)

               DIST2 = TX*TX + TY*TY + TZ*TZ

                 IF(DIST2+1.D-6 .LT. 1.D0) RETURN 

            ENDDO

          ENDDO

          EOLD = 0.D0
          ENEW = 0.D0

          DO K1 = 1, NTOT-1

            DO K2 = K1+1, NTOT

              TX = X(K1)-X(K2)
              TY = Y(K1)-Y(K2)
              TZ = Z(K1)-Z(K2)

              D2 = TX*TX + TY*TY + TZ*TZ
              D6 = D2*D2*D2

              EOLD = EOLD - 1.D0/D6

              TX = XN(K1)-XN(K2)
              TY = YN(K1)-YN(K2)
              TZ = ZN(K1)-ZN(K2)

              D2 = TX*TX + TY*TY + TZ*TZ
              D6 = D2*D2*D2

              ENEW = ENEW - 1.D0/D6

            ENDDO

          ENDDO

           EOLD = EOLD*ASP
           ENEW = ENEW*ASP

C          print*," "
C          print*,"ASP",ENEW,EOLD

C             HALX_NEW = 0.5D0*AL_NEW_X
C             HALY_NEW = 0.5D0*AL_NEW_Y
C             HALZ_NEW = 0.5D0*AL_NEW_Z

            TT1 = 0.D0
            TT2 = 0.D0
            EWOLD = 0.D0
            EWNEW = 0.D0

          DO K1 = 1 , NTOT

            TX1 = (X(K1)-HALX-0.5D0)**6
            TX2 = (X(K1)+HALX+0.5D0)**6

            TY1 = (Y(K1)-HALY-0.5D0)**6
            TY2 = (Y(K1)+HALY+0.5D0)**6

            TZ1 = (Z(K1)-HALZ-0.5D0)**6
            TZ2 = (Z(K1)+HALZ+0.5D0)**6

            TT1 =1.D0/TX1 +1.D0/TX2 +1.D0/TY1+1.D0/TY2
            TT1 = TT1 + 1.D0/TZ1 + 1.D0/TZ2
         
           EWOLD = EWOLD - TT1

            TX1 = (XN(K1)-HALX_NEW-0.5D0)**6
            TX2 = (XN(K1)+HALX_NEW+0.5D0)**6

            TY1 = (YN(K1)-HALY_NEW-0.5D0)**6
            TY2 = (YN(K1)+HALY_NEW+0.5D0)**6

            TZ1 = (ZN(K1)-HALZ_NEW-0.5D0)**6
            TZ2 = (ZN(K1)+HALZ_NEW+0.5D0)**6

            TT2 =1.D0/TX1 +1.D0/TX2 +1.D0/TY1+1.D0/TY2
            TT2 = TT2 + 1.D0/TZ1 + 1.D0/TZ2
         
           EWNEW = EWNEW - TT2

          ENDDO

             EOLD = EOLD + EWOLD*AW
             ENEW = ENEW + EWNEW*AW

C           print*,"AW",EWNEW*AW,EWOLD*AW

         VOLD = (ALX+1.D0)*(ALY+1.D0)*(ALZ+1.D0)
         VNEW = (AL_NEW_X+1.D0)*(AL_NEW_Y+1.D0)*(AL_NEW_Z+1.D0)

C        print*,"VOLD,VNEW-1",VNEW/VOLD, -NTOT*DLOG(VNEW/VOLD)
C        print*,"ENEW,EOLD",ENEW-EOLD, ENEW, EOLD
C        print*," "

         DH = PS*(VNEW-VOLD) - NTOT*DLOG(VNEW/VOLD)+(ENEW-EOLD)

C        print*,"DH",DH,PS*(VOLD-VNEW),NTOT*DLOG(VNEW/VOLD),ENEW-EOLD
C        print*,"exp(-DH)",DEXP(-DH)
C        print*, " "

         IF(RAND() .GT. DEXP(-DH)) RETURN

         ISUC = 0

         ALX = AL_NEW_X
         ALY = AL_NEW_Y
         ALZ = AL_NEW_Z

         HALX = 0.5D0*ALX
         HALY = 0.5D0*ALY
         HALZ = 0.5D0*ALZ

         DO J = 1 , NTOT

            X(J) = XN(J)
            Y(J) = YN(J)
            Z(J) = ZN(J)

         ENDDO
            

       RETURN
       END 

      SUBROUTINE OMEGA()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       DIMENSION X(1000),Y(1000),Z(1000)
       COMMON/POS/X,Y,Z
       COMMON/BOX/ALX,ALY,ALZ
       COMMON/HBOX/HALX,HALY,HALZ
       COMMON/TOT/NTOT
       COMMON/P123/PG(250),PC(250)


         DO K1 = 1 , NTOT-1

           DO K2 = K1 + 1, NTOT

             TX = X(K1)-X(K2)
             TY = Y(K1)-Y(K2)
             TZ = Z(K1)-Z(K2)

C            TX = TX1 - ALX*DNINT(TX1/ALX)
C            TY = TY1 - ALY*DNINT(TY1/ALY)
C            TZ = TZ1 - ALZ*DNINT(TZ1/ALZ)

             RT = DSQRT(TX*TX + TY*TY + TZ*TZ)

             I = INT(RT/0.04D0) + 1


             IF(I.GT.0 .AND. I.LE.250) THEN

                PG(I) = PG(I) + 2.D0

             ENDIF

            ENDDO
          ENDDO

          DO K1 = 1 , NTOT

            D2 = X(K1)*X(K1) + Y(K1)*Y(K1) + Z(K1)*Z(K1)

            D2 = DSQRT(D2)

            I = INT(D2/0.04D0) + 1

            PC(I) = PC(I) + 1.D0

          ENDDO

         RETURN
         END   

