       program pendulum2_see_formula2_in Celletti_2022
       implicit double precision (a-h ,o-z)
       dimension x(2), v1(2), v2(2), adata(41), fli(1), outdat(3,3)
       dimension rot(2), birave(2)
       integer num, denum
       logical isnres
       
       ! This program computes the Poincare map, FLIs and the rotation number for the Hamiltonian described by equation (2) in [1] (see below).
       ! Modifying the subroutines FUN and FUNG, one can easy adapt this program to Hamiltonians given by the formulas (1), (3) or (6) in [1].
       ! The rotation number is computed both by using the classical definition (see formula (7) in [1] ) and the weighted Birkhoff averages (see formula (8) in [1] or formula (6) in [2]).


       ! References:
       ! [1]   Celletti et al. (2022): Classification of regular and chaotic motions in Hamiltonian systems with deep learning, Scientific Reports
       ! [2]   Sander  and Meiss (2021): Birkhoff averages and rotational invariant circles for area-preserving maps, Physica D, 411.
       
      ! Inputs (are read from the file "input_pen_grid.txt"):
           ! ttime - total time of the integration
           ! h0 - step size of the integration
           ! eps - parameter entering the Hamiltonian (2) (see [1]).
      

      ! Outputs: the time series of Poincare type, FLIs, rotation number and index (see [1]) for a grid of (Mpoint+1)x(Npoint+1) initial conditions equally spaced in the (x,y) plane.
                ! Here, Npoint=Mpoint=100, so the grid has 10201 initial  conditions.
         ! The results are included in the files: "Poincare.plt", "fli.plt", "brot_rotnum2.plt", "index_fli_rotnum_birav.plt"
              ! The file "Poincare.plt": the time series of Poincare type for 10201 initial conditions. The blank lines in the Poincare.plt file sepatates the time series. Column 1= time, column 2=x, column 2=y.
              ! "fli.plt" - first column=x, second column=y, third column=fli (the file is used to represent graphically the data in gnuplot.)
              ! "brot_rotnu2.plt" - first column=x, second column=y, third column=rotation number (computed with the formula 8 in [1]), fourth column rotation number (computed with the formula 7 in [1]).
              ! "index_fli_rotnum_birav.plt" - So, in conclusion this new file: index_fli_rotnum_birav_birrefined_v3.plt contains the following information:
                  !a) Columns 1, 2, 3: x, y and FLI.
                  !b) Column 4: rotation number (frequency) computed using the formula (8) in [1] (the weighted Birkhoff averages);
                  !c) Columns 5 and 6: the index is computed by using the FLI (with and without the class of uncertain motions (index -1)).
                  !d) Columns 7 and 8: index computed by using the frequency given by the formulae (7) and (8) in [1] and the property that the rotation number is monotone with respect to the action y.
                  !e) Column 9: index  computed by using the indexing method II: convergence.
                  !f) Column 10: index  computed by using the indexing method II: convergence and the values of the rotation number (in case of convergence). Here, the following coding is used:
                      ! -1=uncertain 0=librational with the rotation number equal to 0; 0.25=librational with the rotation number equal to 1/4; 0.33=librational with the rotation number equal to 1/3;
                      ! 0.5=librational with the rotation number equal to 1/2; 0.67=librational with the rotation number equal to 2/3; 0.75=librational with the rotation number equal to 3/4;
                      ! 1=librational with the rotation number equal to 1; 5=chaotic ; 6=nonrotational.
              
      ! Locals: x(2)- the action angles variables
                !v1(2), v2(2) -  the tangent vectors
                ! adata(41) - the parameters of the differential equations
                ! fli(1) - the FLI
                ! outdat (3,3) - the action y, the roatational number (formula (7) in [1]), the rotational number (formula (8) in [1]) for distinct values of the action y. These are needed to compute the index by using the motonicity property.
                ! rotnum, rotnu2, brot - the rotation number computed in various ways
                ! findex, f2index, aindex, a2index, bindex,  b2index - index computed by using fli, frequency analysis (monotonicity), frequency analysis (convergence)
                !  rot(2), birave(2) - local variables used for storing or computing the rotational number
                !  num, denum - nominator and denominator returned by the subroutine SMALDE (see below)
                ! isnres - used in indexing rorational motions very close to separatrix
                ! q0, p0 - the initial angle x and the initial action y
                ! tsid - the time t
                
       
       j80=18
       j60=12
       j100=14
       j811=30



     
      open(j811,file='brot_rotnu2.plt',
     !status='unknown', form='formatted')
     
       open(j60,file='fli.plt',
     !status='unknown', form='formatted')
        open(j80,file='Poincare.plt',
     !status='unknown', form='formatted')
        open(j100,file='index_fli_rotnum_birav.plt',
     !status='unknown', form='formatted')



        pi=4.d0*datan(1.d0)
        pi2=8.d0*datan(1.d0)
        degree=pi/180.d0




*************************************************************************
*      INPUT DATA from file input.txt
*************************************************************************

        open (6, file='input_pen_grid.txt', form='formatted') !modify input

        read(6,*)
        read(6,*)
        read(6,*)
        read(6,*)
        read(6,*) ttime
        read(6,*) h0
        read(6,*)
        read(6,*)
        read(6,*) eps

        h=pi*2.d0*h0


        N=int(ttime)

        count=0.d0

        Npoint=100
        Mpoint=100
        


*       do iM=1,Mpoint+1
       do iM=1,Mpoint+1
       q0=(0.d0+(dble(iM-1)*360.d0)/dble(Mpoint))*degree

       write(*,*)'iM=',iM
       
        do ia=1,Npoint+1


         p0=-0.5d0 +(dble(ia-1)*2.7d0)/dble(Npoint)



        x(1)=p0
        x(2)=q0

        v1(1)=1.d0
        v1(2)=0.d0

        v2(1)=0.d0
        v2(2)=1.d0

*       pnorm=0.d0
*       do ii=1,2
*         pnorm=pnorm+v1(ii)**2
*       enddo
*       pnorm = dsqrt(pnorm)
*       do ii=1,2
*         v1(ii)=v1(ii)/pnorm
*       enddo

       fli(1) = -10.d0

*********************************************************************
       isnres=.FALSE.
       rotnum=0.d0  ! the rotation number (frequency) computed by the formula (7) from the draft
       rotnu2=0.d0  ! the rotation number (frequency) computed numerically
       brot=0.d0   ! the rotation number (frequeny) is computed using the formula (8) from the draft
       
       ssand=ssander(N/2)  ! constant S see Sander and Meiss, formula (7)
       ssand2=ssander(N)  ! constant S see Sander and Meiss, formula (7)
       birave(1)=0.d0   ! WB_T(h)(z)       see formula (6) in Sander and Meiss
       birave(2)=0.d0   !! WB_T(h)(f^T(z))   see formula (9)  Sander and Meiss


       findex=-1.d0    ! index computed by using the FLI:  -1=uncertain motion, 0=chaotic orbit, 1=rotational motion, 2=librational motion
       aindex=-1.d0
       a2index=-1.d0    ! index computed by using the frequency analysis, method 1, based on rotnum
       bindex=-1.d0    ! index computed by using the frequency analysis, method 2, based on the ideas from Sander and Meiss
       b2index=-1.d0
       num=100
       denum=100

************************************************************************
************************************************************************
*  Integrating the equations of motion and the variational equations by using the Runge-Kutta method

       tsid=0.d0
        do i=1,N
*           count= dble(i-1)

        count=dmod(tsid, pi2)
        if ((count .lt. h) .and. (count .gt. -h)) then
        qq=dmod(x(2),pi*2.d0)
        if (qq .lt. 0.d0) then
         qq=qq+pi*2.d0
         end if
         write(j80,*) tsid, qq/degree, x(1)
        end if
          q2temp=x(2)
          call rukugf(x, v1, v2, tsid, eps, h)

         if (p0 .lt. 0.5d0) then
         qtemp=dmod(x(2), pi2)
         else
         qtemp=dmod(x(2)-tsid, pi2)
         end if
         
         if (qtemp .lt. 0.d0) then
         qtemp=qtemp+pi2
         end if
       if ((qtemp .gt. 160.d0*degree) .and.
     -     (qtemp .lt. 200.d0*degree)) then
        isnres=.TRUE.
        end if
        

        
         tsid=tsid+h
         rotnum=rotnum+x(1)
         rotnu2=rotnu2+(x(2)-q2temp)/h
         brot=brot+ (gsander(dble(i)/(dble(N)))/ssand2)*x(1)
         if(i .le. N/2) then
*         birave(1)=birave(1)+ (gsander(dble(i*2)/(dble(N)))/ssand)*
*     - cos(x(2))*cos(x(1)*pi2)
         birave(1)=birave(1)+ (gsander(dble(i*2)/(dble(N)))/ssand)*x(1)

        else
       birave(2)=birave(2)+ (gsander(dble(i*2-N)/(dble(N)))/ssand)*x(1)
       end if
         
          call flicom(fli, v1, v2)



          
          if (i .eq. int(dble(N)/3.d0)) then
          rot(1)=rotnum/dble(i)                        ! rot(1)= will store the min value of the frequency
          rot(2)=rotnum/dble(i)                        ! rot(2)= will store the max value of the frequency
          end if
          if (i .gt. int(dble(N)/3.d0)) then
          call minmax(rotnum/dble(i), rot)
          end if
          
          
          
           if ( fli(1) .gt. 12) then
           fli(1)=12.d0
*          goto  200
           end if
          


          if (i .eq. 5000) then
          fli5=fli(1)
          end if


          


         end do  ! over i=1,N
*************************************************************************
*************************************************************************


          rotnum=rotnum/dble(N)
          rotnu2=rotnu2/dble(N)

*************************************************************************
*************************************************************************
*       Associate a value for "findex"
*       We use 4 cutoff values to classify the types of motion
*       These values are inferred empirically

*************************************************************************
       fccha=7.d0        ! fccha= cutoff value used to classify chaotic orbits (if fli >= cochaos then the orbit is chaotic)
       fcrot1=3.d0
       fcrot2=3.9d0        ! fcrot1 and fcrot2 are the bounds for rotational motion (if corot1<= fli <= corot2 then the orbit is rotational or nonresonant)
       fclib=2.9d0         ! fclib =cutoff value used to classify the librational motions (if colib<= fli then the orbit is resonant)

       if (fli(1) .ge. fccha) then
       findex=0.d0
       else if ((fli(1) .gt. fcrot2) .and. (fli(1) .lt. fccha)) then
       findex=-1.d0
       else if ((isnres) .and. ((fli(1) .ge. fcrot1) .and.
     -  (fli(1) .le. fcrot2)))  then
       findex= 1.d0
       else ! if  (fli(1) .le. fclib) then
       findex=2.d0
       end if
       
       
        if (findex .lt. -0.5d0) then
         f2index=0.d0
         else
         f2index=findex
         end if

************************************************************************
************************************************************************
*************************************************************************
*      Associate a value for "bindex" using the ideas presented in Sander and Meiss
*      We use the following arguments. To check if the weighted Birkhoff average of the rotation number is poorly convergent (and to infer the chaotic nature of the orbit)
*      we compute and then compare the average over the interval [0,T] and then over [T, 2T] (here T=N/2).
*      To distinguish between rotational and librational motions we check the value of the rotation number. If this is an irrational number then the orbit is rotational, while
*      for rational numbers we are dealing with a librational orbit.

*      We use one cutoff value to distinguish poorly convergent sequences (chaotic orbits)
*      and two for distinguishing rotational and librational motions
************************************************************************

!            bccha=0.0115d0        ! accha= cutoff value used to classify chaotic orbits: if the difference between
             bccha=0.0105d0                    !  the maximum and minimum values of (rotnum/dble(i)), i=N/3, N/3+1, ...N is larger than accha then the orbit is chaotic
                                 
                                  ! delta1 and delta2 are some tolerance values the length of the intervals where we are looking for the


          if (dabs(rot(2)-rot(1)) .gt. bccha) then
          bindex=0.d0
          b2index=5.d0
          end if
          

          if ( bindex .lt. -0.5d0) then
          delta1=1.d-4                          ! 2.95d-3
          call smalde(brot, delta1,  num, denum)
          devom=dabs(dble(denum))
          if (devom .lt. 5.d0)  then                     ! rational number with the denumerator 1,2,3 or 4
          bindex=2.d0
          b2index=dble(num)/dble(denum)
          end if
          end if

          if (bindex .lt. -0.5d0) then
          delta2=1.d-6
          call smalde(brot, delta2,  num, denum)
          devom=dabs(dble(denum))
          if (devom .gt. 50.d0)  then                    ! irrational number
          bindex=1.d0
          b2index=6.d0
          end if
          end if
          
*********************************************************************
*************************************************************************
*************************************************************************
*      Associate a value for "aindex" using the value of rotation number, if it exists.
*      We use the following arguments.The limit of the sequence (rotnum/dble(i), see above) providing the rotatation number does not exist, or
*      the "finite" sequence (rotnum/dble(i)), i=1,2,...N is poorly convergent for chaotic orbits.
*      To distinguish between rotational and librational motions we use the monotonicity property of the rotation number.
*      More precisely, if we are inside a libration region, then the rotation number viewed as a function of the action is constant,
*      while outside the libration island, the rotation number increases if the value of the action increases

*      We use two cutoff values, one to distinguish poorly convergent sequences (chaotic orbits)
*      and one to distinguish constant from increasing rotation numbers (with respect to the action)


*************************************************************************
!        accha=0.0115d0        ! accha= cutoff value used to classify chaotic orbits: if the difference between
                                 !  the maximum and minimum values of (rotnum/dble(i)), i=N/3, N/3+1, ...N is larger than accha then the orbit is chaotic
!        acrot=0.039d0       ! acrot= cutoff value used to classify the rotational motions
        aclib=0.003d0
        a2clib=0.00005d0       ! aclib =cutoff value used to classify the librational motions (if the difference between the rotation number for the current value of the action
                            ! and the rotation number for the previous value of the action is less than aclib then the orbit is resonant)

!          if (dabs(rot(2)-rot(1)) .gt. accha) then
!          aindex=0.d0
!          end if


******************************
          if (ia .eq. 1) then
          outdat(1,1)=p0
          outdat(2,1)=rotnum
          outdat(3,1)= brot


          aindex=1.d0
          a2index=1.d0
          




        write(j100,568) q0/degree, p0,
     -   fli(1), outdat(3,1), findex, f2index,
     -   aindex, a2index, bindex, b2index
  568   format (' ', (F10.5), (F10.5), (F8.3), (F12.7), 6(F7.2), ' ')




          end if
*****************************

          if (ia .eq. 2) then
          outdat(1,2)=p0
          outdat(2,2)=rotnum
          outdat(3,2)= brot


          finxold=findex
          f2inxold=f2index
          binxold=bindex
          b2inxold=b2index
          fliold=fli(1)

          

          end if

************************************
           if ((ia .gt. 2) .and. (ia .lt. Npoint+1)) then
          outdat(1,3)=p0
          outdat(2,3)=rotnum
          outdat(3,3)= brot

          if ((outdat(2,1) .le. outdat(2,2)+aclib) .and.
     -      (outdat(2,2) .le. outdat(2,3)+aclib)) then

          if ( (dabs(outdat(2,2)- outdat(2,1)) .le. aclib) .and.
     -      (dabs(outdat(2,3)- outdat(2,2)) .le. aclib))   then
           aindex=2.d0
           else
           aindex=1.d0
         end if
         else
         aindex=0.d0
         end if
         
         
         
         if ((outdat(3,1) .le. outdat(3,2)+a2clib) .and.
     -      (outdat(3,2) .le. outdat(3,3)+a2clib)) then

          if ( (dabs(outdat(3,2)- outdat(3,1)) .le. a2clib) .and.
     -      (dabs(outdat(3,3)- outdat(3,2)) .le. a2clib))   then
           a2index=2.d0
           else
           a2index=1.d0
         end if
         else
         a2index=0.d0
         end if
         
         
         




        write(j100,548) q0/degree, outdat(1,2),
     -   fliold, outdat(3,2), finxold, f2inxold,  aindex, a2index,
     -    binxold, b2inxold
  548   format (' ', (F10.5), (F10.5), (F8.3), (F12.7), 6(F7.2), ' ')


          finxold=findex
          f2inxold=f2index
          binxold=bindex
          b2inxold=b2index
          fliold=fli(1)

          outdat(1,1)= outdat(1,2)
          outdat(2,1)= outdat(2,2)
          outdat(3,1)= outdat(3,2)

          outdat(1,2)= outdat(1,3)
          outdat(2,2)= outdat(2,3)
          outdat(3,2)= outdat(3,3)

          end if

*****************************
!          if (aindex .lt. -0.5d0) then
          if (ia .eq. Npoint+1) then
          outdat(1,3)=p0
          outdat(2,3)=rotnum
          outdat(3,3)=brot

            if ((outdat(2,1) .le. outdat(2,2)+aclib) .and.
     -      (outdat(2,2) .le. outdat(2,3)+aclib)) then

          if ( (dabs(outdat(2,2)- outdat(2,1)) .le. aclib) .and.
     -      (dabs(outdat(2,3)- outdat(2,2)) .le. aclib))   then
           aindex=2.d0
           else
           aindex=1.d0
         end if
         else
         aindex=0.d0
         end if
         
         
          if ((outdat(3,1) .le. outdat(3,2)+a2clib) .and.
     -      (outdat(3,2) .le. outdat(3,3)+a2clib)) then

          if ( (dabs(outdat(3,2)- outdat(3,1)) .le. a2clib) .and.
     -      (dabs(outdat(3,3)- outdat(3,2)) .le. a2clib))   then
           a2index=2.d0
           else
           a2index=1.d0
         end if
         else
         a2index=0.d0
         end if
         




        write(j100,578) q0/degree, outdat(1,2),
     -   fliold, outdat(3,2), finxold, f2inxold,
     -  aindex, a2index, binxold, b2inxold
  578   format (' ', (F10.5), (F10.5), (F8.3), (F12.7), 6(F7.2), ' ')



          finxold=findex
          f2inxold=f2index
          binxold=bindex
          b2inxold=b2index
          fliold=fli(1)

          outdat(1,1)= outdat(1,2)
          outdat(2,1)= outdat(2,2)
          outdat(3,1)= outdat(3,2)

          outdat(1,2)= outdat(1,3)
          outdat(2,2)= outdat(2,3)
          outdat(3,2)= outdat(3,3)
          
          
          
          
          aindex=1.d0
          a2index=1.d0
          
         if((aindex .gt. -0.5d0) .and. (aindex .lt. 0.5d0)) then
          write(j81,929) q0/degree,outdat(1,3), 5.d0
  929   format (' ', (F10.4), (F10.4), (F14.7) ' ')
          else
          write(j81,928) q0/degree,outdat(1,3), outdat(2,3)
  928   format (' ', (F10.4), (F10.4), (F14.7) ' ')
          end if



        write(j100,528) q0/degree, outdat(1,3),
     -   fli(1), outdat(3,3), finxold, f2inxold,  aindex, a2index,
     -    binxold, b2inxold
  528   format (' ', (F10.5), (F10.5), (F8.3), (F12.7), 6(F7.2), ' ')


          
          end if
*******************************************************************************
          





         if((bindex .gt. -0.5d0) .and. (bindex .lt. 0.5d0 )) then
          write(j811,997) q0/degree, p0, 5.d0, 5.d0
  997   format (' ', (F10.4), (F10.4), 2(F14.7) ' ')
          else
          write(j811,996) q0/degree, p0,  brot, rotnu2
  996   format (' ', (F10.4), (F10.4), 2(F14.7) ' ')
          end if




         write(j60,*) q0/degree, p0 ,  fli
         write(j80,*)
         

         
         

         
         

         end do
         write(j60,*)
         write(j811,*)
         end do
        end




*********************************************************************

*This subroutine implements the Runge-Kutta method of order 4 for the
*variational equations associated to our problem
********************************************************************

        subroutine rukugf(x, v1, v2, tsid, eps, h)
        implicit double precision(a-h,o-z)
        dimension x(2), v1(2), v2(2)
        dimension z(2)
        dimension cone(2), ctwo(2), cthree(2), cfour(2), ctemp(2)
        dimension bone(2), btwo(2), bthree(2), bfour(2), b2temp(2)
        dimension b3temp(2), b4temp(2) , adata(41)
         do ii=1,2
         z(ii)=x(ii)
         end do

        call fun(z, tsid, eps, bone)

          do ii=1,2
          b2temp(ii)=z(ii)+(h/2.d0)*bone(ii)
          end do
        call fun(b2temp ,  tsid+h/2.d0, eps,  btwo)
         
          do ii=1,2
          b3temp(ii)=z(ii)+(h/2.d0)*btwo(ii)
          end do
        call fun(b3temp , tsid+h/2.d0, eps,  bthree)

          do ii=1,2
          b4temp(ii)=z(ii)+h*bthree(ii)
          end do
         call fun(b4temp ,  tsid+h, eps,  bfour)
  
       do ii=1,2
       x(ii)=z(ii)+(h*(bone(ii)+2.d0*btwo(ii)+ 
     - 2.d0*bthree(ii)+bfour(ii)))/6.d0
       end do

       call fung(z, v1, tsid, eps, cone)
       do ii=1,2
       ctemp(ii)=v1(ii)+ h*(cone(ii)/2.d0)
       end do

       call fung(b2temp, ctemp , tsid+h/2.d0, eps, ctwo)
       do ii=1,2
       ctemp(ii)=v1(ii)+ h*(ctwo(ii)/2.d0)
       end do

       call fung(b3temp, ctemp , tsid+h/2.d0, eps,  cthree)
       do ii=1,2
       ctemp(ii)=v1(ii)+ h*cthree(ii)
       end do

       call fung(b4temp, ctemp ,  tsid+h, eps, cfour)
       do ii=1,2
       v1(ii)=v1(ii)+(h*(cone(ii)+2.d0*ctwo(ii)+
     - 2.d0*cthree(ii)+cfour(ii)))/6.d0
       end do


       call fung(z, v2, tsid, eps, cone)
       do ii=1,2
       ctemp(ii)=v2(ii)+ h*(cone(ii)/2.d0)
       end do

       call fung(b2temp, ctemp , tsid+h/2.d0, eps, ctwo)
       do ii=1,2
       ctemp(ii)=v2(ii)+ h*(ctwo(ii)/2.d0)
       end do

       call fung(b3temp, ctemp , tsid+h/2.d0, eps,  cthree)
       do ii=1,2
       ctemp(ii)=v2(ii)+ h*cthree(ii)
       end do

       call fung(b4temp, ctemp ,  tsid+h, eps, cfour)
       do ii=1,2
       v2(ii)=v2(ii)+(h*(cone(ii)+2.d0*ctwo(ii)+
     - 2.d0*cthree(ii)+cfour(ii)))/6.d0
       end do


       return
       end

********************************************************************
*     The function f of the equations of motion
*********************************************************************

        subroutine fun(x, th, eps, f)
	implicit double precision(a-h,o-z)
	dimension x(2), f(2)

		p=x(1)
		q=x(2)

        f(1)=  -(eps*sin(q)+eps*sin(q-th)+eps*sin(2.d0*q-3.d0*th));

       f(2)= p

        return
        end

*********************************************************************
*     The function g of the variational equations
********************************************************************

       subroutine fung(x, v,  th, eps, ga)
       implicit double precision(a-h, o-z)
       dimension x(2), v(2),  ga(2)

	  p=x(1)
	  q=x(2)

	  df1dp=0.d0


         df1dq=-(eps*cos(q)+eps*cos(q-th)+2.d0*eps*cos(2.d0*q-3.d0*th))



        df2dp= 1.d0

	df2dq=0.d0



       ga(1)=(df1dp*v(1)+df1dq*v(2))
       ga(2)=(df2dp*v(1)+df2dq*v(2))

       return
       end
       
**********************************************************************

       subroutine minmax(val, vec)
       implicit double precision(a-h, o-z)
       dimension vec(2)
       
       if(val .lt. vec(1)) then
       vec(1)=val
       end if
       if(val .gt. vec(2)) then
       vec(2)=val
       end if
       
       return
       end
       
*********************************************************************
*      The computation of FLI
**********************************************************************
       subroutine flicom(fli, v1, v2)
       implicit double precision(a-h, o-z)
       dimension v1(2), v2(2),  fli(1)

          do ifli=1,11
       	  pnorm=0.d0
       	  alpha=dble(ifli-1)/10.d0
	  do ii=1,2
	     pnorm=pnorm+(alpha*v1(ii)+dsqrt(1.d0-alpha**2)*v2(ii))**2
	  enddo
	  pnorm=0.5d0*dlog10(pnorm)
	  fli(1) = max(dabs(pnorm),fli(1))
	  enddo
	  
	  do ifli=1,11
       	  pnorm=0.d0
       	  alpha=dble(ifli-1)/10.d0
	  do ii=1,2
	     pnorm=pnorm+(alpha*v1(ii)-dsqrt(1.d0-alpha**2)*v2(ii))**2
	  enddo
	  pnorm=0.5d0*dlog10(pnorm)
	  fli(1) = max(dabs(pnorm),fli(1))
	  enddo
	  

       return
       end
       
**************************************************************************
*     The subroutine returns the closest rational number num/denum to xre, with the minimum denum (see Sander and Meiss for technical details)
**************************************************************************
       subroutine smalde(xre, delta,  num, denum)
       implicit double precision(a-h, o-z)
       integer num, denum, pl, pr, ql, qr

       sgn=1.d0
       if (xre .lt.0.d0) then
       sgn=-1.d0
       end if

       xreal=dabs(xre)

       num=0
       denum=1
       pl=0
       ql=1
       pr=1
       qr=0


       do ii=1,20000000
         if ( dabs(xreal-dble(num)/dble(denum)) .ge. delta) then
         num=pl+pr
         denum=ql+qr
         if (xreal .lt. dble(num)/dble(denum)) then
         pr=num
         qr=denum
         else
         pl=num
         ql=denum
         end if
         else
         goto 973
        end if
        end do
 973    continue

        num=sgn*num
       return
       end


**************************************************************************
*      The exponential bump function taken from Sander and Meiss (see that precede eq. (6).)
**************************************************************************
       function gsander(t)
       implicit double precision (a-h ,o-z)
       if ((t .gt. 0.d0) .and. (t .lt. 1.d0)) then
       gsander=exp(-1.d0/(t*(1.d0-t)))
       else
       gsander=0.d0
       end if
       return
       end function
************************************************************************
*      The function computes the quantity S provided by formula (7) in Sander and Meiss
****************************************************************************
       function ssander(N)
       implicit double precision (a-h ,o-z)
       ssander=0.d0

       do ii=1,N

       ssander=ssander+gsander(dble(ii)/(dble(N)))
       end do

       return
       end function
************************************************************************


