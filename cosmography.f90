!!! Module to cosmomc for estimation of cosmographic parameters 
!!! as described in arXiv:160X.XXXXX
!!! Alejandro Aviles (avilescervantes@gmail.com)
!!! 
!!! Declare the cosmographic (eis) parameters as usual. For example as type CMB: CMB%sf_E1,...
!!! Call the luminosity distance function: 
!!!                    dLcg(z,CMB%sf_E1,CMB%sf_E2,CMB%sf_E3,CMB%sf_H0) 
!!! in the supernovae.f90 file (or BAO.f90 or ...) instead of the usual angular distance function ( times (1+z)**2 ).
!!!
!!! You can get the statefinders (q_0,j_0,s_0) directly from CosmoMC by adding them as derived parameters
!!! to subroutine TP_CalcDerivedParams   (or CalcDerivedParams for old versions of CosmoMC )
!!! For example: 
!!!      q_0 =  -1  + CMB%sf_E1
!!!      derived(7)  = q0 
!!! An alternative is to perform the transformations directly to the chains when the fit is done.
!!!
!!! It is expected that the distributions are skewed and have long tails. Then, it is important 
!!! to set  MPI_Check_Limit_Converge = T in the .ini file 


module cosmography
Implicit none
real :: divcut, zlowcut,zhighcut
logical :: use_Eis, use_SC, use_hybrid1, use_hybrid2, force_LCDM

contains


subroutine cosmography_init
implicit none

    !Use just one of the following 4
	use_SC = .false.
	use_Eis = .false.
	use_hybrid1 = .false.
	use_hybrid2 = .true.
	
	
	force_LCDM =.true.  !This is to let vary only E1 and fix E2 and E3 using the degeneracies of LCDM 
	divcut = 0.4
	zlowcut = 0.05
	zhighcut = 0.9

    !Or better: set them in your ini file and add them through Ini_Read_* functions.


	If (use_Eis) write(*,*) 'Using Eis cosmography approach.'	
	If (use_hybrid1) write(*,*) 'Using hybrid1 cosmography approach.'	
	If (use_hybrid2) write(*,*) 'Using hybrid2 cosmography approach.'	
	If (use_SC) write(*,*) 'Using standard cosmography approach.'

	If (use_Eis .or. use_hybrid1) then
	write(*,*) 'zlow  = ', zlowcut, ', zmid  = ', divcut,  ', zhigh = ', zhighcut
	end if	 

	If (use_hybrid2) then
	write(*,*) 'zlow  = ', zlowcut, ', zmid  = ', divcut, ', zhigh=zmid'
	end if	 


end subroutine cosmography_init



function dLSC(z,E1,E2,E3)
implicit none
real dLSC
double precision z
real :: E1, E2, E3

dLSC = ( &
          z &
        + z**2*(2. - E1)/2. &
        + z**3*(-3.*E1 + 2.*E1**2 - E2 )/6. &
	    + z**4*(8.*E1**2 - 6*E1**3 - 4*E2 + 6.*E1*E2 - E3)/24. &
	   )

end function dLSC


function dLEis(z,E1,E2,E3)
implicit none
real dLEis, intint, znew 
double precision z
real :: E1, E2, E3
integer j

dLEis = 0
           if (z<zlowcut) then    !MODELO Eis
		    intint = z/100.
			Do j=1,100         ! The curves are quite smooth. Riemann integration is enough over the sampled domains.
				znew = znew +intint
				dLEis = dLEis + &
					   intint*(1. / (1. + E1*znew) )
	        end do 
			dLEis = (1+z)*dLEis            
			else if(z<divcut .and. z>= zlowcut) then       
		    intint = z/200.
			Do j=1,200
				znew = znew +intint
				dLEis = dLEis + &
					   intint*(1. / (1. + E1*znew + E2*znew**2/2.) )
	        end do 
			dLEis = (1+z)*dLEis 
			else if (z>=divcut .and. z<zhighcut) then
            intint = z/400.
			Do j=1,400
				znew = znew +intint
				dLEis = dLEis + &
					   intint*(1. / (1. + E1*znew + E2*znew**2/2. + E3*znew**3/6.) )
	        end do  
	        dLEis = (1+z)*dLEis
			else if (z>=zhighcut) then 
				dLEis = ( &                      
				  z     & 
				+ z**2* (- E1)/2. &
				+ z**3* (2*E1**2 - E2)/6.   & 
				+ z**4* (-6*E1**3 + 6*E1*E2 - E3) / 24.      &
!~ 				!+ z**5* (12*E1**4 - 18*E1**2*E2 + 3*E2**2 + 4*E1*E3) /60.  & !Maybe you are interested in expanding to higher order
				)   			
			dLEis = (1+z)*dLEis  
			end if

end function dLEis



function dLhybrid1(z,E1,E2,E3)
implicit none
real dLhybrid1
double precision z
real :: E1, E2, E3

        if (z<zlowcut) then     !SC upto q0
        dLhybrid1 = ( &
          z &
        + z**2*(2 - E1)/2. &    
        )  
        else if ((z >= zlowcut) .and. (z< divcut) ) then         !SC upto j0       
        dLhybrid1 = ( &
          z &
        + z**2*(2. - E1)/2. &
        + z**3*(-3. + 3.*(-1 + E1)**2 + 3.*E1 - E1**2 - E2)/6. &     
        )  
        else if ((z >= divcut) .and. (z< zhighcut) ) then
        dLhybrid1 = ( &
          z     & 
        + z**2* (- E1)/2. &
        + z**3* (2.*E1**2 - E2)/6.   & 
        + z**4* (-6.*E1**3 + 6.*E1*E2 - E3) / 24.      &
        )  
        dLhybrid1 = (1. + z)*dLhybrid1 
                
        else 
        dLhybrid1 = ( &
          z     & 
        + z**2* (- E1)/2. &
        + z**3* (2.*E1**2 - E2)/6.   & 
        + z**4* (-6.*E1**3 + 6.*E1*E2 - E3) / 24.      &
        + z**5* (12.*E1**4 - 18.*E1**2*E2 + 3.*E2**2 + 4.*E1*E3) /60.  &
        )  
        dLhybrid1 = (1. + z)*dLhybrid1 
        end if  

end function dLhybrid1



function dLhybrid2(z,E1,E2,E3)
implicit none
real dLhybrid2
double precision z
real :: E1, E2, E3

        if (z<zlowcut) then     !SC upto q0
        dLhybrid2 = ( &
          z &
        + z**2*(2 - E1)/2. &    
        )  
        else if ((z >= zlowcut) .and. (z< divcut) ) then         !SC upto j0       
        dLhybrid2 = ( &
          z &
        + z**2*(2. - E1)/2. &
        + z**3*(-3. + 3.*(-1 + E1)**2 + 3.*E1 - E1**2 - E2)/6. &     
        )  
        else 
        dLhybrid2 = ( &
          z     & 
        + z**2* (- E1)/2. &
        + z**3* (2.*E1**2 - E2)/6.   & 
        + z**4* (-6.*E1**3 + 6.*E1*E2 - E3) / 24.      &
        + z**5* (12.*E1**4 - 18.*E1**2*E2 + 3.*E2**2 + 4.*E1*E3) /60.  &
        )  
        dLhybrid2 = (1. + z)*dLhybrid2 
        end if  

end function dLhybrid2



function dLcg(z,E1,E2,E3,H0)
implicit none
real dLcg
double precision, intent(in) :: z
real, parameter :: clight = 2.99792458e5
logical, save :: do_cosmography_init = .true.
real, intent(in) :: E1, E2, E3, H0
real :: pE1, pE2, pE3

	if (do_cosmography_init) then 
		call cosmography_init
		do_cosmography_init = .false.
	end if
		
	if (force_LCDM) then 
		if (use_Eis) then
			dLcg = dLEis(z,E1,-E1*(E1-2.),- E1 +  3.*(E1-1.)**2 + 3.*(E1-1.)**3) / H0 * clight
		else if (use_SC) then
			dLcg = dLSC(z,E1,-E1*(E1-2.),- E1 +  3.*(E1-1.)**2 + 3.*(E1-1.)**3) / H0 * clight
		else if (use_hybrid1) then
			dLcg = dLhybrid1(z,E1,-E1*(E1-2.),- E1 +  3.*(E1-1.)**2 + 3.*(E1-1.)**3) / H0 * clight	
		else if (use_hybrid2) then
			dLcg = dLhybrid2(z,E1,-E1*(E1-2.),- E1 +  3.*(E1-1.)**2 + 3.*(E1-1.)**3) / H0 * clight				
		end if
	else 
		if (use_Eis) then
			dLcg = dLEis(z,E1,E2,E3) / H0 * clight
		else if (use_SC) then
			dLcg = dLSC(z,E1,E2,E3) / H0 * clight
		else if (use_hybrid1) then
			dLcg = dLhybrid1(z,E1,E2,E3) / H0 * clight	
		else if (use_hybrid2) then
			dLcg = dLhybrid2(z,E1,E2,E3) / H0 * clight				
		end if
	end if


end function dLcg


end module
