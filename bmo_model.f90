  module bmo_model



    use parameters

    use parameters_bmo



    implicit none



    character (len=100) :: out_file_bmo

    double precision :: rho_cen_b, Pc_b, cp_m

    double precision :: k0_bm, k1_bm, k2_bm

    double precision :: Tm0_m, Tm1_m, Tm2_m, Tm3_m, Tm4_m

    double precision :: M_m, M_om, M_bm

    double precision :: rho0_om, rho1_om, rho2_om, rho3_om, rho0_bm, rho2_bm

    double precision :: T_cen_b,  t1_m,  t2_m,  t3_m

    double precision :: ds0_m, ds1_m, ds2_m, ds3_m, ds4_m

    double precision :: A, A2, A3, B2, D, D2, D4, L, L2, L3

    double precision :: pr_m(0:Nm)   , rhor_m(0:Nm) , gr_m(0:Nm) , psir_m(0:Nm)

    double precision :: Tar_m(0:Nm)  , Tm_c(0:Nm)   , ds_m(0:Nm)

    double precision :: dTadr_m(0:Nm), dTmdP_m(0:Nm), kr_m(0:Nm)

    double precision ::  r_m(0:Nm),r2_m(0:Nm),r3_m(0:Nm),r4_m(0:Nm),r5_m(0:Nm)

    double precision :: r6_m(0:Nm),r7_m(0:Nm),r8_m(0:Nm),r9_m(0:Nm),r10_m(0:Nm),r11_m(0:Nm)

    double precision :: rb0, rr, h0_m, h

    double precision :: dt_b, ttot_b

    double precision, allocatable :: xflux(:), time_b_file(:), dt_b_file(:)

    double precision, allocatable :: rb_b_file(:), Qr_b_file(:)

    double precision, allocatable :: dridt_b_file(:),dr_b_file(:)

    integer          :: ntb, out_file_bmo_len, sol, davies_or_nimmo, ah_b, hor_b

    integer          :: iteration_b



  contains



!******************************************************************************

! read_input: reads the input file from command

!

! Inputs: Qx  - CMB heat flux

!         Tx - CMb temperature

!         rb - Inner core radius

!         EJ_b - Ohmic heating

!

! Comments in the input file are lines starting with an *

! If the solution method sol=0 then Qx is read in from a file

!******************************************************************************

  subroutine read_input_bmo(Qx,Tx,rb,EJ_b,cbar_elm1,cbar_elm2,cbar_elm3)



    double precision      :: cbar_elm1, cbar_elm2, cbar_elm3

    double precision      :: tmp, Tx, rb, Qx, EJ_b

    character (len=100)   :: line



21    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 21

      out_file_bmo     = adjustl(line)

      out_file_bmo_len = index(out_file_bmo, ' ')-1

      out_file_bmo    = out_file_bmo(1:out_file_bmo_len)



22    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 22

      read(line, *) davies_or_nimmo, ah_b, hor_b



23    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 23

      read(line, *) rho_cen_b, Pc_b, cp_m



24    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 24

      read(line, *) rho0_om, rho1_om, rho2_om, rho3_om, rho0_bm, rho2_bm



25    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 25

      read(line, *) ds0_m, ds1_m, ds2_m, ds3_m, ds4_m



26    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 26

      read(line, *) Tm0_m, Tm1_m, Tm2_m, Tm3_m, Tm4_m



27    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 27

      read(line, *) T_cen_b, t1_m, t2_m, t3_m



28    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 28

      read(line, *) k0_bm, k1_bm, k2_bm



29    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 29

      read(line, *) cbar_elm1, cbar_elm2, cbar_elm3



30    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 30

      read(line, *) Tx, rb, rr, h0_m

      h = h0_m



31    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 31

      read(line, *) dt_b,ttot_b



      Qx  = 0.d0

!      Ej = 0.d0

      EJ_b = 0.d0

32    continue

      read(5, 80) line

      if(line(1:1) .eq. '*') goto 32

      read(line, *) sol, tmp



      if(sol==1 .or. sol==3 .or. sol==4) then

        Qx  = tmp

      elseif(sol==2) then

        EJ_b = tmp

      endif



!     If read in Qx, ntb to be determined

      ntb = 0

      if(sol==0) then



33      continue

        read(5, 80) line

        if(line(1:1) .eq. '*') goto 33



        call read_layerhf_bmo(line)

      else

        ntb = ttot_b/dt_b



        allocate(xflux(0:ntb))



        xflux = 0.d0

      endif



      write(*,*) '# timepts = ', ntb



      close(5)



 80   FORMAT(A)



  end subroutine read_input_bmo



!******************************************************************************

! read_layerhf_bmo: reads the heat flux Qx out of BMO from a file

!

! Inputs: Qcmb_file - name of file containing Qx

!                     format is 2 columns: t, Qx

! LZ 6/5/2015 - file format is t, radius to top of bmo, Qx, Qr

!   **just to reproduce labrosse, will eventually calculate Qr directly

!

! n is a parameter set arbitrarily large (hopefully larger than the # rows!)

!******************************************************************************

  subroutine read_layerhf_bmo(Qcmb_file)



    integer         , parameter   :: n=1e6

    double precision              :: tmp(n), t(n)

    double precision              :: tmp2(n), tmp3(n)

    character (len=100) :: Qcmb_file

    integer :: i



    open(99,file=Qcmb_file)



    t  = 0.d0; tmp = 0.d0; tmp2 = 0.d0; tmp3 = 0.d0

    ntb = 0



    do i = 0, n

      read(99,*,end=100) t(i), tmp2(i), tmp(i), tmp3(i)

      ntb = ntb + 1

    enddo



100 continue



    ntb = ntb - 1



    allocate(xflux(0:ntb))

    allocate(time_b_file(0:ntb))

    allocate(dt_b_file(0:ntb))



    allocate(Qr_b_file(0:ntb))            !force radiogen Q and layer height from file

    allocate(rb_b_file(0:ntb))



    allocate(dr_b_file(0:ntb))

    allocate(dridt_b_file(0:ntb))



!   Set up time and Qx so they run backwards.

    dt_b_file = 0.d0

    time_b_file = 0.d0

    xflux = 0.d0



    Qr_b_file = 0.d0

    rb_b_file = 0.d0

    dr_b_file = 0.d0

    dridt_b_file = 0.d0



    do i = 0, ntb

      time_b_file(i) = t(ntb-i)*1e3

      xflux(i)   = tmp(ntb-i)*1e12                                  !Qx in TW

      Qr_b_file(i)   = tmp3(ntb-i)*1e12                                  !Qx in TW

      rb_b_file(i)  =tmp2(ntb-i)			!radius to top of bmo in meters

      if(i > 0) dt_b_file(i) = dabs(time_b_file(i) - time_b_file(i-1))*1e6 !dt_b in yr

      if(i > 0) dr_b_file(i) = dabs(rb_b_file(i) - rb_b_file(i-1)) !dr_b in meters

      if(i > 0) dridt_b_file(i) = dr_b_file(i)/(secinyr*dt_b_file(i)) !rough dridt in m/yr

!      print *, time_b_file(i), dt_b_file(i-1), xflux(i), Qr_b_file(i), rb_b_file(i), &

!      dr_b_file(i-1), dridt_b_file(i-1)

!      write(546,*) time_b_file(i), dt_b_file(i-1), xflux(i), Qr_b_file(i)

      write(546,*) time_b_file(i), dt_b_file(i-1), xflux(i), Qr_b_file(i), &

         rb_b_file(i), dr_b_file(i-1), dridt_b_file(i-1)

    enddo

    dt_b_file(0) = dt_b_file(1)

    dridt_b_file(0) = dridt_b_file(1)



  end subroutine read_layerhf_bmo



!******************************************************************************

! rpts: sets up the radial grid

!

! Inputs: rb - dimensional radius of the BMO

!         rr - formerly dimensional radius of the CMB, now meant as

!	       reference radius, maximum BMO gets/could be

!

! Outputs: r - array of radius points (global variable)

!

! Nb points are always allocated to the inner core; Nr to outer core

!******************************************************************************

  subroutine rpts_bmo(rb,rt)


! DRS 4/25/18
! in Leah's version (v1) this grided pts from cmb (rc) to moving interface (rb) and the second argument that was input
! was rr and was not used within the subroutine but passed in as a remnent from it's previous usage

    double precision, intent(in) :: rb, rt

    double precision :: dric, droc, drs      !dr, ic, oc, strat layer

    integer          :: i,bottom,top,oclower !where to start and stop in core



        r_m(0) = rc                              !First point is cmb

        r2_m(0) = rc*rc

        r3_m(0) = rc*rc*rc



        if(rb .ne. rc) then                      !Check there is an inner core

          dric =     (rb-rc) /real(Nb)           !Divide the points equally in space

          do i = 1, Nb

            r_m(i)  = r_m(i-1) + dric

            r2_m(i) = r_m(i)*r_m(i)

            r3_m(i) = r_m(i)*r_m(i)*r_m(i)

            r4_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)

            r5_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

            r6_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

            r7_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

            r8_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

            r9_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

            r10_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

            r11_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

          enddo

          bottom = Nb+1

        else                                     !If not, OC starts from the centre

          bottom = 1

        endif



        if(rb .ne. rc) then                    !Check there is an inner core

          droc = (rt-rb)/real(Nm-bottom)

          r_m(bottom) = r_m(Nb)

          r2_m(bottom) = r_m(Nb)*r_m(Nb)

          r3_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)

          r4_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)

          r5_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)

          r6_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)

          r7_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)

          r8_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)

          r9_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)

          r10_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)

          r11_m(bottom) = r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)*r_m(Nb)



          oclower = bottom + 1

        else

          droc = (rt-rb)/real(Nm-bottom+1)

          oclower = bottom

        endif



        do i = oclower, Nm                       !Convecting core

          r_m(i)  = r_m(i-1) + droc

          r2_m(i) = r_m(i)*r_m(i)

          r3_m(i) = r_m(i)*r_m(i)*r_m(i)

          r4_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)

          r5_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

          r6_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

          r7_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

          r8_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

          r9_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

          r10_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

          r11_m(i) = r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)*r_m(i)

        enddo



  end subroutine rpts_bmo



!******************************************************************************

! mole2massconc: convert the inputted molar concentrations of O, S and Si into

!                mass concentrations that are evolved in time.

!

! Inputs: Molar concentrations (bars; global variables)

! Outputs: Mass concentrations c_elm1, c_elm2, c_elm3

!

! Mole fraction (#atoms of solute/total #atoms) c' = A'/AL * c, c = mass fraction

!

! Only O is evolved in time because the core chemistry model used here assumes

! that only O is responsible for the compositional part of the ICB density jump

!

! Molecular weights are set as parameters in parameters.f90

!******************************************************************************

  subroutine mole2massconc_bmo(cbar_elm1, cbar_elm2, cbar_elm3, c_elm1, c_elm2, c_elm3)



    double precision, intent(in)  :: cbar_elm1, cbar_elm2, cbar_elm3

    double precision, intent(out) :: c_elm1, c_elm2, c_elm3

    double precision              :: Abar



    Abar = cbar_elm1*Aelm1 + cbar_elm2*Aelm2 + cbar_elm3*Aelm3 + (1.d0-cbar_elm1-cbar_elm2-cbar_elm3)*AXx



    c_elm1  = cbar_elm1 *((Aelm1*1.d0)/Abar)

    c_elm2  = cbar_elm2 *((Aelm2*1.d0)/Abar)

    c_elm3 = cbar_elm3*((Aelm3*1.d0)/Abar)



    write(*,'(A16,E16.6)') 'c_elm1=',c_elm1

    write(*,'(A16,E16.6)') 'c_elm2=',c_elm2

    write(*,'(A16,E16.6)') 'c_elm3=',c_elm3



  end subroutine mole2massconc_bmo



  subroutine mass2moleconc_bmo(c_elm1, c_elm2, c_elm3, cbar_elm1, cbar_elm2, cbar_elm3)



    double precision, intent(in)    :: c_elm1, c_elm2, c_elm3

    double precision, intent(inout) :: cbar_elm1, cbar_elm2, cbar_elm3

    double precision                :: Abar



    Abar = 1.d0 / (c_elm1/Aelm1 + c_elm2/Aelm2 + c_elm3/Aelm3 + (1.d0-c_elm1-c_elm2-c_elm3)/AXx)



    cbar_elm1  = c_elm1 / ((Aelm1*1.d0)/Abar)

    cbar_elm2  = c_elm2 / ((Aelm2*1.d0)/Abar)

    cbar_elm3 = c_elm3 /((Aelm3*1.d0)/Abar)



    write(*,'(A16,E16.6)') 'cbar_elm1=',cbar_elm1

    write(*,'(A16,E16.6)') 'cbar_elm2=',cbar_elm2

    write(*,'(A16,E16.6)') 'cbar_elm3=',cbar_elm3



  end subroutine mass2moleconc_bmo



!******************************************************************************

! adiabat_poly: Calculates the adiabatic temperature using equation (14) of

!               Davies (PEPI, 2014, submitted)

!

! Inputs: Nm - index for CMB

!         Tx - Temperature at CMB

!

! Outputs: Qa - Adiabatic heat flux

!          Tar_m - Adiabatic temperature profile (global variable)

!          dTadr_m - Adiabatic temperature gradient (global variable)

!

! The correction procedure for Ta is the same as described in adiabat()

!******************************************************************************

  subroutine adiabat_poly_bmo(Nm, Tx, Qa)



    integer         , intent(in) :: Nm

    double precision, intent(in) :: Tx

    double precision :: Qa, dTadp_m(0:Nm)

    double precision :: delta

    integer:: i



    T_cen_b  = Tx /      (1.d0 + t1_m*(pr_m(Nm)/1e9) + t2_m*(pr_m(Nm)*pr_m(Nm))/1e18 + t3_m*(pr_m(Nm)*pr_m(Nm)*pr_m(Nm)/1e27))

    Tar_m    = T_cen_b * (1.d0 + t1_m*(pr_m    /1e9) + t2_m*(pr_m*pr_m)        /1e18 + t3_m*(pr_m*pr_m*pr_m            /1e27))

    dTadp_m  = T_cen_b * (       t1_m                + 2.d0*t2_m*pr_m          /1e9  + 3.d0*t3_m*(pr_m*pr_m            /1e18))

    dTadP_m  = dTadP_m / 1e9


    dTadr_m  = -dTadp_m*rhor_m*gr_m



!    T_cen_b  = Tx / (1.d0 + t1_m*(r_m(Nm)/1e3) + t2_m*(r2_m(Nm)/1e6) + t3_m*(r3_m(Nm)/1e9))

!    Tar_m  = T_cen_b * (1.d0 + t1_m*(r_m/1e3) + t2_m*(r2_m/1e6)    + t3_m*(r3_m/1e9))



!    Tar_m   = T_cen_b * (1.d0 + t1_m*(r_m/1e3) + t2_m*(r2_m/1e6) + t3_m*(r3_m/1e9))

!    delta = Tx - Tar_m(Nm)              !Compare Ta @ CMB to Tx



!    T_cen_b  = T_cen_b + delta

!    Tar_m   = T_cen_b * (1.d0 + t1_m*(r_m/1e3) + t2_m*(r2_m/1e6) + t3_m*(r3_m/1e9))



!    dTadr_m = T_cen_b * (t1_m + 2.d0*t2_m*(r_m/1e3) + 3.d0*t3_m*(r2_m/1e6))

!    dTadr_m = dTadr_m / 1e3



    Qa    = -4.d0*pi*r2_m(Nm)*kr_m(Nm)*dTadr_m(Nm)



  end subroutine adiabat_poly_bmo



!******************************************************************************

! density_poly: Calculate core density using equation (7) of Davies (2014, PEPI)

!

! rhor_m(Nb)*conc*alphac_c is the part of the density jump due to Oxygen.

! Presently this is not used because the density jump input to the code does not

! have to equal that in the density profile.

!******************************************************************************

  subroutine density_poly_bmo(Nm, conc, alphac)



    integer, intent(in)          :: Nm

    double precision, intent(in) :: conc, alphac

    double precision             :: drho_comp



    rhor_m       = 0.d0

    rhor_m(0:Nm) = rho0_om + rho1_om*(r_m(0:Nm)/1e3) + rho2_om*(r2_m(0:Nm)/1e6) + rho3_om*(r3_m(0:Nm)/1e9)

    drho_comp = rhor_m(Nb)*conc*alphac



  end subroutine density_poly_bmo



  !******************************************************************************

  ! gravity: Calculates gravity numerically

  !******************************************************************************

    subroutine gravity_num_bmo(Nm)



      double precision, parameter :: M_c=0.15499e24

      integer, intent(in) :: Nm

      double precision  :: integrand(0:Nm)

      integer           :: i



      do i = 0, Nm

        integrand(i) = rhor_m(i)*r_m(i)*r_m(i)

      enddo

      call splat(integrand,r_m,Nm+1)



      gr_m = 4.d0*pi*G* (integrand + M_c) / r2_m



    end subroutine gravity_num_bmo



!******************************************************************************

! gravity_poly: Calculates gravity using equations (10) and (11) of Davies

!               (2014, PEPI, submitted)

!******************************************************************************

  subroutine gravity_poly_bmo(Nb)



    integer, intent(in) :: Nb

    gr_m =  10.0000



  end subroutine gravity_poly_bmo



  subroutine grav_pot_num_bmo(Nm)



      integer :: i, Nm

      double precision  :: integrand(0:Nm)



      do i = 0, Nm

        integrand(i)=-gr_m(i)

      enddo

      call splat(integrand,r_m,Nm+1)

      do i = 0, Nm

        psir_m(i)=integrand(Nm-i)

      enddo



    end subroutine grav_pot_num_bmo



!******************************************************************************

! grav_pot_poly: Calculates the gravitational potential using equation (12) of

!                Davies (2014, PEPI, submitted)

!

! Inputs: Nb - Index for ICB

!         Nm - Index for CMB

!******************************************************************************

  subroutine grav_pot_poly_bmo(Nb, Nm)



    integer, intent(in) :: Nb, Nm

    double precision    :: gpot_cmb, gpot_icb, gpot_icbi

    integer :: i



    psir_m = 0.d0; gpot_cmb=0.d0; gpot_icb=0.d0; gpot_icbi=0.d0



    gpot_cmb = rho0_om*(r2_m(Nm)/1e6)/6.d0    + rho1_om*(r3_m(Nm)/1e9)/12.d0 + &

               rho2_om*(r4_m(Nm)/1e12)/20.d0  + rho3_om*(r5_m(Nm)/1e15)/30.d0



    gpot_icb = rho0_om*(r2_m(Nb+1)/1e6)/6.d0   + rho1_om*(r3_m(Nb+1)/1e9)/12.d0 + &

               rho2_om*(r4_m(Nb+1)/1e12)/20.d0 + rho3_om*(r5_m(Nb+1)/1e15)/30.d0



    gpot_icbi= rho0_bm*(r2_m(Nb)/1e6)/6.d0     + &

               rho2_bm*(r4_m(Nb)/1e12)/20.d0



    do i = 0, Nb

      psir_m(i) = rho0_bm*(r2_m(i)/1e6)/6.d0 + rho2_bm*(r4_m(i)/1e12)/20.d0 - gpot_cmb + gpot_icb - gpot_icbi

    enddo



    do i = Nb+1, Nm

      psir_m(i) = rho0_om*(r2_m(i)/1e6)/6.d0   + rho1_om*(r3_m(i)/1e9)/12.d0 + &

                rho2_om*(r4_m(i)/1e12)/20.d0 + rho3_om*(r5_m(i)/1e15)/30.d0 - gpot_cmb !-gpot_icb + &

!                rho0_bm*(r2_m(Nb)/1e6)/6.d0 + rho2_bm*(r4_m(Nb)/1e12)/20.d0

    enddo



    psir_m = 4.d0 * pi * G * psir_m * 1e6



  end subroutine grav_pot_poly_bmo



  subroutine pressure_num_bmo(Nm,Pc_b)



      integer :: i, Nm

      double precision  :: integrand(0:Nm), Pc_b



      do i = 0, Nm

        integrand(i)=-rhor_m(i)*gr_m(i)

      enddo

      call splat(integrand,r_m,Nm+1)



      pr_m = integrand + Pc_b



    end subroutine pressure_num_bmo



!******************************************************************************

! pressure_poly: Calculates pressure using equation (13) of Davies (2014, PEPI)

!

! Inputs: Nb - Index for ICB

!         Nm  - Index for CMB

!         Pc_b - Pressure at centre of Earth

!******************************************************************************

  subroutine pressure_poly_bmo(Nb, Nm, Pc_b)



    double precision, intent(in)  :: Pc_b

    integer         , intent(in)  :: Nm, Nb

    double precision              :: one_coeff, two_coeff, three_coeff, four_coeff, five_coeff, six_coeff, seven_coeff

    double precision              :: P_icb, P_cmb, P_icbi



    pr_m = 0.d0



    one_coeff    =       rho0_om*rho0_om/6.d0

    two_coeff    = 7.d0 *rho0_om*rho1_om/36.d0

    three_coeff  = 2.d0 *rho0_om*rho2_om/15.d0  +      rho1_om*rho1_om/16.d0

    four_coeff   =       rho0_om*rho3_om/10.d0  + 9.d0*rho1_om*rho2_om/100.d0

    five_coeff   = 5.d0 *rho1_om*rho3_om/72.d0  +      rho2_om*rho2_om/30.d0

    six_coeff    = 11.d0*rho2_om*rho3_om/210.d0

    seven_coeff  =       rho3_om*rho3_om/42.d0



    P_icb        =  4.d6 * pi * G * (one_coeff  *(r2_m(Nb+1)/1e6)  + two_coeff *(r3_m(Nb+1)/1e9)  + &

                   three_coeff*(r4_m(Nb+1)/1e12) + four_coeff *(r5_m(Nb+1)/1e15) + five_coeff*(r6_m(Nb+1)/1e18) + &

                   six_coeff  *(r7_m(Nb+1)/1e21) + seven_coeff*(r8_m(Nb+1)/1e24))



    P_cmb        =  4.d6 * pi * G * (one_coeff  *(r2_m(Nm)/1e6)      + two_coeff *(r3_m(Nm)/1e9) +   &

                   three_coeff*(r4_m(Nm)/1e12) + &

                   four_coeff *(r5_m(Nm)/1e15)    + five_coeff*(r6_m(Nm)/1e18)    + six_coeff  *(r7_m(Nm)/1e21)  + &

                   seven_coeff*(r8_m(Nm)/1e24))

    if(Nb > 0)  then

!  LZ this is off by a minus sign from the paper/notes, but makes more sense here

!  Also, Pc_b is NOT the pressure at center of earth as infile says, it's pressure at CMB from overlying mantle

      pr_m(Nb+1:Nm)   = Pc_b + P_cmb -4.d6 * pi * G * &

                     (one_coeff  *(r2_m(Nb+1:Nm)/1e6)  + two_coeff *(r3_m(Nb+1:Nm)/1e9)  + three_coeff*(r4_m(Nb+1:Nm)/1e12) + &

                     four_coeff *(r5_m(Nb+1:Nm)/1e15) + five_coeff*(r6_m(Nb+1:Nm)/1e18) + six_coeff  *(r7_m(Nb+1:Nm)/1e21) + &

                     seven_coeff*(r8_m(Nb+1:Nm)/1e24))



      one_coeff    =      rho0_bm*rho0_bm/6.d0

      two_coeff    = 8.d0*rho0_bm*rho2_bm/60.d0

      three_coeff  =      rho2_bm*rho2_bm/30.d0



      P_icbi       = 4.d6 * pi * G * (one_coeff*(r2_m(Nb)/1e6) + two_coeff*(r4_m(Nb)/1e12) + three_coeff*(r6_m(Nb)/1e18))



      pr_m(0:Nb)     = Pc_b + P_cmb - P_icb + P_icbi -4.d6 * pi * G * (one_coeff*(r2_m(0:Nb)/1e6) + two_coeff*(r4_m(0:Nb)/1e12) + &

                                                               three_coeff*(r6_m(0:Nb)/1e18))



    else

!  LZ this is off by a minus sign from the paper/notes, but makes more sense here

!  Also, Pc_b is NOT the pressure at center of earth as infile says, it's pressure at CMB from overlying mantle

!      pr_m(0:Nm)   = Pc_b + P_cmb -4.d6 * pi * G * &

!                     (one_coeff  *(r2_m(0:Nm)/1e6)  + two_coeff *(r3_m(0:Nm)/1e9)  + three_coeff*(r4_m(0:Nm)/1e12) + &

!                     four_coeff *(r5_m(0:Nm)/1e15) + five_coeff*(r6_m(0:Nm)/1e18) + six_coeff  *(r7_m(0:Nm)/1e21) + &

!                     seven_coeff*(r8_m(0:Nm)/1e24))

!  For BMO:  always only considering the liquid layer.

!  since density of melt off (assuming constant to start), reformulate as P = Pcmb (139GPa) - Pressure due to

!  material between rcmb and r...

!  Picb is actually P evaluated at rc, cmb radius

!  P_r = Pc_b (pressure at CMB present day from overlying mantle) -P(r) + P_icb (which is P(r=rcmb))

!      pr_m(0:Nm)   = Pc_b + P_icb -4.d6 * pi * G * &

!                     (one_coeff  *(r2_m(0:Nm)/1e6)  + two_coeff *(r3_m(0:Nm)/1e9)  + three_coeff*(r4_m(0:Nm)/1e12) + &

!                     four_coeff *(r5_m(0:Nm)/1e15) + five_coeff*(r6_m(0:Nm)/1e18) + six_coeff  *(r7_m(0:Nm)/1e21) + &

!                     seven_coeff*(r8_m(0:Nm)/1e24))

! slope in PREM is 54700 Pa/m

   pr_m(0:Nm) = Pc_b + 54700.d0 * (rc - r_m(0:Nm))

    endif



  end subroutine pressure_poly_bmo



!******************************************************************************

! mass: Calculates the mass of the inner core (Mic) and outer core (M_om) using

!       equations (8) and (9) of Davies (2014, PEPI, submitted)

!

! Inputs: Nb - Index for ICB

!         Nr - Index for CMB

!******************************************************************************

  subroutine mass_poly_bmo(Nb,Nr)



    integer, intent(in) :: Nr,Nb

    double precision    :: M_om_cmb, M_om_icb, Mic_icb



    M_om_cmb = rho0_om*(r3_m(Nm) /1e9)/3.d0    + rho1_om*(r4_m(Nm) /1e12)/4.d0   + &

               rho2_om*(r5_m(Nm) /1e15)/5.d0   + rho3_om*(r6_m(Nm) /1e18)/6.d0

    M_om_icb = rho0_om*(r3_m(Nb+1)/1e9)/3.d0  + rho1_om*(r4_m(Nb+1)/1e12)/4.d0 + &

               rho2_om*(r5_m(Nb+1)/1e15)/5.d0 + rho3_om*(r6_m(Nb+1)/1e18)/6.d0

    M_om     = 4.d0 * pi * (M_om_cmb - M_om_icb) * 1e9



    Mic_icb = rho0_bm*(r3_m(Nb)/1e9)/3.d0 + rho2_bm*(r5_m(Nb)/1e15)/5.d0

    M_bm    = 4.d0 * pi * (Mic_icb) * 1e9



    M_m = M_om + M_bm



    write(*,'(A16,E16.9)') 'M_om=',M_om

    write(*,'(A16,E16.9)') 'M_m =',M_m



  end subroutine mass_poly_bmo



!******************************************************************************

! mass_poly_correct_bmo: Correct for changes in the core mass over time.

!                 Assume mass of IC is constant.

!

! Inputs: Nb - Index for ICB

!         Nr - Index for CMB

!******************************************************************************



  subroutine mass_poly_correct_bmo(rb)



    double precision, intent(in) :: rb

    double precision             :: rho_new, fac, ic, oc1, oc2



    if(rb == rb0) return



    fac =  3.d0 / ( (rc*rc*rc/1e9) - (rb*rb*rb/1e9) )



    oc1 = rho0_om*( (rc*rc*rc/1e9) - (rb0*rb0*rb0)/1e9)/3.d0



    oc2 = rho1_om*( (rb0*rb0*rb0*rb0/1e12)         - (rb*rb*rb*rb/1e12) )/4.d0 + &

          rho2_om*( (rb0*rb0*rb0*rb0*rb0/1e15)     - (rb*rb*rb*rb*rb)/1e15 )/5.d0 + &

          rho3_om*( (rb0*rb0*rb0*rb0*rb0*rb0/1e18) - (rb*rb*rb*rb*rb*rb)/1e18 )/6.d0



    ic  = rho0_bm*( (rb0*rb0*rb0/1e9)          - (rb*rb*rb)/1e9  )/3.d0 + &

          rho2_bm*( (rb0*rb0*rb0*rb0*rb0/1e15) - (rb*rb*rb*rb*rb)/1e15 )/5.d0



    rho_new = (oc1 - oc2 + ic) * fac



    print *, ic,oc1, oc2, fac

    print *, 'rho_new, old = ', rho_new, rho0_om



    rho0_om = rho_new



  end subroutine mass_poly_correct_bmo



!

!******************************************************************************

! melting_poly: The melting curve using equation (16) of Davies (2014, submitted

!               to PEPI)

!******************************************************************************

  subroutine entropy_melting_bmo()



    ds_m = 0.d0



    ds_m    = ds0_m + ds1_m*(pr_m/1e9) + ds2_m*(pr_m*pr_m)/1e18 + ds3_m*(pr_m*pr_m*pr_m)/1e27 + ds4_m*(pr_m*pr_m*pr_m*pr_m)/1e36



  end subroutine entropy_melting_bmo



!******************************************************************************

! melting_poly: The melting curve using equation (16) of Davies (2014, submitted

!               to PEPI)

!******************************************************************************





  subroutine melting_poly_bmo()

    Tm_c = 0.d0; dTmdP_m = 0.d0

! LZ June 30, 2015:  temporarily commenting out the polyfit and using instead
! an ad hoc log fit to the zigsteg2013 'melting curve' made from labrosse
! methods.
!  The log fit was done in matlab after converting radius of zigsteg2013
!  bmo layer height into Pressure using PREM. Pressure calculated here may
!  deviate slightly from that so that is a minor inconsistancy
!
!DS 04/26/18 3 lines below are for polynomial approach with input params of 1967.23038731663 0.0229434063905606 -7.45249905371733e-05 0.0 0.0
!DS cont. and gives a melting curve for ~5300K at CMB to ~ 4700K at CMB+1000km depth
!    Tm_c    = Tm0_m*(1.d0+Tm1_m*(pr_m/1e9)+Tm2_m*(pr_m*pr_m)/1e18+Tm3_m*(pr_m*pr_m*pr_m)/1e27 + Tm4_m*(pr_m*pr_m*pr_m*pr_m)/1e36)
!    dTmdP_m = Tm0_m * (       Tm1_m       + 2.d0*Tm2_m*pr_m/1e9    + 3.d0*Tm3_m*pr_m*pr_m/1e18 + 4.d0*Tm4_m*(pr_m*pr_m*pr_m)/1e27)
!    dTmdP_m = dTmdP_m / 1e9

!DS 04/26/18 5 lines below are from Leah's implementation of zigsteg that reproduced Labrosse's model using params: 8676.00 -40.00 -0.d0 0.d0  0.d0 0.d0
      Tm_c =    243.d0 * (1/log(10.d0)) * (log(135.75 - pr_m/1e9)) + 4538.2
      dTmdP_m = 243.d0 * (1/log(10.d0)) * (-1/(135.75 - pr_m/1e9))
      dTmdP_m = dTmdP_m / 1e9
! these values were overriden in the original version V1 correpsonding to how rm was grided in rpts
! specifically the first two values in rm = rc 
!     Tm_c(0) = Tm_c(2)
!     Tm_c(1) = Tm_c(2)

! trying this to be consistent with how rm gets grided in rpts in this version
      Tm_c(0) = Tm_c(1)

  end subroutine melting_poly_bmo



!******************************************************************************

! Solid_conc: get c_elm2 from cL

!  cl read in in input file, so specified by user

!******************************************************************************

  subroutine solid_conc_bmo(cl, lambdal, lambdas, dmu, cs)



    double precision, intent(in)  :: cl, lambdal, lambdas, dmu

    double precision, intent(out) :: cs

    double precision              :: xmin, xmax, TmFe, dsri



    xmin = 0.1d-10

    xmax = 1.2d0

    TmFe = Tm_c(Nb+1)

    dsri = ds_m(Nb+1)



    cs = 0.d0

    if(cl==0.d0) then

      cs = 0.d0

      return

    endif



    cs = rtbis(xmin, xmax, 1d-10, lambdal, lambdas, dmu ,cl ,TmFe, dsri)



  end subroutine solid_conc_bmo



  double precision function f(x,lambdal,lambdas,dmu,cl,TmFe,dSFe)

    implicit none



    double precision, intent(in) :: x, TmFe, dSFe, cl, lambdas, lambdal, dmu



    f  = dmu + lambdal*cl - lambdas*x  - kb*TmFe*dlog(x/cl)*( dSFe + (x-cl) ) /dSFe



  end function f



  double precision function rtbis( x1, x2, xacc, lambdal,lambdas,dmu,cl,TmFe,dSFe)



    integer , parameter :: n=200

    integer             :: i

    double precision    :: x1, x2, xacc

    double precision    :: fst, fmid, dx, xmid

    double precision, intent(in) :: TmFe, dSFe, cl, lambdas, lambdal, dmu



    fmid = f(x2,lambdal,lambdas,dmu,cl,TmFe,dSFe)

    fst  = f(x1,lambdal,lambdas,dmu,cl,TmFe,dSFe)



    if(fst*fmid .ge. 0.d0) stop 'Root must be bracketed for bisection'



    if(fst .le. 0.d0) then

      rtbis = x1

      dx   = x2-x1

    else

      rtbis = x2

     dx   = x1-x2

    endif



    do i = 1, n

      dx = dx * 0.5d0

      xmid = rtbis + dx

      fmid = f(xmid,lambdal,lambdas,dmu,cl,TmFe,dSFe)



      if(fmid .le. 0.d0) rtbis=xmid

      if(dabs(dx) .lt. xacc .or. fmid .eq. 0.d0) return

    enddo



  end function rtbis



!******************************************************************************

! melting_pt_dep_bmo: The melting point depression for each light element

! LZ:  See J Chem Ph. Vol 116 Alfe et al, equation 12

! ds_c is in units here of K_b

!******************************************************************************

  subroutine melting_pt_dep_bmo(Nm, cl, cs, dTm)



    integer, intent(in) :: Nm

    double precision, intent(in)  :: cs, cl

    double precision, intent(out) :: dTm(0:Nm)

    double precision :: ds(0:Nm)

    integer :: i



    ds = ds_m / (kB * Ev * Na * 1000.d0  / AXx)



!   CHK - Assumed kb/ds = 1

    dTm = 0.d0

    dTm = (Tm_c) *  (cs - cl)



  end subroutine melting_pt_dep_bmo



!******************************************************************************

! conductivity_poly: Thermal conductivity using equation (17)

!       of Davies (2014, submitted to PEPI)

!******************************************************************************

  subroutine conductivity_poly_bmo()



    integer :: i



    do i = 0, Nm

      kr_m(i) =  k2_bm*r_m(i)*r_m(i)/1e6 + k1_bm*r_m(i)/1e3 + k0_bm

    enddo



    write(*,'(A16,E16.6)') 'k(ro)=',kr_m(Nm)



  end subroutine conductivity_poly_bmo



!******************************************************************************

! ricalc_bmo: find the layers boundary from the intersection of Ta and Tm

!

! Inputs: Nm - total number of gridpoints

!

! Outputs: rb - Interface radius (bottom point of top layer)

!          lflag - layer(s) scenario index, one or two layers only

!

!--------------

! LZ Update 4/2015: generalize the grid search to allow crystalization to

!                  progress from top down or bottom up

!******************************************************************************

  subroutine ricalc_bmo(Nm,rb,lflag)



    implicit none



    integer, intent(in)           :: Nm

    integer, intent(out)           :: lflag

    double precision, intent(out) :: rb

    double precision              :: Tdiff,Tdiffn,Tdd

    integer                       :: i, ript, test



    Tdiff  = 0.d0

    Tdiffn = 0.d0



!   The two T profiles are on the same grid at each time pt

!   so we can just compare then pointwise



    ript = 0

    test = 0

    lflag = 0



! lflag :: 1 = one liquid layer, 2 = liquid bottom layer + solid upper layer,

!          3 = one solid layer, 4 = solid bottom + liquid top

! ript will remain 0 if there is no interface (i e. one layer only)

! ok, june 22, need to update ript, since running now as one layer bmo

! changing so if only one layer, ript is Nm



!note, want rb to be top of the bottom layer, it was one index off in

!C. davies code



    a: do i = 0, Nm-1

         Tdiff = Tar_m(i) - Tm_c(i)

         Tdiffn = Tar_m(i+1) - Tm_c(i+1)

         if(Tdiff .ge. 0.d0 ) then         !Found liquid

            lflag = 1

            if(Tdiffn .le. 0.d0) then      !Found solid

               ript   = i+1

               test   = 1

               lflag  = 2

               exit a

            else                           !Found more liquid

               ript   = i+1

            endif

         endif

         if(Tdiff .le. 0.d0 ) then

            lflag = 3

            if(Tdiffn .ge. 0.d0) then

               ript   = i+1

               test   = 1

               lflag  = 4

               exit a

            else

               ript   = i+1

            endif

         endif

       enddo a



    rb = r_m(ript)



    if(rb < rc) stop "rb < rc!"



    write(*,'(A16,E16.6,I16,I16)') 'rb =',rb,ript,lflag



  end subroutine ricalc_bmo



!******************************************************************************

! crfac_bmo: Calculate Cr using equation (50)

!       of Francis Nimmo's revised treatise art`icle on Energetics of the core

!

! Inputs: Nb - Index for ICB radius

!         Tx - CMb temperature

!         rb - Current ICB radius

!

! Outputs: Cr

!******************************************************************************

  subroutine crfac_bmo(Nb,Tx,rb,Cr)



    integer, intent(in)           :: Nb

    double precision, intent(in)  :: Tx,rb

    double precision, intent(out) :: Cr

    double precision :: dTadP



!   If the temperature at the centre of the Earth is not below the

!   melting temperature at the same pressure then we have no IC

    if(rb==0.d0) then

      dTadP = 0.d0

      Cr = 0.d0

    else

      dTadP = -dTadr_m(Nb+1)/(rhor_m(Nb+1)*gr_m(Nb+1)) !Adiabat at ICB

      Cr    = -(Tar_m(Nb+1)/Tx)*(1.d0/(rhor_m(Nb+1)*gr_m(Nb+1)))*(1.d0/(dTmdP_m(Nb+1) - dTadP))

    endif



    write(*,'(A16,3E16.6)') 'Cr =',Cr, Tar_m(Nb+1), Tx



    write(14,'(10E16.6)') dTadP, dTmdP_m(Nb+1),rhor_m(Nb+1),Tar_m(Nb+1),gr_m(Nb+1),pr_m(Nb+1),Cr,dTadr_m(Nb+1), &

                         -dTmdP_m(Nb+1)*rhor_m(Nb+1)*gr_m(Nb+1), T_cen_b



  end subroutine crfac_bmo



!******************************************************************************

! ccfac_bmo: Calculate Cc using equation (51)

!        of Francis Nimmo's revised treatise article on Energetics of the core

!

! Inputs: Nb - Index for ICB

!         conc - O concentration

!

! Outputs: Cc

!******************************************************************************

  subroutine ccfac_bmo(Nb, cl, cs, Cc)



    integer, intent(in) :: Nb

    double precision, intent(in)  :: cl, cs

    double precision, intent(out) :: Cc



    Cc = 0.d0

    Cc = 4.d0*pi*r2_m(Nb+1)*rhor_m(Nb+1) * (cl  - cs) / M_om



    write(*,'(A16,3E16.6)') 'Cc =',Cc, cl, cs



  end subroutine ccfac_bmo



!******************************************************************************

! dcdt_bmo: Calculate the time change of concentration, dcdt_bmo, using equation (51)

!       of Francis Nimmo's revised treatise article on Energetics of the core

!

! Inputs: Cc - (see ccfac_bmo)

!         Cr - (see crfac_bmo)

!         dTxdt - CMB cooling rate

!         conc - O ocncentration

!         forback - Going forwards (=1) or backwards (=2) in time

!******************************************************************************

  subroutine dcdt_bmo(Cc, Cr, dTxdt, conc, forback)



    double precision, intent(in)  :: Cc, Cr, dTxdt

    double precision, intent(out) :: conc

    double precision              :: dconcdt

    integer                       :: forback



    dconcdt = Cc * Cr * dTxdt * secinyr



    if(forback==1) conc    = conc + dconcdt * dt_b

    if(forback==2) conc    = conc - dconcdt * dt_b



    write(*,'(A16, E16.6)') 'Dc/Dt = ', dconcdt



  end subroutine dcdt_bmo



  subroutine secular_num(Nt,Nb,Tx,Qs,Es)



    integer, intent(in)           :: Nt, Nb

    double precision, intent(in)  :: Tx

    double precision, intent(out) :: Qs, Es

    double precision  :: integrand1(0:Nt), integrand2(0:Nt)

    integer           :: i



    Qs=0.d0; Es=0.d0

    integrand1=0.d0; integrand2=0.d0

    do i = 0, Nt

      integrand1(i) = cp_m * r_m(i)*r_m(i) * rhor_m(i) * Tar_m(i)

      integrand2(i) = cp_m * r_m(i)*r_m(i) * rhor_m(i) * (Tar_m(i)/Tx - 1.d0)

    enddo

    call splat(integrand1,r_m,Nt+1)

    call splat(integrand2,r_m,Nt+1)



    Qs = -4.d0 * pi * integrand1(Nt) / Tx

    Es = -4.d0 * pi * integrand2(Nt) / Tx



  end subroutine secular_num



  subroutine secular_poly_bmo(Nm,Nb,Tx,Qs,Es)



    integer, intent(in)           :: Nm, Nb

    double precision, intent(in)  :: Tx

    double precision, intent(out) :: Qs, Es

    double precision              :: one_coeff, two_coeff, three_coeff, four_coeff, five_coeff, six_coeff, seven_coeff

    double precision              :: Soc_cmb, Soc_icb, Sic_icb, int_rhoT_dV



    one_coeff=0.d0; two_coeff=0.d0; three_coeff=0.d0; four_coeff=0.d0; five_coeff=0.d0; six_coeff=0.d0; seven_coeff=0.d0



    one_coeff    = rho0_om*T_cen_b

    two_coeff    = rho0_om*T_cen_b*t1_m   + rho1_om*T_cen_b

    three_coeff  = rho2_om*T_cen_b      + rho1_om*T_cen_b*t1_m   + rho0_om*T_cen_b*t2_m

    four_coeff   = rho2_om*T_cen_b*t1_m   + rho1_om*T_cen_b*t2_m   + rho3_om*T_cen_b    + rho0_om*T_cen_b*t3_m

    five_coeff   = rho2_om*T_cen_b*t2_m   + rho3_om*T_cen_b*t1_m   + rho1_om*T_cen_b*t3_m

    six_coeff    = rho3_om*T_cen_b*t2_m   + rho2_om*T_cen_b*t3_m

    seven_coeff  = rho3_om*T_cen_b*t3_m



    print *, Nb,Nm



    Soc_cmb = one_coeff*(r3_m(Nm)/1e9)/3.d0  +  two_coeff*(r4_m(Nm)/1e12)/4.d0 + three_coeff*(r5_m(Nm)/1e15)/5.d0 + &

              four_coeff*(r6_m(Nm)/1e18)/6.d0 + five_coeff*(r7_m(Nm)/1e21)/7.d0 + six_coeff*(r8_m(Nm)/1e24)/8.d0 + &

              seven_coeff*(r9_m(Nm)/1e27)/9.d0

    Soc_icb = one_coeff*(r3_m(Nb)/1e9)/3.d0  +  two_coeff*(r4_m(Nb)/1e12)/4.d0 + three_coeff*(r5_m(Nb)/1e15)/5.d0 + &

              four_coeff*(r6_m(Nb)/1e18)/6.d0 + five_coeff*(r7_m(Nb)/1e21)/7.d0 + six_coeff*(r8_m(Nb)/1e24)/8.d0 + &

              seven_coeff*(r9_m(Nb)/1e27)/9.d0



    one_coeff=0.d0; two_coeff=0.d0; three_coeff=0.d0; four_coeff=0.d0; five_coeff=0.d0; six_coeff=0.d0; seven_coeff=0.d0



    one_coeff    = rho0_bm*T_cen_b

    two_coeff    = rho0_bm*T_cen_b*t1_m

    three_coeff  = rho2_bm*T_cen_b      + rho0_bm*T_cen_b*t2_m

    four_coeff   = rho2_bm*T_cen_b*t1_m   + rho0_bm*T_cen_b*t3_m

    five_coeff   = rho2_bm*T_cen_b*t2_m

    six_coeff    = rho2_bm*T_cen_b*t3_m



    Sic_icb = one_coeff*(r3_m(Nb)/1e9)/3.d0  +  two_coeff*(r4_m(Nb)/1e12)/4.d0 + three_coeff*(r5_m(Nb)/1e15)/5.d0 + &

              four_coeff*(r6_m(Nb)/1e18)/6.d0 + five_coeff*(r7_m(Nb)/1e21)/7.d0 + six_coeff*(r8_m(Nb)/1e24)/8.d0



    int_rhoT_dV = 4.d0 * pi * (Soc_cmb - Soc_icb + Sic_icb) * 1e9



    Qs = -cp_m * int_rhoT_dV / Tx



    Es =  cp_m * (M_m - int_rhoT_dV/Tx) /Tx



  end subroutine secular_poly_bmo



  subroutine radiogenic_decay_bmo(time_b, h0_m, h, forback)


    integer               :: forback

    double precision, intent(in)    :: time_b, h0_m

    double precision, intent(inout) :: h

    double precision                :: t,h1,h2,h3,h4

    double precision                :: thuratio,kuratio,urat



     print*,forback,t,'forback,t before'
! If forback = 1, forward; 2=backward
      if(forback==1) then
         t = ttot_b-time_b*1.e6
      elseif(forback==2) then
         t = ttot_b - time_b*1.e6
     endif

     print*,forback,t,ttot_b,'forback,t,ttot_b  after'


!    h = h0_m * 2**(t/hl)

!   where h0_m is the ppb Uranium



     urat = 21.0 / 20.0    ! labrosse factor, not sure why

     thuratio = 1.0/20.0 * 75.0  ! 3.75 labrosse, we might prefer 4

     kuratio = 1.0/20.0 * 2.7e5  ! 13500 labrosse, we might prefer 10000



     h1 = urat * 0.0071 * h0_m * 1e-9 * HU235 * exp( t * log(2.d0) / TU235)

     h2 = urat * 0.9928 * h0_m * 1e-9 * HU238 * exp( t * log(2.d0) / TU238)

     h3 = thuratio * h0_m * 1e-9 * HTH * exp( t * log(2.d0) / TTH)

     h4 = kuratio * 1.19e-4 * h0_m * 1e-9 * HK40 * exp( t * log(2.d0) / TK40)

     h = h1 + h2+ h3 + h4



  end subroutine radiogenic_decay_bmo



  subroutine radiogenic_poly_bmo(Nt,Tx,h,Qr,Er)



    integer, intent(in)           :: Nt

    double precision, intent(in)  :: Tx, h

    double precision, intent(out) :: Qr, Er

    double precision              :: integrand(0:Nm),M_mtmp

    integer :: i



    M_mtmp = 4e24

    Qr = M_mtmp*h





    do i = 0, Nt

      integrand(i) = r_m(i)*r_m(i)*rhor_m(i)*(1.d0/Tar_m(Nt) - 1.d0/Tar_m(i))

    enddo

    call splat(integrand,r_m,Nt+1)

    Er=4.d0*pi*integrand(Nt)*h!1.d12



  end subroutine radiogenic_poly_bmo



!I see that you have put the latent heat in the table back to 0.75, are you sure about that ?

!with an entropy of melting of 1.05 and a temperature of 6350 K it seems to me that T\Delta S

!should be 0.99.



  subroutine lh_coeff_bmo(N,L)



    integer         , intent(in) :: N

    double precision, intent(out) :: L



    L = 0.d0

    L = Tm_c(N) * ds_m(N) !* 1.3806488e-23 * Na * 1000.d0  / AXx



  end subroutine lh_coeff_bmo



  subroutine latent_bmo(Nb,L,Tx,dridt,Ql,El)



    integer         , intent(in)  :: Nb

    double precision, intent(in)  :: Tx, dridt, L

    double precision, intent(out) :: Ql, El



! NOTE:  for BMO, adding a (-) sign.  dridt is oppos. sign than core case,

! since it varies the opposite with T (e.g. as T increases, r of bmo increases, but

! r of inner core would decrease).  I've not changed Cr calc, it's all fine,

! but we need to factor in that latent heat is produced as r shrinks

!    Ql = 4.d0*pi*r2_m(Nb+1)*L*rhor_m(Nb+1)*dridt

! CD - CHK

    Ql = -4.d0*pi*r2_m(Nb+1)*L*rhor_m(Nb+1)*dridt



 write(*,*) '--Ql terms (r2_m), L, rhor_m, dridt:', r2_m(Nb+1),L,rhor_m(Nb+1),dridt

    El = Ql*(Tar_m(Nb+1)-Tx)/(Tar_m(Nb+1)*Tx)



  end subroutine latent_bmo



  subroutine gravitational_num_bmo(N,Ni,Tc,Cr,Cc,alphac,Qg,Eg)



    integer, intent(in)          :: N, Ni

    double precision, intent(in) :: Tc, Cr, Cc, alphac

    double precision             :: integrand(0:N)

    double precision, intent(out):: Qg, Eg

    integer                      :: i



    do i = 0, N

      integrand(i)=psir_m(i)*rhor_m(i)*r_m(i)*r_m(i)

    enddo

    call splat(integrand,r_m,N+1)

    Qg = 4.d0 * pi * (integrand(N)-integrand(N-Ni+1))

    Qg = (Qg - M_m * psir_m(Nb)) * (Cr*Cc) * alphac

    Eg = Qg/Tc



    write(*,*) '**Grv Test**'

    write(*,*) 'Qg, Eg = ', Qg, Eg

    write(*,*) r_m(N), r_m(Ni), integrand(N), integrand(N-Ni+1)



  end subroutine gravitational_num_bmo



! ASSUMES - Only O contributes to the density jump.

  subroutine gravitational_poly_bmo(Nm, Nb, Tx, Cr, Cc, alphac, Qg, Eg)



    integer, intent(in)           :: Nm,Nb

    double precision, intent(in)  :: Tx, Cr, Cc, alphac

    double precision, intent(out) :: Qg, Eg



    double precision    :: gpot_cmb,gpot_icb,gpot_icbi,rho_psicmb_cmb,rho_psicmb_icb,rho_psi_cmb, rho_psi_icb

    double precision    :: five_coeff,six_coeff,seven_coeff,eight_coeff,nine_coeff,ten_coeff,eleven_coeff

    double precision    :: int_rhopsi_dV



    gpot_cmb=0.d0; rho_psicmb_cmb=0.d0; rho_psicmb_icb=0.d0; rho_psi_cmb=0.d0; rho_psi_icb=0.d0



    gpot_cmb = -(rho0_om*(r2_m(Nm)/1e6)/6.d0      + rho1_om*(r3_m(Nm)/1e9)/12.d0 + &

                 rho2_om*(r4_m(Nm)/1e12)/20.d0    + rho3_om*(r5_m(Nm)/1e15)/30.d0)! - &



    rho_psicmb_cmb = gpot_cmb*(rho0_om*(r3_m(Nm)/1e9)/3.d0   + rho1_om*(r4_m(Nm)/1e12)/4.d0  + &

                               rho2_om*(r5_m(Nm)/1e15)/5.d0  + rho3_om*(r6_m(Nm)/1e18)/6.d0)

    rho_psicmb_icb = gpot_cmb*(rho0_om*(r3_m(Nb)/1e9)/3.d0  + rho1_om*(r4_m(Nb)/1e12)/4.d0 + &

                               rho2_om*(r5_m(Nb)/1e15)/5.d0 + rho3_om*(r6_m(Nb)/1e18)/6.d0)



    five_coeff   = rho0_om*rho0_om/30.d0

    six_coeff    = rho0_om*rho1_om/24.d0

    seven_coeff  = rho1_om*rho1_om/84.d0  + rho0_om*rho2_om*13.d0/420.d0

    eight_coeff  = rho1_om*rho2_om/60.d0  + rho0_om*rho3_om/40.d0

    nine_coeff   = rho2_om*rho2_om/180.d0 + 7.d0*rho1_om*rho3_om/540.d0

    ten_coeff    = rho2_om*rho3_om/120.d0

    eleven_coeff = rho3_om*rho3_om/330.d0



    rho_psi_cmb = five_coeff*(r5_m(Nm)/1e15) + six_coeff*(r6_m(Nm)/1e18)  + seven_coeff*(r7_m(Nm)/1e21) + &

                  eight_coeff*(r8_m(Nm)/1e24)+ nine_coeff*(r9_m(Nm)/1e27) + ten_coeff*(r10_m(Nm)/1e30)  + &

                  eleven_coeff*(r11_m(Nm)/1e33)



    rho_psi_icb = five_coeff*(r5_m(Nb)/1e15) + six_coeff*(r6_m(Nb)/1e18)  + seven_coeff*(r7_m(Nb)/1e21) + &

                  eight_coeff*(r8_m(Nb)/1e24)+ nine_coeff*(r9_m(Nb)/1e27) + ten_coeff*(r10_m(Nb)/1e30)  + &

                  eleven_coeff*(r11_m(Nb)/1e33)



    int_rhopsi_dV = 16.d0 * pi * pi * G * (rho_psi_cmb - rho_psi_icb + rho_psicmb_cmb - rho_psicmb_icb) *1e15



    Qg = Cc * Cr  * alphac * (int_rhopsi_dV - M_om*psir_m(Nb))



    Eg = Qg/Tx



  end subroutine gravitational_poly_bmo



  subroutine heatofreaction_poly_bmo(Nm,Nb,dridt,Cc,Eh)



    integer         , intent(in)  :: Nm, Nb

    double precision, intent(in)  :: dridt,Cc

    double precision              :: integrand(Nm-Nb), roc(Nm-Nb)

    double precision, intent(out) :: Eh

    integer                       :: i



    integrand = 0.d0



    do i = Nb+1, Nm

      integrand(i-Nb) = r_m(i)*r_m(i)*rhor_m(i)*(1.d0/Tar_m(i))

      roc(i-Nb)       = r_m(i)

    enddo

    call splat(integrand,roc,(Nm-Nb))



    Eh = Rh * ( (4.d0 * pi * integrand(Nm-Nb)) - (M_om/Tar_m(Nb+1))) * Cc * dridt



  end subroutine heatofreaction_poly_bmo



  subroutine heatofreaction_poly_correct_bmo(Nm,Nb,dridt,Cc,Eh)



    integer         , intent(in)  :: Nm, Nb

    double precision, intent(in)  :: dridt,Cc

    double precision, intent(out) :: Eh

    double precision              :: integrand(Nm+1), roc(Nm-Nb), dmudT(0:Nm)

    integer                       :: i



    dmudT = dmudT0 + dmudT1*r_m/1e3

    dmudT = dmudT * Ev * Na * 1000.d0 / (Aelm1*1.d0)



!print *, dmudT(0), dmudT(Nm)



    integrand = 0.d0

    do i = 0, Nm

      integrand(i+1) = r_m(i)*r_m(i)*rhor_m(i)*dmudT(i)

!      roc(i-Nb)       = r_m(i)

    enddo

    call splat(integrand,r_m,Nm+1)



    Eh = -4.d0 * pi * Cc * dridt * ( integrand(Nm) - integrand(Nb) )



  end subroutine heatofreaction_poly_correct_bmo



! ASSUMES - alpha_c and alpha_D are constants.

!           core is well-mixed and thermodiffusion is negligible (grad(c) & K_T = 0 in eq 15 of Gub etal 14)

!           IS THERE AN ERROR IN EQN 32 OF GUBBINS ET AL 2004?

  subroutine baro_entropy_bmo(Nbc, cbar_elm1, cbar_elm2, cbar_elm3, Ealpha)



    integer         , intent(in):: Nbc

    double precision, intent(in):: cbar_elm1, cbar_elm2, cbar_elm3

    double precision            :: alphaDelm1, alphaDelm2, alphaDelm3

    double precision            :: dmubar_dcbar_elm1, dmubar_dcbar_elm2, dmubar_dcbar_elm3

    double precision            :: dmudc_elm1       , dmudc_elm2       , dmudc_elm3

    double precision            :: Tbar, rhobar, Abar

    double precision            :: integrand(Nm+1), roc(Nm)

    double precision            :: Ealpha_elm1, Ealpha_elm2, Ealpha_elm3, Ealpha

    integer                     :: i



    Abar = cbar_elm1*Aelm1 + cbar_elm2*Aelm2 + cbar_elm3*Aelm3 + (1.d0-cbar_elm1-cbar_elm2-cbar_elm3)*AXx



    Tbar   = (Tar_m(Nm) + Tar_m(0)) / 2

    rhobar = (rhor_m(Nm) + rhor_m(0)) / 2



    dmubar_dcbar_elm1  = kb*Tbar/cbar_elm1  + lambda_elm1_om

    dmubar_dcbar_elm2  = kb*Tbar/cbar_elm2  + lambda_elm2_om

    dmubar_dcbar_elm3 = kb*Tbar/cbar_elm3 + lambda_elm3_om



    dmudc_elm1  = dmubar_dcbar_elm1  * (Abar/Aelm1)  * Ev * Na * (1000.d0/Aelm1)

    dmudc_elm2  = dmubar_dcbar_elm2  * (Abar/Aelm2)  * Ev * Na * (1000.d0/Aelm2)

    dmudc_elm3 = dmubar_dcbar_elm3 * (Abar/Aelm3) * Ev * Na * (1000.d0/Aelm3)



    alphaDelm1  = rhobar * Delm1  / dmudc_elm1

    alphaDelm2  = rhobar * Delm2  / dmudc_elm2

    alphaDelm3 = rhobar * Delm3 / dmudc_elm3



    integrand = 0.d0

    roc = 0.d0

    do i = 0, Nm

      integrand(i+1) = r_m(i)*r_m(i)*gr_m(i)*gr_m(i)*(1.d0/Tar_m(i))

      !roc(i-Nbc)       = r_m(i)

    enddo

    call splat(integrand,r_m,Nm+1)



    Ealpha_elm1  = 4.d0* pi * alphac_c_elm1 *alphac_c_elm1 *alphaDelm1  * (integrand(Nm)-integrand(Nbc))

    Ealpha_elm2  = 4.d0* pi * alphac_c_elm2 *alphac_c_elm2 *alphaDelm2  * (integrand(Nm)-integrand(Nbc))

    Ealpha_elm3 = 4.d0* pi * alphac_c_elm3*alphac_c_elm3*alphaDelm3 * (integrand(Nm)-integrand(Nbc))



    Ealpha   = Ealpha_elm1 + Ealpha_elm2 + Ealpha_elm3



    write(*,'(A16,2E16.6)') 'Ealpha_elm1=',Ealpha_elm1, alphaDelm1

    write(*,'(A16,2E16.6)') 'Ealpha_elm2=',Ealpha_elm2, alphaDelm2

    write(*,'(A16,2E16.6)') 'Ealpha_elm3=',Ealpha_elm3, alphaDelm3



  end subroutine baro_entropy_bmo



  subroutine cond_entropy_poly_bmo(Nm,Ea)



    integer, intent(in)           :: Nm

    double precision, intent(out) :: Ea

    double precision              :: integrand(Nm)

    integer                       :: i



    integrand = 0.d0



    do i = 1,Nm

        integrand(i)=kr_m(i)*(r_m(i)*dTadr_m(i)/Tar_m(i))**2

    enddo

    call splat(integrand,r_m,Nm)



    Ea=4.d0*pi*integrand(Nm)



  end subroutine cond_entropy_poly_bmo



  subroutine open_files_bmo(char)



    character (len=4)   :: char

    character (len=100) :: energy_file, entropy_file, diag_file, pro_file, icb_file, conc_file

    character (len=100) :: pro1k_file



    energy_file  = out_file_bmo(1:out_file_bmo_len)//"_"//char//'_energy'

    entropy_file = out_file_bmo(1:out_file_bmo_len)//"_"//char//'_entropy'

    diag_file    = out_file_bmo(1:out_file_bmo_len)//"_"//char//'_diagnostics'

    pro_file     = out_file_bmo(1:out_file_bmo_len)//"_"//char//'_profiles'

    icb_file     = out_file_bmo(1:out_file_bmo_len)//"_"//char//'_icb'

    conc_file    = out_file_bmo(1:out_file_bmo_len)//"_"//char//'_conc'

    pro1k_file    = out_file_bmo(1:out_file_bmo_len)//"_"//char//'_1kprofile'



    open(10,file=energy_file)

    open(11,file=entropy_file)

    open(12,file=diag_file)

    open(13,file=pro_file)

    open(14,file=icb_file)

    open(15,file=conc_file)

    open(16,file=pro1k_file)

    write(10,'(10A16)') 'Time (Myr)', 'Qs', 'Qg', 'Ql', 'Qr', 'Qx', 'Qa', 'dTdr_c(CMB)', 'dTrsdr', 'dTadr_m(rr)'

    write(11,'(9A16)')  'Time (Myr)', 'Es', 'Eg', 'El', 'Er', 'EJ_b', 'Ea', 'Eh', 'Ealpha'

    write(12,'(8A16)')  'Time (Myr)', 'Tx', 'dTxdt', 'rb', 'rr', 'dridt', 'M_m', 'M_om'

    write(13,'(9A16)')  'r', 'T', 'dTadr_m', 'Tm', 'p', 'rho', 'g', 'psi', 'k'

    write(14,'(10A16)') 'dTadP(rb)','dTmdP(rb)','rho(rb)','T(rb)','g(rb)','p(rb)','Cr','dTadr(rb)','dTmdr(rb)','T_cen_b'

    write(15,'(14A16)') 'Time (Myr)','rb','clO','clS','clSi','csO','csS','csSi','QgO','QgS','QgSi','dTmO', 'dTmS', 'dTmSi'

    write(16,'(7A16)')  'r', 'T', 'dTadr_m', 'Tm', 'p','g','rho'



  end subroutine open_files_bmo



  subroutine close_files_bmo()



    close(10)

    close(11)

    close(12)

    close(13)

    close(14)

    close(15)

    close(16)



  end subroutine close_files_bmo



  subroutine write_profiles_bmo()



    implicit none



    integer :: j



    do j = 0, Nm

      write(13,'(9E16.8)') r_m(j), Tar_m(j), dTadr_m(j), Tm_c(j), pr_m(j), rhor_m(j), gr_m(j), psir_m(j), kr_m(j)

    enddo

    write(13,*)

    write(13,*)



  end subroutine write_profiles_bmo



end module bmo_model

