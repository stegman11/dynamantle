  program thermal_history

    use bmo_evolution

    implicit none

    double precision :: Qx, Tx, dTxdt, rb, time_b, EJ_b
    double precision :: cbarO_bm, cbarS_bm, cbarSi_bm, c_elm1, c_elm2, c_elm3
    integer          :: forback
    character(len=4) :: forbackc

    call read_input_bmo(Qx,Tx,rb,EJ_b,cbarO_bm, cbarS_bm, cbarSi_bm) ! Read the input file

!    forback=1                                              ! Tells code to go backwards in time (forwards =1)
!    forbackc="forw"
     forback=2                                              ! Tells code to go backwards in time (forwards =1)
     forbackc="back"
    call open_files_bmo(forbackc)                         ! Open output files

    rb0 = rb
    Nbc = Nb
    iteration_b=0

    if(forback==2) then

    do i = 0, ntb

      call get_time_bmo(i,time_b,dt_b,forback)
      write(*,*) '----------------------', i

      if(sol==0) then
        time_b = time_b_file(i)
        dt_b   = dt_b_file(i)
        Qx    = xflux(i)
 !       Qr_b = Qr_b_file(i)
 !       rb = rb_b_file(i)
 !       dridt_b = dridt_b_file(i)
      endif


      call rpts_bmo(rt,rb)                                       ! Create grid from cmb to cmb+1000k so Tm and Ta can be obtained
!DS commenting out line below which is a variation from what the previous version used and the reason for the undocumented change is unknown
!      call rpts_bmo(rb,rt)                                       ! Create grid from cmb to cmb+1000k so Tm and Ta can be obtained
      call mole2massconc_bmo(cbarO_bm, cbarS_bm, cbarSi_bm, c_elm1, c_elm2, c_elm3)
      call bmo_evol_calc(Qx,Tx,time_b,dTxdt,rb,EJ_b,c_elm1,c_elm2,c_elm3,forback)
      call mass2moleconc_bmo(c_elm1, c_elm2, c_elm3, cbarO_bm, cbarS_bm, cbarSi_bm)

      if(i==0) iteration_b = 1
      if(sol .ne. 0) xflux(ntb-i) = Qx

      rb0 = r_m(Nbc)

    enddo


    Tx = Tx + dTxdt*secinyr*dt_b                             ! Undo last increment of Tx
!    h = h0_m * 2**(((ttot_b-dt_b)/1e6)/hl)
    write(*,*) '---END OF BACKWARD CALC---'

    elseif(forback==1) then

    do i = 0, ntb

      call get_time_bmo(i,time_b,dt_b,forback)
      write(*,*) '----------------------', i

       if(sol==0) then
         time_b = time_b_file(ntb-i)
         dt_b   = dt_b_file(ntb-i)
         Qx    = xflux(ntb-i)
       endif

      call rpts_bmo(rt,rb)                                       ! Create grid from cmb to cmb+1000k so Tm and Ta can be obtained
      call mole2massconc_bmo(cbarO_bm, cbarS_bm, cbarSi_bm, c_elm1, c_elm2, c_elm3)
      call bmo_evol_calc(Qx,Tx,time_b,dTxdt,rb,EJ_b,c_elm1,c_elm2,c_elm3,forback)   !   arg1=forward(1) or backward(2)
      call mass2moleconc_bmo(c_elm1, c_elm2, c_elm3, cbarO_bm, cbarS_bm, cbarSi_bm)

      if(i==0) iteration_b = 1
      if(sol .ne. 0) xflux(ntb-i) = Qx

      rb0 = r_m(Nbc)

     enddo

    Tx = Tx + dTxdt*secinyr*dt_b                             ! Undo last increment of Tx
    write(*,*) '---END OF FORWARD  CALC---'

    endif

    call close_files_bmo()

    write(*,'(A16,F16.6)') 'Tcen = ', T_cen_b
    write(*,'(A16,F16.6)') 'Tx   = ', Tx
    write(*,'(A16,F16.6)') 'rb   = ', rb
    write(*,'(A16,E16.6)') 'Qcmb = ', Qx
    write(*,'(A16,E16.6)') 'h    = ', h



    deallocate(xflux)

  end program thermal_history
