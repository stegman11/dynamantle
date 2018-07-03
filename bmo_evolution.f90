  module bmo_evolution



    use bmo_model

    use parameters



    implicit none



    double precision, parameter :: Mc_labrosse=2e24, Cpc_labrosse=860e0

    integer          :: Nbc, i, j, z

    double precision :: Qa_b, Qst_b, Qr_b, Qlt_b, Qgt_b, Qtt_b

    double precision :: Ealpha_b, Ea_b, Est_b, Er_b, Elt_b, Egt_b, Ett_b, Eht_b

    double precision :: dTplusdr_b, dcplusdr_b, dTcmbdr_b, dTrsdr_b, Cr_b, dridt_b



  contains



  subroutine bmo_evol_calc(Qx,Tx,time_b,dTxdt,rb,EJ_b,cl_elm1,cl_elm2,cl_elm3,forback)



    double precision :: Qx, Tx, time_b, dTxdt, rb, EJ_b

    double precision :: Qg_elm1, Eg_elm1, Qg_elm2, Eg_elm2, Qg_elm3, Eg_elm3

    double precision :: LH_m

    double precision :: Cc_elm1, Cc_elm2, Cc_elm3, dTm_elm1(0:Nm), dTm_elm2(0:Nm), dTm_elm3(0:Nm)

    double precision :: cl_elm1, cl_elm2, cl_elm3, cs_elm1, cs_elm2, cs_elm3

    double precision :: clbar_elm1, clbar_elm2, clbar_elm3, csbar_elm1, csbar_elm2, csbar_elm3

    integer          :: forback, i, lflag



      call mass2moleconc_bmo(cl_elm1, cl_elm2, cl_elm3, clbar_elm1, clbar_elm2, clbar_elm3)

      call conductivity_poly_bmo()                        ! Get conductivity profile

      call density_poly_bmo(Nm, cl_elm1, alphac_c_elm1)

      call gravity_num_bmo(Nm)

      call grav_pot_num_bmo(Nm)

      call pressure_num_bmo(Nm,Pc_b)

      call adiabat_poly_bmo(Nm,Tx,Qa_b)                   !Find Ta on this grid

      call melting_poly_bmo()                             !Find Tm(P)

      call entropy_melting_bmo()

      csbar_elm1=0.d0; csbar_elm2=clbar_elm2; csbar_elm3=clbar_elm3

      call melting_pt_dep_bmo(Nm, clbar_elm1, csbar_elm1, dTm_elm1)    !Melting pt depression due to _elm1

      call melting_pt_dep_bmo(Nm, clbar_elm2, csbar_elm2, dTm_elm2)    !Melting pt depression due to _elm2

      call melting_pt_dep_bmo(Nm, clbar_elm3, csbar_elm3, dTm_elm3) !Melting pt depression due to _elm3

      call mole2massconc_bmo(csbar_elm1, csbar_elm2, csbar_elm3, cs_elm1, cs_elm2, cs_elm3)

      Tm_c = Tm_c + dTm_elm1 + dTm_elm2 + dTm_elm3

      do i = 0, Nm

        write(555,*) r_m(i), rhor_m(i), gr_m(i), psir_m(i), pr_m(i), Tm_c(i), Tar_m(i)

      enddo



      call ricalc_bmo(Nm, rb,lflag)            !Find inner core radius using Tm & Ta on old grid



      Nbc=Nb                                     !Nb is total number of pts used for IC

      if(rb == rc) Nbc=0                        !If no IC, put no pts in IC



      call rpts_bmo(rb,rt)                       !Set new grid based on rb



      call lh_coeff_bmo(Nbc, LH_m)

      write(*,'(A16,E16.6)') 'L = ', LH_m

      call conductivity_poly_bmo()                             ! Get conductivity profile

      call mass_poly_bmo(Nbc,Nm)                !Now compute everything on new grid

      call density_poly_bmo(Nm, cl_elm1, alphac_c_elm1)

      call gravity_num_bmo(Nm)

      call grav_pot_num_bmo(Nm)

      call pressure_num_bmo(Nm,Pc_b)

      call adiabat_poly_bmo(Nm,Tx,Qa_b)

      call melting_poly_bmo()

      call entropy_melting_bmo()

      call melting_pt_dep_bmo(Nm, clbar_elm1, csbar_elm1, dTm_elm1)    !Melting pt depression due to _elm1

      call melting_pt_dep_bmo(Nm, clbar_elm2, csbar_elm2, dTm_elm2)    !Melting pt depression due to _elm2

      call melting_pt_dep_bmo(Nm, clbar_elm3, csbar_elm3, dTm_elm3)    !Melting pt depression due to _elm3

      Tm_c = Tm_c + dTm_elm1 + dTm_elm2 + dTm_elm3



      call write_profiles_bmo()



      call crfac_bmo(Nbc,Tx,rb,Cr_b)                                                    !dridt_b

      call ccfac_bmo(Nbc, cl_elm1, cs_elm1,  Cc_elm1)

      call ccfac_bmo(Nbc, cl_elm2, cs_elm2,  Cc_elm2)

      call ccfac_bmo(Nbc, cl_elm3, cs_elm3, Cc_elm3)                                    !Dc/dt -> dTc/dt



      call latent_bmo(Nbc,LH_m,Tx,Cr_b,Qlt_b,Elt_b)                                     !Latent heat



      call radiogenic_decay_bmo(time_b, h0_m, h, forback)

      call radiogenic_poly_bmo(Nbc,Tx,h,Qr_b,Er_b)                                        !Radiogenic



      Qgt_b = 0.d0; Egt_b = 0.d0

      call secular_num(Nbc, 0, Tx, Qst_b, Est_b)

      call gravitational_num_bmo(Nbc,0,Tx,Cr_b,Cc_elm1 ,alphac_c_elm1 ,Qg_elm1,Eg_elm1) !Gravitational energy

      Qgt_b=Qg_elm1; Egt_b=Eg_elm1

      call gravitational_num_bmo(Nbc,0,Tx,Cr_b,Cc_elm2 ,alphac_c_elm2 ,Qg_elm2,Eg_elm2) !Gravitational energy

      Qgt_b=Qgt_b+Qg_elm2; Egt_b=Egt_b+Eg_elm2

      call gravitational_num_bmo(Nbc,0,Tx,Cr_b,Cc_elm3,alphac_c_elm3,Qg_elm3,Eg_elm3)   !Gravitational energy

      Qgt_b=Qgt_b+Qg_elm3; Egt_b=Egt_b+Eg_elm3

      call heatofreaction_poly_correct_bmo(Nm,Nbc,Cr_b,Cc_elm1,Eht_b)                    !Heat of reaction

      call cond_entropy_poly_bmo(Nm,Ea_b)                                                !Thermal conduction

      call baro_entropy_bmo(Nbc, clbar_elm1, clbar_elm2, clbar_elm3, Ealpha_b)           !Barodiffusion



!     Need to chk this

      Qgt_b = -1.d0*Qgt_b; Egt_b = -1.d0*Egt_b



      if(hor_b .ne. 1) Eht_b = 0.d0



      call bmo_en_ent_calc(Qx,time_b,dTxdt,rb,EJ_b,forback) !Relate EJ_b and Qx



      dridt_b = Cr_b*dTxdt                                  !IC growth rate

      dTcmbdr_b = -Qx/(4.d0*pi*rb*rb*kr_m(Nm))              !dTdr @ CMB



      write(10,'(11E16.6)') time_b, Qst_b*dTxdt, Qgt_b*dTxdt, Qlt_b*dTxdt, Qr_b, &

                            Qx, Qa_b, dTcmbdr_b, dTplusdr_b, dTadr_m(Nm), &

                           -Mc_labrosse*Cpc_labrosse*dTxdt

      write(11,'(9E16.6)')  time_b, Est_b*dTxdt, Egt_b*dTxdt, Elt_b*dTxdt, Er_b, EJ_b, Ea_b, Eht_b*dTxdt, &

                            Ealpha_b

      write(12,'(8E16.8)')  time_b, Tx, dTxdt*secingyr,rb, rr, dridt_b*secingyr, M_m ,M_om

      write(15,'(14E16.6)') time_b, rb, cl_elm1, cl_elm2, cl_elm3, cs_elm1, cs_elm2, cs_elm3, Qg_elm1*dTxdt, &

                            Qg_elm2*dTxdt, Qg_elm3*dTxdt, dTm_elm1(Nb+1), dTm_elm2(Nb+1), dTm_elm3(Nb+1)



      call dcdt_bmo(Cc_elm1 , Cr_b, dTxdt, cl_elm1 , forback)

      call dcdt_bmo(Cc_elm2 , Cr_b, dTxdt, cl_elm2 , forback)

      call dcdt_bmo(Cc_elm3, Cr_b, dTxdt, cl_elm3, forback)

      call radiogenic_decay_bmo(time_b, h0_m, h, forback)



      print *,'Tx, dTxdt, and dt_b are: ',Tx,dTxdt,dt_b

      if(forback==1) Tx  = Tx + dTxdt*secinyr*dt_b        !New CMB temperature

      if(forback==2) Tx  = Tx - dTxdt*secinyr*dt_b        !New CMB temperature



      if(Tx < 0.d0) then

        write(*,*) '********************'

        write(*,*) 'CMB temperature < 0!'

        write(*,*) 'time_b, Tx = ', time_b, Tx

        write(*,*) '********************'

        stop

      endif



  end subroutine bmo_evol_calc



  subroutine bmo_en_ent_calc(Qx,time_b,dTxdt,rb,EJ_b,forback)



    double precision, parameter :: Mc_labrosse=2e24, Cpc_labrosse=860e0

    double precision,save :: Qlast, Efirst

    double precision      :: time_b, rb, Qx, dTxdt, EJ_b

    integer               :: forback

    integer, save         :: qq



      Qtt_b = Qst_b + Qgt_b + Qlt_b - Mc_labrosse*Cpc_labrosse

      Ett_b = Est_b + Egt_b + Elt_b + Eht_b



      write(*,'(A16,E16.6)') 'Qst_l=',-Mc_labrosse*Cpc_labrosse

      write(*,'(A16,E16.6)') 'Qst_b=',Qst_b

      write(*,'(A16,E16.6)') 'Qgt_b=',Qgt_b

      write(*,'(A16,E16.6)') 'Qlt_b=',Qlt_b

      write(*,'(A16,E16.6)') 'Qtt_b=',Qtt_b

      write(*,'(A16,E16.6)') 'Est_b=',Est_b

      write(*,'(A16,E16.6)') 'Egt_b=',Egt_b

      write(*,'(A16,E16.6)') 'Elt_b=',Elt_b

      write(*,'(A16,E16.6)') 'Eht_b=',Eht_b

      write(*,'(A16,E16.6)') 'Ett_b=',Ett_b



      if(forback==1) then

        dTxdt = (Qx - Qr_b)/Qtt_b

!        dTxdt = (Qx - Qr_b - Qlt_b)/Qtt_b

        EJ_b    = Ett_b*dTxdt + Er_b - Ea_b - Ealpha_b

        write(*,'(A16, E16.6, I16)') 'time_b (Qx fixed) = ', time_b, i

      elseif(forback==2) then



        if((sol==0) .or. (sol==1)) then              !Qx specified

          dTxdt = (Qx - Qr_b)/Qtt_b

!           dTxdt = (Qx - Qr_b - Qlt_b)/Qtt_b

          EJ_b    = Ett_b*dTxdt + Er_b - Ea_b - Ealpha_b

          write(*,'(A16, E16.6, I16)') 'time_b (Qx fixed) = ', time_b, i



        elseif(sol==2) then                          !OH specified

          dTxdt = (EJ_b + Ea_b + Ealpha_b - Er_b)/Ett_b

          Qx     = Qtt_b*dTxdt + Qr_b

          write(*,'(A16, E16.6, I16)') 'time_b (E fixed) = ', time_b, i



        endif



      endif



  end subroutine bmo_en_ent_calc



! If forback = 1, forward; 2=backward

  subroutine get_time_bmo(i,time_b,dt_b,forback)



    integer          :: i, forback

    double precision :: dt_b, time_b



    if(forback==1)  time_b  = real(i)*dt_b/1e6

    if(forback==2)  time_b  = (ttot_b/1e6) - real(i)*dt_b/1e6



  end subroutine get_time_bmo



  end module bmo_evolution

