  module parameters_bmo

    implicit none

    integer         , parameter :: Nr=1000, Nb=10001, Nm=Nr+Nb  !Number of pts in the outer, basal mantle & total
    double precision, parameter :: rc=3480e3     !CMB radius
    double precision, parameter :: rt=4480e3     !1000km above CMB radius, top of height to compare Ta and Tm
    double precision, parameter :: hl=1.248e3    !Half life of K40
    double precision, parameter :: dmudT0=-0.00023, dmudT1=-6.32184e-08
    double precision, parameter :: g_cmb = 10.935211   ! gravity at top of the core from c. davies outputs when whole core was molten

    double precision, parameter :: rho_zeroP = 7.9e3 !Density, zero pressure
    double precision, parameter :: K_zeroP   = 500e9   !Bulk modulus, zero pressure

    integer, parameter          :: Aelm1         = 16    , Aelm2         = 32    , Aelm3         = 28,   AXx=56
    double precision, parameter :: lambda_elm1_om = 3.25d0, lambda_elm2_om = 6.15d0, lambda_elm3_om = 3.6d0
    double precision, parameter :: lambda_elm1_bm = 0.00d0, lambda_elm2_bm = 5.9d0 , lambda_elm3_bm = 2.7d0
    double precision, parameter :: dmu_elm1       =-2.6d0 , dmu_elm2       =-0.25d0, dmu_elm3       =-0.05d0
    double precision, parameter :: Delm1         = 1e-8  , Delm2         = 5e-9  , Delm3         = 5e-9
    double precision, parameter :: alphac_c_elm1  = 1.1d0 , alphac_c_elm2  = 0.64  , alphac_c_elm3  = 0.87, alphaT_c=1.25d-5
    double precision, parameter :: Rh=-27.7d6

  contains

  end module parameters_bmo
