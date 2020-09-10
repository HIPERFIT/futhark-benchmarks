-- Rodinia's Myocyte benchmark translated to Futhark.
-- Code and comments based on
-- https://github.com/kkushagra/rodinia/blob/master/openmp/myocyte/
--
-- ==
--
-- input @ data/small.in
-- output @ data/small.out
-- input @ data/medium.in

----------------------
----- Ecc MODULE -----
----------------------

let pow(x: f32, y: f32): f32 = x ** y
let log10(x: f32): f32 = f32.log(x) / f32.log(10.0f32)

let fmod(a: f32, b: f32): f32 = ( a - b * r32(t32(a / b)) )

let fabs(a: f32): f32 = if (a < 0.0f32) then -a else a

let ecc [equs][pars] (timeinst:  f32, initvalu: [equs]f32, initvalu_offset: i32,
                      parameter: [pars]f32, parameter_offset: i32, finavalu: *[equs]f32 ): *[equs]f32 =
  #[unsafe]
  -- variable references
  let offset_1  = initvalu_offset
  let offset_2  = initvalu_offset+1
  let offset_3  = initvalu_offset+2
  let offset_4  = initvalu_offset+3
  let offset_5  = initvalu_offset+4
  let offset_6  = initvalu_offset+5
  let offset_7  = initvalu_offset+6
  let offset_8  = initvalu_offset+7
  let offset_9  = initvalu_offset+8
  let offset_10 = initvalu_offset+9
  let offset_11 = initvalu_offset+10
  let offset_12 = initvalu_offset+11
  let offset_13 = initvalu_offset+12
  let offset_14 = initvalu_offset+13
  let offset_15 = initvalu_offset+14
  let offset_16 = initvalu_offset+15
  let offset_17 = initvalu_offset+16
  let offset_18 = initvalu_offset+17
  let offset_19 = initvalu_offset+18
  let offset_20 = initvalu_offset+19
  let offset_21 = initvalu_offset+20
  let offset_22 = initvalu_offset+21
  let offset_23 = initvalu_offset+22
  let offset_24 = initvalu_offset+23
  let offset_25 = initvalu_offset+24
  let offset_26 = initvalu_offset+25
  let offset_27 = initvalu_offset+26
  let offset_28 = initvalu_offset+27
  let offset_29 = initvalu_offset+28
  let offset_30 = initvalu_offset+29
  let offset_31 = initvalu_offset+30
  let offset_32 = initvalu_offset+31
  let offset_33 = initvalu_offset+32
  let offset_34 = initvalu_offset+33
  let offset_35 = initvalu_offset+34
  let offset_36 = initvalu_offset+35
  let offset_37 = initvalu_offset+36
  let offset_38 = initvalu_offset+37
  let offset_39 = initvalu_offset+38
  let offset_40 = initvalu_offset+39
  let offset_41 = initvalu_offset+40
  let offset_42 = initvalu_offset+41
  let offset_43 = initvalu_offset+42
  let offset_44 = initvalu_offset+43
  let offset_45 = initvalu_offset+44
  let offset_46 = initvalu_offset+45

  -- variable references
  let parameter_offset_1  = parameter_offset

  -- decoded input initial data
  let initvalu_1  = initvalu[offset_1 ]
  let initvalu_2  = initvalu[offset_2 ]
  let initvalu_3  = initvalu[offset_3 ]
  let initvalu_4  = initvalu[offset_4 ]
  let initvalu_5  = initvalu[offset_5 ]
  let initvalu_6  = initvalu[offset_6 ]
  let initvalu_7  = initvalu[offset_7 ]
  let initvalu_8  = initvalu[offset_8 ]
  let initvalu_9  = initvalu[offset_9 ]
  let initvalu_10 = initvalu[offset_10]
  let initvalu_11 = initvalu[offset_11]
  let initvalu_12 = initvalu[offset_12]
  let initvalu_13 = initvalu[offset_13]
  let initvalu_14 = initvalu[offset_14]
  let initvalu_15 = initvalu[offset_15]
  let initvalu_16 = initvalu[offset_16]
  let initvalu_17 = initvalu[offset_17]
  let initvalu_18 = initvalu[offset_18]
  let initvalu_19 = initvalu[offset_19]
  let initvalu_20 = initvalu[offset_20]
  let initvalu_21 = initvalu[offset_21]
  let initvalu_23 = initvalu[offset_23]
  let initvalu_24 = initvalu[offset_24]
  let initvalu_25 = initvalu[offset_25]
  let initvalu_26 = initvalu[offset_26]
  let initvalu_27 = initvalu[offset_27]
  let initvalu_28 = initvalu[offset_28]
  let initvalu_29 = initvalu[offset_29]
  let initvalu_30 = initvalu[offset_30]
  let initvalu_31 = initvalu[offset_31]
  let initvalu_32 = initvalu[offset_32]
  let initvalu_33 = initvalu[offset_33]
  let initvalu_34 = initvalu[offset_34]
  let initvalu_35 = initvalu[offset_35]
  let initvalu_36 = initvalu[offset_36]
  let initvalu_37 = initvalu[offset_37]
  let initvalu_38 = initvalu[offset_38]
  let initvalu_39 = initvalu[offset_39]
  let initvalu_40 = initvalu[offset_40]

  -- decoded input parameters
  let parameter_1 = parameter[parameter_offset_1]

  -- constants
  let pi   = 3.1416f32
  let r    = 8314.0f32     in -- [J/kmol*K]
  let frdy = 96485.0f32    in -- [C/mol]
  let temp = 310.0f32      in -- [K] 310
  let foRT = frdy/r/temp
  let cmem = 1.3810e-10f32 in -- [F] membrane capacitance
  let qpow = (temp-310.0f32)/10.0f32

  -- Cell geometry
  let cellLength = 100.0f32   in  -- cell length []um
  let cellRadius = 10.25f32   in  -- cell radius []um
  let vcell = pi*pow(cellRadius,2.0f32)*cellLength*1.0e-15f32 in -- [L]
  let vmyo = 0.65f32*vcell
  let vsr = 0.035f32*vcell
  let vsl = 0.02f32 *vcell
  let vjunc = 0.0539f32 * 0.01f32 * vcell
  let j_ca_juncsl = 1.0f32/1.2134e12f32  in    -- [L/msec]
  let j_ca_slmyo  = 1.0f32/2.68510e11f32 in    -- [L/msec]
  let j_na_juncsl = 1.0f32/(1.6382e12f32 / 3.0f32 * 100.0f32) in -- [L/msec]
  let j_na_slmyo  = 1.0f32/(1.8308e10f32 / 3.0f32 * 100.0f32) in -- [L/msec]

  -- Fractional currents in compartments
  let fjunc = 0.11f32
  let fsl = 1.0f32 - fjunc
  let fjunc_CaL = 0.9f32
  let fsl_CaL = 1.0f32 - fjunc_CaL

  -- Fixed ion concentrations
  let cli = 15.0f32  in -- Intracellular Cl  []mM
  let clo = 150.0f32 in -- Extracellular Cl  []mM
  let ko = 5.4f32    in -- Extracellular K   []mM
  let nao = 140.0f32 in -- Extracellular Na  []mM
  let cao = 1.8f32   in -- Extracellular Ca  []mM
  let mgi = 1.0f32

  -- Nernst Potentials
  let ena_junc = (1.0f32/foRT)*f32.log(nao/initvalu_32)        in -- []mV
  let ena_sl = (1.0f32/foRT)*f32.log(nao/initvalu_33)          in -- []mV
  let ek = (1.0f32/foRT)*f32.log(ko/initvalu_35)               in -- []mV
  let eca_junc = (1.0f32/foRT/2.0f32)*f32.log(cao/initvalu_36) in -- []mV
  let eca_sl = (1.0f32/foRT/2.0f32)*f32.log(cao/initvalu_37)   in -- []mV
  let ecl = (1.0f32/foRT)*f32.log(cli/clo)                     in -- []mV

  -- Na transport parameters
  let gna =  16.0f32       in -- [#mS/uF]
  let gnaB = 0.297e-3f32   in -- [#mS/uF]
  let iBarNaK = 1.90719f32 in -- [#uA/uF]
  let kmNaip = 11.0f32     in -- []mM
  let kmko = 1.5f32        in -- []mM

  -- K current parameters
  let pNaK = 0.01833f32
  let dToSlow = 0.06f32 in -- [#mS/uF]
  let gToFast = 0.02f32 in -- [#mS/uF]
  let gkp = 0.001f32

  -- Cl current parameters
  let gClCa = 0.109625f32 in -- [#mS/uF]
  let gClB = 9.0e-3f32      in -- [#mS/uF]
  let kdClCa = 100.0e-3f32     -- []mM

  -- i_Ca parameters
  let pNa = 1.5e-8f32  in -- [#cm/sec]
  let pCa = 5.4e-4f32  in -- [#cm/sec]
  let pK = 2.7e-7f32   in -- [#cm/sec]
  let q10CaL = 1.8f32

  -- Ca transport parameters
  let ibarNCX = 9.0f32   in -- [#uA/uF]
  let kmCai = 3.59e-3f32 in -- []mM
  let kmcao = 1.3f32     in -- []mM
  let kmNai = 12.29f32   in -- []mM
  let kmnao = 87.5f32    in -- []mM
  let ksat = 0.27f32     in -- []none
  let nu = 0.35f32       in -- []none
  let kdact = 0.256e-3f32 in -- []mM
  let q10NCX = 1.57f32   in -- []none
  let ibarSLCaP = 0.0673f32 in -- [#uA/uF]
  let kmPCa = 0.5e-3f32  in -- []mM
  let gCaB = 2.513e-4f32 in -- [#uA/uF]
  let q10SLCaP = 2.35f32 in -- []none

  -- SR flux parameters
  let q10SRCaP = 2.6f32 in -- []none
  let vmax_SRCaP = 2.86e-4f32 in -- [#mM/msec] (mmol/L cytosol/msec)
  let kmf = 0.246e-3f32    in -- []mM
  let kmr = 1.7f32         in -- []mML cytosol
  let hillSRCaP = 1.787f32 in -- []mM
  let ks   = 25.0f32 in -- [1/ms]
  let koCa = 10.0f32 in -- [#mM^-2 1/ms]
  let kom  = 0.06f32 in -- [1/ms]
  let kiCa = 0.5f32  in -- [1/mM/ms]
  let kim = 0.005f32 in -- [1/ms]
  let ec50SR = 0.45f32  in -- []mM

  -- Buffering parameters
  let bmax_Naj = 7.561f32 in -- []mM
  let bmax_Nasl = 1.65f32 in -- []mM
  let koff_na = 1.0e-3f32   in -- [1/ms]
  let kon_na = 0.1e-3f32  in -- [1/mM/ms]
  let bmax_Tnclow = 70e-3f32 in -- []mM, TnC low affinity
  let koff_tncl = 19.6e-3f32 in -- [1/ms]
  let kon_tncl = 32.7f32     in -- [1/mM/ms]
  let bmax_TnChigh = 140.0e-3f32  in -- []mM, TnC high affinity
  let koff_tnchca = 0.032e-3f32 in -- [1/ms]
  let kon_tnchca = 2.37f32      in -- [1/mM/ms]
  let koff_tnchmg = 3.33e-3f32  in -- [1/ms]
  let kon_tnchmg = 3.0e-3f32      in -- [1/mM/ms]
  let bmax_myosin = 140.0e-3f32   in -- []mM, Myosin buffering
  let koff_myoca = 0.46e-3f32     in -- [1/ms]
  let kon_myoca = 13.8f32         in -- [1/mM/ms]
  let koff_myomg = 0.057e-3f32    in -- [1/ms]
  let kon_myomg = 0.0157f32       in -- [1/mM/ms]
  let bmax_SR = 19.0f32 * 0.9e-3f32 in -- []mM
  let koff_sr = 60.0e-3f32        in -- [1/ms]
  let kon_sr = 100.0f32           in -- [1/mM/ms]
  let bmax_SLlowsl = 37.38e-3f32 * vmyo / vsl in -- []mM, SL buffering
  let bmax_SLlowj = 4.62e-3f32   * vmyo / vjunc * 0.1f32 in -- []mM
  let koff_sll = 1300.0e-3f32     in -- [1/ms]
  let kon_sll = 100.0f32          in -- [1/mM/ms]
  let bmax_SLhighsl = 13.35e-3f32 * vmyo / vsl in -- []mM
  let bmax_SLhighj = 1.65e-3f32 * vmyo / vjunc * 0.1f32 in -- []mM
  let koff_slh = 30e-3f32 in -- [1/ms]
  let kon_slh = 100.0f32  in -- [1/mM/ms]
  let bmax_Csqn = 2.7f32  in -- 140e-3*vmyo/Vsr; []mM
  let koff_csqn = 65.0f32 in -- [1/ms]
  let kon_csqn = 100.0f32 in -- [1/mM/ms]

  -- i_Na: Fast Na Current
  let am = 0.32f32 * (initvalu_39 + 47.13f32) / (1.0f32 - f32.exp(-0.1f32*(initvalu_39+47.13f32)))
  let bm = 0.08f32 * f32.exp(-initvalu_39/11.0f32)
  let ( ah, aj, bh, bj ) =  if (initvalu_39 >= -40.0f32)
                            then ( 0.0f32,
                                   0.0f32,
                                   1.0f32/(0.13f32*(1.0f32+f32.exp(-(initvalu_39+10.66f32)/11.1f32))),
                                   0.3f32*f32.exp(-2.535e-7f32*initvalu_39)/(1.0f32+f32.exp(-0.1f32*(initvalu_39+32.0f32)))
                                 )
                            else ( 0.135f32*f32.exp((80.0f32+initvalu_39)/(-6.8f32)),
                                   ( -127140.0f32*f32.exp(0.2444f32*initvalu_39) - 3.474e-5f32*f32.exp(-0.04391f32*initvalu_39) ) *
                                        (initvalu_39+37.78f32)/(1.0f32+f32.exp(0.311f32*(initvalu_39+79.23f32))),
                                   3.56f32*f32.exp(0.079f32*initvalu_39)+3.1e5f32*f32.exp(0.35f32*initvalu_39),
                                   0.1212f32 * f32.exp(-0.01052f32*initvalu_39) / (1.0f32 + f32.exp(-0.1378f32 * (initvalu_39+40.14f32)))
                                 )

  let finavalu[offset_1] = am*(1.0f32-initvalu_1) - bm*initvalu_1
  let finavalu[offset_2] = ah*(1.0f32-initvalu_2) - bh*initvalu_2
  let finavalu[offset_3] = aj*(1.0f32-initvalu_3) - bj*initvalu_3
  let i_Na_junc = fjunc*gna*pow(initvalu_1,3.0f32)*initvalu_2*initvalu_3*(initvalu_39-ena_junc)
  let i_Na_sl = fsl*gna*pow(initvalu_1,3.0f32)*initvalu_2*initvalu_3*(initvalu_39-ena_sl)

  -- i_nabk: Na Background Current
  let i_nabk_junc = fjunc*gnaB*(initvalu_39-ena_junc)
  let i_nabk_sl   = fsl  *gnaB*(initvalu_39-ena_sl  )

  -- i_nak: Na/K Pump Current
  let sigma = (f32.exp(nao/67.3f32)-1.0f32)/7.0f32
  let fnak = 1.0f32/(1.0f32+0.1245f32*f32.exp(-0.1f32*initvalu_39*foRT)+0.0365f32*sigma*f32.exp(-initvalu_39*foRT))
  let i_nak_junc = fjunc*iBarNaK*fnak*ko / (1.0f32+pow((kmNaip/initvalu_32),4.0f32)) /(ko+kmko)
  let i_nak_sl = fsl*iBarNaK*fnak*ko /(1.0f32+pow((kmNaip/initvalu_33),4.0f32)) /(ko+kmko)
  let i_nak = i_nak_junc + i_nak_sl

  -- i_kr: Rapidly Activating K Current
  let gkr = 0.03f32 * f32.sqrt(ko/5.4f32)
  let xrss = 1.0f32 / (1.0f32 + f32.exp(-(initvalu_39+50.0f32)/7.5f32))
  let tauxr = 1.0f32/(0.00138f32*(initvalu_39+7.0f32) / (1.0f32-f32.exp(-0.123f32*(initvalu_39+7.0f32))) +
                6.1e-4f32*(initvalu_39+10.0f32)/(f32.exp(0.145f32*(initvalu_39+10.0f32))-1.0f32))
  let finavalu[offset_12] = (xrss-initvalu_12)/tauxr
  let rkr = 1.0f32 / (1.0f32 + f32.exp((initvalu_39+33.0f32)/22.4f32))
  let i_kr = gkr*initvalu_12*rkr*(initvalu_39-ek)

  -- i_ks: Slowly Activating K Current
  let pcaks_junc = -log10(initvalu_36)+3.0f32
  let pcaks_sl = -log10(initvalu_37)+3.0f32
  let gks_junc = 0.07f32*(0.057f32 +0.19f32/(1.0f32+ f32.exp((-7.2f32+pcaks_junc)/0.6f32)))
  let gks_sl = 0.07f32*(0.057f32 +0.19f32/(1.0f32+ f32.exp((-7.2f32+pcaks_sl)/0.6f32)))
  let eks = (1.0f32/foRT)*f32.log((ko+pNaK*nao)/(initvalu_35+pNaK*initvalu_34))
  let xsss = 1.0f32/(1.0f32+f32.exp(-(initvalu_39-1.5f32)/16.7f32))
  let tauxs = 1.0f32/(7.19e-5f32*(initvalu_39+30.0f32)/(1.0f32-f32.exp(-0.148f32*(initvalu_39+30.0f32))) +
                1.31e-4f32*(initvalu_39+30.0f32)/(f32.exp(0.0687f32*(initvalu_39+30.0f32))-1.0f32))
  let finavalu[offset_13] = (xsss-initvalu_13) / tauxs
  let i_ks_junc = fjunc*gks_junc*pow(initvalu_12,2.0f32)*(initvalu_39-eks)
  let i_ks_sl = fsl*gks_sl*pow(initvalu_13,2.0f32)*(initvalu_39-eks)
  let i_ks = i_ks_junc+i_ks_sl

  -- i_kp: Plateau K current
  let kp_kp = 1.0f32/(1.0f32+f32.exp(7.488f32-initvalu_39/5.98f32))
  let i_kp_junc = fjunc*gkp*kp_kp*(initvalu_39-ek)
  let i_kp_sl = fsl*gkp*kp_kp*(initvalu_39-ek)
  let i_kp = i_kp_junc+i_kp_sl

  -- i_to: Transient Outward K Current (slow and fast components)
  let xtoss = 1.0f32/(1.0f32+f32.exp(-(initvalu_39+3.0f32)/15.0f32))
  let ytoss = 1.0f32/(1.0f32+f32.exp((initvalu_39+33.5f32)/10.0f32))
  let rtoss = 1.0f32/(1.0f32+f32.exp((initvalu_39+33.5f32)/10.0f32))
  let tauxtos = 9.0f32/(1.0f32+f32.exp((initvalu_39+3.0f32)/15.0f32))+0.5f32
  let tauytos = 3.0e3f32/(1.0f32+f32.exp((initvalu_39+60.0f32)/10.0f32))+30.0f32
  let taurtos = 2800.0f32/(1.0f32+f32.exp((initvalu_39+60.0f32)/10.0f32))+220.0f32
  let finavalu[offset_8] = (xtoss-initvalu_8)/tauxtos
  let finavalu[offset_9] = (ytoss-initvalu_9)/tauytos
  let finavalu[offset_40]= (rtoss-initvalu_40)/taurtos
  let i_tos = dToSlow*initvalu_8*(initvalu_9+0.5f32*initvalu_40)*(initvalu_39-ek) -- [#uA/uF]

  --
  let tauxtof = 3.5f32*f32.exp(-initvalu_39*initvalu_39/30.0f32/30.0f32)+1.5f32
  let tauytof = 20.0f32/(1.0f32+f32.exp((initvalu_39+33.5f32)/10.0f32))+20.0f32
  let finavalu[offset_10] = (xtoss-initvalu_10)/tauxtof
  let finavalu[offset_11] = (ytoss-initvalu_11)/tauytof
  let i_tof = gToFast*initvalu_10*initvalu_11*(initvalu_39-ek)
  let i_to = i_tos + i_tof

  -- i_ki: Time-independent K Current
  let aki = 1.02f32/(1.0f32+f32.exp(0.2385f32*(initvalu_39-ek-59.215f32)))
  let bki =(0.49124f32*f32.exp(0.08032f32*(initvalu_39+5.476f32-ek)) + f32.exp(0.06175f32*(initvalu_39-ek-594.31f32))) /
                (1.0f32 + f32.exp(-0.5143f32*(initvalu_39-ek+4.753f32)))
  let kiss = aki/(aki+bki)
  let i_ki = 0.9f32*f32.sqrt(ko/5.4f32)*kiss*(initvalu_39-ek)

  -- i_ClCa: Ca-activated Cl Current, i_Clbk: background Cl Current
  let i_ClCa_junc = fjunc*gClCa/(1.0f32+kdClCa/initvalu_36)*(initvalu_39-ecl)
  let i_ClCa_sl = fsl*gClCa/(1.0f32+kdClCa/initvalu_37)*(initvalu_39-ecl)
  let i_ClCa = i_ClCa_junc+i_ClCa_sl
  let i_Clbk = gClB*(initvalu_39-ecl)

  -- i_Ca: L-type Calcium Current
  let dss = 1.0f32/(1.0f32+f32.exp(-(initvalu_39+14.5f32)/6.0f32))
  let taud = dss*(1.0f32-f32.exp(-(initvalu_39+14.5f32)/6.0f32))/(0.035f32*(initvalu_39+14.5f32))
  let fss = 1.0f32/(1.0f32+f32.exp((initvalu_39+35.06f32)/3.6f32))+0.6f32/(1.0f32+f32.exp((50.0f32-initvalu_39)/20.0f32))
  let tauf = 1.0f32/(0.0197f32*f32.exp(-pow(0.0337f32*(initvalu_39+14.5f32),2.0f32))+0.02f32)
  let finavalu[offset_4] = (dss-initvalu_4)/taud
  let finavalu[offset_5] = (fss-initvalu_5)/tauf
  let finavalu[offset_6] = 1.7f32*initvalu_36*(1.0f32-initvalu_6)-11.9e-3f32*initvalu_6 in -- fCa_junc
  let finavalu[offset_7] = 1.7f32*initvalu_37*(1.0f32-initvalu_7)-11.9e-3f32*initvalu_7    -- fCa_sl

  --
  let ibarca_j = pCa*4.0f32*(initvalu_39*frdy*foRT) * (0.341f32*initvalu_36*f32.exp(2.0f32*initvalu_39*foRT)-0.341f32*cao) /
                    (f32.exp(2.0f32*initvalu_39*foRT)-1.0f32)
  let ibarca_sl= pCa*4.0f32*(initvalu_39*frdy*foRT) * (0.341f32*initvalu_37*f32.exp(2.0f32*initvalu_39*foRT)-0.341f32*cao) /
                    (f32.exp(2.0f32*initvalu_39*foRT)-1.0f32)
  let ibark = pK*(initvalu_39*frdy*foRT)*(0.75f32*initvalu_35*f32.exp(initvalu_39*foRT)-0.75f32*ko) /
                    (f32.exp(initvalu_39*foRT)-1.0f32)
  let ibarna_j = pNa*(initvalu_39*frdy*foRT) *(0.75f32*initvalu_32*f32.exp(initvalu_39*foRT)-0.75f32*nao) /
                    (f32.exp(initvalu_39*foRT)-1.0f32)
  let ibarna_sl= pNa*(initvalu_39*frdy*foRT) *(0.75f32*initvalu_33*f32.exp(initvalu_39*foRT)-0.75f32*nao) /
                    (f32.exp(initvalu_39*foRT)-1.0f32)
  let i_Ca_junc = (fjunc_CaL*ibarca_j*initvalu_4*initvalu_5*(1.0f32-initvalu_6)*pow(q10CaL,qpow))*0.45f32
  let i_Ca_sl   = (fsl_CaL*ibarca_sl*initvalu_4*initvalu_5*(1.0f32-initvalu_7)*pow(q10CaL,qpow))*0.45f32
  let i_Ca = i_Ca_junc+i_Ca_sl
  let finavalu[offset_43] = -i_Ca*cmem/(vmyo*2.0f32*frdy)*1e3f32
  let i_CaK = (ibark*initvalu_4*initvalu_5*(fjunc_CaL*(1.0f32-initvalu_6)+fsl_CaL*(1.0f32-initvalu_7))*pow(q10CaL,qpow))*0.45f32
  let i_CaNa_junc = (fjunc_CaL*ibarna_j*initvalu_4*initvalu_5*(1.0f32-initvalu_6)*pow(q10CaL,qpow))*0.45f32
  let i_CaNa_sl = (fsl_CaL*ibarna_sl*initvalu_4*initvalu_5*(1.0f32-initvalu_7)*pow(q10CaL,qpow))*0.45f32

  -- i_ncx: Na/Ca Exchanger flux
  let ka_junc = 1.0f32/(1.0f32+pow((kdact/initvalu_36),3.0f32))
  let ka_sl = 1.0f32/(1.0f32+pow((kdact/initvalu_37),3.0f32))
  let s1_junc = f32.exp(nu*initvalu_39*foRT)*pow(initvalu_32,3.0f32)*cao
  let s1_sl = f32.exp(nu*initvalu_39*foRT)*pow(initvalu_33,3.0f32)*cao
  let s2_junc = f32.exp((nu-1.0f32)*initvalu_39*foRT)*pow(nao,3.0f32)*initvalu_36
  let s3_junc = (kmCai*pow(nao,3.0f32)*(1.0f32+pow((initvalu_32/kmNai),3.0f32))+pow(kmnao,3.0f32)*initvalu_36 +
                    pow(kmNai,3.0f32)*cao*(1.0f32+initvalu_36/kmCai)+kmcao*pow(initvalu_32,3.0f32)+pow(initvalu_32,3.0f32)*cao +
                    pow(nao,3.0f32)*initvalu_36)*(1.0f32+ksat*f32.exp((nu-1.0f32)*initvalu_39*foRT))
  let s2_sl = f32.exp((nu-1.0f32)*initvalu_39*foRT)*pow(nao,3.0f32)*initvalu_37
  let s3_sl = (kmCai*pow(nao,3.0f32)*(1.0f32+pow((initvalu_33/kmNai),3.0f32)) +
                pow(kmnao,3.0f32)*initvalu_37+pow(kmNai,3.0f32)*cao*(1.0f32+initvalu_37/kmCai) +
                kmcao*pow(initvalu_33,3.0f32)+pow(initvalu_33,3.0f32)*cao+pow(nao,3.0f32)*initvalu_37)*
                (1.0f32+ksat*f32.exp((nu-1.0f32)*initvalu_39*foRT))
  let i_ncx_junc = fjunc*ibarNCX*pow(q10NCX,qpow)*ka_junc*(s1_junc-s2_junc)/s3_junc
  let i_ncx_sl = fsl*ibarNCX*pow(q10NCX,qpow)*ka_sl*(s1_sl-s2_sl)/s3_sl
  let i_ncx = i_ncx_junc+i_ncx_sl
  let finavalu[offset_45] = 2.0f32*i_ncx*cmem/(vmyo*2.0f32*frdy)*1e3f32

  -- i_pca: Sarcolemmal Ca Pump Current
  let i_pca_junc =  fjunc *
                    pow(q10SLCaP,qpow) *
                    ibarSLCaP *
                    pow(initvalu_36,1.6f32) /
                    (pow(kmPCa,1.6f32) + pow(initvalu_36,1.6f32))
  let i_pca_sl = fsl *
                pow(q10SLCaP,qpow) *
                ibarSLCaP *
                pow(initvalu_37,1.6f32) /
                (pow(kmPCa,1.6f32) + pow(initvalu_37,1.6f32))
  let i_pca = i_pca_junc + i_pca_sl
  let finavalu[offset_44] = -i_pca*cmem/(vmyo*2.0f32*frdy)*1e3f32

  -- i_cabk: Ca Background Current
  let i_cabk_junc = fjunc*gCaB*(initvalu_39-eca_junc)
  let i_cabk_sl = fsl*gCaB*(initvalu_39-eca_sl)
  let i_cabk = i_cabk_junc+i_cabk_sl
  let finavalu[offset_46] = -i_cabk*cmem/(vmyo*2.0f32*frdy)*1e3f32

  -- SR fluxes: Calcium Release, SR Ca pump, SR Ca leak
  let maxSR = 15.0f32
  let minSR = 1.0f32
  let kCaSR = maxSR - (maxSR-minSR)/(1.0f32+pow(ec50SR/initvalu_31,2.5f32))
  let koSRCa = koCa/kCaSR
  let kiSRCa = kiCa*kCaSR
  let ri = 1.0f32 - initvalu_14 - initvalu_15 - initvalu_16
  let finavalu[offset_14] = (kim*ri-kiSRCa*initvalu_36*initvalu_14) -
                            (koSRCa*pow(initvalu_36,2.0f32) * initvalu_14-kom*initvalu_15) in -- R
  let finavalu[offset_15] = (koSRCa*pow(initvalu_36,2.0f32) * initvalu_14-kom*initvalu_15) -
                            (kiSRCa*initvalu_36*initvalu_15-kim*initvalu_16) in -- O
  let finavalu[offset_16] = (kiSRCa*initvalu_36*initvalu_15-kim*initvalu_16) -
                            (kom*initvalu_16-koSRCa*pow(initvalu_36,2.0f32)*ri) in -- i
  let j_SRCarel = ks*initvalu_15*(initvalu_31-initvalu_36) in -- [#mM/ms]
  let j_serca = pow(q10SRCaP,qpow)*vmax_SRCaP*(pow((initvalu_38/kmf),hillSRCaP)-pow((initvalu_31/kmr),hillSRCaP))
                                        /(1.0f32+pow((initvalu_38/kmf),hillSRCaP)+pow((initvalu_31/kmr),hillSRCaP))
  let j_SRleak = 5.348e-6f32*(initvalu_31-initvalu_36) -- [#mM/ms]

  -- Sodium and Calcium Buffering
  let finavalu[offset_17] = kon_na*initvalu_32*(bmax_Naj-initvalu_17)-koff_na*initvalu_17   in -- NaBj  [#mM/ms]
  let finavalu[offset_18] = kon_na*initvalu_33*(bmax_Nasl-initvalu_18)-koff_na*initvalu_18  in -- NaBsl [#mM/ms]
  -- Cytosolic Ca Buffers
  let finavalu[offset_19] = kon_tncl*initvalu_38*(bmax_Tnclow-initvalu_19)-koff_tncl*initvalu_19 in -- TnCL [#mM/ms]
  let finavalu[offset_20] = kon_tnchca*initvalu_38*(bmax_TnChigh-initvalu_20-initvalu_21)-koff_tnchca*initvalu_20 in -- TnCHc [#mM/ms]
  let finavalu[offset_21] = kon_tnchmg*mgi*(bmax_TnChigh-initvalu_20-initvalu_21)-koff_tnchmg*initvalu_21 in -- TnCHm [#mM/ms]
  let finavalu[offset_22] = 0.0f32 in -- CaM [#mM/ms]
  let finavalu[offset_23] = kon_myoca*initvalu_38*(bmax_myosin-initvalu_23-initvalu_24)-koff_myoca*initvalu_23 in -- Myosin_ca [#mM/ms]
  let finavalu[offset_24] = kon_myomg*mgi*(bmax_myosin-initvalu_23-initvalu_24)-koff_myomg*initvalu_24 in -- Myosin_mg [#mM/ms]
  let finavalu[offset_25] = kon_sr*initvalu_38*(bmax_SR-initvalu_25)-koff_sr*initvalu_25 in -- SRB [#mM/ms]
  let j_CaB_cytosol =   finavalu[offset_19] + finavalu[offset_20] + finavalu[offset_21] +
                        finavalu[offset_22] + finavalu[offset_23] + finavalu[offset_24] + finavalu[offset_25]

  -- Junctional and SL Ca Buffers
  let finavalu[offset_26] = kon_sll*initvalu_36*(bmax_SLlowj-initvalu_26)-koff_sll*initvalu_26  in -- SLLj  [#mM/ms]
  let finavalu[offset_27] = kon_sll*initvalu_37*(bmax_SLlowsl-initvalu_27)-koff_sll*initvalu_27 in -- SLLsl [#mM/ms]
  let finavalu[offset_28] = kon_slh*initvalu_36*(bmax_SLhighj-initvalu_28)-koff_slh*initvalu_28 in -- SLHj  [#mM/ms]
  let finavalu[offset_29] = kon_slh*initvalu_37*(bmax_SLhighsl-initvalu_29)-koff_slh*initvalu_29 in-- SLHsl [#mM/ms]
  let j_CaB_junction = finavalu[offset_26]+finavalu[offset_28]
  let j_CaB_sl = finavalu[offset_27]+finavalu[offset_29]

  -- SR Ca Concentrations
  let finavalu[offset_30] = kon_csqn*initvalu_31*(bmax_Csqn-initvalu_30)-koff_csqn*initvalu_30 in -- Csqn [#mM/ms]
  let oneovervsr = 1.0f32 / vsr
  let finavalu[offset_31] = j_serca*vmyo*oneovervsr-(j_SRleak*vmyo*oneovervsr+j_SRCarel)-finavalu[offset_30] -- Ca_sr [#mM/ms]

  -- Sodium Concentrations
  let i_Na_tot_junc = i_Na_junc+i_nabk_junc+3.0f32*i_ncx_junc+3.0f32*i_nak_junc+i_CaNa_junc in -- [#uA/uF]
  let i_Na_tot_sl = i_Na_sl+i_nabk_sl+3.0f32*i_ncx_sl+3.0f32*i_nak_sl+i_CaNa_sl in -- [#uA/uF]
  let finavalu[offset_32] = -i_Na_tot_junc*cmem/(vjunc*frdy)+j_na_juncsl/vjunc*(initvalu_33-initvalu_32)-finavalu[offset_17]
  let oneovervsl = 1.0f32 / vsl
  let finavalu[offset_33] = -i_Na_tot_sl*cmem*oneovervsl/frdy+j_na_juncsl*oneovervsl*(initvalu_32-initvalu_33) +
                            j_na_slmyo*oneovervsl*(initvalu_34-initvalu_33)-finavalu[offset_18]
  let finavalu[offset_34] = j_na_slmyo/vmyo*(initvalu_33-initvalu_34) -- [#mM/msec]

  -- Potassium Concentration
  let i_K_tot = i_to+i_kr+i_ks+i_ki-2.0f32*i_nak+i_CaK+i_kp in -- [#uA/uF]
  let finavalu[offset_35] = 0.0f32 -- [#mM/msec]

  -- Calcium Concentrations
  let i_Ca_tot_junc = i_Ca_junc+i_cabk_junc+i_pca_junc-2.0f32*i_ncx_junc in -- [#uA/uF]
  let i_Ca_tot_sl = i_Ca_sl+i_cabk_sl+i_pca_sl-2.0f32*i_ncx_sl in -- [#uA/uF]
  let finavalu[offset_36] = -i_Ca_tot_junc*cmem/(vjunc*2.0f32*frdy)+
                            j_ca_juncsl/vjunc*(initvalu_37-initvalu_36) -
                            j_CaB_junction+(j_SRCarel)*vsr/vjunc+j_SRleak*vmyo/vjunc in -- Ca_j
  let finavalu[offset_37] = -i_Ca_tot_sl*cmem/(vsl*2.0f32*frdy) +
                            j_ca_juncsl/vsl*(initvalu_36-initvalu_37) +
                            j_ca_slmyo/vsl*(initvalu_38-initvalu_37)-j_CaB_sl in -- Ca_sl
  let finavalu[offset_38] = -j_serca-j_CaB_cytosol +j_ca_slmyo/vmyo*(initvalu_37-initvalu_38)

  let i_app = if(fmod(timeinst,parameter_1) <= 5.0f32)
              then 9.5f32 else 0.0f32

  -- Membrane Potential
  let i_Na_tot = i_Na_tot_junc + i_Na_tot_sl in -- [#uA/uF]
  let i_Cl_tot = i_ClCa+i_Clbk in -- [#uA/uF]
  let i_Ca_tot = i_Ca_tot_junc + i_Ca_tot_sl
  let i_tot = i_Na_tot+i_Cl_tot+i_Ca_tot+i_K_tot
  let finavalu[offset_39] = -(i_tot-i_app)

  -- Set unused output values to 0 (MATLAB does it by default)
  let finavalu[offset_41] = 0.0f32
  let finavalu[offset_42] = 0.0f32
  in finavalu


----------------------
----- Cam MODULE -----
----------------------

let cam [equs][pars] (_timeinst:  f32, initvalu: [equs]f32,
                      initvalu_offset: i32,
                      parameter: [pars]f32, parameter_offset: i32,
                      finavalu: *[equs]f32, ca: f32 ): (f32, *[equs]f32) =
  #[unsafe]

  -- input data and output data variable references
  let offset_1  = initvalu_offset
  let offset_2  = initvalu_offset+1
  let offset_3  = initvalu_offset+2
  let offset_4  = initvalu_offset+3
  let offset_5  = initvalu_offset+4
  let offset_6  = initvalu_offset+5
  let offset_7  = initvalu_offset+6
  let offset_8  = initvalu_offset+7
  let offset_9  = initvalu_offset+8
  let offset_10 = initvalu_offset+9
  let offset_11 = initvalu_offset+10
  let offset_12 = initvalu_offset+11
  let offset_13 = initvalu_offset+12
  let offset_14 = initvalu_offset+13
  let offset_15 = initvalu_offset+14

  -- input parameters variable references
  let parameter_offset_2  = parameter_offset+1
  let parameter_offset_3  = parameter_offset+2
  let parameter_offset_4  = parameter_offset+3
  let parameter_offset_5  = parameter_offset+4

  -- decoding input array
  let caM   = initvalu[offset_1]
  let ca2caM  = initvalu[offset_2]
  let ca4caM  = initvalu[offset_3]
  let caMB    = initvalu[offset_4]
  let ca2caMB = initvalu[offset_5]
  let ca4caMB = initvalu[offset_6]
  let pb2     = initvalu[offset_7]
  let pb      = initvalu[offset_8]
  let pt      = initvalu[offset_9]
  let pt2     = initvalu[offset_10]
  let pa      = initvalu[offset_11]
  let ca4caN  = initvalu[offset_12]
  let caMca4caN = initvalu[offset_13]
  let ca2caMca4caN = initvalu[offset_14]
  let ca4caMca4caN = initvalu[offset_15]

  -- decoding input parameters
  let btot      = parameter[parameter_offset_2]
  let caMKiitot = parameter[parameter_offset_3]
  let caNtot    = parameter[parameter_offset_4]
  let pp1tot    = parameter[parameter_offset_5]

  -- values [CONSTANTS FOR ALL THREADS]
  let k = 135.0f32
  let mg = 1.0f32

  -- ca/caM parameters
  let (kd02, kd24) = if (mg <= 1.0f32)
                     then ( 0.0025f32*(1.0f32+k/0.94f32-mg/0.012f32)*(1.0f32+k/8.1f32+mg/0.022f32),
                            0.128f32*(1.0f32+k/0.64f32+mg/0.0014f32)*(1.0f32+k/13.0f32-mg/0.153f32)
                          )
                     else ( 0.0025f32*(1.0f32+k/0.94f32-1.0f32/0.012f32+(mg-1.0f32)/0.060f32) *
                                (1.0f32+k/8.1f32+1.0f32/0.022f32+(mg-1.0f32)/0.068f32),
                            0.128f32*(1.0f32+k/0.64f32+1.0f32/0.0014f32+(mg-1.0f32)/0.005f32) *
                                (1.0f32+k/13.0f32-1.0f32/0.153f32+(mg-1.0f32)/0.150f32)
                          )

  let k20 = 10.0f32  in -- [#s^-1]
  let k02 = k20/kd02 in -- [#uM^-2 s^-1]
  let k42 = 500.0f32 in -- [#s^-1]
  let k24 = k42/kd24    -- [#uM^-2 s^-1]

  -- caM buffering (B) parameters
  let k0Boff = 0.0014f32 in -- [#s^-1]
  let k0Bon = k0Boff/0.2f32 in -- [#uM^-1 s^-1] kon = koff/kd
  let k2Boff = k0Boff/100.0f32 in -- [#s^-1]
  let k2Bon = k0Bon in -- [#uM^-1 s^-1]
  let k4Boff = k2Boff in -- [#s^-1]
  let k4Bon = k0Bon      -- [#uM^-1 s^-1]

  -- using thermodynamic constraints
  let k20B = k20/100.0f32 in -- [#s^-1] thermo constraint on loop 1
  let k02B = k02 in -- [#uM^-2 s^-1]
  let k42B = k42 in -- [#s^-1] thermo constraint on loop 2
  let k24B = k24    -- [#uM^-2 s^-1]

  -- Wi Wa Wt Wp
  let kbi = 2.2f32 in -- [#s^-1] (ca4caM dissocation from Wb)
  let kib = kbi/33.5e-3f32 in -- [#uM^-1 s^-1]
  let kpp1 = 1.72f32 in -- [#s^-1] (pp1-dep dephosphorylation rates)
  let kmpp1 = 11.5f32 in -- []uM
  let kib2 = kib
  let kb2i = kib2*5.0f32
  let kb24 = k24
  let kb42 = k42*33.5e-3f32/5.0f32
  let kta = kbi/1000.0f32 in -- [#s^-1] (ca4caM dissociation from Wt)
  let kat = kib in -- [#uM^-1 s^-1] (ca4caM reassociation with Wa)
  let kt42 = k42*33.5e-6f32/5.0f32
  let kt24 = k24
  let kat2 = kib
  let kt2a = kib*5.0f32

  -- caN parameters
  let kcancaoff = 1.0f32 in -- [#s^-1]
  let kcancaon = kcancaoff/0.5f32 in -- [#uM^-1 s^-1]
  let kcancaM4on = 46.0f32 in -- [#uM^-1 s^-1]
  let kcancaM4off = 0.0013f32 in -- [#s^-1]
  let kcancaM2on = kcancaM4on
  let kcancaM2off = 2508.0f32*kcancaM4off
  let kcancaM0on = kcancaM4on
  let kcancaM0off = 165.0f32*kcancaM2off
  let k02can = k02
  let k20can = k20/165.0f32
  let k24can = k24
  let k42can = k20/2508.0f32

  -- caM Reaction fluxes
  let rcn02 = k02*pow(ca,2.0f32)*caM - k20*ca2caM
  let rcn24 = k24*pow(ca,2.0f32)*ca2caM - k42*ca4caM

  -- caM buffer fluxes
  let b = btot - caMB - ca2caMB - ca4caMB
  let rcn02B = k02B*pow(ca,2.0f32)*caMB - k20B*ca2caMB
  let rcn24B = k24B*pow(ca,2.0f32)*ca2caMB - k42B*ca4caMB
  let rcn0B = k0Bon*caM*b - k0Boff*caMB
  let rcn2B = k2Bon*ca2caM*b - k2Boff*ca2caMB
  let rcn4B = k4Bon*ca4caM*b - k4Boff*ca4caMB

  -- caN reaction fluxes
  let ca2caN = caNtot - ca4caN - caMca4caN - ca2caMca4caN - ca4caMca4caN
  let rcnca4caN = kcancaon*pow(ca,2.0f32)*ca2caN - kcancaoff*ca4caN
  let rcn02caN = k02can*pow(ca,2.0f32)*caMca4caN - k20can*ca2caMca4caN
  let rcn24caN = k24can*pow(ca,2.0f32)*ca2caMca4caN - k42can*ca4caMca4caN
  let rcn0caN = kcancaM0on*caM*ca4caN - kcancaM0off*caMca4caN
  let rcn2caN = kcancaM2on*ca2caM*ca4caN - kcancaM2off*ca2caMca4caN
  let rcn4caN = kcancaM4on*ca4caM*ca4caN - kcancaM4off*ca4caMca4caN

  -- caMkii reaction fluxes
  let pix = 1.0f32 - pb2 - pb - pt - pt2 - pa
  let rcnCkib2 = kib2*ca2caM*pix - kb2i*pb2
  let rcnCkb2b = kb24*pow(ca,2.0f32)*pb2 - kb42*pb
  let rcnCkib = kib*ca4caM*pix - kbi*pb
  let t = pb + pt + pt2 + pa
  let kbt = 0.055f32 * t + 0.0074f32 * pow(t,2.0f32) + 0.015f32*pow(t,3.0f32)
  let rcnCkbt = kbt*pb - kpp1*pp1tot*pt/(kmpp1+caMKiitot*pt)
  let rcnCktt2 = kt42*pt - kt24*pow(ca,2.0f32)*pt2
  let rcnCkta = kta*pt - kat*ca4caM*pa
  let rcnCkt2a = kt2a*pt2 - kat2*ca2caM*pa
  let rcnCkt2b2 = kpp1*pp1tot*pt2/(kmpp1+caMKiitot*pt2)
  let rcnCkai = kpp1*pp1tot*pa/(kmpp1+caMKiitot*pa)

  -- caM equations
  let dcaM = 1e-3f32*(-rcn02 - rcn0B - rcn0caN)
  let dca2caM = 1e-3f32*(rcn02 - rcn24 - rcn2B - rcn2caN + caMKiitot*(-rcnCkib2 + rcnCkt2a) )
  let dca4caM = 1e-3f32*(rcn24 - rcn4B - rcn4caN + caMKiitot*(-rcnCkib+rcnCkta) )
  let dcaMB = 1e-3f32*(rcn0B-rcn02B)
  let dca2caMB = 1e-3f32*(rcn02B + rcn2B - rcn24B)
  let dca4caMB = 1e-3f32*(rcn24B + rcn4B)

  -- caMkii equations
  let dpb2 = 1e-3f32*(rcnCkib2 - rcnCkb2b + rcnCkt2b2) in -- pb2
  let dpb = 1e-3f32*(rcnCkib + rcnCkb2b - rcnCkbt) in -- pb
  let dpt = 1e-3f32*(rcnCkbt-rcnCkta-rcnCktt2) in -- pt
  let dpt2 = 1e-3f32*(rcnCktt2-rcnCkt2a-rcnCkt2b2) in -- pt2
  let dpa = 1e-3f32*(rcnCkta+rcnCkt2a-rcnCkai)        -- pa

  -- caN equations
  let dca4caN = 1e-3f32*(rcnca4caN - rcn0caN - rcn2caN - rcn4caN) in -- ca4caN
  let dcaMca4caN = 1e-3f32*(rcn0caN - rcn02caN) in           -- caMca4caN
  let dca2caMca4caN = 1e-3f32*(rcn2caN+rcn02caN-rcn24caN) in -- ca2caMca4caN
  let dca4caMca4caN = 1e-3f32*(rcn4caN+rcn24caN)             -- ca4caMca4caN

  -- encode output array
  let finavalu[offset_1] = dcaM
  let finavalu[offset_2] = dca2caM
  let finavalu[offset_3] = dca4caM
  let finavalu[offset_4] = dcaMB
  let finavalu[offset_5] = dca2caMB
  let finavalu[offset_6] = dca4caMB
  let finavalu[offset_7] = dpb2
  let finavalu[offset_8] = dpb
  let finavalu[offset_9] = dpt
  let finavalu[offset_10] = dpt2
  let finavalu[offset_11] = dpa
  let finavalu[offset_12] = dca4caN
  let finavalu[offset_13] = dcaMca4caN
  let finavalu[offset_14] = dca2caMca4caN
  let finavalu[offset_15] = dca4caMca4caN

  -- write to global variables for adjusting ca buffering in EC coupling model
  let jca = 1e-3f32*(2.0f32*caMKiitot*(rcnCktt2-rcnCkb2b) -
            2.0f32*(rcn02+rcn24+rcn02B+rcn24B+rcnca4caN+rcn02caN+rcn24caN)) -- [#uM/msec]
  in ( jca, finavalu)

----------------------
----- Fin MODULE -----
----------------------

let fin [equs][pars]
       (initvalu:     [equs]f32, initvalu_offset_ecc: i32,
        initvalu_offset_Dyad: i32, initvalu_offset_SL: i32,
        initvalu_offset_Cyt: i32, parameter: [pars]f32,
        finavalu: *[equs]f32, jcaDyad: f32, jcaSL: f32, jcaCyt: f32 ): *[equs]f32 =
  #[unsafe]

  let btotDyad      = parameter[2]
  let caMKiitotDyad = parameter[3]

  -- set variables
  let vmyo = 2.1454e-11f32
  let vdyad = 1.7790e-14f32
  let vsl = 6.6013e-13f32
  let kSLmyo = 8.587e-15f32
  let k0Boff = 0.0014f32
  let k0Bon = k0Boff/0.2f32
  let k2Boff = k0Boff/100.0f32
  let k2Bon = k0Bon
  let k4Bon = k0Bon

  -- ADJUST ECC incorporate ca buffering from caM, convert jcaCyt from uM/msec to mM/msec
  let finavalu[initvalu_offset_ecc+35] = finavalu[initvalu_offset_ecc+35] + 1e-3f32*jcaDyad
  let finavalu[initvalu_offset_ecc+36] = finavalu[initvalu_offset_ecc+36] + 1e-3f32*jcaSL
  let finavalu[initvalu_offset_ecc+37] = finavalu[initvalu_offset_ecc+37] + 1e-3f32*jcaCyt

  -- incorporate caM diffusion between compartments
  let caMtotDyad =    initvalu[initvalu_offset_Dyad+0]
                    + initvalu[initvalu_offset_Dyad+1]
                    + initvalu[initvalu_offset_Dyad+2]
                    + initvalu[initvalu_offset_Dyad+3]
                    + initvalu[initvalu_offset_Dyad+4]
                    + initvalu[initvalu_offset_Dyad+5]
                    + caMKiitotDyad * (  initvalu[initvalu_offset_Dyad+6]
                                    + initvalu[initvalu_offset_Dyad+7]
                                    + initvalu[initvalu_offset_Dyad+8]
                                    + initvalu[initvalu_offset_Dyad+9])
                    + initvalu[initvalu_offset_Dyad+12]
                    + initvalu[initvalu_offset_Dyad+13]
                    + initvalu[initvalu_offset_Dyad+14]

  let bdyad = btotDyad - caMtotDyad in -- [#uM dyad]
  let j_cam_dyadSL = 1e-3f32 * (  k0Boff*initvalu[initvalu_offset_Dyad+0] -
                        k0Bon*bdyad*initvalu[initvalu_offset_SL+0]) in -- [#uM/msec dyad]
  let j_ca2cam_dyadSL = 1e-3f32 * (  k2Boff*initvalu[initvalu_offset_Dyad+1] -
                        k2Bon*bdyad*initvalu[initvalu_offset_SL+1]) in -- [#uM/msec dyad]
  let j_ca4cam_dyadSL = 1e-3f32 * (  k2Boff*initvalu[initvalu_offset_Dyad+2] -
                        k4Bon*bdyad*initvalu[initvalu_offset_SL+2])    -- [#uM/msec dyad]

  let j_cam_SLmyo    = kSLmyo * (  initvalu[initvalu_offset_SL+0] - initvalu[initvalu_offset_Cyt+0]) in -- [#umol/msec]
  let j_ca2cam_SLmyo = kSLmyo * (  initvalu[initvalu_offset_SL+1] - initvalu[initvalu_offset_Cyt+1]) in -- [#umol/msec]
  let j_ca4cam_SLmyo = kSLmyo * (  initvalu[initvalu_offset_SL+2] - initvalu[initvalu_offset_Cyt+2])    -- [#umol/msec]

  -- ADJUST CAM Dyad
  let finavalu[initvalu_offset_Dyad+0] = finavalu[initvalu_offset_Dyad+0] - j_cam_dyadSL
  let finavalu[initvalu_offset_Dyad+1] = finavalu[initvalu_offset_Dyad+1] - j_ca2cam_dyadSL
  let finavalu[initvalu_offset_Dyad+2] = finavalu[initvalu_offset_Dyad+2] - j_ca4cam_dyadSL

  -- ADJUST CAM Sl
  let finavalu[initvalu_offset_SL+0] = finavalu[initvalu_offset_SL+0] + j_cam_dyadSL*vdyad/vsl - j_cam_SLmyo/vsl
  let finavalu[initvalu_offset_SL+1] = finavalu[initvalu_offset_SL+1] + j_ca2cam_dyadSL*vdyad/vsl - j_ca2cam_SLmyo/vsl
  let finavalu[initvalu_offset_SL+2] = finavalu[initvalu_offset_SL+2] + j_ca4cam_dyadSL*vdyad/vsl - j_ca4cam_SLmyo/vsl

  -- ADJUST CAM Cyt
  let finavalu[initvalu_offset_Cyt+0] = finavalu[initvalu_offset_Cyt+0] + j_cam_SLmyo/vmyo
  let finavalu[initvalu_offset_Cyt+1] = finavalu[initvalu_offset_Cyt+1] + j_ca2cam_SLmyo/vmyo
  let finavalu[initvalu_offset_Cyt+2] = finavalu[initvalu_offset_Cyt+2] + j_ca4cam_SLmyo/vmyo
  in finavalu

-------------------------
----- Master MODULE -----
-------------------------

let master [equs][pars] (timeinst:  f32, initvalu: [equs]f32, parameter: [pars]f32 ): *[equs]f32 =
  let finavalu = replicate equs 0.0f32
  -- ecc function
  let initvalu_offset_ecc  = 0
  let parameter_offset_ecc = 0
  let finavalu = ecc( timeinst, initvalu, initvalu_offset_ecc, parameter, parameter_offset_ecc, finavalu)

  let jcaDyad = 0.0f32
  let jcaSL   = 0.0f32
  let jcaCyt  = 0.0f32
  let (jcaDyad, jcaSL, jcaCyt, finavalu) =
    loop ((jcaDyad, jcaSL, jcaCyt, finavalu)) for ii < 3i32 do
        let (initvalu_offset, parameter_offset, ind) =
            if      (ii == 0) then -- cam function for Dyad
                ( 46, 1, 35 )
            else if (ii == 1) then -- cam function for SL
                ( 61, 6, 36 )
            else -- if (ii == 2) then -- cam function for Cyt
                ( 76, 11, 37 )

        let inp_val = #[unsafe] (initvalu[ind]*1e3f32)
        let (res_val, finavalu) = cam(  timeinst, initvalu, initvalu_offset, parameter,
                                        parameter_offset, finavalu, inp_val )

        let ( jcaDyad, jcaSL, jcaCyt ) =
            if      (ii == 0) then ( res_val, jcaSL,   jcaCyt )
            else if (ii == 1) then ( jcaDyad, res_val, jcaCyt )
            else                   ( jcaDyad, jcaSL,   res_val)
        in  ( jcaDyad, jcaSL, jcaCyt, finavalu)

  -- final adjustments
  let finavalu = fin(  initvalu, initvalu_offset_ecc, 46, 61, 76,
                       parameter, finavalu, jcaDyad, jcaSL, jcaCyt )

  in map (\(x: f32): f32  ->
            if ( f32.isnan(x) || f32.isinf(x) )
            then 0.0001f32
            else x
         ) finavalu


----------------------------------------
----- embedded_fehlberg_7_8 MODULE -----
----------------------------------------

let embedded_fehlberg_7_8 [equs][pars]
                         (timeinst:   f32, h: f32,
                          initvalu: [equs]f32,
                          parameter: [pars]f32 ): (*[equs]f32, *[equs]f32) =
  let c_1_11 = 41.0f32 / 840.0f32
  let c6 = 34.0f32 / 105.0f32
  let c_7_8= 9.0f32 / 35.0f32
  let c_9_10 = 9.0f32 / 280.0f32

  let a2 = 2.0f32 / 27.0f32
  let a3 = 1.0f32 / 9.0f32
  let a4 = 1.0f32 / 6.0f32
  let a5 = 5.0f32 / 12.0f32
  let a6 = 1.0f32 / 2.0f32
  let a7 = 5.0f32 / 6.0f32
  let a8 = 1.0f32 / 6.0f32
  let a9 = 2.0f32 / 3.0f32
  let a10 = 1.0f32 / 3.0f32

  let b31 = 1.0f32 / 36.0f32
  let b32 = 3.0f32 / 36.0f32
  let b41 = 1.0f32 / 24.0f32
  let b43 = 3.0f32 / 24.0f32
  let b51 = 20.0f32 / 48.0f32
  let b53 = -75.0f32 / 48.0f32
  let b54 = 75.0f32 / 48.0f32
  let b61 = 1.0f32 / 20.0f32
  let b64 = 5.0f32 / 20.0f32
  let b65 = 4.0f32 / 20.0f32
  let b71 = -25.0f32 / 108.0f32
  let b74 =  125.0f32 / 108.0f32
  let b75 = -260.0f32 / 108.0f32
  let b76 =  250.0f32 / 108.0f32
  let b81 = 31.0f32/300.0f32
  let b85 = 61.0f32/225.0f32
  let b86 = -2.0f32/9.0f32
  let b87 = 13.0f32/900.0f32
  let b91 = 2.0f32
  let b94 = -53.0f32/6.0f32
  let b95 = 704.0f32 / 45.0f32
  let b96 = -107.0f32 / 9.0f32
  let b97 = 67.0f32 / 90.0f32
  let b98 = 3.0f32
  let b10_1 = -91.0f32 / 108.0f32
  let b10_4 = 23.0f32 / 108.0f32
  let b10_5 = -976.0f32 / 135.0f32
  let b10_6 = 311.0f32 / 54.0f32
  let b10_7 = -19.0f32 / 60.0f32
  let b10_8 = 17.0f32 / 6.0f32
  let b10_9 = -1.0f32 / 12.0f32
  let b11_1 = 2383.0f32 / 4100.0f32
  let b11_4 = -341.0f32 / 164.0f32
  let b11_5 = 4496.0f32 / 1025.0f32
  let b11_6 = -301.0f32 / 82.0f32
  let b11_7 = 2133.0f32 / 4100.0f32
  let b11_8 = 45.0f32 / 82.0f32
  let b11_9 = 45.0f32 / 164.0f32
  let b11_10 = 18.0f32 / 41.0f32
  let b12_1 = 3.0f32 / 205.0f32
  let b12_6 = - 6.0f32 / 41.0f32
  let b12_7 = - 3.0f32 / 205.0f32
  let b12_8 = - 3.0f32 / 41.0f32
  let b12_9 = 3.0f32 / 41.0f32
  let b12_10 = 6.0f32 / 41.0f32
  let b13_1 = -1777.0f32 / 4100.0f32
  let b13_4 = -341.0f32 / 164.0f32
  let b13_5 = 4496.0f32 / 1025.0f32
  let b13_6 = -289.0f32 / 82.0f32
  let b13_7 = 2193.0f32 / 4100.0f32
  let b13_8 = 51.0f32 / 82.0f32
  let b13_9 = 33.0f32 / 164.0f32
  let b13_10 = 12.0f32 / 41.0f32

  let err_factor  = -41.0f32 / 840.0f32
  let h2_7 = a2 * h

  -- initvalu_temp[equations]
  -- finavalu_temp[13, equations]
  let finavalu_temp = replicate 13 (replicate equs 0.0f32)

  let finavalu_temp = loop (finavalu_temp) for ii < 13 do
    let initvalu_temp =
      map (\i -> if ii == 0 then initvalu[i]

                 else if ii == 1
                 then let x = initvalu[i]
                      let y = finavalu_temp[0, i]
                      in x + h2_7 * y

                 else if ii == 2
                 then let x = initvalu[i]
                      let y1 = finavalu_temp[0,i]
                      let y2 = finavalu_temp[1,i]
                      in x + h * ( b31*y1 + b32*y2)

                 else if ii == 3
                 then let x = initvalu[i]
                      let y1 = finavalu_temp[0,i]
                      let y2 = finavalu_temp[2,i]
                      in x + h * ( b41*y1 + b43*y2 )

                 else if ii == 4
                 then initvalu[i] + h * ( b51*finavalu_temp[0,i] +
                                          b53*finavalu_temp[2,i] + b54*finavalu_temp[3,i])

                 else if ii == 5
                 then initvalu[i] + h * ( b61*finavalu_temp[0,i] +
                                          b64*finavalu_temp[3,i] + b65*finavalu_temp[4,i] )

                 else if ii == 6
                 then initvalu[i] + h * ( b71*finavalu_temp[0,i] +
                                          b74*finavalu_temp[3,i] + b75*finavalu_temp[4,i] +
                                          b76*finavalu_temp[5,i] )

                 else if ii == 7
                 then initvalu[i] + h * ( b81*finavalu_temp[0,i] +
                                          b85*finavalu_temp[4,i] + b86*finavalu_temp[5,i] +
                                          b87*finavalu_temp[6,i] )

                 else if ii == 8
                 then initvalu[i] + h * ( b91*finavalu_temp[0,i] +
                                          b94*finavalu_temp[3,i] + b95*finavalu_temp[4,i] +
                                          b96*finavalu_temp[5,i] + b97*finavalu_temp[6,i] +
                                          b98*finavalu_temp[7,i] )

                 else if ii == 9
                 then initvalu[i] + h * ( b10_1*finavalu_temp[0,i] +
                                          b10_4*finavalu_temp[3,i] + b10_5*finavalu_temp[4,i] +
                                          b10_6*finavalu_temp[5,i] + b10_7*finavalu_temp[6,i] +
                                          b10_8*finavalu_temp[7,i] + b10_9*finavalu_temp[8,i] )

                 else if ii == 10
                 then initvalu[i] + h * ( b11_1*finavalu_temp[0,i] +
                                          b11_4*finavalu_temp[3,i] + b11_5*finavalu_temp[4,i] +
                                          b11_6*finavalu_temp[5,i] + b11_7*finavalu_temp[6,i] +
                                          b11_8*finavalu_temp[7,i] + b11_9*finavalu_temp[8,i] +
                                          b11_10 * finavalu_temp[9,i])

                 else if ii == 11
                 then initvalu[i] + h * ( b12_1*finavalu_temp[0,i] +
                                          b12_6*finavalu_temp[5,i] + b12_7*finavalu_temp[6,i] +
                                          b12_8*finavalu_temp[7,i] + b12_9*finavalu_temp[8,i] +
                                          b12_10 * finavalu_temp[9,i] )

                 else -- if (ii == 12) then
                   initvalu[i] +
                   h * ( b13_1*finavalu_temp[0,i] +
                         b13_4*finavalu_temp[3,i] + b13_5*finavalu_temp[4,i] +
                         b13_6*finavalu_temp[5,i] + b13_7*finavalu_temp[6,i] +
                         b13_8*finavalu_temp[7,i] + b13_9*finavalu_temp[8,i] +
                         b13_10*finavalu_temp[9,i] + finavalu_temp[11,i]))
           (iota equs)
    let timeinst_temp =
        if      (ii == 0) then timeinst
        else if (ii == 1) then timeinst + h2_7
        else if (ii == 2) then timeinst + a3*h
        else if (ii == 3) then timeinst + a4*h
        else if (ii == 4) then timeinst + a5*h
        else if (ii == 5) then timeinst + a6*h
        else if (ii == 6) then timeinst + a7 * h
        else if (ii == 7) then timeinst + a8 * h
        else if (ii == 8) then timeinst + a9*h
        else if (ii == 9) then timeinst + a10*h
        else if (ii == 10) then timeinst + h
        else if (ii == 11) then timeinst
        else -- if (ii == 12) then
          timeinst + h
    in #[unsafe] let finavalu_temp[ii] = master( timeinst_temp, initvalu_temp, parameter )
              in finavalu_temp
  --------------
  -- end loop --
  --------------

  let finavalu = map (\i  ->
                        initvalu[i] +  h * (c_1_11 * (finavalu_temp[0,i] + finavalu_temp[10,i]) +
                            c6 * finavalu_temp[5,i] + c_7_8 * (finavalu_temp[6,i] + finavalu_temp[7,i]) +
                            c_9_10 * (finavalu_temp[8,i] + finavalu_temp[9,i]) )
                    ) (iota(equs) )

  let error = map (\i  ->
                        fabs(err_factor * (finavalu_temp[0,i] + finavalu_temp[10,i] - finavalu_temp[11,i] - finavalu_temp[12,i]))
                 ) (iota(equs) )
  in ( finavalu, error )


-------------------------
----- SOLVER MODULE -----
-------------------------

let min_scale_factor(): f32 = 0.125f32
let max_scale_factor(): f32 = 4.0f32
let attempts(): i32         = 12

let solver [pars][equs] (xmax: i32, params: [pars]f32, y0: [equs]f32): (bool,[equs]f32) =
  let err_exponent  = 1.0f32 / 7.0f32
  let h_init = 1.0f32
  let h = h_init
  let xmin = 0
  let tolerance = 10.0f32 / f32.i32 (xmax-xmin)
  let y_km1  = y0 in

  if xmax < xmin || h <= 0.0f32 then (false, y_km1)
  else if xmax == xmin then (true, y_km1)
  else

  -- initialize return and loop-variant params
  let failed = false
  let km1    = 0
  let (_, failed, y_km1) =
    loop((km1, failed, y_km1)) while ( (!failed) && (km1 < xmax) ) do -- for km1 < xmax do
    -- reinitialize variables
    let h          = h_init
    let scale_fina = 1.0f32
    let y_k = replicate equs 0.0f32

    -- make attempts to minimize error
    let breakLoop  = false
    let j = 0
    let (_,_,y_k,breakLoop,_) =
      loop (j,h,_,breakLoop,_) =
           (j,h,y_k,breakLoop,scale_fina) while ( (!breakLoop) && (j < attempts()) ) do
      -- EVALUATE ALL equations
      let (y_k: [equs]f32, err: [equs]f32) = embedded_fehlberg_7_8( f32.i32 (km1), h, y_km1, params)

      -- iF THERE WAS NO ERROR FOR ANY OF equations, SET SCALE AND LEAVE THE LOOP
      let errs = map (\(e: f32): bool  -> if e > 0.0f32 then true else false) err
      let error= or errs

      in
      if (!error)
      then (j+1, h, y_k, true, max_scale_factor())
      else
      -- FiGURE OUT SCALE AS THE MiNiMUM OF COMPONENT SCALES
      let yy = map (\(x: f32): f32  ->
                        if (x == 0.0f32) then tolerance else fabs(x)
                  ) (y_km1 )

      let scale = map  (\(yy_erri: (f32,f32)): f32  ->
                          let (yyi, erri) = yy_erri
                          in 0.8f32 * pow( tolerance * yyi / erri , err_exponent )
                      ) (zip yy err )

      let scale_min = f32.minimum scale
      let scale_fina = f32.min(f32.max scale_min (min_scale_factor())) (max_scale_factor())

      -- iF WiTHiN TOLERANCE, FiNiSH attempts...
      let tmps =map  (\(err_yyi: (f32,f32)): bool  ->
                        let (erri, yyi) = err_yyi
                        in erri <= ( tolerance * yyi )
                    ) (zip err yy )
      let breakLoop = and tmps

      -- ...OTHERWiSE, ADJUST STEP FOR NEXT ATTEMPT
      -- scale next step in a default way
      let h = h * scale_fina

      -- limit step to 0.9, because when it gets close to 1, it no longer
      -- makes sense, as 1 is already the next time instance
      -- (added to original algorithm)
      let h = if (h >= 0.9f32) then 0.9f32 else h

      -- if instance+step exceeds range limit, limit to that range
      let h = if ( f32.i32(km1) + h > f32.i32(xmax) ) then f32.i32(xmax - km1)
              else if ( f32.i32(km1) + h + 0.5f32 * h > f32.i32(xmax) )
                   then 0.5f32 * h else h
      in (j+1, h, y_k, breakLoop, scale_fina)
    in ( km1+1, !breakLoop, y_k )
  in (!failed,y_km1)

-----------------------
----- MAiN MODULE -----
-----------------------

let equations (): i32 = 91
let parameters(): i32 = 16

let main (repeat: i32) (eps: f32) (workload: i64)
         (xmax: i32) (y0: [91]f32) (params: [16]f32): (bool, [workload][91]f32) =
  let (oks, y_res) =
    unzip <|
    #[sequential_inner]
    map (\i  ->
            let add_fact = f32.i32(i32.i64 i % repeat)*eps
            let y_row = map (+add_fact) y0
            in solver(xmax, params, y_row))
    (iota workload)

  in ( and (opaque oks), y_res )
