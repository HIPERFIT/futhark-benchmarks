entry main [m'][N'] (_trend: i32) (k': i32) (n': i32) (_freq: f32)
                  (hfrac: f32) (_lam: f32)
                  (_mappingindices : [N']i32)
                  (_images : [m'][N']f32) =
    let k = u64.i32 k'
    let n = u64.i32 n'
    let m = u64.i64 m'
    let N = u64.i64 N'

    let K = 2*k+2

    -- 1. mkX
    let mkX_flops = 8*K*N
    let mkX_words = N*(K+1)

    -- 2. Xsqr = map (matmul_filt Xh Xth) Yh
    let bmmm_flops = 4*m*n*K*K
    let bmmm_words = 2*m*n*K

    -- 3. batch matrix inversion; Xinv = map mat_inv Xsqr
    let bminv_flops = 3*m*K*K
    let bminv_words = 4*m*K*K

    -- 4. 1st mat-mat-mult; beta0 : [m][2k+2]f32 = map (matvecmul_row_filt Xh) Yh
    let mmm1_flops = 3*m*n*K
    let mmm1_words = 2*m*n

    -- 5. 2nd mat-mat-mult; beta : [m][2k+2]f32 = map2 matvecmul_row Xinv beta0
    let mmm2_flops = 2*m*K*K
    let mmm2_words = 2*m*K*K

    -- 6. 3rd mat-mat-mult; y_preds : [m][N]f32 = map (matvecmul_row Xt) beta
    let mmm3_flops = 2*m*N*K
    let mmm3_words = 2*m*N

    -- 7. (Nss, y_errors, val_indss) = map2 (...filterPadWithKeys...)
    let filtpad_flops = 5*m*N
    let filtpad_words = 4*m*N

    -- 8. (hs, nss, sigmas) = map2 (... redomap o redomap ...)
    let sgmred_flops = 8*m*n -- 4*m*n + 4*m
    let sgmred_words = m * (2*n+3)

    -- 9. moving-sums-fst; MO_fsts = ...  I AM HERE!!!
    let hmax = u64.f32 <| (f32.u64 n) * hfrac
    let mvsfst_flops = m * (1 + hmax)
    let mvsfst_words = m * (4 + hmax)

    -- 10. BOUND
    let bound_flops = 4*(N-n)
    let bound_words = 3*(N-n)

    -- 11. final computation and extracting result; (_,_,breaks,means) = ...
    let res_flops = 9*m*(N-n)
    let res_words = m*(7 + 3*(N-n))

    -- Summing up all components
    let total_flops = mkX_flops + bmmm_flops + bminv_flops + mmm1_flops + mmm2_flops + mmm3_flops +
                      filtpad_flops + sgmred_flops + mvsfst_flops + bound_flops + res_flops
    let total_bytes = (mkX_words + bmmm_words + bminv_words + mmm1_words + mmm2_words + mmm3_words +
                       filtpad_words + sgmred_words + mvsfst_words + bound_words + res_words) * 4

    in  (total_flops, total_bytes, m, N, n, K, hmax)
