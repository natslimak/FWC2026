

# # ========================================
# # ----- JONSWAP Wave Spectrum Function -----
# # ========================================

# def JONSWAP(ws, Hs, Tp, Gamma=None):
#     """
#     JONSWAP (Hasselmann) wave spectrum
    
#     Parameters:
#     -----------
#     ws : float or array
#         Wave frequencies [rad/s]
#     Hs : float
#         Significant wave height [m]
#     Tp : float
#         Peak wave period [s]
#     Gamma : float, optional
#         Peak shape parameter (default: auto from IEC 61400-3)
        
#     Returns:
#     --------
#     S : array
#         Wave spectral density [m^2 s/rad]
#     """
#     # If peak shape parameter gamma is not specified, use recommendation from IEC 61400-3
#     if Gamma is None:
#         TpOvrSqrtHs = Tp / np.sqrt(Hs)
#         if TpOvrSqrtHs <= 3.6:
#             Gamma = 5.0
#         elif TpOvrSqrtHs >= 5.0:
#             Gamma = 1.0
#         else:
#             Gamma = np.exp(5.75 - 1.15 * TpOvrSqrtHs)
    
#     # Handle both scalar and array inputs
#     if isinstance(ws, (list, tuple, np.ndarray)):
#         ws = np.array(ws)
#     else:
#         ws = np.array([ws])
    
#     # Initialize output
#     S = np.zeros(len(ws))
    
#     # Wave frequency in Hz and common terms
#     f = 0.5 / np.pi * ws                        # wave frequencies in Hz
#     fpOvrf4 = pow((Tp * f), -4.0)               # (fp/f)^4 = (Tp*f)^(-4)
#     C = 1.0 - (0.287 * np.log(Gamma))           # normalizing factor
#     Sigma = 0.07 * (f <= 1.0/Tp) + 0.09 * (f > 1.0/Tp)  # peakedness scaling
#     Alpha = np.exp(-0.5 * ((f * Tp - 1.0) / Sigma) ** 2)
    
#     S = 0.5/np.pi * C * 0.3125 * Hs**2 * fpOvrf4 / f * np.exp(-1.25 * fpOvrf4) * Gamma ** Alpha
#     return S


# # ========================================================================
# # -----  CASE 3: Dynamic Frequency-Domain Analysis (Wave Excitation) -----
# # ========================================================================

# print("\n\n=== Dynamic Line Tension Analysis (Frequency Domain) ===")

# # Wave environment parameters (shallow water case)
# Hs = 0.5                                    # Significant wave height [m] (small waves for sheltered site)
# Tp = 3.0                                    # Peak period [s] (typical for wind-driven waves in shallow water)

# # Wave frequency range for analysis [rad/s]
# # For Tp=3s, the peak frequency is omega_p = 2*pi/Tp = 2.094 rad/s
# # We'll analyze well below and above peak frequency to capture full response shape
# w_min = 0.5
# w_max = 4.0
# n_freqs = 60
# omegas = np.linspace(w_min, w_max, n_freqs)
# omega_peak = 2.0 * np.pi / Tp

# # Generate JONSWAP wave spectrum at these frequencies
# Sw = JONSWAP(omegas, Hs=Hs, Tp=Tp)

# print(f"Wave environment: Hs={Hs} m, Tp={Tp} s")
# print(f"Frequency range: {w_min:.3f} to {w_max:.3f} rad/s ({n_freqs} points)")
# print(f"Spectral peak frequency: omega_p={omega_peak:.3f} rad/s")

# # Response Amplitude Operator (RAO) for the floater
# # Format: RAO_fl[freq_index, [surge, sway, heave]]
# # For a simple small floater in shallow water, use moderate RAO values
# # Surge and sway RAO typically range from 0.5 to 2.0 m/m
# # Heave RAO typically ranges from 0.2 to 1.0 m/m
# RAO_surge = 0.8 * np.ones(n_freqs)          # 0.8 m/m surge response
# RAO_sway = 0.0 * np.ones(n_freqs)           # No sway (symmetric system)
# RAO_heave = 0.2 * np.ones(n_freqs)          # 0.2 m/m heave response
# RAO_fl = np.column_stack([RAO_surge, RAO_sway, RAO_heave])  # shape: (n_freqs, 3)

# print(f"Floater RAO: surge={RAO_surge[0]:.2f} m/m, heave={RAO_heave[0]:.2f} m/m")

# # Seabed parameters
# kbot = 3E+06                                # seabed stiffness [N/m]
# cbot = 3E+05                                # seabed damping [N·s/m]

# # ===== Perform dynamic solve on all lines =====
# fig_dyn, axes_dyn = plt.subplots(2, 1, figsize=(12, 8))

# # Store results for all lines
# chain_fairlead_psd = []
# rope_fairlead_psd = []
# chain_node_std = []
# rope_node_std = []

# # Chain lines (first line of each leg)
# for leg_idx in range(len(angles)):
#     line_idx = 2 * leg_idx
#     chain_line = ms.lineList[line_idx]
    
#     print(f"\n--- Chain Line (Leg {leg_idx+1}) ---")
#     print(f"  Length: {chain_line.L:.2f} m, Nodes: {chain_line.nNodes}")
    
#     try:
#         # Run frequency-domain dynamic solver
#         # Returns: T_nodes_amp, T_nodes_psd, T_nodes_std, s, r_static, r_dynamic, r_total, X
#         (T_nodes_amp, T_nodes_psd, T_nodes_std, s, 
#          r_static, r_dynamic, r_total, X) = chain_line.dynamicSolve(
#             omegas, Sw,
#             RAO_A=0,                         # No motion of anchor (fixed point)
#             RAO_B=RAO_fl,                    # Motion of fairlead (floater motion)
#             depth=ms.depth,                  # Water depth [m]
#             kbot=kbot,
#             cbot=cbot,
#             seabed_tol=1e-4,
#             tol=0.01,
#             iters=100,
#             w=0.8
#         )
        
#         # Extract fairlead tension PSD (last node of the line)
#         T_fairlead_psd = T_nodes_psd[:, -1]
#         chain_fairlead_psd.append(T_fairlead_psd)
#         chain_node_std.append(T_nodes_std)
        
#         print(f"  Fairlead tension PSD range: {T_fairlead_psd.min():.2f} to {T_fairlead_psd.max():.2f} N^2 s/rad")
#         print(f"  Tension std dev at fairlead: {T_nodes_std[-1]:.2f} N")
        
#         # Plot
#         axes_dyn[0].semilogy(omegas, T_fairlead_psd, '-o', linewidth=2, markersize=5, 
#                              label=f'Leg {leg_idx+1}')
        
#     except Exception as e:
#         print(f"  ERROR: {e}")
#         chain_fairlead_psd.append(np.zeros(n_freqs))
#         chain_node_std.append(np.zeros(chain_line.nNodes))

# # Format chain plot
# axes_dyn[0].set_xlabel('Wave frequency [rad/s]')
# axes_dyn[0].set_ylabel('Fairlead Tension PSD [N^2 s/rad]')
# axes_dyn[0].set_title(f'Chain Line Dynamic Tension (Wave: Hs={Hs}m, Tp={Tp}s)')
# axes_dyn[0].legend(loc='best')
# axes_dyn[0].grid(True, alpha=0.3, which='both')
# axes_dyn[0].axvline(omega_peak, color='k', linestyle='--', linewidth=1.2, alpha=0.7)
# axes_dyn[0].text(omega_peak + 0.03, axes_dyn[0].get_ylim()[1] / 5, 'omega_p', fontsize=9)

# # Rope lines (second line of each leg)
# for leg_idx in range(len(angles)):
#     line_idx = 2 * leg_idx + 1
#     rope_line = ms.lineList[line_idx]
    
#     print(f"\n--- Rope Line (Leg {leg_idx+1}) ---")
#     print(f"  Length: {rope_line.L:.2f} m, Nodes: {rope_line.nNodes}")
    
#     try:
#         # Run frequency-domain dynamic solver
#         (T_nodes_amp, T_nodes_psd, T_nodes_std, s, 
#          r_static, r_dynamic, r_total, X) = rope_line.dynamicSolve(
#             omegas, Sw,
#             RAO_A=0,
#             RAO_B=RAO_fl,
#             depth=ms.depth,
#             kbot=kbot,
#             cbot=cbot,
#             seabed_tol=1e-4,
#             tol=0.01,
#             iters=100,
#             w=0.8
#         )
        
#         # Extract fairlead tension PSD (last node)
#         T_fairlead_psd = T_nodes_psd[:, -1]
#         rope_fairlead_psd.append(T_fairlead_psd)
#         rope_node_std.append(T_nodes_std)
        
#         print(f"  Fairlead tension PSD range: {T_fairlead_psd.min():.2f} to {T_fairlead_psd.max():.2f} N^2 s/rad")
#         print(f"  Tension std dev at fairlead: {T_nodes_std[-1]:.2f} N")
        
#         # Plot
#         axes_dyn[1].semilogy(omegas, T_fairlead_psd, '-s', linewidth=2, markersize=5, 
#                              label=f'Leg {leg_idx+1}')
        
#     except Exception as e:
#         print(f"  ERROR: {e}")
#         rope_fairlead_psd.append(np.zeros(n_freqs))
#         rope_node_std.append(np.zeros(rope_line.nNodes))

# # Format rope plot
# axes_dyn[1].set_xlabel('Wave frequency [rad/s]')
# axes_dyn[1].set_ylabel('Fairlead Tension PSD [N^2 s/rad]')
# axes_dyn[1].set_title(f'Rope Line Dynamic Tension (Wave: Hs={Hs}m, Tp={Tp}s)')
# axes_dyn[1].legend(loc='best')
# axes_dyn[1].grid(True, alpha=0.3, which='both')
# axes_dyn[1].axvline(omega_peak, color='k', linestyle='--', linewidth=1.2, alpha=0.7)
# axes_dyn[1].text(omega_peak + 0.03, axes_dyn[1].get_ylim()[1] / 5, 'omega_p', fontsize=9)

# plt.tight_layout()
# plt.show()

# # ===== Dynamic design criteria checks =====
# print("\n=== Dynamic Design Criteria Check ===")

# # 1) Peak/characteristic tension margin to MBL
# n_sigma = 3.0
# chain_MBL = float(ms.lineTypes[chainType].get("MBL", np.nan))
# rope_MBL = float(ms.lineTypes[ropeType].get("MBL", np.nan))

# print(f"\n1) Ultimate strength margin using characteristic tension T_char = T_static + {n_sigma:.1f}*sigma")
# print("LineType | Leg | T_static_end(N) | sigma_end(N) | T_char(N) | MBL(N) | Utilization(%) | Margin(%)")

# for i in range(len(angles)):
#     ch_line = ms.lineList[2 * i]
#     rp_line = ms.lineList[2 * i + 1]

#     # Use larger of end tensions as conservative static end tension
#     ch_static_end = float(max(ch_line.TA, ch_line.TB))
#     rp_static_end = float(max(rp_line.TA, rp_line.TB))

#     # Use fairlead-end sigma (last node)
#     ch_sigma_end = float(chain_node_std[i][-1]) if len(chain_node_std[i]) > 0 else 0.0
#     rp_sigma_end = float(rope_node_std[i][-1]) if len(rope_node_std[i]) > 0 else 0.0

#     ch_T_char = ch_static_end + n_sigma * ch_sigma_end
#     rp_T_char = rp_static_end + n_sigma * rp_sigma_end

#     ch_util = 100.0 * ch_T_char / chain_MBL if chain_MBL > 0 else np.nan
#     rp_util = 100.0 * rp_T_char / rope_MBL if rope_MBL > 0 else np.nan

#     ch_margin = 100.0 - ch_util if np.isfinite(ch_util) else np.nan
#     rp_margin = 100.0 - rp_util if np.isfinite(rp_util) else np.nan

#     print(f"Chain    | {i+1:>3} | {ch_static_end:14.2f} | {ch_sigma_end:11.2f} | {ch_T_char:9.2f} | {chain_MBL:7.2f} | {ch_util:13.3f} | {ch_margin:8.3f}")
#     print(f"Rope     | {i+1:>3} | {rp_static_end:14.2f} | {rp_sigma_end:11.2f} | {rp_T_char:9.2f} | {rope_MBL:7.2f} | {rp_util:13.3f} | {rp_margin:8.3f}")

# # 2) Anchor uplift margin under dynamic amplification
# print("\n2) Anchor uplift with dynamic amplification (chain anchor end)")
# print("Leg | Static_TA(N) | sigma_anchor(N) | T_dyn_char(N) | Angle(deg) | T_vertical_dyn(N) | Anchor_W(N) | Margin(N) | Status")

# for i in range(len(angles)):
#     ch_line = ms.lineList[2 * i]

#     # Geometric angle at anchor from current equilibrium
#     r_anchor = ch_line.rA
#     r_connector = ch_line.rB
#     dz = r_connector[2] - r_anchor[2]
#     dx_dy = np.sqrt((r_connector[0] - r_anchor[0])**2 + (r_connector[1] - r_anchor[1])**2)
#     theta_rad = np.arctan2(dz, dx_dy)
#     theta_deg = np.degrees(theta_rad)

#     T_static_anchor = float(ch_line.TA)
#     sigma_anchor = float(chain_node_std[i][0]) if len(chain_node_std[i]) > 0 else 0.0
#     T_dyn_char_anchor = T_static_anchor + n_sigma * sigma_anchor

#     T_vertical_dyn = T_dyn_char_anchor * np.sin(theta_rad)
#     margin_N = anchor_weight_N - T_vertical_dyn
#     status = "OK" if margin_N > 0 else "FAIL"

#     print(f"{i+1:>3} | {T_static_anchor:12.2f} | {sigma_anchor:15.2f} | {T_dyn_char_anchor:11.2f} | {theta_deg:10.2f} | {T_vertical_dyn:16.2f} | {anchor_weight_N:10.2f} | {margin_N:8.2f} | {status}")

# # 3) Fatigue relevance based on PSD concentration near frequent sea-state band
# print("\n3) Fatigue relevance: PSD concentration near peak sea-state band")
# band_low = 0.8 * omega_peak
# band_high = 1.2 * omega_peak
# band_mask = (omegas >= band_low) & (omegas <= band_high)

# print(f"Peak band used for fatigue relevance: [{band_low:.3f}, {band_high:.3f}] rad/s")
# print("LineType | Leg | PSD_peak_freq(rad/s) | PeakBandEnergyFrac(%) | Indicator")

# for i in range(len(angles)):
#     ch_psd = np.asarray(chain_fairlead_psd[i], dtype=float)
#     rp_psd = np.asarray(rope_fairlead_psd[i], dtype=float)

#     # Chain
#     ch_peak_idx = int(np.argmax(ch_psd)) if ch_psd.size > 0 else 0
#     ch_peak_w = float(omegas[ch_peak_idx])
#     ch_total = float(simpson(ch_psd, omegas)) if ch_psd.size > 1 else 0.0
#     ch_band = float(simpson(ch_psd[band_mask], omegas[band_mask])) if np.count_nonzero(band_mask) > 1 else 0.0
#     ch_frac = 100.0 * ch_band / ch_total if ch_total > 0 else 0.0
#     ch_indicator = "HIGH" if ch_frac >= 40.0 else "MODERATE" if ch_frac >= 20.0 else "LOW"

#     # Rope
#     rp_peak_idx = int(np.argmax(rp_psd)) if rp_psd.size > 0 else 0
#     rp_peak_w = float(omegas[rp_peak_idx])
#     rp_total = float(simpson(rp_psd, omegas)) if rp_psd.size > 1 else 0.0
#     rp_band = float(simpson(rp_psd[band_mask], omegas[band_mask])) if np.count_nonzero(band_mask) > 1 else 0.0
#     rp_frac = 100.0 * rp_band / rp_total if rp_total > 0 else 0.0
#     rp_indicator = "HIGH" if rp_frac >= 40.0 else "MODERATE" if rp_frac >= 20.0 else "LOW"

#     print(f"Chain    | {i+1:>3} | {ch_peak_w:19.3f} | {ch_frac:21.2f} | {ch_indicator}")
#     print(f"Rope     | {i+1:>3} | {rp_peak_w:19.3f} | {rp_frac:21.2f} | {rp_indicator}")

# print("\nInterpretation: high utilization or negative margins require redesign (line type, length, anchor mass, or load assumptions).")
# print("\nDynamic analysis complete.")

