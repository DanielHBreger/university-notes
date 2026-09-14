"""Regression tests against independent equations and deliberately bad data."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'chua')]
from diode_analysis import fit_segments, slope_ranges, block_average
from scope_data import read_scope
from rpot import rpot
from sweeplib import list_csvs, load_rpot, sibling, folder_labels
from lorenz_map import maxima, period_samples
from lyapunov import lyapunov, local_fit, mean_ln_slope
from rosenstein import nearest_recurrent, rosenstein, embed
from simulate import make_g, equilibria, integrate, rhs, midpoint, metrics, stability, R0, C1, C2, L
import integration as kern
from identify import (loop_c1, inductor_loop, rayleigh_fit, kennedy_from_pwl, cycle_markers,
                      averaged_cycle, derivative, fit_tank, admittance_rows)


def write_scope(path, t, *channels):
    pd.DataFrame(np.column_stack([t,*channels]),
                 columns=['Time(s)']+[f'CH{i+1}(V)' for i in range(len(channels))]).to_csv(path,index=False)
    return path


def test_divider_recovers_resistance_despite_offset_and_midpoint_gain(tmp_path):
    t=np.arange(4000)*1e-5
    v1=2*np.sin(230*t)+0.5; v2=.4*np.cos(317*t)-.1
    vm=midpoint(v1,v2,673)*1.07+.09
    p=write_scope(tmp_path/'trace.csv',t,v1,v2,vm)
    value,residual=rpot(p,R0,1)
    assert value == pytest.approx(673,rel=1e-10)
    assert residual < 1e-12


def test_divider_gain_ratio(tmp_path):
    t=np.arange(1000)*1e-5
    a,b=np.sin(t*700),np.cos(t*1900)
    p=write_scope(tmp_path/'trace.csv',t,a*1.02,b*.97,midpoint(a,b,510))
    assert rpot(p,R0,.97/1.02)[0] == pytest.approx(510,rel=1e-10)


@pytest.mark.parametrize('kind',['flat','collinear','negative'])
def test_divider_rejects_unidentifiable_and_nonphysical_data(tmp_path,kind):
    t=np.arange(1000)*1e-5; a=np.sin(230*t); b=np.cos(720*t)
    if kind=='flat': a=np.zeros_like(t)
    if kind=='collinear': b=2*a+4
    vm=.5*a+.5*b if kind!='negative' else -.2*a+.8*b
    p=write_scope(tmp_path/'trace.csv',t,a,b,vm)
    with pytest.raises(ValueError): rpot(p,990,1)


@pytest.mark.parametrize('bad',['nan','duplicate','gap','backwards'])
def test_reader_rejects_broken_time_series(tmp_path,bad):
    t=np.arange(1000)*1e-5; x=np.sin(t)
    if bad=='nan': x[10]=np.nan
    if bad=='duplicate': t[10]=t[9]
    if bad=='gap': t[10:]+=1e-3
    if bad=='backwards': t=t[::-1]
    p=write_scope(tmp_path/'trace.csv',t,x,x)
    with pytest.raises(ValueError): read_scope(p)


def test_reader_keeps_named_channel_order_and_float64_time(tmp_path):
    t=1+np.arange(1000)*1e-8
    p=write_scope(tmp_path/'trace.csv',t,t*2,t*3)
    tt,y=read_scope(p,('CH2','CH1'))
    assert np.all(np.diff(tt)>0)
    np.testing.assert_allclose(y[:,0],3*t)


def test_peak_amplitude_period_and_return_time(tmp_path):
    t=np.arange(20000)*1e-5
    x=2*np.sin(2*np.pi*100*t)
    p=write_scope(tmp_path/'trace.csv',t,x,x)
    M,info=maxima(p)
    assert len(M)==20
    np.testing.assert_allclose(M,2,atol=1e-3)
    assert info['f0_hz']==pytest.approx(100,rel=.02)
    np.testing.assert_allclose(np.diff(info['t_peaks']),.01,atol=1e-5)


def test_peaks_do_not_crash_for_short_flat_and_large_window(tmp_path):
    assert period_samples(np.ones(3))==0
    p=write_scope(tmp_path/'trace.csv',np.arange(100)*1e-5,np.ones(100),np.ones(100))
    M,_=maxima(p,period=100000)
    assert not len(M)


def test_piecewise_fit_recovers_known_slopes_without_duplicate_samples():
    x=np.repeat(np.arange(-6,6.01,.1),5)
    breaks=np.array([-4,-2,1,4])
    labels=np.searchsorted(breaks,x)
    slopes=np.array([2,-.5,-1,-.4,2.5]); intercepts=np.array([8,-2,-3,-3.6,-15.2])
    y=slopes[labels]*x+intercepts[labels]
    fits,r2=fit_segments(x,y,min_points=20)
    np.testing.assert_allclose([f['slope_S'] for f in fits],slopes,atol=1e-8)
    assert sum(f['n'] for f in fits)==len(x)
    assert r2==pytest.approx(1,abs=1e-10)


def test_ranges_do_not_hide_invalid_slope_order():
    f=[{'slope_S':v} for v in [1,-2,-1,-.5,1]]
    result=slope_ranges(f,0)
    assert result[0]['valid'] is False
    assert result[1]['total_lo_ohm']==1
    assert result[1]['total_hi_ohm']==2


def test_local_regression_handles_large_coordinate_offset():
    x=1e6+np.linspace(-1,1,1000); y=3*(x-1e6)+2e6
    slopes,pred,se=local_fit(x,y,.2,10)
    np.testing.assert_allclose(slopes,3,atol=1e-7)
    np.testing.assert_allclose(pred,y,atol=1e-7)


def test_logistic_map_exponent_against_ln2():
    x=np.empty(12000); x[0]=.123456789
    for i in range(1,len(x)): x[i]=4*x[i-1]*(1-x[i-1])
    x=x[1000:]; t=np.arange(len(x),dtype=float)
    result=lyapunov(x,t,width=.06,min_width=.01,min_spread=.01,min_r2=.8,min_points=10)
    assert result['status']=='ok'
    assert result['lam']==pytest.approx(np.log(2),abs=.08)
    shuffled=np.random.default_rng(4).permutation(x)
    assert lyapunov(shuffled,t)['status']!='ok'


def test_return_mean_not_median_and_nonmonotonic_times_rejected():
    M=np.ones(40); t=np.r_[0,np.cumsum(np.tile([1,1,4],13))]
    result=lyapunov(M,t)
    assert result['mean_T']==2
    assert result['median_T']==1
    t[10]=t[9]
    assert 'increasing' in lyapunov(M,t)['status']


def test_exact_zero_derivative_not_silently_discarded():
    assert np.isneginf(mean_ln_slope(np.array([0.,2.,2.]),np.zeros(3))[0])


def test_width_uncertainty_changes_window_even_when_floor_dominates():
    x=np.empty(8000);x[0]=.314159
    for i in range(1,len(x)): x[i]=4*x[i-1]*(1-x[i-1])
    result=lyapunov(x,np.arange(len(x),dtype=float),width=.001,min_width=.15,
                   min_spread=.01)
    assert result['status']=='ok'
    assert result['lam_sys']>0.001


def test_nearest_recurrent_matches_exhaustive_search_even_k1():
    X=np.random.default_rng(4).normal(size=(71,3)); last=70
    i,j=nearest_recurrent(X,last,8,1,1,1)
    for a,b in zip(i,j):
        valid=np.flatnonzero(np.abs(np.arange(last)-a)>8)
        assert abs(a-b)>8
        assert np.linalg.norm(X[a]-X[b])==pytest.approx(np.min(np.linalg.norm(X[valid]-X[a],axis=1)))


@pytest.mark.parametrize('args',[{'dt':0},{'period':0},{'fit':(2,1)},{'fit':(0,4)},{'stride':0},{'tau':0}])
def test_rosenstein_invalid_options(args):
    settings=dict(dt=1e-5,period=100); settings.update(args)
    with pytest.raises(ValueError): rosenstein(np.ones(1000),np.ones(1000),**settings)


def test_equilibria_include_exact_origin_and_satisfy_three_equations():
    g=make_g(i0=-.412974865e-3)
    # Use exact fitted offset for a root exactly on the origin.
    from simulate import fitted_offset
    g=make_g(i0=fitted_offset())
    roots=equilibria(g,1650,30)
    assert len(roots)==3
    assert min(abs(np.array(roots))) < 1e-10
    for v1 in roots:
        i=v1/1680; v2=30*i
        np.testing.assert_allclose(rhs(np.array([[v1],[v2],[i]]),np.array([1650]),30,g),0,atol=1e-8)


def test_rk4_matches_exact_linear_circuit_and_time_grid():
    G=.001; g=lambda v:G*v
    t,y=integrate([660],1,.000103,2e-7,20,g,t_skip=.00002)
    Rt=R0+660
    A=np.array([[(-1/Rt-G)/C1,1/(Rt*C1),0],[1/(Rt*C2),-1/(Rt*C2),-1/C2],[0,1/L,-20/L]])
    exact=expm(A*t[-1])@np.array([1,0,0])
    np.testing.assert_allclose(y[-1,:,0],exact,atol=1e-8,rtol=1e-6)
    assert t[0]==pytest.approx(.00002)
    assert t[-1]==pytest.approx(.000103)


def test_outer_segments_extrapolate_linearly_and_offset_is_removed():
    g=make_g()
    assert float(g(30)-g(20))==pytest.approx(10*g.segments[-1][2],rel=1e-9)
    assert float(g(-20)-g(-30))==pytest.approx(10*g.segments[0][2],rel=1e-9)
    assert float(g(0.0))==pytest.approx(0,abs=1e-12)


def test_csv_discovery_and_sidecar_validation(tmp_path):
    t=np.arange(10,dtype=float)
    for name in ['trace10.csv','trace2.csv']: write_scope(tmp_path/name,t,t,t)
    (tmp_path/'summary.csv').write_text('value,result\n1,2\n')
    assert list_csvs(tmp_path,recursive=True)==['trace2.csv','trace10.csv']
    Path(sibling(str(tmp_path),'_rpot.csv')).write_text('filename,rpot_ohm,residual_pct\na,nan,1\nb,500,nan\nc,500,1\n')
    assert load_rpot(tmp_path)=={'c':(500.,1.)}
    assert sibling(str(tmp_path)+'/', '_x')==sibling(str(tmp_path),'_x')
    assert len(set(folder_labels(['/a/forward','/b/forward'])))==2


def test_simulation_skips_sidecars_without_a_sweep_column(tmp_path):
    from simulate import load_measured_bifurcation
    (tmp_path/'sweep').mkdir()
    (tmp_path/'sweep'/'trace1.csv').touch()
    rows=['sweep,filename,rpot_ohm,max_v','sweep,trace1.csv,600,2','sweep,deleted.csv,610,-1',
          'other,trace1.csv,610,-1','sweep,trace1.csv,nan,0','']
    (tmp_path/'sweep_bifurcation_points.csv').write_text('\n'.join(rows))
    (tmp_path/'simulated_bifurcation_points.csv').write_text('rpot_ohm,initial_v1_V,max_v\n600,1,3\n')
    assert load_measured_bifurcation(str(tmp_path))==[(600.,2.,'sweep')]


def test_continuation_carries_the_actual_previous_state(monkeypatch):
    import simulate as s
    monkeypatch.setattr(s, 'sweep_starts', lambda g,r,rl:(np.zeros((3,1)),np.ones((3,1))))
    def hold(y,rt,*args):
        y=y+rt[0]
        return y,[[float(y[0,0])]]
    monkeypatch.setattr(s,'_hold',hold)
    r=np.array([30.,10.,20.])
    got=s.sweep_continuation(None,r,0,1,1,1,'up')
    expected=np.cumsum(R0+np.sort(r))+1
    assert got==[[expected[2]],[expected[0]],[expected[1]]]


def test_compiled_kernel_matches_matrix_exponential():
    from integration import trajectory
    G=.001; Rt=1650.; rl=20.; dt=2e-7
    vk=np.array([-40.,40.]);ik=G*vk
    y=trajectory(np.array([1.,0.,0.]),Rt,rl,C1,C2,L,vk,ik,dt,500,0)
    A=np.array([[(-1/Rt-G)/C1,1/(Rt*C1),0],[1/(Rt*C2),-1/(Rt*C2),-1/C2],[0,1/L,-rl/L]])
    np.testing.assert_allclose(y[-1],expm(A*dt*500)@np.array([1.,0.,0.]),atol=1e-8)


def test_variational_exponent_agrees_with_linear_eigenvalues():
    from benchmark_lyapunov import simulate
    alpha,beta,m=9.,14.286,.2
    _,exponent=simulate(alpha,beta,m,m,5000,.01,transient=10000)
    A=np.array([[-alpha*(1+m),alpha,0],[1,-1,1],[0,-beta,0]])
    expected=np.linalg.eigvals(A).real.max()
    assert exponent==pytest.approx(expected,abs=.025)


def test_pwl_knots_are_continuous_with_the_given_slopes():
    bp = [-6.7, -1.0, 0.9, 6.0]; G = [3.9e-3, -0.42e-3, -0.77e-3, -0.42e-3, 3.9e-3]
    vk, ik = kern.pwl_knots(bp, G)
    assert ik[np.searchsorted(vk, 0.0) - 1] == pytest.approx(G[2] * vk[np.searchsorted(vk, 0.0) - 1])
    np.testing.assert_allclose(np.diff(ik) / np.diff(vk), G, rtol=1e-12)
    assert np.interp(0.0, vk, ik) == 0.0


def test_kennedy_mapping_reproduces_the_segments():
    bp = [-6.7, -1.0, 0.9, 6.0]; G = [3.9e-3, -0.42e-3, -0.77e-3, -0.42e-3, 3.9e-3]
    k = kennedy_from_pwl(G[2], G[1], G[0], bp[1], bp[2], bp[0], bp[3])
    P = kern.bench_params(1e-8, 1e-7, 1e-2, 0, 0, 0, k['RA_ohm'], k['AA'], k['RB_ohm'], k['AB'],
                          k['VpA_V'], k['VnA_V'], k['VpB_V'], k['VnB_V'], 1e6, 1e-7, 1e-6)
    bp2, G2 = kern.static_pwl(P)
    np.testing.assert_allclose(bp2, bp, rtol=1e-9)
    np.testing.assert_allclose(G2, G, rtol=1e-9)


def test_bench_kernel_matches_matrix_exponential_in_the_linear_limit():
    # rails far away, no Rayleigh terms: a five-state linear system
    c1, c2, l0, r0 = 11e-9, 89e-9, 19e-3, 3.0
    RA, AA, RB, AB, sr, tauA, tauB = 260.0, 1.12, 22e3, 7.6, 1e12, 0.05e-6, 1.2e-6
    P = kern.bench_params(c1, c2, l0, 0.0, r0, 0.0, RA, AA, RB, AB, 1e3, -1e3, 1e3, -1e3, sr, tauA, tauB)
    Rt, dt, n = 1800.0, 0.02e-6, 400
    y0 = np.array([0.3, 0.1, 0.0002, 0.2, 1.0])
    y = kern.bench_trajectory(y0, Rt, dt, n, 0, P)
    A = np.array([[-1 / (Rt * c1), 1 / (Rt * c1), 0, 1 / (RA * c1), 1 / (RB * c1)],
                  [1 / (Rt * c2), -1 / (Rt * c2), -1 / c2, 0, 0],
                  [0, 1 / l0, -r0 / l0, 0, 0],
                  [AA / tauA, 0, 0, -1 / tauA, 0],
                  [AB / tauB, 0, 0, 0, -1 / tauB]])
    # the diode current -(1/RA + 1/RB) v1 / c1 sits on the v1 row
    A[0, 0] -= (1 / RA + 1 / RB) / c1
    np.testing.assert_allclose(y[-1], expm(A * dt * n) @ y0, rtol=1e-6, atol=1e-9)


def test_bench_op_amps_saturate_and_hold_at_the_rails():
    P = kern.bench_params(11e-9, 89e-9, 19e-3, 0.8, 2.8, 4700, 260.0, 1.12, 22e3, 7.6, 6.6, -7.5, 6.6, -7.5,
                          0.5e6, 0.05e-6, 1.2e-6)
    y = kern.bench_state(3.0, 0.0, 0.0, P)
    assert y[4] == 6.6 and y[3] == pytest.approx(3.36)
    d = kern.bench_derivative(3.0, 0.0, 0.0, 3.36, 6.6, 1800.0, P)
    assert d[4] == 0.0                       # railed and pushed outward: held
    d = kern.bench_derivative(-3.0, 0.0, 0.0, -3.36, 6.6, 1800.0, P)
    assert d[4] == pytest.approx(-0.5e6)     # slew limited on the way down


def test_loop_integral_c1_and_flux_loop_recover_a_linear_circuit():
    # a sinusoidal cycle on node 1 with a static linear element g = G v1 and C1 = 11 nF
    C1t, G, Rt, T = 11e-9, -0.42e-3, 1800.0, 340e-6
    n = 2048; th = np.arange(n) / n * T; w = 2 * np.pi / T
    v1 = -3.0 + 1.5 * np.sin(w * th)
    dv1 = 1.5 * w * np.cos(w * th)
    v2 = Rt * (C1t * dv1 + G * v1) + v1        # solves C1 v1' = (v2 - v1)/Rt - G v1
    assert loop_c1(v1, v2, dv1, Rt) == pytest.approx(C1t, rel=1e-9)
    # a tank of L, r with C2 known: v2 = L iL' + r iL
    C2t, Lt, rt = 89e-9, 20e-3, 11.0
    # several harmonics, so that C2 and L separate in the admittance fit
    iL = 2e-3 * np.sin(w * th) + 0.5e-3 * np.sin(3 * w * th) + 0.2e-3 * np.sin(5 * w * th)
    diL = w * (2e-3 * np.cos(w * th) + 1.5e-3 * np.cos(3 * w * th) + 1e-3 * np.cos(5 * w * th))
    v2b = Lt * diL + rt * iL
    dv2 = derivative(v2b, T)
    v1b = v2b + Rt * (C2t * dv2 + iL)
    _, _, I, L_eff, r_eff = inductor_loop(v1b, v2b, dv2, Rt, C2t, T)
    assert I == pytest.approx(0.5 * np.ptp(iL), rel=1e-3)
    assert r_eff == pytest.approx(rt, rel=2e-3)
    rows = admittance_rows(v1b, v2b, Rt, T)
    (c2f, lf, rf), resid = fit_tank(rows)
    assert (c2f, lf, rf) == pytest.approx((C2t * 1e9, Lt * 1e3, rt), rel=0.02)


def test_rayleigh_fit_and_cycle_markers():
    I = np.array([1, 5, 10]) * 1e-3
    out = rayleigh_fit(I, 0.019 + 0.4 * I, 2.8 + (8 / (3 * np.pi)) * 4700 * I)
    assert out['L0_H'] == pytest.approx(0.019) and out['nu_H_per_A'] == pytest.approx(0.8)
    assert out['r0_ohm'] == pytest.approx(2.8) and out['rho_ohm_per_A'] == pytest.approx(4700)
    t = np.arange(20000) * 1e-6
    v = np.sin(2 * np.pi * 3000 * t) + 0.02 * np.sin(2 * np.pi * 30000 * t)
    tc = cycle_markers(t, v, 1 / 3000)
    np.testing.assert_allclose(np.diff(tc), 1 / 3000, rtol=1e-3)
    cyc, T, n = averaged_cycle(t, np.column_stack([v, v]), tc, nph=256)
    assert T == pytest.approx(1 / 3000, rel=1e-3) and n >= 50
    np.testing.assert_allclose(cyc[:, 0], np.sin(2 * np.pi * np.arange(256) / 256) + 0.02 * np.sin(2 * np.pi * 10 * np.arange(256) / 256), atol=0.02)


def test_regime_metrics_read_a_synthetic_sweep():
    rs = np.arange(300, 901, 5.0)
    def down(r):
        if r > 850: return [-3.0]
        if r > 800: return [-0.6]
        if r > 780: return [-0.7, -0.5]
        if r > 700: return list(np.linspace(-1.0, 0.0, 30))
        if r > 400: return list(np.linspace(-1.5, 2.5, 40))
        return [6.5]
    def up(r):
        return [6.5] if r < 680 else down(r)
    m = metrics(rs, [down(r) for r in rs], [up(r) for r in rs])
    assert m['first doubling R1 (down)'] == 800
    assert m['chaos onset (down)'] == 780
    assert m['double scroll from (down)'] == 700
    assert m['double scroll ends (down)'] == 405
    assert m['large cycle from (down)'] == 400
    assert m['large cycle survives to (up)'] == 675
    assert m['double scroll from (up)'] == 680


def test_shilnikov_rows_classify_the_double_scroll_geometry():
    from shilnikov import eigen_rows, kind
    g = make_g()
    rows = eigen_rows(g, R0 + 700.0, 0.0, C1, C2, L)
    assert len(rows) == 3
    origin = rows[1]
    assert abs(origin['v1']) < 1e-9 and origin['gamma'] > 0 > origin['sigma']
    assert kind(origin) == 'saddle-focus, 1-D unstable'
    for outer in (rows[0], rows[2]):
        assert outer['gamma'] < 0 < outer['sigma'] and outer['ratio'] < 1
    # the eigenvalues are those of the Jacobian used by stability()
    ev, _ = stability(g, origin['v1'], R0 + 700.0, 0.0)
    assert origin['gamma'] == pytest.approx(ev.real[np.abs(ev.imag) < 1].max())
