using Revise

using CSV
using JSON
using DataFrames, DataFramesMeta
using PyPlot
using Setfield

using Oiler


save_results = false

midas_vel_file = "../geod/midas_no_ak.geojson"
gsrm_vel_file = "../geod/gsrm_na_rel_no_ak.geojson"
all_vels_file = "../data/vels_consolidated.geojson"
elliott_vels_file = "../data/elliott_freymueller_2020_vels.geojson"
blocks_file = "../data/cascadia_blocks.geojson"
idl_block_file = "../../../geodesy/global_blocks/old/data/idl_blocks.geojson"
s_us_block_file = "../../../us_faults/s_us_faults/s_us_blocks.geojson"
s_us_fault_file = "../../../us_faults/s_us_faults/s_us_faults.geojson"
cascadia_geol_slip_rates_file = "../data/cascadia_geol_slip_rate_pts.geojson"
cali_geol_slip_rates_file = "../../../us_faults/s_us_faults/ca_geol_slip_rates.geojson"
new_us_geol_rates_file = "../../../us_faults/s_us_faults/new_us_faults_geol_slip_rates.geojson"
#aleut_tris_file = "../../../geodesy/global_block_comps/subduction/sub_tri_meshes/alu_tris.geojson"
aleut_tris_file = "../../../geodesy/global_block_comps/subduction/sub_tri_meshes/alu_tris_slab2.geojson"

#trench_faults_file = "../data/cascadia_trench_faults.geojson"

# tris_file = "../data/cascadia_subduction_tris.geojson"
tris_file = "../data/graham_cascadia_subduction_tris.geojson"
faults_file = "../data/cascadia_block_faults.geojson"

midas_df = Oiler.IO.gis_vec_file_to_df(midas_vel_file)
gsrm_df = Oiler.IO.gis_vec_file_to_df(gsrm_vel_file)
elliott_df = Oiler.IO.gis_vec_file_to_df(elliott_vels_file)
gnss_df = Oiler.IO.gis_vec_file_to_df(all_vels_file)
fault_df = Oiler.IO.gis_vec_file_to_df(faults_file)
s_us_fault_df = Oiler.IO.gis_vec_file_to_df(s_us_fault_file)

#trench_fault_df = Oiler.IO.gis_vec_file_to_df(trench_faults_file)

#fault_df = vcat(fault_df, s_us_fault_df, trench_fault_df)
fault_df = vcat(fault_df, s_us_fault_df)

block_df = Oiler.IO.gis_vec_file_to_df(blocks_file)
idl_block_df = Oiler.IO.gis_vec_file_to_df(idl_block_file)
s_us_block_df = Oiler.IO.gis_vec_file_to_df(s_us_block_file)
block_df = vcat(block_df, idl_block_df, s_us_block_df)
#block_df = vcat(block_df, idl_block_df)

casc_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cascadia_geol_slip_rates_file)
geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(new_us_geol_rates_file)
cali_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cali_geol_slip_rates_file)

geol_slip_rate_df = vcat(geol_slip_rate_df, casc_geol_slip_rate_df, cali_geol_slip_rate_df)

tri_json = JSON.parsefile(tris_file)

aleut_tri_json = JSON.parsefile(aleut_tris_file)


@info "culling blocks"
#println("n blocks before ", size(block_df, 1))
#bound_df = Oiler.IO.gis_vec_file_to_df("../data/cascadia_qua_cascadia_boundary.geojson")
#block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df; epsg=2991)
println("n blocks after ", size(block_df, 1))

# load GNSS data

@info "doing GNSS"
@time gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df, block_df;
                                                           fix="na",
                                                           ve=:e_vel,
                                                           vn=:n_vel,
                                                           ee=:e_err,
                                                           en=:n_err,
                                                           name=:station,
                                                           epsg=2991)

println("n gnss vels: ", length(gnss_vels))

# Fake JDF vel points
jdf_pt_file = "../data/jdf_vel_pts.csv"
jdf_pts = CSV.read(jdf_pt_file, DataFrame)

# From DeMets et al 2010, pre-converted to Cartesian
jdf_na_pole = Oiler.PoleCart(
    x=5.915996643479488e-9,
    y=1.486624308775784e-8,
    z=-9.997991640995085e-9,
    ex=sqrt(103.05e-8),
    ey=sqrt(186.44e-8),
    ez=sqrt(259.51e-8),
    cxy=135.05e-8,
    cxz=-155.74e-8,
    cyz=-203.56e-8,
    fix="na",
    mov="c006"
)

jdf_vels = Oiler.BlockRotations.predict_block_vels(jdf_pts[:,:lon], 
    jdf_pts[:,:lat], jdf_na_pole)


jdf_vels = [Oiler.VelocityVectorSphere(vel; vel_type="fault") for vel in jdf_vels]
#jdf_vels = jdf_vels


exp_pac_pole = Oiler.PoleSphere(
    lat=53.99, lon=120.04, rotrate=3.3, 
    elat=0.25, elon=0.25, erotrate=0.33,
    fix="c024", mov="c112")

exp_pac_pole2 = Oiler.PoleSphere(
    lat=54.80, lon=-116.62, rotrate=3.1, 
    elat=0.25 * 0.5, elon=0.25 * 0.5, erotrate=0.31 * 0.1,
    fix="c024", mov="c112")

exp_pt_file = "../data/explorer_vel_pts.csv"
exp_pts = CSV.read(exp_pt_file, DataFrame)

exp_vels = Oiler.BlockRotations.predict_block_vels(exp_pts[:,:lon],
                                                   exp_pts[:,:lat],
                                                   exp_pac_pole2)

exp_vels = [Oiler.VelocityVectorSphere(vel; vel_type="fault") for vel in exp_vels]

cascadia_tris = Oiler.IO.tris_from_geojson(tri_json)


cascadia_tris = Oiler.Utils.tri_priors_from_pole(cascadia_tris, jdf_na_pole,
                                                 locking_fraction=0.5,
                                                 depth_adjust=true,
                                                 err_coeff=1.)
aleut_tris = Oiler.IO.tris_from_geojson(aleut_tri_json)

# Pac-NA pole for Aleut priors
pac_na_pole = Oiler.PoleCart(
  x = -2.4885570284209037e-9,
  y = 8.602351086527437e-9,
  z = -1.0424548581912783e-8,
  ex = 5.4017557235402144e-11,
  ey = 3.829966190320154e-11,
  ez = 4.100348456652813e-11,
  cxy = 1.0584250692981392e-21,
  cxz = -9.478742893924508e-22,
  cyz = -7.463051012669244e-22,
  fix = "na",
  mov = "c024",
)

aleut_tris = Oiler.Utils.tri_priors_from_pole(aleut_tris, pac_na_pole,
                                              locking_fraction=0.25,
                                              depth_adjust=true,
                                              depth_max=100.,
                                              err_coeff=1e6)

tris = vcat(cascadia_tris, aleut_tris)

#function set_tri_rates(tri; ds=20., de=20., ss=0., se=10.)
#    tri = @set tri.dip_slip_rate = ds
#    tri = @set tri.dip_slip_err = de
#    tri = @set tri.strike_slip_rate = ss
#    tri = @set tri.strike_slip_err = se
#    tri
#end
#
#tris = map(set_tri_rates, tris)
#tris = []

# faults

fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        faults_file, 
                                                        s_us_fault_file;
                                                        block_df=block_df,
                                                        subset_in_bounds=true,
                                                        usd=:upper_seis_depth,
                                                        lsd=:lower_seis_depth)
# filter ridge vels which immobilize JdF plate
jdf_ridge_vels = filter( x -> x.mov == "c006", fault_vels)
fault_vels = filter( x -> x.mov != "c006", fault_vels)

println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))


@info "doing geol slip rates"
geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                        geol_slip_rate_df, 
                                                        fault_df;
                                                        usd="upper_seis_depth",
                                                        lsd="lower_seis_depth")
println("n fault slip rate vels: ", length(geol_slip_rate_vels))


vels = vcat(fault_vels, 
            gnss_vels, 
            jdf_vels, 
            geol_slip_rate_vels, 
            exp_vels
           )

vel_groups = Oiler.group_vels_by_fix_mov(vels)


"""
SOLVE
"""
# poles, tri_rates = 
results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
     tris=tris,
     regularize_tris=true,
     #tris=[],
     faults=faults,
     tri_distance_weight=20.,
     tri_priors=false,
     weighted=true,
     sparse_lhs=true,
     predict_vels=true,
     check_closures=true,
     pred_se=false,
     check_nans=true,
     constraint_method="kkt_sym",
     factorization="lu",
    );


if save_results
    Oiler.IO.write_fault_results_to_gj(results, 
    "../results/w_us_cascadia_fault_results.geojson";
    name="Western North America faults")

    Oiler.IO.write_tri_results_to_gj(tris, results,
    "../results/cascadia_aleut_tri_results.geojson";
    name="Cascadia/Aleut tri rates")
end

Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)
Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, geol_slip_rate_vels,
                               fault_df, results, usd=:upper_seis_depth,
                               lsd=:lower_seis_depth)

#na_rel_poles = [Oiler.Utils.get_path_euler_pole(pole_arr, "na",
#                                                string(block_df[i, :fid]))
#                for i in 1:size(block_df, 1)]
#if save_results
#    CSV.write("../results/na_rel_poles.csv", 
#            Oiler.IO.poles_to_df(na_rel_poles, convert_to_sphere=true))
#end

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 directory="../web_viewer", ref_pole="na")
println("done")
