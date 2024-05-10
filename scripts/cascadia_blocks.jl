using Revise

using CSV
using JSON
using DataFrames, DataFramesMeta
using PyPlot
using Setfield

using Oiler

geol_slip_rate_weight = 2.
save_results = false


# cascadia blocks
cascadia_blocks_file = "../data/cascadia_blocks.geojson"
cascadia_faults_file = "../data/cascadia_block_faults.geojson"
cascadia_geol_slip_rates_file = "../data/cascadia_geol_slip_rate_pts.geojson"


# s us blocks
s_us_block_file = "../../s_us_faults/s_us_blocks.geojson"
s_us_fault_file = "../../s_us_faults/s_us_faults.geojson"
cali_geol_slip_rates_file = "../../s_us_faults/ca_geol_slip_rates.geojson"
new_us_geol_rates_file = "../../s_us_faults/new_us_faults_geol_slip_rates.geojson"

# global blocks
global_block_file = "../../global_scale_plates/global_scale_plates.geojson"
global_fault_file = "../../global_scale_plates/global_scale_faults.geojson"
global_slip_rate_file = "../../global_scale_plates/global_scale_slip_rates.geojson"

# nea blocks
nea_block_file = "../../ne_asia_blocks/ne_asia_blocks.geojson"
nea_fault_file = "../../ne_asia_blocks/ne_asia_faults.geojson"

# geodesy
midas_vel_file = "../geod/midas_no_ak.geojson"
gsrm_vel_file = "../geod/gsrm_na_rel_no_ak.geojson"
all_vels_file = "../data/vels_consolidated.geojson"
elliott_vels_file = "../data/elliott_freymueller_2020_vels.geojson"


# subduction
aleut_tris_file = "../../subduction/sub_tri_meshes/alu_tris_slab2_updated.geojson"
cascadia_tris_file = "../data/jdf_explorer_interface.geojson"
#cascadia_tris_file = "/home/itchy/Desktop/cascadia_tris_graham_priors.geojson"
#tris_file = "../data/graham_cascadia_subduction_tris.geojson"


@info "joining blocks"
cascadia_block_df = Oiler.IO.gis_vec_file_to_df(cascadia_blocks_file)
s_us_block_df = Oiler.IO.gis_vec_file_to_df(s_us_block_file)
nea_block_df = Oiler.IO.gis_vec_file_to_df(nea_block_file)

global_block_df = Oiler.IO.gis_vec_file_to_df(global_block_file,
                                              fid_drop=["ant"])

block_df = vcat(cascadia_block_df, 
                s_us_block_df,
                nea_block_df,
                global_block_df, 
                cols=:union)


@info "culling blocks"
println("n blocks before ", size(block_df, 1))
bound_df = Oiler.IO.gis_vec_file_to_df("../data/nw_nam_boundary.geojson")
#bound_df = Oiler.IO.gis_vec_file_to_df("../data/cascadia_qua_cascadia_boundary.geojson")
#bound_df = Oiler.IO.gis_vec_file_to_df("../data/just_casc.geojson")
#bound_df = Oiler.IO.gis_vec_file_to_df("../data/ak_bounds.geojson")
block_df = Oiler.IO.get_blocks_in_bounds!(block_df, 
                                          bound_df)
println("n blocks after ", size(block_df, 1))

# load GNSS data

@info "doing GNSS"
#midas_df = Oiler.IO.gis_vec_file_to_df(midas_vel_file)
#gsrm_df = Oiler.IO.gis_vec_file_to_df(gsrm_vel_file)
#elliott_df = Oiler.IO.gis_vec_file_to_df(elliott_vels_file)
gnss_df = Oiler.IO.gis_vec_file_to_df(all_vels_file)

@time gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df, block_df;
                                                           fix="na",
                                                           ve=:e_vel,
                                                           vn=:n_vel,
                                                           ee=:e_err,
                                                           en=:n_err,
                                                           name=:station,
                                                           #epsg=2991,
                                                          )

println("n gnss vels: ", length(gnss_vels))

@info "doing JdF vels"
# Fake JDF vel points
jdf_pt_file = "../data/jdf_vel_pts.csv"
jdf_pts = CSV.read(jdf_pt_file, DataFrame)

# From DeMets et al 2010, pre-converted to Cartesian
jdf_na_pole = Oiler.PoleSphere(lon=66.8, lat=-32.5, rotrate=1.099,
                                elon=2., elat=2., erotrate=0.1,
                               fix="na", mov="c006")


jdf_vels = Oiler.BlockRotations.predict_block_vels(jdf_pts[:,:lon], 
    jdf_pts[:,:lat], jdf_na_pole)


jdf_vels = [Oiler.VelocityVectorSphere(vel; vel_type="plate") for vel in jdf_vels]
#jdf_vels = jdf_vels

# eplorer poles from Braunmiller and Nabelek, 2002
exp_pac_pole = Oiler.PoleSphere(
    lat=53.99, lon=-120.04, rotrate=3.3, 
    elat=0.25, elon=0.25, erotrate=0.33,
    fix="c024", mov="c112")

exp_pac_pole2 = Oiler.PoleSphere(
    lat=54.80, lon=-116.62, rotrate=3.1, 
    elat=0.25 * 0.5, elon=0.25 * 0.5, erotrate=0.31 * 0.1,
    fix="c024", mov="c112")


exp_nam_pole = Oiler.PoleSphere(
    lat=52.67, lon=-131.90, rotrate=2.65,
    elat=0.25, elon=0.35, erotrate=0.265,
    fix="na", mov="c112")


exp_pt_file = "../data/explorer_vel_pts.csv"
exp_pts = CSV.read(exp_pt_file, DataFrame)

exp_vels = Oiler.BlockRotations.predict_block_vels(exp_pts[:,:lon],
                                                   exp_pts[:,:lat],
                                                   exp_nam_pole)

exp_vels = [Oiler.VelocityVectorSphere(vel; vel_type="plate") for vel in exp_vels]



@info "doing tris"

cascadia_tri_json = JSON.parsefile(cascadia_tris_file)
aleut_tri_json = JSON.parsefile(aleut_tris_file)

cascadia_tris = Oiler.IO.tris_from_geojson(cascadia_tri_json)
jdf_tris = filter(t -> t.fw == "jdf", cascadia_tris)
jdf_tris = Oiler.Utils.tri_priors_from_pole(jdf_tris, jdf_na_pole,
                                                 locking_fraction=0.5,
                                                 depth_adjust=true,
                                                 err_coeff=4.0)

exp_tris = filter(t -> t.fw == "explorer", cascadia_tris)
exp_tris = Oiler.Utils.tri_priors_from_pole(exp_tris, exp_nam_pole,
                                                 locking_fraction=0.5,
                                                 depth_adjust=true,
                                                 err_coeff=8.0)

cascadia_tris = vcat(jdf_tris, exp_tris)

function set_tri_rates(tri)
   #tri = @set tri.dip_slip_rate = 20.
   tri = @set tri.dip_slip_err = 2.e-2
   #tri = @set tri.strike_slip_rate = 0.
   tri = @set tri.strike_slip_err = 2.e-2
   tri
end

#cascadia_tris = map(set_tri_rates, cascadia_tris)

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
                                              locking_fraction=0.7,
                                              depth_adjust=true,
                                              depth_max=50.,
                                              err_coeff=0.5)

tris = cascadia_tris
tris = aleut_tris
tris = vcat(cascadia_tris, aleut_tris)
#tris = cascadia_tris


# faults
@info "doing faults"
fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        cascadia_faults_file, 
                                                        s_us_fault_file;
                                                        block_df=block_df,
                                                        #fid_drop=["cf197"],
                                                        subset_in_bounds=true,
                                                        check_blocks=true,
                                                        adjust_err_by_dip=true,
                                                        usd=:upper_seis_depth,
                                                        lsd=:lower_seis_depth,
                                                       )
# filter ridge vels which immobilize JdF plate
#jdf_ridge_vels = filter( x -> x.mov == "c006", fault_vels)
fault_vels = filter( x -> x.mov != "c006", fault_vels)
fault_vels = filter( x -> x.mov != "c112", fault_vels)

println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))

@info "doing non-fault block boundaries"
@time non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(b->Oiler.Boundaries.boundary_to_vels(b, ee=0.5, en=0.5), 
                      non_fault_bounds)...)

bound_vels = filter(x->x.fix != "c006", bound_vels)
bound_vels = filter(x->x.mov != "c006", bound_vels)
bound_vels = filter(x->x.fix != "c112", bound_vels)
bound_vels = filter(x->x.mov != "c112", bound_vels)
bound_vels = filter(x->x.fix != "c024", bound_vels)
bound_vels = filter(x->x.mov != "c024", bound_vels)

println("n non-fault-bound vels: ", length(bound_vels))

@info "doing geol slip rates"
casc_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cascadia_geol_slip_rates_file)
imw_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(new_us_geol_rates_file)
cali_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cali_geol_slip_rates_file)
glob_geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(global_slip_rate_file)

geol_slip_rate_df = vcat(imw_geol_slip_rate_df, 
                         casc_geol_slip_rate_df, 
                         cali_geol_slip_rate_df,
                         glob_geol_slip_rate_df,
                        )

geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                                        geol_slip_rate_df, 
                                                        fault_df;
                                                        weight = geol_slip_rate_weight,
                                                        usd="upper_seis_depth",
                                                        lsd="lower_seis_depth")
println("n fault slip rate vels: ", length(geol_slip_rate_vels))


vels = vcat(fault_vels, 
            gnss_vels, 
            #jdf_vels, 
            geol_slip_rate_vels, 
            exp_vels,
            bound_vels,
           )

vel_groups = Oiler.group_vels_by_fix_mov(vels)


"""
SOLVE
"""
# poles, tri_rates = 
results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
     tris=tris,
     regularize_tris=true,
     faults=faults,
     tri_distance_weight=15.,
     tri_priors=true,
     weighted=true,
     sparse_lhs=true,
     predict_vels=true,
     check_closures=false,
     pred_se=false,
     check_nans=true,
     constraint_method="kkt_sym",
     factorization="lu",
    );

#Oiler.ResultsAnalysis.get_block_centroid_vels(results, block_df; fix="na")
#Oiler.ResultsAnalysis.compare_data_results(results=results,
#                                           vel_groups=vel_groups,
#                                           geol_slip_rate_df=geol_slip_rate_df,
#                                           geol_slip_rate_vels=geol_slip_rate_vels,
#                                           fault_df=fault_df,
#                                           usd=:upper_seis_depth,
#                                           lsd=:lower_seis_depth,
#                                          )
#Oiler.ResultsAnalysis.calculate_resid_block_strain_rates(results)
#
#println(results["stats_info"])

if save_results
    Oiler.IO.write_fault_results_to_gj(results, 
    "../results/cascadia_fault_results.geojson";
    name="Western North America faults")

    Oiler.IO.write_tri_results_to_gj(tris, results,
    "../results/cascadia_aleut_tri_results.geojson";
    name="Cascadia/Aleut tri rates")

    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups,
                        name="../results/cascadia_gnss_results.csv")
end

#Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)
#Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, geol_slip_rate_vels,
#                               fault_df, results, usd=:upper_seis_depth,
#                               lsd=:lower_seis_depth)
#Oiler.Plots.plot_tri_prior_post(tris, results)

#na_rel_poles = [Oiler.Utils.get_path_euler_pole(pole_arr, "na",
#                                                string(block_df[i, :fid]))
#                for i in 1:size(block_df, 1)]
if save_results
    #CSV.write("../results/na_rel_poles.csv", 
    #        Oiler.IO.poles_to_df(na_rel_poles, convert_to_sphere=true))
end

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 directory="../web_viewer", ref_pole="na")
println("done")
