using Revise

using CSV
using JSON
using DataFrames, DataFramesMeta
using PyPlot
using Setfield

using Oiler

save_results = false

locking_priors = false
locking_frac_prior = 0.5
geol_slip_rate_weight = 2.
tri_distance_weight = 2.


all_vels_file = "../data/vels_consolidated.geojson"
#blocks_file = "../data/cascadia_blocks_big_forearc.geojson"
blocks_file = "../data/cascadia_blocks.geojson"
idl_block_file = "../../../geodesy/global_blocks/old/data/idl_blocks.geojson"
s_us_block_file = "../../../us_faults/s_us_faults/s_us_blocks.geojson"
s_us_fault_file = "../../../us_faults/s_us_faults/s_us_faults.geojson"
cascadia_geol_slip_rates_file = "../data/cascadia_geol_slip_rate_pts.geojson"
new_us_geol_rates_file = "../../../us_faults/s_us_faults/new_us_faults_geol_slip_rates.geojson"
aleut_tris_file = "../../../geodesy/global_block_comps/subduction/sub_tri_meshes/alu_tris.geojson"

#trench_faults_file = "../data/cascadia_trench_faults.geojson"

# tris_file = "../data/cascadia_subduction_tris.geojson"
tris_file = "../data/graham_cascadia_subduction_tris.geojson"
faults_file = "../data/cascadia_faults_big_forearc.geojson"
faults_file = "../data/cascadia_block_faults.geojson"


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

geol_slip_rate_df = vcat(geol_slip_rate_df, casc_geol_slip_rate_df)
geol_slip_rate_df[:, :dextral_err] /= geol_slip_rate_weight
geol_slip_rate_df[:, :extension_err] /= geol_slip_rate_weight

tri_json = JSON.parsefile(tris_file)

aleut_tri_json = JSON.parsefile(aleut_tris_file)


@info "culling blocks"
println("n blocks before ", size(block_df, 1))
bound_df = Oiler.IO.gis_vec_file_to_df("../data/cascadia_qua_cascadia_boundary.geojson")
block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df; epsg=2991)
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
#gnss_vels = [g for g in gnss_vels if sqrt(g.ee^2 + g.en^2) < 2.5 ]
gnss_vels = filter(g->sqrt(g.ee^2+g.en^2) < 2.5, gnss_vels)
println("n gnss vels: ", length(gnss_vels))


@info "doing Explorer and JdF plates"
# Fake JDF vel points
jdf_pt_file = "../data/jdf_vel_pts.csv"
jdf_pts = CSV.read(jdf_pt_file, DataFrame)

jdf_na_pole = Oiler.PoleSphere(
    lon=68.3,
    lat=-32.0,
    rotrate=1.081,
    elat=0.5,
    elon=0.5,
    erotrate=0.15,
    fix="na",
    mov="c006")


jdf_vels = Oiler.BlockRotations.predict_block_vels(jdf_pts[:,:lon], 
    jdf_pts[:,:lat], jdf_na_pole)


jdf_vels = [Oiler.VelocityVectorSphere(vel; vel_type="fault") for vel in jdf_vels]
#jdf_vels = jdf_vels


exp_pac_pole = Oiler.PoleSphere(
    lat=54.80, lon=-116.62, rotrate=3.1, 
    elat=0.25, elon=0.25, erotrate=0.31,
    fix="c024", mov="c112")

exp_pt_file = "../data/explorer_vel_pts.csv"
exp_pts = CSV.read(exp_pt_file, DataFrame)

exp_vels = Oiler.BlockRotations.predict_block_vels(exp_pts[:,:lon],
                                                   exp_pts[:,:lat],
                                                   exp_pac_pole)

exp_vels = [Oiler.VelocityVectorSphere(vel; vel_type="fault") for vel in exp_vels]

@info "doing tris"

tris = Oiler.IO.tris_from_geojson(tri_json)

#@info " getting hanging walls"
#@time tris = map(x->Oiler.Tris.get_tri_hanging_wall(x,block_df, epsg=2991), tris)


function set_tri_rates(tri; ds=25., de=10., ss=0., se=5., cds=0.)
    tri = @set tri.dip_slip_rate = ds
    tri = @set tri.dip_slip_err = de
    tri = @set tri.strike_slip_rate = ss
    tri = @set tri.strike_slip_err = se
    tri = @set tri.cds = cds
    tri
end

function get_tri_priors(tri, pole; ds_only=false, locking_frac=0.9)
    ds, de, ss, se, cds = Oiler.Tris.get_tri_rate_from_pole(tri, pole)
    ds *= locking_frac
    de *= locking_frac
    ss *= locking_frac
    se *= locking_frac

    if ds_only == true
        tri = set_tri_rates(tri; ds=ds, de=de, ss=0.)
    else
        tri = set_tri_rates(tri; ds=ds, de=de, ss=ss, se=se, cds)
    end
    tri
end

if locking_priors
    @info " making priors, locking frac: $locking_frac_prior"
    @time tris = map(x->get_tri_priors(x, jdf_na_pole; locking_frac=locking_frac_prior), tris)
end


#Oiler.Utils.get_path_euler_pole(results["poles"], "c024", "c_21")
#PoleCart
#  x: Float64 2.1524081602284653e-9
#  y: Float64 -8.04901767143253e-9
#  z: Float64 9.010117425991496e-9
#  ex: Float64 1.4605686611857602e-9
#  ey: Float64 8.742507288560317e-10
#  ez: Float64 3.500260387614162e-9
#  cxy: Float64 1.2575269093784705e-18
#  cxz: Float64 -5.0873372724404186e-18
#  cyz: Float64 -3.0234291652265446e-18
#  fix: String "c024"
#  mov: String "c_21"


#aleut_tris = Oiler.IO.tris_from_geojson(aleut_tri_json)
#tris = vcat(tris, aleut_tris)

#tris = []

@info "doing faults"

fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        faults_file, 
                                                        s_us_fault_file;
                                                        block_df=block_df,
                                                        subset_in_bounds=true,
                                                        usd=:upper_seis_depth,
                                                        lsd=:lower_seis_depth)


#faults = filter( x -> x.fw != "c006", faults)
jdf_ridge_vels = filter( x -> x.mov == "c006", fault_vels)
fault_vels = filter( x -> x.mov != "c006", fault_vels)
fault_vels = vcat([jdf_ridge_vels[1]], fault_vels)

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
            geol_slip_rate_vels, 
            jdf_vels, 
            exp_vels
           )

vel_groups = Oiler.group_vels_by_fix_mov(vels)


"""
SOLVE
"""
# poles, tri_rates = 
results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
     tris=tris,
     #tris=[],
     faults=faults,
     tri_distance_weight=tri_distance_weight,
     regularize_tris=true,
     tri_priors=false,
     weighted=true,
     sparse_lhs=true,
     predict_vels=true,
     check_closures=true,
     pred_se=true,
     constraint_method="kkt_sym",
     factorization="lu",
    );

Oiler.ResultsAnalysis.compare_data_results(results=results,
                                           vel_groups=vel_groups,
                                           geol_slip_rate_df=geol_slip_rate_df,
                                           geol_slip_rate_vels=geol_slip_rate_vels,
                                           fault_df=fault_df,
                                           usd=:upper_seis_depth,
                                           lsd=:lower_seis_depth)

tris_ = Oiler.IO.make_tris_from_results(tris, results)
#tris_ = map(x->Oiler.ResultsAnalysis.get_tri_locking_rate(x, 
#                                                         results["poles"],
#                                                         set=true), tris_)

println(results["stats_info"])

if save_results
    Oiler.IO.write_fault_results_to_gj(results, 
    "../results/cascadia_fault_results_reg_20_no_priors.geojson";
    name="NW faults")

    Oiler.IO.write_tri_results_to_gj(tris_, results,
    "../results/cascadia_tri_results_reg_20_no_priors.geojson";
    name="Cascadia")

    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups;
                     name="../results/cascadia_gnss_vels_reg_20_no_priors.csv")

    Oiler.IO.write_solution_poles(
                     "../results/cascadia_na_poles_reg_20_no_priors.csv",
                     results, block_df, "na")

    Oiler.IO.write_results_stats(results, 
                                 "../results/cascadia_reg_20_no_priors.json",
                                 description="standard run, no priors, reg 40")
end



map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, geol_slip_rate_vels,
                                           fault_df, results; usd=:upper_seis_depth,
                                           lsd=:lower_seis_depth)


Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df, 
                                 directory="../web_viewer", ref_pole="na")

println("done")
