using Revise

using CSV
using JSON
using DataFrames, DataFramesMeta
using PyPlot
using Setfield

using Oiler


save_results = true

midas_vel_file = "../geod/midas_no_ak.geojson"
gsrm_vel_file = "../geod/gsrm_na_rel_no_ak.geojson"
all_vels_file = "../data/vels_consolidated.geojson"
elliott_vels_file = "../data/elliott_freymueller_2020_vels.geojson"
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

geol_slip_rate_df = vcat(geol_slip_rate_df, casc_geol_slip_rate_df)

tri_json = JSON.parsefile(tris_file)

aleut_tri_json = JSON.parsefile(aleut_tris_file)


@info "culling blocks"
#println("n blocks before ", size(block_df, 1))
#bound_df = Oiler.IO.gis_vec_file_to_df("../data/cascadia_qua_cascadia_boundary.geojson")
#block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df; epsg=2991)
println("n blocks after ", size(block_df, 1))

# load GNSS data

#@info "doing GSRM GNSS"
#@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_df, block_df;
##@info "doing GNSS"
##@time gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df, block_df;
#                                                           fix="na",
#                                                           ve=:e_vel,
#                                                           vn=:n_vel,
#                                                           ee=:e_err,
#                                                           en=:n_err,
#                                                           name=:station,
#                                                           epsg=2991)
#
#@info "doing MIDAS"
#@time midas_vels = Oiler.IO.make_vels_from_gnss_and_blocks(midas_df, block_df;
##@time gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df, block_df;
#                                                           fix="na",
#                                                           ve=:e_vel,
#                                                           vn=:n_vel,
#                                                           ee=:e_err,
#                                                           en=:n_err,
#                                                           name=:station,
#                                                           epsg=2991)
#
#@info "doing Elliott"
#@time elliott_vels = Oiler.IO.make_vels_from_gnss_and_blocks(elliott_df, block_df;
##@time gnss_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_df, block_df;
#                                                           fix="na",
#                                                           ve=:e_vel,
#                                                           vn=:n_vel,
#                                                           ee=:e_err,
#                                                           en=:n_err,
#                                                           name=:station,
#                                                           epsg=2991)
#gnss_vels = vcat(midas_vels, gsrm_vels, elliott_vels)

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
    elat=0.25, elon=0.25, erotrate=0.31,
    fix="c024", mov="c112")

exp_pt_file = "../data/explorer_vel_pts.csv"
exp_pts = CSV.read(exp_pt_file, DataFrame)

exp_vels = Oiler.BlockRotations.predict_block_vels(exp_pts[:,:lon],
                                                   exp_pts[:,:lat],
                                                   exp_pac_pole2)

exp_vels = [Oiler.VelocityVectorSphere(vel; vel_type="fault") for vel in exp_vels]

tris = Oiler.IO.tris_from_geojson(tri_json)

aleut_tris = Oiler.IO.tris_from_geojson(aleut_tri_json)

tris = vcat(tris, aleut_tris)

function set_tri_rates(tri; ds=25., de=5., ss=0., se=5.)
    tri = @set tri.dip_slip_rate = ds
    tri = @set tri.dip_slip_err = de
    tri = @set tri.strike_slip_rate = ss
    tri = @set tri.strike_slip_err = se
    tri
end

tris = map(set_tri_rates, tris)
#tris = []

# faults

fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        faults_file, 
                                                        s_us_fault_file;
                                                        block_df=block_df,
                                                        subset_in_bounds=true,
                                                        usd=:upper_seis_depth,
                                                        lsd=:lower_seis_depth)

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
     #tris=[],
     faults=faults,
     tri_distance_weight=10.,
     tri_priors=true,
     weighted=true,
     sparse_lhs=true,
     predict_vels=true,
     check_closures=true,
     pred_se=true,
     constraint_method="kkt",
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

pole_arr = collect(values(results["poles"]))
pole_arr = [pole for pole in pole_arr if typeof(pole) == Oiler.PoleCart]

obs_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(vel_groups)]
pred_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(results["predicted_vels"])]

glons, glats = Oiler.Utils.get_coords_from_vel_array(obs_vels)
ove = [v.ve for v in obs_vels]
ovn = [v.vn for v in obs_vels]

pve = [v.ve for v in pred_vels]
pvn = [v.vn for v in pred_vels]

figure(figsize=(14, 14))

cm = get_cmap(:viridis)

function get_tri_total_rate(tri)
    ds = results["tri_slip_rates"][tri.name]["dip_slip"]
    ss = results["tri_slip_rates"][tri.name]["strike_slip"]
    total_rate = sqrt(ds^2 + ss^2)
end

if length(tris) > 0
    tri_rates = [get_tri_total_rate(tri) for tri in tris]
    tri_rate_min = minimum(tri_rates)
    tri_rate_max = maximum(tri_rates)
end

function plot_tri(tri; vmin=tri_rate_min, vmax=tri_rate_max)
    lons = [tri.p1[1], tri.p2[1], tri.p3[1], tri.p1[1]]
    lats = [tri.p1[2], tri.p2[2], tri.p3[2], tri.p1[2]]
    total_rate = get_tri_total_rate(tri)
    rate_frac = (total_rate - vmin) / (vmax - vmin)
    color = cm(rate_frac)
    # color = [0., 0., rate_frac]
    fill(lons, lats, color=color, alpha=0.25, zorder=0)
end

for tri in tris
    plot_tri(tri)
end


for fault in faults
    plot(fault.trace[:,1], fault.trace[:,2], "k-", lw=0.3)
end

quiver(glons, glats, ove, ovn, color="b", scale=300)
quiver(glons, glats, pve, pvn, color="r", scale=300)
axis("equal")

pred_geol_slip_rates = []
for (i, rate) in enumerate(geol_slip_rate_vels)
    #fault_idx = parse(Int, rate.name)
    fault_idx = rate.name
    fault_row = @where(fault_df, :fid .== fault_idx)[1,:]
    fault = Oiler.IO.row_to_fault(fault_row; lsd="lower_seis_depth", 
                                             usd="upper_seis_depth")
    
    if haskey(results["poles"], (rate.fix, rate.mov))
        pred_rate = Oiler.Faults.get_fault_slip_rate_from_pole(fault, 
                        results["poles"][(rate.fix, rate.mov)]; 
                        lon=rate.lon, lat=rate.lat)
    else
        pred_rate = Oiler.Faults.get_fault_slip_rate_from_pole(fault,
                        results["poles"][(rate.mov, rate.fix)]; 
                        lon=rate.lon, lat=rate.lat)
    end
    
    push!(pred_geol_slip_rates, pred_rate)
end


inc = map(!, iszero.(geol_slip_rate_df[!,:include]))

function check_missing(val)
    if val == ""
        return false
    elseif val == 0.
        return false
    elseif ismissing(val) | isnothing(val)
        return false
    else
        return true
    end
end

dex_inc = [check_missing(d) for d in geol_slip_rate_df[!,:dextral_rate]]
ext_inc = [check_missing(d) for d in geol_slip_rate_df[!,:extension_rate]]

dex_geol_obs = geol_slip_rate_df[!,:dextral_rate][dex_inc .* inc]
#dex_geol_err = parse.(Float64, geol_slip_rate_df[!,:dextral_err][dex_inc .* inc])
dex_geol_err = geol_slip_rate_df[!,:dextral_err][dex_inc .* inc];
dex_geol_pred = [p[1] for p in pred_geol_slip_rates[dex_inc[inc]]];
dex_geol_pred_err = [p[3] for p in pred_geol_slip_rates[dex_inc[inc]]];

ext_geol_obs = geol_slip_rate_df[!,:extension_rate][ext_inc .* inc];
ext_geol_err = geol_slip_rate_df[!,:extension_err][ext_inc .* inc];
ext_geol_pred =  [p[2] for p in pred_geol_slip_rates[ext_inc[inc]]];
ext_geol_pred_err = [p[4] for p in pred_geol_slip_rates[ext_inc[inc]]];

all_obs = vcat(ove, ovn, dex_geol_obs, ext_geol_obs);
all_errs = vcat([v.ee for v in obs_vels], [v.en for v in obs_vels], 
                dex_geol_err, ext_geol_err);

all_pred = vcat(pve, pvn, dex_geol_pred, ext_geol_pred);

chi_sq = sum(((all_obs .- all_pred).^2 ) ./ all_errs.^2);

n_param = length(keys(vel_groups)) / 3.;

red_chi_sq = chi_sq / (length(all_obs) - n_param);


figure(figsize=(5, 9));
# suptitle("Observed vs. Modeled Quaternary Slip Rates")

subplot(2,1,1);
data_max = maximum([maximum(dex_geol_obs), maximum(dex_geol_pred)]);
data_min = minimum([minimum(dex_geol_obs), minimum(dex_geol_pred)]);

plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5);

axis("equal");
errorbar(dex_geol_obs, dex_geol_pred, xerr=dex_geol_err, yerr=dex_geol_pred_err,
         fmt=",", elinewidth=0.3);

autoscale(false);

fill_between([data_min, 0., -data_min], 
             [data_min, 0., data_min], 
             [data_min, data_min, data_min],
             color="grey",
             lw=0.,
             alpha=0.1);

fill_between([data_min, 0., -data_min], 
             [-data_min, -data_min, -data_min],
             [-data_min, 0., -data_min], 
             color="grey",
             lw=0.,
             alpha=0.1);

xlabel("observed");
ylabel("modeled");
title("dextral");
subplot(2,1,2);

data_max = maximum([maximum(ext_geol_obs), maximum(ext_geol_pred)]);
data_min = minimum([minimum(ext_geol_obs), minimum(ext_geol_pred)]);
plot([data_min, data_max], [data_min, data_max], "C1--", lw=0.5);

axis("equal")
errorbar(ext_geol_obs, ext_geol_pred, xerr=ext_geol_err, yerr=ext_geol_pred_err, 
         fmt=",", elinewidth=0.3);

autoscale(false);

fill_between([data_min, 0., -data_min], 
             [data_min, 0., data_min], 
             [data_min, data_min, data_min],
             color="grey",
             lw=0.,
             alpha=0.1);

fill_between([data_min, 0., -data_min], 
             [-data_min, -data_min, -data_min],
             [-data_min, 0., -data_min], 
             color="grey",
             lw=0.,
             alpha=0.1);


xlabel("observed");
ylabel("modeled");
title("extension");

tight_layout();


show();


#fault_pts_out = [Oiler.Faults.fault_to_vel_point(fault) 
#                 for fault in results["predicted_slip_rates"]]

#fault_out_df = DataFrame(fault_pts_out[1])
#for i in 2:length(fault_pts_out)
#    push!(fault_out_df, fault_pts_out[i])
#end

#CSV.write("/home/itchy/research/cascadia/cascadia_blocks/results/fault_slip_rate_pts.csv",
#          fault_out_df)

na_rel_poles = [Oiler.Utils.get_path_euler_pole(pole_arr, "na",
                                                string(block_df[i, :fid]))
                for i in 1:size(block_df, 1)]
if save_results
    CSV.write("../results/na_rel_poles.csv", 
            Oiler.IO.poles_to_df(na_rel_poles, convert_to_sphere=true))
end

println("done")
