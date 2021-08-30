using Revise

using CSV
using JSON
using DataFrames, DataFramesMeta
using PyPlot
using Setfield

using Oiler

midas_vel_file = "../data/midas_vels.geojson"
gsrm_vel_file = "../data/gsrm_na_rel.geojson"
blocks_file = "../data/cascadia_blocks.geojson"
idl_block_file = "../../../geodesy/global_blocks/data/idl_blocks.geojson"
s_us_block_file = "../../../us_faults/s_us_faults/s_us_blocks.geojson"
s_us_fault_file = "../../../us_faults/s_us_faults/s_us_faults.geojson"
cascadia_geol_slip_rates_file = "../data/cascadia_geol_slip_rates.geojson"
new_us_geol_rates_file = "../../../us_faults/s_us_faults/new_us_faults_geol_slip_rates.geojson"

trench_faults_file = "../data/cascadia_trench_faults.geojson"

# tris_file = "../data/cascadia_subduction_tris.geojson"
tris_file = "../data/graham_cascadia_subduction_tris.geojson"
faults_file = "../data/cascadia_block_faults.geojson"

midas_df = Oiler.IO.gis_vec_file_to_df(midas_vel_file)
gsrm_df = Oiler.IO.gis_vec_file_to_df(gsrm_vel_file)
fault_df = Oiler.IO.gis_vec_file_to_df(faults_file)
s_us_fault_df = Oiler.IO.gis_vec_file_to_df(s_us_fault_file)

trench_fault_df = Oiler.IO.gis_vec_file_to_df(trench_faults_file)

fault_df = vcat(fault_df, s_us_fault_df, trench_fault_df)

block_df = Oiler.IO.gis_vec_file_to_df(blocks_file)
idl_block_df = Oiler.IO.gis_vec_file_to_df(idl_block_file)
s_us_block_df = Oiler.IO.gis_vec_file_to_df(s_us_block_file)
block_df = vcat(block_df, idl_block_df, s_us_block_df)

#geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cascadia_geol_slip_rates_file)
geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(new_us_geol_rates_file)

tri_json = JSON.parsefile(tris_file)


# load GNSS data

@info "doing GSRM GNSS"
@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_df, block_df;
                                                           fix="na",
                                                           ve=:e_vel,
                                                           vn=:n_vel,
                                                           ee=:e_err,
                                                           en=:n_err,
                                                           name=:station,
                                                           epsg=2991)

@info "doing MIDAS GNSS"
@time midas_vels = Oiler.IO.make_vels_from_gnss_and_blocks(midas_df, block_df;
                                                            fix="na",
                                                            ve=:na_ve_mm_yr,
                                                            vn=:na_vn_mm_yr,
                                                            ee=:na_ee_mm_yr,
                                                            en=:na_en_mm_yr,
                                                            name=:site,
                                                            epsg=2991)

gnss_vels = vcat(gsrm_vels, midas_vels)


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


jdf_vels = [Oiler.VelocityVectorSphere(vel; vel_type="GNSS") for vel in jdf_vels]
jdf_vels = jdf_vels



tris = Oiler.IO.tris_from_geojson(tri_json)

function set_tri_rates(tri; ds=25., de=1., ss=0., se=0.5)
    tri = @set tri.dip_slip_rate = ds
    tri = @set tri.dip_slip_err = de
    tri = @set tri.strike_slip_rate = ss
    tri = @set tri.strike_slip_err = se
    tri
end

tris = map(set_tri_rates, tris)

# faults

dropmissing!(fault_df, :hw)
dropmissing!(fault_df, :fw)

fault_df = filter(row -> !(row.hw == ""), fault_df)
fault_df = filter(row -> !(row.fw == ""), fault_df)


fault_weight = 1.

function vel_nothing_fix(vel; return_val=0.)
    if vel == ""
        return return_val
    elseif isnothing(vel) | ismissing(vel)
        return return_val
    else
        if typeof(vel) == String
            return parse(Float64, vel)
        else
            return vel
        end
    end
end

function err_nothing_fix(err; return_val=2.)
    if err == ""
        return return_val
    elseif isnothing(err) | ismissing(err) | (err == 0.)
        return return_val
    else
        if typeof(err) == String
            err = parse(Float64, err)
            if iszero(err)
                return return_val
            else
                return err
            end
        else
            return err
        end
    end
end


function row_to_fault(row)
    trace = Oiler.IO.get_coords_from_geom(row[:geometry])

    Oiler.Fault(trace=trace,
        dip_dir=row[:dip_dir],
        extension_rate=vel_nothing_fix(row[:v_ex]),
        extension_err=err_nothing_fix(row[:e_ex], return_val=0.5),
        dextral_rate=vel_nothing_fix(row[:v_rl]),
        dextral_err=err_nothing_fix(row[:e_rl], return_val=0.5),
        dip=row[:dip],
        fid=row[:fid],
        hw=row[:hw],
        fw=row[:fw],
        #name=row[:name],
        # usd=row[:usd],
        lsd=vel_nothing_fix(row[:lower_seis_depth]; return_val=15.),
        )
end

faults = [row_to_fault(fault_df[i,:]) for i in 1:size(fault_df, 1)]
# faults = [Oiler.IO.row_to_fault(fault_df[i,:]; usd=:upper_seis_depth, lsd=:lower_seis_depth) for i in 1:size(fault_df, 1)]

# faults = []
# for i in 1:size(fault_df, 1)
#    try
#        push!(faults, Oiler.IO.row_to_fault(fault_df[i,:]; 
#                                            usd=:upper_seis_depth, 
#                                            lsd=:lower_seis_depth)
#             )
#    catch
#        println(fault_df[i,:])
#    end
# end


#faults = [fault for fault in faults if fault.fw != "c006"]
faults = [fault for fault in faults if fault.hw != "c_07"]
# faults = [fault for fault in faults if (fault.name == "cf197") ⊻ (fault.fw != "c006")]

# faults = map(feat_to_Fault, fault_json["features"]);
# fault_vels_ = map(Oiler.fault_to_vels, faults);

fault_vels = []
for fault in faults
    try
        vel = Oiler.fault_to_vel(fault)
        push!(fault_vels, vel)
    catch
        println(fault.name)
    end
end

fault_vels = convert(Array{VelocityVectorSphere}, fault_vels)

# fault_vels_ = map(Oiler.fault_to_vel, faults);
# fault_vels = reduce(vcat, fault_vels_)
println("n faults: ", length(faults))

# geol slip rates
function make_vel_from_slip_rate(slip_rate_row, fault_df)
    fault_seg = slip_rate_row[:fault_seg]
    fault_idx = fault_seg # parse(Int, fault_seg)
    fault_row = @where(fault_df, :fid .== fault_idx)[1,:]
    fault = row_to_fault(fault_row)

    extension_rate = vel_nothing_fix(slip_rate_row[:extension_rate])
    extension_err = err_nothing_fix(slip_rate_row[:extension_err])# ; return_val=5.)
    dextral_rate = vel_nothing_fix(slip_rate_row[:dextral_rate])
    dextral_err = err_nothing_fix(slip_rate_row[:dextral_err])# ; return_val=5.)

    ve, vn = Oiler.Faults.fault_slip_rate_to_ve_vn(dextral_rate, 
                                                   extension_rate,
                                                   fault.strike)

    ee, en, cen = Oiler.Faults.fault_slip_rate_err_to_ee_en(dextral_err, 
                                                            extension_err,
                                                            fault.strike)

    pt = Oiler.IO.get_coords_from_geom(slip_rate_row[:geometry])
    lon = pt[1]
    lat = pt[2]
    
    VelocityVectorSphere(lon=lon, lat=lat, ve=ve, vn=vn, fix=fault.hw,
                         ee=ee, en=en, cen=cen,
                         mov=fault.fw, vel_type="fault", name=fault_seg)

end


geol_slip_rate_vels = []

for i in 1:size(geol_slip_rate_df, 1)
    slip_rate_row = geol_slip_rate_df[i,:]
    if (slip_rate_row[:include] == true) | (slip_rate_row[:include] == "1")
        try
            push!(geol_slip_rate_vels, make_vel_from_slip_rate(slip_rate_row, 
                                                           fault_df))
        catch
            println(slip_rate_row[:fault_seg])
        end
    end
end

geol_slip_rate_vels = convert(Array{VelocityVectorSphere}, geol_slip_rate_vels)

println("n fault slip rate vels: ", length(geol_slip_rate_vels))


vels = vcat(fault_vels, gnss_vels, jdf_vels, geol_slip_rate_vels)

vels = [vel for vel in vels if (vel.fix != "")]
vels = [vel for vel in vels if (vel.mov != "")]


vel_groups = Oiler.group_vels_by_fix_mov(vels)



"""
SOLVE
"""
# poles, tri_rates = 
results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
     #tris=tris,
     tris=[],
     faults=faults,
     tri_distance_weight=50.,
     tri_priors=true,
     weighted=true,
     sparse_lhs=true,
     predict_vels=true,
     check_closures=false,
     pred_se=true,
     constraint_method="kkt",
    );


Oiler.IO.write_fault_results_to_gj(results, 
    "../results/w_us_cascadia_fault_results.geojson";
    name="Western North America faults")

#Oiler.IO.write_tri_results_to_gj(tris, results,
#    "../results/cascadia_tri_results_sm.geojson";
#    name="Cascadia tri rates")

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

#cm = get_cmap(:viridis)
#
#function get_tri_total_rate(tri)
#    ds = results["tri_slip_rates"][tri.name]["dip_slip"]
#    ss = results["tri_slip_rates"][tri.name]["strike_slip"]
#    total_rate = sqrt(ds^2 + ss^2)
#end
#
#if length(tris) > 0
#    tri_rates = [get_tri_total_rate(tri) for tri in tris]
#    tri_rate_min = minimum(tri_rates)
#    tri_rate_max = maximum(tri_rates)
#end
#
#function plot_tri(tri; vmin=tri_rate_min, vmax=tri_rate_max)
#    lons = [tri.p1[1], tri.p2[1], tri.p3[1], tri.p1[1]]
#    lats = [tri.p1[2], tri.p2[2], tri.p3[2], tri.p1[2]]
#    total_rate = get_tri_total_rate(tri)
#    rate_frac = (total_rate - vmin) / (vmax - vmin)
#    color = cm(rate_frac)
#    # color = [0., 0., rate_frac]
#    fill(lons, lats, color=color, alpha=0.25, zorder=0)
#end
#
#for tri in tris
#    plot_tri(tri)
#end


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
    #fault = Oiler.IO.row_to_fault(fault_row)
    fault = row_to_fault(fault_row)
    
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
CSV.write("../results/na_rel_poles.csv", 
            Oiler.IO.poles_to_df(na_rel_poles, convert_to_sphere=true))

println("done")
