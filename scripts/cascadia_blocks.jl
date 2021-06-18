using Revise

using CSV
using JSON
using DataFrames, DataFramesMeta
using PyPlot

using Oiler

midas_vel_file = "../data/midas_vels.geojson"
gsrm_vel_file = "../data/gsrm_na_rel.geojson"
blocks_file = "../data/cascadia_blocks.geojson"
idl_block_file = "../../../geodesy/global_blocks/data/idl_blocks.geojson"
s_us_block_file = "../../../us_faults/s_us_faults/s_us_blocks.geojson"
s_us_fault_file = "../../../us_faults/s_us_faults/s_us_faults.geojson"
cascadia_geol_slip_rates_file = "../data/cascadia_geol_slip_rates.geojson"

# tris_file = "../data/cascadia_subduction_tris.geojson"
tris_file = "../data/graham_cascadia_subduction_tris.geojson"
faults_file = "../data/cascadia_block_faults.geojson"

midas_df = Oiler.IO.gis_vec_file_to_df(midas_vel_file)
gsrm_df = Oiler.IO.gis_vec_file_to_df(gsrm_vel_file)
fault_df = Oiler.IO.gis_vec_file_to_df(faults_file)
s_us_fault_df = Oiler.IO.gis_vec_file_to_df(s_us_fault_file)
fault_df = vcat(fault_df, s_us_fault_df)

block_df = Oiler.IO.gis_vec_file_to_df(blocks_file)
idl_block_df = Oiler.IO.gis_vec_file_to_df(idl_block_file)
s_us_block_df = Oiler.IO.gis_vec_file_to_df(s_us_block_file)
block_df = vcat(block_df, idl_block_df, s_us_block_df)

geol_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cascadia_geol_slip_rates_file)

tri_json = JSON.parsefile(tris_file)


# load GNSS data
@info "getting GSRM idxs"
@time gsrm_block_idx = Oiler.IO.get_block_idx_for_points(gsrm_df, block_df, 2991)
@info "getting MIDAS idxs"
@time midas_block_idx = Oiler.IO.get_block_idx_for_points(midas_df, block_df, 2991)

function gsrm_vel_from_row(row, block)
    pt = Oiler.IO.get_coords_from_geom(row[:geometry])
    lon = pt[1]
    lat = pt[2]
    Oiler.VelocityVectorSphere(lon=lon, lat=lat, ve=row.e_vel,
        vn=row.n_vel, ee=row.e_err, en=row.n_err, name=row.station,
        fix="na", mov=string(block), vel_type="GNSS")
end

function midas_vel_from_row(row, block)
    pt = Oiler.IO.get_coords_from_geom(row[:geometry])
    lon = pt[1]
    lat = pt[2]
    Oiler.VelocityVectorSphere(lon=lon, lat=lat, 
        ve=row.na_ve_mm_yr,
        vn=row.na_vn_mm_yr,
        ee=row.na_ee_mm_yr,
        en=row.na_en_mm_yr,
        name=row.site,
        fix="na", mov=string(block), vel_type="GNSS")
end

gnss_vels = []

for (i, block) in enumerate(gsrm_block_idx)
    if !ismissing(block)
        gv = gsrm_vel_from_row(gsrm_df[i,:], block) 
        push!(gnss_vels, gv)
    end
end

for (i, block) in enumerate(midas_block_idx)
    if !ismissing(block)
        gv = midas_vel_from_row(midas_df[i,:], block) 
        push!(gnss_vels, gv)
    end
end
gnss_vels = convert(Array{VelocityVectorSphere}, gnss_vels)


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
jdf_vels = jdf_vels[1:10]



# do tris
function tri_from_feature(feat)
    tri = Oiler.Tris.Tri(
       p1=Float64.(feat["geometry"]["coordinates"][1][1]),
       p2=Float64.(feat["geometry"]["coordinates"][1][2]),
       p3=Float64.(feat["geometry"]["coordinates"][1][3]),
       name=string(feat["properties"]["fid"])
       )
end

tris = map(tri_from_feature, tri_json["features"])


# faults

dropmissing!(fault_df, :hw)
dropmissing!(fault_df, :fw)

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

function err_nothing_fix(err, return_val=2.)
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
        extension_err=err_nothing_fix(row[:e_ex]),
        dextral_rate=vel_nothing_fix(row[:v_rl]),
        dextral_err=err_nothing_fix(row[:e_rl]),
        dip=row[:dip],
        name=row[:fid],
        hw=row[:hw],
        fw=row[:fw],
        #usd=row[:usd],
        lsd=vel_nothing_fix(row[:lower_seis_depth]; return_val=25.),
        )
end

faults = [row_to_fault(fault_df[i,:]) for i in 1:size(fault_df, 1)]

faults = [fault for fault in faults if fault.hw != ""]
faults = [fault for fault in faults if fault.fw != ""]
faults = [fault for fault in faults if fault.fw != "c006"]
#faults = [fault for fault in faults if (fault.name == "cf197") ⊻ (fault.fw != "c006")]

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
    fault_idx = fault_seg #parse(Int, fault_seg)
    fault_row = @where(fault_df, :fid .== fault_idx)[1,:]
    fault = row_to_fault(fault_row)

    extension_rate = vel_nothing_fix(slip_rate_row[:extension_rate])
    extension_err = err_nothing_fix(slip_rate_row[:extension_err])#; return_val=5.)
    dextral_rate =vel_nothing_fix(slip_rate_row[:dextral_rate])
    dextral_err = err_nothing_fix(slip_rate_row[:dextral_err])#; return_val=5.)

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
        push!(geol_slip_rate_vels, make_vel_from_slip_rate(slip_rate_row, 
                                                           fault_df))
    end
end

geol_slip_rate_vels = convert(Array{VelocityVectorSphere}, geol_slip_rate_vels)

println("n fault slip rate vels: ", length(geol_slip_rate_vels))


vels = vcat(fault_vels, gnss_vels, jdf_vels, geol_slip_rate_vels)


vel_groups = Oiler.group_vels_by_fix_mov(vels)



"""
SOLVE
"""
# poles, tri_rates = 
results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
     tris=tris,
     #tris=[],
     faults=faults,
     weighted=true,
     sparse_lhs=true,
     predict_vels=true,
     check_closures=false,
     pred_se=false
    );

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

tri_rates = [get_tri_total_rate(tri) for tri in tris]
tri_rate_min = minimum(tri_rates)
tri_rate_max = maximum(tri_rates)

function plot_tri(tri; vmin=tri_rate_min, vmax=tri_rate_max)
    lons = [tri.p1[1], tri.p2[1], tri.p3[1], tri.p1[1]]
    lats = [tri.p1[2], tri.p2[2], tri.p3[2], tri.p1[2]]
    total_rate = get_tri_total_rate(tri)
    rate_frac = (total_rate - vmin) / (vmax - vmin)
    color = cm(rate_frac)
    #color = [0., 0., rate_frac]
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
show()

fault_pts_out = [Oiler.Faults.fault_to_vel_point(fault) 
                 for fault in results["predicted_slip_rates"]]

fault_out_df = DataFrame(fault_pts_out[1])
for i in 2:length(fault_pts_out)
    push!(fault_out_df, fault_pts_out[i])
end

CSV.write("/home/itchy/research/cascadia/cascadia_blocks/results/fault_slip_rate_pts.csv",
          fault_out_df)

na_rel_poles = [Oiler.Utils.get_path_euler_pole(pole_arr, "na",
                                                string(block_df[i, :fid]))
                for i in 1:size(block_df,1)]
CSV.write("../results/na_rel_poles.csv", 
            Oiler.IO.poles_to_df(na_rel_poles, convert_to_sphere=true))

println("done")
