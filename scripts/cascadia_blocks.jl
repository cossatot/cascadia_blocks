using Revise

using CSV
using JSON
using DataFrames
using PyPlot

using Oiler

midas_vel_file = "../data/midas_vels.geojson"
gsrm_vel_file = "../data/gsrm_na_rel.geojson"
blocks_file = "../data/cascadia_blocks.geojson"
tris_file = "../data/cascadia_subduction_tris.geojson"
faults_file = "../data/cascadia_block_faults.geojson"

midas_df = Oiler.IO.gis_vec_file_to_df(midas_vel_file)
gsrm_df = Oiler.IO.gis_vec_file_to_df(gsrm_vel_file)
block_df = Oiler.IO.gis_vec_file_to_df(blocks_file)
fault_df = Oiler.IO.gis_vec_file_to_df(faults_file)

tri_json = JSON.parsefile(tris_file)


# load GNSS data
gsrm_block_idx = Oiler.IO.get_block_idx_for_points(gsrm_df, block_df)
midas_block_idx = Oiler.IO.get_block_idx_for_points(midas_df, block_df)

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
jdf_pts = CSV.read(jdf_pt_file)

# From DeMets et al 2010, pre-converted to Cartesian
jdf_na_pole = Oiler.PoleCart(
    x=5.915996643479488e-9,
    y=1.486624308775784e-8,
    z=-9.997991640995085e-9,
    # ex=sqrt(103.05e-8),
    # ey=sqrt(186.44e-8),
    # ez=sqrt(259.51e-8),
    # cxy=135.05e-8,
    # cxz=-155.74e-8,
    # cyz=-203.56e-8,
    fix="na",
    mov="c006"
)

jdf_vels = Oiler.BlockRotations.predict_block_vels(jdf_pts[:,:lon], 
    jdf_pts[:,:lat], jdf_na_pole)


jdf_vels = [Oiler.VelocityVectorSphere(vel; vel_type="GNSS") for vel in jdf_vels]

# do tris
function tri_from_feature(feat)
    tri = Oiler.Tris.Tri(
       p1=Float64.(feat["geometry"]["coordinates"][1][1]),
       p2=Float64.(feat["geometry"]["coordinates"][1][2]),
       p3=Float64.(feat["geometry"]["coordinates"][1][3]),
       name=feat["properties"]["fid"]
       )
end

tris = map(tri_from_feature, tri_json["features"])


# faults

dropmissing!(fault_df, :hw)
dropmissing!(fault_df, :fw)

fault_weight = 1.


function vel_nothing_fix(vel)
    if vel == ""
        return 0.
    elseif isnothing(vel) | ismissing(vel)
        return 0.
    else
        if typeof(vel) == String
            return parse(Float64, vel)
        else
            return vel
        end
    end
end

function err_nothing_fix(err)
    if err == ""
        return 20.
    elseif isnothing(err) | ismissing(err)
        return 20.
    else
        if typeof(err) == String
            err = parse(Float64, err)
            if iszero(err)
                return 20.
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
        # usd=row[:usd],
        # lsd=row[:lsd],
        )
end

faults = [row_to_fault(fault_df[i,:]) for i in 1:size(fault_df, 1)]

faults = [fault for fault in faults if fault.hw != ""]
faults = [fault for fault in faults if fault.fw != ""]

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

vels = vcat(fault_vels, gnss_vels, jdf_vels)


vel_groups = Oiler.group_vels_by_fix_mov(vels)



"""
SOLVE
"""
# poles, tri_rates = 
results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
     tris=tris,
     faults=faults,
     weighted=true,
     sparse_lhs=true,
     predict_vels=true,
     check_closures=false,
     pred_se=true
    );


obs_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(vel_groups)]
pred_vels = [v["vel"] for v in Oiler.Utils.get_gnss_vels(results["predicted_vels"])]

glons, glats = Oiler.Utils.get_coords_from_vel_array(obs_vels)
ove = [v.ve for v in obs_vels]
ovn = [v.vn for v in obs_vels]

pve = [v.ve for v in pred_vels]
pvn = [v.vn for v in pred_vels]

figure(figsize=(14, 14))

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
