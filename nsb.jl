using LibPQ
using DataFrames

# Connect to NSB
function nsbConnect(username, password, inMfn=false)
    host = inMfn ? "192.168.101.168" : "212.201.100.111"
    con = LibPQ.Connection("dbname=nsb host=$host port=5432 user=$username password=$password")
    return con
end

# Compute age(s) at given depth on a given site.
function findAge(con, hole_id, x)
    nam = DataFrame(execute(con, "SELECT a.depth_mbsf, a.age_ma FROM neptune_age_model as a, neptune_age_model_history as b, neptune_hole_summary as c WHERE a.site_hole=b.site_hole AND a.site_hole=c.site_hole AND a.revision_no=b.revision_no AND b.current_flag='Y' AND c.hole_id='$hole_id' ORDER BY a.depth_mbsf, a.age_ma;"))
    depth = map(Float64,nam[!,1])
    age = map(Float64,nam[!,2])
    res = [NaN for i=1:length(x)]
    for i = 1:length(x)
        if (x[i] .<= depth[end]) & (x[i].>=depth[1]) & (length(depth) > 0)
            x1 = depth[depth .<= x[i]][end]
            x2 = depth[depth .>= x[i]][1]
            y1 = age[depth .<= x[i]][end]
            y2 = age[depth .>= x[i]][1]
            res[i] = y1 + (x[i] - x1)*(y2 - y1)/(x2 - x1)
        end
    end
    return res
end

# Linear Sedimentation Rates
function lsr(con, hole_id, mbsf=nothing)
    nam = DataFrame(execute(con,"SELECT b.depth_mbsf, b.age_ma FROM neptune_hole_summary as a, neptune_age_model as b, neptune_age_model_history as c WHERE a.site_hole=b.site_hole AND b.site_hole=c.site_hole AND b.revision_no=c.revision_no AND c.current_flag='Y' AND a.hole_id='$hole_id' ORDER BY b.depth_mbsf, b.age_ma;"))
    depth = map(Float64,nam.depth_mbsf)
    age = map(Float64,nam.age_ma)
    from = depth[1:(end-1)]
    to = depth[2:end]
    b = (to-from)./(age[2:end]-age[1:(end-1)])/10
    df = DataFrame(hole_id=hole_id, from=from, to=to, lsr=b)
    df = df[df.lsr .> 0.00,:]
    if mbsf==nothing
        return df
    elseif mbsf > depth[end] || mbsf < depth[1]
        return missing
    else
        return df.lsr[(mbsf .>= df.from) .& (mbsf .<= df.to)]
    end
end

##Test
# Connect as user "guest"
nsb = nsbConnect("guest", "arm_aber_sexy")

# Compute age for site 748B at depths 10 to 60
depths = [10 20 30 40 50 60]
age = findAge(nsb, "120_748B", depths)

#Grab the table of linear sedimentation rates implied by the age model
lsr_table = lsr(nsb,"120_748B")

# Compute the LSR at depth 60
lsr_at_60 = lsr(nsb,"120_748B",60)