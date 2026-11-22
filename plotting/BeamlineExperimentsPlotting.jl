module BeamlineExperimentsPlotting

using PALSParserJ, AcceleratorLatticePlot, Makie
using GeometryBasics: Polygon

export interactive_apertures, floor_view, default_lattice

default_lattice() = normpath(joinpath(@__DIR__, "..", "src", "bta_lattice", "bta.pals.yaml"))

const ALP = AcceleratorLatticePlot

function ribbon_geometry(tab, records, scale)
function edges(i, f, lo, hi)
    placement = place(tab, i, f)
    c = (placement.r[3], placement.r[1])
    normal = (placement.W.m31, placement.W.m11)
    a = Point2f(c[1] + scale*lo*normal[1], c[2] + scale*lo*normal[2])
    b = Point2f(c[1] + scale*hi*normal[1], c[2] + scale*hi*normal[2])
    # Verify transverse scaling against the unrounded construction.
    @assert isapprox(hypot(normal...), 1.; atol=1e-6)
    @assert isapprox(hypot(Float64(b[1]-a[1]), Float64(b[2]-a[2])), scale*(hi-lo); atol=1e-5)
    return a, b
end
nanpoint = Point2f(NaN, NaN)
lowedge, highedge, checks = Point2f[], Point2f[], Point2f[]
polygons = Polygon{2,Float32}[]
profile_polys = Polygon{2,Float32}[]
profile_lo, profile_hi = Point2f[], Point2f[]
for r in records
    i, lo, hi = r.i, r.lo, r.hi
    # Fine samples also follow short bends smoothly; do not join across patches.
    fs = range(0., 1.; length=max(3, ceil(Int, abs(tab.angle[i])/0.005)+1))
    pairs = [edges(i, f, lo, hi) for f in fs]
    lower, upper = first.(pairs), last.(pairs)
    if tab.length[i] > 0
        push!(polygons, Polygon(vcat(lower, reverse(upper))))
        append!(lowedge, lower); push!(lowedge, nanpoint)
        append!(highedge, upper); push!(highedge, nanpoint)
        s0, s1 = tab.s[i], tab.s[i]+tab.length[i]
        push!(profile_polys, Polygon(Point2f[(s0,lo*1e3),(s1,lo*1e3),(s1,hi*1e3),(s0,hi*1e3)]))
        append!(profile_lo, [Point2f(s0,lo*1e3),Point2f(s1,lo*1e3),nanpoint])
        append!(profile_hi, [Point2f(s0,hi*1e3),Point2f(s1,hi*1e3),nanpoint])
    end
    locations = r.location == "BOTH_ENDS" ? (0.,1.) :
                r.location == "EXIT_END" ? (1.,) :
                r.location == "CENTER" ? (0.5,) :
                r.location == "EVERYWHERE" ? (0.,0.5,1.) : (0.,)
    for f in locations
        a,b = edges(i,f,lo,hi)
        append!(checks, [a,b])
    end
end

return (; polygons, lowedge, highedge, checks, profile_polys, profile_lo, profile_hi)
end


"""
    interactive_apertures(input; scale=30, show=true)

Pan/zoom the nominal aperture ribbons and adjust transverse exaggeration.
The input is a PALS lattice file. Returns figure, axes, and slider controls.
"""
function interactive_apertures(input=default_lattice(); scale=30, show=true)
    1 <= scale <= 100 || error("scale must be between 1 and 100")
    lat = parse_and_expand_pals(String(input); problems=:none)
    isempty(lat.problems) || error("PALS expansion problems: $(lat.problems)")
    tab = element_table(lat)
num(n, key, default) = haskey(n, key) ? parse(Float64, String(n[key])) : default
str(n, key, default) = haskey(n, key) ? String(n[key]) : default
records = NamedTuple[]
for i in eachindex(tab.name)
    n = tab.node[i]
    haskey(n, "ApertureP") || continue
    ap = n["ApertureP"]
    lowercase(str(ap, "aperture_active", "true")) == "false" && continue
    location = str(ap, "location", "ENTRANCE_END")
    location == "NOWHERE" && continue
    lo = num(ap, "x_min", -Inf)
    hi = num(ap, "x_max", Inf)
    if haskey(ap, "x_width")
        width, center = num(ap, "x_width", Inf), num(ap, "x_center", 0.)
        lo, hi = center - width/2, center + width/2
    end
    isfinite(lo) && isfinite(hi) || continue
    lo <= hi || error("Invalid horizontal aperture on $(tab.name[i])")
    push!(records, (; i, lo, hi, location))
end
isempty(records) && error("No finite active horizontal apertures found")


    fig = Figure(size=(1500,1000), fontsize=17)
    Label(fig[0,1], "BeamlineExperiments · Interactive horizontal apertures", fontsize=27, font=:bold, halign=:left)
    Label(fig[1,1], "Nominal limits around the curved reference frame; body misalignments omitted.",
          fontsize=15, color=:gray35, halign=:left)
    ax = Axis(fig[2,1], aspect=DataAspect(), xlabel="Global z [m]", ylabel="Global x [m]")
    control = GridLayout(fig[3,1])
    Label(control[1,1], "Transverse exaggeration")
    slider = Slider(control[1,2], range=1:100, startvalue=round(Int,scale))
    Label(control[1,3], lift(v -> "$(round(Int,v))×", slider.value))
    reset = Button(control[1,4], label="Reset view")
    data = lift(s -> ribbon_geometry(tab, records, s), slider.value)
    poly!(ax, lift(d -> d.polygons, data); color=(:teal,0.18), strokewidth=0)
    lines!(ax, lift(d -> d.lowedge, data); color=:slateblue, linewidth=2)
    lines!(ax, lift(d -> d.highedge, data); color=:teal, linewidth=2)
    scatter!(ax, lift(d -> d.checks, data); color=:gray25, markersize=4)
    centerline = build_geometry(tab, ShapeMap([ele_shape("*", :none, :gray30; label=:none)]); view="zx", arc_tol=0.005)
    lines!(ax, centerline.ref_pts; color=:gray20, linestyle=:dash, linewidth=1.5)
    Legend(fig[4,1], [LineElement(color=:gray20,linestyle=:dash),
        LineElement(color=:teal,linewidth=3), LineElement(color=:slateblue,linewidth=3),
        MarkerElement(color=:gray25,marker=:circle)],
        ["Reference centerline", "+x edge", "−x edge", "Aperture-check positions"];
        orientation=:horizontal, framevisible=false, tellwidth=false)
    profile = Axis(fig[5,1], xlabel="Distance along reference path s [m]", ylabel="Local x [mm]",
        title="Actual limits · unchanged by the exaggeration slider", titlealign=:left)
    initial = data[]
    poly!(profile, initial.profile_polys; color=(:teal,0.15), strokewidth=0)
    lines!(profile, initial.profile_lo; color=:slateblue, linewidth=1.5)
    lines!(profile, initial.profile_hi; color=:teal, linewidth=1.5)
    hlines!(profile,[0.]; color=:gray35, linestyle=:dash)
    xlims!(profile,0,maximum(tab.s))
    rowsize!(fig.layout,5,180)
    Label(fig[6,1], "Scroll: zoom   •   Right-drag: pan   •   Left-drag: rectangle zoom   •   Reset view: fit current apertures",
        fontsize=15, halign=:left)
    Label(fig[7,1], "Bands extend stated limits along each element for illustration; missing apertures remain gaps. Dots mark check positions.",
        fontsize=14, color=:gray35, halign=:left)
    function fitview!()
        margin = maximum(max(abs(r.lo),abs(r.hi)) for r in records)*slider.value[]+1
        xlims!(ax,minimum(tab.z)-margin,maximum(tab.z)+margin)
        ylims!(ax,minimum(tab.x)-margin,maximum(tab.x)+margin)
    end
    on(_ -> fitview!(), reset.clicks)
    fitview!()
    show && display(fig)
    return (; figure=fig, axis=ax, profile, slider, data, reset, records, tab)
end


"""Create the regular machine floor plan from a PALS file using the active Makie backend."""
function floor_view(input=default_lattice(); kwargs...)
    lat = parse_and_expand_pals(String(input); problems=:none)
    isempty(lat.problems) || error("PALS expansion problems: $(lat.problems)")
    return floor_plot(lat; kwargs...)
end

end # module
