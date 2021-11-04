using LinearAlgebra

export construct_oneloop_curve

function construct_oneloop_curve(x0,y0,ds,length,gap,angle)
    x, y = Float64[], Float64[]

    xt1 = x0 + length/2
    xt2 = x0 - length/2
    yt1 = y0 + gap/2
    yt2 = y0 - gap/2


    xi,yi = _straight_line(xt1,yt1,xt2,yt1,ds)
    append!(x,xi);append!(y,yi)
    xi,yi = _semi_circle(xt2,yt1,xt2,yt2,ds)
    append!(x,xi);append!(y,yi)
    xi,yi = _straight_line(xt2,yt2,xt1,yt2,ds)
    append!(x,xi);append!(y,yi)
    xi,yi = _semi_circle(xt1,yt2,xt1,yt1,ds)
    append!(x,xi);append!(y,yi)

    x , y, x[1], y[1]
end

function _straight_line(x1,y1,x2,y2,ds)
    dis = norm([x1-x2, y1-y2],2)
    num_pts = Int(cld(dis,ds))

    x = Array(LinRange(x1, x2, num_pts))
    y = Array(LinRange(y1, y2, num_pts))

    x,y
end

function _semi_circle(x1,y1,x2,y2,ds)
    radius = norm([x1-x2, y1-y2],2)/2
    center = ([x1,y1] + [x2,y2])   /2
    num_pts = Int(cld(radius*pi,ds))

    vector = [x1,y1] - center

    θ = Array(LinRange(0, pi, num_pts))


    x = center[1] .+ vector[1] * cos.(θ) .- vector[2] * sin.(θ)
    y = center[2] .+ vector[1] * sin.(θ) .+ vector[2] * cos.(θ)

    x,y
end
