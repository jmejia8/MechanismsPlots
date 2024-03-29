function generateTrayectory(p, nframes=30)
    # Save mechanism dynamics 
    if length(p) == 6
        r1,r2,r3,r4,rcx,rcy = p
        θ0,x0,y0 = 0.0, 0.0, 0.0
    else
        r1,r2,r3,r4,rcx,rcy,θ0,x0,y0 = p
    end
 
    θ1 = 0

    X0 = zeros(nframes, 2)
    X1 = zeros(nframes, 2)
    X2 = zeros(nframes, 2)
    X3 = zeros(nframes, 2)
    
    C = zeros(nframes, 2)
    i = 1
    for θ2 = range(0, 2π, length = nframes)

        A1 = 2p[3] * (p[2] * cos(θ2) - p[1]*cos(θ1))
        B1 = 2p[3] * (p[2] * sin(θ2) - p[1]*sin(θ1))
        C1 = p[1]^2+p[2]^2+p[3]^2-p[4]^2- 2p[1]*p[2]*cos(θ2 - θ1)

        d = B1^2+A1^2-C1^2

        if d < 0
            error("It is not a four-bar mechanism.")
        end

        θ3 = 2atan((-B1+sqrt(d))/(C1-A1))


        Cxr = r2*cos(θ2)+p[5]*cos(θ3) - rcy*sin(θ3)
        Cyr = r2*sin(θ2)+p[5]*sin(θ3) + rcy*cos(θ3)

        C[i, 1] = Cxr*cos(θ0)-Cyr*sin(θ0) + x0;
        C[i, 2] = Cxr*sin(θ0)+Cyr*cos(θ0) + y0;
        
        X0[i, 1] = x0
        X0[i, 2] = y0

        X1[i, 1] = x0 + r1*cos(θ0)
        X1[i, 2] = y0 + r1*sin(θ0)

        X2[i, 1] = x0 + r2*cos(θ2+θ0)
        X2[i, 2] = y0 + r2*sin(θ2+θ0)

        X3[i, 1] = X2[i, 1] + r3*cos(θ3+θ0)
        X3[i, 2] = X2[i, 2] + r3*sin(θ3+θ0)
        
        i += 1
    end

    return C, X0, X1, X2, X3

end

function animate(p; precision_points=nothing,
                                    xlimits=nothing,
                                    ylimits=nothing,
                                    title="Four-Bar Mechanism",
                                    nframes=50)

    C, X0, X1, X2, X3 = generateTrayectory(p, nframes)

    if xlimits == nothing
        xmin = min( minimum(C[:, 1]), minimum(X0[:, 1]), minimum(X1[:, 1]), minimum(X2[:, 1]), minimum(X3[:, 1]), )
        xmax = max( maximum(C[:, 1]), maximum(X0[:, 1]), maximum(X1[:, 1]), maximum(X2[:, 1]), maximum(X3[:, 1]), )
    else
        xmin, xmax = xlimits
    end

    if ylimits == nothing
        ymin = min( minimum(C[:, 2]), minimum(X0[:, 2]), minimum(X1[:, 2]), minimum(X2[:, 2]), minimum(X3[:, 2]), )
        ymax = max( maximum(C[:, 2]), maximum(X0[:, 2]), maximum(X1[:, 2]), maximum(X2[:, 2]), maximum(X3[:, 2]), )

    else
        ymin, ymax = ylimits
    end

    nx = xmax - xmin
    ny = ymax - ymin

    if xlimits == nothing && nx < ny
        offset = (ny - nx) / 2.0
        xmin -= offset
        xmax += offset 
    end
    if ylimits == nothing && nx > ny
        offset = abs(ny - nx) / 2.0
        ymin -= offset
        ymax += offset 
    end

    xlimits = [xmin, xmax]
    ylimits = [ymin, ymax]


    anim = @animate for t = 1:nframes
        # plot()
        plot(C[:,1], C[:,2], linestyle=:dot, linecolor=:green, title=title, xlimits=xlimits, ylimits=ylimits)

        if precision_points != nothing
            scatter!(precision_points[:,1], precision_points[:,2], markercolor=:black, markershape=:star5)
        end
        
        scatter!(X0[t, 1:1], X0[t, 2:2])
        scatter!(X1[t, 1:1], X1[t, 2:2])
        scatter!(X2[t, 1:1], X2[t, 2:2])
        scatter!(X3[t, 1:1], X3[t, 2:2])
        
        scatter!(C[t, 1:1], C[t, 2:2], markercolor=:red)
       

        plot!([X0[t, 1] , X2[t, 1] ], [X0[t, 2] , X2[t, 2] ])
        plot!([X1[t, 1] , X0[t, 1] ], [X1[t, 2] , X0[t, 2] ])
        plot!([X1[t, 1] , X3[t, 1] ], [X1[t, 2] , X3[t, 2] ])
        
        plot!([X2[t, 1] , X3[t, 1] ], [X2[t, 2] , X3[t, 2] ], linecolor=:red)
        plot!([X2[t, 1] , C[t, 1]], [X2[t, 2] , C[t, 2]], linecolor=:red)
        plot!([X3[t, 1] , C[t, 1]], [X3[t, 2] , C[t, 2]], linecolor=:red) 
    end

    anim
end
