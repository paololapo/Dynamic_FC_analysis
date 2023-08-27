using NPZ
using Plots
using LaTeXStrings
using LinearAlgebra

#Define plot
default(fontfamily="Computer Modern",
        linewidth=1, framestyle=:box, label=nothing, grid=true)
#scalefontsizes(0.5)

#Import table
#table = npzread("table_mid.npy")
classification = npzread("./files/classification.npy")
#sigma = range(0, stop=0.7, length=10)


function plot_table(table, classification)
    sigma = range(0, stop=0.7, length=length(table[:, 1]))
    p = plot(guidefontsize=12, legendfontsize=12, dpi=500)
    for i in 1:96
        if classification[i] == 0
            plot!(p, sigma, table[:, i], color="green")
        else
            plot!(p, sigma, table[:, i], color="lightblue")
        end
    end

    plot!(p, sigma, table[:, 1], color="green", label="Limit cycle")
    plot!(p, sigma, table[:, 3], color="lightblue", label="Fixed point")

    plot!(ylabel=L"\mathrm{KS}(V_{\mathrm{MINDy}}, V)", xlabel=L"\sigma")

    return p
end


#savefig(p, "./plots/sigma_julia.png")
#display(p)

p1 = plot_table(npzread("./files/table_short.npy"), classification)
title!(p1, "Short")
p2 = plot_table(npzread("./files/table_mid.npy"), classification)
title!(p2, "Mid")
p3 = plot_table(npzread("./files/table_long.npy"), classification)
title!(p3, "Long")

plt = plot(p1, p2, p3, layout=(1, 3), size=(2000, 600))
plot!(plt, left_margin=13Plots.mm, bottom_margin=10Plots.mm, top_margin=5Plots.mm)
#savefig(plt, "./plots/sigma_dFC.png")
