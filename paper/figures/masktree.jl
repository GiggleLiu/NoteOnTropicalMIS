using UnitDiskMapping.TikzGraph

function fig1()  # mask tree
    filename = joinpath(@__DIR__, "masktree.tex")
    graph = canvas(props=Dict("scale"=>"1.5")) do c
        node!(x, y) = Node(x, y; minimum_size="0.3cm") >> c
        dx1 = 0.8
        dx2 = 0.5
        dx3 = 0.3
        locs = [(0, 0), (-dx1, 0.6), (dx1, 0.6), (-dx1-dx2, 1.2), (-dx1+dx2, 1.2), (dx1+dx2, 1.2), (dx1-dx2, 1.2), (-dx1-dx2-dx3, 1.8), (-dx1-dx2+dx3, 1.8)]
        dx = dy = 0.0
        nodes = []
        for (x,y) in locs
            push!(nodes, node!(x+dx, -y-dy))
        end
        for (i, j) in [(1, 2), (1,3), (2,4), (2,5), (3,6), (3,7), (4,8), (4,9)]
            Line(nodes[i], nodes[j], arrow="<-", line_width="0.5pt") >> c
        end
        annotate(nodes[1], "\$\\alpha(G)\$", offsety=0.3) >> c
        annotate(nodes[4], "\$A\$", offsetx=-0.3) >> c
        annotate(nodes[5], "\$B\$", offsetx=0.0, offsety=-0.3) >> c
        annotate(nodes[2], "\$C\$", offsety=0.3) >> c
        annotate(nodes[2], "\$*\$", offsety=-0.3) >> c
        PlainText(dx, -2.2, "\\large (a)") >> c

        dx += 4
        nodes = []
        for (x,y) in locs
            push!(nodes, node!(x+dx, -y-dy))
        end
        for (i, j) in [(1, 2), (1,3), (2,4), (2,5), (3,6), (3,7), (4,8), (4,9)]
            Line(nodes[i], nodes[j], arrow="->", line_width="0.5pt", draw="red") >> c
        end
        annotate(nodes[1], "\$\\overline{\\alpha(G)}\$", offsety=0.3) >> c
        annotate(nodes[4], "\$\\overline{A}\$", offsetx=-0.3) >> c
        annotate(nodes[5], "\$\\overline{B}\$", offsetx=0.0, offsety=-0.3) >> c
        annotate(nodes[2], "\$\\overline{C}\$", offsety=0.3) >> c
        annotate(nodes[2], "\\small Eq. 6.9", offsety=-0.4) >> c
        PlainText(dx, -2.2, "\\large (b)") >> c

        dx += 4
        nodes = []
        for (x,y) in locs
            push!(nodes, node!(x+dx, -y-dy))
        end
        for (i, j) in [(1, 2), (1,3), (2,4), (2,5), (3,6), (3,7), (4,8), (4,9)]
            Line(nodes[i], nodes[j], arrow="<-", line_width="0.5pt", draw="black") >> c
        end
        annotate(nodes[1], "\$s_{\\alpha(G)}\\infty^{\\alpha(G)}\$", offsety=0.3) >> c
        annotate(nodes[4], "\$A'\\circ\\overline{A}\$", offsetx=-0.45) >> c
        annotate(nodes[5], "\$B'\\circ\\overline{B}\$", offsetx=0.0,offsety=-0.3) >> c
        annotate(nodes[2], "\$C'\\circ\\overline{C}\$", offsety=0.3) >> c
        annotate(nodes[2], "\$*\$", offsety=-0.3) >> c
        PlainText(dx, -2.2, "\\large (c)") >> c
    end
    writepdf(filename, graph)
end

fig1()