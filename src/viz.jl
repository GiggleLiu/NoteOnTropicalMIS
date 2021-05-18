using Viznet
export viztnet, vizconfig
using Compose

vizconfig(ud::UnitdiskWidget, config=zeros(Int, length(ud.nodes))) = vizconfig(ud.nodes, edges(ud), config)
viztnet(ud::UnitdiskWidget, config=zeros(Int, length(ud.nodes))) = viztnet(ud.nodes, edges(ud), config)

function vizconfig(nodes, edges, config=zeros(Int, length(nodes)), unit=1.0, graphsize=12cm)
	tb = textstyle(:default, fill("white"))
	nb = nodestyle(:circle, fill("black"), r=0.02*unit)
	nb2 = nodestyle(:circle, fill("red"),r=0.02*unit)
	eb = bondstyle(:default)
	img = canvas() do
		for (i, (t, p)) in enumerate(nodes)
			(config[i]==1 ? nb2 : nb) >> p
			tb >> (p, t)
		end
		for (i,j) in edges
			eb >> (nodes[i].second, nodes[j].second)
		end
	end
    XMIN = minimum(x->x.second[1], nodes)
    YMIN = minimum(x->x.second[2], nodes)
    XMAX = maximum(x->x.second[1], nodes)
    YMAX = maximum(x->x.second[2], nodes)
    zoom_into(img, XMIN, XMAX, YMIN, YMAX; graphsize=graphsize)
end

function zoom_into(img, XMIN, XMAX, YMIN, YMAX; graphsize)
    sx = 0.8/(XMAX-XMIN)
    sy = 0.8/(YMAX-YMIN)
	Compose.set_default_graphic_size(graphsize*sy/sx, graphsize)
    Compose.compose(context(-((XMIN+XMAX)/2*sx-0.5), -((YMIN+YMAX)/2*sy-0.5), sx, sy), img)
end

function viztnet(nodes, edges, config=zeros(Int, length(nodes)), unit=1.0, graphsize=12cm)
    XMIN = minimum(x->x.second[1], nodes)
    YMIN = minimum(x->x.second[2], nodes)
    XMAX = maximum(x->x.second[1], nodes)
    YMAX = maximum(x->x.second[2], nodes)
	tb = textstyle(:default, fill("blue"))
	nb = nodestyle(:circle, fill("black"),r=0.005*unit, stroke("black"))
	bt = nodestyle(:square, fill("black"),r=0.015*unit)
	nb2 = nodestyle(:circle, fill("red"),r=0.005*unit)
	eb = bondstyle(:default)
	img = canvas() do
		for (i, (t, p)) in enumerate(nodes)
			(config[i]==1 ? nb2 : nb) >> p
			tb >> (p, t)
		end
		for (i,j) in edges
			eb >> (nodes[i].second, nodes[j].second)
            bt >> ((nodes[i].second .+ nodes[j].second) ./ 2)
		end
	end
    zoom_into(img, XMIN, XMAX, YMIN, YMAX; graphsize=graphsize)
end