module OctTree
using StaticArrays
using Healpix
const NSIDE = 1
const NPIX = 12*NSIDE^2


export Node, TreeGather
export buildtree, get_scatter_ngb_tree, get_gather_ngb_tree, treewalk, nearest
export NSIDE, NPIX

mutable struct PartData{N,T}
	pos::SVector{N,T}
	idx::Int64
	hsml::T
	mass::T
	mass_H2::T
	mass_CO::T
end

PartData{N,T}() where {N,T} = PartData{N,T}(zero(SVector{N,T}),0,0,0,0,0)

mutable struct NodeData{N,T} #auxiliary data carried by treenodes
	pos_c::SVector{N,T} #center of mass
	max_hsml::T  #for scatter ngb search
	mass::T
	mass_H2::T
	mass_CO::T
	#idxs::Vector{Int64}
end

NodeData{N,T}() where {N,T} = NodeData{N,T}(zero(SVector{N,T}),0,0,0,0)

mutable struct Node{N,T}
    center::SVector{N,T}
    length::SVector{N,T}
	p::Union{PartData{N,T}, Nothing}           #initilize when an empty leaf becomes nonempty (0->1)
	child::Union{Array{Node{N,T},1}, Nothing}  #initilize when a nonempty leaf becomes a node (1->2)
	n::Union{NodeData{N,T}, Nothing}           #initilize when a nonempty leaf becomes a node (1->2)
end

Node{N,T}(center::SVector{N,T}, length::SVector{N,T}) where {N,T} =
Node{N,T}(center::SVector{N,T}, length::SVector{N,T}, nothing, nothing, nothing)

function isLeaf(node::Node{N,T}) where{N,T}
    return node.child == nothing ? true : false
end

function getChildIndex(pos::SVector{N,T}, node::Node{N,T}) where {N,T}
    idx = 0
    @inbounds for i in 1:N
        if pos[i] > node.center[i]
            idx |= 2^(i-1)
        end
        #println(i, 2^(i-1), idx)
    end
    return idx+1
end

function getoffset(N::T) where {T}
	offset_child = zeros(T,N,2^N)
	for i in 1:2^N
	    a = bitstring(i-1)
	    #offset = SVector(2*parse(Int, a[end])-1, 2*parse(Int, a[end-1])-1, 2*parse(Int, a[end-2])-1)
	    offset = T[]
	    @inbounds for j in 1:N
	        push!(offset,2*parse(T, a[end-(j-1)])-1)
	    end
	    #println(offset)
		offset_child[:,i] = offset
	end
	SMatrix{N,2^N,T}(offset_child)
end

const offset = getoffset(3) #precomputed offset array



function insertpart!(p::PartData{N,T}, node::Node{N,T}) where {N,T}
    if isLeaf(node)
        if node.p == nothing
            #println("at an empty leaf, insert particle and we're done")
			@assert node.n == nothing
			node.p = p
        else
            #println("at a leaf with a preexisting particle, so we split it")
			#initializing the child nodes and calculate their centers
            node.child = Array{Node{N,T},1}(undef,2^N)
            @inbounds for i in 1:2^N
                #println(getoffset(i,N))
                childcenter = node.center + offset[:,i] * 0.25 .* node.length
                #println("childcenter=", childcenter, "  getoffset(i,N)=", getoffset(i,N))
				node.child[i] = Node{N,T}(childcenter, 0.5 .* node.length)
                #println("center = ", childcenter, "  edge_min = ", childcenter - 0.25*node.length,
                #        "  edge_max = ", childcenter + 0.25*node.length)
            end
			#println("insert back the preexsisting particle...")
			insertpart!(node.p, node.child[getChildIndex(node.p.pos, node)])
			#since it's not a leaf node anymore, initialize node data
			node.n = NodeData{N,T}()
			#push!(node.n.idxs, node.p.idx, p.idx)
			masssum = node.p.mass + p.mass
			inv_masssum = 1.0 / masssum
			node.n.pos_c = (node.p.mass .* node.p.pos .+ p.mass .* p.pos) .* inv_masssum
            node.n.mass = masssum #this line has to be put after insertpart!, otherwise we're inserting back the wrong node.mass!!!
			node.n.mass_H2 = node.p.mass_H2 + p.mass_H2
			node.n.mass_CO = node.p.mass_CO + p.mass_CO
			node.n.max_hsml = node.p.hsml > p.hsml ? node.p.hsml : p.hsml
            #println("insert the new particle...")
            insertpart!(p, node.child[getChildIndex(p.pos, node)])
        end
    else
        #println("open node")
		#push!(node.n.idxs, p.idx)
		masssum = node.n.mass + p.mass
		inv_masssum = 1.0 / masssum
		node.n.pos_c = (node.n.mass .* node.n.pos_c .+ p.mass .* p.pos) .* inv_masssum
        node.n.mass = masssum
		node.n.mass_H2 += p.mass_H2
		node.n.mass_CO += p.mass_CO
		node.n.max_hsml = node.n.max_hsml > p.hsml ? node.n.max_hsml : p.hsml
        insertpart!(p, node.child[getChildIndex(p.pos,node)])
    end
end

#PBC wrap: typically used as nearest.(dx) where dx is a separation vector connecting two points
function nearest(x::T, boxsize::T) where {T}
	#can be negative as this is just a vector component
    return (x > 0.5*boxsize) ? (x - boxsize) : ( (x < -0.5*boxsize) ? (x + boxsize) : x )
end

function nearest_plain(x::T) where {T}
    if x < -boxHalf_X
        x = x + boxHalf_X*2
    else
        x = x
    end
    if x > boxHalf_X
        x = x - boxHalf_X*2
    else
        x = x
    end
    return x
end

function get_distance2(a::SVector{N,T}, b::SVector{N,T}, boxsizes::SVector{N,T}, periodic::Bool) where {N,T}
    c = a - b
    if periodic
        c = nearest.(c, boxsizes)
    end
    return sum(c .^ 2)
end

res = Healpix.Resolution(NSIDE)

function vec2pix(p)
    theta,phi = vec2ang(p[1],p[2],p[3])
    ang2pixRing(res, theta,phi)
end

mutable struct TreeGather{T}
    mass::T
    column_all::Vector{T}
    column_H2::Vector{T}
    column_CO::Vector{T}
    #nodecenters::Vector{Vector{T}} #for debugging/visualization
    #nodelengths::Vector{Vector{T}} #for debugging/visualization
end

TreeGather{T}() where {T} = TreeGather{T}(0,zeros(NPIX),zeros(NPIX),zeros(NPIX))

#when changing the argument of treewalk, remember to also change the recursive treewalk call inside the function body!
function treewalk(ga::TreeGather{T}, p::SVector{N,T}, node::Node{N,T},
	openingangle::T, shieldinglength::T, boxsizes::SVector{N,T}) where {N,T}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p != nothing
            #println("nonempty leaf")
            ga.mass += node.p.mass
            dx = nearest.(node.p.pos - p, boxsizes)
			dist2 = sum(dx.^2)
			#@show dist2, node.p.pos, p
			if 0.0 < dist2 < shieldinglength^2  #exclude self contribution (dist2=0)
	            ipix = vec2pix(dx)
	            #ga.column_all[ipix] += node.mass
	            area = 4*pi/NPIX * dist2
				ga.column_all[ipix] += node.p.mass / area
				ga.column_H2[ipix] += node.p.mass_H2 / area
	            ga.column_CO[ipix] += node.p.mass_CO / area
	            #push!(ga.nodecenters, node.p.pos)
	            #push!(ga.nodelengths, zeros(N))
			end
            #@show "it's a particle..." node.center, node.length
            #dist2 = get_distance2(node.part, p, boxsizes, true)
        end
    else
        #println("This is a node... check the opening angle")
		#dist2 = get_distance2(pos_c, p, boxsizes, true)
		#don't use get_distance2 as we need both dist2 and dx (for vec2pix)
		dx = nearest.(node.n.pos_c - p, boxsizes)
		dist2 = sum(dx.^2)
		if dist2 > (node.length[1] / openingangle)^2
			#println("skip node ", i)
			if dist2 < shieldinglength^2
				ga.mass += node.n.mass
				#dx = nearest.(node.n.pos_c - p, boxsizes)
				ipix = vec2pix(dx)
				#ga.column_all[ipix] += node.mass
				area = 4*pi/NPIX * dist2
				ga.column_all[ipix] += node.n.mass / area
				ga.column_H2[ipix] += node.n.mass_H2 / area
				ga.column_CO[ipix] += node.n.mass_CO / area
				#push!(ga.nodecenters, node.center)
				#push!(ga.nodelengths, node.length)
			end
			#@show "use this node!" node.center, node.length
		else
	        @inbounds for i in 1:2^N
	            #println("open this node")
	            treewalk(ga, p, node.child[i], openingangle, shieldinglength, boxsizes)
	        end
		end
    end
end

#ngb treewalk: recursively find gather ngbs for particle p within the searching radius
function get_ptl_in_box_gather(radius::T, p::SVector{N,T}, node::Node{N,T}, boxsizes::SVector{N,T}, idx_ngbs::Vector{Int64}) where {N,T}
    if isLeaf(node)
        #println("in a leaf node")
        if node.p != nothing
            #println("nonempty leaf")
	    	dist2 = get_distance2(node.p.pos, p, boxsizes, true)
            if dist2 < radius^2
                #println("push ", node.part)
				push!(idx_ngbs, node.p.idx)
            end
        end
    else
        #println("This is a node... check if node and radius overlap")
		@inbounds for j in 1:N
			if abs(nearest(node.center[j] - p[j], boxsizes[j])) > 0.5*node.length[j] + radius
				return  #no overlap, so node can be skipped!
			end
		end
		#=
		if get_distance2(node.center, p, boxsizes, true) + 0.75 * node.length[1]^2 < radius^2
			append!(idx_ngbs, node.n.idxs)
			return
		end
		=#
		#println("node and radius overlap, so we need to open the node and check its children")
        @inbounds for i in 1:2^N
			get_ptl_in_box_gather(radius, p, node.child[i], boxsizes, idx_ngbs)
        end
    end
end

#ngb treewalk: recursively find scatter ngbs for particle p
function get_ptl_in_box_scatter(p::SVector{N,T}, node::Node{N,T}, boxsizes::SVector{N,T}, idx_ngbs::Vector{Int64}) where {N,T}
    if isLeaf(node)
        #println("in a leaf node")
		if node.p != nothing
            #println("nonempty leaf")
	    	dist2 = get_distance2(node.p.pos, p, boxsizes, true)
            if dist2 < node.p.hsml^2
                #println("push ", node.part)
				push!(idx_ngbs, node.p.idx)
            end
        end
    else
		#println("This is a node... check if node and radius overlap")
		@inbounds for j in 1:N
			if abs(nearest(node.center[j] - p[j], boxsizes[j])) > 0.5*node.length[j] + node.n.max_hsml
				return  #no overlap, so node can be skipped!
			end
		end
		#println("node and radius overlap, so we need to open the node and check its children")
        @inbounds for i in 1:2^N
			get_ptl_in_box_scatter(p, node.child[i], boxsizes, idx_ngbs)
        end
    end
end

#driver function to get gather ngbs for particle x within the searching radius h
function get_gather_ngb_tree(x::SVector{N,T}, h::T, node::Node{N,T}, boxsizes::SVector{N,T}) where {N,T}
    idx_ngbs = Int64[]
    get_ptl_in_box_gather(h, x, node, boxsizes, idx_ngbs)
    return idx_ngbs
end

#driver function to get scatter ngbs for particle x; no need to specify a searching radius
function get_scatter_ngb_tree(x::SVector{N,T}, node::Node{N,T}, boxsizes::SVector{N,T}) where {N,T}
    idx_ngbs = Int64[]
    get_ptl_in_box_scatter(x, node, boxsizes, idx_ngbs)
    return idx_ngbs
end

#build a tree for given particle locations X
#we distinguish between topnode_length (for the tree) and boxsizes (for PBC wrap)
function buildtree(X::Vector{SVector{N,T}}, hsml::Vector{T}, mass::Vector{T},
	mass_H2::Vector{T}, mass_CO::Vector{T}, center::SVector{N,T}, topnode_length::SVector{N,T}) where {N,T}

    #construct the tree for ngb search
	tree = Node{N,T}(center, topnode_length)

    @inbounds for i in eachindex(X)
		part = PartData(SVector(X[i]), i, hsml[i], mass[i], mass_H2[i], mass_CO[i])
		insertpart!(part, tree)
    end
    return tree
end

end #module OctTree
